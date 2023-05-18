from math import cos, isclose, radians, sin, sqrt
from typing import Any, Dict, List

from pydantic import BaseModel, root_validator

class C(BaseModel):
    """ Complex number as re + i*im """
    re: float # real part
    im: float # imaginary part
    
    @property
    def mag(self) -> float:
        """ Magnitude """
        return sqrt(self.re**2 + self.im**2)

    @property
    def prob(self) -> float:
        """ Probability of colapsing to this state (mag**2). """
        return self.mag**2

    def __str__(self) -> str:
        if self.im == 0.0:
            return str(self.re)
        elif self.re == 0.0:
            return f'{self.im}i'
        else:
            return f'({self.re} + {self.im}i)'

    def __add__(self, x: 'C') -> 'C':
        if isinstance(x, C):
            return C(re=self.re+x.re, im=self.im+x.im)
        else:
            raise ValueError('huet addition: {x}')

    def __sub__(self, x: 'C') -> 'C':
        if isinstance(x, C):
            return C(re=self.re-x.re, im=self.im-x.im)
        else:
            raise ValueError('huet addition: {x}')

    def __mul__(self, x: [int, float, 'C']) -> 'C':
        if isinstance(x, C):
            return C(
                re=self.re*x.re + self.im*x.im,
                im=self.re*x.im + self.im*x.re,
            )
        else:
            return C(re=x*self.re, im=x*self.im)


class Q(BaseModel):
    """ Qubit as amp0*|0> + amp1*|1>  """
    amp0: C # Amplitude of 0
    amp1: C # Amplitude of 1

    def __str__(self) -> str:
        return f'{self.amp0}*|0> + {self.amp1}*|1>'

    @root_validator()
    def _has_prob_1(cls, values: Dict[str, Any]) -> Dict[str, Any]:
        """ Assert the probability of collapsing to 0 or 1 == 1. """
        prob = values['amp0'].prob + values['amp1'].prob
        if not isclose(prob, 1.0):
            raise ValueError(f'Probability is not 1: {values} {prob}')
        return values


# Constants
MS2 = 1/sqrt(2) # Square root of 2 to the minus 1
c0 = C(re=0,im=0) # Complex 0
c1 = C(re=1,im=0) # Complex 1
q0 = Q(amp0=c1, amp1=c0) # |0> aka [1, 0]T
q1 = Q(amp0=c0, amp1=c1) # |1> aka [0, 1]T


def _assert_0_or_1(qubits: List[Q]) -> None:
    for q in qubits:
        if not (q == q0 or q == q1):
            raise ValueError(f'Qubit must be |0> or |1>: {q}')


def _round_qubit(q: Q) -> None:
    """ Bring values close to 0 or 1 to 0 or 1. """
    def _round(f: float) -> float:
        if isclose(f, 0, abs_tol=1e-9):
            return 0
        elif isclose(f, 1):
            return 1
        else:
            return f
    q.amp0.re = _round(q.amp0.re)
    q.amp0.im = _round(q.amp0.im)
    q.amp1.re = _round(q.amp1.re)
    q.amp1.im = _round(q.amp1.im)

def _apply_UT(q: Q, ut) -> None:
    """ Apply a 2x2 unitary transformation to qubit q. """
    print(q)
    amp0 = q.amp0 * ut[0][0] + q.amp1 * ut[0][1]
    amp1 = q.amp0 * ut[1][0] + q.amp1 * ut[1][1]
    q.amp0 = amp0
    q.amp1 = amp1
    _round_qubit(q)
    print(q)


# GATES
def g_NOT(q: Q) -> None:
    """ Reverse the amplitudes. """
    amp0 = q.amp0.copy()
    q.amp0 = q.amp1.copy()
    q.amp1 = amp0


def g_HADAMARD(q: Q) -> None:
    amp0 = (q.amp0 + q.amp1) * MS2
    amp1 = (q.amp0 - q.amp1) * MS2
    q.amp0 = amp0
    q.amp1 = amp1


def g_CNOT(control: Q, target: Q) -> None:
    """ Controlled NOT """
    _assert_0_or_1([control, target])
    if control == q1:
        g_NOT(target)


def g_CCNOT(cc: Q, c: Q, t: Q) -> None:
    """ Controlled CNOT """
    _assert_0_or_1([cc, c, t])
    if cc == q1:
        g_CNOT(c, t)


def g_SWAP(a: Q, b: Q) -> None:
    _assert_0_or_1([a, b])
    bak = a.copy()
    a.amp0 = b.amp0
    a.amp1 = b.amp1
    b.amp0 = bak.amp0
    b.amp1 = bak.amp1
    

def g_CSWAP(c: Q, a: Q, b: Q) -> None:
    _assert_0_or_1([c, a, b])
    if c == q1:
        g_SWAP(a, b)


def g_ROT(q: Q, phi: float) -> None:
    """ Rotate by phi degrees counter-clockwise. """
    rad_phi = radians(phi)
    ut = [
        [cos(rad_phi), -sin(rad_phi)],
        [sin(rad_phi), cos(rad_phi)],
    ]
    _apply_UT(q, ut)


def q_S(q: Q) -> None:
    """ Multiply amp1 by i """
    ut = [
        [1, 0],
        [0, C(re=0, im=1)],
    ]
    _apply_UT(q, ut)


def g_T(q: Q) -> None:
    ut = [
        [1, 0],
        [0, C(re=MS2, im=MS2)],
    ]
    _apply_UT(q, ut)


def g_T_dgr(q: Q) -> None:
    """ T dagger gate """
    ut = [
        [1, 0],
        [0, C(re=MS2, im=-MS2)],
    ]
    _apply_UT(q, ut)


def g_Z(q: Q) -> None:
    """ Multiply amp1 by -1 """
    ut = [
        [1, 0],
        [0, -1],
    ]
    _apply_UT(q, ut)
    

