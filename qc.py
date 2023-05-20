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
        """ Probability of colapsing to this state (= mag**2). """
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

    @property
    def conj(self) -> 'C':
        """ Get conjugate. """
        return C(re=self.re, im=self.im * -1)


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

    @property
    def vec(self) -> List[C]:
        """ Qbit to 2-dim vector. """
        return [self.amp0, self.amp1]


# Constants
MS2 = 1/sqrt(2) # Square root of 2 to the minus 1
C0 = C(re=0,im=0) # Complex 0
C1 = C(re=1,im=0) # Complex 1
Q0 = Q(amp0=C1.copy(), amp1=C0.copy()) # |0> aka [1, 0]T
Q1 = Q(amp0=C0.copy(), amp1=C1.copy()) # |1> aka [0, 1]T


def _assert_0_or_1(qubits: List[Q]) -> None:
    for q in qubits:
        if not (q == Q0 or q == Q1):
            raise ValueError(f'Qubit must be |0> or |1>: {q}')


def gen_eye(n: int) -> List[List[float]]:
    """ Make identity matrix. """
    eye = [[C0.copy() for _ in range(n)] for _ in range(n)]
    for i in range(n):
        eye[i][i] = C1.copy()
    return eye
            

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


class SuPos:
    """ Super position """
    def __init__(self, tensor: List[C]) -> None:
        self.reg = tuple(t for t in tensor)


def scale_mat(mat: 'Matrix', s: float) -> 'Matrix':
    """ Scale a matrix by a scalar <s>. """
    for row in mat:
        for j in range(len(row)):
            row[j] *= s
    return mat


def mat_dagger(mat: 'Matrix') -> 'Matrix':
    """ Get hermitian (conjugate) transpose of mat. """
    dagger = []
    for col_ix in range(len(mat[0])):
        new_row = [row[col_ix] for row in mat]
        print(new_row)
        dagger.append(new_row)
    return dagger
    

def join_mats(mats: List['Matrix'], axis: str) -> 'Matrix':
    ''' Join matrices horizontally or vertically. '''
    if axis == 'h':
        row_nums = set(len(mat) for mat in mats)
        if len(row_nums) > 1:
            raise ValueError('All matrices must have the same numnber of rows.')

        num_rows = list(row_nums)[0]
        joined = [[] for _ in range(num_rows)]
        for mat in mats:
            for row_ix in range(num_rows):
                joined[row_ix].extend(mat[row_ix])
        return joined
    elif axis == 'v':
        #V join is the same as: dagger --> H join --> dagger.
        daggered = join_mats([mat_dagger(mat) for mat in mats], 'h') 
        return mat_dagger(daggered)
    else:
        raise ValueError('Unsupported axis: {axis}')


def tensor_mats(mats: List['Matrix']) -> 'Matrix':
    """ Make a tensor product of matrices. """


def tensor_vecs(vectors: List['Vector']) -> 'Vector':
    """ Make a tensor product of vectors. """
    if len(vectors) < 2:
        return vectors
    elif len(vectors) > 2:
        sub = tensor_vecs(vectors[1:])
        return tensor_vecs([vectors[0], sub])
    else:
        tensor = []
        for a0 in vectors[0]:
            for a1 in vectors[1]:
                tensor.append(a0 * a1)
        return tensor







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
    if control == Q1:
        g_NOT(target)


def g_CCNOT(cc: Q, c: Q, t: Q) -> None:
    """ Controlled CNOT """
    _assert_0_or_1([cc, c, t])
    if cc == Q1:
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
    if c == Q1:
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
    

