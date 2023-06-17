import cirq

from typing import List

NBITS = 3

def oracle(qubits: List[cirq.LineQubit], ancillas: List[cirq.NamedQubit]):
   XPRIME = [0,1,0]
   q0, q1, q2 = qubits
   a0, a1 = ancillas
   yield (cirq.X(q) for (q,bit) in zip(qubits, XPRIME) if not bit)
   yield cirq.CCNOT(q0, q1, a0)
   yield cirq.CCNOT(a0, q2, a1)
   yield cirq.Z(a1)
   yield cirq.CCNOT(q0, q1, a0)
   yield cirq.CCNOT(a0, q2, a1)
   yield (cirq.X(q) for (q,bit) in zip(qubits, XPRIME) if not bit)


def grover_operator(qubits: List[cirq.LineQubit]):
   yield cirq.H.on_each(*qubits)
   yield cirq.X.on_each(*qubits)
   yield cirq.CCZ(*qubits)
   yield cirq.X.on_each(*qubits)
   yield cirq.H.on_each(*qubits)


def grover():
   qubits = cirq.LineQubit.range(NBITS)
   ancillas = [cirq.NamedQubit('Ancilla0'), cirq.NamedQubit('Ancilla1')]
   circuit = cirq.Circuit()

   circuit.append(cirq.H.on_each(*qubits))
   for _ in range(2):
      circuit.append(oracle(qubits, ancillas))
      circuit.append(grover_operator(qubits))
   circuit.append(cirq.measure(*qubits, key="result"))

   # Sample from the circuit a couple times.
   simulator = cirq.Simulator()
   result = simulator.run(circuit, repetitions=100)

   # Look at the sampled bitstrings.
   frequencies = result.histogram(key="result", fold_func= lambda bits: "".join(str(int(b)) for b in bits))
   print('Sampled results:\n{}'.format(frequencies))

   # Check if we actually found the secret value.
   most_common_bitstring = frequencies.most_common(1)[0][0]
   print("\nMost common bitstring: {}".format(most_common_bitstring))

