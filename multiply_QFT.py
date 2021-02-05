from qiskit import *
import matplotlib.pyplot as mpl
import numpy as np
from numpy import pi
import math
from qiskit.visualization import plot_histogram
#初始化
def init_reg(num1,num2):
    for index in range(l1):
        if num1[index]=='1':
            qc.x(qra[l1-1-index])
    for index in range(l2):
        if num2[index]=='1':
            qc.x(qrb[l2-1-index])



#量子傅里叶变换
def qft_rotations(circuit, n):
    """Performs qft on the first n qubits in circuit (without swaps)"""
    if n==l3-l1-1:
        return circuit
    n -=1
    circuit.h(qrsum[n])
    for qubit in range(n-1,l3-l1-2,-1):
        circuit.cu1(pi/2**(n-qubit), qrsum[qubit], qrsum[n])

    qft_rotations(circuit, n)
    return circuit


#加法
def phase_rotations(circuit,l2_i):
    for qubitl3 in range(l3-1,l3-l1-2,-1):
        k=(qubitl3-(l3-l1-1))+1
        for qubitl1 in range(l1):
            if qubitl1>(qubitl3-(l3-l1-1)):
                break
            circuit.ccx(qra[qubitl1],qrb[l2_i],d[0])
            circuit.cu1(pi / 2 ** (k - 1), d[0], qrsum[qubitl3])
            circuit.reset(d[0])
            k=k-1
            # circuit.cu1(pi / 2 ** (k - 1), qra[qubitl1], qrsum[qubitl3])
    return circuit

#逆QFT变换
def qft_dagger(circuit, n):
    '''for qubit in range(n//2):
        circ.swap(qubit,n-qubit-1)'''
    for j in range(l3-l1-1,l3):
        for m in range(l3-l1-1,j):
            circuit.cu1(-math.pi/float(2**(j-m)),qrsum[m],qrsum[j])
        circuit.h(qrsum[j])
    return circuit

#右移
def move_right(circuit,l3):
    for j in range(l3):
        if j==0:
            circuit.cx(qrsum[0],d[0])
            circuit.cx(d[0],qrsum[0])
        else:
            circuit.cx(qrsum[j],qrsum[j-1])
            circuit.cx(qrsum[j-1],qrsum[j])
    circuit.reset(d[0])
    return circuit

try:
    a = int(input("Input binary value: "), 2)
    print("num (decimal format):", a)
except ValueError:
    print("Please input only binary value...")
try:
    b = int(input("Input binary value: "), 2)
    print("num (decimal format):", b)
except ValueError:
    print("Please input only binary value...")


if b>a:
    a,b = b,a
num1=bin(a)[2:]
num2=bin(b)[2:]
l1=len(num1)
l2=len(num2)
l3=l1+l2
qra=QuantumRegister(l1,name='a')   #被乘数
qrb=QuantumRegister(l2,name='b')   #乘数
qrsum=QuantumRegister(l3,name="sum")
d=QuantumRegister(1)
test=ClassicalRegister(l3,name='test')
qc=QuantumCircuit(qra,qrb,qrsum,d,test)
init_reg(num1,num2)
qc.barrier()
for i in range(l2):
    qft_rotations(qc,l3)
    qc.barrier()
    phase_rotations(qc,i)
    qc.barrier()
    qft_dagger(qc,l3)
    qc.barrier()
    if i!=l2-1:
        move_right(qc,l3)
        qc.barrier()

for qubit in range(l3):
    qc.measure(qrsum[qubit],test[qubit])
print(qc)
backend = Aer.get_backend('qasm_simulator')
shots=1024
results = execute(qc, backend=backend, shots=shots, optimization_level=0).result()
counts = results.get_counts()

plot_histogram(counts)
mpl.show()