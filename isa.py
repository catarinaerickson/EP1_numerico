#item b EP1 Numérico

import numpy as np
def k(x): 
    k = a + b*x
    return(k)

m = int(input("Qual a massa, em Kg? "))
n = int(input("Quantas massas o sistema apresenta? "))
a = int(input("O termo livre da função k: "))
b = int(input("O termo acompanhado de i da função k:")) 
A = []

#Criação da matriz A
A = []
for i in range(n):
   A.append( [0] * n )
i=0
j=0
for i in range (n):
    for j in range (n):
        if i==j:
            A[i][j]= (i+1)+(j+1)
        if j==i+1:
            A[i][j]=-(j+1)
        if j==i-1:
            A[i][j]=-(i+1)
        j = j+1 
    j=0
    i=i+1
    
#Falta multiplicar a matriz A por 1/m 
#Fazer a decomposição QR para A e encontrar V (matriz com autovetores de A)
#Nesses itens, Q^T é a matriz de autovetores de A

X0 = np.array([[-2],[-3],[-1],[-3],[-1]])
Qt = V.T
Y0 = np.dot (Qt,X0)





