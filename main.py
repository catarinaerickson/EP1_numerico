##############################################################################
#                               EP1 - MAP3121                                #
#       Autovalores e Autovetores de Matrizes Tridiagonais Simétricas        #
##############################################################################
#        Alunas:                                                             #
#                Catarina Rodrigues Erickson   (NUSP: 11258742)              #
#                Isabella Mulet e Souza        (NUSP: 11259184)              #
#                                                                            #
##############################################################################

#Importando bibliotecas para trabalhar com aritmética de vetores e gráficos
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Cabeçalho
    print("EP1 - MAP3121")
    print("Autovalores e Autovetores de Matrizes Tridiagonais Simétricas")
    print("Catarina Erickson e Isabella Mulet")
    print()

# Calcula u de acordo com a heurística de Wilkinson
def calc_heuristica(alpha_arr, beta_arr, n):
    d = (alpha_arr[n-1] - alpha_arr[n])/2
    if d>=0:
        sgn = 1
    else:
        sgn = -1
    u = alpha_arr[n] + d - sgn*((d**2 + beta_arr[n-1]**2)**(1/2))
    return u

# Calcula o seno e o cosseno da matriz Q de transformação de Givens
def calc_theta(alpha, beta):
    raiz = (alpha**2+beta**2)**(1/2)
    cos = alpha/raiz
    sen = -beta/raiz
    return (cos, sen)

# Algoritmo QR
# Recebe uma lista da matriz diagonal principal, uma lista da matriz subdiagonal 
# e True ou False para o deslocamento espectral.
# Retorna os autovalores, a mariz de autovetores e o númeor de iterações
def qr(alpha_list, beta_list, deslocamento):
    alpha_arr = np.array(alpha_list, dtype = float)
    beta_arr = np.array(beta_list, dtype = float)

    # Define a dimensão inicial da matriz Anxn
    dim = alpha_arr.shape[0]

    # Define V inicial (matriz identidade)
    v_arr = np.empty((dim,dim))
    for i in range(dim):
        v = np.zeros(dim)
        v[i] = 1
        v_arr[i] = v

    # Inicializa vetor para armazenar os autovalores encontrados
    autovalor_arr = np.array([])
    k=0
    # Itera para m = dim, dim-1, ..., 2
    for m in range(dim, 1, -1):
        
        # Itera enquanto o beta de referência da iteração m for maior do que 10e-6
        while abs(beta_arr[m-2])>=0.000001:
            # Calcula e subtrai ukI de A(k) a partir da iteração k caso haja deslocamento
            if k>0 and deslocamento:
                u = calc_heuristica(alpha_arr, beta_arr, m-1)
                u_arr = u*np.ones(alpha_arr.shape)
                alpha_arr = alpha_arr - u_arr

            # Inicializa os vetores para armazenar os valores de seno e cosseno da iteração k
            c_arr = np.array([])
            s_arr = np.array([])


            # Calcula a matriz R
            beta_fixed = beta_arr.copy()
            for i in range(m-1):
                # Calcula o seno e cosseno da iteração i e guarda nos respectivos vetores
                (cos, sen) = calc_theta(alpha_arr[i],beta_fixed[i])
                c_arr = np.append(c_arr, cos)
                s_arr = np.append(s_arr, sen)

                # Cria a diagonal principal e a subdiagonal da matriz Qi*Qi-1*...*Q1*A
                alpha_R = alpha_arr.copy()
                beta_R = beta_arr.copy()

                # Modifica a linha i
                alpha_R[i] = cos*alpha_arr[i] - sen*beta_fixed[i]
                beta_R[i] = cos*beta_arr[i] - sen*alpha_arr[i+1]
                # Modifica a linha i+1
                alpha_R[i+1] = sen*beta_arr[i] + cos*alpha_arr[i+1]
                # Modifica a linha i+1 de beta apenas se não for a última iteração i=m-2
                if i < beta_arr.shape[0]-1:
                    beta_R[i+1] = cos*beta_arr[i+1]

                # Define as diagonais da matriz Qi*Qi-1*...*Q1*A como ponto de partida para a próxima iteração
                alpha_arr = alpha_R.copy()
                beta_arr = beta_R.copy()

            # Calcula A(k+1)
            for i in range(m-1):
                # Cria a diangonal principal e a subdiagonal da matriz R*Q1T*...*QiT
                alpha_next = alpha_arr.copy()
                beta_next = beta_arr.copy()

                # Modifica a coluna i
                alpha_next[i] = c_arr[i]*alpha_arr[i] - s_arr[i]*beta_arr[i]
                # Modifica a coluna i+1
                beta_next[i] = -s_arr[i]*alpha_arr[i+1]
                alpha_next[i+1] = c_arr[i]*alpha_arr[i+1]

                # Define as diagonais da matriz R*Q1T*...*QiT como ponto de partida para a próxima iteração
                alpha_arr = alpha_next.copy()
                beta_arr = beta_next.copy()

            # Soma A(k+1) e ukI a partir da iteração k caso haja deslocamento
            if k>0 and deslocamento:
                alpha_arr = alpha_arr + u_arr 

            # Calcula V
            for i in range(m-1):
                v_next = v_arr.copy()
                # Altera a colunas i e i+1 para cada linha j da matriz V
                for j in range(dim):
                    v_next[j][i] = v_arr[j][i]*c_arr[i] - v_arr[j][i+1]*s_arr[i]
                    v_next[j][i+1] = v_arr[j][i]*s_arr[i] + v_arr[j][i+1]*c_arr[i]

                # Define as diagonais da matriz V*Q1T*...*QiT como ponto de partida para a próxima iteração
                v_arr = v_next.copy()
            k+=1   
        
        # Guarda o autovalor encontrado na iteração m no vetor de autovalores
        autovalor_arr = np.append(autovalor_arr, alpha_arr[m-1])

        # Define as diagonais da submatriz para a iteração m-1
        alpha_arr = np.delete(alpha_arr, m-1)
        beta_arr = np.delete(beta_arr, m-2)

    # Guarda o autovetor que sobrou depois da última iteração m-1
    autovalor_arr = np.append(autovalor_arr, alpha_arr[0])

    # Inverte o vetor dos autovalores para que o autovalor da posição i corresponda ao autovetor da coluna i da matriz V
    autovalor_arr = np.array(list(reversed(autovalor_arr)))

    # Retorna o vetor com os autovalores, a matriz V e o número de iterações, nessa ordem
    return autovalor_arr, v_arr, k

# Aplica a função qr() em uma matriz tridiagonal simétrica de ordem n
# com diagonal principal constante igual a 2 e subdiagonal igual a -1
def teste_qr(n, deslocamento):
    alpha_arr = 2*np.ones((n,))
    beta_arr = -1*np.ones((n-1,))   
    return qr(alpha_arr, beta_arr, deslocamento)

# Recebe uma matriz e retorna a matriz transposta correspondente
def transpose(mat):
    lin = mat.shape[0]
    col = mat.shape[1]

    mat_t = np.zeros((col,lin))

    for i in range(lin):
        for j in range(col):
            mat_t[j][i] = mat[i][j]

    return mat_t

def teste_analitico(n):    
    autovalor_analitico = np.array([])
    autovetor_transp = np.zeros((n,n))
    for i in range(1,n+1):
        lmb = 2*(1-np.cos(i*np.pi/(n+1)))
        v = 0
        norma = 0
        for j in range (1,n+1):
            v = np.sin(j*i*np.pi/(n+1))
            norma += (v)**2
            autovetor_transp[i-1][j-1] = v
        autovetor_transp[i-1] = autovetor_transp[i-1]/(norma**(1/2))
        autovalor_analitico = np.append(autovalor_analitico, lmb)
    autovetor_analitico = transpose(autovetor_transp)
    return autovalor_analitico, autovetor_analitico

# Cria a matriz A tridiagonal simétrica do sistema massa-mola com n massas
def create_A(n):
    A_principal = np.array([])
    A_subdiagonal = np.array([])

    for i in range(1,n+1):
        # Calcula os elementos da diagonal principal com a fórmula dada para 5 ou 10 massas
        if n==5:
            p = ((40+2*i)+(40+2*(i+1)))/2
        else:
            p = ((40+2*(-1)**i)+(40+2*(-1)**(i+1)))/2
        A_principal = np.append(A_principal, p)

        # Calcula os elementos da subdiagonal com a fórmula dada para 5 ou 10 massas
        if i>1:
            if n==5:
                s = -(40+2*i)/2
            else:
                s = -(40+2*(-1)**i)/2
            A_subdiagonal = np.append(A_subdiagonal, s)

    return A_principal, A_subdiagonal

# Calcula a coordenada(deslocamento) para uma abscissa(tempo) da função xi(t) de cada mola do sistema,
# incluindo a mudança de variável de Y(t) para X(t)
def calc_coord(Y_zero, Q_autovetor, A_autovalor, lin, absc):
    coord = float(0)
    for col in range(Q_autovetor.shape[1]):
        coord += Q_autovetor[lin][col]*Y_zero[col]*np.cos(absc*(A_autovalor[col])**(1/2))
    return coord

# Plota o gráfico dos deslocamento das massas de um sistema massa mola, dado o vetor de deslocamentos
# iniciais X(0) e a diagonal principal e a subdiagonal da matriz A  
def plot_desl(X_inicial, A_principal, A_subdiagonal):
    # Define o vetor de deslocamentos iniciais X(0)
    X_zero = np.array(X_inicial, dtype = float)

    # Define os autovalores e a matriz de autovetores de A
    A_autovalor, Q_autovetor, k = qr(A_principal, A_subdiagonal, True)

    # Realiza a mudança de variável de X(0) para Y(0), em que Y(t)=QT*X(t)
    Y_zero = np.dot(transpose(Q_autovetor), X_zero)

    # Define as abscissas do gráfico de 0 a 10 segundos a cada 0.025 segundos
    absc = np.linspace(0, 10, 401)
    # Calcula as coordenadas da função xi(t) de cada mola do sistema
    for i in range(A_autovalor.shape[0]):
        plt.plot(absc, calc_coord(Y_zero, Q_autovetor, A_autovalor, i, absc), label=f"x{i+1}")

    # Define os nomes do eixos, o título e o posicionamento da legenda no gráfico
    plt.xlabel('tempo (s)')
    plt.ylabel('deslocamento')
    plt.title("Sistema massa mola: deslocamento ao longo do tempo")
    plt.legend(loc=(1.04,0))
    plt.tight_layout() 

    plt.show()

# Plota o gráfico dos deslocamento das massas de um sistema massa mola com os deslocamentos iniciais
# correspondentes ao modo de maior frequência
def plot_desl_freq(A_principal, A_subdiagonal):
    # Define os autovalores e a matriz de autovetores de A
    A_autovalor, Q_autovetor, k = qr(A_principal, A_subdiagonal, True)

    # Calcula X_zero como o autovetor associado ao maior autovalor lambda
    lmb, X_zero = select_X_zero(A_autovalor, Q_autovetor)
    
    # Realiza a mudança de variável de X(0) para Y(0), em que Y(t)=QT*X(t)
    Y_zero = np.dot(transpose(Q_autovetor), X_zero)

    # Define o vetor com o valor da maior frequência ao quadrado (freq = raiz de lambda)
    lmb_arr = lmb*np.ones((A_autovalor.shape))
    # Define as abscissas do gráfico de 0 a 10 segundos a cada 0.025 segundos
    absc = np.linspace(0, 10, 401)
    # Calcula as coordenadas da função xi(t) de cada mola do sistema
    for i in range(A_autovalor.shape[0]):
        plt.plot(absc, calc_coord(Y_zero, Q_autovetor, lmb_arr, i, absc), label=f'x{i+1}')

    # Define os nomes do eixos, o título e o posicionamento da legenda no gráfico
    plt.xlabel('tempo (s)')
    plt.ylabel('deslocamento')
    plt.title("Sistema massa mola: modo de maior frequência")
    plt.legend(loc=(1.04,0))
    plt.tight_layout()

    plt.show()

# Retorna o maior autovalor e seu autovetor correspondente, dados o vetor com os autovalores e 
# a matriz de autovetores
def select_X_zero(mat_autovalor, mat_autovetor):
    lmb = mat_autovalor[0]
    X_zero = np.array([])
    pos = 0
    
    # Define o maior autovalor e sua posição
    for i in range(1, mat_autovalor.shape[0]):
        if mat_autovalor[i]>lmb:
            lmb = mat_autovalor[0]
            pos = i

    # Define o autovetor correspondente
    for i in range(mat_autovalor.shape[0]):
        X_zero = np.append(X_zero, mat_autovetor[i][pos])

    print(X_zero)

    return lmb, X_zero

#Chama a função main()
if __name__ == "__main__":
    main()