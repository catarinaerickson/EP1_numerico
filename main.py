import numpy as np

print('EP1 - MAP3121 \nAutovalores e Autovetores de Matrizes Tridiagonais Simétricas')
print('Alunas: Catarina e Isabella')

#Calcula u de acordo com a heurística de Wilkinson
def calc_heuristica(alpha_arr, beta_arr, n):
    d = (alpha_arr[n-1] - alpha_arr[n])/2
    if d>=0:
        sgn = 1
    else:
        sgn = -1
    u = alpha_arr[n] + d - sgn*(d**2 + beta_arr[n-1]**2)**(1/2)
    return u

#Calcula o seno e o cosseno da matriz Q de transformação de Givens
def calc_theta(alpha, beta):
    raiz = (alpha**2+beta**2)**(1/2)
    cos = alpha/raiz
    sen = -beta/raiz
    return (cos, sen)

#Algoritmo QR
#Recebe uma lista da matriz diagonal principal, uma lista da matriz subdiagonal 
#e True ou False para o deslocamento espectral
def qr(alpha_list, beta_list, deslocamento):
    alpha_arr = np.array(alpha_list, dtype = float)
    beta_arr = np.array(beta_list, dtype = float)

    #Define a dimensão inicial da matriz Anxn
    dim = alpha_arr.shape[0]

    #Define V inicial (matriz identidade)
    v_arr = np.empty((dim,dim))
    for i in range(dim):
        v = np.zeros(dim)
        v[i] = 1
        v_arr[i] = v

    #Inicializa vetor para armazenar os autovalores encontrados
    autovalor_arr = np.array([])

    #Itera para m = dim, dim-1, ..., 2
    for m in range(dim, 1, -1):
        k=0
        #Itera enquanto o beta de referência da iteração m for maior do que 10e-6
        while abs(beta_arr[m-2])>=0.000001:
            #Calcula e subtrai ukI de A(k) a partir da iteração k caso haja deslocamento
            if k>0 and deslocamento:
                u = calc_heuristica(alpha_arr, beta_arr, m-1)
                u_arr = u*np.ones(alpha_arr.shape)
                alpha_arr = alpha_arr - u_arr

            #inicializa os vetores para armazenar os valores de seno e cosseno da iteração k
            c_arr = np.array([])
            s_arr = np.array([])


            #Calcula a matriz R
            beta_fixed = beta_arr.copy()
            for i in range(m-1):
                #Calcula o seno e cosseno da iteração i e guarda nos respectivos vetores
                (cos, sen) = calc_theta(alpha_arr[i],beta_fixed[i])
                c_arr = np.append(c_arr, cos)
                s_arr = np.append(s_arr, sen)

                #Cria a diagonal principal e a subdiagonal da matriz Qi*Qi-1*...*Q1*A
                alpha_R = alpha_arr.copy()
                beta_R = beta_arr.copy()

                #Modifica a linha i
                alpha_R[i] = cos*alpha_arr[i] - sen*beta_fixed[i]
                beta_R[i] = cos*beta_arr[i] - sen*alpha_arr[i+1]
                #Modifica a linha i+1
                alpha_R[i+1] = sen*beta_arr[i] + cos*alpha_arr[i+1]
                #Modifica a linha i+1 de beta apenas se não for a última iteração i=m-2
                if i < beta_arr.shape[0]-1:
                    beta_R[i+1] = cos*beta_arr[i+1]

                #Define as diagonais da matriz Qi*Qi-1*...*Q1*A como ponto de partida para a próxima iteração
                alpha_arr = alpha_R.copy()
                beta_arr = beta_R.copy()

            #Calcula A(k+1)
            for i in range(m-1):
                #Cria a diangonal principal e a subdiagonal da matriz R*Q1T*...*QiT
                alpha_next = alpha_arr.copy()
                beta_next = beta_arr.copy()

                #Modifica a coluna i
                alpha_next[i] = c_arr[i]*alpha_arr[i] - s_arr[i]*beta_arr[i]
                beta_next[i] = -s_arr[i]*alpha_arr[i+1]
                #Modifica a coluna i+1
                alpha_next[i+1] = c_arr[i]*alpha_arr[i+1]

                #Define as diagonais da matriz R*Q1T*...*QiT como ponto de partida para a próxima iteração
                alpha_arr = alpha_next.copy()
                beta_arr = beta_next.copy()

            #Soma A(k+1) e ukI a partir da iteração k caso haja deslocamento
            if k>0 and deslocamento:
                alpha_arr = alpha_arr + u_arr


            #Calcula V
            for i in range(m-1):
                v_next = v_arr.copy()
                #Altera a colunas i e i+1 para cada linha j da matriz V
                for j in range(dim):
                    v_next[j][i] = v_arr[j][i]*c_arr[i] - v_arr[j][i+1]*s_arr[i]
                    v_next[j][i+1] = v_arr[j][i]*s_arr[i] + v_arr[j][i+1]*c_arr[i]

                #Define as diagonais da matriz V*Q1T*...*QiT como ponto de partida para a próxima iteração
                v_arr = v_next.copy()
            k+=1
        
        #Guarda o autovalor encontrado na iteração m no vetor de autovalores
        autovalor_arr = np.append(autovalor_arr, alpha_arr[m-1])

        #Define as diagonais da submatriz para a iteração m-1
        alpha_arr = np.delete(alpha_arr, m-1)
        beta_arr = np.delete(beta_arr, m-2)

    #Guarda o autovetor que sobrou depois da última iteração m-1
    autovalor_arr = np.append(autovalor_arr, alpha_arr[0])

    #Inverte o vetor dos autovalores para que o autovalor da posição i corresponda ao autovetor da coluna i da matriz V
    autovalor_arr = np.array(list(reversed(autovalor_arr)))

    #Retorna o vetor com os autovalores e a matriz V, nessa ordem
    return autovalor_arr, v_arr
