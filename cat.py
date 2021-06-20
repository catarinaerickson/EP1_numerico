import numpy as np;

alpha_arr = np.array([2., 2., 2., 2.])
beta_arr = np.array([-1., -1., -1])
v_arr = np.array([[1., 0., 0., 0.], [0., 1., 0., 0.],[0., 0., 1., 0.],[0., 0., 0., 1.]])
m = alpha_arr.shape[0]

def calc_heuristica(alpha_arr, beta_arr, n):
    d = (alpha_arr[n-1] - alpha_arr[n])/2
    if d>=0:
        sgn = 1
    else:
        sgn = -1
    u = alpha_arr[n] + d - sgn*(d**2 + beta_arr[n-1]**2)**(1/2)
    return u


def calc_theta(alpha, beta):
    raiz = (alpha**2+beta**2)**(1/2)
    cos = alpha/raiz
    sen = -beta/raiz
    return (cos, sen)


autovalor_arr = np.array([])

dim_fixed = alpha_arr.shape[0]
while m>=2:
    k=0
    dim = alpha_arr.shape[0]
    while abs(beta_arr[m-2])>=0.000001:
        if k>0:
            u = calc_heuristica(alpha_arr, beta_arr, m-1)
            u_arr = u*np.ones(alpha_arr.shape)
            alpha_arr = alpha_arr - u_arr

        #inicializa os vetores para armazenar os valores de sen e cos
        c_arr = np.array([])
        s_arr = np.array([])


        #cálculo da matriz R
        beta_ref = beta_arr.copy()
        for i in range(dim-1):
            (cos, sen) = calc_theta(alpha_arr[i],beta_ref[i])
            c_arr = np.append(c_arr, cos)
            s_arr = np.append(s_arr, sen)
            alpha_R = alpha_arr.copy()
            beta_R = beta_arr.copy()

            #linha i
            alpha_R[i] = cos*alpha_arr[i] - sen*beta_ref[i]
            beta_R[i] = cos*beta_arr[i] - sen*alpha_arr[i+1]
            #linha i+1
            alpha_R[i+1] = sen*beta_arr[i] + cos*alpha_arr[i+1]
            if i < beta_arr.shape[0]-1:
                beta_R[i+1] = cos*beta_arr[i+1]

            alpha_arr = alpha_R.copy()
            beta_arr = beta_R.copy()

        #calculo de A(k+1)
        for i in range(dim-1):
            alpha_next = alpha_arr.copy()
            beta_next = beta_arr.copy()

            #coluna i
            alpha_next[i] = c_arr[i]*alpha_arr[i] - s_arr[i]*beta_arr[i]
            beta_next[i] = -s_arr[i]*alpha_arr[i+1]
            #coluna i+1
            alpha_next[i+1] = c_arr[i]*alpha_arr[i+1]

            alpha_arr = alpha_next.copy()
            beta_arr = beta_next.copy()
        #soma u
        if k>0:
            alpha_arr = alpha_arr + u_arr


        # cálculo de V
        for i in range(dim-1):
            v_next = v_arr.copy()
            #altera colunas i e i+1 para cada linha j
            for j in range(4):
                v_next[j][i] = v_arr[j][i]*c_arr[i] - v_arr[j][i+1]*s_arr[i]
                v_next[j][i+1] = v_arr[j][i]*s_arr[i] + v_arr[j][i+1]*c_arr[i]

            v_arr = v_next.copy()

        k+=1
    autovalor_arr = np.append(autovalor_arr, alpha_arr[m-1])
    alpha_arr = np.delete(alpha_arr, m-1)
    beta_arr = np.delete(beta_arr, m-2)
    m=m-1

autovalor_arr = np.append(autovalor_arr, alpha_arr[m-1])

print('autovalores:', autovalor_arr)
print('autovetores:', v_arr)


autovalor_analitico = np.array([])
autovetor_analitico = np.array([])
for j in range(1,5):
    lmb = 2*(1-np.cos(j*np.pi/5))
    v = np.array([])
    for k in range (1,5):
        v = np.append(v, np.sin(k*j*np.pi/5))
    autovetor_analitico = np.append(autovetor_analitico, v)
    autovalor_analitico = np.append(autovalor_analitico, lmb)

print('autovalor analitico:', autovalor_analitico)
print('autovetor analitico',autovetor_analitico)
