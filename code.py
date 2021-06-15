# -*- coding: utf-8 -*-
"""
INTERPOLAÇÃO POLINOMIAL - MÉTODO DE NEWTON
"""

import numpy as np

# Ajusta o numero de casas decimais que devem ser exibidas
np.set_printoptions(precision=5)

"""
o usuario deve alterar a definicao da funcao a ser trabalhada abaixo
"""
def exact(x):
    '''
    Define-se a função exata
    '''
    # return np.exp(3*x)/(x**2-1)
    return np.sin(x+1)

def build_finite_difference_matrix(x,y,N):
    '''
    Constroi a tabela de diferenças divididas para N pontos
    '''
    DD = np.zeros((N,N))
    DD[:,0] = y
    
    for j in range(1,N):
        for i in range(N-j):
            DD[i,j] = (DD[i+1,j-1] - DD[i,j-1])/(x[i+j] - x[i])
    return DD


def poly_Newton(DD, x, xp, degree):
    '''
    Retorna:
        1. Valor interpolado pelo método de Newton
        2. Estimativa de erro
        3. Erro exato
    '''
    
    n = degree
    L = len(DD)
    if n>=L:
        print('\nErro: A tabela de pontos fornecida não é capaz de gerar o polinomio exigido.')
    elif n<=0:
        print('\nErro: Grau do polinomio deve ser maior ou igual a 1.')
    elif xp < x[0] or xp > x[-1]:
        print('\nErro: Forneça um valor para x tal que:\n x_min = %.4e < x < x_max= %.4e' % (x[0],x[-1]))
    else:
        # Encontra o indice i da tabela para a partir deste construir 
        #o polinomio interpolador
        pos = 0
        for i in range(L):
            if xp > x[i]:
                pos = i
            else:
                break
        # Ajusta o indice caso este cai fora dos limites
        pos1 = pos
        pos2 = pos1 + n
        while pos2 >= L:
            pos2 -=1
            pos1 -=1
        
        # inicia a construcao do polinomio
        Px = DD[pos1, 0]
        prod = 1
        for i in range(n):            
            prod *= (xp - x[pos1+i])
            Px += DD[pos1, i+1]*prod      
        
        # inicia a contrucao do estimador de erro
        prod = 1
        while pos2 >= L:
            pos2 -=1
            pos1 -=1
        for i in range(pos1, pos2+1):
            prod *= abs(xp - x[i])
        
        err_approx = prod*max( abs(DD[:,n+1]) )
        err_exact = abs(exact(xp)-Px)
        
        return (Px, err_approx, err_exact)

"""
o usuario deve alterar as definicoes da funcao abaixo como quiser
"""
def main():
    
    # grau maximo para o polinomio interpolador
    n = 5
    
    # Numero de pontos da tabela
    N = n + 1
    
    # vetoriza a funcao. Evita o uso do comando 'for' para calcular
    # os valores do vetor 'x'.
    vfunc = np.vectorize(exact)
    
    # pontos tabelados
    x = np.array([0, 0.62831, 1.25663, 1.8849, 2.51327, 3.14159])
    y = vfunc(x)

    # Constroi a tabela de diferencas divididas
    DD = build_finite_difference_matrix(x,y,N)
    print('Tabela de diferencas divididas:\n', DD)
    
    # grau do polinomio interpolador
    degree = 3
    
    # valor para ser interpolado
    xp = 0.5236
    
    # Calculo do valor interpolado
    Px, err_appox, err_exact = poly_Newton(DD, x, xp, degree)
    print(' P(%f) = %f\n err_appox: %f\n err_exact: %f' % (xp, Px, err_appox, err_exact))
    
    
# Inicia o programa pela chamada da funcao 'main'
main()
