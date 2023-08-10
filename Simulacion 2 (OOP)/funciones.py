import numpy as np
import matplotlib.pyplot as plt
import time


def hil(x, n, k):
    """
    Funcion que recibe los parametros y calcula el valor respectivo
    de la funcion de Hill

    Parameters
    ----------
    x: float
        El numero de proteinas
    n: int
        El coeficiente de Hill que determina el tipo de retroalimentacion
    k: float
        La constante de activacion

    Returns
    -------
    float
        La fraccion de proteina unida al ligando
    """
    return 1/(1 + (x/k)**n)

 
def pob(f, c, num):
    """
    Genera un array para una poblacion (num) de celulas modeladas por
    la simulacion (f), con constantes (c)

    Parameters
    ----------
    f : class function
        Funcion que representa el tipo de simulacion
    c : dict
        Diccionario que contiene todas las constantes a usar en la simulacion
    num : int
        El tamaño de la poblacion a simular

    Returns
    -------
    array
        Un arreglo donde cada elemento es 
    """
    return np.array(list(map(f, [c for i in range(num)])))


def prom(data):
    l, ld = len(data[0][0]), len(data)
    return [np.array([np.sum([i[j] for i in data], axis=0) for j in range(3)])/ld]


def est(data, tipo):
    mean = sum(data)/len(data)
    cov = np.sqrt(sum(list(map(lambda x: (x-mean)**2, data)))/len(data))
    return 'Las '+tipo+' tienen un promedio de '+str(mean)+' y un ruido de '+str(cov/mean)+'.'


def tm(f, c, num):
    i = time.time()
    pob(f, c, num)
    return 'El tiempo de simulacion para '+str(num)+' celulas es: '+str(time.time()-i)+' segundos.'


def graf(data, c):
    fig, axs =plt.subplots(2)
    for i in data:
        axs[0].step(i[0], i[1])
        axs[1].step(i[0], i[2])
    axs[0].set_title('mRNA ($\kappa_r =$'+str(c['kr'])+' $min^{-1}$,$\gamma_r =$'+str(c['gr'])+'$min^{-1}$)')
    axs[1].set_title('Proteína ($\kappa_p =$'+str(c['kp'])+' $min^{-1}$,$\gamma_p =$'+str(c['gp'])+'$min^{-1}$)')
    axs[0].get_xaxis().set_visible(False)
    fig.supxlabel('t (min)')
    fig.supylabel('Concentración')
    fig.suptitle('Simulación estocástica primitiva')
    #plt.savefig(dir+'')
    plt.show()


