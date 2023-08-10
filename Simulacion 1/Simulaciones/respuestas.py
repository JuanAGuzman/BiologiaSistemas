import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from .simulaciones import *
from .funciones import *

# Inicializacion de las constantes
const = {'kr': 1,
         'gr': 1/5,
         'kp': 30,
         'gp': 1/30,
         't': 100,
         'dt': 0.1,
         'h': 2,
         'l': 10000,
         'k': 5,
         'r0': 0,
         'p0': 0
}
const['t'] = int(10/const['gp'])
const['k'] = const['gp']/const['kp']

# Direcciones guardado de respuestas
g = os.getcwd()+'/Respuestas/Graficas/'

# Definicion de funciones
def uno(x, n, k):
    return 1

def hil(x, n, k):
    return 1/(1 + (x/k)**n)

def dis_exp(t, r):
    return np.random.exponential(r, t)

def dis_log(u, r):
    return np.log(1/u)/r

def ajuste(f, x):
    n, bins, patches = plt.hist(x)
    bin_center = (bins[:-1] + np.diff(bins)) / 2
    return curve_fit(f, bins[0:-1], n)

def s1():
    #R = open(r, 'r+')
    # b
    data= pob(sim_prim, uno, const, 500)
    p = prom(data1)
    graf(data,  const)
    print(est(p1[0][1], 'mRNA'))
    print(est(p1[0][2], 'proteinas'))

    # c
    print(tm(sim_prim, const, 100))

    # d
    data2 = pob(sim_prim_neg, const, 500)
    p2 = prom(data2)
    graf(p2, const)
    print(est(p2[0][1], 'mRNA'))
    print(est(p2[0][2], 'proteinas'))
    #R.close()


    # b
    cel = celula(const)
    pob = poblacion(cel, 500, const)
    pob.sim_prim(uno)
    pob.p()
    print(pob.est_r())
    print(pob.est_p())
    pob.graf(g, 'Simulacion estocastica primitiva')

    # c
    print(pob.tm_prim(cel, 100, uno))

    # d
    cel = celula(const)
    pob = poblacion(cel, 500, const)
    pob.sim_prim(hil)
    pob.p()
    print(pob.est_r())
    print(pob.est_p())
    pob.graf(g, 'Simulacion estocastica primitiva retroalimentacion negativa')
    

def s2():
    # a
    
    plt.hist(dis_log(np.random.uniform(0,1,100), 1), label='$\\tau = \\frac{1}{r}ln(\\frac{1}{u})$')
    plt.hist(dis_exp(100, 1), label='$r e^{-r x}$')
    plt.legend()
    plt.savefig(g+'/Ilustracion_grafica.png')
    plt.show()

    # b

    plt.hist(dis_log(np.random.uniform(0,1,500), 10), label='r=10')
    plt.savefig(g+'/500datos_r10.png')
    plt.show()

    # c
    A = dis_log(np.random.uniform(size=100), 10)
    B = dis_log(np.random.uniform(size=100), 5)
    ra, c = ajuste(dis_exp, A)
    rb, c = ajuste(dis_exp, B)
    x = np.linspace(0,4, 50)
    plt.hist(A, label='r=10')
    plt.plot(x, dis_exp(x, ra[0]), label='$\\lambda =$'+str(ra[0]-10))
    plt.legend()
    plt.savefig(g+'/A_r10.png')
    plt.show()
    x = np.linspace(0,4, 50)
    plt.hist(B, label='r=5')
    plt.plot(x, dis_exp(x, rb[0]), label='$\\lambda =$'+str(rb[0]-10))
    plt.legend()
    plt.savefig(g+'/B_r5.png')
    plt.show()
    plt.hist(A, label='r=10')
    plt.hist(B, label='r=5')
    plt.legend()
    plt.savefig(g+'/AB.png')
    plt.show()

    # d

    A = np.random.uniform(size=100)
    for i in range(len(A)):
        if np.random.uniform()<= 1/3:
            A[i] = dis_log(A[i], 5)
        else:
            A[i] = dis_log(A[i], 15)
    plt.hist(A)
    plt.savefig(g+'/CambioTipoReaccion.png')
    plt.show()

def s3():
    # b
    cel = celula(const)
    pob = poblacion(cel, 1000, const)
    pob.sim_gil(uno)
    pob.p()
    print(pob.est_r())
    print(pob.est_p())
    pob.graf(g, 'Simulacion estocastica con algoritmo de Gillespie')

    # c
    print(pob.tm_gil(cel, 100, uno))

def s4():
    # a
    cel = celula(const)
    pob = poblacion(cel, 1000, const)
    pob.sim_gil(hil)
    pob.p()
    print(pob.est_r())
    print(pob.est_p())
    pob.graf(g, 'Retroalimentacion negativa')

    # b
    

    # c
    const['h'] = 4
    cel = celula(const)
    pob = poblacion(cel, 1000, const)
    pob.sim_gil(hil)
    pob.p()
    print(pob.est_r())
    print(pob.est_p())

    # d
    const['h'] = -2
    cel = celula(const)
    pob = poblacion(cel, 1000, const)
    pob.sim_gil(hil)
    pob.p()
    print(pob.est_r())
    print(pob.est_p())
    pob.graf(g, 'Retroalimentacion positiva')

    # e
    const['p0'] = (const['kp']*const['kr'])/(const['gp']*const['gr'])
    cel = celula(const)
    pob = poblacion(cel, 1000, const)
    pob.sim_gil(hil)
    pob.p()
    pob.graf(g, 'Retroalimentacion positiva con $p_0 = $'+str(const['p0']))


def s5():
    const['p0'] = 0
    cel = celula()
    pob = poblacion(cel, 1000, const)
    pob.sim_esc_tem(uno)
    pob.p()
    pob.graf(g, 'Aproximacion por escalas temporales')

def gen_res():
    res = open(os.getcwd()+'\\Respuestas\\respuestas.txt', 'w')

    res.write('Simulacion estocastica primitiva\n')
    res.write('\n')
    # b
    res.write('b)\n')
    data= pob(sim_prim, uno, const, 500)
    p = prom(data)
    res.write(est(p[1], 'mRNA')+'\n')
    res.write(est(p[2], 'proteina')+'\n')
    graf(p, const, g, 'Simulacion estocastica primitiva')

    # c
    res.write('\n')
    res.write('c)\n')
    res.write('El delta de tiempo no tiene un efecto en el cambio de los resultados diferente del costo computacional, esto debido a que la definicion de los s tienen este factor, por lo tanto al dividirlos con s_t este factor se anula.')
    res.write(tm(sim_prim, uno, const, 100)+'\n')

    # d
    res.write('\n')
    res.write('d)\n')
    data= pob(sim_prim, hil, const, 500)
    p = prom(data)
    res.write(est(p[1], 'mRNA')+'\n')
    res.write(est(p[2], 'proteina')+'\n')
    graf(p, const, g, 'Simulacion estocastica primitiva retroalimentacion negativa')
    res.write('\n')

    res.write('Fundamentos del algoritmo de Gillespie\n')
    res.write('\n')


    res.write('\n')
    res.write('Simulacion estocastica con algoritmo de Gillespie\n')
    res.write('\n')
    # b
    res.write('b)\n')
    data= pob(sim_gil, uno, const, 500)
    p = prom(data)
    res.write(est(p[1], 'mRNA')+'\n')
    res.write(est(p[2], 'proteina')+'\n')
    graf(p, const, g, 'Simulacion estocastica con algoritmo de Gillespie')
    # c
    res.write('\n')
    res.write('c)\n')
    res.write(tm(sim_gil, uno, const, 100)+'\n')
    res.write('\n')

    res.write('Retroalimentacion\n')
    res.write('\n')
    # a
    res.write('a)\n')
    data= pob(sim_gil, hil, const, 500)
    p = prom(data)
    res.write(est(p[1], 'mRNA')+'\n')
    res.write(est(p[2], 'proteina')+'\n')
    graf(p, const, g, 'Retroalimentacion negativa')

    # b
    res.write('\n')
    

    # c
    res.write('\n')
    res.write('c)\n')
    const['h'] = 4
    data= pob(sim_gil, hil, const, 500)
    p = prom(data)
    res.write(est(p[1], 'mRNA')+'\n')
    res.write(est(p[2], 'proteina')+'\n')

    # d
    res.write('\n')
    res.write('d)\n')
    const['h'] = -2
    const['p0'] = 0.1
    data= pob(sim_gil, hil, const, 500)
    p = prom(data)
    res.write(est(p[1], 'mRNA')+'\n')
    res.write(est(p[2], 'proteina')+'\n')
    graf(p, const, g, 'Retroalimentacion positiva')

    # e
    const['p0'] = (const['kp']*const['kr'])/(const['gp']*const['gr'])
    data = pob(sim_gil, hil, const, 500)
    p = prom(data)
    graf(p, const, g, 'Retroalimentacion positiva con p_0 = '+str(const['p0']))
    res.write('\n')

    res.write('Aproximacion por escalas temporales\n')
    res.write('\n')
    # b
    res.write('b)\n')
    const['p0'] = 0
    data= pob(sim_esc_tem, uno, const, 500)
    p = prom(data)
    res.write(est(p[2], 'proteina')+'\n')
    graf(p, const, g, 'Aproximacion por escalas temporales')

    res.close()




