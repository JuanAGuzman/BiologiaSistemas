import numpy as np
import matplotlib.pyplot as plt
import time

def dis_exp(t, r):
    return np.random.exponential(r, t)

def dis_log(u, r):
    return np.log(1/u)/r
 
def pob(f, h, c, num):
    return np.array(list(map(f, [[c, h] for i in range(num)])))


def prom(data):
    l, ld = len(data[0][0]), len(data)
    return np.array([np.sum([i[j] for i in data], axis=0) for j in range(len(data[0]))])/ld


def est(data, tipo):
    mean = sum(data)/len(data)
    cov = np.sqrt(sum(list(map(lambda x: (x-mean)**2, data)))/len(data))
    return 'Las '+tipo+' tienen un promedio de '+str(mean)+' y un ruido de '+str(cov/mean)+'.'


def tm(f, h, c, num):
    i = time.time()
    pob(f, h, c, num)
    return 'El tiempo de simulacion para '+str(num)+' celulas es: '+str(time.time()-i)+' segundos.'


def graf(data, c, dir, title):
    fig, axs =plt.subplots(2)
    axs[0].step(data[0], data[1])
    axs[1].step(data[0], data[2])
    axs[0].set_title('mRNA ($\kappa_r =$'+str(c['kr'])+' $min^{-1}$,$\gamma_r =$'+str(c['gr'])+'$min^{-1}$)')
    axs[1].set_title('Proteína ($\kappa_p =$'+str(c['kp'])+' $min^{-1}$,$\gamma_p =$'+str(c['gp'])+'$min^{-1}$)')
    axs[0].get_xaxis().set_visible(False)
    fig.supxlabel('t (min)')
    fig.supylabel('Concentración')
    fig.suptitle(title)
    plt.savefig(dir+title+'.png')
    plt.show()


