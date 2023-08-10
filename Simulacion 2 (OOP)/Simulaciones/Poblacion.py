import time
import numpy as np
import matplotlib.pyplot as plt


class poblacion():
    def __init__(self, celula, n, c):
        self.n, self.c = n, c
        self.pob = [celula]*n
        self.tp, self.rp, self.pp = [], [], []

    def sim_prim(self, f):
        [celula.sim_prim(f) for celula in self.pob]

    def sim_gil(self, f):
        [celula.sim_gil(f) for celula in self.pob]

    def sim_esc_tem(self, f):
        [celula.sim_esc_tem(f) for celula in self.pob]

    def p(self):
        l = self.pob[0].getl()
        self.tp, self.rp, self.pp = np.zeros(l), np.zeros(l), np.zeros(l)
        for i in range(l):
            for j in range(self.n):
                cel = self.pob[j].get()
                self.tp[i] += cel[0][i]/l
                self.rp[i] += cel[1][i]/l
                self.pp[i] += cel[2][i]/l

    def est_r(self):
        mean = np.mean(self.rp)
        noise = np.std(self.rp)/mean
        return 'Los mRNA tienen un promedio de '+str(mean)+' y un ruido de '+str(noise)+'.'

    def est_p(self):
        mean = np.mean(self.pp)
        noise = np.std(self.pp)/mean
        return 'Las proteinas tienen un promedio de '+str(mean)+' y un ruido de '+str(noise)+'.'

    def tm_prim(self, celula, n, f):
        p = [celula]*n
        i = time.time()
        [celula.sim_prim(f) for celula in p]
        return 'El tiempo de simulacion para '+str(n)+' celulas es: '+str(time.time()-i)+' segundos.'

    def tm_gil(self, celula, n, f):
        p = [celula]*n
        i = time.time()
        [celula.sim_gil(f) for celula in p]
        return 'El tiempo de simulacion para '+str(n)+' celulas es: '+str(time.time()-i)+' segundos.'

    def graf(self, dir, title):
        fig, axs =plt.subplots(2)
        axs[0].step(self.tp, self.rp)
        axs[1].step(self.tp, self.pp)
        axs[0].set_title('mRNA ($\kappa_r =$'+str(self.c['kr'])+' $min^{-1}$,$\gamma_r =$'+str(self.c['gr'])+'$min^{-1}$)')
        axs[1].set_title('Proteína ($\kappa_p =$'+str(self.c['kp'])+' $min^{-1}$,$\gamma_p =$'+str(self.c['gp'])+'$min^{-1}$)')
        axs[0].get_xaxis().set_visible(False)
        fig.supxlabel('t (min)')
        fig.supylabel('Concentración')
        fig.suptitle(title)
        plt.savefig(dir+title+'.png')
        plt.show()
        