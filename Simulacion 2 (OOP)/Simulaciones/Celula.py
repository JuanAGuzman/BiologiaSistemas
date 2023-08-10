import numpy as np
import matplotlib.pyplot as plt


class celula():
    def __init__(self, const):
        self.c = const
        self.t, self.r, self.p = [], [], []

    def sim_prim(self, f):
        t, r, p = [0], [self.c['r0']], [self.c['p0']]
        for i in range(self.c['l']):
            s = [self.c['kr']*self.c['dt']*f(p[i], self.c['h'], self.c['k']),
                 self.c['gr']*r[i]*self.c['dt'],
                 self.c['kp']*r[i]*self.c['dt'],
                 self.c['gp']*p[i]*self.c['dt']]
            st = sum(s)
            m = np.random.uniform()
            t.append(t[i]+self.c['dt'])
            if m<= s[0]/st:
                r.append(r[i]+1)
                p.append(p[i])
            elif m<= (s[0]+s[1])/st:
                r.append(r[i]-1)
                p.append(p[i])
            elif m<=(s[0]+s[1]+s[2])/st:
                r.append(r[i])
                p.append(p[i]+1)
            else:
                r.append(r[i])
                p.append(p[i]-1)
        self.t, self.r, self.p = t, r, p

    def sim_gil(self, f):
        t, r, p = [0], [self.c['r0']], [self.c['p0']]
        for i in range(self.c['l']):
            s = [self.c['kr']*f(p[i], self.c['h'], self.c['k']),
                 self.c['gr']*r[i],
                 self.c['kp']*r[i],
                 self.c['gp']*p[i]]
            st = sum(s)
            m = np.random.uniform()
            t.append(t[i]- np.log(np.random.uniform())/st)
            if m<= s[0]/st:
                r.append(r[i]+1)
                p.append(p[i])
            elif m<= (s[0]+s[1])/st:
                r.append(r[i]-1)
                p.append(p[i])
            elif m<=(s[0]+s[1]+s[2])/st:
                r.append(r[i])
                p.append(p[i]+1)
            else:
                r.append(r[i])
                p.append(p[i]-1)
        self.t, self.r, self.p = t, r, p

    def sim_esc_tem(self, f):
        t, p = [0], [self.c['p0']]
        for i in range(self.c['l']):
            s = [self.c['kr']*f(p[i], self.c['h'], self.c['k']),
                 self.c['gp']*p[i]]
            st = sum(s)
            m = np.random.uniform()
            t.append(t[i]- np.log(m)/st)
            if m<= s[0]/st:
                p.append(p[i] + (self.c['kp']/self.c['gr']))
            else:
                p.append(p[i]-1)
        self.t, self.p = t, p

    def get(self):
        return np.array([self.t, self.r, self.p])
    
    def getl(self):
        return len(self.t)

        

