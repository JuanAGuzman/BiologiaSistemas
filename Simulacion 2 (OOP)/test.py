import numpy as np

def sim_prim(c):
    t = np.arange(0, c['t'], c['dt'])
    r = np.zeros(len(t))
    p = np.zeros(len(t))
    r[0], p[0] = 0, 0
    for i in range(len(t)-1):
        u = np.random.uniform()
        m = np.random.uniform()
        if u<=c['kr']*c['dt']:
            r[i+1]=r[i]+1
        elif u<=c['gr']*c['dt']*r[i]:
            r[i+1]=r[i]-1
        if m<=c['kp']*c['dt']*r[i]:
            p[i+1] = p[i]+1
        elif m<=c['gp']*c['dt']*p[i]:
            p[i+1] = p[i]-1
    return np.array([t, r, p])


def sim_prim_neg(c):
    t = np.arange(0, c['t'], c['dt'])
    r = np.zeros(len(t))
    p = np.zeros(len(t))
    r[0], p[0] = 0, 0
    for i in range(len(t)-1):
        u = np.random.uniform()
        m = np.random.uniform()
        if u<=c['kr']*hil(p[i], 2, c['gr']/c['kr'])*c['dt']:
            r[i+1]=r[i]+1
        elif u<=c['gr']*c['dt']*r[i]:
            r[i+1]=r[i]-1
        if m<=c['kp']*c['dt']*r[i]:
            p[i+1] = p[i]+1
        elif m<=c['gp']*c['dt']*p[i]:
            p[i+1] = p[i]-1
    return np.array([t, r, p])

a = [5]*5
a[2] = 4
print(a)

