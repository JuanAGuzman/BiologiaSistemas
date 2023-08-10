import numpy as np
import matplotlib.pyplot as plt


def sim_prim(a):
    t = np.arange(0, a[0]['t'], a[0]['dt'])
    r = np.zeros(len(t))
    p = np.zeros(len(t))
    r[0], p[0] = a[0]['r0'], a[0]['p0']
    for i in range(len(t)-1):
        s = [a[0]['kr']*a[0]['dt']*a[1](p[i], a[0]['h'], a[0]['gr']/a[0]['kr']), a[0]['gr']*r[i]*a[0]['dt'], a[0]['kp']*r[i]*a[0]['dt'], a[0]['gp']*p[i]*a[0]['dt']]
        st = sum(s)
        m = np.random.uniform()
        t[i+1] = t[i]+a[0]['dt']
        if m<= s[0]/st:
            r[i+1] = r[i]+1
            p[i+1] = p[i]
        elif m<= (s[0]+s[1])/st:
            r[i+1] = r[i]-1
            p[i+1] = p[i]
        elif m<=(s[0]+s[1]+s[2])/st:
            r[i+1] = r[i]
            p[i+1] = p[i]+1
        else:
            r[i+1] = r[i]
            p[i+1] = p[i]-1
    return np.array([t, r, p])

def sim_gil(a):
    t = np.zeros(a[0]['l'])
    r = np.zeros(len(t))
    p = np.zeros(len(t))
    r[0], p[0] = a[0]['r0'], a[0]['p0']
    for i in range(a[0]['l']-1):
        s = [a[0]['kr']*a[1](p[i], a[0]['h'], a[0]['k']), a[0]['gr']*r[i], a[0]['kp']*r[i], a[0]['gp']*p[i]]
        st = sum(s)
        m =  np.random.uniform()
        t[i+1] = t[i] - np.log(np.random.uniform())/st
        if m<= s[0]/st:
            r[i+1] = r[i]+1
            p[i+1] = p[i]
        elif m<= (s[0]+s[1])/st:
            r[i+1] = r[i]-1
            p[i+1] = p[i]
        elif m<=(s[0]+s[1]+s[2])/st:
            r[i+1] = r[i]
            p[i+1] = p[i]+1
        else:
            r[i+1] = r[i]
            p[i+1] = p[i]-1
    return np.array([t, r, p])

def sim_esc_tem(a):
    t = np.zeros(a[0]['l'])
    r = np.zeros(len(t))
    p = np.zeros(len(t))
    r[0], p[0] = a[0]['r0'], a[0]['p0']
    for i in range(a[0]['l']-1):
        s = [a[0]['kr']*a[1](p[i], a[0]['h'], a[0]['k']),
            a[0]['gp']*p[i]]
        st = sum(s)
        m = np.random.uniform()
        t[i+1] = t[i]- np.log(np.random.uniform())/st
        if m<= s[0]/st:
            p[i+1] = p[i] + (a[0]['kp']/a[0]['gr'])
        else:
            p[i+1] = p[i]-1
    return np.array([t, r, p])


