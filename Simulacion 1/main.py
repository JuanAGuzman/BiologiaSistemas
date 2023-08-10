import matplotlib.pyplot as plt
import numpy as np
import sys,os
sys.path.append(os.getcwd())
from Simulaciones.respuestas import *


if __name__ == '__main__':
    punto = input('Que ejercicio del taller desea simular? (1/2/3/4/5/6): ')
    if int(punto)==1:
        s1()
    elif int(punto)==2:
        s2()
    elif int(punto)==3:
        s3()
    elif int(punto)==4:
        s4()
    elif int(punto)==5:
        s5()
    elif int(punto)==6:
        gen_res()
