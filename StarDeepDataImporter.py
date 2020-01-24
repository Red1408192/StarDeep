# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 16:42:05 2020

@author: Red
"""
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
#Read Element and Isotopes in consideration from the csv file:
f = open("Elements&Isotopes.csv", 'r')
header = f.readline().strip().split('\t')
elIso = []
for line in f:
    a = line.strip().split('\t')
    a.append(header)
    elIso.append(a)
f.close()


#read the Conversion table and create a dictionary of mettalicity percentecy-key and a dictionary
#as value, the internal dictionary return an interpolative function that return of each element
#the yield by star mass 
f = open("ConversionTable.csv", 'r')
cT, c, g, b, d, masses = {}, {}, {}, [], [], []
h, maxZetaKeys = [], []
for line in f:
    if line.strip().split('\t')[0][0] == "<":
        if len(b) != 0:
            for li in b:
                masses.append(float(li[0]))
                for i in range(len(li)-1):
                    c[d[i]].append(float(li[i+1]))
            for k in c.keys():
                g.update({k:interp1d(masses, c[k], kind='cubic', fill_value="extrapolate")})
            c, b, masses = {}, [], []
            h.append(g.copy())
            g = {}
        maxZetaKeys.append(float(line.strip().split('\t')[0].strip('<')))
        d = line.strip().split('\t')[1:]
        for u in d:
            c.update({u:[]})
        continue
    a = line.strip().split('\t')
    b.append(a)
if len(b) != 0:
    for li in b:
        masses.append(float(li[0]))
        for i in range(len(li)-1):
            c[d[i]].append(float(li[i+1]))
    for k in c.keys():
        g.update({k:interp1d(masses, c[k], kind='cubic', fill_value="extrapolate")})
    c, b, masses = {}, [], []
    h.append(g.copy())
    g = {}
f.close()
for z, i in zip(maxZetaKeys, range(len(maxZetaKeys))):
    cT.update({z:h[i]})


def getStarYield(metallicity, starMass):
    """
    Get star yield, a function that take a metallicity precentage and the star mass of star and
    return a dictionary of element+renmant keys, for each key a the value of hydrogen converted
    in solar masses. If the metallicity is superior to any key, return the last one
    """
    elements = []
    sYield = []
    for k in cT.keys():
        if metallicity < k:
            currentTable = cT[k]
            break
        currentTable = cT[k]
    for el in currentTable.keys():
        elements.append(el)
        sYield.append(float(currentTable[el](starMass)))
    return {k:v for k,v in zip(elements, sYield)}
    
def getStarAllInterpolatedYield(metallicity):
    """
    ***
    """
    x = [i for i in range(100)]
    xi = [i for i in range(200)]
    for k in cT.keys():
        if metallicity < k:
            currentTable = cT[k]
            break
        currentTable = cT[k]
    for el in currentTable.keys():
        f = currentTable[el]
        yi = f(xi)
        plt.plot(xi, yi, label = "intr.Data")
        plt.suptitle(el)
        y = f(x)
        plt.plot(x, y, "r--", label = "data")
        plt.legend(loc="best")
        plt.show()