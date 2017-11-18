#!/usr/bin/python3
from matplotlib.pyplot import *

x = list()
y = list()

with open("polynome.pol", encoding = "utf-8") as f:
    for ligne in f:
        donnees = ligne.split(",")
        d = donnees[0]
        x.append(float(d))
        d = donnees[1]
        y.append(float(d))

axvline(0, color = 'k')
axhline(0, color = 'k')
plot(x, y)
show()