#!/usr/bin/python3
from matplotlib.pyplot import *

x = list()
y = list()

with open("polynome.pol", encoding = "utf-8") as f:
    for ligne in f:
        donnees = ligne.split(",")
        d = donnees[1]
        x.append(float(d))
        d = donnees[2]
        y.append(float(d))

print(help(matplotlib))