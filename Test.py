# -*- coding: utf-8 -*-
"""
ChE 373 - Homework 11

Daniel Hill (dhill25)

Date: 9 Feb 2017
Python3
"""
import ChETools as ChE
P = 22          # bar
Pc = 42.480     # bar
T = 423         # K
Tc = 369.83     # K
R = 8.314       # J/mol K

Pr = P / Pc
Tr = T / Tc

Z = ChE.eosRK(Pr, Tr)
print(Z)

H1, S1 = ChE.residRK(Pr, Tr)

H1 = H1*R*T
S1 = S1*R
print(H, S)
