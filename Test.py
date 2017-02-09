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

Pr = P / Pc
Tr = T / Tc

Z = ChE.eosRK(Pr, Tr)
print(Z)
