# -*- coding: utf-8 -*-
"""
Python Chemical Engineering Tools

A collection of various tools for chemical engineering applications

@author: Daniel Hill (dhill25)

Tools:
    itp1 - single interpolation
    itp2 - double interpolation
    alphaSRK - alpha for SRK EOS
    alphaPR - alpha for PR EOS
    cubicEOS - cuboc EOS, used for iterative solving
"""
from  scipy.optimize import fsolve
import numpy as np

def itp1(a, a1, valA1, a2, valA2):
    """
    Single linear interpolation tool
    Input:
        a       point to interpolate at
        a1      lower known point
        valA1   value at lower known point
        a2      upper known point
        valA2   value at upper known point
    Return:
        linearly interpolated value at a
    """
    return ((a2-a)/(a2-a1))*valA1 + ((a-a1)/(a2-a1))*valA2

def itp2(x, y, x1, y1, valX1Y1, x2, y2, valX2Y2, valX1Y2, valX2Y1):
    """
    Double linear interpolation tool
    Input:
        x       x coord of interpolation point
        y       y coord of interpolation point
        x1      x coord of lower known point
        y1      y coord of lower known point
        valX1Y1 value at x1, y1
        x2      x coord of upper known point
        y2      y coord of upper known point
        valX2Y2 value at x2, y2
        valX1Y2 value at x1, y2
        valX2Y1 value at x2, y1
    return:
        double linearly interpolated value at x, y
    """
    partA = ((((x2-x)/(x2-x1))*valX1Y1) + (((x-x1)/(x2-x1))*valX2Y1))*((y2-2)/(y2-y1))
    partB = ((((x2-x)/(x2-x1))*valX1Y2) + (((x-x1)/(x2-x1))*valX2Y2))*((y-y1)/(y2-y1))
    return (partA + partB)

def alphaSRK(Tr, w):
    return (1+(0.480+1.574*w - 0.176*w**2)*(1-Tr**(1/2)))**2

def alphaPR(Tr, w):
    return (1+(0.37464+1.54226*w-0.26992**2)*(1-Tr**(1/2)))**2

def eosCubic(z, beta, q, epsilon, sigma):
    return z-1 + beta - q*beta*((z-beta)/((z+epsilon*beta)*(z+sigma*beta)))

def iterSolve(Z, Pr, Tr, alpha, sigma, epsilon, Ohm, Psi, Zc):
    # Used for iteratively solving the EOS by other functions, only internal
    beta = Ohm * Pr / Tr
    q = Psi * alpha / (Ohm * Tr)
    return (1+beta-q*beta*((Z - beta)/((Z+epsilon*beta)*(Z+sigma*beta)))-Z)


def eosRK(Pr, Tr):
    """
    Redlich-Kwong EOS
    Input:
        Pr      Reduced pressure
        Tr      Reduced Temperature

    Output:
        Z       Compression factor
    """
    alpha = Tr**(-1/2)
    sigma = 1
    epsilon = 0
    Ohm = 0.08664
    Psi = 0.42748
    Zc = 1/3
    return fsolve(iterSolve, 0.5, (Pr, Tr, alpha, sigma, epsilon, Ohm, Psi, Zc))

def residRK(Pr, Tr):
    """
    Calculates the enthalpy and entropy residuals
    Inputs:
        Pr      Reduced pressure
        Tr      Reduced Temperature
    Output:
        Enthalpy residual divided by RT
        Entropy residual divided by R
    """
    alpha = Tr**(-1/2)
    sigma = 1
    epsilon = 0
    Ohm = 0.08664
    Psi = 0.42748
    Zc = 1/3
    beta = Ohm * Pr / Tr
    q = Psi * alpha / (Ohm * Tr)
    Z = eosRK(Pr, Tr)
    resid = np.zeros(2)

    if (sigma == epsilon):
        I = beta / (Z+epsilon*beta)
    else:
        I = 1/(sigma - epsilon)*np.log((Z+sigma*beta)/(Z+epsilon*beta))

    scaryThing = -1/2
    resid[0] = Z - 1 +(scaryThing-1)*q*I
    resid[1] = np.log(Z-beta)+scaryThing*q*I
    return (resid)
