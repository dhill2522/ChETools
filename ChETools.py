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
    alpha = Tr**(-1/2)
    sigma = 1
    epsilon = 0
    Ohm = 0.08664
    Psi = 0.42748
    Zc = 1/3
    Z = eosRK(Pr, Tr)
    if (sigma = epsilon):
        I = beta / (Z+epsilon*beta)
    else:
        I = 1/(sigma - epsilon)*np.log((Z+sigma*beta)/(Z+epsilon*beta))
    scaryThing = -1/2
    residH = Z - 1 +(scaryThing-1)*q*I
    residS = np.log(Z-beta)+scaryThing*q*I





# OLD - MAY NOT BE COMPLETE OR RIGOROUS
# This returns P* according to the DIPPR equation when given an array containing the
def dipprPstar(Array, T):
    Y = np.exp(Array[0]+(Array[1]/T)+Array[2]*np.log(T)+Array[3]*T**Array[4])
    return Y        # returns the vapor pressure in Pascals

def findHeatCap1 (array, T):
    C_p = array[0]+array[1]*T+array[2]*T**2+array[3]*T**3
    return C_p

def findHeatCap2 (array, T):
    C_p = array[0]+array[1]*T+array[2]*T**-2
    return C_p


# Calculating heat capacities of compounds
C5H12g_array = [1.15E-01,3.41E-04,-1.90E-07,4.23E-11]
C5H12l_array = [1.55E-01,4.37E-04,0,0]
CaOs_array = [4.18E-02,2.03E-05,-4.52E+02,0]
CH3OHg_array = [4.29E-02,8.30E-05,-1.87E-08,-8.03E-12]
CH3OHl_array = [7.59E-02,1.68E-04,0,0]
C3H8g_array = [6.80E-02,2.26E-04,-1.31E-07,3.17E-11]
COg_array = [2.90E-02,4.11E-06,3.55E-09,-2.22E-12]
CO2g_array = [3.61E-02,4.23E-05,-2.89E-08,7.46E-12]
H2Og_array = [3.35E-02,6.88E-06,7.60E-09,-3.59E-12]
N2g_array = [2.90E-02,2.20E-06,5.72E-09,-2.87E-12]
NH3g_array = [3.52E-02,2.95E-05,4.42E-09,-6.69E-12]
NOg_array = [2.95E-02,8.19E-06,-2.93E-09,3.65E-13]
O2g_array = [2.91E-02,1.16E-05,-6.08E-09,1.31E-12]
H2g_array = [2.88E-02,7.65E-08,3.29E-09,-8.70E-13]
H2Ol_array = [7.54E-02,0,0,0]
SO2g_array = [3.89E-02,3.90E-05,-3.10E-08,8.61E-12]
C2H4g_array = [4.08E-02,1.15E-04,-6.89E-08,1.77E-11]
CH4g_array = [3.43E-02,5.47E-05,3.66E-09,-1.10E-11]
C2H6g_array = [4.94E-02,1.39E-04,-5.82E-08,7.28E-12]

def calcCp(compound, T):
    # Calculates the heat capacity of the compound at the given temperature
    if (compound == "CaOs"):
        Cp = findHeatCap2(CaOs_array, T)
    else:
        Cp = findHeatCap1(compound+"_array",T)
    return Cp
