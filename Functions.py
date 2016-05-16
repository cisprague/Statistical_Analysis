# -*- coding: utf-8 -*-
"""
Created on Mon May 16 10:06:21 2016

@author: cisprague
"""
import numpy as np
from tabulate import tabulate
from scipy import stats
import matplotlib.pyplot as plt

IJ_List = np.array(
[[200, 226, 240, 261],
[278, 312, 330, 381],
[369, 416, 462, 517],
[500, 575, 645, 733]])


def ANOVA_2(IJ_List):
    '''
    IJ_List is a (I,J) shaped numpy array
    | x_i,j    x_i,j+1    ...  x_i,J   |
    | x_i+1,j  x_i+1,j+1  ...  x_i+1,J |
    | :        :          ...  :       |
    | x_I,j    x_I,j+1    ...  x_I,J   |
    '''
    IJ_List = IJ_List.astype(np.float64)
    I = np.shape(IJ_List)[0]
    J = np.shape(IJ_List)[1]
    dfi = I - 1
    dfj = J - 1
    dfe = dfi * dfj
    dft = I * J -1
    xibar = [np.mean(IJ_List[i,:]) for i in range(I)]
    xjbar = [np.mean(IJ_List[:,j]) for j in range(J)]
    xbar = np.mean(IJ_List)
    SSA = J * np.sum([(xibar[i] - xbar)**2 for i in range(I)])
    SSB = I * np.sum([(xjbar[j] - xbar)**2 for j in range(J)])
    SSE = np.sum([np.sum([
        (IJ_List[i,j] - xibar[i] - xjbar[j] + xbar)**2 for j in range(J)
        ]) for i in range(I)])
    SST = np.sum([np.sum([
        (IJ_List[i,j] - xbar)**2 for j in range(J)
        ]) for i in range(I)])
    MSA = SSA / dfi
    MSB = SSB / dfj
    MSE = SSE / dfe
    fA = MSA / MSE
    fB = MSB / MSE
    results = (
    dfi, dfj, dfe, dft,
    SSA, SSB, SSE, SST,
    MSA, MSB, MSE,
    fA, fB,
    xibar, xjbar, xbar)
    return results

def ANOVA_2_Table(IJ_List):
    dfi,dfj,dfe,dft,SSA,SSB,SSE,SST,MSA,MSB,MSE,fA,fB,xibar,xjbar,xbar = ANOVA_2(IJ_List)
    headers = [
    'Source of Variation', 'df', 'Sum of Squares', 'Mean Square', 'f'
    ]
    table = [
    ['Factor A', dfi, SSA, MSA, fA],
    ['Factor B', dfj, SSB, MSB, fB],
    ['Error',    dfe, SSE, MSE, ''],
    ['Total',    dft, SST, '',  '']
    ]
    print tabulate(table, headers = headers)

def ANOVA_2_Plot(IJ_List):
    results = ANOVA_2(IJ_List)
    xibar = results[-2]
    xjbar = results[-3]
    plt.plot(xibar)
    plt.plot(xjbar)
    plt.show()

ANOVA_2_Plot(IJ_List)
