# -*- coding: utf-8 -*-
"""
Created on Mon May 16 10:06:21 2016

@author: cisprague
"""
import numpy as np
from tabulate import tabulate
from scipy import stats
import matplotlib.pyplot as plt

IJ_List = np.array([[200, 226, 240, 261],
                    [278, 312, 330, 381],
                    [369, 416, 462, 517],
                    [500, 575, 645, 733]])

def ANOVA_2(IJ_List):
    '''
    Ouputs values of interest from factor
    ANOVA analysis with 1 observation per
    treatment.

    IJ_List is a (I,J) shaped numpy array:
    | x_i,j    x_i,j+1    ...  x_i,J   |
    | x_i+1,j  x_i+1,j+1  ...  x_i+1,J |
    | :        :          ...  :       |
    | x_I,j    x_I,j+1    ...  x_I,J   |
    '''
    IJ_List = IJ_List.astype(np.float64)
    I       = np.shape(IJ_List)[0]
    J       = np.shape(IJ_List)[1]
    dfi     = I - 1
    dfj     = J - 1
    dfe     = dfi * dfj
    dft     = I * J -1
    xibar   = [np.mean(IJ_List[i,:])
            for i in range(I)]
    xjbar   = [np.mean(IJ_List[:,j])
            for j in range(J)]
    xbar    = np.mean(IJ_List)
    SSA     = J * np.sum([(xibar[i] - xbar)**2
            for i in range(I)])
    SSB     = I * np.sum([(xjbar[j] - xbar)**2
            for j in range(J)])
    SSE     = np.sum([np.sum([
            (IJ_List[i,j] - xibar[i] - xjbar[j] + xbar)**2
            for j in range(J)])
            for i in range(I)])
    SST     = np.sum([np.sum([
            (IJ_List[i,j] - xbar)**2
            for j in range(J)])
            for i in range(I)])
    MSA     = SSA / dfi
    MSB     = SSB / dfj
    MSE     = SSE / dfe
    fA      = MSA / MSE
    fB      = MSB / MSE
    results = (
            dfi,   dfj,   dfe,   dft,
            SSA,   SSB,   SSE,   SST,
            MSA,   MSB,   MSE,
            fA,    fB,
            xibar, xjbar, xbar)
    return results

def ANOVA_2_Table(IJ_List):
    (
    dfi,   dfj,   dfe,   dft,
    SSA,   SSB,   SSE,   SST,
    MSA,   MSB,   MSE,
    fA,    fB,
    xibar, xjbar, xbar
    )       = ANOVA_2(IJ_List)
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

IJK_List = np.array([[[ 61., 63.], [ 69., 69.], [ 67., 69.]],
                     [[ 67., 69.], [ 69., 74.], [ 69., 74.]],
                     [[ 65., 74.], [ 74., 72.], [ 74., 74.]]])

def ANOVA_IJK(IJK_List):
    '''
    Outputs values of interst from a 2-factor
    ANOVA analysis with K number of observations
    per treatment.

    I       = Number of factor A levels
    J       = Number of factor B levels
    K       = Number of observations at level (i,j)
    IJK_List is a numpy array of shape (I,J,K)
    '''
    shape   = np.shape(IJK_List)
    I       = shape[0]
    J       = shape[1]
    K       = shape[2]
    dfi     = I - 1
    dfj     = J - 1
    dfij    = dfi * dfj
    dfe     = I * J * (K - 1)
    dft     = I * J * K - 1
    xibar   = np.array([np.mean(IJK_List[i,:,:])
            for i in range(I)])
    xjbar   = np.array([np.mean(IJK_List[:,j,:])
            for j in range(J)])
    xijbar  = np.asarray([[np.mean(IJK_List[i,j,:])
            for j in range(J)]
            for i in range(I)])
    xbar    = np.mean(IJK_List)
    SSA     = np.sum([np.sum([np.sum([
            (xibar[i] - xbar)**2
            for k in range(K)])
            for j in range(J)])
            for i in range(I)])
    SSB     = np.sum([np.sum([np.sum([
            (xjbar[j] - xbar)**2
            for k in range(K)])
            for j in range(J)])
            for i in range(I)])
    SSAB    = np.sum([np.sum([np.sum([
            (xijbar[i,j] - xibar[i] - xjbar[j] + xbar)**2
            for k in range(K)])
            for j in range(J)])
            for i in range(I)])
    SSE     = np.sum([np.sum([np.sum([
            (IJK_List[i,j,k] - xijbar[i,j])**2
            for k in range(K)])
            for j in range(J)])
            for i in range(I)])
    SST     = np.sum([np.sum([np.sum([
            (IJK_List[i,j,k] - xbar)**2
            for k in range(K)])
            for j in range(J)])
            for i in range(I)])
    MSA     =  SSA  / float(dfi)
    MSB     =  SSB  / float(dfj)
    MSAB    =  SSAB / float(dfij)
    MSE     =  SSE  / float(dfe)
    fA      =  MSA  / MSE
    fB      =  MSB  / MSE
    fAB     =  MSAB / MSE
    results = (dfi,   dfj,   dfij,   dfe,   dft,
               SSA,   SSB,   SSAB,   SSE,   SST,
               MSA,   MSB,   MSAB,   MSE,
               fA,    fB,    fAB,
               xibar, xjbar, xbar)
    return results

def ANOVA_IJK_Table(IJK_List):
    (
    dfi,    dfj,    dfij,    dfe,    dft,
    SSA,    SSB,    SSAB,    SSE,    SST,
    MSA,    MSB,    MSAB,    MSE,
    fA,     fB,     fAB,
    xibar,  xjbar,  xbar
    ) = ANOVA_IJK(IJK_List)
    headers = [
    'Source of Variation', 'df', 'Sum of Squares', 'Mean Square', 'f'
    ]
    table = [['Factor A',  dfi,  SSA,  MSA,  fA ],
             ['Factor B',  dfj,  SSB,  MSB,  fB ],
             ['Factor AB', dfij, SSAB, MSAB, fAB],
             ['Error',     dfe,  SSE,  MSE,  '' ],
             ['Total',     dft,  SST,  '',   '' ]]
    print tabulate(table, headers = headers)

x_list = np.array([99.0 , 101.1, 102.7, 103.0, 105.4,
                   107.0, 108.7, 110.8, 112.1, 112.4,
                   113.6, 113.8, 115.1, 115.4, 120.0])
y_list = np.array([28.8 , 27.9 , 27.0 , 25.2 , 22.8,
                   21.5 , 20.9 , 19.6 , 17.1 , 18.9,
                   16.0 , 16.7 , 13.0 , 13.6 , 10.8 ])

def Linear_Regression_Eq(x_list, y_list, verbose = False):
    '''
    x_list and y_list are cooresponding arrays
    of 1 dimensional shape and common length N.

    x_list = [x_n, x_n+1, ..., x_N]
    y_list = [y_n, y_n+1, ..., y_N]
    '''
    (slope, intercept,
    r_value, p_value,
    std_err) = stats.linregress(x_list, y_list)
    if verbose == False:
        return slope, intercept
    else:
        print('y = ' + str(intercept) + ' ' + str(slope) + ' x')

def Linear_Regression_Plot(x_list, y_list):
    plt.plot(x_list, y_list, 'o')
    x = np.linspace(min(x_list), max(x_list), 200)
    slope, intercept = Linear_Regression_Eq(x_list, y_list)
    y = intercept + slope * x
    plt.plot(x, y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

Linear_Regression_Plot(x_list, y_list)
