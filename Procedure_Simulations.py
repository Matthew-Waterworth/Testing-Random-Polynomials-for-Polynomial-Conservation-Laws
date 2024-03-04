#!/usr/bin/env python
# coding: utf-8

##Code that generates the syzygy module basis

import sympy as sp
from sympy import QQ
import random
from operator import add

#This function takes a list of polynomials, and two lists of symbols
def Syzcalc(Polys, var, base=[]):
    
    #to get this code to work, I had to set the ks as variables,
    #but throughout the code I'm vary careful to not let them influence anything
    varlen = len(var)
    var += base
    
    #This is the part of the code that does the calculations to find the 
    #syzygy module basis
    R = QQ.old_poly_ring(*var).free_module(1)
    Syz = R.submodule(*[[f] for f in Polys]).syzygy_module()
    
    #in order to manipulate the syzygy module basis, it needs to be a non-sympy specific type
    #unfortunately, I haven't found a better way to do this
    Syzstr = '{}'.format(Syz)
    Syzstr = Syzstr.replace(">", "").replace("<", "")
    Syzstr = Syzstr.strip('][').split(', ')
    for j in range(len(Syzstr)):
        Syzstr[j] = Syzstr[j].strip('][')
    for i in range(len(Syzstr)):
        Syzstr[i] = sp.sympify(Syzstr[i])
    Syzbase = list(Syzstr)
    
    n = len(Polys)
    SyzList = [Syzstr[m:m+n] for m in range(0, len(Syzstr), n)]
    
    return SyzList, varlen, var

##This code determines random rational linear combinations of the elements of Syz

#This function takes a list of vectors as input
def ranlin(S, var, n_x):
    Gs = []
    for i in S:
        g = []
        
        #This code creates a polynomial
        nterms = 1
        termlist = []
        for j in range(nterms):
            p = (-1)**random.randint(1,2)
            C = 255*random.random()
            degs = [random.randint(0,255) for k in range(n_x)]
            factors = [var[l]**degs[l] for l in range(n_x)]
            term = 1
            for m in range(len(factors)):
                term *= factors[m]
            termlist.append(term*p*C)
        Pol = sum(termlist)
        
        #python wont let me multiply the whole vector by a scalar 
        #so I'll multiply each element
        for j in i:
            j *= Pol
            g.append(j)
        Gs.append(g)
        
    #making sure that the vectors in Gs are the same dimension
    for i in range(0,len(Gs)):
        lendiff = len(Gs[i])-len(Gs[i-1])
        if lendiff > 0:
            for j in range(lendiff):
                Gs[i-1].append(0)
        elif lendiff < 0:
            for j in range(-1*lendiff):
                Gs[i].append(0)
            
    #Adding together the 
    SumGs = Gs[0]
    for i in range(1,len(Gs)):
        SumGs = (map(add, SumGs, Gs[i]))
        
    return list(SumGs)

##This code calculates the curl of the polynomials

#This function takes a list of vectors, a list of variables, and a number as inputs
def curlzero(Gs, var, varlen):
    rlzeros = []
    for i in range(len(Gs)):
        for j in range(len(Gs[i])):
            Zero = True
            for k in range(varlen):
                calc = 0
                if j == k:
                    continue
                    
                lendiff = len(Gs[i]) - varlen
                if lendiff < 0:
                    for l in range(-1*lendiff):
                        Gs[i].append(0)
                    
                calc = sp.diff(Gs[i][k], var[j]) - sp.diff(Gs[i][j], var[k])
                if calc != 0:
                    Zero = False
                    break
                    
            #This if statement is necessary so that we dont pick up the wrong vectors
            if Zero == False:
                break
                
        if Zero == True:
            rlzeros.append(Gs[i])
                
    return rlzeros

##This code calculates the potential of any found conservative vector fields 

#This function takes as input a list of vectors, a list of variables, and a number
def inte(rlzeros, var, varlen):
    Fs = []
    c = sp.symbols('c')
    for f in rlzeros:
        F = sp.integrate(f[0], var[0])
        C = c 
        #SymPy doesnt give a constant of integration, so we have to find one ourselves
        for i in range(1, len(f)):            
            newF = F
            dF = sp.diff(F, var[i])
            rem = f[i] - dF
            I = sp.integrate(rem, var[i])
            C += I    
            
            newF += sp.expand(C)
        Fs.append(newF)
    return Fs

##This code joins together the other pieces of code to make a 
#comprehensive calculation of polynomial conservation laws

def genpolconslaw(F, n_x, n_k, kvals):
    
    for i in range(len(F)):
        for j in range(1, len(n_k)):
           F[i] = F[i].subs(n_k[j], kvals[j-1])
    
    Ss = Syzcalc(F, n_x)
    Bases, varlen, var = Ss[0], Ss[1], Ss[2]
    
    print("The {} Syzygy Module Bases: ".format(len(Bases)), Bases, "\n")
    
    GsList = [x for x in Bases]
    for i in range(500):
        Gs = ranlin(Bases, var, varlen)
        GsList.append(Gs)
        
    rlzeros = curlzero(GsList, var, varlen)
    
    Fs = inte(rlzeros, var, varlen)
    
    return Fs, len(Fs)
    