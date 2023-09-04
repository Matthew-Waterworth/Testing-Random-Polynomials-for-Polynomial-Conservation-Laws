#!/usr/bin/env python
# coding: utf-8

# In[1]:


##Code that generates the syzygy module basis

import sympy as sp
from sympy import QQ

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


# In[3]:


##This code determines random rational linear combinations of the elements of Syz
import random
from operator import add


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
            
    #Adding together the elements of Gs
    SumGs = Gs[0]
    for i in range(1,len(Gs)):
        SumGs = (map(add, SumGs, Gs[i]))
        
    return list(SumGs)


# In[4]:


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


# In[5]:


##This code calculates the potential of any conservative vector fields found earlier

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


# In[6]:


##This code joins together the other pieces of code to make a 
#comprehensive calculation of polynomial conservation laws

def genpolconslaw(F, n_x, n_k, kvals):
    
    #Giving the different ks their respective values
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


# In[7]:


## Simulation 1 (genbm92)

#Setting up the appropriate amount of xs and ks based on the model

n_x = 5
n_k = 6

x = [y for y in sp.symbols('x0:{}'.format(n_x))]
k = [y for y in sp.symbols('k0:{}'.format(n_k))]

F = [k[4]*x[4] - k[2]*x[1] - k[3]*x[1]*x[2], \
    k[2]*x[1] + k[4]*x[4] + 2*k[5]*x[4] - k[3]*x[1]*x[2], \
    k[2]*x[1] + k[5]*x[4], \
    -k[4]*x[4] - k[5]*x[4] + k[3]*x[1]*x[2]]

kvals = [1, 1/250, 1000, 21/100000, 27/50000, 3/125000, 3/125000]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Conservation Laws:".format(ret[1]), ret[0])


# In[8]:


## Simulation 2 (genbm69)

n_x = 11
n_k = 22

x = [y for y in sp.symbols('x0:{}'.format(n_x))]
k = [y for y in sp.symbols('k0:{}'.format(n_k))]

F = [(k[1]*x[2]*x[5]) + (k[4]*k[7]*x[4]) + (-1*k[11]*k[2]*x[1]) + (-1*k[17]*k[2]*x[1]*x[8]),
    (k[7]*x[3]) + (k[11]*k[2]*x[1]) + (-1*k[1]*x[2]*x[5]) + (-1*k[10]*k[3]*x[2]) + (-1*k[15]*k[3]*(x[2]**2)) + (k[17]*k[2]*x[1]*x[8]) + (-1*k[14]*k[3]*x[2]*x[3]) + (-1*k[16]*k[3]*x[2]*x[4]), \
    (-1*k[7]*x[3]) + (k[10]*k[3]*x[2]) + (k[11]*k[2]*x[4]) + (k[15]*k[3]*(x[2]**2)) + (-1*k[1]*x[3]*x[5]) + (k[14]*k[3]*x[2]*x[3]) + (k[16]*k[3]*x[2]*x[4]) + (k[17]*k[2]*x[4]*x[8]), \
    (k[1]*x[3]*x[5]) + (-1*k[11]*k[2]*x[4]) + (-1*k[4]*k[7]*x[4]) + (-1*k[17]*k[2]*x[4]*x[8]), \
    (-1*k[13]*x[5]) + (k[12]*x[10]*x[6]), \
    (k[13]*x[5]) + (-1*k[12]*x[10]*x[6]), \
    (k[8]*x[8]) + (-1*k[9]*x[7]) + (-1*k[10]*k[5]*x[7]) + (-1*k[14]*k[5]*x[3]*x[7]) + (-1*k[15]*k[5]*x[2]*x[7]) + (-1*k[16]*k[5]*x[4]*x[7]), \
    (k[9]*x[7]) + (-1*k[8]*x[8]) + (k[10]*k[5]*x[7]) + (k[14]*k[5]*x[3]*x[7]) + (k[15]*k[5]*x[2]*x[7]) + (k[16]*k[5]*x[4]*x[7]), \
    (-1*k[10]*k[6]*x[9]) + (-1*k[14]*k[6]*x[3]*x[9]) + (-1*k[15]*k[6]*x[2]*x[9]) + (-1*k[16]*k[6]*x[4]*x[9]), \
    (k[13]*x[5]) + (k[10]*k[6]*x[9]) + (-1*k[12]*x[10]*x[6]) + (k[14]*k[6]*x[3]*x[9]) + (k[15]*k[6]*x[2]*x[9]) + (k[16]*k[6]*x[4]*x[9])] 

kvals = [1, 4/5, 1, 10, 1, 1, 1/20, 3/20, 7/200, 1/10000, 0, 1/10, 1/100, 1, 0, 1, 1, 0, 0, 1/10000, 1, 1, 1, 1, 1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[9]:


## Simulation 3 (genbm72)

n_x = 8
n_k = 10

x = [y for y in sp.symbols('x0:{}'.format(n_x))]
k = [y for y in sp.symbols('k0:{}'.format(n_k))]

F = [(k[3]*x[7]) + (-1*k[2]*x[1]*x[2]),
    (k[4] + (k[3]*x[7]) + (-1*k[5]*x[2]) + (-1*k[2]*x[1]*x[2])),
    (k[6]*x[4]*x[5]) + (-1*k[7]*x[3]*x[7]),
    (k[7]*x[3]*x[7]) + (-1*k[6]*x[4]*x[5]),
    (k[9]*x[6]) + (-1*k[6]*x[4]*x[5]),
    (-1*k[9]*x[6]) + (k[7]*x[3]*x[7]),
    (-1*k[3]*x[7]) + (-1*k[8]*x[7] + (k[2]*x[1]*x[2]))]

kvals = [1, 83/25000000000000000000, 1/100, 4, 1/2500, 1, 1/100000, 1/250, 11/100, 1000, 1000]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[10]:


## Simulation 4 (genbm80)

n_x = 11
n_k = 12

x = [y for y in sp.symbols('x0:{}'.format(n_x))]
k = [y for y in sp.symbols('k0:{}'.format(n_k))]

F = [(k[3]*x[2]) + (-1*k[2]*x[1]*x[10]),
    (k[10]*x[7]) + (k[5]*x[3]) + (-1*k[3]*x[2]) + (k[2]*x[1]*x[10]) + (-1*k[4]*x[2]*x[4]),
    (-1*k[5]*x[3]) + (-1*k[6]*x[3]) + (k[4]*x[2]*x[4]) + (k[7]*x[5]*x[6]),
    (k[11]*x[9]) + (k[5]*x[3]) + (-1*k[4]*x[2]*x[4]),
    (k[6]*x[3]) + (k[9]*x[7]) + (-1*k[7]*x[5]*x[6]) + (-1*k[8]*x[5]*x[8]),
    (k[6]*x[3]) + (-1*k[7]*x[5]*x[6]),
    (-1*k[10]*x[7]) + (-1*k[9]*x[7]) + (k[8]*x[5]*x[8]),
    (k[9]*x[7]) + (-1*k[8]*x[5]*x[8]),
    (k[10]*x[7]) + (-1*k[11]*x[9]),
    (k[3]*x[2]) + (-1*k[2]*x[1]*x[10])]

kvals = [1, 5000000, 10, 100000000, 1/10, 5, 100000, 5000000, 55, 1, 2, 31/1000000, 41/1000000, 1/1000000000, 100001/10000000000, 1/1000000, 11/1000000]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[11]:


## Simulation 5 (genbm82)

n_k = 12
n_x = 11

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(k[3]*x[2]) + (-1*k[2]*x[1]*x[8]),
    (k[10]*x[6]) + (k[5]*x[3]) + (-1*k[3]*x[2]) + (k[2]*x[1]*x[8]) + (-1*k[4]*x[2]*x[9]),
    (-1*k[5]*x[3]) + (-1*k[6]*x[3]) + (k[4]*x[2]*x[9]) + (k[7]*x[4]*x[5]),
    (k[6]*x[3]) + (k[9]*x[6]) + (-1*k[7]*x[4]*x[5]) + (-1*k[8]*x[4]*x[7]),
    (k[6]*x[3]) + (-1*k[7]*x[4]*x[5]),
    (-1*k[10]*x[6]) + (-1*k[9]*x[6]) + (k[8]*x[4]*x[7]),
    (k[9]*x[6]) + (-1*k[8]*x[4]*x[7]),
    (k[3]*x[2]) + (-1*k[2]*x[1]*x[8]),
    (k[11]*x[10]) + (k[5]*x[3]) + (-1*k[4]*x[2]*x[9]),
    (k[10]*x[6]) + (-1*k[11]*x[10])]

kvals = [1, 500000, 1/2, 100000000, 1/10, 1/10, 100000, 10000000, 1/10, 1/20, 1/10, 1/100000000, 1001/100000000, 1/1000000000, 10001/1000000000, 1/1000000000, 10001/1000000000]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[12]:


## Simulation 6 (genbm99)

n_k = 17
n_x = 8

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(k[11]*x[5]) + (-1*k[12]*x[1]),
    (k[9]*x[5]) + (-1*k[10]*x[2]*x[4]),
    (k[3]*x[2]) + (-1*k[4]*x[3]),
    (k[7] + (-1*k[8]*x[4]*x[7])),
    (k[1]*x[7]) + (-1*k[2]*x[5]),
    (k[13]*x[1]) + (-1*k[14]*x[3]*x[6]),
    (k[5]*x[6]) + (-1*k[6]*x[3]*x[7])]

kvals = [7/5, 9/10, 5/2, 3/2, 3/5, 4/5, 2, 13/10, 29/100, 1, 3/5, 31/10, 33, 9/2, 1, 1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[13]:


## Simulation 7 (genbm156)

n_k = 9
n_x = 4

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(k[1]*k[7]*x[1]) + (-1*k[3]*x[1]*x[2]),
    (k[6]*x[3]) + (-1*k[5]*x[2]),
    (-1*k[6]*x[3]) + (k[1]*k[4]*x[1])]

kvals = [1, 0, 37/10, 3/2, 9/10, 11/10, 2, 1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[14]:


## Simulation 8 (genbm159)

n_k = 9
n_x = 4

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(k[1]*k[2]) + (-1*k[3]*x[1]) + (-1*k[4]*x[1]*x[2]),
    (k[7]*x[3]) + (-1*k[6]*x[2]),
    (-1*k[7]*x[3]) + (k[2]*k[5]*x[1])]

kvals = [3/10, 1, 0, 16/5, 2/5, 1/10, 1/10, 1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[15]:


## Simulation 10 (genbm252)

n_k = 10
n_x = 5

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(k[1] + (k[5]*x[4])) + (k[8]*x[4]) + (-1*k[2]*x[1]) + (-1*k[6]*x[1]*x[3]),
    (k[3]*(x[1]**2)) + (-1*k[7]*x[2]),
    (k[4]*x[2]) + (k[5]*x[4]) + (k[9]*x[4]) + (-1*k[8]*x[3]) + (-1*k[6]*x[1]*x[3]),
    (-1*k[5]*x[4]) + (-1*k[8]*x[4]) + (-1*k[9]*x[4]) + (k[6]*x[1]*x[3])]

kvals = [1000, 1/10, 3/100, 7/5, 7200, 5000, 3/5, 1/5, 11, 1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[16]:


## Simulation 11 (genbm283)

n_k = 4
n_x = 5

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(k[1]*x[3]) + (-1*x[1]*x[2]),
    (k[1]*x[3]) + (k[2]*x[3]) + (-1*x[1]*x[2]),
    (x[1]*x[2]) + (-1*k[1]*x[3]) + (-1*k[2]*x[3]),
    (k[2]*x[3])]

kvals = [0, 1/2, 1, 8, 1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[17]:


## Simulation 12 (genbm357)

n_k = 14
n_x = 10

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(k[10]*x[5]) + (k[11]*x[7]) + (k[12]*x[9]) + (k[2]*x[2]) + (k[4]*x[5]) + (k[6]*x[7]) + (k[8]*x[9]) + (k[9]*x[2]) + (-1*k[1]*x[1]*x[3]) + (-1*k[3]*x[1]*x[4]) + (-1*k[5]*x[1]*x[3]) + (-1*k[7]*x[1]*x[8]),
    (-1*k[2]*x[2]) + (-1*k[9]*x[2]) + (k[1]*x[1]*x[3]),
    (k[11]*x[7]) + (k[9]*x[2]) + (-1*k[1]*x[1]*x[3]) + (-1*k[5]*x[1]*x[3]),
    (k[10]*x[5]) + (k[2]*x[2]) + (-1*k[3]*x[1]*x[4]),
    (-1*k[10]*x[5]) + (-1*k[4]*x[5]) + (k[3]*x[1]*x[4]),
    (k[4]*x[5]) + (k[8]*x[9]),
    (-1*k[11]*x[7]) + (-1*k[6]*x[7]) + (k[5]*x[1]*x[3]),
    (k[12]*x[9]) + (k[6]*x[7]) + (-1*k[7]*x[1]*x[8]),
    (-1*k[12]*x[9]) + (-1*k[8]*x[9]) + (k[7]*x[1]*x[8])]

kvals = [459/5, 412/5, 303/2, 2099/10, 129/25, 323/10, 47/10, 213/5, 167/5, 37/200, 109/5, 133/5000000, 1, 3/20000, 1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[18]:


## Simulation 13 (genbm957)

n_k = 5
n_x = 5

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(-1*k[1]*x[1]*x[2]),
    (-1*k[2]*x[2]) + (-1*k[3]*x[2]) + (k[1]*x[1]*x[2]),
    (k[3]*x[2]),
    (k[2]*x[2])]

kvals = [209/1000000000, 909/1000, 1/10, 1, 5999815]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[19]:


## Simulation 14 (genbm359)

n_k = 17
n_x = 10

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(k[10]*x[8]) + (k[3]*x[3]) + (k[5]*x[4]) + (-1*k[2]*x[1]*x[2]) + (-1*k[6]*x[1]*x[5]) + (-1*k[9]*x[1]*x[7]),
    (k[3]*x[3]) + (k[15]*x[3]*x[7]) + (-1*k[16]*x[2]*x[9]) + (-1*k[2]*x[1]*x[2]),
    (-1*k[3]*x[3]) + (-1*k[4]*x[3]) + (k[16]*x[2]*x[9]) + (k[2]*x[1]*x[2]) + (-1*k[15]*x[3]*x[7]),
    (k[12]*x[9]) + (k[4]*x[3]) + (-1*k[5]*x[4]) + (k[6]*x[1]*x[5]) + (-1*k[11]*x[4]*x[6]),
    (k[5]*x[4]) + (k[8]*x[7]) + (-1*k[6]*x[1]*x[5]) + (-1*k[7]*x[5]*x[6]),
    (k[12]*x[9]) + (k[8]*x[7]) + (-1*k[11]*x[4]*x[6]) + (-1*k[7]*x[5]*x[6]),
    (k[10]*x[8]) + (-1*k[8]*x[7]) + (k[16]*x[2]*x[9]) + (k[7]*x[5]*x[6]) + (-1*k[15]*x[3]*x[7]) + (-1*k[9]*x[1]*x[7]),
    (k[13]*x[9]) + (-1*k[10]*x[8]) + (-1*k[14]*x[8]) + (k[9]*x[1]*x[7]),
    (k[14]*x[8]) + (-1*k[12]*x[9]) + (-1*k[13]*x[9]) + (k[11]*x[4]*x[6]) + (k[15]*x[3]*x[7]) + (-1*k[16]*x[2]*x[9])]

kvals = [1,5,770,420,770,5,27/500,1/50,11/25,0,6,1/50,0,0,20,0,1,170,12/5]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[20]:


## Simulation 15 (genbm360)

n_k = 17
n_x = 10

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(k[10]*x[8]) + (k[13]*x[9]) + (k[3]*x[3]) + (k[5]*x[4]) + (-1*k[14]*x[1]*x[7]) + (-1*k[2]*x[1]*x[2]) + (-1*k[6]*x[1]*x[5]) + (-1*k[9]*x[1]*x[7]),
    (k[3]*x[3]) + (-1*k[2]*x[1]*x[2]),
    (-1*k[3]*x[3]) + (-1*k[4]*x[3]) + (k[2]*x[1]*x[2]),
    (k[12]*x[9]) + (k[4]*x[3]) + (-1*k[5]*x[4]) + (k[6]*x[1]*x[5]) + (-1*k[11]*x[4]*x[6]),
    (k[5]*x[4]) + (k[8]*x[7]) + (-1*k[6]*x[1]*x[5]) + (-1*k[7]*x[5]*x[6]),
    (k[12]*x[9]) + (k[8]*x[7]) + (-1*k[11]*x[4]*x[6]) + (-1*k[7]*x[5]*x[6]),
    (k[10]*x[8]) + (k[13]*x[9]) + (-1*k[8]*x[7]) + (k[7]*x[5]*x[6]) + (-1*k[14]*x[1]*x[7]) + (-1*k[9]*x[1]*x[7]),
    (k[15]*x[9]) + (-1*k[10]*x[8]) + (-1*k[16]*x[8]) + (k[9]*x[1]*x[7]),
    (k[16]*x[8]) + (-1*k[12]*x[9]) + (-1*k[13]*x[9]) + (-1*k[15]*x[9]) + (k[11]*x[4]*x[6]) + (k[14]*x[1]*x[7])]

kvals = [1,5,770,420,770,5,27/500,1/50,11/25,33/500,10,0,0,0,0,0,9999997/10000000, 1699999/10000, 2399999/1000000]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[21]:


## Simuluation 16 (genbm361)

n_k = 11
n_x = 9

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(k[10]*x[8]) + (k[3]*x[3]) + (k[5]*x[4]) + (-1*k[2]*x[1]*x[2]) + (-1*k[6]*x[1]*x[5]) + (-1*k[9]*x[1]*x[7]),
    (k[3]*x[3]) + (-1*k[2]*x[1]*x[2]),
    (-1*k[3]*x[3]) + (-1*k[4]*x[3]) + (k[2]*x[1]*x[2]),
    (k[4]*x[3]) + (-1*k[5]*x[4]) + (k[6]*x[1]*x[5]),
    (k[5]*x[4]) + (k[8]*x[7]) + (-1*k[6]*x[1]*x[5]) + (-1*k[7]*x[5]*x[6]),
    (k[8]*x[7]) + (-1*k[7]*x[5]*x[6]),
    (k[10]*x[8]) + (-1*k[8]*x[7]) + (k[7]*x[5]*x[6]) + (-1*k[9]*x[1]*x[7]),
    (-1*k[10]*x[8]) + (k[9]*x[1]*x[7])]
                       
kvals = [1,5,770,420,770,5,27/500,1/50,11/25,33/500,9999997/10000000,1699999/10000,2399999/1000000]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[22]:


## Simulation 17 (genbm552)

n_k = 6
n_x = 3

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(-1*k[1]*x[1]*x[2]),
    (k[2] + (-1*k[3]*x[1]) + (-1*k[4]*x[2]))]

kvals = [7/1000,33/100,21/5000,1/100,1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[23]:


## Simulation 18 (genbm553)

n_k = 6
n_x = 3

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(-1*k[1]*x[1]*x[2]),
    (k[2] + (-1*k[3]*x[1]) + (-1*k[4]*x[2]))]

kvals = [7/1000,33/100,21/5000,1/100,1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[24]:


## Simulation 19 (genbm707)

n_k = 12
n_x = 6

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(k[2] + (-1*k[1]*x[1]) + (-1*k[3]*x[1]*x[2])),
    ((k[7]*x[3]) + (-1*k[8]*x[2])),
    ((-1*k[4]*x[3]) + (k[3]*x[1]*x[2]) + (-1*k[5]*x[3]*x[4])),
    ((k[9]*x[5])) + (-1*k[10]*x[4]),
    ((-1*k[6]*x[5]) + (k[5]*x[3]*x[4]))]

kvals = [1/100,2,1/250,33/100,1/250,2,50,2,2000,2,1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[25]:


## Simulation 20 (genbm815)

n_k = 8
n_x = 3

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [((-1/32*(x[1]**2)) + (1/16*x[1]) + (-1/8*x[1]*x[2])),
    ((-1/32*(x[2]**2)) + (1/32*x[2]) + (-11/128*x[1]*x[2]))]

kvals = [0,0,0,0,0,0,0,0,0,0,0,0]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[26]:


## Simulation 21 (genbm922)

n_k = 10
n_x = 4

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [((k[1]*k[2]) + (-1*k[3]*x[1]) + (-1*k[4]*x[1])),
    ((k[3]*x[1]) + (-1*k[5]*(x[2]**2)) + (-1*k[6]*x[2]) + (-1*k[9]*x[2])),
    ((k[6]*x[2]) + (-1*k[7]*x[3]) + (-1*k[8]*x[3]))]

kvals = [20,30,361/1000,1/20,1/5000,67/500,171/500,1/400,501/10000,1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[27]:


## Simulation 22 (genbm1024)

n_k = 6
n_x = 3

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [((-1*k[2]*x[1]) + (k[1]*x[1]*x[2])),
    ((k[3]*x[2]) + (-1*k[4]*x[1]*x[2]))]

kvals = [9/5000,4/25,3/20,3/20000,1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[28]:


## Simulation 23 (genbm1031)

n_k = 6
n_x = 4

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [((k[1]*x[1]) + (-1*x[1]*x[3])),
    ((x[1]*x[3]) + (-1*k[5]*x[2])),
    ((k[5]*x[2]) + (-1*k[2]*x[3]) + (-1*k[3]*x[1]*x[3]))]

kvals = [0,3/20,1/50,1,1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[29]:


## Simulation 24 (genbm1035)

n_k = 13
n_x = 5

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [((k[1]*x[1]) + (-1*k[1]*(x[1]**2)) + (-1*k[4]*x[1]) + (-1*k[1]*x[1]*x[2]) + (-1*k[2]*x[1]*x[3]) + (-1*k[3]*x[1]*x[4])),
    ((-1*x[2]) + (k[2]*x[1]*x[3]) + (-1*k[5]*x[2]*x[4])),
    ((k[11]*x[2]) + (-1*k[7]*x[3]) + (-1*k[2]*x[1]*x[3]) + (-1*k[6]*x[3]*x[4])),
    ((-1*k[10]*x[4]) + (k[8]*x[2]*x[4]) + (k[9]*x[1]*x[4]))]

kvals = [9/25,1/2,9/2,639/5000,12/25,4/25,1/5,3/5,29/100,4/25,2,1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[30]:


## Simulation 25 (genbm1037)

n_k = 8
n_x = 3

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [((k[1]*x[1]) + (-1*k[1]*k[2]*(x[1]**2)) + (-1*k[3]*x[1]*x[2])),
    ((k[4]*x[2]) + (k[6]*x[1]*x[2]) + (-1*k[4]*k[5]*(x[2]**2)))]

kvals = [539/1250,299/100000000, 4657,5000,2213/5000,2/5,11891/10000]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[31]:


## Simulation 26 (genbm1038)

n_k = 12
n_x = 4

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [((k[1]*x[1]) + (k[9]*x[1]*x[3]) + (-1*k[1]*k[2]*(x[1]**2)) + (-1*k[3]*x[1]*x[2])),
    ((k[4]*x[2]) + (k[6]*x[1]*x[2]) + (-1*k[10]*x[2]*x[3]) + (-1*k[4]*k[5]*(x[2]**2)))]

kvals = [539/1250, 299/100000000, 9817/10000, 2213/5000, 2/5, 2291/10000, 561/625, 9611/10000, 443/2000, 199/400, 1]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[32]:


## Simulation 27 (genbm827)

n_k = 16
n_x = 11

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [((k[5]*x[5]) + (-1*k[2]*(x[1]**2)) + (2*k[1]*x[2]) + (-1*k[7]*x[1]*x[3])),
    ((1/2*k[2]*(x[1]**2)) + (-1*k[1]*x[2])),
    ((k[5]*x[5]) + (-1*k[4]*(x[3]**2)) + (2*k[3]*x[4]) + (-1*k[7]*x[1]*x[3])),
    ((k[14]*x[10]) + (1/2*k[4]*(x[3]**2)) + (-1*k[3]*x[4]) + (-1*k[11]*x[4]*x[8])),
    ((k[13]*x[9]) + (-1*k[5]*x[5]) + (k[7]*x[1]*x[3]) + (-1*k[10]*x[5]*x[8])),
    ((-1*k[8]*(x[6]**2)) + (2*k[6]*x[7])),
    ((k[12]*x[8]) + (1/2*k[8]*(x[6]**2)) + (-1*k[6]*x[7]) + (-1*k[9]*x[7])),
    ((k[9]*x[7]) + (-1*k[12]*x[8])),
    ((-1*k[13]*x[9]) + (k[10]*x[5]*x[8])),
    ((-1*k[14]*x[10]) + (k[11]*x[4]*x[8]))]

kvals = [31/25,23000000000,27/25,1900000000000,1,1,240000000000,2600000000000,9/200,10000000000,10000000000,7/250,7/250,7/250,1,9/2500000000000,11/200000000000000,37/5000000000000]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[33]:


## Simulation 28 (genbm755)

n_k = 9
n_x = 10

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(-1*k[2]*x[1]*x[2]),
    (-1*k[2]*x[1]*x[2]),
    ((k[5]*x[6]) + (k[2]*x[1]*x[2]) + (-1*k[4]*x[3]*x[4])),
    ((-1*k[4]*x[3]*x[4]) + (-1*k[3]*x[1]*x[2]*x[4])),
    ((-1*k[8]*x[5]) + (k[6]*x[3]*x[7]) + (k[3]*x[1]*x[2]*x[4])),
    ((-1*k[5]*x[6]) + (k[4]*x[3]*x[4])),
    ((k[5]*x[6]) + (-1*k[7]*x[7]) + (-1*k[6]*x[3]*x[7])),
    (k[7]*x[7]),
    (k[8]*x[5])]

kvals = [1,121267/1000,94929/200000000000000000000, 2569840000000, 387701/5, 69679400000,472749/100000000, 201671/10000000, 1/1000000000, 1/6250000,7/5000000]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[34]:


## Simulation 9 (genbm243)

n_k = 19
n_x = 24

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

kvals = [1, 19957/156250000, 1673329/2500000, 1/100000, 5946569/10000000000, 9999999/10000000,
        8875063/10000000, 4022189/5000000000, 2249759/1000000000, 602629/5000000, 2891451/100000000, 
        751457/5000000, 7204261/10000000000, 28033/78125, 1842081/500000, 278739/12500000, 32091/5000000, 1, 5083923/1000000]

F = [(-1*k[1]*x[1]*x[16]),
    (k[1]*x[1]*x[16]) + (-1*k[2]*x[2]*x[8]) + (-1*k[3]*x[17]*x[2]) + (-1*k[4]*x[18]*x[2]),
    (k[2]*x[2]*x[8]) + (-1*k[5]*x[3]*x[8]) + (-1*k[6]*x[17]*x[3]) + (-1*k[7]*x[18]*x[3]),
    (k[3]*x[17]*x[2]) + (-1*k[5]*x[4]*x[8]) + (-1*k[6]*x[17]*x[4]) + (-1*k[7]*x[18]*x[4]),
    (k[4]*x[18]*x[2]) + (-1*k[5]*x[5]*x[8]) + (-1*k[6]*x[17]*x[5]) + (-1*k[7]*x[18]*x[5]),
    (-1*k[8]*x[6]**2) + (k[10]*x[10]*x[8]) + (k[5]*x[3]*x[8]),
    (-1*k[9]*x[7]*x[9]),
    (-1*k[10]*x[10]*x[8]) + (-1*k[2]*x[2]*x[8]) + (-1*k[5]*x[3]*x[8]) + (-1*k[5]*x[4]*x[8]) + (-1*k[5]*x[5]*x[8]),
    (k[8]*x[6]**2) + (-1*k[11]*x[9]),
    (-1*k[12]*x[10]) + (k[9]*x[7]*x[9]),
    (k[5]*x[4]*x[8]) + (k[6]*x[17]*x[3]) + (-1*k[13]*x[11]*x[19]),
    (-1*k[14]*x[12]*x[14]),
    (-1*k[15]*x[13]) + (k[14]*x[12]*x[14]),
    (-1*k[16]*x[14]) + (k[13]*x[11]*x[19]),
    (k[15]*x[13]) + (-1*k[17]*x[15]),
    (-1*k[1]*x[1]*x[16]),
    (-1*k[3]*x[17]*x[2]) + (-1*k[6]*x[17]*x[3]) + (-1*k[6]*x[17]*x[4]) + (-1*k[6]*x[17]*x[5]),
    (-1*k[4]*x[18]*x[2]) + (-1*k[7]*x[18]*x[3]) + (-1*k[7]*x[18]*x[4]) + (-1*k[7]*x[18]*x[5]),
    (-1*k[13]*x[11]*x[19]),
    (k[5]*x[5]*x[8]) + (k[7]*x[18]*x[3]),
    (k[6]*x[17]*x[4]),
    (k[6]*x[17]*x[5]) + (k[7]*x[18]*x[4]),
    (k[7]*x[18]*x[5])]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[35]:


## Simulation 29 (genbm747)

n_k = 46
n_x = 34

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [((k[2]*x[2]) + (-1*k[1]*x[1]*x[3])),
    ((-1*k[2]*x[2]) + (-1*k[3]*x[2]) + (k[1]*x[1]*x[3])),
    ((k[2]*x[2]) + (k[3]*x[2]) + (k[5]*x[5]) + (k[6]*x[5]) + (-1*k[1]*x[1]*x[3]) + (-1*k[4]*x[3]*x[6])),
    ((k[14]*x[13]) + (k[3]*x[2]) + (-1*k[33]*x[4]) + (-1*k[13]*x[10]*x[4])),
    ((-1*k[5]*x[5]) + (-1*k[6]*x[5]) + (k[4]*x[3]*x[6])),
    ((k[16]*x[14]) + (k[5]*x[5]) + (-1*k[15]*x[13]*x[6]) + (-1*k[4]*x[3]*x[6])),
    ((3*x[29]) + (k[17]*x[14]) + (k[19]*x[16]) + (k[20]*x[16]) + (k[25]*x[19]) + (k[6]*x[5]) + (k[8]*x[9]) + (k[9]*x[9]) + (-1*k[34]*x[7]) + (-1*k[18]*x[15]*x[7]) + (-1*k[24]*x[17]*x[7]) + (-1*k[36]*x[30]*x[7]) + (-1*k[7]*x[7]*x[8])),
    ((k[11]*x[12]) + (k[8]*x[9]) + (-1*k[10]*x[11]*x[8]) + (-1*k[7]*x[7]*x[8])),
    ((-1*k[8]*x[9]) + (-1*k[9]*x[9]) + (k[7]*x[7]*x[8])),
    ((k[12]*x[12]) + (k[14]*x[13]) + (k[9]*x[9]) + (-1*k[32]*x[10]) + (-1*k[13]*x[10]*x[4])),
    ((3*x[32]) + (k[11]*x[12]) + (k[12]*x[12]) + (k[22]*x[18]) + (k[23]*x[18]) + (k[28]*x[21]) + (k[30]*x[23]) + (k[31]*x[23]) + (-1*k[35]*x[11]) + (-1*k[10]*x[11]*x[8]) + (-1*k[21]*x[11]*x[15]) + (-1*k[29]*x[11]*x[22]) + (-1*k[42]*x[11]*x[33])),
    ((-1*k[11]*x[12]) + (-1*k[12]*x[12]) + (k[10]*x[11]*x[8])),
    ((k[16]*x[14]) + (k[17]*x[14]) + (-1*k[14]*x[13]) + (k[13]*x[10]*x[4]) + (-1*k[15]*x[13]*x[6])),
    ((-1*k[16]*x[14]) + (-1*k[17]*x[14]) + (k[15]*x[13]*x[6])),
    ((k[19]*x[16]) + (k[22]*x[18]) + (-1*k[18]*x[15]*x[7]) + (-1*k[21]*x[11]*x[15])),
    ((-1*k[19]*x[16]) + (-1*k[20]*x[16]) + (k[18]*x[15]*x[7])),
    ((k[20]*x[16]) + (k[23]*x[18]) + (k[25]*x[19]) + (k[39]*x[31]) + (-1*k[24]*x[17]*x[7]) + (-1*k[38]*x[17]*x[29])),
    ((-1*k[22]*x[18]) + (-1*k[23]*x[18]) + (k[21]*x[11]*x[15])),
    ((3*x[31]) + (k[27]*x[21]) + (k[28]*x[21]) + (-1*k[25]*x[19]) + (k[24]*x[17]*x[7]) + (-1*k[26]*x[19]*x[20]) + (-1*k[40]*x[19]*x[30])),
    ((k[27]*x[21]) + (-1*k[26]*x[19]*x[20])),
    ((-1*k[27]*x[21]) + (-1*k[28]*x[21]) + (k[26]*x[19]*x[20])),
    ((k[30]*x[23]) + (-1*k[29]*x[11]*x[22])),
    ((-1*k[30]*x[23]) + (-1*k[31]*x[23]) + (k[29]*x[11]*x[22])),
    (k[31]*x[23]),
    (k[32]*x[10]),
    (k[33]*x[4]),
    (k[34]*x[7]),
    (k[35]*x[11]),
    ((-3*x[29]) + (k[39]*x[31]) + (k[36]*x[30]*x[7]) + (-1*k[38]*x[17]*x[29])),
    ((3*x[29]) + (3*x[31]) + (-1*k[36]*x[30]*x[7]) + (-1*k[40]*x[19]*x[30])),
    ((-3*x[31]) + (-1*k[39]*x[31]) + (k[38]*x[17]*x[29]) + (k[40]*x[19]*x[30])),
    ((-3*x[32]) + (k[42]*x[11]*x[33])),
    ((3*x[32]) + (-1*k[42]*x[11]*x[33]))]

kvals = [1/10,11/5,47/100,1/10,11/2,7/5,1/10,21/10,23/1000,1/10,15,9/10,1/10,17/100,1/10,19,29,1/10,1,43/1000, 
        1/10,36/5,13/50,1/10,1/10,1/10,100,35,1/10,720,84,11/10000,17/10000,11/1000,3/125,1/10,3,1/10,1/10,1/10,3,1/10,
        3,0,1,90,1/200,170,7/10,1400,20,7900,0,0]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[36]:


## Simulation 30 (genbm475)

n_k = 28
n_x = 24

k = [y for y in sp.symbols('k0:{}'.format(n_k))]
x = [y for y in sp.symbols('x0:{}'.format(n_x))]

F = [(-1*k[3]*x[1]*x[2]),
    (-1*k[3]*x[1]*x[2]),
    ((k[8]*x[9]) + (k[3]*x[1]*x[2]) + (-1*k[7]*x[3]*x[4])),
    ((k[12]*x[11]) + (k[4]*(x[5]**2)) + (k[8]*x[9]) + (-1*k[5]*x[4]) + (-1*k[7]*x[3]*x[4])),
    ((-2*k[4]*(x[5]**2)) + (2*k[5]*x[4])),
    ((k[11]*x[10]) + (-1*k[6]*x[6]*x[8])),
    ((k[10]*x[10]) + (k[6]*x[6]*x[8]) + (-1*k[9]*x[7]*x[9])),
    ((k[25]*x[12]) + (2*k[26]*x[19]) + (3*k[27]*x[22]) + (-1*k[13]*x[18]*x[8]) + (-1*k[6]*x[6]*x[8])),
    ((k[10]*x[10]) + (-1*k[8]*x[9]) + (k[7]*x[3]*x[4]) + (-1*k[9]*x[7]*x[9])),
    ((-1*k[10]*x[10]) + (-1*k[11]*x[10]) + (k[9]*x[7]*x[9])),
    ((k[11]*x[10]) + (-1*k[12]*x[11])),
    ((k[12]*x[11]) + (k[15]*x[14]) + (-1*k[25]*x[12]) + (-1*k[14]*x[12]*x[13])),
    ((k[15]*x[14]) + (k[19]*x[17]) + (k[24]*x[21]) + (-1*k[14]*x[12]*x[13]) + (-1*k[20]*x[13]*x[19])),
    ((k[17]*x[16]) + (-1*k[15]*x[14]) + (k[14]*x[12]*x[13]) + (-1*k[16]*x[14]*x[15])),
    ((k[17]*x[16]) + (k[22]*x[20]) + (k[13]*x[18]*x[8]) + (-1*k[16]*x[14]*x[15]) + (-1*k[21]*x[15]*x[17])),
    ((-1*k[17]*x[16]) + (-1*k[18]*x[16]) + (k[16]*x[14]*x[15])),
    ((k[18]*x[16]) + (k[22]*x[20]) + (-1*k[19]*x[17]) + (k[20]*x[13]*x[19]) + (-1*k[21]*x[15]*x[17])),
    ((k[18]*x[16]) + (k[23]*x[20]) + (-1*k[13]*x[18]*x[8])),
    ((k[19]*x[17]) + (-1*k[26]*x[19]) + (-1*k[20]*x[13]*x[19])),
    ((-1*k[22]*x[20]) + (-1*k[23]*x[20]) + (k[21]*x[15]*x[17])),
    ((k[23]*x[20]) + (-1*k[24]*x[21])),
    ((k[24]*x[21]) + (-1*k[27]*x[22])),
    ((k[25]*x[12]) + (k[26]*x[19]) + (k[27]*x[22]))]

kvals = [0,1,3/200000000,1/100,1000,1/4000000,100000,1000,351/10000,1/100,1/100,1,100000,1/200000,1/200,39/500,1/10000000000, 1/20,3/400000,1/200000,39/500,1/10000000000, 1/200,1/200,1/1250,1/200,1001,7480,206,194,8698,1520,193]

ret = genpolconslaw(F, x, k, kvals)
print("The {} Polynomial Covservation Laws:".format(ret[1]), ret[0])


# In[ ]:




