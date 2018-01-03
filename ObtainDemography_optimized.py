import sys, os, time
import numpy as np
from scipy import linalg
from scipy.sparse import spdiags
"""
This code is to set the demography of stratified worm burden (SWB) system in 
terms of fraction of the population of each age in each strata. It takes the following input:
    (1)Pop0: population size
    (2)Nstrata: number of worm-burden strata
    (3)Nage: number of age groups
    (4)Morrate: mortality rates for each age group
    (5)Matrate: maturation rate for each age group
For questions contact the developer at ramzi.alsallaq@gmail.com 
"""

def updated_Dpara(**kwargs):
    params={}
    params['Pop0']=23321.0
    params['Nstrata']=40
    params['Nage']=4
    #was 0-4, 5-8, 9-12, 13+ 
    #Now 0-4, 5-14, 15-29, 30+
    params['Mortrate']=np.array([0.074, 0.007, 0.025, 0.06])
    params['Matrate']= 1./np.array([5., 10., 15., 70.])
    params.update(**kwargs)
    return params


def ObtainDemographic(params=updated_Dpara()):
    """ lamb, gamma, Mortrate, and Matrate are arrays of size
        equal to the number of age groups """
    Pop0=params['Pop0']
    Nage = params['Nage']
    AgeMort = params['Mortrate']
    Maturation = params['Matrate']
    Nstrata = params['Nstrata']
    AgeFrac=np.ones(Nage)*1./Nage
    Sources = np.zeros(Nage*Nstrata)
    """-----------------Obtain Infection Free equilibrium populations-------------------------------------"""
    lamb=np.zeros(Nage, dtype=float)
    gamma=np.zeros(Nage, dtype=float)

    AgeMort_vector = np.array([elm*np.ones(Nstrata) for i,elm in enumerate(AgeMort)]).flatten() 
    hstrata0 =  np.array([np.append(elm*Pop0,np.zeros(Nstrata-1)) for i,elm in enumerate(AgeFrac)]).flatten()
    Sources.itemset(0, np.sum(AgeMort_vector*hstrata0)) 
    dhstrata = -Sources #this is b

    """ Getting Amatrix= matrix of coefficients to solve for x in Ax=b"""
    Amatrix = np.zeros((Nage*Nstrata, Nage*Nstrata))  #matrix of matrices
    L_diag = np.zeros(Nstrata)
    #CH
    i =0
    L_diag.fill(lamb[i/Nstrata])
    Main_diag = [-(lamb[i/Nstrata]+Maturation[i/Nstrata]+AgeMort[i/Nstrata]+float(j)*gamma[i/Nstrata]) for j in xrange(Nstrata)]
    Main_diag = np.asarray(Main_diag) 
    U_diag = [float(j)*gamma[i] for j in xrange(Nstrata)]
    U_diag = np.asarray(U_diag)
    data=np.vstack((L_diag, Main_diag, U_diag))
    diags = np.array([-1, 0, 1])
    Amatrix[i:i+Nstrata, i:i+Nstrata] = spdiags(data, diags, Nstrata, Nstrata).toarray()

    #SA and Y or any ages between children and adults
    for i in xrange(Nstrata, (Nage-1)*Nstrata, Nstrata):
            L_diag.fill(lamb[i/Nstrata])
            Main_diag = [-(lamb[i/Nstrata]+Maturation[i/Nstrata]+AgeMort[i/Nstrata]+float(j)*gamma[i/Nstrata]) for j in xrange(Nstrata)]
            Main_diag = np.asarray(Main_diag) # -lamb_vector[i/Nstrata] -VRate_t[i/Nstrata]   #vaccine updates
            U_diag = [float(j)*gamma[i/Nstrata] for j in xrange(Nstrata)]
            U_diag = np.asarray(U_diag)
            data=np.vstack((L_diag, Main_diag, U_diag))
            diags = np.array([-1, 0, 1])
            Amatrix[i:i+Nstrata, i:i+Nstrata] = spdiags(data, diags, Nstrata, Nstrata).toarray()
            Amatrix[i:i+Nstrata, i-Nstrata:i] = Maturation[i/Nstrata-1]*np.identity(Nstrata)

    #O
    i=(Nage-1)*Nstrata
    L_diag.fill(lamb[i/Nstrata])
    Main_diag = [-(lamb[i/Nstrata]+AgeMort[i/Nstrata]+float(j)*gamma[i/Nstrata]) for j in xrange(Nstrata)]
    Main_diag = np.asarray(Main_diag) #-lamb_vector[i/Nstrata] -VRate_t[i/Nstrata]   #vaccine updates
    U_diag = [float(j)*gamma[i/Nstrata] for j in xrange(Nstrata)]
    U_diag = np.asarray(U_diag)
    data=np.vstack((L_diag, Main_diag, U_diag))
    diags = np.array([-1, 0, 1])
    Amatrix[i:i+Nstrata, i:i+Nstrata] = spdiags(data, diags, Nstrata, Nstrata).toarray()
    Amatrix[i:i+Nstrata, i-Nstrata:i] = Maturation[i/Nstrata-1]*np.identity(Nstrata)
    """ --------------- END of Amatrix -------------------------"""

    equi_h = linalg.solve(Amatrix,dhstrata)
    nn=np.arange(0, Nstrata, dtype=float)
    nn_vector = np.array([nn for i,elm in enumerate(AgeMort)]).flatten() 
    Normal_equi_h = equi_h.copy()
    AgeFrac = np.zeros(Nage)
    #get normalized numbers
    for i in xrange(0, Nage*Nstrata, Nstrata):
        Number = equi_h[i:i+Nstrata].sum()
        Normal_equi_h[i:i+Nstrata] = equi_h[i:i+Nstrata]/Number  
        AgeFrac[i/Nstrata] =equi_h[i:i+Nstrata].sum()/equi_h.sum() 

    return AgeFrac,  equi_h
