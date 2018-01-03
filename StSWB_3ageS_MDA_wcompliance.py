import sys, os, time
from scipy.misc import comb
from scipy import integrate # odeint is part of the submodule integrate
from scipy.sparse import spdiags
import numpy as np
from scipy.interpolate import interp1d
"""
This code is developed to simulate the evolution of a dynamic system of
ordinary differential equations (ODE) representing the 
transmission-contamination cycle of schistosomiasis. It is composed
of ODEs for disease acquisition and progression in humans stratified 
by age and worm-burden and ODEs for disease acquistion and progression in 
snails as intermediate hosts for schistosomiasis. Further, it simulates the 
effect of providing preventive chemotherapy for a given age group that
can be specified along with coverage among that age group, frequency and efficacy of treatment on killing (or clearing) worms within treated individuals.  
For questions please contact the developer at ramzi.alsallaq@gmail.com
"""

def updated_para(**kwargs):
    params={}
    params['Pop0']=2000.0
    params['Nage']=4
    params['zstate0']=0.0
    params['ystate0']=0.9
    params['sigma']=0.06 
    params['epsi']=52./9. #nine weeks 
    params['deldblu']=2.
    params['Mortrate']= np.array([0.127, 0.007, 0.007, 0.06])
    params['Matrate']=1./np.array([5.,4., 4., 70.])
    params.update(**kwargs)
    params['gamma']=1./5.7*np.ones(params['Nage'])
    VEff_wmort =(-np.log(1.-0.863)/(28./365.))/params['gamma'][0]-1. #from (1+VEm)*gamma=surviving fraction of worms over a given period 
    x=np.linspace(0,100,100*52) #100 years of weeks
    y=np.array([VEff_wmort]*4 + (100*52-4)*[0.0]) #only the fisrt 4 weeks the efficacy is on
    params['Treat_Eff'] = interp1d(x,y) #a step function for the efficacy versus time from when coverage starts
    return params


def dState(State, time, Sources, Amatrix_fixed, Nage, Nstrata, alpha, beta, rho0, sigma, gamma, params, phi_k, Events_time, Events, delta_t):
    # prepare matrices
    #print time
    output = np.zeros(Nage*Nstrata+2)
    hstrata = State[0:Nage*Nstrata]
    ystate = State[Nage*Nstrata:Nage*Nstrata+2]
    dhstrata = np.zeros(Nage*Nstrata)
    dystate = np.zeros(2)
    Lambda = np.zeros(Nage)
    lamb = np.zeros(Nage)
    # get some useful parameters
    deldblu = params['deldblu']
    epsi=params['epsi']
    #Treatment params
    strata_num = np.arange(Nstrata)
    lamb[0:Nage] = alpha*ystate[1]
    for j in xrange(Nage):
        #No density dependence 
        Lambda[j] = beta[j]*np.sum(phi_k*rho0[j]*hstrata[j*Nstrata:j*Nstrata+Nstrata])
    # make a copy of the fixed matrix, so that it stays intact **important not to change the fixed part


    """ ------ Update Amatrix with the time varying elements ----"""
    Amatrix = Amatrix_fixed.copy() 
    if Events is not None:
        Amatrix_T = Amatrix_fixed.copy() 
        Amatrix_N = Amatrix_fixed.copy() 
        Treat_Eff = params['Treat_Eff'](time-Events_time) # a function of time from when coverage starts
        Groups = Events['groups']
        coverage = Events['coverage']
        compliance = Events['compliance']
        TrueCoverage = coverage*compliance
        gamma_T = gamma*(1.+Treat_Eff)
        for i in xrange(0, Nage*Nstrata, Nstrata):
            if TrueCoverage[i/Nstrata] >0. and Treat_Eff>0.:
                #not treated
                np.fill_diagonal(Amatrix_N[i+1:i+Nstrata, i:i+Nstrata-1],  lamb[i/Nstrata]) #lower diagonal
                np.fill_diagonal(Amatrix_N[i:i+Nstrata, i:i+Nstrata], np.diag(Amatrix_N[i:i+Nstrata, i:i+Nstrata])-lamb[i/Nstrata]-np.asarray([float(j)*gamma[i/Nstrata] for j in xrange(Nstrata)]))
                np.fill_diagonal(Amatrix_N[i:i+Nstrata-1, i+1:i+Nstrata], np.asarray([float(j+1)*gamma[i/Nstrata] for j in xrange(Nstrata)])) #upper diagonal 
    
                #treated
                np.fill_diagonal(Amatrix_T[i+1:i+Nstrata, i:i+Nstrata-1],  lamb[i/Nstrata]) #lower diagonal
                np.fill_diagonal(Amatrix_T[i:i+Nstrata, i:i+Nstrata], np.diag(Amatrix_T[i:i+Nstrata, i:i+Nstrata])-lamb[i/Nstrata]-np.asarray([float(j)*gamma_T[i/Nstrata] for j in xrange(Nstrata)]))
                np.fill_diagonal(Amatrix_T[i:i+Nstrata-1, i+1:i+Nstrata], np.asarray([float(j+1)*gamma_T[i/Nstrata] for j in xrange(Nstrata)])) #upper diagonal 

                Amatrix[i:i+Nstrata,i:i+Nstrata] = TrueCoverage[i/Nstrata]*Amatrix_T[i:i+Nstrata,i:i+Nstrata] + (1.-TrueCoverage[i/Nstrata])*Amatrix_N[i:i+Nstrata,i:i+Nstrata]    
            else:
                #not treated
                np.fill_diagonal(Amatrix_N[i+1:i+Nstrata, i:i+Nstrata-1],  lamb[i/Nstrata]) #lower diagonal
                np.fill_diagonal(Amatrix_N[i:i+Nstrata, i:i+Nstrata], np.diag(Amatrix_N[i:i+Nstrata, i:i+Nstrata])-lamb[i/Nstrata]-np.asarray([float(j)*gamma[i/Nstrata] for j in xrange(Nstrata)]))
                np.fill_diagonal(Amatrix_N[i:i+Nstrata-1, i+1:i+Nstrata], np.asarray([float(j+1)*gamma[i/Nstrata] for j in xrange(Nstrata)])) #upper diagonal 
                Amatrix[i:i+Nstrata,i:i+Nstrata] = Amatrix_N[i:i+Nstrata,i:i+Nstrata]    
         


    else:
        for i in xrange(0, Nage*Nstrata, Nstrata):
            np.fill_diagonal(Amatrix[i+1:i+Nstrata, i:i+Nstrata-1],  lamb[i/Nstrata]) #lower diagonal
            np.fill_diagonal(Amatrix[i:i+Nstrata, i:i+Nstrata], np.diag(Amatrix[i:i+Nstrata, i:i+Nstrata])-lamb[i/Nstrata]-np.asarray([float(j)*gamma[i/Nstrata] for j in xrange(Nstrata)]))
            np.fill_diagonal(Amatrix[i:i+Nstrata-1, i+1:i+Nstrata], np.asarray([float(j+1)*gamma[i/Nstrata] for j in xrange(Nstrata)])) #upper diagonal 


    dhstrata = np.dot(Amatrix,hstrata) + Sources
    dystate[0] =Lambda.sum()*(1.-ystate.sum()) - sigma*ystate[0] 
    dystate[1] = sigma*ystate[0] - epsi*ystate[1] 
    #if dystate is restricted to positive values then decoupling between human and snails will occur; never restrain dystate by dystate[dystate<0] = 1.e-40
    output = np.append(dhstrata, [dystate[0], dystate[1]])

    return output


#dummy values
alpha=np.array([ 1.12 ,  0.21, 1.0 , 0.85])*180.
beta=alpha*1.e-10
Nstrata=120
print Nstrata
rho0=np.ones(4)*36.
sigma=0.997


def RunSim(equi_h, AgeFrac, Events, alpha=alpha,beta=beta, sigma=sigma,rho0=rho0, Nstrata=Nstrata, params=updated_para()):
    """ lamb, gamma, Mortrate, and Matrate are arrays of size
        equal to the number of age groups """
    Pop0=params['Pop0']
    deldblu=params['deldblu']
    gamma = params['gamma']
    Nage = params['Nage']
    AgeMort = params['Mortrate']
    Maturation = params['Matrate']
    ystate0=params['ystate0']
    zstate0=params['zstate0']
    Sources = np.zeros(Nage*Nstrata)
    #

    """ ----------Define fixed part of the matrix of Coefficients-------------------------------"""
    Amatrix_fixed = np.zeros((Nage*Nstrata, Nage*Nstrata))  #matrix of matrices
    L_diag = np.zeros(Nstrata)

    #CH
    i = 0
    #Main_diag = [-(Maturation[i/Nstrata]+AgeMort[i/Nstrata]+float(j)*gamma[i/Nstrata]) for j in xrange(Nstrata)]
    Main_diag = -(Maturation[i/Nstrata]+AgeMort[i/Nstrata])*np.ones(Nstrata)
    #Main_diag = np.asarray(Main_diag) 
    #U_diag = [float(j)*gamma[i/Nstrata] for j in xrange(Nstrata)]
    U_diag = np.zeros(Nstrata) 
    #U_diag = np.asarray(U_diag)
    data=np.vstack((L_diag, Main_diag, U_diag))
    diags = np.array([-1, 0, 1])
    Amatrix_fixed[i:i+Nstrata, i:i+Nstrata] = spdiags(data, diags, Nstrata, Nstrata).toarray()

    # in between age groups
    for i in xrange(Nstrata, (Nage-1)*Nstrata, Nstrata):
        #Main_diag = [-(Maturation[i/Nstrata]+AgeMort[i/Nstrata]+float(j)*gamma[i/Nstrata]) for j in xrange(Nstrata)]
        Main_diag = -(Maturation[i/Nstrata]+AgeMort[i/Nstrata])*np.ones(Nstrata)
        #Main_diag = np.asarray(Main_diag) # -lamb[i/Nstrata] -VRate_t[i/Nstrata]   #vaccine updates
        #U_diag = [float(j)*gamma[i/Nstrata] for j in xrange(Nstrata)]
        U_diag = np.zeros(Nstrata) 
        #U_diag = np.asarray(U_diag)
        data=np.vstack((L_diag, Main_diag, U_diag))
        diags = np.array([-1, 0, 1])
        Amatrix_fixed[i:i+Nstrata, i:i+Nstrata] = spdiags(data, diags, Nstrata, Nstrata).toarray()
        Amatrix_fixed[i:i+Nstrata, i-Nstrata:i] = Maturation[i/Nstrata-1]*np.identity(Nstrata)

    # old
    i=(Nage-1)*Nstrata
    #Main_diag = [-(AgeMort[i/Nstrata]+float(j)*gamma[i/Nstrata]) for j in xrange(Nstrata)]
    Main_diag = -(AgeMort[i/Nstrata])*np.ones(Nstrata)
    #Main_diag = np.asarray(Main_diag) #-lamb[i/Nstrata] -VRate_t[i/Nstrata]   #vaccine updates
    #U_diag = [float(j)*gamma[i/Nstrata] for j in xrange(Nstrata)]
    U_diag = np.zeros(Nstrata) 
    #U_diag = np.asarray(U_diag)
    data=np.vstack((L_diag, Main_diag, U_diag))
    diags = np.array([-1, 0, 1])
    Amatrix_fixed[i:i+Nstrata, i:i+Nstrata] = spdiags(data, diags, Nstrata, Nstrata).toarray()
    Amatrix_fixed[i:i+Nstrata, i-Nstrata:i] = Maturation[i/Nstrata-1]*np.identity(Nstrata)

    """-------------------END of Amatrix_fixed definition -------------------------------"""

    """ -------------- initial values ---------------------------------------------------------------------"""
    AgeMort_vector = np.array([elm*np.ones(Nstrata) for i,elm in enumerate(AgeMort)]).flatten()
    hstrata0 =  np.array([np.append(elm*Pop0,np.zeros(Nstrata-1)) for i,elm in enumerate(AgeFrac)]).flatten()
    Sources.itemset(0,  np.sum(AgeMort_vector*hstrata0))
    #calculate phi for one SWB:
    strata_num = np.arange(Nstrata)
    dblu = deldblu*strata_num
    temp = np.floor(dblu/2.)
    calc_comb = comb(dblu, temp)
    calc_comb[calc_comb==np.inf]=sys.float_info.max
    phi_k = temp*(1.-np.power(2,-dblu)*calc_comb)

    State0 = np.append(equi_h, np.array([ystate0, zstate0]))
    #stabilizing baseline
    burntime = np.linspace(0, 300, 300*10)
    delta_t = burntime[-1] - burntime[-2]
    #print State0.shape
    State_burnBL = integrate.odeint(dState, State0, burntime, args=(Sources, Amatrix_fixed, Nage, Nstrata, alpha, beta, rho0, sigma, gamma, params, phi_k, None, None, delta_t))
    State_burnBL=State_burnBL[-1,:]
    #sort events according to time
    TreatmentTimes = np.sort(Events.keys())
    #up to first time nothing happens
    time1 = np.linspace(0,TreatmentTimes[0],52)
    delta_t = time1[-1]-time1[-2]
    State1 = integrate.odeint(dState, State_burnBL, time1, args=(Sources, Amatrix_fixed, Nage, Nstrata, alpha, beta, rho0, sigma, gamma, params, phi_k, None, None, delta_t))
    State = State1.copy() 
    State1 = State1[-1,:]
    time = time1.copy()
    for tind in xrange(TreatmentTimes.shape[0]-1):
        time_dummy = np.linspace(TreatmentTimes[tind],TreatmentTimes[tind+1], 100.0)
        delta_t = time_dummy[-1] - time_dummy[-2]
        time = np.append(time, time_dummy)
        Events_dummy = Events[TreatmentTimes[tind]]
        State_dummy = integrate.odeint(dState, State1, time_dummy, args=(Sources, Amatrix_fixed, Nage, Nstrata, alpha, beta, rho0, sigma, gamma, params, phi_k, TreatmentTimes[tind], Events_dummy, delta_t))
        State = np.vstack((State, State_dummy)) 
        State1 = State_dummy[-1,:]

    timef = np.linspace(TreatmentTimes[-1], TreatmentTimes[-1]+2, 200.0)
    delta_t = timef[-1]-timef[-2]
    time = np.append(time, timef)
    Statef = integrate.odeint(dState, State1, timef, args=(Sources, Amatrix_fixed, Nage, Nstrata, alpha, beta, rho0,sigma, gamma, params, phi_k, TreatmentTimes[-1], Events[TreatmentTimes[-1]], delta_t))
    State = np.vstack((State, Statef)) 
    hstate = State[:, 0:Nage*Nstrata]
    sstate = State[:, Nage*Nstrata:Nage*Nstrata+2]
    nn=np.arange(0, Nstrata, dtype=float)
    nn_vector = np.array([nn for i,elm in enumerate(AgeMort)]).flatten()

    return time, hstate, sstate, nn_vector 

