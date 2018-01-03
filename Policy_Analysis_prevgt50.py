import sys, os
import time as timeM
from scipy.misc import comb
from scipy import stats
import numpy as np
from scipy.stats.mstats import mquantiles
import ObtainDemography_optimized as DemCode
reload(DemCode)
import StSWB_3ageS_MDA_wcompliance as MDACode
reload(MDACode)
import cPickle as pickle
import pandas as pd
from scipy.interpolate import interp1d

"""
This code is developed for part of the analyses in the paper:
Jaspreet Toor, Ramzi Alsallaq, James Truscott, Hugo Turner, Marleen Werkman, et al: Are we on our way to achieving the 2020 goals for schistosomiasis morbidity control using current WHO guidelines?, to appear at Clinical Infectious Diseases 2017.
For more information refer to README file
For questions please contact Ramzi Alsallaq (rxa313@case.edu or ramzi.alsallaq@gmail.com)
"""

class MDAModelRunSetUp:
    def __init__(self, **kwargs):
        self.params = {}

        self.params['Pop0']= 2000.0 
        self.params['Nage']=4
        self.params['zstate0']=0.0
        self.params['ystate0']=0.9
        self.params['epsi']= 52./9. #nine weeks 
        self.params['deldblu']=2.
        #was 0-4, 5-8, 9-12, 13+ in fitting module 
        #Now 0-4, 5-14, 15-29, 30+
        self.params['Mortrate']= np.array([0.074, 0.007, 0.025, 0.06])
        self.params['Matrate']=1./np.array([5., 10., 15., 70.])
    
        #adjusted to reflect prevalence and heavy prevalence of data points
        #the corresponding relative coefficients are np.array([0.38, 1.0, 0.02, 0.02])
        self.params['RelLamb'] = np.array([0.2, 0.53, 0.010486887, 0.010486887])
        self.params['RelRho'] = np.array([1., 1., 1., 1.0])
    
        self.params.update(**kwargs)
        self.params['gamma']=(1./5.7)*np.ones(self.params['Nage'])
        VEff_wmort =(-np.log(1.-0.863)/(28./365.))/self.params['gamma'][0]-1. #from (1+VEm)*gamma=surviving fraction of worms over a given period 
        x=np.linspace(0,100,100*52) #100 years of weeks
        y=np.array([VEff_wmort]*4 + (100*52-4)*[0.0]) #only the fisrt 4 weeks the efficacy is on
        self.params['Treat_Eff'] = interp1d(x,y) #a step function for the efficacy versus time from when coverage starts


    def Nstrata(self):
        #print self.alpha(), self.params['gamma']
        return 126

    def GetDemography(self):
        params=self.params.copy()
        params.update({'Nstrata':self.Nstrata()})
        print 'Nstrata=',params['Nstrata']
        self.AgeFrac, self.equi_h = DemCode.ObtainDemographic(params=params)

    def nwroms_npairs(self):
        nworms = np.arange(0, self.params['deldblu']*self.Nstrata())
        npairs = np.arange(0,int(nworms.max()/2.)+1)
        return nworms, npairs

    def phi(self):
        #calculate phi for one SWB:
        strata_num = np.arange(self.Nstrata())
        dblu = self.params['deldblu']*strata_num
        temp = np.floor(dblu/2.)
        calc_comb = comb(dblu, temp)
        calc_comb[calc_comb==np.inf]=sys.float_info.max
        phi_k = temp*(1.-np.power(2,-dblu)*calc_comb)
        return phi_k

def Prob_pair_nworms():
    """ This gives the number of pairs given the number of worms
        in a matrix form where number of pairs are on the rows and
        the number of worms are on the columns shape = npairsxnworms
        for a given column index (nworms), the sum of all values=1 """
    return np.load("Prob_pair_Given_nworms_deldblu2_Nstrata126.npy")


def main(argv=None):
    if argv is None:
        argv = sys.argv
    #the indices of posterior samples as they are split according to previously
    # calculated SAC prevalence at Y6
    ind_Prev_atY6ge50 = np.array([174, 175, 241, 242, 243, 244, 245, 246, 247, 248, 258, 261, 263,
       265, 266, 267, 274, 275, 276, 277, 278, 279, 280, 281, 285, 294,
       316, 317, 318, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335,
       336, 337, 338, 347, 348, 355, 356, 357, 360, 361, 366, 367, 368,
       369, 370, 371, 372, 377, 378, 379, 380, 381, 382, 383, 390, 391,
       392, 393, 394, 400, 401, 402, 404, 409, 413, 414, 424, 425, 426,
       427, 428, 429, 431, 432, 433, 434, 435])
    df_scen = pd.DataFrame({'scenN':ind_Prev_atY6ge50})
    df_full = pd.DataFrame({'scenN':range(436)})

    ScenN=int(argv[1])

    Events = {}
   
    if ScenN==1: #annual PCT for SAC @75% for scenarios with 10%<SACprev<50% at year 6  
        id_samples = df_full[~df_full['scenN'].isin(df_scen['scenN'])]['scenN'].values #apply annual to them from year 6 to year 9 because their Prev<50%
        #id_samples = df_full['scenN'].values #apply annual to them from year 6 to year 9 because their Prev<50%
        alpha_factor = 1.0
        beta_factor = 1.0
        TotCovSAC = np.array(30*[0.75])
        TotCovAdults = np.array(30*[0.0])
        for i in xrange(1,22):
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[i-1]]+2*[TotCovAdults[i-1]]), 'compliance':np.array(4*[1.0])}
        timepoints = np.arange(0,21,0.05)
    elif ScenN==2: #Twice A year PCT for SAC @75% for scenarios with SACprev>=50% at year 6 
        id_samples = ind_Prev_atY6ge50 #apply twice a year from year 6 to year 9 because their Prev>=50%
        alpha_factor = 1.0
        beta_factor = 1.0
        TotCovSAC = np.array(30*[0.75])
        TotCovAdults = np.array(30*[0.0])
        for i in np.arange(1,7): #from Y0-Y5
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[i-1]]+2*[TotCovAdults[i-1]]), 'compliance':np.array(4*[1.0])}
        for i in np.arange(7,22,0.5): #from Y6-Y20
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[int(i-1)]]+2*[TotCovAdults[int(i-1)]]), 'compliance':np.array(4*[1.0])}
        timepoints = np.arange(0,21,0.05)
    elif ScenN==3: ##Annual Y0-Y9, then Y10- +40% adults for scenarios with 10%<SACPrev<50% at Y6 
        id_samples = df_full[~df_full['scenN'].isin(df_scen['scenN'])]['scenN'].values #apply annual to them from year 6 to year 9 because their Prev<50%
        alpha_factor = 1.0
        beta_factor = 1.0
        TotCovSAC = np.array(10*[0.75]+19*[0.75])
        TotCovAdults = np.array(10*[0.0]+19*[0.40])
        for i in xrange(1,22):
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[i-1]]+2*[TotCovAdults[i-1]]), 'compliance':np.array(4*[1.0])}
        timepoints = np.arange(0,21,0.05)
    elif ScenN==4: ## twice a year form year6 , then Y10- +40% adults for scenarios with SACprev>=50% at Y6
        id_samples = ind_Prev_atY6ge50 #apply twice a year from year 6 to year 9 because their Prev>=50%
        alpha_factor = 1.0
        beta_factor = 1.0
        TotCovSAC = np.array(10*[0.75]+19*[0.75])
        TotCovAdults = np.array(10*[0.0]+19*[0.40])
        for i in np.arange(1,7): #from Y0-Y5 Annual
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[i-1]]+2*[TotCovAdults[i-1]]), 'compliance':np.array(4*[1.0])}
        for i in np.arange(7,22,0.5): #from Y6-Y10 twice a year then from y10 twice a year+40% for adults 
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[int(i-1)]]+2*[TotCovAdults[int(i-1)]]), 'compliance':np.array(4*[1.0])}
        timepoints = np.arange(0,21,0.05)
    elif ScenN==5: #twice a year from y10 
        id_samples = df_full[~df_full['scenN'].isin(df_scen['scenN'])]['scenN'].values #apply annual to them from year 6 to year 9 because their Prev<50%
        alpha_factor = 1.0
        beta_factor = 1.0
        TotCovSAC = np.array(20*[0.75]+19*[0.75])
        TotCovAdults = np.array(20*[0.0]+19*[0.0])
        for i in xrange(1,11): #Y0-9 annual
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[i-1]]+2*[TotCovAdults[i-1]]), 'compliance':np.array(4*[1.0])}
        for i in np.arange(11,22,0.5): #from Y10-Y20 twice a year
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[int(i-1)]]+2*[TotCovAdults[int(i-1)]]), 'compliance':np.array(4*[1.0])}
        timepoints = np.arange(0,21,0.05)

    elif ScenN==6: ##Annual Y0-Y9 @75%, then Y10- annual at 85% + 40% adults 
        id_samples = df_full[~df_full['scenN'].isin(df_scen['scenN'])]['scenN'].values #apply annual to them from year 6 to year 9 because their Prev<50%
        alpha_factor = 1.0
        beta_factor = 1.0
        TotCovSAC = np.array(10*[0.75]+19*[0.85])
        TotCovAdults = np.array(10*[0.0]+19*[0.40])
        for i in xrange(1,22):
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[i-1]]+2*[TotCovAdults[i-1]]), 'compliance':np.array(4*[1.0])}
        timepoints = np.arange(0,21,0.05)
    elif ScenN==7: ## twice a year from Y6 @75%, then Y10- twice a year at 85% +40% adults
        id_samples = ind_Prev_atY6ge50 #apply twice a year from year 6 to year 9 because their Prev>=50%
        alpha_factor = 1.0
        beta_factor = 1.0
        TotCovSAC = np.array(10*[0.75]+19*[0.85])
        TotCovAdults = np.array(10*[0.0]+19*[0.40])
        for i in np.arange(1,7): #from Y0-Y5 Annual
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[i-1]]+2*[TotCovAdults[i-1]]), 'compliance':np.array(4*[1.0])}
        for i in np.arange(7,22,0.5): #from Y6-Y10 twice a year then from y10 twice a year+40% for adults 
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[int(i-1)]]+2*[TotCovAdults[int(i-1)]]), 'compliance':np.array(4*[1.0])}
        timepoints = np.arange(0,21,0.05)

    else:
        id_samples = df_full['scenN'].values 
        alpha_factor = 1.0
        beta_factor = 1.0
        TotCovSAC = np.array(24*[0.])
        TotCovAdults = np.array(6*[0.0,0.0,0.0,0.0])
        for i in xrange(1,10):
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array([0.]+[TotCovSAC[i-1]]+[TotCovAdults[i-1]]+[0.]), 'compliance':np.array(4*[1.0])}
        timepoints = np.arange(0,9,0.05)

    setup = MDAModelRunSetUp()
    setup.GetDemography()
    AgeFrac = setup.AgeFrac
    equi_h = setup.equi_h
    phi_k = setup.phi()
    nworms, npairs = setup.nwroms_npairs()
    params=setup.params
    params.update({'Events':Events})
    Pop0 = params['Pop0']
    deldblu = params['deldblu']
    Nage = params['Nage']
    RelLamb = params['RelLamb']
    Mortrate = params['Mortrate']
    gamma = params['gamma']    
    epsi = params['epsi']
    Nstrata=setup.Nstrata()
    run_params = params.copy()
    run_params['equi_h'] = equi_h
    run_params['AgeFrac'] = AgeFrac
    run_params['Nstrata'] = Nstrata
    run_params['Pop0'] = Pop0
    run_params['Nage'] = Nage
    
    """ specify params """
    #Parameter specification against SCORE experience in high prevalence villages
    pfile=pickle.load(open("MCMC.p","rb"))
    rho0 = pfile["rho0"][0][600:]
    #adjusted to reflect prevalence and heavy prevalence of data points
    clpara = pfile["clpara"][0][600:]*7.0 
    #I do not read alphas and betas directly because these were generated for a different age grouping
    #data points are for age groups 5-8, 9-12 and 13+, but here we have 0-4, 5-14, 15-29, and 30+
    logbeta0 = pfile["logbeta0"][0][600:]
    logalphaSAC = pfile["logalphaSAC"][0][600:]
    alphas = np.power(10,logalphaSAC)[:,None]*RelLamb[None,:] 
    betas = alphas*np.power(10,logbeta0)[:,None]
    sigma = pfile["sigma"][0][600:]
    zmu = pfile["zmu"][0][600:]
    zstate0=np.exp(zmu)/(1.+np.exp(zmu))
    LProb_eggs_age = pfile["Prob_eggs_age"][0][600:]
    Prob_eggs_age = np.exp(LProb_eggs_age)/(1.+np.exp(LProb_eggs_age))
    
    Eq_strata = pfile["Eq_strata"][0][600:]

    #Calculating Rnot
    #Averaging alpha and beta for one age group
    alpha = (alphas*AgeFrac).sum(axis=1)
    beta = (betas*AgeFrac).sum(axis=1)
    mu = (Mortrate*AgeFrac).sum()
    gamma = gamma[0]
    #Rnot has the same indicies as alphas and betas
    Rnot = np.sqrt(Pop0*alpha*beta*rho0*phi_k[1]/(epsi*(gamma+mu)))
    print "full number of samples=", Rnot.shape[0]
    SACPrev0 = 1.-Prob_eggs_age[:,0,2]#here index 2 was referring to 9-12 
    #samples up to Rnot<6.0
    Rpairs = zip(np.arange(1,5,0.5),np.arange(1.5,6,0.5))
    i_samples=np.array([], dtype=int)
    for Rmn, Rmx in Rpairs:
        temp = np.where((Rnot>=Rmn)&(Rnot<Rmx))[0][-80:] 
        i_samples = np.append(i_samples, temp)
    nsamples = i_samples.shape[0]
    #selected samples
    print "number of selected samples=",nsamples

    #to calculate prevalence
    intdeldblu = int(deldblu)
    outshape1 = Nstrata*intdeldblu
    #probability of having certain number of pairs up to Nworms/2 given number of worms Nworms
    Prob_pair_nworms = np.load("Prob_pair_Given_nworms_deldblu2_Nstrata126.npy") 
    pop = lambda x:(x[:,0:Nstrata].sum(axis=1),x[:,Nstrata:2*Nstrata].sum(axis=1),x[:,2*Nstrata:3*Nstrata].sum(axis=1),x[:,3*Nstrata:4*Nstrata].sum(axis=1))

    nbinomcdf = stats.distributions.nbinom.cdf
    #id_samples are the ids of the samples defined according to scenario
    # and connected to Events, they are split according to previously
    # calculated SAC prevalence at Y6
    SAC_eqwdist = 0.0
    for s in i_samples[id_samples][0:2]:
        alpha_i = alphas[s]*alpha_factor
        beta_i = betas[s]*beta_factor
        sigma_i = sigma[s]
        rho0_i=rho0[s]
        clpara_i=clpara[s]

        run_params['sigma'] = sigma_i
        run_params['rho0'] = rho0_i*np.ones(Nage)
        run_params['clpara'] = clpara_i
        run_params['ystate0'] = 0.9 
        run_params['zstate0'] = zstate0[s]
        run_params.update({'sigma':sigma})

        time, hstate, sstate, nn_vector = MDACode.RunSim(equi_h=equi_h, AgeFrac=AgeFrac, Events=Events, alpha=alpha_i,beta=beta_i, sigma=sigma_i, rho0=rho0_i*np.ones(Nage), Nstrata=Nstrata, params=run_params)

        #To calculate prev
        NB_n = 1./clpara_i
        NB_p =  NB_n/(rho0_i*npairs + NB_n)
        nb0 = np.power( NB_n/(rho0_i*npairs + NB_n) , NB_n ) #for zero prevalence
        nbgt400 = 1.0-nbinomcdf(400.,n=NB_n, p=NB_p) #for heavy prevalence
        ind_mat = np.array([np.where(time>=timepoints[j])[0][0] for j in xrange(len(timepoints))])
        times = time[ind_mat]
        hstate = hstate[ind_mat]
        sstate = sstate[ind_mat]
        Incidence =sstate[:,1]*deldblu*np.sum([hstate[:,jj*Nstrata:jj*Nstrata+Nstrata-1].sum(axis=1)*alpha_i[jj] for jj in xrange(Nage)], axis=0 )/hstate.sum(axis=1)

        Pop_tupe = pop(hstate) 
        #Normalized states each with size timepointsxNstrata 
        Inf= hstate[:,0:Nstrata]/Pop_tupe[0][:,None]
        SAC= hstate[:,Nstrata:2*Nstrata]/Pop_tupe[1][:,None]
        YAdults= hstate[:,2*Nstrata:3*Nstrata]/Pop_tupe[2][:,None]
        Adults= hstate[:,3*Nstrata:4*Nstrata]/Pop_tupe[3][:,None]
        #as a measure of spread in worm distribution
        Mean_worm_SAC = np.sum(hstate[-1,Nstrata:2*Nstrata]*np.arange(Nstrata)*deldblu*0.5)/Pop_tupe[1][-1] #mean number of worms
        var_worm_SAC = np.sum( hstate[-1,Nstrata:2*Nstrata]*np.power(np.arange(Nstrata)*deldblu*0.5,2))/Pop_tupe[1][-1] - np.power(Mean_worm_SAC,2) #variance of the number of worms in a strata from the mean

        Aggregation_SWB_SAC = np.sqrt(var_worm_SAC)/Mean_worm_SAC 
        SAC_eqwdist = hstate[-1,Nstrata:2*Nstrata] 

        #Prob_age_nworms for different age groups    
        Prob_age_nworms_Inf,Prob_age_nworms_SAC,Prob_age_nworms_YAdults,Prob_age_nworms_Adults = np.zeros((ind_mat.shape[0], outshape1)), np.zeros((ind_mat.shape[0], outshape1)),np.zeros((ind_mat.shape[0], outshape1)), np.zeros((ind_mat.shape[0], outshape1))

        for i in xrange(Nstrata):
            Prob_age_nworms_Inf[:,i*intdeldblu:(i+1)*intdeldblu] = Inf[:,i][:,None]
            Prob_age_nworms_SAC[:,i*intdeldblu:(i+1)*intdeldblu] = SAC[:,i][:,None]
            Prob_age_nworms_YAdults[:,i*intdeldblu:(i+1)*intdeldblu] = YAdults[:,i][:,None]
            Prob_age_nworms_Adults[:,i*intdeldblu:(i+1)*intdeldblu] = Adults[:,i][:,None]
        Prob_age_nworms_Inf = Prob_age_nworms_Inf.reshape((Inf.shape[0],1,Nstrata*intdeldblu))/deldblu #division by deldblu to normalize
        Prob_age_nworms_SAC = Prob_age_nworms_SAC.reshape((SAC.shape[0],1,Nstrata*intdeldblu))/deldblu #division by deldblu to normalize
        Prob_age_nworms_YAdults = Prob_age_nworms_YAdults.reshape((YAdults.shape[0],1,Nstrata*intdeldblu))/deldblu #division by deldblu to normalize
        Prob_age_nworms_Adults = Prob_age_nworms_Adults.reshape((Adults.shape[0],1,Nstrata*intdeldblu))/deldblu #division by deldblu to normal

        #From Prob_eggs_age: the full deVlas modified model:
        Prev_Inf = np.zeros((ind_mat.shape[0],1))
        Prev_SAC = np.zeros((ind_mat.shape[0],1))
        Prev_YAdults = np.zeros((ind_mat.shape[0],1))
        Prev_Adults = np.zeros((ind_mat.shape[0],1))
        #nb0=Prob_eggs_wormpair_I[0,:] #size=npairs
        HPrev_Inf = np.zeros((ind_mat.shape[0],1))
        HPrev_SAC = np.zeros((ind_mat.shape[0],1))
        HPrev_YAdults = np.zeros((ind_mat.shape[0],1))
        HPrev_Adults = np.zeros((ind_mat.shape[0],1))
        MeanInt_Inf = np.zeros((ind_mat.shape[0],1))
        MeanInt_SAC = np.zeros((ind_mat.shape[0],1))
        MeanInt_YAdults = np.zeros((ind_mat.shape[0],1))
        MeanInt_Adults = np.zeros((ind_mat.shape[0],1))
        for t in xrange(ind_mat.shape[0]):
            #Prevalence and heavy prevalence both are obtained using deVlas model
            Prev_Inf[t,:]=1.-np.dot( nb0, np.dot( Prob_pair_nworms,Prob_age_nworms_Inf[t,:,:].T) )
            Prev_SAC[t,:]=1.-np.dot( nb0, np.dot( Prob_pair_nworms,Prob_age_nworms_SAC[t,:,:].T) )
            Prev_YAdults[t,:]=1.-np.dot( nb0, np.dot( Prob_pair_nworms,Prob_age_nworms_YAdults[t,:,:].T) )
            Prev_Adults[t,:]=1.-np.dot( nb0, np.dot( Prob_pair_nworms,Prob_age_nworms_Adults[t,:,:].T) )
    
            HPrev_Inf[t,:]=np.dot( nbgt400, np.dot( Prob_pair_nworms,Prob_age_nworms_Inf[t,:,:].T) )
            HPrev_SAC[t,:]=np.dot( nbgt400, np.dot( Prob_pair_nworms,Prob_age_nworms_SAC[t,:,:].T) )
            HPrev_YAdults[t,:]=np.dot( nbgt400, np.dot( Prob_pair_nworms,Prob_age_nworms_YAdults[t,:,:].T) )
            HPrev_Adults[t,:]=np.dot( nbgt400, np.dot( Prob_pair_nworms,Prob_age_nworms_Adults[t,:,:].T) )
    
            MeanInt_Inf[t,:]=rho0_i*(phi_k*Inf[t,:]).sum()
            MeanInt_SAC[t,:]=rho0_i*(phi_k*SAC[t,:]).sum()
            MeanInt_YAdults[t,:]=rho0_i*(phi_k*YAdults[t,:]).sum()
            MeanInt_Adults[t,:]=rho0_i*(phi_k*Adults[t,:]).sum()

        if s==i_samples[id_samples[0]]:
            Prev_Inf_samples = Prev_Inf
            Prev_SAC_samples = Prev_SAC
            Prev_YAdults_samples = Prev_YAdults
            Prev_Adults_samples = Prev_Adults
    
            HPrev_Inf_samples = HPrev_Inf
            HPrev_SAC_samples = HPrev_SAC
            HPrev_YAdults_samples = HPrev_YAdults
            HPrev_Adults_samples = HPrev_Adults
    
            MeanInt_Inf_samples = MeanInt_Inf
            MeanInt_SAC_samples = MeanInt_SAC
            MeanInt_YAdults_samples = MeanInt_YAdults
            MeanInt_Adults_samples = MeanInt_Adults
    
            ystate_samples = sstate[:,0]
            zstate_samples = sstate[:,1]

            Incidence_samples = Incidence
            Mean_worm_SAC_samples = Mean_worm_SAC
            var_worm_SAC_samples = var_worm_SAC 
            Aggregation_SWB_SAC_samples = Aggregation_SWB_SAC 
            SAC_eqwdist_samples = SAC_eqwdist 
        else:
            #print Prev_Inf_samples.shape
            Prev_Inf_samples = np.concatenate((Prev_Inf_samples,Prev_Inf), axis=1)
            Prev_SAC_samples = np.concatenate((Prev_SAC_samples,Prev_SAC), axis=1)
            Prev_YAdults_samples = np.concatenate((Prev_YAdults_samples,Prev_YAdults), axis=1)
            Prev_Adults_samples = np.concatenate((Prev_Adults_samples,Prev_Adults), axis=1)

            HPrev_Inf_samples = np.concatenate((HPrev_Inf_samples,HPrev_Inf), axis=1)
            HPrev_SAC_samples = np.concatenate((HPrev_SAC_samples,HPrev_SAC), axis=1)
            HPrev_YAdults_samples = np.concatenate((HPrev_YAdults_samples,HPrev_YAdults), axis=1)
            HPrev_Adults_samples = np.concatenate((HPrev_Adults_samples,HPrev_Adults), axis=1)

            MeanInt_Inf_samples = np.concatenate((MeanInt_Inf_samples,MeanInt_Inf), axis=1)
            MeanInt_SAC_samples = np.concatenate((MeanInt_SAC_samples,MeanInt_SAC), axis=1)
            MeanInt_YAdults_samples = np.concatenate((MeanInt_YAdults_samples,MeanInt_YAdults), axis=1)
            MeanInt_Adults_samples = np.concatenate((MeanInt_Adults_samples,MeanInt_Adults), axis=1)
            Aggregation_SWB_SAC_samples = np.append(Aggregation_SWB_SAC_samples, Aggregation_SWB_SAC) 
            Mean_worm_SAC_samples = np.append(Mean_worm_SAC_samples,Mean_worm_SAC)
            var_worm_SAC_samples = np.append(var_worm_SAC_samples,var_worm_SAC)  

            if np.ndim(ystate_samples)==1:
                ystate_samples = np.concatenate((ystate_samples[:,None], sstate[:,0][:,None]), axis=1)
                zstate_samples = np.concatenate((zstate_samples[:,None], sstate[:,1][:,None]), axis=1)
                Incidence_samples = np.concatenate((Incidence_samples[:,None], Incidence[:,None]), axis=1)
                SAC_eqwdist_samples = np.concatenate((SAC_eqwdist_samples[:,None], SAC_eqwdist[:,None]), axis=1) 
            else:
                ystate_samples = np.concatenate((ystate_samples, sstate[:,0][:,None]), axis=1)
                zstate_samples = np.concatenate((zstate_samples, sstate[:,1][:,None]), axis=1)
                Incidence_samples = np.concatenate((Incidence_samples, Incidence[:,None]), axis=1)
                SAC_eqwdist_samples = np.concatenate((SAC_eqwdist_samples, SAC_eqwdist[:,None]), axis=1) 
      
        
        
    TobeSavedDic= {'times':times, \
                   'SAC_eqwdist_samples':SAC_eqwdist_samples, \
                   'Aggregation_SWB_SAC_samples':Aggregation_SWB_SAC_samples, \
                   'Mean_worm_SAC_samples': Mean_worm_SAC_samples, \
                   'var_worm_SAC_samples':var_worm_SAC_samples, \
                   'Incidence_samples':Incidence_samples, \
                   'ystate_samples':ystate_samples, \
                   'zstate_samples':zstate_samples, \
                   'Prev_Inf_samples':Prev_Inf_samples, \
                   'Prev_SAC_samples':Prev_SAC_samples, \
                   'Prev_YAdults_samples':Prev_YAdults_samples, \
                   'Prev_Adults_samples':Prev_Adults_samples, \
                   'HPrev_Inf_samples':HPrev_Inf_samples, \
                   'HPrev_SAC_samples':HPrev_SAC_samples, \
                   'HPrev_YAdults_samples':HPrev_YAdults_samples, \
                   'HPrev_Adults_samples':HPrev_Adults_samples, \
                   'MeanInt_Inf_samples':MeanInt_Inf_samples, \
                   'MeanInt_SAC_samples':MeanInt_SAC_samples, \
                   'MeanInt_YAdults_samples':MeanInt_YAdults_samples, \
                   'MeanInt_Adults_samples':MeanInt_Adults_samples}


    filename = "Results_Prevgt50_Scen_"+str(ScenN)+".p"
    outfile = open(filename, 'wb')
    pickle.dump(TobeSavedDic, outfile)
    outfile.close()
    return TobeSavedDic

#allow this code to be run as a script
if __name__ == "__main__":
    pfile=main()
