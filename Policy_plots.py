import sys, os
import time as timeM
from scipy.misc import comb
from scipy import stats
import numpy as np
import matplotlib as mpl
""" set up fonts here before importing matplotlib.pylab """
parms = {'font.family': 'serif', 'font.serif': 'Palatino','font.size':23,'xtick.major.size':6,'ytick.major.size':6}
mpl.rcParams.update(parms)
from matplotlib import pyplot as pl
from scipy.stats.mstats import mquantiles
import ObtainDemography_optimized as DemCode
reload(DemCode)
import cPickle as pickle
from scipy.stats.mstats import mquantiles
import pandas as pd
"""
This code is developed for part of the analyses in the paper:
Jaspreet Toor, Ramzi Alsallaq, James Truscott, Hugo Turner, Marleen Werkman, et al: Are we on our way to achieving the 2020 goals for schistosomiasis morbidity control using current WHO guidelines?, to appear at Clinical Infectious Diseases 2017.
For questions please contact Ramzi Alsallaq (rxa313@case.edu or ramzi.alsallaq@gmail.com)
Refer to README for more information on how to run this code

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
        #was 0-4, 5-8, 9-12, 13+ 
        #Now 0-4, 5-14, 15-29, 30+
        self.params['Mortrate']= np.array([0.074, 0.007, 0.025, 0.06])
        self.params['Matrate']=1./np.array([5., 10., 15., 70.])
        self.params['Treat_Eff'] = 0.863 #as per the meta-analysis by Zwang et.al 2014 (egg reduction rate)
    
        self.params['RelLamb'] = np.array([0.2, 0.53, 0.010486887, 0.010486887])
        self.params['RelRho'] = np.array([1., 1., 1., 1.0])
        Events = {}
        TotCov = np.array([0.75]*10)
        for i in xrange(1,21):
            Events[i] = {'groups':np.array([0,1,2,3]),'coverage':np.array(4*[TotCov[0]]), 'compliance':np.array(4*[1.0])}
        self.params.update({'Events':Events})
    
        self.params.update(**kwargs)
        self.params['gamma']=(1./5.7)*np.ones(self.params['Nage'])

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


def PlotScenFig(pfile=None, ScenN=None, PrevTh=None):
    times = pfile['times']
    Prev_SAC_samples = pfile['Prev_SAC_samples']
    ind, = np.where((Prev_SAC_samples[0]<PrevTh[1])&(Prev_SAC_samples[0]>=PrevTh[0]))
    Prev_SAC_samples = Prev_SAC_samples[:,ind]
    HPrev_SAC_samples = pfile['HPrev_SAC_samples'][:,ind]
    Incidence_samples = pfile['Incidence_samples'][:,ind]
    fig1,ax=pl.subplots()
    #ax.set_title("Infection levels among SAC")
    SACPrev_mean =Prev_SAC_samples.mean(axis=1)*100.
    SACPrev_qs = mquantiles(Prev_SAC_samples,[0.025,0.975], axis=1)*100.
    SACHPrev_mean =HPrev_SAC_samples.mean(axis=1)*100.
    SACHPrev_qs = mquantiles(HPrev_SAC_samples,[0.025,0.975], axis=1)*100.
    #incidence
    PopInc_mean =Incidence_samples.mean(axis=1)*100.
    #print PopInc_mean
    PopInc_qs = mquantiles(Incidence_samples,[0.025,0.975], axis=1)*100.

    if PrevTh[0]==0.5:
        fh = "gt50p"
        ar1_h=(6,67); ar2_h=(0,67)
        textat =(2.,77)
        #Processed data from SCORE villages
        try:
            DATA=pd.read_table("Prev_HPrev_villages.txt",sep=" ")
            DATA_villages = DATA.Village_ID.values 
            DATA_SACprev = DATA.SACPrev.values*100. 
            DATA_SACHprev = DATA.SACHeavyPrev*100. 
            DATA_time = np.array([-0.5]*DATA_villages.shape[0])
            ax.scatter(DATA_time[:,None], DATA_SACprev[:,None], s=100, c='k', marker="*",alpha=0.6, label="Village data")
            ax.scatter(DATA_time[:,None], DATA_SACHprev[:,None], s=100,c='k', marker="*", alpha=0.6)
        except:
            pass

    if PrevTh[0]==0.1:
        fh = "gt10p"
        ar1_h=(6,33); ar2_h=(0,50)
        textat =(2.,77)
    if PrevTh[0]==0.:
        fh = "lt10p"
        ar1_h=(6,5); ar2_h=(0,10)
        textat =(2.,77)
    """save to excel"""
    filename="output_"+fh+".xlsx"
    df = pd.DataFrame({'Year':times-1,'SACPrev-Mean (%)':SACPrev_mean,'SACPrev-0.025q (%)':SACPrev_qs.T[0], 'SACPrev-0.975q (%)':SACPrev_qs.T[1], 'SACHPrev-Mean (%)':SACHPrev_mean, 'SACHPrev-0.025q (%)':SACHPrev_qs.T[0],'SACHPrev-0.975q (%)':SACHPrev_qs.T[1], 'PopIncidence-Mean (w/100 p-y)':PopInc_mean, 'PopIncidence-0.025q (w/100 p-y)': PopInc_qs.T[0], 'PopIncidence-0.975q (w/100 p-y)': PopInc_qs.T[1]})
    writer = pd.ExcelWriter(filename)
    df.to_excel(writer,'DataPrev'+fh, columns=['Year','SACPrev-Mean (%)','SACPrev-0.025q (%)', 'SACPrev-0.975q (%)', 'SACHPrev-Mean (%)', 'SACHPrev-0.025q (%)','SACHPrev-0.975q (%)', 'PopIncidence-Mean (w/100 p-y)', 'PopIncidence-0.025q (w/100 p-y)','PopIncidence-0.975q (w/100 p-y)'])
    writer.save()

    ax.plot(times-1.0,SACPrev_mean, linewidth=3, label="Mean prevalence")
    ax.fill_between(times-1.0, *SACPrev_qs.T , color='k', alpha=0.3, label="95% CI")
    ax.plot(times-1.0,SACHPrev_mean, linewidth=3, label="Mean heavy-intensity prevalence")
    ax.fill_between(times-1.0, *SACHPrev_qs.T , color='k', alpha=0.3)
    ax.plot([times[0]-1.0,times[-1]-1.0],[5.,5.],"--k",linewidth=1, label="2020 target")
    ax.plot([times[0]-1.0,times[-1]-1.0],[1.,1.],":k",linewidth=1, label="2025 target")
    ax.annotate('decision points',xy=ar1_h,xytext=textat, arrowprops=dict(facecolor='k',shrink=0.05, width=0.5),)
    ax.annotate('',xy=ar2_h,xytext=textat, arrowprops=dict(facecolor='k',shrink=0.05, width=0.5),)
    ax.set_xlabel("Years since intervention starts")
    ax.set_ylabel("Percentage among SAC")
    ax.set_ylim(0,100)
    ax.set_xlim(-1,16)
    ax.legend(loc=0)

    #incidence
    fig2,ax=pl.subplots()
    ax.fill_between(times-1.0, *PopInc_qs.T , color='k', alpha=0.3, label="95% CI")
    ax.plot(times-1, PopInc_mean, linewidth=3, label="Mean")
    ax.plot([times.min()-6, times.max()+20],[1.0, 1.0], "--k", linewidth=3, label='transmission interruption threshold')
    ax.set_ylabel("Population incidence (worm/100 person-year)", fontsize=20)
    ax.set_xlabel("Years since intervention starts")
    ax.set_xlim(-1,16)
    ax.set_ylim(0,260)
    ax.legend(loc=0)

    return fig1, fig2

def PrintSummariesAtYear(indY, indScen, time, pfile):
    print "-----------------","Year"+str(time[indY])+"-----------"
    print time[indY]
    MY_SACPrev=pfile['Prev_SAC_samples'][indY,indScen].mean()
    QsY_SACPrev=mquantiles(pfile['Prev_SAC_samples'][indY,indScen], [0.025,0.975])
    print "mean Y sac prev=",MY_SACPrev
    print "95% CI Y sac prev=",(QsY_SACPrev[0], QsY_SACPrev[1])
    MY_SACHPrev=pfile['HPrev_SAC_samples'][indY,indScen].mean()
    QsY_SACHPrev=mquantiles(pfile['HPrev_SAC_samples'][indY,indScen], [0.025,0.975])
    print "mean Y sac Hprev=",MY_SACHPrev
    print "95% CI Y sac Hprev=",(QsY_SACHPrev[0], QsY_SACHPrev[1])
    PercSACPrevLt1p = (pfile['Prev_SAC_samples'][indY,indScen]<0.01).sum()/float(pfile['Prev_SAC_samples'][indY,indScen].shape[0])
    PercSACPrevBt1pA10p = ((pfile['Prev_SAC_samples'][indY,indScen]>=0.01)&(pfile['Prev_SAC_samples'][indY,indScen]<0.1)).sum()/float(pfile['Prev_SAC_samples'][indY,indScen].shape[0])
    PercSACPrevBt10pA50p = ((pfile['Prev_SAC_samples'][indY,indScen]>=0.1)&(pfile['Prev_SAC_samples'][indY,indScen]<0.5)).sum()/float(pfile['Prev_SAC_samples'][indY,indScen].shape[0])
    PercSACPrevGe50p = (pfile['Prev_SAC_samples'][indY,indScen]>=0.5).sum()/float(pfile['Prev_SAC_samples'][indY,indScen].shape[0])
    #Prob to achieve goals in Y
    PercSACHPrevLt5p = (pfile['HPrev_SAC_samples'][indY,indScen]<0.05).sum()/float(pfile['HPrev_SAC_samples'][indY,indScen].shape[0])

    #Prob to achieve goals in Y
    PercSACHPrevLt5p = (pfile['HPrev_SAC_samples'][indY,indScen]<0.05).sum()/float(pfile['HPrev_SAC_samples'][indY,indScen].shape[0])
    print "Probability to achieve goal 1 in Y = ", PercSACHPrevLt5p
    PercSACHPrevLt1p = (pfile['HPrev_SAC_samples'][indY,indScen]<0.01).sum()/float(pfile['HPrev_SAC_samples'][indY,indScen].shape[0])
    print "Probability to achieve goal 2 in Y = ", PercSACHPrevLt1p

    PopIncY=pfile['Incidence_samples'][indY,indScen]
    print "Probability of elimination (incidence<0.01) in Y=",(PopIncY<0.01).sum()/float(PopIncY.shape[0])
    print "--------------------------------------------------------------------"
    return PercSACPrevLt1p, PercSACPrevBt1pA10p, PercSACPrevBt10pA50p, PercSACPrevGe50p



def main(argv=None):
    if argv is None:
        argv = sys.argv

    filename1 = argv[1]
    filename2 = argv[2]
    ScenN = int(argv[3])
    print ScenN
    PrevTh = (float(argv[4]),float(argv[5]))

    setup = MDAModelRunSetUp()
    setup.GetDemography()
    AgeFrac = setup.AgeFrac
    equi_h = setup.equi_h
    phi_k = setup.phi()
    nworms, npairs = setup.nwroms_npairs()
    params=setup.params
    Pop0 = params['Pop0']
    deldblu = params['deldblu']
    Nage = params['Nage']
    RelLamb = params['RelLamb']
    Mortrate = params['Mortrate']
    gamma = params['gamma']    
    epsi = params['epsi']
    Nstrata=setup.Nstrata()
    
    """ specify params """
    pfile=pickle.load(open("MCMC.p","rb"))
    rho0 = pfile["rho0"][0][600:]
    #I do not read alphas and betas directly because these were generated for a different age grouping
    logbeta0 = pfile["logbeta0"][0][600:]
    logalphaSAC = pfile["logalphaSAC"][0][600:]
    alphas = np.power(10,logalphaSAC)[:,None]*RelLamb[None,:] 
    betas = alphas*np.power(10,logbeta0)[:,None]

    #Calculating Rnot
    #Averaging alpha and beta for one age group
    alpha = (alphas*AgeFrac).sum(axis=1)
    beta = (betas*AgeFrac).sum(axis=1)
    mu = (Mortrate*AgeFrac).sum()
    gamma = gamma[0]
    Rnot = np.sqrt(Pop0*alpha*beta*rho0*phi_k[1]/(epsi*(gamma+mu)))
    print "full number of samples=", Rnot.shape[0]
    Rpairs = zip(np.arange(1,5,0.5),np.arange(1.5,6,0.5))
    i_samples=np.array([], dtype=int)
    for Rmn, Rmx in Rpairs:
        temp = np.where((Rnot>=Rmn)&(Rnot<Rmx))[0][-80:] 
        i_samples = np.append(i_samples, temp)
    nsamples = i_samples.shape[0]
    print "total number of samples", nsamples 

    #read file
    outfile = open(filename1, 'rb')
    pfile1 = pickle.load(outfile)
    outfile.close()
    outfile = open(filename2, 'rb')
    pfile2 = pickle.load(outfile)
    outfile.close()
    pfile={}
    for key in pfile1:
        pfile[key] = np.hstack((pfile1[key],pfile2[key]))
    #make sure time is not dublicated
    pfile['times']=pfile1['times']

    PlotScenFig(pfile=pfile, ScenN=ScenN, PrevTh=PrevTh)

    #Getting an idea about the aggregation in the equilibrium SWB system:
    Aggregation_SWB_SAC_samples = pfile['Aggregation_SWB_SAC_samples']
    print "Aggregation of SWB eqm distribution for SAC=",mquantiles(Aggregation_SWB_SAC_samples,[0.025,0.5,0.975])
    time = pfile['times']-1.0 #shifted time to start from zero
    ind0 = np.where(time<=0)[0][0:10]
    ind, =np.where((pfile['Prev_SAC_samples'][0]>=PrevTh[0]) & (pfile['Prev_SAC_samples'][0]<PrevTh[1]))
    print "number of scenarios=",pfile['Prev_SAC_samples'][:,ind].shape[1]
    MBase_SACPrev=pfile['Prev_SAC_samples'][ind0][:,ind].mean(axis=1)[0]
    QsBase_SACPrev=mquantiles(pfile['Prev_SAC_samples'][ind0][:,ind], [0.025,0.975],axis=1).T
    print "mean BL sac prev=",MBase_SACPrev
    print "95% CI BL sac prev=",(QsBase_SACPrev[0].mean(), QsBase_SACPrev[1].mean())
    MBase_SACHPrev=pfile['HPrev_SAC_samples'][ind0][:,ind].mean(axis=1)[0]
    QsBase_SACHPrev=mquantiles(pfile['HPrev_SAC_samples'][ind0][:,ind], [0.025,0.975],axis=1).T
    print "mean BL sac Hprev=",MBase_SACHPrev
    print "95% CI BL sac Hprev=",(QsBase_SACHPrev[0].mean(), QsBase_SACHPrev[1].mean())

    #Y6 analyses
    ind6 = np.where(time<6)[0][-1]
    print time[ind6]
    PercSACPrevLt1p, PercSACPrevBt1pA10p, PercSACPrevBt10pA50p,PercSACPrevGe50p =PrintSummariesAtYear(ind6, ind, time, pfile)
    #Decision 2
    ProbNextDecision=0.0
    if PercSACPrevLt1p>ProbNextDecision:
        ProbNextDecision = PercSACPrevLt1p
        printable = "Conduct serology"
    if PercSACPrevBt1pA10p>ProbNextDecision:
        ProbNextDecision = PercSACPrevBt1pA10p
        printable = "PCT once every two years"
    if PercSACPrevBt10pA50p>ProbNextDecision:
        ProbNextDecision = PercSACPrevBt10pA50p
        printable = "PCT at previous frequency"
    if PercSACPrevGe50p>ProbNextDecision:
        ProbNextDecision = PercSACPrevGe50p 
        printable = "PCT two times a year"
    
    print "Decision is:", printable, "with probability:",ProbNextDecision  

    #Y7 end point analyses
    ind7 = np.where(time<7)[0][-1]
    print time[ind7]
    PercSACPrevLt1p, PercSACPrevBt1pA10p, PercSACPrevBt10pA50p,PercSACPrevGe50p =PrintSummariesAtYear(ind7, ind, time, pfile)

    #Y10 end point analyses
    ind10 = np.where(time<10)[0][-1]
    print time[ind10]
    PercSACPrevLt1p, PercSACPrevBt1pA10p, PercSACPrevBt10pA50p,PercSACPrevGe50p =PrintSummariesAtYear(ind10, ind, time, pfile)


    #Y11 end point analyses
    ind11 = np.where(time<11)[0][-1]
    print time[ind11]
    PercSACPrevLt1p, PercSACPrevBt1pA10p, PercSACPrevBt10pA50p,PercSACPrevGe50p =PrintSummariesAtYear(ind11, ind, time, pfile)

    #Y15 end point analyses
    ind15 = np.where(time<15)[0][-1]
    print time[ind15]
    PercSACPrevLt1p, PercSACPrevBt1pA10p, PercSACPrevBt10pA50p,PercSACPrevGe50p =PrintSummariesAtYear(ind15, ind, time, pfile)
    

    #Y16 end point analyses
    ind16 = np.where(time<16)[0][-1]
    print time[ind16]
    PercSACPrevLt1p, PercSACPrevBt1pA10p, PercSACPrevBt10pA50p,PercSACPrevGe50p =PrintSummariesAtYear(ind16, ind, time, pfile)

    #end year goal2 acheived?
    indf = np.where(time<time[-1])[0][-1]
    print time[indf]
    PercSACPrevLt1p, PercSACPrevBt1pA10p, PercSACPrevBt10pA50p,PercSACPrevGe50p =PrintSummariesAtYear(indf, ind, time, pfile)


#allow this code to be run as a script
if __name__ == "__main__":
    try:
        main()
    except:
        print "Usage: %run Policy_plots.py picklefile1 picklefile2 scenN prevth1 prevth2"
        print "Usage example: %run Policy_plots.py Results_Prevgt50_D1_Scen_1.p Results_Prevgt50_D1_Scen_2.p 1 0.5 1.0"
