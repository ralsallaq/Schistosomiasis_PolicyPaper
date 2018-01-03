# Schistosomiasis_PolicyPaper
Code For High Prevalence Communities

Schistosomiasis is a parasitic disease with the parasitic life cycle completed in infected humans and in specific species of snails. Three species of schistosomiasis are of clinical importance: mansoni, haematobium, and japonicum. Infected persons with S.mansoni shed eggs in stools. Schistosomiasis can result in anaemia chronic pain, diarrhea and malnutrition. The WHO has drafted guidelines in 2012/2013 for controlling the morbidity from schistosomiasis. 
The analyses are for evaluating whether the WHO guidelines concerning control of schistosomiasis-induced morbidities would be sufficient to attain WHO morbidity-control targets. The guidelines specify a coverage of 75% of preventive chemotherapy among school-aged children (SAC) applied at years 0 for 5-6 years with frequency that increases with baseline prevalence among SAC. Then in year 6 the frequency is adapted according to prevalence attained with higher frequency if prevalence is still particularly high. Prevalence decision thresholds are high prevalence (>=50%) moderate prevalence (>=10% and <50%) and low prevalence (<10%) among SAC. The WHO morbidity-control targets pertains to levels of heavy prevalence among SAC (proportion with 400 eggs/gram of stool or more for S.mansoni) are to reach 5% or less in 6 years and 1% or less in 10 years.  
The steps are: 
1. Calibrate the dynamical model to reflect baseline prevalence in these settings and to account for the influence of the associated uncertainties in specifying the epidemiological parameters.
2. Apply WHO guidelines accordingly on the ensemble obtained in 1 and study the probabilities of attaining WHO targets.
3. Assess the adequacy of WHO guidelines in high prevalence settings.
Files:
* The file MCMC.p contains the parameters specified by calibrating the dynamical model (human coupled to snails) against SCORE experience in high prevalence villages. The calibration was done by a combination of Markov Chain fitting to data from a selected village of high prevalence and heavy prevalence among SAC followed by adjusting both relative rates of transmission of age groups and the aggregation parameter specifying the variation in egg release to reflect prevalence and heavy prevalence among SAC in a number of SCORE villages with prevalence and heavy prevalence among SAC in these villages measured to be >50% and >10%  (see file Prev_HPrev_villages.txt for processed data for these villages). Adjusting these two parameters is justified by being the least certain as reflected by experimenting with the MCMC fitting.

Please note that the raw data for the villages cannot be shared without data sharing agreement. The raw data was processed using python pandas to get age-specific profiles for prevalence and heavy prevalence.  

* The code Policy_Analysis_prevgt50.py is tasked of obtaining an ensemble of time series for high prevalence communities and for different scenarios. The following output files can be obtained using this code (D1=decision at year 6, D2=decision at year 10):

Results_Prevgt50_Scen_1.p : 81% of simulations D1: annual, D2: annual
To get, do from 
python Policy_Analysis_prevgt50.py 1
Results_Prevgt50_Scen_2.p : 19% of simulations D1: annual, D2: twice per year
To get, do
python Policy_Analysis_prevgt50.py 2

19% of simulations from the ensemble have prevalence >50% in year 6 and are thus subjected to higher frequency.

* The code Policy_plots.py is tasked of analyzing the ensembles of different scenarios to get probabilities of attaining WHO goals for a given prevalence setting at different years. Also the code saves Excel files for the time series. Given the pickle files (the files Results_Prev.. ending with .p above), the following few lines explain how to run this code to obtain results:

To get a combined figure and an Excel file for high prev at WHO guidelines, do
python Policy_plots.py Results_Prevgt50_D1_Scen_1.p Results_Prevgt50_D1_Scen_2.p 1 0.5 1.0
* The Excel file output_gt50p_uptoY20.xlsx contains the results and plots for high prevalence settings.
