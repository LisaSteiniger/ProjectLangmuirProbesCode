""" This project deals with the shortening of the pop-up Langmuir Probes of Wendelstein 7-X in OP2.2 and OP2.3
    It enables automatic detection of plunges that have been subject to shortening by investigating the range of voltage
    and the probe line resistance. 
    Additionally, the operational and plasma parameters of measurements and discharges which experienced those shortenings
    are tested for any kind of systematics including same magnetic field configuration, similar iota, P_rad, ...
    This file is the actual application file. It calls the neccessary functions that are saved in the src folder. 
    General settings are stored in settings.py, input Files in the inputFile folder and files containing results
    in the result folder"""

import numpy as np
import pandas as pd

import os
import itertools
import matplotlib

import settings
import src.readData as read
import src.plotData as plot
import src.processData as process 

#avoids pop up of plot windows
matplotlib.use('agg')

#######################################################################################################################################################################        
#parameters for running program:

#list of all released configurations (date: 01. Dec. 2025) from W7X-info
#-> alternativly, list of strings can be provided with a set/single configuration, e.g. ['EIM000+2520'] 
configurations = pd.read_csv('inputFiles/configurationListWithSettings.csv', sep=';')['configuration']
#configurations.remove('DBM000+2520')
#configurations.remove('FTM004+2520')
#configurations.remove('KJM008+2520')
#configurations.remove('EIM000+2520')
#configurations = ['DBM000+2520', 'FTM004+2520', 'KJM008+2520']

#list of campaigns to be looked at, 'OP223' means OP2.2 and OP2.3
campaigns = ['OP22' ,'OP23']    #['OP223', 'OP22', 'OP23']

#set filter options for the discharges here (configuration filter is applied later)
filterSelected = settings.q1 + settings.q2 + settings.q3

#Langmuir Probes that should be investigated
LP_list = ['50209', '51222'] #['50209']

#discharges to be looked at
#leave blank if list is given by campaign and configuration
dischargeIDs = [] #['20241127.010', '20241105.067']


#limits of voltage range and averaged probe line resistance in (V) and (Ohm)
V_min_ideal=-180
V_max_ideal=20
R_limit=50

#tolerance for voltage range and outliers of voltage extrema and average probe line resistance
V_tolerance = 0.05          #voltage extrema are ok if in [V_min_ideal + (V_max_ideal - V_min_ideal)*V_tolerance, V_max_ideal - (V_max_ideal - V_min_ideal)*V_tolerance]
outliers_tolerance = 0.03   #percentage of voltage extrema and average values that might lie outside of the V_tolarance/R_limit withuot being interpreted as potential shortening

#are the discharge lists per configuration/Langmuir Probe already saved (at least partially, missing files will be created anyways)? 
#-> inputFiles/configurations/dischargeList*.csv and results/LP_*/*_dischargePlungeList_FailureIndicators.csv
#set ConfigurationExist to False only if you want to reset your whole data set (e.g. when having changed the filter for the discharges)
#set LangmuirProbeExist to False if list of discharge IDs has changed
filesDischargesPerConfigurationExist = True                   
filesDischargesPerLangmuirProbeExist = False                  

#should the measurement values be read out again?
#set True only if there was a problem with the reading routine
#if False, then already read out data wont be downloaded again
reReadData = False

#should probe voltage, probe current, and probe line resistance be plotted?
plottingRawData = False

#should the main program run or is just some other testing going on? -> only commands above "if not run:" will be executed
run = True

#######################################################################################################################################################################        
#RUNNING PROGRAM -> DO NOT CHANGE 

if not run:
    exit()

#read dischargeIDs, get corresponding configurations and campaigns
dischargeIDs, dischargeConfigurations, dischargeCampaigns = process.getDischargeIDsAndAttributes(dischargeIDs, campaigns, configurations, filterSelected, filesDischargesPerConfigurationExist)

#get Langmuir Probe data and test it for indicators of shortening, results are saved as .csv file for each Langmuir Probe
#-> saved under f"results/LP_{LP}/{LP}_dischargePlungeList_FailureIndicators.csv"
process.processLangmuirProbeData(LP_list, dischargeIDs, dischargeCampaigns, dischargeConfigurations, R_limit, V_min_ideal, V_max_ideal, V_tolerance, fetch=True, plottingRawData=plottingRawData, filesExist=filesDischargesPerLangmuirProbeExist)

#test indicators for shortening, and select discharges and plunges that actually might have been subject to shortening while discarding working data 
#results are saved in .csv file under f'results/LP_{LP}/{LP}failed_dischargePlungeList_FailureIndicators.csv' 
process.filterShorteningCandidatePlunges(LP_list, outliers_tolerance)

'''
failureDataFrame = pd.read_csv(f'results/LP_50209/50209failed_dischargePlungeList_FailureIndicators.csv', sep=';', dtype={'dischargeID': str}) 

failureDischargeIDs = np.unique(np.array(failureDataFrame['dischargeID']))
maxFailureNumbersDischargeIDs = []

for failureDischargeID in failureDischargeIDs:
    failureFilter = np.array([x == failureDischargeID for x in failureDataFrame['dischargeID']])
    maxFailureNumbersDischargeIDs.append(np.max(np.array(failureDataFrame['numberOfFailureIndicators'])[failureFilter]))

'''

plot.plotFailuresInDependanceOfX('configuration', '50209', 'OP223')
plot.plotFailuresInDependanceOfX('configuration', '51222', 'OP223')
plot.plotFailuresInDependanceOfX('configurationShort', '50209', 'OP223')
plot.plotFailuresInDependanceOfX('configurationShort', '51222', 'OP223')