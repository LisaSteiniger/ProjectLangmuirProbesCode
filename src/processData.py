''' This file contains functions for processin Langmuir Probe data with the purpose of identifying shortening. 
    It uses functions from readData to fetch data from the archive and applies detection routines on it.'''

import itertools
import os

import numpy as np
import pandas as pd

import src.readData as read

#######################################################################################################################################################################        
def getDischargeIDsAndAttributes(dischargeIDs: list[str],
                                 campaigns: list[str] =['OP22', 'OP23'],
                                 configurations: list[str] =['EIM000+2520'],
                                 filterSelected: str ='',
                                 filesExist: bool =False,
                                 ) -> list[list[str], list[str], list[str]]:
    ''' This function creates a 1D list with all discharge IDs, and second and third list of identical shape with the corresponding campaigns and configurations
        -> dischargeIDs, corresponding configuartions, corresponding campaigns
        If the discharge ID list that is given is not empty, only campaigns are determined, configurations remain unknown
        -> the other input parameters are not important in that case, default values can be used
        If the given list with IDs is empty, then discharges are actually read out with configuration and campaign
        
        "dischargeIDs" is   either a list with all discharge IDs of interest that will be returned
                            or an empty list, indicating that all dischrges of certain configurations in certain campaigns should be used
        "campaigns" is a list with the campaigns from which discharges should be considered
        "configurations" is a list with the configurations from which discharges should be considered
        "filterSelected" is a string that provides the criteria for which discharges should be filtered from the logbook
        -> this is besides campaign and configuration criteria
        "filesExist" holds information over the reading of already downloaded data
        -> if it is true, existing files created by "read.readAllShotNumbersFromLogbook" are read out while missing files are created
            -> this only works if the filters have not be changed, otherwise the result will be a discharge list filtered by the old filter
        -> if it is false, possibly existing files are overwritten as all discharges are read again from the logbook'''
    
    if dischargeIDs == []:
        dischargeConfigurations = []
        dischargeCampaigns = []

        for campaign in campaigns:      
            for configuration in configurations:     
                #find all discharge IDs according to the filters activated in "filterSelected" 
                #usually filters by !conditioning, !gas valve tests, !sniffer tests, and configuration (internal filter of "read.readAllShotNumbersFromLogbook")
                print('Logbook search: reading discharges for ' + configuration)
                dischargeIDs.append(read.readAllShotNumbersFromLogbook(configuration, filterSelected, q_add=campaign, filesExist=filesExist)['dischargeID'])
                dischargeConfigurations.append([configuration]*len(dischargeIDs[-1]))
            
            dischargeCampaigns.append([campaign]*len(list(itertools.chain.from_iterable(dischargeIDs))))
        
        dischargeIDs = list(itertools.chain.from_iterable(dischargeIDs))
        dischargeConfigurations = list(itertools.chain.from_iterable(dischargeConfigurations))
        dischargeCampaigns = list(itertools.chain.from_iterable(dischargeCampaigns))

    else:
        dischargeConfigurations = ['unknown']*len(dischargeIDs)
        dischargeCampaigns = ['OP2.2' if x.startswith('2024') else 'OP2.3' for x in dischargeIDs]

    return dischargeIDs, dischargeConfigurations, dischargeCampaigns

#######################################################################################################################################################################        
def processLangmuirProbeData(LP_list: list[str],
                             dischargeIDs: list[str|float],
                             dischargeCampaigns: list[str],
                             dischargeConfigurations: list[str],
                             R_limit: int|float= 300,
                             V_min_ideal: int|float= -180,
                             V_max_ideal: int|float= 20,
                             V_tolerance: float|int= 0.05,
                             fetch: bool= True,
                             plottingRawData: bool= False,
                             safePlot: str= 'results/LP_',
                             filesExist: bool= False
                             ) -> None:
    ''' This function is the framework for fetching Langmuir Probe data from the archive and the application failure detection routines
        For each discharge, discharge ID, configuration, V_limit_failure, V_limit_change, and, R_limit failure are collected
        -> Results are written to a .csv file "results/LP_{LP}/{LP}_dischargePlungeList_FailureIndicators.csv"
        -> V_limit_failure gives the percentage of voltage extrema that is too high/too low for V_limit with tolerance during hold time
        -> V_limit_change is indicating a stron change in voltage limits during the plunge
        -> R_limit_failure gives the percentage of average resistance values of the probe line that is too high/too low for R_limit with tolerance during hold time
        
        
        "LP_list" holds the internal IDs of the Langmuir Probes e.g. '51222'
        -> see settings to get full set of Langmuir Probe IDs
        "dischargeIDs" is a list of IDs of the discharges e.g. "20241127.010"
        "dischargeCampaigns" is a list with the corresponding OP name for every discharge (same len as "dischargeIDs")
        "dischargeConfigurations" is a list with the corresponding configuration for every discharge (same len as "dischargeIDs")
        
        "R_limit" is the limit for averaged probe line resistance in (Ohm)
        -> if the averaged resistance drops below, shortening is assumed
        "V_min_ideal" is the lower limit of the ideal voltage range in (V)
        "V_max_ideal" is the upper limit of the ideal voltage range in (V)
        "V_tolerance" is the tolerance on the ideal voltage range
        -> e.g. 0.05 tolerance means that voltages < -171 V and > 19 V are ok
        "fetch" is True - only the data from the start of a plunge + 100 ms is fetched
                  is False, the whole discharge data is fetched
        "plottingRawData" decides if the raw data (probe voltage, probe current, and probe line resistance) are plotted
        "safePlot" is defalut path structure for saving a plot
        "filesExist" being False means that possibly existing .csv files (which were created by this function) are going to be overwritten'''

    for LP in LP_list:
        #if filesExist is True, existing files are read out, and only missing ones are created
        #if filesExist is False, existing files are overwritted, and missing ones are created
        if filesExist:  
            if os.path.isfile(f"results/LP_{LP}/{LP}_dischargePlungeList_FailureIndicators.csv"):
                print(f'"results/LP_{LP}/{LP}_dischargePlungeList_FailureIndicators.csv" exists')
                continue
                
        #lists with configurations and campaigns of the discharges and their plunges
        campaignListDict = []
        configurationListDict = []

        #lists with dischargeIDs and their plunges as lists -> each dischargeID is added x times with x being the number of plunges
        dischargeListDict = []
        plungeListDict = []

        #different indicators of shortening -> each discharge gets one entry (list with values for plunges)
        #if no shortening occurs, [0, False, 0] are entries of the three lists
        VlimitListDict = []    #percentage of voltage extrema that lie outside the ideal voltage range with tolerance (e.g. minima above V_min_ideal*V_tolerance)
        VchangeListDict = []   #do the extrema change significantly during a plunge? e.g. minima getting higher (=shortening occurs) or lower (=shortening disappears)
        RlimitListDict = []    #percentage of averaged R_probeline values that lie below the threshold of R_limit

        #values for electron density, electron temperature and their standard deviations
        neListDict = []
        TeListDict = []
        sneListDict = []
        sTeListDict = []

        for counter, dischargeID, dischargeCampaign, dischargeConfiguration in zip(range(len(dischargeIDs)), dischargeIDs, dischargeCampaigns, dischargeConfigurations):
            dischargeID = str(dischargeID)
            print(dischargeID)
    
            results = read.readLangmuirProbeOperationalParameters(dischargeID, LP, R_limit, V_min_ideal, V_max_ideal, V_tolerance, fetch, plottingRawData, safePlot)

            dischargeListDict.append(results[0])
            plungeListDict.append(results[1])
            VlimitListDict.append(results[2])
            VchangeListDict.append(results[3])
            RlimitListDict.append(results[4])
            neListDict.append(results[5])
            TeListDict.append(results[6])
            sneListDict.append(results[7])
            sTeListDict.append(results[8])
        
            campaignListDict.append([dischargeCampaign]*len(plungeListDict[-1]))
            configurationListDict.append([dischargeConfiguration]*len(plungeListDict[-1]))

            #in here, so that the intermediate results are saved after each 100 discharges
            if counter in range(0, 10001, 100):
                LP_Dict = pd.DataFrame({'LP':[LP]*len(list(itertools.chain.from_iterable(dischargeListDict))),
                                        'campaign': list(itertools.chain.from_iterable(campaignListDict)),
                                        'configuration': list(itertools.chain.from_iterable(configurationListDict)),
                                        'dischargeID': list(itertools.chain.from_iterable(dischargeListDict)),
                                        'plunge': list(itertools.chain.from_iterable(plungeListDict)),
                                        'V_limit_failure': list(itertools.chain.from_iterable(VlimitListDict)),
                                        'V_limit_change': list(itertools.chain.from_iterable(VchangeListDict)),
                                        'R_limit_failure': list(itertools.chain.from_iterable(RlimitListDict)),
                                        'ne': list(itertools.chain.from_iterable(neListDict)),
                                        'std_ne': list(itertools.chain.from_iterable(sneListDict)),
                                        'Te': list(itertools.chain.from_iterable(TeListDict)),
                                        'std_Te': list(itertools.chain.from_iterable(sTeListDict))})
                
                LP_Dict.to_csv(f'results/LP_{LP}/{LP}_dischargePlungeList_FailureIndicators.csv', sep=';')

#######################################################################################################################################################################        
def filterShorteningCandidatePlunges(LP_list: list[str],
                                     failureTolerance: int|float
                                     ) -> list[list[list[str],list[int]]]:
    ''' This function filters the complete list of discharges and plunges with their failure indicators for real candidates of shortening
        -> if V_limit_change is True or V_limit_failure or R_limit_failure percentage is outside of failure tolerance
        -> number of applying indicators is counted and addes to the dataframe
        -> filtered discharges and plunges are written to dataframe and saved as .csv file under f'results/LP_{LP}/{LP}failed_dischargePlungeList_FailureIndicators.csv' 
        returns list of discharge IDs that are probably subject to shortening and the max. number of indicators for all plunges of that discharge
        -> for every Langmuir Probe x, y, z,..., a tuple LPx/y/z = [IDs, indicatorNumber] is returned in [LPx, LPy, LPz, ...]

        "LP_list" is the list with Langmuir Probes that should be investigated (internal probe ID such as '50209')
        "failureTolerance" is the percentage of allowed outliers of voltage extrema and averaged probe line resistances without declaring the plunge a candidate for shortening'''

    failures = []
    for LP in LP_list:
        failureIndices = []
        numberOfFailureIndicatorsList = []
        V_limit_failure_yn_List = []
        R_limit_failure_yn_List = []
        #note that all columns are interpreted as str values
        shorteningIndicatorList = pd.read_csv(f"results/LP_{LP}/{LP}_dischargePlungeList_FailureIndicators.csv", sep=';', dtype={'dischargeID': str})

        for i in range(len(shorteningIndicatorList['dischargeID'])):
            numberOfFailureIndicators = 0
            V_limit_failure_yn = False
            R_limit_failure_yn = False

            #if no data was read out or no extrema were found, there is no valid data on the failure indicator columns
            if shorteningIndicatorList['V_limit_failure'][i] == 'no extrema':
                continue
            if np.isnan(float(shorteningIndicatorList['V_limit_failure'][i])):
                continue

            if float(shorteningIndicatorList['V_limit_failure'][i]) > failureTolerance:
                numberOfFailureIndicators += 1
                V_limit_failure_yn = True
            if shorteningIndicatorList['V_limit_change'][i] == 'True':
                numberOfFailureIndicators += 1
            if float(shorteningIndicatorList['R_limit_failure'][i]) > failureTolerance:
                numberOfFailureIndicators += 1
                R_limit_failure_yn = True
            if numberOfFailureIndicators != 0:
                failureIndices.append(i)
                numberOfFailureIndicatorsList.append(numberOfFailureIndicators)
                V_limit_failure_yn_List.append(V_limit_failure_yn)
                R_limit_failure_yn_List.append(R_limit_failure_yn)
            
        failureIndices = np.array(failureIndices)

        failureDataFrame = pd.DataFrame({'LP': np.array(shorteningIndicatorList['LP'])[failureIndices],
                                        'campaign': np.array(shorteningIndicatorList['campaign'])[failureIndices],
                                        'configuration': np.array(shorteningIndicatorList['configuration'])[failureIndices],
                                        'dischargeID': np.array(shorteningIndicatorList['dischargeID'])[failureIndices],
                                        'plunge': np.array(shorteningIndicatorList['plunge'])[failureIndices],
                                        'V_limit_failure': np.array(shorteningIndicatorList['V_limit_failure'])[failureIndices],
                                        'V_limit_failure_yn': V_limit_failure_yn_List,
                                        'V_limit_change': np.array(shorteningIndicatorList['V_limit_change'])[failureIndices],
                                        'R_limit_failure': np.array(shorteningIndicatorList['R_limit_failure'])[failureIndices],
                                        'R_limit_failure_yn': R_limit_failure_yn_List,
                                        'numberOfFailureIndicators': numberOfFailureIndicatorsList,
                                        'ne': np.array(shorteningIndicatorList['ne'])[failureIndices],
                                        'std_ne': np.array(shorteningIndicatorList['std_ne'])[failureIndices],
                                        'Te': np.array(shorteningIndicatorList['Te'])[failureIndices],
                                        'std_Te': np.array(shorteningIndicatorList['std_Te'])[failureIndices]})
    
        failureDataFrame = failureDataFrame.sort_values(by=["numberOfFailureIndicators", 'V_limit_change', 'R_limit_failure', 'V_limit_failure'], ascending=False)
        failureDataFrame.to_csv(f'results/LP_{LP}/{LP}failed_dischargePlungeList_FailureIndicators.csv', sep=';') 

        failureDischargeIDs = np.unique(np.array(failureDataFrame['dischargeID']))
        maxFailureNumbersDischargeIDs = []

        for failureDischargeID in failureDischargeIDs:
            failureFilter = np.array([x == failureDischargeID for x in failureDataFrame['dischargeID']])
            maxFailureNumbersDischargeIDs.append(np.max(np.array(failureDataFrame['numberOfFailureIndicators'])[failureFilter]))
        
        failures.append([failureDischargeIDs, maxFailureNumbersDischargeIDs])
    
    return failures