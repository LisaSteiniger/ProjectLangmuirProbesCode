''' This file contains plotting fuctions for visualization of possilbe systematics in shortening'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from settings import legendsize, scalesize, labelsize

#######################################################################################################################################################################
def plotFailuresInDependanceOfX(X: str,
                                LP: str,
                                OP: str,
                                pathAllDischarges: str ='',
                                pathFailureCandidates: str ='',
                                pathFailureDischarges: str =''
                                ) -> None:
    ''' This function plots the relative number of failures in voltage limit, voltage limit change, and probe line resistance limit
        -> in dependence of a characteristic X for a certain campaign
        
        resulting plot is saved under "results/LP_{LP}/failuresInDependanceOf{X}.png}
        
        "X" is the characteristic that should be investigated 
        -> "configurationShort" e.g. "EIM"
        -> "configuration" e.g. "EIM000+2520"

        "LP" is the investigated Langmuir Probes internal ID, e.g. "50209"
        "OP" is the campaign that should be investigated
        -> 'OP22', 'OP23', 'OP223'
        "pathAllDischarges" is the path to the overview file with the Langmuir Probe failure indicators for all discharges
        "pathFailureCandidates" is the path to the overview file with the Langmuir Probe failure indicators for failure candidates
        "pathFailureDischarges" is the path to the overview file with the discharges for which the Langmuir Probe data showed at least one failure indicator
        -> contains max number of positive failure indicators for all plunges and the information if the discharge showed a failed plunge'''
    
    if pathAllDischarges == '':
        pathAllDischarges = f'results/LP_{LP}/{LP}_dischargePlungeList_FailureIndicators.csv'
    
    if pathFailureCandidates == '':
        pathFailureCandidates = f'results/LP_{LP}/{LP}failed_dischargePlungeList_FailureIndicators.csv'
    
    if pathFailureDischarges == '':
        pathFailureDischarges = f'results/LP_{LP}/{LP}failed_dischargeList_maxFailureIndicators.csv'

    allDischargesDataFrame = pd.read_csv(pathAllDischarges, sep=';')
    failureCandidatesDataFrame = pd.read_csv(pathFailureCandidates, sep=';', dtype={'V_limit_failure_yn': bool, 'R_limit_failure_yn': bool, 'V_limit_change': bool})
    failureDischargesDataFrame = pd.read_csv(pathFailureDischarges, sep=';', dtype={'numberOfFailedPlunges': int, 'maxNumberOfFailureIndicators': int})

    if OP == 'OP22':
        allDischargeFilter = np.array([x == 'OP22' for x in allDischargesDataFrame['campaign']])
        failureCandidatesFilter = np.array([x == 'OP22' for x in failureCandidatesDataFrame['campaign']])
        failureDischargesFilter = np.array([x == 'OP22' for x in failureDischargesDataFrame['campaign']])
    elif OP == 'OP23':
        allDischargeFilter = np.array([x == 'OP23' for x in allDischargesDataFrame['campaign']])
        failureCandidatesFilter = np.array([x == 'OP23' for x in failureCandidatesDataFrame['campaign']])
        failureDischargesFilter = np.array([x == 'OP23' for x in failureDischargesDataFrame['campaign']])
    elif OP == 'OP223':
        allDischargeFilter = np.array([x == 'OP22' or x == 'OP23' for x in allDischargesDataFrame['campaign']])
        failureCandidatesFilter = np.array([x == 'OP22' or x == 'OP23' for x in failureCandidatesDataFrame['campaign']])
        failureDischargesFilter = np.array([x == 'OP22' or x == 'OP23' for x in failureDischargesDataFrame['campaign']])
    else:
        raise ValueError(f'"{OP}" is not one of the allowed campaigns [OP22, OP23, OP223]')

    allDischargeFilterNan = np.array([not np.isnan(x) and not np.isnan(y) for x, y in zip(allDischargesDataFrame['ne'], allDischargesDataFrame['Te'])])
    failureCandidatesFilterNan = np.array([not np.isnan(x) and not np.isnan(y) for x, y in zip(failureCandidatesDataFrame['ne'], failureCandidatesDataFrame['Te'])])
    #failureDischargesFilterNan = np.array([np.isnan(x) or np.isnan(y) for x, y in zip(failureDischargesDataFrame['ne'], failureDischargesDataFrame['Te'])])
    
    allDischargeFilter = np.logical_and(allDischargeFilter, allDischargeFilterNan)
    failureCandidatesFilter = np.logical_and(failureCandidatesFilter, failureCandidatesFilterNan)

    configurationTotalPlungeNumber = []
    configurationFailVLimitPlungeNumber = []
    configurationFailVChangePlungeNumber = []
    configurationFailRLimitPlungeNumber = []
    configurationFailedPlungeNumber = []

    configurationTotalDischargeNumber = []
    configurationFailedDischargeNumber = []

    if X == 'configuration':
        configurationList = pd.read_csv('inputFiles/configurationListWithSettings.csv', sep=';')['configuration']
        for configuration in configurationList:
            configurationTotalPlungeNumber.append(sum([x == configuration for x in np.array(allDischargesDataFrame['configuration'])[allDischargeFilter]]))
            configurationFailVLimitPlungeNumber.append(sum([x == configuration and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['V_limit_failure_yn'])[failureCandidatesFilter])]))
            configurationFailRLimitPlungeNumber.append(sum([x == configuration and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['R_limit_failure_yn'])[failureCandidatesFilter])]))
            configurationFailVChangePlungeNumber.append(sum([x == configuration and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['V_limit_change'])[failureCandidatesFilter])]))
            configurationFailedPlungeNumber.append(sum([x == configuration and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['Failed'])[failureCandidatesFilter])]))
            
            dischargeFilter = np.array([x == configuration for x in np.array(allDischargesDataFrame['configuration'])[allDischargeFilter]])
            configurationTotalDischargeNumber.append(len(np.unique(np.array(allDischargesDataFrame['dischargeID'])[allDischargeFilter][dischargeFilter])))
            configurationFailedDischargeNumber.append(sum([x == configuration and y > 0 for x,y in zip(np.array(failureDischargesDataFrame['configuration'])[failureDischargesFilter], np.array(failureDischargesDataFrame['numberOfFailedPlunges'])[failureDischargesFilter])]))

    elif X == 'configurationShort':
        configurationList = pd.read_csv('inputFiles/configurationListWithSettings.csv', sep=';')['configuration']
        configurationList = list(map(lambda x: x[:3], configurationList))
        configurationList = np.unique(np.array(configurationList))
        
        for configuration in configurationList:	
            configurationTotalPlungeNumber.append(sum([x.startswith(configuration) for x in np.array(allDischargesDataFrame['configuration'])[allDischargeFilter]]))    
            configurationFailVLimitPlungeNumber.append(sum([x.startswith(configuration) and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['V_limit_failure_yn'])[failureCandidatesFilter])]))            
            configurationFailRLimitPlungeNumber.append(sum([x.startswith(configuration) and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['R_limit_failure_yn'])[failureCandidatesFilter])]))
            configurationFailVChangePlungeNumber.append(sum([x.startswith(configuration) and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['V_limit_change'])[failureCandidatesFilter])]))
            configurationFailedPlungeNumber.append(sum([x.startswith(configuration) and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['Failed'])[failureCandidatesFilter])]))

            dischargeFilter = np.array([x.startswith(configuration) for x in np.array(allDischargesDataFrame['configuration'])[allDischargeFilter]])
            configurationTotalDischargeNumber.append(len(np.unique(np.array(allDischargesDataFrame['dischargeID'])[allDischargeFilter][dischargeFilter])))
            configurationFailedDischargeNumber.append(sum([x.startswith(configuration) and y > 0 for x,y in zip(np.array(failureDischargesDataFrame['configuration'])[failureDischargesFilter], np.array(failureDischargesDataFrame['numberOfFailedPlunges'])[failureDischargesFilter])]))

    else:
        raise ValueError(f'The characteristic "{X}" is not implemented for analization')
    
    configurationFailVLimitPlungeNumber = np.array(configurationFailVLimitPlungeNumber)
    configurationTotalPlungeNumber = np.array(configurationTotalPlungeNumber)
    configurationFailRLimitPlungeNumber = np.array(configurationFailRLimitPlungeNumber)
    configurationFailVChangePlungeNumber = np.array(configurationFailVChangePlungeNumber)
    configurationFailedPlungeNumber = np.array(configurationFailedPlungeNumber)

    configurationFailedDischargeNumber = np.array(configurationFailedDischargeNumber)
    configurationTotalDischargeNumber = np.array(configurationTotalDischargeNumber)

    configFilter = np.array([x!=0 for x in configurationTotalPlungeNumber])    
    runs = configurationList[configFilter]
    bars = pd.DataFrame({'V_limit_failure': configurationFailVLimitPlungeNumber[configFilter]/configurationTotalPlungeNumber[configFilter],
                         'R_limit_failure': configurationFailRLimitPlungeNumber[configFilter]/configurationTotalPlungeNumber[configFilter],
                         'V_limit_change': configurationFailVChangePlungeNumber[configFilter]/configurationTotalPlungeNumber[configFilter],
                         'failed_plunges': configurationFailedPlungeNumber[configFilter]/configurationTotalPlungeNumber[configFilter]})
    
    bars2 = pd.DataFrame({'total_plunges': configurationTotalPlungeNumber[configFilter],
                          'failed_plunges': configurationFailedPlungeNumber[configFilter]
                          #'V_limit_failure': configurationFailVLimitPlungeNumber[configFilter],
                          #'R_limit_failure': configurationFailRLimitPlungeNumber[configFilter],
                          #'V_limit_change': configurationFailVChangePlungeNumber[configFilter]
                          })
          
    bars3 = pd.DataFrame({'total_discharges': configurationTotalDischargeNumber[configFilter],
                          'failed_discharges': configurationFailedDischargeNumber[configFilter]})
    
    bars4 = pd.DataFrame({'failed_discharges': configurationFailedDischargeNumber[configFilter]/configurationTotalDischargeNumber[configFilter]})
    
    x = np.arange(len(runs))  # the label locations
    width = 1/(len(bars.keys()) + 1)  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained', figsize = (10, 8))

    for combination, difference in bars.items():
        offset = width * multiplier
        ax.bar(x + offset, difference, width, label=combination)
        multiplier += 1

    ax.grid(True)
    ax.set_xticks(x + width, runs)
    ax.tick_params(axis='x', labelrotation = 90, labelsize = scalesize)
    ax.tick_params(axis='y', labelsize = scalesize)
    #ax.set_xlabel('configuration', fontsize = labelsize)
    ax.set_ylabel('relative number of plunges that are characterized as\nshortening based on the shown criteria', fontsize = labelsize)
    ax.legend(loc='upper left', prop={'size': legendsize})

    fig.savefig(f'results/LP_{LP}/compareFailuresDependingOn{X}.png', bbox_inches='tight')
    plt.show()
    plt.close()

          
    x = np.arange(len(runs))  # the label locations
    width = 1/(len(bars2.keys()) + 1)  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained', figsize = (10, 8))

    for combination, difference in bars2.items():
        offset = width * multiplier
        ax.bar(x + offset, difference, width, label=combination)
        multiplier += 1

    ax.grid(True)
    ax.set_xticks(x + width, runs)
    ax.tick_params(axis='x', labelrotation = 90, labelsize = scalesize)
    ax.tick_params(axis='y', labelsize = scalesize)
    #ax.set_xlabel('configuration', fontsize = labelsize)
    ax.set_yscale('log')
    ax.set_ylabel('absolute number of plunges that are characterized as\nshortening based on the shown criteria', fontsize = labelsize)
    ax.legend(loc='upper left', prop={'size': legendsize})

    fig.savefig(f'results/LP_{LP}/compareFailuresDependingOn{X}TotalNo.png', bbox_inches='tight')
    plt.show()
    plt.close()

    x = np.arange(len(runs))  # the label locations
    width = 1/(len(bars3.keys()) + 1)  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained', figsize = (10, 8))

    for combination, difference in bars3.items():
        offset = width * multiplier
        ax.bar(x + offset, difference, width, label=combination)
        multiplier += 1

    ax.grid(True)
    ax.set_xticks(x + width, runs)
    ax.tick_params(axis='x', labelrotation = 90, labelsize = scalesize)
    ax.tick_params(axis='y', labelsize = scalesize)
    ax.set_yscale('log')
    #ax.set_xlabel('configuration', fontsize = labelsize)
    ax.set_ylabel('absolute number of discharges that are characterized\nas shortening', fontsize = labelsize)
    ax.legend(loc='upper left', prop={'size': legendsize})

    fig.savefig(f'results/LP_{LP}/compareFailuresDependingOn{X}DischargesTotalNo.png', bbox_inches='tight')
    plt.show()
    plt.close()

    x = np.arange(len(runs))  # the label locations
    width = 1/(len(bars4.keys()) + 1)  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained', figsize = (10, 8))

    for combination, difference in bars4.items():
        offset = width * multiplier
        ax.bar(x + offset, difference, width, label=combination)
        multiplier += 1

    ax.grid(True)
    ax.set_xticks(x + width, runs)
    ax.tick_params(axis='x', labelrotation = 90, labelsize = scalesize)
    ax.tick_params(axis='y', labelsize = scalesize)
    #ax.set_xlabel('configuration', fontsize = labelsize)
    ax.set_ylabel('relative number of discharges that are characterized\nas shortening', fontsize = labelsize)
    ax.legend(loc='upper left', prop={'size': legendsize})

    fig.savefig(f'results/LP_{LP}/compareFailuresDependingOn{X}Discharges.png', bbox_inches='tight')
    plt.show()
    plt.close()
