''' This file contains plotting fuctions for visualization of possilbe systematics in shortening'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plotFailuresInDependanceOfX(X: str,
                                LP: str,
                                OP: str,
                                pathAllDischarges: str ='',
                                pathFailureCandidates: str =''
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
        "pathFailureCandidates" is the path to the overview file with the Langmuir Probe failure indicators for failure candidates'''
    
    if pathAllDischarges == '':
        pathAllDischarges = f'results/LP_{LP}/{LP}_dischargePlungeList_FailureIndicators.csv'
    
    if pathFailureCandidates == '':
        pathFailureCandidates = f'results/LP_{LP}/{LP}failed_dischargePlungeList_FailureIndicators.csv'

    allDischargesDataFrame = pd.read_csv(pathAllDischarges, sep=';')
    failureCandidatesDataFrame = pd.read_csv(pathFailureCandidates, sep=';', dtype={'V_limit_failure_yn': bool, 'R_limit_failure_yn': bool, 'V_limit_change': bool})

    if OP == 'OP22':
        allDischargeFilter = np.array([x == 'OP22' for x in allDischargesDataFrame['campaign']])
        failureCandidatesFilter = np.array([x == 'OP22' for x in failureCandidatesDataFrame['campaign']])
    elif OP == 'OP23':
        allDischargeFilter = np.array([x == 'OP23' for x in allDischargesDataFrame['campaign']])
        failureCandidatesFilter = np.array([x == 'OP23' for x in failureCandidatesDataFrame['campaign']])
    elif OP == 'OP223':
        allDischargeFilter = np.array([x == 'OP22' or x == 'OP23' for x in allDischargesDataFrame['campaign']])
        failureCandidatesFilter = np.array([x == 'OP22' or x == 'OP23' for x in failureCandidatesDataFrame['campaign']])
    else:
        raise ValueError(f'"{OP}" is not one of the allowed campaigns [OP22, OP23, OP223]')
    
    configurationTotalPlungeNumber = []
    configurationFailVLimitPlungeNumber = []
    configurationFailVChangePlungeNumber = []
    configurationFailRLimitPlungeNumber = []

    if X == 'configuration':
        configurationList = pd.read_csv('inputFiles/configurationListWithSettings.csv', sep=';')['configuration']
        for configuration in configurationList:
            configurationTotalPlungeNumber.append(sum([x == configuration for x in np.array(allDischargesDataFrame['configuration'])[allDischargeFilter]]))
            configurationFailVLimitPlungeNumber.append(sum([x == configuration and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['V_limit_failure_yn'])[failureCandidatesFilter])]))
            configurationFailRLimitPlungeNumber.append(sum([x == configuration and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['R_limit_failure_yn'])[failureCandidatesFilter])]))
            configurationFailVChangePlungeNumber.append(sum([x == configuration and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['V_limit_change'])[failureCandidatesFilter])]))

    elif X == 'configurationShort':
        configurationList = pd.read_csv('inputFiles/configurations/..', sep=';')['configuration']
        configurationList = list(map(x[:3], configurationList))
        configurationList = np.unique(np.array(configurationList))
        
        for configuration in configurationList:	
            configurationTotalPlungeNumber.append(sum([x.startswith('configuration') for x in np.array(allDischargesDataFrame['configuration'])[allDischargeFilter]]))    
            configurationFailVLimitPlungeNumber.append(sum([x.startswith(configuration) and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['V_limit_failure_yn'])[failureCandidatesFilter])]))            
            configurationFailRLimitPlungeNumber.append(sum([x.startswith(configuration) and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['R_limit_failure_yn'])[failureCandidatesFilter])]))
            configurationFailVChangePlungeNumber.append(sum([x.startswith(configuration) and y for x,y in zip(np.array(failureCandidatesDataFrame['configuration'])[failureCandidatesFilter], np.array(failureCandidatesDataFrame['V_limit_change'])[failureCandidatesFilter])]))

    else:
        raise ValueError(f'The characteristic "{X}" is not implemented for analization')
    
    configurationFailVLimitPlungeNumber = np.array(configurationFailVLimitPlungeNumber)
    configurationTotalPlungeNumber = np.array(configurationTotalPlungeNumber)
    configurationFailRLimitPlungeNumber = np.array(configurationFailRLimitPlungeNumber)
    configurationFailVChangePlungeNumber = np.array(configurationFailVChangePlungeNumber)

    configFilter = np.array([x!=0 for x in configurationTotalPlungeNumber])    
    runs = configurationList[configFilter]
    bars = pd.DataFrame({'V_limit_failure': configurationFailVLimitPlungeNumber[configFilter]/configurationTotalPlungeNumber[configFilter],
                         'R_limit_failure': configurationFailRLimitPlungeNumber[configFilter]/configurationTotalPlungeNumber[configFilter],
                         'V_limit_change': configurationFailVChangePlungeNumber[configFilter]/configurationTotalPlungeNumber[configFilter]})
        
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
    ax.tick_params(axis='x', labelrotation = 90)#, labelsize = scalesize)
    ax.tick_params(axis='y')#, labelsize = scalesize)
    #ax.set_xlabel('configuration', fontsize = labelsize)
    ax.set_ylabel('relative number of plunges that are characterized as shortening based on the shown criteria')#, fontsize = labelsize)
    ax.legend(loc='upper left')#, prop={'size': legendsize})

    fig.savefig(f'results/LP_{LP}/compareFailuresDependingOn{X}.png', bbox_inches='tight')
    plt.show()
    plt.close()