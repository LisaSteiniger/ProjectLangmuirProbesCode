""" This file contains program routines to read data from input files, xDrive and dbArchive"""

import requests
import os
import itertools
import w7xarchive  

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import settings as settings
from src.dlp_data import extract_divertor_probe_data as extract
from src.dlp_data import fetch_data_dlp as fetch 
from scipy.signal import argrelextrema

#######################################################################################################################################################################        
def readLangmuirProbeOperationalParameters(dischargeID: str,
                                           LP: str,
                                           R_limit: int|float= 300,
                                           V_min_ideal: int|float= -180,
                                           V_max_ideal: int|float= 20,
                                           V_tolerance: float|int= 0.05,
                                           fetched: bool= True,
                                           plotting: bool= False,
                                           safe: str= 'results/LP_'
                                           ) -> list[list[str], list[float], list[int|float], list[bool], list[int|float]]:
    ''' This function reads operational parameters of a given Langmuir Probe
        -> this includes input voltage, probe potential (u_probe-u_bridge), and current
        Tests each plunge on shortening 
        -> by comparing range of probe potential to desired range of [-180, 20] V with some tolerance to it
        -> by comparing probe potential before and after insertion of the probe
        -> by calculating the probe line resistance and comparing it to the expected value
            -> R = V_input/I should not be equal/smaller than 5 Ohm
            
        Returns five lists of same length: 
        -> [dischargeID]*len(plunges)
        -> plunge times
        -> number of voltage extrema that didn't reach the required limit per plunge
        -> true/false for voltage extrema changes being in the allowed range or not for each plunge
        -> number of R_probeLine values that didn't reach the required limit per plunge
        
        "dischargeID" is the ID of the discharge e.g. "20241127.010"
        "LP" is the internal ID of the Langmuir Probe e.g. '51222'
        -> see settings to get full set of Langmuir Probe IDs
        "R_limit" is the limit for averaged probe line resistance in (Ohm)
        -> if the averaged resistance drops below, shortening is assumed
        "V_min_ideal" is the lower limit of the ideal voltage range in (V)
        "V_max_ideal" is the upper limit of the ideal voltage range in (V)
        "V_tolerance" is the tolerance on the ideal voltage range
        -> e.g. 0.05 tolerance means that voltages < -171 V and > 19 V are ok
        "fetched" is True - only the data from the start of a plunge + 100 ms is fetched
                  is False, the whole discharge data is fetched
        "plotting" decides if the raw data (probe voltage, probe current, and probe line resistance) are plotted
        "safe" is defalut path structure for saving a plot'''
    
    #different indicators of shortening -> each plunge gets one entry
    #if no shortening occurs, [0, False, 0] are added to the three lists
    V_limit_failure = []    #percentage of voltage extrema that lie outside the ideal voltage range with tolerance (e.g. minima above V_min_ideal*V_tolerance)
    V_limit_change = []     #do the extrema change significantly during a plunge? e.g. minima getting higher (=shortening occurs) or lower (=shortening disappears)
    R_limit_failure = []    #percentage of averaged R_probeline values that lie below the threshold of R_limit

    try:
        probe = fetch.Probes(probe = LP, shot = dischargeID, piece_fetching = fetched) 
    except ValueError:
        print(f'value error in fetch.Probes for {LP} in {dischargeID}')
        return [dischargeID], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan]
    except TypeError:
        print(f'type error in fetch.Probes for {LP} in {dischargeID}, most probably due to triggers')
        return [dischargeID], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan]

    #quantities: time, probe voltage and probe current
    #-> stored as probe.time, probe.probe_voltage and probe.probe_current
    try:
        probe.get_data()
    except ValueError:
        print(f'value error in get_data for {LP} in {dischargeID}')
        return [dischargeID], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan] 

    #quantity: input voltage is strored as probe.input_voltage
    probe.get_input_voltage()

    #quantity: reference voltage is strored as probe.reference_voltage
    probe.get_reference_voltage()

    #get the plunge times as probe.pop_up_time (set time instances for the diagnostic) 
    #-> has a 1 sec delay, real plunging in time is (probe.pop_up_time - 1) sec
    plunges = np.array(probe.pop_up_time) - 1
    
    #control if all arrays are of same size
    if len(probe.probe_voltage) == len(probe.probe_current)\
        and len(probe.probe_voltage)== len(probe.time)\
             and len(probe.probe_voltage) == len(probe.input_voltage):
        pass

    else:
        return [dischargeID], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan] 

    #compensation for insertion delay depends on divertor unit
    #lower divertor unit Langmuir Probes have 0 as 2. digit in their ID, upper DU probes have 1
    if LP[1] == '0':
        t_compensationDU = 0.01
    if LP[1] == '1':
        t_compensationDU = 0.015

    #determine local maximal and minimal voltages for each plunge
    for counter, plunge in enumerate(plunges):
        #filter the timeinterval for a plunge
        if counter == len(plunges) - 1:
            t_in = plunge
            t_out = probe.time[-1]
        else:
            t_in = plunge
            t_out = plunges[counter + 1]

        t_hold_start = t_in + probe.insertion_time + t_compensationDU
        t_hold_stop = t_in + probe.insertion_time + probe.hold_time + t_compensationDU
        timefilter = np.array([t_in < x < t_out for x in probe.time])
        timefilter_hold = np.array([t_hold_start < x < t_hold_stop for x in probe.time])
        voltage_timeHold_filter = np.array([(y < -0.5 or 0.5 < y) and t_hold_start < x < t_hold_stop for x, y in zip(probe.time, probe.input_voltage)])

        #find indices local extrema (also flat ones)
        V_max_indices = argrelextrema(probe.probe_voltage[timefilter_hold], np.greater_equal, order=20, mode='wrap')
        V_min_indices = argrelextrema(probe.probe_voltage[timefilter_hold], np.less_equal, order=20, mode='wrap')

        V_max_indices_plunge = argrelextrema(probe.probe_voltage[timefilter], np.greater_equal, order=20, mode='wrap')
        V_min_indices_plunge = argrelextrema(probe.probe_voltage[timefilter], np.less_equal, order=20, mode='wrap')

        #if no extrema are found:
        for orderGrade in range(19, 9):
            if len(V_max_indices[0]) == 0\
                or len(V_min_indices[0]) == 0\
                    or len(V_min_indices_plunge[0]) == 0\
                        or len(V_max_indices_plunge[0]) == 0:
                V_max_indices = argrelextrema(probe.probe_voltage[timefilter_hold], np.greater_equal, order=orderGrade, mode='wrap')
                V_min_indices = argrelextrema(probe.probe_voltage[timefilter_hold], np.less_equal, order=orderGrade, mode='wrap')

                V_max_indices_plunge = argrelextrema(probe.probe_voltage[timefilter], np.greater_equal, order=orderGrade, mode='wrap')
                V_min_indices_plunge = argrelextrema(probe.probe_voltage[timefilter], np.less_equal, order=orderGrade, mode='wrap')
            else:
                break
        
        if len(V_max_indices[0]) == 0\
            or len(V_min_indices[0]) == 0\
                or len(V_min_indices_plunge[0]) == 0\
                    or len(V_max_indices_plunge[0]) == 0:
            V_limit_failure.append('no extrema') 
            V_limit_change.append('no extrema') 
            R_limit_failure.append('no extrema') 
            continue

        #find local extrema values
        V_max = probe.probe_voltage[timefilter_hold][V_max_indices]
        V_min = probe.probe_voltage[timefilter_hold][V_min_indices]

        V_max_plunge = probe.probe_voltage[timefilter][V_max_indices_plunge]
        V_min_plunge = probe.probe_voltage[timefilter][V_min_indices_plunge]

        #test for too high minima and too high maxima
        V_max_tooLow = [x < V_max_ideal - (abs(V_max_ideal) + abs(V_min_ideal)) * V_tolerance for x in V_max]
        V_min_tooHigh = [x > V_min_ideal + (abs(V_max_ideal) + abs(V_min_ideal)) * V_tolerance for x in V_min]

        V_limit_failure.append((sum(V_min_tooHigh) + sum(V_max_tooLow))/(len(V_min_tooHigh) + len(V_max_tooLow)))
        if sum(V_max_tooLow) > 0 or sum(V_min_tooHigh) > 0:
            print(f'Shortening in {dischargeID} plunge {counter}?: {str((sum(V_min_tooHigh)+sum(V_max_tooLow))/(len(V_min_tooHigh) + len(V_max_tooLow)))} of voltage extrema outside of [{str(V_min_ideal)}, {str(V_max_ideal)}] V with tolerance {str(V_tolerance*100)}%')
        
        #test for changes of voltage minima/maxima
        V_max_min = min(V_max_plunge)
        V_max_max = max(V_max_plunge)
        
        V_min_min = min(V_min_plunge)
        V_min_max = max(V_min_plunge)

        V_limit_change.append(abs(V_min_max - V_min_min) > abs(V_max_max - V_min_min) * V_tolerance\
                              or abs(V_max_max - V_max_min) > abs(V_max_max - V_min_min) * V_tolerance)
        if abs(V_min_max - V_min_min) > abs(V_max_max - V_min_min) * V_tolerance\
           or abs(V_max_max - V_max_min) > abs(V_max_max - V_min_min) * V_tolerance:
            print(f'Shortening in {dischargeID} plunge {counter}?: voltage extrema range changes more than the allowed {str(abs(V_max_max - V_min_min)*V_tolerance)} V')

        #investigate probe line resistance
        if V_limit_failure[-1] > 0.03 or V_limit_change[-1]:
            R_probeLine = abs(probe.input_voltage[voltage_timeHold_filter]/probe.probe_current[voltage_timeHold_filter])
        else: 
            R_probeLine = abs(50/(probe.reference_voltage[voltage_timeHold_filter]/probe.probe_voltage[voltage_timeHold_filter] - 1))

        if len(R_probeLine) > 50:
            R_probeLine_av = list(itertools.chain.from_iterable([[np.mean(x)] * len(x) for x in np.array_split(R_probeLine, len(R_probeLine)/50)]))
        else:
            print(f'{dischargeID} plunge {counter}: only {str(len(R_probeLine))} R_probeLine values -> averaged to get one R_probeLine_av')
            R_probeLine_av = [np.mean(R_probeLine)]

        shortening_R = np.array([x < R_limit for x in np.array(R_probeLine_av)])

        R_limit_failure.append(sum(shortening_R)/len(shortening_R))
        if sum(shortening_R) > 0:
            print(f'Shortening in {dischargeID} plunge {counter}?: {str(sum(shortening_R/len(shortening_R)))} of R_probeLine values below {str(R_limit)} Ohm')
        
        #plotting
        if plotting:
            fig, [ax1, ax2, ax3] = plt.subplots(3, 1, figsize=(10, 15), layout='constrained', sharex=True)

            ax1.axhline(-180, c='b', linestyle=':')  
            ax1.axhline(20, c='r', linestyle=':')
            ax1.plot(probe.time[timefilter_hold][V_max_indices], probe.probe_voltage[timefilter_hold][V_max_indices], 'rx')
            ax1.plot(probe.time[timefilter_hold][V_min_indices], probe.probe_voltage[timefilter_hold][V_min_indices], 'bx')
            ax1.plot(probe.time[voltage_timeHold_filter][shortening_R], probe.probe_voltage[voltage_timeHold_filter][shortening_R], 'cx')
            ax1.plot(probe.time[timefilter], probe.probe_voltage[timefilter], 'k-')
            ax1.axvspan(t_hold_start, t_hold_stop, color='darkgrey')
            ax1.set_ylabel('probe voltage\nin (V)', fontsize=settings.labelsize)
            #ax1.tick_params(axis='x', labelsize = settings.scalesize)
            ax1.tick_params(axis='y', labelsize = settings.scalesize)

            ax2.plot(probe.time[voltage_timeHold_filter][shortening_R], probe.probe_current[voltage_timeHold_filter][shortening_R], 'cx')
            ax2.plot(probe.time[timefilter], probe.probe_current[timefilter], 'k-')
            ax2.axvspan(t_hold_start, t_hold_stop, color='darkgrey')
            ax2.set_ylabel('probe current\nin (A)', fontsize=settings.labelsize)
            #ax2.tick_params(axis='x', labelsize = settings.scalesize)
            ax2.tick_params(axis='y', labelsize = settings.scalesize)

            ax3.plot(probe.time[voltage_timeHold_filter][shortening_R], R_probeLine[shortening_R], 'cx')
            ax3.plot(probe.time[voltage_timeHold_filter], R_probeLine, 'k-')
            ax3.plot(probe.time[voltage_timeHold_filter], R_probeLine_av, 'r-')
            ax3.axvspan(t_hold_start, t_hold_stop, color='darkgrey')
            ax3.set_ylabel('probe line\nresistance\nin (Ohm)', fontsize=settings.labelsize)
            ax3.set_xlabel('time in (s) since beginning of discharge', fontsize=settings.labelsize)
            ax3.tick_params(axis='x', labelsize = settings.scalesize)
            ax3.tick_params(axis='y', labelsize = settings.scalesize)

            plt.figtext(0.5, 1.01, f'Langmuir Probe {LP}: {dischargeID} plunge {str(counter)}', horizontalalignment='center', fontsize=settings.labelsize)
        
            fig.savefig(f'{safe}{LP}/{dischargeID}_plunge{str(counter)}.png', bbox_inches='tight')
            plt.show()
            plt.close()
    
    if LP[1] == '0':
        LPs_lower = [LP]
        LPs_upper = []
        result_index = 0
    else:
        LPs_lower = []
        LPs_upper = [LP]
        result_index = 1
    
    results = readLangmuirProbeDataFromXdrive(dischargeID, LPs_lower, LPs_upper)
    
    if type(results) == str or (len(results[0+result_index][0]) == 1 and np.isnan(results[0+result_index][0][0])):
        ne = [np.nan]*len(plunges)
        sne = [np.nan]*len(plunges)
        Te = [np.nan]*len(plunges)
        sTe = [np.nan]*len(plunges)
    else:
        if len(results[0 + result_index][0]) == len(plunges):
            ne = results[0 + result_index][0] 
            Te = results[2 + result_index][0]
            sne = results[8 + result_index][0]
            sTe = results[10 + result_index][0]
        else:
            print(f'problem with number of plunges: measurement for {str(len(results[0 + result_index][0]))}, plunge times for {str(len(plunges))}')
            if len(results[0 + result_index][0]) < len(plunges):
                plunges = plunges[:len(results[0 + result_index][0])]
                V_limit_change = V_limit_change[:len(results[0 + result_index][0])]
                V_limit_failure = V_limit_failure[:len(results[0 + result_index][0])]
                R_limit_failure = R_limit_failure[:len(results[0 + result_index][0])]
                ne = results[0 + result_index][0] 
                Te = results[2 + result_index][0]
                sne = results[8 + result_index][0]
                sTe = results[10 + result_index][0]
            else:
                ne = results[0 + result_index][0][:len(plunges)]
                Te = results[2 + result_index][0][:len(plunges)]
                sne = results[8 + result_index][0][:len(plunges)]
                sTe = results[10 + result_index][0][:len(plunges)]

    return [dischargeID]*len(plunges), plunges, V_limit_failure, V_limit_change, R_limit_failure, ne, Te, sne, sTe

#######################################################################################################################################################################
def readLangmuirProbeDataFromXdrive(dischargeID: str, 
                                    LPs_lower: list[str],
                                    LPs_upper: list[str],
                                    xdrive_directory: str= "//x-drive/Diagnostic-logbooks/QRP-LangmuirProbes/QRP02-Divertor Langmuir Probes/Analysis/OP2/"
                                    ) -> str|list[list]:
    ''' This functions reads the electron density and electron temperature measured by the pop-up Langmuir Probes for a given discharge from xdrive
        Langmuir Probes are numbered as follows: TM2h07 holds probe 0 to 5, TM3h01 6 to 13, TM8h01 14 to 17
        -> internal order: with increasing distance from pumping gap 
        -> index_divertorUnits applies same scheme
        
        If data is available:
        -> Returns electron temperature Te in [eV], electron density ne in [1/m^3], corresponding measurement times in [s], and standard deviations
            -> for upper and lower divertor unit
            -> each of them is nested list of ndim=2 with each line representing measurements over time at one Langmuir Probe positions 
            -> (=each column represents measurements at all positions at one time)
        -> Additionally returns information about which Langmuir Probes were active 
            -> index_upper and index_lower provide the indices for the active probes
        Returns list of these lists [ne_lower, ne_upper, Te_lower, Te_upper, t_lower, t_upper, index_lower, index_upper, sne_lower, sne_upper, sTe_lower, sTe_upper]
        -> for the required probes

        If no data is available or wrong units are used: 
        Returns string
        
        "dischargeID" is the ID of the discharge in question and should be something like "20241127.010"
        "LPs_*" are the internal IDs of divertor Langmuir Probes on lower/upper divertor unit, e.g. "50209" (LDU -> 0 at 2. pos.) or "51222" (UDU -> 1 at 2. pos.)
        "xdrive_directory" is the path to the xdrive folder with Langmuir Probe data for all available discharges
        -> default is given for OP2 data'''
    
    ###########################################################
    #test if langmuir probe data is available for that discharge
    
    #add dischargeID to path (as a subfolder)
    pathLP = f"{xdrive_directory}/{dischargeID}/"

    if not os.path.exists(pathLP):
        print(f'No LP data available for {dischargeID}')
        return f'No LP data available for {dischargeID}'
    ###########################################################
    #setup of probe names and index lists for active Langmuir Probes

    #internal naming of probes in xdrive, saved in settings for each target module
    probes_lower = list(itertools.chain.from_iterable([settings.OP2_TM2_IDlower, settings.OP2_TM3_IDlower, settings.OP2_TM8_IDlower]))
    probes_upper = list(itertools.chain.from_iterable([settings.OP2_TM2_IDupper, settings.OP2_TM3_IDupper, settings.OP2_TM8_IDupper]))
    
    #convert internal naming to indices
    LP_indices_lower = []
    LP_indices_upper = []
    for LP in LPs_lower:
        if LP in probes_lower:
            LP_indices_lower.append(probes_lower.index(LP))
        else:
            raise ValueError(f'{LP} is no lower divertor Langmuir Probe')
    for LP in LPs_upper:
        if LP in probes_upper:
            LP_indices_upper.append(probes_upper.index(LP))
        else:
            raise ValueError(f'{LP} is no upper divertor Langmuir Probe')

    #index lists for each divertor unit in case that all LPs were active
    #indices of inactive probes will be removed later
    index_lower = list(range(18))
    index_upper = list(range(18))
    
    ###########################################################
    #read data and remove indices of inactive Langmuir Probes
    
    #read data from xdrive
    #data_DU stores ne, Te, t
    #fails_DU record missing LP data
    data_lower, data_upper, fails_lower, fails_upper = extract.fetch_xdrive_data(shot = dischargeID)
    
    #remove indices of inactive/missing LPs
    for fail in fails_lower:
        if not fail in LPs_lower:
            index_lower.remove(probes_lower.index(fail))
        else:
            index_lower[index_lower.index(probes_lower.index(fail))] = f'failed_{probes_lower.index(fail)}'
    for fail in fails_upper:
        if not fail in LPs_upper:
            index_upper.remove(probes_upper.index(fail))
        else:
            index_upper[index_upper.index(probes_upper.index(fail))] = f'failed_{int(probes_upper.index(fail))}'

    #needs at least one active Langmuir Probe to continue
    #test_index must contain that active Langmuir Probe as it will be used to test for correct units
    #test_data is the corresponing data set
    if len(index_lower) != 0 and sum([type(x)==int for x in index_lower]) != 0:
        failed = np.array([type(x)==int for x in index_lower])
        test_index = int(np.array(index_lower)[failed][0])
        test_data = data_lower
    elif len(index_upper) != 0 and sum([type(x)==int for x in index_upper]) != 0:
        failed = np.array([type(x)==int for x in index_upper])
        test_index = int(np.array(index_upper)[failed][0])
        test_data = data_upper
    else:
        return f'No LP data available for {dischargeID}'
    
    ###########################################################
    #test units to be (s), (1E+18 1/m^3) and (eV) 
    #adds data for all active Langmuir Probes if units are correct
    
    if test_data[test_index].units['time'] == 's' and test_data[test_index].units['ne'] == '10$^{18}$m$^{-3}$' and test_data[test_index].units['Te'] == 'eV':
        #lists will hold the values for all active Langmuir Probes
        ne_lower, ne_upper, Te_lower, Te_upper, t_lower, t_upper = [], [], [], [], [], []
        
        #lists will hold the standard deviations of ne and Te for all active Langmuir Probes
        sne_lower, sne_upper, sTe_lower, sTe_upper = [], [], [], []

        for counter, i in enumerate(index_lower):  
            if type(i)==str:
                ne_lower.append([np.nan])          
                sne_lower.append([np.nan])          
                Te_lower.append([np.nan])          
                sTe_lower.append([np.nan])          
                t_lower.append([np.nan])
                index_lower[counter] = int(i.split('_')[-1])
                continue

            #filter out measurements that are nonexisting (fake measurement value for t = 0 is inserted)
            filter_lower = np.array([j == 0 for j in data_lower[i].time])

            #values are given as ne (1e+18 1/m^3), convert to (1/m^3)
            ne_lower.append(list(np.array(data_lower[i].ne)[~filter_lower]*1e+18))  
            sne_lower.append(list(np.array(data_lower[i].sne)[~filter_lower]*1e+18))

            Te_lower.append(list(np.array(data_lower[i].Te)[~filter_lower]))
            sTe_lower.append(list(np.array(data_lower[i].sTe)[~filter_lower]))

            t_lower.append(list(np.array(data_lower[i].time)[~filter_lower]))
        
        for counter, i in enumerate(index_upper):
            if type(i)==str:
                ne_upper.append([np.nan])          
                sne_upper.append([np.nan])          
                Te_upper.append([np.nan])          
                sTe_upper.append([np.nan])          
                t_upper.append([np.nan])
                index_upper[counter] = int(i.split('_')[-1])
                continue
            
            #filter out measurements that are nonexisting (file exists but fake measurement value for t = 0 is inserted)
            filter_upper = np.array([j == 0 for j in data_upper[i].time])

            #values are given as ne (1e+18 1/m^3), convert to (1/m^3)
            ne_upper.append(list(np.array(data_upper[i].ne)[~filter_upper]*1e+18))  
            sne_upper.append(list(np.array(data_upper[i].sne)[~filter_upper]*1e+18))  

            Te_upper.append(list(np.array(data_upper[i].Te)[~filter_upper]))
            sTe_upper.append(list(np.array(data_upper[i].sTe)[~filter_upper]))
            
            t_upper.append(list(np.array(data_upper[i].time)[~filter_upper]))

        filter_indices_lower = []
        for LP_index_lower in LP_indices_lower:
            filter_indices_lower.append(index_lower.index(LP_index_lower))
        if len(filter_indices_lower) != 0:
            ne_lower = [ne_lower[i] for i in filter_indices_lower]
            sne_lower = [sne_lower[i] for i in filter_indices_lower]
            Te_lower = [Te_lower[i] for i in filter_indices_lower]
            sTe_lower = [sTe_lower[i] for i in filter_indices_lower]
            t_lower = [t_lower[i] for i in filter_indices_lower]
            #ne_lower = np.array(ne_lower)[np.array(filter_indices_lower)]
            #sne_lower = np.array(sne_lower)[np.array(filter_indices_lower)]
            #Te_lower = np.array(Te_lower)[np.array(filter_indices_lower)]
            #sTe_lower = np.array(sTe_lower)[np.array(filter_indices_lower)]
            #t_lower = np.array(t_lower)[np.array(filter_indices_lower)]
        else:
            ne_lower = np.array([[np.nan]])          
            sne_lower = np.array([[np.nan]])          
            Te_lower = np.array([[np.nan]])          
            sTe_lower = np.array([[np.nan]])          
            t_lower = np.array([[np.nan]])

        filter_indices_upper = []
        for LP_index_upper in LP_indices_upper:
            filter_indices_upper.append(index_upper.index(LP_index_upper))
        if len(filter_indices_upper) != 0:
            ne_upper = [ne_upper[i] for i in filter_indices_upper]
            sne_upper = [sne_upper[i] for i in filter_indices_upper]
            Te_upper = [Te_upper[i] for i in filter_indices_upper]
            sTe_upper = [sTe_upper[i] for i in filter_indices_upper]
            t_upper = [t_upper[i] for i in filter_indices_upper]
            #ne_upper = np.array(ne_upper)[np.array(filter_indices_upper)]
            #sne_upper = np.array(sne_upper)[np.array(filter_indices_upper)]
            #Te_upper = np.array(Te_upper)[np.array(filter_indices_upper)]
            #sTe_upper = np.array(sTe_upper)[np.array(filter_indices_upper)]
            #t_upper = np.array(t_upper)[np.array(filter_indices_upper)]
        else:
            ne_upper = np.array([[np.nan]])          
            sne_upper = np.array([[np.nan]])          
            Te_upper = np.array([[np.nan]])          
            sTe_upper = np.array([[np.nan]])          
            t_upper = np.array([[np.nan]])

        #careful, not all subarrays have the same shape
        #len(Te_upper[0]) == len(ne_upper[0]), but not neccessarily len(Te_upper[0]) == len(Te_upper[1])
        return [ne_lower, ne_upper, Te_lower, Te_upper, t_lower, t_upper, index_lower, index_upper, sne_lower, sne_upper, sTe_lower, sTe_upper]
    
    else:
        #for units not matching the required units
        return f'Wrong units for Langmuir Probe measurement in {dischargeID}'

#######################################################################################################################################################################        
def readAllShotNumbersFromLogbook(config: str, 
                                  filterSelected: str, 
                                  q_add: str,
                                  filesExist: bool =False, 
                                  safe: str ='inputFiles/configurations/dischargeList', 
                                  overviewTableLink: str ='results/calculationTables/results_',
                                  url: str= 'https://w7x-logbook.ipp-hgw.mpg.de/api/_search',
                                  urlTriggerBase: str= 'http://archive-webapi.ipp-hgw.mpg.de/programs.json?from='
                                  ) -> pd.DataFrame|str:
  
    ''' This function is responsible for finding all discharges in the Logbook according to the given filter for OP2.2 and OP2.3
        
        Returns string if no discharges are found, and pd.DataFrame if discharges are found
            -> dataFrame keys: "configuration", "dischargeID", "duration", "duration_planned", "durationHeating", "overviewTable"
                -> "duration" is time difference between trigger t1 and t4
                -> "duration_planned" is planned ECRH duration
                -> "durationHeating" is summed duration of ECRH and NBI
                -> "overviewTable" is default safe for calculation tables
            -> DataFrame is saved as .csv file under modified f"{safe}_{q_add}_{config}.csv"        
        
        "config" is magnetic field configuration such as "EIM000+2520" 
            -> applies internal configuration filter and determines path for saveing .csv file
        "filterSelected" is externally given filter (no conditioning, no gas valve tests, no sniffer tests,...)
            -> such as '!"Conditioning" AND !"gas valve tests" AND !"sniffer tests"'
        "q_add" is additional time filter for the year in string format
            -> 'OP22' means 2024, 'OP23' means 2025, and 'OP223' means 2024-2025
        "filesExist" determines, if possibly existing files with discharges for that configuration are overwritten or read out and returned
        "safe" is the basic structure of the path where the created .csv file is saved
        "overviewTableLink" is the basic structure of the path where latter calculation tables for a discharge are going to be saved
        "url" is where to read logbook data from for discharges  
        "urlTriggerBase" is the link to the archive for reading the internal triggers'''
    
    #if filesExist is True, existing files are read out, and only missing ones are created
    #if filesExist is False, existing files are overwritted, and missing ones are created
    if filesExist:  
        if os.path.isfile(f"{safe}_{q_add}_{config}.csv"):
            print(f'Confiuration discharge file for {config} in {q_add} exists and is read out')
            return pd.read_csv(f"{safe}_{q_add}_{config}.csv", sep=';', dtype={'dischargeID': str})

    ###########################################################
    #create complete filter

    #filter for configurations (config)
    q7 = f'tags.value:"{config}"'

    #combine neccessary filters   
    q = filterSelected + q7 

    if q_add == 'OP22':
        p = {'time':'[2024 TO 2024]', 'size':'9999', 'q' : q }   # OP2.2
    elif q_add == 'OP23':
        p = {'time':'[2025 TO 2025]', 'size':'9999', 'q' : q }   # OP2.3
    elif q_add == 'OP223':
        p = {'time':'[2024 TO 2025]', 'size':'9999', 'q' : q }   # OP2.2 and OP2.3

    ###########################################################
    #read discharge attributes (durations, frad, HHeRatio)

    dischargeIDs, duration_Trigger, duration_planned, duration_ECRH_NBI = [], [], [], []
    frad, HHeRatio, overviewTables = [], [], []

    #read logbook
    res = requests.get(url, params=p).json()

    if res['hits']['total']==0:
        print(f'no discharges found for {config} in {q_add}')
        return f'no discharges found for {config} in {q_add}'

    #transform dischargeID to format with three digits after '.' -> e.g. '20241127.010'
    for discharge in res['hits']['hits']:

        dischargeID = discharge['_id'].split('.')
        if len(dischargeID[1]) == 3:
            dischargeID = f'{dischargeID[0][3:]}.{dischargeID[1]}'
        elif len(dischargeID[1]) == 2:
            dischargeID = f'{dischargeID[0][3:]}.0{dischargeID[1]}'
        else:
            dischargeID = f'{dischargeID[0][3:]}.00{dischargeID[1]}'

        ###########################################################
        #read trigger times of t1 and t4 if they both exist (otherwise append np.nan to duration_Trigger)
        urlTrigger = urlTriggerBase + dischargeID + '#'
        resTrigger = requests.get(urlTrigger).json()

        t1, t4 = np.nan, np.nan
        trigger = np.nan
        if 'programs' in resTrigger.keys():
            if type(resTrigger['programs']) == list:
                if 'trigger' in resTrigger['programs'][0].keys():
                    if type(resTrigger['programs'][0]['trigger']) == dict:
                        if '1' in resTrigger['programs'][0]['trigger'].keys() and '4' in resTrigger['programs'][0]['trigger'].keys():
                            if type(resTrigger['programs'][0]['trigger']['4']) == list and type(resTrigger['programs'][0]['trigger']['1']) == list:
                                if len(resTrigger['programs'][0]['trigger']['4']) > 0 and len(resTrigger['programs'][0]['trigger']['1']) > 0:
                                    trigger = ((resTrigger['programs'][0]['trigger']['4'][0] - resTrigger['programs'][0]['trigger']['1'][0])/1e9)
                                    t4 = resTrigger['programs'][0]['trigger']['4'][0]
                                    t1 = resTrigger['programs'][0]['trigger']['1'][0]
                                else:
                                    print(f'{dischargeID}: Trigger "1" or "4" is an empty list')
                            else:
                                print(f'{dischargeID}: Trigger "1" or "4" is not a list')
                        else:
                            print(f'{dischargeID}: Trigger "1" or "4" does not exist')
                    else:
                        print(f'{dischargeID}: Trigger is not a dictionary')
                else:
                    print(f'{dischargeID}: No trigger in "programs"')
            else:
                print(f'{dischargeID}: "programs" is no list')
        else:
            print(f'{dischargeID}: Key "programs" does not exist')

        ###########################################################
        #try to get duration of the discharge based on heating power with threshhold 0.1 MW
        #filter out discharges with total heating power below 1 MW or electron density below 1E+19 1/m^3 for whole discharge
        #also filter out discharges without datastream for electron density

        if not np.isnan(t1) and not np.isnan(t4):
            heating_data = w7xarchive.get_signal({'ECRH': "Test/raw/W7X/CBG_ECRH/TotalPower_DATASTREAM/V1/0/Ptot_ECRH",
                                                  'NBIS3': "ArchiveDB/codac/W7X/ControlStation.2176/BE000_DATASTREAM/4/s3_Pel/scaled", 
                                                  'NBIS4': "ArchiveDB/codac/W7X/ControlStation.2176/BE000_DATASTREAM/5/s4_Pel/scaled", 
                                                  'NBIS7': "ArchiveDB/codac/W7X/ControlStation.2092/BE000_DATASTREAM/4/s7_Pel/scaled", 
                                                  'NBIS8': "ArchiveDB/codac/W7X/ControlStation.2092/BE000_DATASTREAM/5/s8_Pel/scaled", 
                                                  'ne': "ArchiveDB/raw/W7X/ControlStation.2185/Density_DATASTREAM/0/Density",
                                                  'Prad': "ArchiveDB/raw/W7X/ControlStation.2224/QRB_Divertor_Bolometry_Prad_DATASTREAM/15/Prad QRB weigted sum",
                                                  'HHeRatio': "ArchiveDB/raw/W7XAnalysis/QSK06_PassiveSpectroscopy/FastSpec_USB_HR4000_HHeRatio_DATASTREAM/V1/0/HHeRatio"},
                                                  t1, t4)
            
            #merge heating DataFrames on rounded time axis
            for i, key in enumerate(heating_data.keys()):
                if i == 0:
                    merged_heating = pd.DataFrame({'time': np.round(np.array(heating_data[key][0]), 3), key: heating_data[key][1]})
                else:
                    merged_heating = pd.merge(merged_heating, pd.DataFrame({'time': np.round(np.array(heating_data[key][0]), 3), key: heating_data[key][1]}), 'outer', on='time')
            
            #merged is only created if any data stream is existant
            if len(heating_data.keys()) > 1:
                heating = np.zeros_like(np.array(merged_heating['time']))
                for key in merged_heating.keys():
                    if 'ECRH' in key:
                        heating = heating + np.nan_to_num(np.array(merged_heating[key]))/1000
                    elif 'NBI' in key:
                        heating = heating + np.nan_to_num(np.array(merged_heating[key]))
                merged_heating['heating'] = heating

                filter_heating_discharge = [x > 1 for x in merged_heating['heating']]
                
                if 'ne' in heating_data.keys():
                    filter_ne_discharge = [x > 1e19 for x in merged_heating['ne']]
                    if sum(filter_heating_discharge) == 0 or sum(filter_ne_discharge) == 0:
                        continue
                else:
                    print(f'Is {dischargeID} a real discharge?')
                    continue
                
                filter_heating = [x > 0.1 for x in merged_heating['heating']]
                if sum(filter_heating) > 0:
                    duration_ECRH_NBI.append((merged_heating['time'][len(filter_heating) - 1 - filter_heating[::-1].index(True)] - merged_heating['time'][filter_heating.index(True)])/1e9)
                else:
                    duration_ECRH_NBI.append(np.nan)
        
                ###########################################################
                #calculate frad and read HHe ratio
                if 'Prad' in merged_heating.keys():
                    heating_0 = np.array([j != 0 for j in heating])
                    if sum(heating_0) == 0:
                        frad.append(np.nan)
                    else:
                        x = np.nansum(np.array(merged_heating['Prad'])[heating_0]/heating[heating_0])/len(heating[heating_0])
                        frad.append(x)
                else:
                    frad.append(np.nan)

                if 'HHeRatio' in merged_heating.keys():
                    x = np.nansum(np.array(merged_heating['HHeRatio']))/len(heating)
                    if x == 0:
                        HHeRatio.append(np.nan)
                    else:
                        HHeRatio.append(x)
                else:
                    HHeRatio.append(np.nan)

                ###########################################################

            #if no data stream existed including ne -> discard that discharge
            else:
                print(f'Is {dischargeID} a real discharge?')
                continue

        #if triggers do not exist
        else:
            HHeRatio.append(np.nan)    
            frad.append(np.nan)    
            duration_ECRH_NBI.append(np.nan)

        duration_Trigger.append(trigger)
        dischargeIDs.append(dischargeID) 
        overviewTables.append(f'{overviewTableLink}{dischargeID}.csv')

        print(dischargeID)

        ###########################################################
        #try to get planned ECRH heating durtion (planned != real (program abort) and no inclusion of NBI heating nor ICRH heating)
        for tag in discharge['_source']['tags']: 
            if 'catalog_id' in tag.keys():
                if tag['catalog_id'] == '1#3':
                    duration_planned.append(tag['ECRH duration'])
                    continue
        #if no ECRH duration tag was found
        if len(duration_planned) != len(dischargeIDs):
            duration_planned.append(np.nan)

        if abs(duration_planned[-1] - trigger) > 20:
            print(f'Planned and trigger duration do not match: control {dischargeID}')
    
        ###########################################################
    
    #create and save dataFrame with all discharges and their characteristics 
    dischargeTable = pd.DataFrame({'configuration': [config] * len(dischargeIDs),
                                   'dischargeID': dischargeIDs,
                                   'duration': duration_Trigger,
                                   'duration_planned': duration_planned,
                                   'durationHeating': duration_ECRH_NBI, 
                                   'frad': frad,
                                   'HHeRatio': HHeRatio,
                                   'overviewTable': overviewTables})
    
    dischargeTable.to_csv(f"{safe}_{q_add}_{config}.csv", sep=';')
    
    return dischargeTable