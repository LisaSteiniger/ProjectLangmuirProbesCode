""" This file contains general settings for the whole project"""

#######################################################################################################################################################################
#INPUT THAT MIGHT BE CHANGED

#for plotting
legendsize = 15
scalesize = 15
labelsize = 15

#######################################################################################################################################################################
#INPUT THAT SHOULD NOT BE CHANGED

#IDs of Langmuir Probes 
#naming: campaign_targetModule_divertorUnit
#order: with increasing distance from pumping gap
OP2_TM2_IDlower = ['50201', '50203', '50205', '50207', '50209', '50211']
OP2_TM3_IDlower = ['50218', '50220', '50222', '50224', '50226', '50228', '50230', '50232']
OP2_TM8_IDlower = ['50246', '50248', '50249', '50251']

OP2_TM2_IDupper = ['51201', '51203', '51205', '51207', '51209', '51211']
OP2_TM3_IDupper = ['51218', '51220', '51222', '51224', '51226', '51228', '51230', '51232']
OP2_TM8_IDupper = ['51246', '51248', '51249', '51251']

#distances from pumping gap in (m) for Langmuir probes
#equal for both divertor units    
#naming: campaign_targetModule
#order: with increasing distance from pumping gap
OP2_TM2Distances = [0.10602897, 0.13182066, 0.15761507, 0.18341105, 0.20920802, 0.23500566] #Are those the right six values as TM2h07 is said to have Langmuir Probes closer to pumping gap?
OP2_TM3Distances = [0.32530727, 0.35110524, 0.37690349, 0.40270196, 0.42850063, 0.45429944, 0.48009839, 0.50589745] #Are those the right eight values or are they messed up with TM2h07?
OP2_TM8Distances = [0.09177, 0.11728, 0.13278, 0.15828] 

#filter options to choose from for the discharges (configuration filter is applied later):
#!!! HELIUM AND HYDROGEN DISCHARGES ARE NOT DISTINGUISHED !!!
q1 = '!"Conditioning" AND '
q2 = '!"gas valve tests" AND '
q3 = '!"sniffer tests" AND '
q4 = '!"reference discharge" AND '
q41 = '"ReferenceProgram"'   # in OP2.1 does not work
q5 = '!"Open valves" AND '
q6 = 'id:XP_* AND tags.value:"ok" AND '
q66 = 'id:XP_* AND '
q71 = 'tags.value:"Reference"'  # for OP2.1
q44 = '"reference discharge"'   # for OP2.2, OP2.3
q45 = '"Reference discharge"'   # for OP1.2b
qNBI = 'tags.value:"NBI source 7" OR tags.value:"NBI source 8" AND'

