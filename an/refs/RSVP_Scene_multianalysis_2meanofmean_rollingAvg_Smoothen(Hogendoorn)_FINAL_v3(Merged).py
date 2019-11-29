# This Python Source File Available in python_source/data_visualization/sample_trace_plot.py
import matplotlib.pyplot as plt
import numpy as np

import pylab
import os
from scipy import stats

import random
import copy
from scipy.optimize import curve_fit, differential_evolution
import warnings
import math

import sys
sys.path.append("D:\Program Files\PsychoPy2\Lib\site-packages")
sys.path.append("D:\Program Files\PsychoPy2\Lib\site-packages\psychopy")

from psychopy.iohub.datastore.util import ExperimentDataAccessUtility
from psychopy.iohub import EventConstants

'''
v2 Patch Notes: Re-scale bins to make sure trials are more eventually distributed across.

v3 Patch Notes: Standardised windows.

v4 Patch Notes: With error bars, consolidated via mean error across participants.

v5 Patch Notes: Threshold Analysis added.

'''

NrOfIterations = 10000
HogBin = 60
NrOfPS = 10 #Excluding ignored subjs
IgnoreListG = [] #Participants to exclude from analysis for goggle
IgnoreListN = [] #Participants to exclude from analysis for natural blink
averageRange = 120 #in ms

HogenAnalysis = True; NonDamp = False; Gaussian = False

GroupAnalysis = True; IndAnalysis = False
#GroupAnalysis = False; IndAnalysis = True #Unusable in this analysis.
                                          
#OnsetAnalysis = True; OffsetAnalysis = False
OnsetAnalysis = False; OffsetAnalysis = True

partset = [] 
filenameGogHDF = []
filenameNatHDF = []
    
if IndAnalysis == True:
    filenameNatHDF =  ['eventsNat11_C.hdf5','eventsNatNB11_C.hdf5']
    filenameGogHDF =  ['eventsGog11_C.hdf5','eventsGogNB11_C.hdf5'] #'eventsGogNB07_C_Nat.hdf5', 'eventsGogNB10_C_E.hdf5'
  
#Counters to extract data
GcounterNB = 0; GcounterB0 = 0; GcounterB1 = 0; GcounterB2 = 0; GcounterB3 = 0; GcounterB4 = 0; GcounterB5 = 0
GcounterB6 = 0; GcounterB7 = 0; GcounterB8 = 0; GcounterB9 = 0; GcounterB10 = 0; GcounterB11 = 0; GcounterB12 = 0

counterNB = 0; counterB0 = 0; counterB1 = 0; counterB2 = 0; counterB3 = 0; counterB4 = 0; counterB5 = 0
counterB6 = 0; counterB7 = 0; counterB8 = 0; counterB9 = 0; counterB10 = 0; counterB11 = 0; counterB12 = 0

TcounterNB = 0; TcounterB0 = 0; TcounterB1 = 0; TcounterB2 = 0; TcounterB3 = 0; TcounterB4 = 0; TcounterB5 = 0
TcounterB6 = 0; TcounterB7 = 0; TcounterB8 = 0; TcounterB9 = 0; TcounterB10 = 0; TcounterB11 = 0; TcounterB12 = 0

NBTotal = 0; B0Total = 0; B1Total = 0; B2Total = 0; B3Total = 0; B4Total = 0; B5Total = 0; B6Total = 0
B7Total = 0; B8Total = 0; B9Total = 0; B10Total = 0; B11Total = 0; B12Total = 0

TNBTotal = 0; TB0Total = 0; TB1Total = 0; TB2Total = 0; TB3Total = 0; TB4Total = 0; TB5Total = 0; TB6Total = 0
TB7Total = 0; TB8Total = 0; TB9Total = 0; TB10Total = 0; TB11Total = 0; TB12Total = 0

NBList = []; B0List = []; B1List = []; B2List = []; B3List = []; B4List = []; B5List = []; B6List = []
B7List = []; B8List = []; B9List = []; B10List = []; B11List = []; B12List = []

stderr = []

# Pupil size counter
IndPupilSizeList = []
TrialPupilSizeList = []
TBNormalizedList = []
PupilSizeListGog = []

thisfile = 0
part = 0
Sum2 = 0

ThresholdCounter = 0
GoggleD0Count = 0
Avgcounter = 0
Avgcounter2 = 1
maxLength = 0

PSdata1 = []
PSdata2 = []
averagedList = []
ListOfRaw = []
currentSet1 = []
currentSet2 = []
IndData = []

NextBin = HogBin
TotalTrialsInBin = 0
CorrectTrialsInBin = 0
Bin60msList = []

ThresholdCounter = 0
NoiseCount = 0
RejectCount = 0
TypoCount = 0

def add_np(A, B):
    m = max(len(A), len(B))
    return np.pad(A, (0, m-len(A)), 'constant') + np.pad(B, (0, m-len(B)), 'constant')

#Paste filename to be analysised here. (Hdf)
if GroupAnalysis == True:
    os.chdir('D:\Dropbox\RSVP\Scene Goggle\hdf5')
    readfile = os.listdir('Cleaned')
    
    start1 = 'eventsSGog0'
    start2 = 'eventsSGogNB0'
    num = 1
    ignorecounter = 0
    for readparticipantfile in range(NrOfPS+len(IgnoreListG)):
        if len(IgnoreListG) != 0:
            if num == IgnoreListG[ignorecounter]:
                num += 1
                ignorecounter += 1
                if ignorecounter > len(IgnoreListG)-1:
                    ignorecounter = len(IgnoreListG)-1
                continue
        
        if num == 10:
            start1 = 'eventsSGog'
            start2 = 'eventsSGogNB'
        
        participant = [str('%s%d' % (start1,num)),str('%s%d' % (start2,num))]
        for file in readfile:
            #Append filename to empty list if it ends with .csv
            if file.startswith(participant[0]):
                partset.append(file)
            if file.startswith(participant[1]):
                partset.append(file)
        filenameGogHDF.append(partset)
        num += 1
        partset = []
    divisions = len(filenameGogHDF)   

if GroupAnalysis == True:
    multiples = 2
    divisions = len(filenameGogHDF)
if IndAnalysis == True:
    multiples = 1
    divisions = 1

print filenameGogHDF, len(filenameGogHDF)

for analyse in range(len(filenameGogHDF*multiples)):
    #print part
    #print filename[thisfile][part]
    index1 = 0
    index2 = 0
    IndBlinkList = []
    # Load an ioDataStore file containing 500 Hz sample data from a
    # remote eye tracker that was recording both eyes. In the plotting example
    if GroupAnalysis == True:
        dataAccessUtil=ExperimentDataAccessUtility('D:\Dropbox\RSVP\Scene Goggle\hdf5\Cleaned', filenameGogHDF[thisfile][part], experimentCode=None,sessionCodes=[])

    if IndAnalysis == True:
        dataAccessUtil=ExperimentDataAccessUtility('D:\Dropbox\RSVP\Scene Goggle\hdf5\Cleaned', filenameGogHDF[part], experimentCode=None,sessionCodes=[])  

    ##### STEP A. #####
    # Retrieve a subset of the BINOCULAR_EYE_SAMPLE event attributes, for events that occurred
    # between each time period defined by the TRIAL_START and TRIAL_END trial variables of each entry
    # in the trial_conditions data table.
    #
    exp_data = dataAccessUtil.getConditionVariables(filter=None)

    # No need to keep the hdf5 file open anymore...
    #
    dataAccessUtil.close()

    for trial_index,trial_samples in enumerate(exp_data):
        index1 += 1
        
        if trial_samples.SESSION_ID == 0: 
            RejectCount += 1
            continue
        ## Remove trials where blink duration is abnoramlly low, most likely noise which triggered blink condition.
        if trial_samples.blinkend - trial_samples.blinkstart <= 0.011 and trial_samples.blinkend - trial_samples.blinkstart != 0.0:
            #print trial_samples.blinkstart, trial_samples.blinkend
            NoiseCount += 1
            continue
        if trial_samples.response2 < 0 or trial_samples.response2 > 4:
            TypoCount += 1
            #print filenameGogHDF[thisfile][part]
            continue
        
        if trial_samples.blinked == 1 and trial_samples.delay == 0:
            GoggleD0Count += 1
            continue
        #elif trial_samples.blinked == 1 and trial_samples.target_display - trial_samples.blinkstart <= 0.20:  
        elif trial_samples.blinked == 1 and trial_samples.target_display - trial_samples.blinkstart <= 0.20:                  
            B1Total += 1.
            if trial_samples.response2 == trial_samples.target_category:
                counterB1 += 1.            

        if trial_samples.blinked == 0 and trial_samples.target_category != -1:
            NBTotal += 1.
            if trial_samples.response2 == trial_samples.target_category:
                counterNB += 1.      

        elif trial_samples.blinked == 1 and trial_samples.blinkstart != 0. and trial_samples.blinkend != 0:
            PStargonset = trial_samples.target_display*1000 - trial_samples.blinkend*1000

            if PStargonset != 0:
                if trial_samples.response2 == trial_samples.target_category:
                    PSresp = 1.
                elif trial_samples.response2 != trial_samples.target_category:
                    PSresp = 0.                    
                PSdata1.append(PStargonset)    
                PSdata2.append(PSresp)  
                currentSet1.append(PStargonset)
                currentSet2.append(PSresp)
            
    if part == 1:
        NBList.append(counterNB/NBTotal)
        counterNB = 0; NBTotal = 0
        part = -1            
        thisfile += 1
        listforsorting = []

        ## Sort individual data and have it appended for error analysis just right after loop function.
        for sort in range(len(currentSet1)):
            objectappend = (currentSet1[sort],currentSet2[sort])
            listforsorting.append(objectappend)
        
        sortedlist = sorted(listforsorting, key=lambda x: x[0])
 
        NextBin = HogBin
        
        ## Original data Rolling Average
        for axisrange in range(int(np.max(currentSet1))):
            if axisrange == NextBin or axisrange == 11: #11 because first bin (D1) starts at 11ms.
                tail1 = axisrange-(averageRange/2)
                tail2 = axisrange+(averageRange/2)
            
                for consolidate in range(len(sortedlist)):
                    if sortedlist[consolidate][0] >= tail1 and sortedlist[consolidate][0] < tail2:
                        ListOfRaw.append(sortedlist[consolidate][1]) 
                        
                if axisrange != 11:
                    NextBin += HogBin
                    
                averagedList.append(np.mean(ListOfRaw))
                ListOfRaw = []     
        
        IndData.append(averagedList)
        ## Identify the longest list out of all participants
        thisMax = len(averagedList)
        if thisMax > maxLength:
            maxLength = thisMax
        currentSet1 = []
        currentSet2 = []
        averagedList = []

    part += 1
    
print 'GOGGLE'   

ValList = []
ErrList = []

## Generate error value for each x point based on combined data of all participant at that x point.
for GenError in range(maxLength):
    for participant in range(len(IndData)):
        try:
            if IndData[participant][GenError] == None or IndData[participant][GenError] == '' or math.isnan(IndData[participant][GenError]):
                continue
            ValList.append(IndData[participant][GenError])
        except:
            continue
    currentErr = stats.sem(ValList)
    ErrList.append(currentErr)

    ValList = []

listforsorting = []

## Sort group data and have it appended for error analysis.
for sort in range(len(PSdata1)):
    objectappend = (PSdata1[sort],PSdata2[sort])
    listforsorting.append(objectappend)

sortedlist = sorted(listforsorting, key=lambda x: x[0])

NextBin = HogBin

## Original data Rolling Average
for axisrange in range(int(np.max(PSdata1))):
    if axisrange == NextBin or axisrange == 11: #11 because first bin (D1) starts at 11ms.
        tail1 = axisrange-(averageRange/2)
        tail2 = axisrange+(averageRange/2)
    
        for consolidate in range(len(sortedlist)):
            if sortedlist[consolidate][0] >= tail1 and sortedlist[consolidate][0] < tail2:
                ListOfRaw.append(sortedlist[consolidate][1]) 
                
        if axisrange != 11:
            NextBin += HogBin
            
        averagedList.append(np.mean(ListOfRaw))
        ListOfRaw = []

NextBin = HogBin
fftGrpList = []
replacementGrpList = []

## Do iterations shuffled data before averaging it out at the end
for iterate in range(NrOfIterations): 
    
    listforsorting = []
    DataList = []
    RepDataList = []
    PSdata2_Scramble = copy.deepcopy(PSdata2)

    ## Merge original data's correct/wrong to the original timings
    for sort in range(len(PSdata1)):
        objectappend = (PSdata1[sort],PSdata2_Scramble[sort])
        listforsorting.append(objectappend)
    
    sortedlist = sorted(listforsorting, key=lambda x: x[0])    
    replacelist = []
    
    ## Replacement of original data
    for replace in range(len(sortedlist)):
        replacelist.append(sortedlist[random.randint(0,len(sortedlist)-1)])
     
    ## NOTE ONE IS replacelist and the other is replaceDlist     
    replacedlist = sorted(replacelist, key=lambda x: x[0])       
    
    NextBin = HogBin
    
    ## Original data replacement
    for axisrange in range(int(np.max(PSdata1))):
        if axisrange == NextBin or axisrange == 11: #11 because first bin (D1) starts at 11ms.
            tail1 = axisrange-(averageRange/2)
            tail2 = axisrange+(averageRange/2)
        
            for consolidate in range(len(sortedlist)):
                if replacedlist[consolidate][0] >= tail1 and replacedlist[consolidate][0] < tail2:
                    ListOfRaw.append(replacedlist[consolidate][1]) 
                    
            if axisrange != 11:
                NextBin += HogBin
                
            RepDataList.append(np.mean(ListOfRaw))
            ListOfRaw = []  
            
            
    averagedList2 = []           
    AvgListMean = np.mean(RepDataList)
    for div in RepDataList:
        div = div-AvgListMean
        averagedList2.append(div)
    
    ## Refer to https://plot.ly/matplotlib/fft/ ##
    ## FFT of iterations of original replacement data ##
    Fs = 1./(HogBin/1000.) #sampling rate; 1 divided by the sampling rate of 60ms used in the experiment.
    n = len(averagedList2) # length of the signal
    k = np.arange(n)
    T = np.float(n/Fs)
    y = averagedList2
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range
    
    Y = np.fft.fft(y)/n
    Y = Y[range(n/2)]
    
    replacementGrpList.append(abs(Y))
            
    ## ------------------------- Confidence Interval ------------------------------
    ## Randomise correct and wrong
    random.shuffle(PSdata2_Scramble)
    listforsorting = []
    
    ## Merge randomised correct/wrong to the original timings
    for sort in range(len(PSdata1)):
        objectappend = (PSdata1[sort],PSdata2_Scramble[sort])
        listforsorting.append(objectappend)
    
    sortedlist = sorted(listforsorting, key=lambda x: x[0])
    
    NextBin = HogBin
    
    ## Shuffled data Rolling Average
    for axisrange in range(int(np.max(PSdata1))):
        if axisrange == NextBin or axisrange == 11: #11 because first bin (D1) starts at 11ms.
            tail1 = axisrange-(averageRange/2)
            tail2 = axisrange+(averageRange/2)
        
            for consolidate in range(len(sortedlist)):
                if sortedlist[consolidate][0] >= tail1 and sortedlist[consolidate][0] < tail2:
                    ListOfRaw.append(sortedlist[consolidate][1]) 
                    
            if axisrange != 11:
                NextBin += HogBin
                
            DataList.append(np.mean(ListOfRaw))
            ListOfRaw = []    
                
    averagedList2 = []           
    AvgListMean = np.mean(DataList)
    for div in DataList:
        div = div-AvgListMean
        averagedList2.append(div)
    
    ## Refer to https://plot.ly/matplotlib/fft/ ##
    ## FFT of iterations for 95% confidence interval ##
    Fs = 1./(HogBin/1000.) #sampling rate; 1 divided by the sampling rate of 60ms used in the experiment.
    n = len(averagedList2) # length of the signal
    k = np.arange(n)
    T = np.float(n/Fs)
    y = averagedList2
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range
    
    Y = np.fft.fft(y)/n
    Y = Y[range(n/2)]
    
    fftGrpList.append(abs(Y))

## Percentile calculation for fft   
fftAvgList = []
CountList = []
fftPercentileList = []
lowerTotal = 0. ; upperTotal = 0.

for access in range(len(Y)):
    CountList = []
    
    for accessList in range(len(fftGrpList)):
        CountList.append(fftGrpList[accessList][access])
        
    fftAvgList.append(CountList)
        
for percent in range(len(fftAvgList)):
    lowerTotal = np.percentile(fftAvgList[percent],2.5)
    upperTotal = np.percentile(fftAvgList[percent],97.5)
        
    fftPercentileList.append([lowerTotal,upperTotal])   
    
## Percentile calculation for replacement data   
GrpAvgList = []
CountList = []
GrpPercentileList = []
lowerTotal = 0. ; upperTotal = 0.

for access in range(len(Y)):
    CountList = []
    
    for accessList in range(len(replacementGrpList)):
        CountList.append(replacementGrpList[accessList][access])
        
    GrpAvgList.append(CountList)
        
for percent in range(len(GrpAvgList)):
    lowerTotal = np.percentile(GrpAvgList[percent],2.5)
    upperTotal = np.percentile(GrpAvgList[percent],97.5)
        
    GrpPercentileList.append([lowerTotal,upperTotal])  
    
NBAvg = np.mean(NBList)  
    
'''        
pylab.figure(figsize=(9,5))
plt.axhline(NBAvg, color = 'green')

x = [0, int(np.max(PSdata1))]
NBerr = stats.sem(NBList)
plt.fill_between(x, NBAvg-NBerr, NBAvg+NBerr,
    alpha=0.75, edgecolor='#7EFF99', facecolor='#7EFF99')

plt.rc('figure', titlesize=15)
pylab.rc('axes', labelsize=15)         

# function for genetic algorithm to minimize (sum of squared error)
# bounds on parameters are set in generate_Initial_Parameters() below
def sumOfSquaredError(parameterTuple):
    warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
    return np.sum((y - HogAnalysis(x, *parameterTuple)) ** 2)

if HogenAnalysis == True:
    def generate_Initial_Parameters():    
        global parameterBounds2
        
        parameterBounds = []
        parameterBounds.append([0.0,1.0]) # parameter bounds for c1; 0.01 arbitrary just to give it a small window to form proper lower and upper bound
        parameterBounds.append([0.0,1.0]) # parameter bounds c3
        parameterBounds.append([0.125, 4.0]) # parameter bounds for Lambda; Nyquist limit, can't detect more than 8Hz in current data set. So 1/8 = 0.125. 1/0.25 = 4.
        parameterBounds.append([0.0, 2*np.pi]) # parameter bounds for theta; Phase offset in radians.
        parameterBounds.append([0.0, 3.0]) # parameter bounds for sigma
    
        parameterBounds2 = ((parameterBounds[0][0], parameterBounds[1][0], parameterBounds[2][0],
                                 parameterBounds[3][0], parameterBounds[4][0]), (parameterBounds[0][1],
                                parameterBounds[1][1], parameterBounds[2][1], parameterBounds[3][1],
                                parameterBounds[4][1]))

        # "seed" the numpy random number generator for repeatable results (N.A. as seed as been removed)
        result = differential_evolution(sumOfSquaredError, parameterBounds)
        return result.x        
                
elif NonDamp == True:
    def generate_Initial_Parameters():    
        global parameterBounds2
        
        parameterBounds = []
        parameterBounds.append([0.0,1.0]) # parameter bounds for c1; 0.01 arbitrary just to give it a small window to form proper lower and upper bound
        parameterBounds.append([0.0,1.0]) # parameter bounds c3
        parameterBounds.append([0.125, 4.0]) # parameter bounds for Lambda; Nyquist limit, can't detect more than 8Hz in current data set. So 1/8 = 0.125. 1/0.25 = 4.
        parameterBounds.append([0.0, 2*np.pi]) # parameter bounds for theta; Phase offset in radians.
    
        parameterBounds2 = ((parameterBounds[0][0], parameterBounds[1][0], parameterBounds[2][0],
                                 parameterBounds[3][0]), (parameterBounds[0][1],
                                parameterBounds[1][1], parameterBounds[2][1], parameterBounds[3][1],
                                ))
        # "seed" the numpy random number generator for repeatable results (N.A. as seed as been removed)
        result = differential_evolution(sumOfSquaredError, parameterBounds)
        return result.x        
    
elif Gaussian == True:
    def generate_Initial_Parameters():    
        global parameterBounds2
        
        parameterBounds = []
        parameterBounds.append([np.mean(averagedList) - 0.01, np.mean(averagedList) + 0.01]) # parameter bounds for c1; 0.01 arbitrary just to give it a small window to form proper lower and upper bound
        parameterBounds.append([np.max(averagedList) - np.mean(averagedList), np.max(averagedList) - np.mean(averagedList) + 0.01]) # parameter bounds c3
        parameterBounds.append([0.0, 1.0]) # parameter bounds for x0 mean.
        parameterBounds.append([0.0, 3.0]) # parameter bounds for sigma
    
        parameterBounds2 = ((parameterBounds[0][0], parameterBounds[1][0], parameterBounds[2][0],
                                 parameterBounds[3][0]), (parameterBounds[0][1],
                                parameterBounds[1][1], parameterBounds[2][1], parameterBounds[3][1],
                                ))
        # "seed" the numpy random number generator for repeatable results (N.A. as seed as been removed)
        result = differential_evolution(sumOfSquaredError, parameterBounds)
        return result.x

## HogAnalysis - To smooth out data
## lambda λ = the wave length, i.e. the inverse of the Frequency 
## theta θ = the time value of the highest first peak in the raw data plot
## sigma σ = dampening of the osciallation for times long after the blink. It’s an open parameter to be fit.

if HogenAnalysis == True:
    def HogAnalysis(x, c1, c3, Lambda, theta, sigma):
        x0 = 0.
        return c1 + c3 * np.cos((2*np.pi*(x/Lambda)) - theta) * np.exp(-(x - x0)**2 / (2 * sigma**2))

if NonDamp == True:    
    def HogAnalysis(x, c1, c3, Lambda, theta):
        return c1 + c3 * np.cos((2*np.pi*(x/Lambda)) - theta)

if Gaussian == True:
    def HogAnalysis(x, c1, c3, x0, sigma):
        return c1 + c3 * np.exp(-(x - x0)**2 / (2 * sigma**2))


## Binned Data
x = np.arange(len(averagedList)) #or rawDataList
lastBin = round(np.max(PSdata1))/1000 #last bin start time in seconds
x = x*(lastBin/len(averagedList)) #Convert to seconds; or rawDataList
x[0] = x[0] + 0.011 #1st bin starts at 11ms
y = np.array(averagedList) #or rawDataList

## Rolling Avg Data
#x = np.arange(len(averagedList))
#x = x/1000.
#y = np.array(averagedList)

# generate initial parameter values
initialParameters = generate_Initial_Parameters()

#popt,pcov = curve_fit(HogAnalysis, x, y, p0=[np.mean(averagedList),np.max(averagedList) - np.mean(averagedList), 0.25,0.0,1.5], bounds=((0., 0., 0. ,0., 0.), (1.0, 1,0, 10.0, 10.0, 10.0)), maxfev=10000)
popt,pcov = curve_fit(HogAnalysis, x, y, initialParameters, method=None, bounds=parameterBounds2, maxfev=50000)

## c1 = est mean performance; c3 = peak of curve - c1; Lambda = 0.75; theta = 0.0; sigma = 1.5
## set x vector in seconds. 1ms = 0.001
## rolling average 100-120ms

print popt

if HogenAnalysis == True:
    print 'c1:', popt[0]
    print 'c3:', popt[1]
    print 'Lambda:', popt[2]
    print 'theta:', popt[3]
    print 'sigma:', popt[4]

if NonDamp == True:
    print 'c1:', popt[0]
    print 'c3:', popt[1]
    print 'lambda:', popt[2]
    print 'theta:', popt[3]

if Gaussian == True:
    print 'c1:', popt[0]
    print 'c3:', popt[1]
    print 'x0:', popt[2]
    print 'sigma:', popt[3]
'''
## FFT ##
averagedList2 = []
AvgListMean = np.mean(averagedList) #or rawDataList
for div in averagedList: #or rawDataList
    div = div-AvgListMean
    averagedList2.append(div)

## Refer to https://plot.ly/matplotlib/fft/ ##
Fs = 1./(HogBin/1000.) #150.0;  # sampling rate
n = len(averagedList2) # length of the signal
k = np.arange(n)
T = np.float(n/Fs)
y = averagedList2
frq = k/T # two sides frequency range
frq = frq[range(n/2)] # one side frequency range

Y = np.fft.fft(y)/n
Y = Y[range(n/2)]

## SD calculation for replacement data   
GrpAvgList = []
CountList = []
GrpSDList = []
lowerSD = 0. ; upperSD = 0.

for access in range(len(Y)):
    CountList = []
    
    for accessList in range(len(replacementGrpList)):
        CountList.append(replacementGrpList[accessList][access])
        
    GrpAvgList.append(CountList)
        
for stdev in range(len(GrpAvgList)):   
    ## Plot bootstrap SD to bootstrap data
    #lowerSD = np.percentile(GrpAvgList[stdev],50.0) - np.std(GrpAvgList[stdev]) #np.percentile(GrpAvgList[stdev],50.0)
    #upperSD = np.percentile(GrpAvgList[stdev],50.0) + np.std(GrpAvgList[stdev])
    #median = np.median(GrpAvgList[stdev]) 
    #GrpSDList.append([lowerSD,upperSD,median]) 
    
    ## Plot bootstrap SD to actual data
    lowerSD = abs(Y[stdev]) - np.std(GrpAvgList[stdev]) 
    upperSD = abs(Y[stdev]) + np.std(GrpAvgList[stdev])
    meann = np.mean(GrpAvgList[stdev])
    GrpSDList.append([lowerSD,upperSD,meann])   
    
## Set a minimum for spectral power SD plot, that lower SD does not go below 0.
for correction in range(len(GrpSDList)):
    if GrpSDList[correction][0] < 0.001:
        GrpSDList[correction][0] = 0.001

# Get rid of 0 Hz data point as it means nothing.
Y = Y[1:]
frq = frq[1:]

plt.plot(frq,abs(Y), color = 'red')

#plt.title('Experiment 1 AB Freq distribution')
#plt.xlabel('Frequency (Hz)')
#plt.ylabel('Spectral Power')

lowPercentileList = []
highPercentileList = []

for plot in range(len(fftPercentileList)):
    lowPercentileList.append(fftPercentileList[plot][0])
    highPercentileList.append(fftPercentileList[plot][1])
    
DatalowPercentileList = []
DatahighPercentileList = []

for plot in range(len(GrpPercentileList)):
    DatalowPercentileList.append(GrpPercentileList[plot][0])
    DatahighPercentileList.append(GrpPercentileList[plot][1])

SDlowList = []
SDhighList = []
meanList= []

for plot in range(len(GrpSDList)):
    SDlowList.append(GrpSDList[plot][0])
    SDhighList.append(GrpSDList[plot][1])
    meanList.append(GrpSDList[plot][2])

# Get rid of 0 Hz data point as it means nothing.
del SDlowList[0]; del lowPercentileList[0]
del SDhighList[0]; del highPercentileList[0]

## To fill in areas in plot where the fft value is more than the 95% conf interval
ShadingList = []
previousSig = 0
for identify in range(len(frq)):
    # if current odds is <=0.05
    if abs(Y[identify]) < lowPercentileList[identify] or abs(Y[identify]) > highPercentileList[identify]:
        print 'HIT 1:', abs(Y[identify]), lowPercentileList[identify], highPercentileList[identify]
        try: 
            # if next odds is also <=0.05, make a note and proceed. If not, shade out the accumulated region that is <=0.05
            if abs(Y[identify+1]) < lowPercentileList[identify+1] or abs(Y[identify+1]) > highPercentileList[identify+1]:
                previousSig += 1
                print 'HIT 1a'
            else:
                if previousSig == 0:
                    # divided by 1000 to convert to sec. divided by 2 to get one tail of window. divided by 2 again to avoid overlapping adjacent bins                    
                    ShadingList.append([frq[identify]-0.2,
                                        frq[identify]+0.2])    
                    print 'HIT 1b'
                else:
                    ShadingList.append([frq[identify-previousSig]-0.2,
                                        frq[identify]+0.2])        
                    previousSig = 0          
                    print 'HIT 1c'
        except:
            # To fix the issue of the last bin, if significant, does not get highlighted
            if abs(Y[identify]) < lowPercentileList[identify] or abs(Y[identify]) > highPercentileList[identify] and identify == len(frq)-1:
                ShadingList.append([frq[identify]-0.2,
                                    frq[identify]+0.2])       
                print 'HIT 1d'           
            continue
        
        else:
            print 'HIT 1e'
            continue
        
for fillSig in range(len(ShadingList)):
    plt.fill_between(ShadingList[fillSig], 0.035, 0.034,
        alpha=1.0, facecolor='red') 
    
#plt.plot(frq,lowPercentileList, color = 'magenta')
#plt.plot(frq,highPercentileList, color = 'magenta', label= '95% Confidence')  
#plt.plot(frq,DatalowPercentileList, color = 'purple')
#plt.plot(frq,DatahighPercentileList, color = 'purple',label= 'Bootstrap')  
#plt.plot(frq,SDlowList, linestyle='dashed', color = 'red')
#plt.plot(frq,SDhighList, linestyle='dashed', color = 'red')
#plt.plot(frq,meanList, color = 'green')

# Hide the right and top spines
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.fill_between(frq, SDlowList, SDhighList,
    alpha=0.15, facecolor='red')

plt.ylim(0,np.max(SDhighList) + 0.005)
plt.yticks(np.arange(0, np.max(SDhighList) + 0.005, step=0.005))

#plt.show()

print '-------------------------------'

#----------------------------------------------------------------------------------------------------------------------
#                                                 NATURAL BLINK SECTION
#----------------------------------------------------------------------------------------------------------------------

#Counters to extract data
GcounterNB = 0; GcounterB0 = 0; GcounterB1 = 0; GcounterB2 = 0; GcounterB3 = 0; GcounterB4 = 0; GcounterB5 = 0
GcounterB6 = 0; GcounterB7 = 0; GcounterB8 = 0; GcounterB9 = 0; GcounterB10 = 0; GcounterB11 = 0; GcounterB12 = 0

counterNB = 0; counterB0 = 0; counterB1 = 0; counterB2 = 0; counterB3 = 0; counterB4 = 0; counterB5 = 0
counterB6 = 0; counterB7 = 0; counterB8 = 0; counterB9 = 0; counterB10 = 0; counterB11 = 0; counterB12 = 0

TcounterNB = 0; TcounterB0 = 0; TcounterB1 = 0; TcounterB2 = 0; TcounterB3 = 0; TcounterB4 = 0; TcounterB5 = 0
TcounterB6 = 0; TcounterB7 = 0; TcounterB8 = 0; TcounterB9 = 0; TcounterB10 = 0; TcounterB11 = 0; TcounterB12 = 0

NBTotal = 0; B0Total = 0; B1Total = 0; B2Total = 0; B3Total = 0; B4Total = 0; B5Total = 0; B6Total = 0
B7Total = 0; B8Total = 0; B9Total = 0; B10Total = 0; B11Total = 0; B12Total = 0

TNBTotal = 0; TB0Total = 0; TB1Total = 0; TB2Total = 0; TB3Total = 0; TB4Total = 0; TB5Total = 0; TB6Total = 0
TB7Total = 0; TB8Total = 0; TB9Total = 0; TB10Total = 0; TB11Total = 0; TB12Total = 0

NBList = []; B0List = []; B1List = []; B2List = []; B3List = []; B4List = []; B5List = []; B6List = []
B7List = []; B8List = []; B9List = []; B10List = []; B11List = []; B12List = []

stderr = []

# Pupil Size Counter
IndPupilSizeList = []
TrialPupilSizeList = []
TBNormalizedList = []
PupilSizeListNat = []

thisfile = 0
part = 0
Sum2 = 0

PSdata1 = []
PSdata2 = []
averagedList = []
ListOfRaw = []

NaturalD0Count = 0
ThresholdCounter = 0
maxLength = 0
currentSet1 = []
currentSet2 = []
IndData = []

NextBin = HogBin
TotalTrialsInBin = 0
CorrectTrialsInBin = 0
Bin60msList = []

ThresholdCounter = 0
NoiseCount = 0
RejectCount = 0
TypoCount = 0

#Paste filename to be analysised here. (Hdf)
if GroupAnalysis == True:
    os.chdir('D:\Dropbox\RSVP\Scene Eye Tracker\hdf5')
    readfile = os.listdir('Cleaned')
    
    start1 = 'eventsSNat0'
    start2 = 'eventsSNatNB0'
    num = 1
    ignorecounter = 0
    for readparticipantfile in range(NrOfPS+len(IgnoreListN)):
        if len(IgnoreListN) != 0:
            if num == IgnoreListN[ignorecounter]:
                num += 1
                ignorecounter += 1
                if ignorecounter > len(IgnoreListN)-1:
                    ignorecounter = len(IgnoreListN)-1
                continue
        
        if num == 10:   
            start1 = 'eventsSNat'
            start2 = 'eventsSNatNB'  
        
        participant = [str('%s%d' % (start1,num)),str('%s%d' % (start2,num))]
        for file in readfile:
            #Append filename to empty list if it ends with .csv
            if file.startswith(participant[0]):
                partset.append(file)
            if file.startswith(participant[1]):
                partset.append(file)
        filenameNatHDF.append(partset)
        num += 1          
        partset = []
    divisions = len(filenameNatHDF)      

if GroupAnalysis == True:
    multiples = 2
    divisions = len(filenameNatHDF)
if IndAnalysis == True:
    multiples = 1
    divisions = 1
    
print filenameNatHDF, len(filenameNatHDF)

for analyse in range(len(filenameNatHDF*multiples)):
    #print part
    #print filenameNatHDF[thisfile][part]
    index1 = 0
    index2 = 0
    IndBlinkList = []    
    
    # Load an ioDataStore file containing 500 Hz sample data from a
    # remote eye tracker that was recording both eyes. In the plotting example
    if GroupAnalysis == True:
        dataAccessUtil=ExperimentDataAccessUtility('D:\Dropbox\RSVP\Scene Eye Tracker\hdf5\Cleaned', filenameNatHDF[thisfile][part], experimentCode=None,sessionCodes=[])
    if IndAnalysis == True:
        dataAccessUtil=ExperimentDataAccessUtility('D:\Dropbox\RSVP\Scene Eye Tracker\hdf5\Cleaned', filenameNatHDF[part], experimentCode=None,sessionCodes=[])    
    
    exp_data = dataAccessUtil.getConditionVariables(filter=None)

    # No need to keep the hdf5 file open anymore...
    #
    dataAccessUtil.close()

    for trial_index,trial_samples in enumerate(exp_data):
        index1 += 1
        
        if trial_samples.SESSION_ID == 0: 
            RejectCount += 1
            continue
        ## Remove trials where blink duration is abnoramlly low, most likely noise which triggered blink condition.
        if trial_samples.blinkend - trial_samples.blinkstart <= 0.011 and trial_samples.blinkend - trial_samples.blinkstart != 0.0:
            #print trial_samples.blinkstart, trial_samples.blinkend
            NoiseCount += 1
            continue
        if trial_samples.response2 < 0 or trial_samples.response2 > 4:
            TypoCount += 1
            #print filenameGogHDF[thisfile][part]
            continue
        
        if trial_samples.blinked == 1 and trial_samples.delay == 0:
            GoggleD0Count += 1
            continue
        #elif trial_samples.blinked == 1 and trial_samples.target_display - trial_samples.blinkstart <= 0.20:  
        elif trial_samples.blinked == 1 and trial_samples.target_display - trial_samples.blinkstart <= 0.20:                  
            B1Total += 1.
            if trial_samples.response2 == trial_samples.target_category:
                counterB1 += 1.            

        if trial_samples.blinked == 0 and trial_samples.target_category != -1:
            NBTotal += 1.
            if trial_samples.response2 == trial_samples.target_category:
                counterNB += 1.      

        elif trial_samples.blinked == 1 and trial_samples.blinkstart != 0. and trial_samples.blinkend != 0:
            PStargonset = trial_samples.target_display*1000 - trial_samples.blinkend*1000

            if PStargonset != 0:
                if trial_samples.response2 == trial_samples.target_category:
                    PSresp = 1.
                elif trial_samples.response2 != trial_samples.target_category:
                    PSresp = 0.                    
                PSdata1.append(PStargonset)    
                PSdata2.append(PSresp)  
                currentSet1.append(PStargonset)
                currentSet2.append(PSresp)               

    if part == 1:
        NBList.append(counterNB/NBTotal)
        counterNB = 0; NBTotal = 0
        part = -1            
        thisfile += 1
        listforsorting = []

        ## Sort individual data and have it appended for error analysis just right after loop function.
        for sort in range(len(currentSet1)):
            objectappend = (currentSet1[sort],currentSet2[sort])
            listforsorting.append(objectappend)
        
        sortedlist = sorted(listforsorting, key=lambda x: x[0])
 
        NextBin = HogBin
        
        ## Original data Rolling Average
        for axisrange in range(int(np.max(currentSet1))):
            if axisrange == NextBin or axisrange == 11: #11 because first bin (D1) starts at 11ms.
                tail1 = axisrange-(averageRange/2)
                tail2 = axisrange+(averageRange/2)
            
                for consolidate in range(len(sortedlist)):
                    if sortedlist[consolidate][0] >= tail1 and sortedlist[consolidate][0] < tail2:
                        ListOfRaw.append(sortedlist[consolidate][1]) 
                        
                if axisrange != 11:
                    NextBin += HogBin
                    
                averagedList.append(np.mean(ListOfRaw))
                ListOfRaw = []     
        
        IndData.append(averagedList)
        ## Identify the longest list out of all participants
        thisMax = len(averagedList)
        if thisMax > maxLength:
            maxLength = thisMax
        currentSet1 = []
        currentSet2 = []
        averagedList = []

    part += 1
    
print 'NATURAL BLINK'   

ValList = []
ErrList = []

## Generate error value for each x point based on combined data of all participant at that x point.
for GenError in range(maxLength):
    for participant in range(len(IndData)):
        try:
            if IndData[participant][GenError] == None or IndData[participant][GenError] == '' or math.isnan(IndData[participant][GenError]):
                continue
            ValList.append(IndData[participant][GenError])
        except:
            continue
    currentErr = stats.sem(ValList)
    ErrList.append(currentErr)

    ValList = []

listforsorting = []

## Sort group data and have it appended for error analysis.
for sort in range(len(PSdata1)):
    objectappend = (PSdata1[sort],PSdata2[sort])
    listforsorting.append(objectappend)

sortedlist = sorted(listforsorting, key=lambda x: x[0])

NextBin = HogBin

## Original data Rolling Average
for axisrange in range(int(np.max(PSdata1))):
    if axisrange == NextBin or axisrange == 11: #11 because first bin (D1) starts at 11ms.
        tail1 = axisrange-(averageRange/2)
        tail2 = axisrange+(averageRange/2)
    
        for consolidate in range(len(sortedlist)):
            if sortedlist[consolidate][0] >= tail1 and sortedlist[consolidate][0] < tail2:
                ListOfRaw.append(sortedlist[consolidate][1]) 
                
        if axisrange != 11:
            NextBin += HogBin
            
        averagedList.append(np.mean(ListOfRaw))
        ListOfRaw = []

NextBin = HogBin
fftGrpList = []
replacementGrpList = []

## Do iterations shuffled data before averaging it out at the end
for iterate in range(NrOfIterations): 
    
    listforsorting = []
    DataList = []
    RepDataList = []
    PSdata2_Scramble = copy.deepcopy(PSdata2)

    ## Merge original data's correct/wrong to the original timings
    for sort in range(len(PSdata1)):
        objectappend = (PSdata1[sort],PSdata2_Scramble[sort])
        listforsorting.append(objectappend)
    
    sortedlist = sorted(listforsorting, key=lambda x: x[0])    
    replacelist = []
    
    ## Replacement of original data
    for replace in range(len(sortedlist)):
        replacelist.append(sortedlist[random.randint(0,len(sortedlist)-1)])
    
    ## NOTE ONE IS replacelist and the other is replaceDlist    
    replacedlist = sorted(replacelist, key=lambda x: x[0])       
    
    NextBin = HogBin
    
    ## Original data replacement
    for axisrange in range(int(np.max(PSdata1))):
        if axisrange == NextBin or axisrange == 11: #11 because first bin (D1) starts at 11ms.
            tail1 = axisrange-(averageRange/2)
            tail2 = axisrange+(averageRange/2)
        
            for consolidate in range(len(sortedlist)):
                if replacedlist[consolidate][0] >= tail1 and replacedlist[consolidate][0] < tail2:
                    ListOfRaw.append(replacedlist[consolidate][1]) 
                    
            if axisrange != 11:
                NextBin += HogBin
                
            RepDataList.append(np.mean(ListOfRaw))
            ListOfRaw = []  
            
            
    averagedList2 = []           
    AvgListMean = np.mean(RepDataList)
    for div in RepDataList:
        div = div-AvgListMean
        averagedList2.append(div)
    
    ## Refer to https://plot.ly/matplotlib/fft/ ##
    ## FFT of iterations of original replacement data ##
    Fs = 1./(HogBin/1000.) #sampling rate; 1 divided by the sampling rate of 60ms used in the experiment.
    n = len(averagedList2) # length of the signal
    k = np.arange(n)
    T = np.float(n/Fs)
    y = averagedList2
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range
    
    Y = np.fft.fft(y)/n
    Y = Y[range(n/2)]
    
    replacementGrpList.append(abs(Y))
            
    ## ------------------------- Confidence Interval ------------------------------
    ## Randomise correct and wrong
    random.shuffle(PSdata2_Scramble)
    listforsorting = []
    
    ## Merge randomised correct/wrong to the original timings
    for sort in range(len(PSdata1)):
        objectappend = (PSdata1[sort],PSdata2_Scramble[sort])
        listforsorting.append(objectappend)
    
    sortedlist = sorted(listforsorting, key=lambda x: x[0])
    
    NextBin = HogBin
    
    ## Shuffled data Rolling Average
    for axisrange in range(int(np.max(PSdata1))):
        if axisrange == NextBin or axisrange == 11: #11 because first bin (D1) starts at 11ms.
            tail1 = axisrange-(averageRange/2)
            tail2 = axisrange+(averageRange/2)
        
            for consolidate in range(len(sortedlist)):
                if sortedlist[consolidate][0] >= tail1 and sortedlist[consolidate][0] < tail2:
                    ListOfRaw.append(sortedlist[consolidate][1]) 
                    
            if axisrange != 11:
                NextBin += HogBin
                
            DataList.append(np.mean(ListOfRaw))
            ListOfRaw = []    
                
    averagedList2 = []           
    AvgListMean = np.mean(DataList)
    for div in DataList:
        div = div-AvgListMean
        averagedList2.append(div)
    
    ## Refer to https://plot.ly/matplotlib/fft/ ##
    ## FFT of iterations for 95% confidence interval ##
    Fs = 1./(HogBin/1000.) #sampling rate; 1 divided by the sampling rate of 60ms used in the experiment.
    n = len(averagedList2) # length of the signal
    k = np.arange(n)
    T = np.float(n/Fs)
    y = averagedList2
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range
    
    Y = np.fft.fft(y)/n
    Y = Y[range(n/2)]
    
    fftGrpList.append(abs(Y))

## Percentile calculation for fft   
fftAvgList = []
CountList = []
fftPercentileList = []
lowerCount = 0. ; upperCount = 0.
lowerTotal = 0. ; upperTotal = 0.

for access in range(len(Y)):
    CountList = []
    
    for accessList in range(len(fftGrpList)):
        CountList.append(fftGrpList[accessList][access])
        
    fftAvgList.append(CountList)
        
for percent in range(len(fftAvgList)):
    lowerTotal = np.percentile(fftAvgList[percent],2.5)
    upperTotal = np.percentile(fftAvgList[percent],97.5)
        
    fftPercentileList.append([lowerTotal,upperTotal])   
    
## Percentile calculation for replacement data   
GrpAvgList = []
CountList = []
GrpPercentileList = []
lowerCount = 0. ; upperCount = 0.
lowerTotal = 0. ; upperTotal = 0.

for access in range(len(Y)):
    CountList = []
    
    for accessList in range(len(replacementGrpList)):
        CountList.append(replacementGrpList[accessList][access])
        
    GrpAvgList.append(CountList)
        
for percent in range(len(GrpAvgList)):
    lowerTotal = np.percentile(GrpAvgList[percent],2.5)
    upperTotal = np.percentile(GrpAvgList[percent],97.5)
        
    GrpPercentileList.append([lowerTotal,upperTotal])  
    
NBAvg = np.mean(NBList)  
    
'''        
pylab.figure(figsize=(9,5))
plt.axhline(NBAvg, color = 'green')

x = [0, int(np.max(PSdata1))]
NBerr = stats.sem(NBList)
plt.fill_between(x, NBAvg-NBerr, NBAvg+NBerr,
    alpha=0.75, edgecolor='#7EFF99', facecolor='#7EFF99')

plt.rc('figure', titlesize=15)
pylab.rc('axes', labelsize=15)         

# function for genetic algorithm to minimize (sum of squared error)
# bounds on parameters are set in generate_Initial_Parameters() below
def sumOfSquaredError(parameterTuple):
    warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
    return np.sum((y - HogAnalysis(x, *parameterTuple)) ** 2)

if HogenAnalysis == True:
    def generate_Initial_Parameters():    
        global parameterBounds2
        
        parameterBounds = []
        parameterBounds.append([0.0,1.0]) # parameter bounds for c1; 0.01 arbitrary just to give it a small window to form proper lower and upper bound
        parameterBounds.append([0.0,1.0]) # parameter bounds c3
        parameterBounds.append([0.125, 4.0]) # parameter bounds for Lambda; Nyquist limit, can't detect more than 8Hz in current data set. So 1/8 = 0.125. 1/0.25 = 4.
        parameterBounds.append([0.0, 2*np.pi]) # parameter bounds for theta; Phase offset in radians.
        parameterBounds.append([0.0, 3.0]) # parameter bounds for sigma
    
        parameterBounds2 = ((parameterBounds[0][0], parameterBounds[1][0], parameterBounds[2][0],
                                 parameterBounds[3][0], parameterBounds[4][0]), (parameterBounds[0][1],
                                parameterBounds[1][1], parameterBounds[2][1], parameterBounds[3][1],
                                parameterBounds[4][1]))

        # "seed" the numpy random number generator for repeatable results (N.A. as seed as been removed)
        result = differential_evolution(sumOfSquaredError, parameterBounds)
        return result.x        
                
elif NonDamp == True:
    def generate_Initial_Parameters():    
        global parameterBounds2
        
        parameterBounds = []
        parameterBounds.append([0.0,1.0]) # parameter bounds for c1; 0.01 arbitrary just to give it a small window to form proper lower and upper bound
        parameterBounds.append([0.0,1.0]) # parameter bounds c3
        parameterBounds.append([0.125, 4.0]) # parameter bounds for Lambda; Nyquist limit, can't detect more than 8Hz in current data set. So 1/8 = 0.125. 1/0.25 = 4.
        parameterBounds.append([0.0, 2*np.pi]) # parameter bounds for theta; Phase offset in radians.
    
        parameterBounds2 = ((parameterBounds[0][0], parameterBounds[1][0], parameterBounds[2][0],
                                 parameterBounds[3][0]), (parameterBounds[0][1],
                                parameterBounds[1][1], parameterBounds[2][1], parameterBounds[3][1],
                                ))
        # "seed" the numpy random number generator for repeatable results (N.A. as seed as been removed)
        result = differential_evolution(sumOfSquaredError, parameterBounds)
        return result.x        
    
elif Gaussian == True:
    def generate_Initial_Parameters():    
        global parameterBounds2
        
        parameterBounds = []
        parameterBounds.append([np.mean(averagedList) - 0.01, np.mean(averagedList) + 0.01]) # parameter bounds for c1; 0.01 arbitrary just to give it a small window to form proper lower and upper bound
        parameterBounds.append([np.max(averagedList) - np.mean(averagedList), np.max(averagedList) - np.mean(averagedList) + 0.01]) # parameter bounds c3
        parameterBounds.append([0.0, 1.0]) # parameter bounds for x0 mean.
        parameterBounds.append([0.0, 3.0]) # parameter bounds for sigma
    
        parameterBounds2 = ((parameterBounds[0][0], parameterBounds[1][0], parameterBounds[2][0],
                                 parameterBounds[3][0]), (parameterBounds[0][1],
                                parameterBounds[1][1], parameterBounds[2][1], parameterBounds[3][1],
                                ))
        # "seed" the numpy random number generator for repeatable results (N.A. as seed as been removed)
        result = differential_evolution(sumOfSquaredError, parameterBounds)
        return result.x

## HogAnalysis - To smooth out data
## lambda λ = the wave length, i.e. the inverse of the Frequency 
## theta θ = the time value of the highest first peak in the raw data plot
## sigma σ = dampening of the osciallation for times long after the blink. It’s an open parameter to be fit.

if HogenAnalysis == True:
    def HogAnalysis(x, c1, c3, Lambda, theta, sigma):
        x0 = 0.
        return c1 + c3 * np.cos((2*np.pi*(x/Lambda)) - theta) * np.exp(-(x - x0)**2 / (2 * sigma**2))

elif NonDamp == True:    
    def HogAnalysis(x, c1, c3, Lambda, theta):
        return c1 + c3 * np.cos((2*np.pi*(x/Lambda)) - theta)

elif Gaussian == True:
    def HogAnalysis(x, c1, c3, x0, sigma):
        return c1 + c3 * np.exp(-(x - x0)**2 / (2 * sigma**2))


## Binned Data
x = np.arange(len(averagedList)) #or rawDataList
lastBin = round(np.max(PSdata1))/1000 #last bin start time in seconds
x = x*(lastBin/len(averagedList)) #Convert to seconds; or rawDataList
x[0] = x[0] + 0.011 #1st bin starts at 11ms
y = np.array(averagedList) #or rawDataList

## Rolling Avg Data
#x = np.arange(len(averagedList))
#x = x/1000.
#y = np.array(averagedList)

# generate initial parameter values
initialParameters = generate_Initial_Parameters()

#popt,pcov = curve_fit(HogAnalysis, x, y, p0=[np.mean(averagedList),np.max(averagedList) - np.mean(averagedList), 0.25,0.0,1.5], bounds=((0., 0., 0. ,0., 0.), (1.0, 1,0, 10.0, 10.0, 10.0)), maxfev=10000)
popt,pcov = curve_fit(HogAnalysis, x, y, initialParameters, method=None, bounds=parameterBounds2, maxfev=50000)

## c1 = est mean performance; c3 = peak of curve - c1; Lambda = 0.75; theta = 0.0; sigma = 1.5
## set x vector in seconds. 1ms = 0.001
## rolling average 100-120ms

print popt

if HogenAnalysis == True:
    print 'c1:', popt[0]
    print 'c3:', popt[1]
    print 'Lambda:', popt[2]
    print 'theta:', popt[3]
    print 'sigma:', popt[4]

elif NonDamp == True:
    print 'c1:', popt[0]
    print 'c3:', popt[1]
    print 'lambda:', popt[2]
    print 'theta:', popt[3]

elif Gaussian == True:
    print 'c1:', popt[0]
    print 'c3:', popt[1]
    print 'x0:', popt[2]
    print 'sigma:', popt[3]
'''
## FFT ##
averagedList2 = []
AvgListMean = np.mean(averagedList) #or rawDataList
for div in averagedList: #or rawDataList
    div = div-AvgListMean
    averagedList2.append(div)

## Refer to https://plot.ly/matplotlib/fft/ ##
Fs = 1./(HogBin/1000.) #150.0;  # sampling rate
n = len(averagedList2) # length of the signal
k = np.arange(n)
T = np.float(n/Fs)
y = averagedList2
frq = k/T # two sides frequency range
frq = frq[range(n/2)] # one side frequency range

Y = np.fft.fft(y)/n
Y = Y[range(n/2)]

## SD calculation for replacement data   
GrpAvgList = []
CountList = []
GrpSDList = []
lowerSD = 0. ; upperSD = 0.

for access in range(len(Y)):
    CountList = []
    
    for accessList in range(len(replacementGrpList)):
        CountList.append(replacementGrpList[accessList][access])
        
    GrpAvgList.append(CountList)
        
for stdev in range(len(GrpAvgList)):   
    ## Plot bootstrap SD to bootstrap data
    #lowerSD = np.percentile(GrpAvgList[stdev],50.0) - np.std(GrpAvgList[stdev]) #np.percentile(GrpAvgList[stdev],50.0)
    #upperSD = np.percentile(GrpAvgList[stdev],50.0) + np.std(GrpAvgList[stdev])
    #median = np.median(GrpAvgList[stdev]) 
    #GrpSDList.append([lowerSD,upperSD,median]) 
    
    ## Plot bootstrap SD to actual data
    lowerSD = abs(Y[stdev]) - np.std(GrpAvgList[stdev]) 
    upperSD = abs(Y[stdev]) + np.std(GrpAvgList[stdev])
    meann = np.mean(GrpAvgList[stdev])
    GrpSDList.append([lowerSD,upperSD,meann])   
    
## Set a minimum for spectral power SD plot, that lower SD does not go below 0.
for correction in range(len(GrpSDList)):
    if GrpSDList[correction][0] < 0.001:
        GrpSDList[correction][0] = 0.001

# Get rid of 0 Hz data point as it means nothing.
Y = Y[1:]
frq = frq[1:]

plt.plot(frq,abs(Y), color = 'blue')

#plt.title('Experiment 1 VB Freq distribution')
#plt.xlabel('Frequency (Hz)')
#plt.ylabel('Spectral Power')

lowPercentileList = []
highPercentileList = []

for plot in range(len(fftPercentileList)):
    lowPercentileList.append(fftPercentileList[plot][0])
    highPercentileList.append(fftPercentileList[plot][1])
    
DatalowPercentileList = []
DatahighPercentileList = []

for plot in range(len(GrpPercentileList)):
    DatalowPercentileList.append(GrpPercentileList[plot][0])
    DatahighPercentileList.append(GrpPercentileList[plot][1])

SDlowList = []
SDhighList = []
meanList= []

for plot in range(len(GrpSDList)):
    SDlowList.append(GrpSDList[plot][0])
    SDhighList.append(GrpSDList[plot][1])
    meanList.append(GrpSDList[plot][2])
    
# Get rid of 0 Hz data point as it means nothing.
del SDlowList[0]; del lowPercentileList[0]
del SDhighList[0]; del highPercentileList[0]

## To fill in areas in plot where the fft value is more than the 95% conf interval
ShadingList2 = []
previousSig = 0
for identify in range(len(frq)):
    # if current odds is <=0.05
    if abs(Y[identify]) < lowPercentileList[identify] or abs(Y[identify]) > highPercentileList[identify]:
        print 'HIT 1:', abs(Y[identify]), lowPercentileList[identify], highPercentileList[identify]
        try: 
            # if next odds is also <=0.05, make a note and proceed. If not, shade out the accumulated region that is <=0.05
            if abs(Y[identify+1]) < lowPercentileList[identify+1] or abs(Y[identify+1]) > highPercentileList[identify+1]:
                previousSig += 1
                print 'HIT 1a'
            else:
                if previousSig == 0:
                    # divided by 1000 to convert to sec. divided by 2 to get one tail of window. divided by 2 again to avoid overlapping adjacent bins                    
                    ShadingList2.append([frq[identify]-0.2,
                                        frq[identify]+0.2])    
                    print 'HIT 1b'
                else:
                    ShadingList2.append([frq[identify-previousSig]-0.2,
                                        frq[identify]+0.2])        
                    previousSig = 0          
                    print 'HIT 1c'
        except:
            # To fix the issue of the last bin, if significant, does not get highlighted
            if abs(Y[identify]) < lowPercentileList[identify] or abs(Y[identify]) > highPercentileList[identify] and identify == len(frq)-1:
                ShadingList2.append([frq[identify]-0.2,
                                    frq[identify]+0.2])       
                print 'HIT 1d'           
            continue
        
        else:
            print 'HIT 1e'
            continue
        
for fillSig in range(len(ShadingList2)):
    plt.fill_between(ShadingList2[fillSig], 0.033, 0.032,
        alpha=1.0, facecolor='blue') 
    
#plt.plot(frq,lowPercentileList, color = 'magenta')
#plt.plot(frq,highPercentileList, color = 'magenta', label= '95% Confidence')  
#plt.plot(frq,DatalowPercentileList, color = 'purple')
#plt.plot(frq,DatahighPercentileList, color = 'purple',label= 'Bootstrap')  
#plt.plot(frq,SDlowList, linestyle='dashed', color = 'blue')
#plt.plot(frq,SDhighList, linestyle='dashed', color = 'blue')
#plt.plot(frq,meanList, color = 'green')

plt.xlim(0.4,8)

# Hide the right and top spines
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.fill_between(frq, SDlowList, SDhighList,
    alpha=0.15, facecolor='blue')

#plt.savefig('D:\Dropbox\RSVP\RSVP report\Publication\Fig 3b.pdf', format='pdf', dpi=1000, bbox_inches="tight")
plt.show()
