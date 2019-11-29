# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 12:14:23 2019

@author: Aaron Ang
"""
# NOTE: This is not a working script. Only explanatory.

## ----------------------------------------------------------------------------
## FFT ##
## This section generates the FFT signal of your original data.
## ----------------------------------------------------------------------------

## Appendix
## --------
## averagedList = Your data 

averagedList2 = []
AvgListMean = np.mean(averagedList) # attain the mean score of all your data 

# You now subtract all your data by the mean score to normalise them to zero.
for div in averagedList: 
    div = div-AvgListMean
    averagedList2.append(div)

## Refer to https://plot.ly/matplotlib/fft/ ## UNFORTUNATELY, THEY REMOVED THE PAGE. 
Fs = 1./(HogBin/1000.) # sampling rate. HogBin is 60ms in my case, as my bin sizes were in 60ms.
n = len(averagedList2) # length of the signal(data)
k = np.arange(n)
T = np.float(n/Fs)
y = averagedList2
frq = k/T # two sides frequency range
frq = frq[range(n/2)] # one side frequency range

Y = np.fft.fft(y)/n
Y = Y[range(n/2)]

# Get rid of 0 Hz data point as it means nothing.
Y = Y[1:]
frq = frq[1:]

## ----------------------------------------------------------------------------
## FFT Steps
## - Though the FFT codes look simple, the process of generating it with the bootstrap data is 
##   complicated. It's quite difficult to copy and paste snippets from my code and translate it here.
##   But hopefully, the following step-by-step list of what I did would help you understand by codes better.
##
##   Step 1) Generate FFT of original data
##   Step 2) Do bootstrap with replacement at 10000 iteration of your data
##   Step 3) Generate the FFT plot for each of your 10000 iteration, append it to a list
##   Step 4) Calculate lower (2.5%) and upper (97.5%) percentile of your list of 10000 FFT outcome at each time point. (Refer to "## Percentile calculation for fft" in original script)
##   Step 5) Plot your original FFT data together with the lower and upper percentile
##   Step 6) Any point of your original FFT data that surpasses the upper percentile means there is significant difference.