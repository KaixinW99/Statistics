#import the package
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as SciSig
import scipy.interpolate as SciInt
from matplotlib.patches import Ellipse

# You can minimize all the Plot function

# Read the data
dt = 1/4096
dataH = np.loadtxt('GW150914_H.dat',skiprows=3)
dataL = np.loadtxt('GW150914_L.dat',skiprows=3)
timeH = np.arange(0,len(dataH)*dt,dt)
timeL = np.arange(0,len(dataL)*dt,dt)

def RawTimePlot():
    # Plot of raw time series
    fig, axes = plt.subplots(2,1,figsize=[8,8])
    axes[0].plot(timeH,dataH,label='GW150914_H')
    axes[1].plot(timeL,dataL,label='GW150914_L')
    for ax in axes.flat:
        ax.legend()
        ax.set(xlabel='Time(s)',ylabel='Intensity(a.u.)')
    fig.suptitle('Raw Data: Time Series')
    #plt.show()
    return fig

# Fast Fourier Transform for real input
funcH = np.fft.rfft(dataH)
funcL = np.fft.rfft(dataL)
freqH = np.fft.rfftfreq(dataH.size,dt)
freqL = np.fft.rfftfreq(dataL.size,dt)

def FreqPlot():
    # Plot the frequency series
    fig, axes = plt.subplots(2,2,figsize=[16,8])
    axes[0,0].plot(freqH,np.absolute(funcH),label='GW150914_H')
    axes[1,0].plot(freqH,np.absolute(funcH),label='GW150914_H (LogY)')
    axes[1,0].set(yscale='log')
    axes[0,1].plot(freqL,np.absolute(funcL),label='GW150914_L')
    axes[1,1].plot(freqL,np.absolute(funcL),label='GW150914_L (LogY)')
    axes[1,1].set(yscale='log')
    for ax in axes.flat:
        ax.legend()
        ax.set(xlabel='Frequency(Hz)',ylabel='Intensity(a.u.)')
    fig.suptitle('Fourier Transform: Frequency Series (Plot and LogPlot)')
    #plt.show()
    return fig

# Power spectral density calculation and interpolation
PSDlistH = SciSig.welch(dataH,fs=1/dt,nperseg=4096)
PSDlistL = SciSig.welch(dataL,fs=1/dt,nperseg=4096)

PSDfH = SciInt.interp1d(PSDlistH[0],PSDlistH[1])
PSDfL = SciInt.interp1d(PSDlistL[0],PSDlistL[1])

#When I increase the "nperseg", the final plot will be more similar to the pictures on the paper.
#Although I do not know exactly why, that corresponds to the values of peak in a much larger window.
#The data with a much larger window will keep small frequency in the time series.
#So when I applied smaller window, my final plot was dominated by the high frequency 
#and the smaller ones were filtered.

#My guess is that the "nperseg" fundamentally depends on the distance between the blackhole and us
#due to the doppler effect and Hubble Law.
#If we know the red shift, maybe we can have a much clearer data.

# Filter the freq series with PSD
funcH_F = np.divide(funcH,np.sqrt(PSDfH(freqH)))
funcL_F = np.divide(funcL,np.sqrt(PSDfL(freqL)))

def FilterFuncPSD():
    # Plot the Filtered Function
    fig, axes = plt.subplots(2,2,figsize=[16,8])
    axes[0,0].plot(freqH,1/np.sqrt(PSDfH(freqH)),label='GW150914_H')
    axes[1,0].plot(freqH,1/np.sqrt(PSDfH(freqH)),label='GW150914_H (LogY)')
    axes[1,0].set(yscale='log')
    axes[0,1].plot(freqL,1/np.sqrt(PSDfL(freqL)),label='GW150914_L')
    axes[1,1].plot(freqL,1/np.sqrt(PSDfL(freqL)),label='GW150914_L (LogY)')
    axes[1,1].set(yscale='log')
    for ax in axes.flat:
        ax.legend()
        ax.set(xlabel='Frequency(Hz)',ylabel='Intensity(a.u.)')
    fig.suptitle('1/sqrt(PSD)')
    #plt.show()
    return fig    

def FilterFreqPlot():
    # Plot the Filtered freq series
    fig, axes = plt.subplots(2,2,figsize=[16,8])
    axes[0,0].plot(freqH,np.absolute(funcH_F),label='GW150914_H')
    axes[1,0].plot(freqH,np.absolute(funcH_F),label='GW150914_H (LogY)')
    axes[1,0].set(yscale='log')
    axes[0,1].plot(freqL,np.absolute(funcL_F),label='GW150914_L')
    axes[1,1].plot(freqL,np.absolute(funcL_F),label='GW150914_L (LogY)')
    axes[1,1].set(yscale='log')
    for ax in axes.flat:
        ax.legend()
        ax.set(xlabel='Frequency(Hz)',ylabel='Intensity(a.u.)')
    fig.suptitle('Frequency Series Filtered with PSD')
    #plt.show()
    return fig

# truncate the freq function
funcH_F_C=funcH_F.copy(); funcH_F_C[(freqH < 25) | (freqH > 300)] = 0
funcL_F_C=funcL_F.copy(); funcL_F_C[(freqL < 25) | (freqL > 300)] = 0

def FilterTrunFreqPlot():
    # Plot the Filtered truncated freq series
    fig, axes = plt.subplots(2,2,figsize=[16,8])
    axes[0,0].plot(freqH[(freqH >= 20) & (freqH <= 305)],np.absolute(funcH_F_C)[(freqH >= 20) & (freqH <= 305)],label='GW150914_H')
    axes[1,0].plot(freqH[(freqH >= 20) & (freqH <= 305)],np.absolute(funcH_F_C)[(freqH >= 20) & (freqH <= 305)],label='GW150914_H (LogY)')
    axes[1,0].set(yscale='log')
    axes[0,1].plot(freqL[(freqL >= 20) & (freqL <= 305)],np.absolute(funcL_F_C)[(freqL >= 20) & (freqL <= 305)],label='GW150914_L')
    axes[1,1].plot(freqL[(freqL >= 20) & (freqL <= 305)],np.absolute(funcL_F_C)[(freqL >= 20) & (freqL <= 305)],label='GW150914_L (LogY)')
    axes[1,1].set(yscale='log')
    for ax in axes.flat:
        ax.legend()
        ax.set(xlabel='Frequency(Hz)',ylabel='Intensity(a.u.)')
    fig.suptitle('Truncated Frequency Series Filtered with PSD')
    #plt.show()
    return fig

# Inverse Fast Fourier Transform
newfuncH = np.fft.irfft(funcH_F_C)
newfuncL = np.fft.irfft(funcL_F_C)

def FilterTrunTimePlot():
    # Plot the Filtered truncated time series
    fig, axes = plt.subplots(2,2,figsize=[16,8])
    axes[0,0].plot(timeH,dataH,label='GW150914_H')
    axes[0,1].plot(timeL,dataL,label='GW150914_L')
    axes[1,0].plot(timeH,newfuncH,label='GW150914_H (Filtered/Truncated)')
    axes[1,1].plot(timeL,newfuncL,label='GW150914_L (Filtered/Truncated)')
    for ax in axes.flat:
        ax.set_xlabel('Time(s)')
        ax.set_ylabel('Intensity(a.u.)')
        ax.legend()
    fig.suptitle('Inverse Fourier Transform: Time Series from Truncated/Filtered Frequency Series')
    #Put a red ellipse around the outlier so the viewer can identify it
    E1 = Ellipse(xy=(16.425,0), width=1, height=500, ec='r', fill=False, lw=2,  clip_on=False)
    E2 = Ellipse(xy=(15.425,0), width=1, height=500, ec='r', fill=False, lw=2,  clip_on=False)
    axes[1,0].add_patch(E1)
    axes[1,1].add_patch(E2)
    axes[1,0].annotate('The important peak in gravitational wave',
        xy=(16.4, 300), xycoords='data',
        xytext=(-70, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"))
    axes[1,1].annotate('The important peak in gravitational wave',
        xy=(15.4, 300), xycoords='data',
        xytext=(-70, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"))
    #plt.show()
    return fig

# Calculate the standard deviation of time series
startH = int(newfuncH.size/(len(dataH)*dt))
startL = int(newfuncL.size/(len(dataL)*dt))
endH = int(newfuncH.size*(1-1/(len(dataH)*dt)))
endL = int(newfuncL.size*(1-1/(len(dataL)*dt)))
#print(startH*dt,endH*dt)
#print(startL*dt,endL*dt)
stdH = np.std(newfuncH[startH:endH+1])
stdL = np.std(newfuncL[startL:endL+1])

def FilterTrunTimePlot1to31():
    # PLot Filtered truncated time series from 1 to 31 s
    fig, axes = plt.subplots(2,2,figsize=[16,8])
    axes[0,0].plot(timeH,dataH,label='GW150914_H')
    axes[0,1].plot(timeL,dataL,label='GW150914_L')
    axes[1,0].plot(timeH[startH:endH+1],newfuncH[startH:endH+1],label='GW150914_H (Filtered/Truncated)')
    #axes[1,0].set(xlim=[16.2,16.5],ylim=[-500,500])  #change it to zoom in/out
    axes[1,1].plot(timeL[startL:endL+1],newfuncL[startL:endL+1],label='GW150914_L (Filtered/Truncated)')
    #axes[1,1].set(xlim=[15.3,15.6],ylim=[-450,450])  #change it to zoom in/out
    for ax in axes.flat:
        ax.set_xlabel('Time(s)')
        ax.set_ylabel('Intensity(a.u.)')
        ax.legend()
    fig.suptitle('Time Series from 1 to 31 seconds')
    #plt.show()
    return fig

# change the unit of yaxis from a.u. to Noise standard deviation
FinalFunH = np.divide(newfuncH,stdH)
FinalFunL = np.divide(newfuncL,stdL)

def FilterTrunTimePlotSTD():
    # Plot the Filtered truncated time series (red ellipses are our focus)
    fig, axes = plt.subplots(2,2,figsize=[16,8])
    axes[0,0].plot(timeH,dataH,label='GW150914_H')
    axes[0,1].plot(timeL,dataL,label='GW150914_L')
    axes[1,0].plot(timeH,FinalFunH,label='GW150914_H (Filtered/Truncated)')
    #axes[1,0].set(xlim=[16.2,16.5],ylim=[-500,500])  #change it to zoom in/out
    axes[1,1].plot(timeL,FinalFunL,label='GW150914_L (Filtered/Truncated)')
    #axes[1,1].set(xlim=[15.3,15.6],ylim=[-450,450])  #change it to zoom in/out
    axes[0,0].set_ylabel('Intensity(a.u)')
    axes[0,1].set_ylabel('Intensity(a.u)')
    axes[1,0].set_ylabel('Intensity(Noise $\sigma$)')
    axes[1,1].set_ylabel('Intensity(Noise $\sigma$)')
    for ax in axes.flat:
        ax.set_xlabel('Time(s)')
        ax.legend()
    #Put a red ellipse around the outlier so the viewer can identify it
    E1 = Ellipse(xy=(16.425,0), width=1, height=50, ec='r', fill=False, lw=2,  clip_on=False)
    E2 = Ellipse(xy=(15.425,0), width=1, height=50, ec='r', fill=False, lw=2,  clip_on=False)
    axes[1,0].add_patch(E1)
    axes[1,1].add_patch(E2)
    axes[1,0].annotate('The important peak in gravitational wave',
        xy=(16.4, 30), xycoords='data',
        xytext=(-70, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"))
    axes[1,1].annotate('The important peak in gravitational wave',
        xy=(15.4, 30), xycoords='data',
        xytext=(-70, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"))    
    fig.suptitle('Time series with Standard Deviation Unit')
    #plt.show()
    return fig

def FilterTrunTimePlotSTD1to31():
    # Plot the Filtered truncated time series with Noise std unit from 1 to 31 s
    fig, axes = plt.subplots(2,2,figsize=[16,8])
    axes[0,0].plot(timeH,dataH,label='GW150914_H')
    axes[0,1].plot(timeL,dataL,label='GW150914_L')
    axes[1,0].plot(timeH[startH:endH+1],FinalFunH[startH:endH+1],label='GW150914_H (Filtered/Truncated)')
    #axes[1,0].set(xlim=[16.2,16.5],ylim=[-500,500])  #change it to zoom in/out
    axes[1,1].plot(timeL[startH:endH+1],FinalFunL[startL:endL+1],label='GW150914_L (Filtered/Truncated)')
    #axes[1,1].set(xlim=[15.3,15.6],ylim=[-450,450])  #change it to zoom in/out
    axes[0,0].set_ylabel('Intensity(a.u)')
    axes[0,1].set_ylabel('Intensity(a.u)')
    axes[1,0].set_ylabel('Intensity(Noise $\sigma$)')
    axes[1,1].set_ylabel('Intensity(Noise $\sigma$)')
    for ax in axes.flat:
        ax.set_xlabel('Time(s)')
        ax.legend()
    fig.suptitle('Time Series with Standard Deviation Unit from 1 t 31 seconds')
    #plt.show()
    return fig

def FinalPlot():
    # Plot the Filtered truncated time series with Noise std unit in detail
    fig, axes = plt.subplots(2,2,figsize=[16,8])
    axes[0,0].plot(timeH,dataH,label='GW150914_H')
    axes[0,1].plot(timeL,dataL,label='GW150914_L')
    axes[1,0].plot(timeH,FinalFunH,label='GW150914_H (Filtered/Truncated)')
    axes[1,0].set(xlim=[16.25,16.46],ylim=[-10,10])  #change it to zoom in/out
    axes[1,1].plot(timeL,FinalFunL,label='GW150914_L (Filtered/Truncated)')
    axes[1,1].set(xlim=[15.24,15.45],ylim=[-10,10])  #change it to zoom in/out
    axes[0,0].set_ylabel('Intensity(a.u)')
    axes[0,1].set_ylabel('Intensity(a.u)')
    axes[1,0].set_ylabel('Intensity(Noise $\sigma$)')
    axes[1,1].set_ylabel('Intensity(Noise $\sigma$)')
    for ax in axes.flat:
        ax.set_xlabel('Time(s)')
        ax.legend()
    fig.suptitle('Time Series (zoom in)')
    #plt.show()
    return fig

def FinalComparePlot():
    # Plot the Filtered truncated time series with Noise std unit in detail
    fig, axes = plt.subplots(3,2,figsize=[16,9])
    axes[0,0].plot(timeH,dataH,label='GW150914_H')
    axes[0,1].plot(timeL,dataL,label='GW150914_L')
    axes[1,0].plot(timeH,FinalFunH,label='GW150914_H (Filtered/Truncated)')
    axes[2,0].plot(timeH,FinalFunH,label='GW150914_H (Filtered/Truncated)')
    axes[2,0].set(xlim=[16.25,16.46],ylim=[-10,10])  #change it to zoom in/out
    axes[1,1].plot(timeL,FinalFunL,label='GW150914_L (Filtered/Truncated)')
    axes[2,1].plot(timeL,FinalFunL,label='GW150914_L (Filtered/Truncated)')
    axes[2,1].set(xlim=[15.24,15.45],ylim=[-10,10])  #change it to zoom in/out
    axes[0,0].set_ylabel('Intensity(a.u)')
    axes[0,1].set_ylabel('Intensity(a.u)')
    axes[1,0].set_ylabel('Intensity(Noise $\sigma$)')
    axes[1,1].set_ylabel('Intensity(Noise $\sigma$)')
    axes[2,0].set_ylabel('Intensity(Noise $\sigma$)')
    axes[2,1].set_ylabel('Intensity(Noise $\sigma$)')
    for ax in axes.flat:
        ax.legend()
    E1 = Ellipse(xy=(16.425,0), width=1, height=50, ec='r', fill=False, lw=2,  clip_on=False)
    E2 = Ellipse(xy=(15.425,0), width=1, height=50, ec='r', fill=False, lw=2,  clip_on=False)
    axes[1,0].add_patch(E1)
    axes[1,1].add_patch(E2)
    axes[1,0].annotate('The important peak in gravitational wave (zoom in)',
        xy=(16.4, 30), xycoords='data',
        xytext=(-70, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"))
    axes[1,1].annotate('The important peak in gravitational wave (zoom in)',
        xy=(15.4, 30), xycoords='data',
        xytext=(-70, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"))    
    fig.suptitle('Time Series (Comparsion)')
    fig.supxlabel('Time(s)')
    plt.show()
    return fig

# Plot multiple pictures in one pdf file
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('hw04a.pdf')
pp.savefig(RawTimePlot())                       # Save and show the Raw data vs time
pp.savefig(FreqPlot())                          # Save and show the FT data vs freq
pp.savefig(FilterFuncPSD())                     # Save and show the interpolated PSD vs freq
pp.savefig(FilterFreqPlot())                    # Save and show the Filtered FT data vs freq
pp.savefig(FilterTrunFreqPlot())                # Save and show the Filtered/Truncated FT data vs freq
pp.savefig(FilterTrunTimePlot())                # Save and show the InvFT data vs time
pp.savefig(FilterTrunTimePlot1to31())           # Save and show the InvFT data vs time (t:1~31s)
pp.savefig(FilterTrunTimePlotSTD())             # Save and show the InvFT data vs time (y-unit: StD)
pp.savefig(FilterTrunTimePlotSTD1to31())        # Save and show the InvFT data vs time (y-unit: StD) (t:1~31s)
pp.savefig(FinalPlot())
pp.savefig(FinalComparePlot())
pp.close()