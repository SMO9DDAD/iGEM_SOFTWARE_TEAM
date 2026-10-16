import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit


def model_quadratic(x,a,b,c):
    y = np.zeros(len(x),float)
    for i in range(len(x)):
        y[i] = (a*x[i]+b)*x[i]+c
    return y


def ProcessSpectrum(filename):#,WL_min,WL_max,f_ini,f_fin):

    #load lines
    with open(filename,'r') as file:
        lines = file.readlines()

    #skip header and load wavelength and reflection
    WL0 = np.zeros(len(lines)-17-1,float)
    R0 = np.zeros(len(lines)-17-1,float)
 
    for i in range(17,len(lines)-1):
        aux = lines[i].split()
        WL0[i-17] = float(aux[0])
        R0[i-17]  = float(aux[1])

    #limit the range of the spectrum
    i_wlmin = 0
    i_wlmax = 0
    for i in range(len(WL0)):
        if(WL0[i]<WL_min):
            i_wlmin = i
        if(WL0[i]<WL_max):
            i_wlmax = i

    WL = WL0[i_wlmin:i_wlmax]
    R  = R0[i_wlmin:i_wlmax]

    #calculate DFT for the selected frequency range
    N_f = 100
    Rp = np.zeros(N_f,complex)
    f0 = np.zeros(N_f,complex)

    for j in range(N_f):
        f0[j] = f_ini+(f_fin-f_ini)/N_f*j
        for i in range(len(R)):
            Rp[j] += R[i]*np.exp(-I*2*np.pi/WL[i]*f0[j])

    #we analyze just the modulus of the DFT
    Rp_mod0 = np.zeros(len(Rp),float)
    Rp_mod0 = np.real(Rp*np.conj(Rp))

    # vertex of the data to be
    i_f_max = np.argmax(Rp_mod0)
    Rp_max = Rp_mod0[i_f_max]
    f_max = f0[i_f_max]

    # subset of Rp_mod0 that will be fitted to quadratic function
    # middle pint: vertex, the same number
    N_pad = 10
    #Rp_mod = np.zeros(N_pad*2+1,float)
    #f1 = np.zeros(N_pad*2+1,float)
    #for i in range(N_subset):
    #    Rp_mod[i] = Rp_mod0[i+i_f_max-N_pad]
    #    f1[i] = f[i+i_f_max-N_subset/2]

    Rp_mod1 = Rp_mod0[i_f_max-N_pad:i_f_max+N_pad]
    f1 = f0[i_f_max-N_pad:i_f_max+N_pad]
    
    # second derivative
    a = (Rp_mod1[N_pad-1]+Rp_mod1[N_pad+1]-2*Rp_mod1[N_pad])/(f1[1]-f1[0])/2

    p1_guess = np.zeros(3,float)
    p1_guess[0] = a
    p1_guess[1] = -2*a
    p1_guess[2] = Rp_max+a*f_max**2
        
    #aux = np.log10(np.real(Rp*np.conj(Rp)))
    
    #print(p1_guess)

    popt,pcov = curve_fit(model_quadratic,f1,Rp_mod1,p0=p1_guess)

    #print(popt)

    #aux1 = Rp_mod1
    #aux2 = (popt[0]*f1+popt[1])*f1+popt[2]

    #fig1, ax1 = plt.subplots()
    #ax1.plot(aux1,'.')
    #ax1.plot(aux2,'.')
    #plt.show()

    return -popt[1]/(2*popt[0])


WL_min = 400  #nm
WL_max = 1000 #nm

f_ini = 7500
f_fin = 10000

I = np.complex128(0,1)

folder = '/Users/pepferre/Library/CloudStorage/OneDrive-URV/Prog/iGem/IGEM/20261006/'
subfolder = 'S02VOC/'
list_file = 'S02VOC_List.txt'

filename0 = folder+list_file

#load lines
with open(filename0,'r') as file:
    lines = file.readlines()

#print(lines)

N_meas = len(lines)/10
N_meas = 10
df_EOT = pd.DataFrame(np.zeros(len(lines),float))

df_EOTavg = pd.DataFrame(np.zeros(N_meas,float))
#df_EOTavg = pd.DataFrame()
df_EOTavg[0] = float(0)
df_EOTavg[1] = float(0)
df_EOTavg[2] = float(0)


N_rep = 10
EOT_aux = np.zeros(N_rep,float)
EOT_avg = np.zeros(N_meas,float)
EOT_std = np.zeros(N_meas,float)

i = 0
j = 0
for file in lines:    
    if(i<N_rep):
        print(j)
        print(i)
        filename0 = folder+subfolder+file[0:-1]
        print(filename0)
        EOT_aux[i] = ProcessSpectrum(filename0)
        df_EOT.iloc[i+j*N_rep] = EOT_aux[i]
        i += 1
    if(i==N_rep):
        EOT_avg[j] = np.average(EOT_aux)
        EOT_std[j] = np.std(EOT_aux)
        df_EOTavg.iloc[j,0] = j
        df_EOTavg.iloc[j,1] = np.average(EOT_aux)
        df_EOTavg.iloc[j,2] = np.std(EOT_aux)
        j += 1
        i = 0

print(EOT_avg,EOT_std)

df_EOT.to_excel(folder+file[0:-10]+"_EOT.xlsx", sheet_name='SheetName', index=False, header=None)
df_EOTavg.to_excel(folder+file[0:-10]+"_EOTavg.xlsx", sheet_name='SheetName', index=False, header=None)