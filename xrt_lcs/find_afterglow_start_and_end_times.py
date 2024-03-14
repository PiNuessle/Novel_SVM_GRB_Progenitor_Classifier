# %%
##package imports##
import pandas as pd
import math
import pandas as pd
import pdb
import os
import numpy as np
from astropy.io import ascii
import julian
from astropy.time import Time
import scipy.stats as st
import re
import datetime
from statistics import NormalDist
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pdb
import scipy.integrate as integrate
from scipy.optimize import fsolve
import gc
from astropy.time import Time
directory_path = os.getcwd()
os.chdir("..")
directory_path = os.getcwd()

# %%
##function definitions##
def friendly_loop(data): #This defines the start and end time of XRT data using the rates.
    i=0
    j=1
    if len(np.where(data['BGrate']!=0)[0])!=0:
        for k in range(0, len(data)):
            if data['BGrate'][i] == 0:
                if i == 0:
                    i=i+1
                elif data['BGrate'][i-1] == 0:
                    i=i+1
                else:
                    continue
            else:
                start_time=data["!Time"][i]+data['T_-ve'][i]
        for k in range(0, len(data)):
            if data['BGrate'][-j] == 0:
                if j==1:
                    j=j+1
                elif data['BGrate'][-(j-1)] == 0:
                    j=j+1
                else:
                    continue
            else:
                stop_time=data["!Time"][-j]+data['T_+ve'][-j]
    else:
        start_time=data["!Time"][0]+data['T_-ve'][0]
        stop_time=data["!Time"][-1]+data['T_+ve'][-1]
    return start_time, stop_time

#these next few functions are just power laws with different numbers o fbreaks, as necessitated by the fits in the XRT catalog 
def power_law (variable, index):
    value=variable**index
    return value

def singly_broken_PL(amplitude, break1, power1, power2, t0, t):
    if t<break1 and t>t0:
        value=amplitude*t**(power1)
    elif t>=break1:
        value=amplitude*t**(power2)*break1**(power1-power2)
    else:
        value=np.NaN
    return value

def doubly_broken_PL(amplitude, break1, break2, power1, power2, power3, t0, t):
    if t<break1 and t>t0:
        value=amplitude*t**(power1)
    elif t>=break1 and t<break2:
        value=amplitude*t**(power2)*break1**(power1-power2)
    elif t>=break2:
        value=amplitude*t**(power3)*break1**(power1-power2)*break2**(power2-power3)
    else:
        value=np.NaN
    return value

def triply_broken_PL(amplitude, break1, break2, break3, power1, power2, power3, power4, t0,\
                     t):
    if t<break1 and t>t0:
        value=amplitude*t**(power1)
    elif t>=break1 and t<break2:
        value=amplitude*t**(power2)*break1**(power1-power2)
    elif t>=break2 and t<break3:
        value=amplitude*t**(power3)*break1**(power1-power2)*break2**(power2-power3)
    elif t>=break3:
        value=amplitude*t**(power4)*break1**(power1-power2)*break2**(power2-power3)*\
                              break3**(power3-power4)
    else:
        value=np.NaN
    return value

def quadruply_broken_PL(amplitude, break1, break2, break3, break4, power1, power2, power3,\
                        power4, power5, t0, t):
    if t<break1 and t>t0:
        value=amplitude*t**(power1)
    elif t>=break1 and t<break2:
        value=amplitude*t**(power2)*break1**(power1-power2)
    elif t>=break2 and t<break3:
        value=amplitude*t**(power3)*break1**(power1-power2)*break2**(power2-power3)
    elif t>=break3 and t<break4:
        value=amplitude*t**(power4)*break1**(power1-power2)*break2**(power2-power3)*\
                              break3**(power3-power4)
    elif t>=break4:
        value=amplitude*t**(power5)*break1**(power1-power2)*break2**(power2-power3)*\
                              break3**(power3-power4)*break4**(power4-power5)
    else:
        value=np.NaN
    return value

def quintuply_broken_PL(amplitude, break1, break2, break3, break4, break5, power1, power2,\
                        power3, power4, power5, power6, t0, t):
    if t<break1 and t>t0:
        value=amplitude*t**(power1)
    elif t>=break1 and t<break2:
        value=amplitude*t**(power2)*break1**(power1-power2)
    elif t>=break2 and t<break3:
        value=amplitude*t**(power3)*break1**(power1-power2)*break2**(power2-power3)
    elif t>=break3 and t<break4:
        value=amplitude*t**(power4)*break1**(power1-power2)*break2**(power2-power3)*\
                              break3**(power3-power4)
    elif t>=break4 and t<break5:
        value=amplitude*t**(power5)*break1**(power1-power2)*break2**(power2-power3)*\
                              break3**(power3-power4)*break4**(power4-power5)
    elif t>=break5:
        value=amplitude*t**(power6)*break1**(power1-power2)*break2**(power2-power3)*\
                              break3**(power3-power4)*break4**(power4-power5)*\
                                break5**(power5-power6)
    else:
        value=np.NaN
    return value

#Then, we needed the power law derivatives to do error propogation for the AG fluence

def PL_derivative (variable, index):
    value=index*(variable**(index-1))
    return value

def singly_broken_PL_derivative(amplitude, break1, power1, power2, t0, t):
    if t<break1 and t>t0:
        value=amplitude*power1*t**(power1-1)
    elif t>=break1:
        value=amplitude*power2*t**(power2-1)*break1**(power1-power2)
    else:
        value=np.NaN
    return value

def doubly_broken_PL_derivative(amplitude, break1, break2, power1, power2, power3, t0, t):
    if t<break1 and t>t0:
        value=amplitude*power1*t**(power1-1)
    elif t>=break1 and t<break2:
        value=amplitude*power2*t**(power2-1)*break1**(power1-power2)
    elif t>=break2:
        value=amplitude*power3*t**(power3-1)*break1**(power1-power2)*\
        break2**(power2-power3)
    else:
        value=np.NaN
    return value

def triply_broken_PL_derivative(amplitude, break1, break2, break3, power1, power2,\
                                power3, power4, t0, t):
    if t<break1 and t>t0:
        value=amplitude*power1*t**(power1-1)
    elif t>=break1 and t<break2:
        value=amplitude*power2*t**(power2-1)*break1**(power1-power2)
    elif t>=break2 and t<break3:
        value=amplitude*power3*t**(power3-1)*break1**(power1-power2)\
        *break2**(power2-power3)
    elif t>=break3:
        value=amplitude*power4*t**(power4-1)*break1**(power1-power2)*break2**(power2-power3)\
        *break3**(power3-power4)
    else:
        value=np.NaN
    return value

def quadruply_broken_PL_derivative(amplitude, break1, break2, break3, break4, power1,\
                                   power2, power3, power4, power5, t0, t):
    if t<break1 and t>t0:
        value=amplitude*power3*t**(power1-1)
    elif t>=break1 and t<break2:
        value=amplitude*power2*t**(power2-1)*break1**(power1-power2)
    elif t>=break2 and t<break3:
        value=amplitude*power3*t**(power3-1)*break1**(power1-power2)*break2**(power2-power3)
    elif t>=break3 and t<break4:
        value=amplitude*power4*t**(power4-1)*break1**(power1-power2)*break2**(power2-power3)\
        *break3**(power3-power4)
    elif t>=break4:
        value=amplitude*power5*t**(power5-1)*break1**(power1-power2)*break2**(power2-power3)\
        *break3**(power3-power4)*break4**(power4-power5)
    else:
        value=np.NaN
    return value

def quintuply_broken_PL_derivative(amplitude, break1, break2, break3, break4, break5, \
                                   power1, power2, power3, power4, power5, power6, t0, t):
    if t<break1 and t>t0:
        value=amplitude*power1*t**(power1-1)
    elif t>=break1 and t<break2:
        value=amplitude*power2*t**(power2-1)*break1**(power1-power2)
    elif t>=break2 and t<break3:
        value=amplitude*power3*t**(power3-1)*break1**(power1-power2)*break2**(power2-power3)
    elif t>=break3 and t<break4:
        value=amplitude*power4*t**(power4-1)*break1**(power1-power2)*break2**(power2-power3)\
        *break3**(power3-power4)
    elif t>=break4 and t<break5:
        value=amplitude*power5*t**(power5-1)*break1**(power1-power2)*break2**(power2-power3)\
        *break3**(power3-power4)*break4**(power4-power5)
    elif t>=break5:
        value=amplitude*power6*t**(power6-1)*break1**(power1-power2)*break2**(power2-power3)\
        *break3**(power3-power4)*break4**(power4-power5)*break5**(power5-power6)
    else:
        value=np.NaN
    return value

#I define the flux here to depend on where you are in the power law, and how many breaks it has

def defining_the_flux(tarjit, temporal_indices):
    flux_catch=' Obs Flux_{} (pc) '.format(tarjit+1)
    tarjete=tarjit+1
    xrt_flux=more_xrt_data[flux_catch][xrt_match]
    state=0
    if xrt_flux==' N/A ' or xrt_flux=='N/A' or xrt_flux=='NaN':
        flux_catch=' Obs Flux_{} (wt) '.format(tarjit+1)
        xrt_flux=more_xrt_data[flux_catch][xrt_match]
        state=1
        #weird case where sometimes XRT is still repointing or something, IDK.
        #It nets us 220101A, which is what I wanted.
    if tarjit==0:
        if xrt_flux==' N/A ' or xrt_flux=='N/A' or xrt_flux=='NaN':
            flux_catch=' Obs Flux_{} (pc) '.format(tarjit+2)
            xrt_flux=more_xrt_data[flux_catch][xrt_match]
            tarjete=tarjit+2
            state=0
            if xrt_flux==' N/A ' or xrt_flux=='N/A' or xrt_flux=='NaN':
                flux_catch=' Obs Flux_{} (wt) '.format(tarjit+2)
                xrt_flux=more_xrt_data[flux_catch][xrt_match]
                state=1
                #weird case where if it's at the very beginning, XRT
                #might not be measuring flux yet
    return flux_catch, xrt_flux, tarjete, state

#whereas the hardness ration depends on very little

def HR_err_func(S_top, S_bot, sig_S_top, sig_S_bot):
    HR_err=np.sqrt(np.power((sig_S_top/S_bot), 2)+np.power((S_top*sig_S_bot)/\
                                                          np.power(S_bot,2), 2))
    return HR_err

#unless you can actually *derive* a well-approximated spectrum like GBM can, at which point you might as well just use the peak energy as the main analysis does.

def Fermi_HR_func(relevant_fermi_data, spectral_model, entry):
    low_e_range=[15, 150]
    high_e_range=[150, 1500]
    if spectral_model=='flnc_plaw              ' or spectral_model==1:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_plaw_index'])
        low_HR=integrate.quad(lambda x: power_law(x, p_1), low_e_range[0], \
                              low_e_range[1]) #PL
        high_HR=integrate.quad(lambda x: power_law(x, p_1), high_e_range[0], \
                              high_e_range[1]) #PL
        HR=high_HR[0]/low_HR[0]
        
    elif spectral_model=='flnc_comp              ' or spectral_model==2:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_comp_index'])
        E_break=float(relevant_fermi_data.at[entry, 'flnc_comp_epeak'])
        low_HR=integrate.quad(lambda x:  Compton_PL(x, p_1, E_break), \
                       low_e_range[0], low_e_range[1]) #CPL
        high_HR=integrate.quad(lambda x: Compton_PL(x, p_1, E_break), \
                       high_e_range[0], high_e_range[1]) #CPL
        HR=high_HR[0]/low_HR[0]
        
    elif spectral_model=='flnc_band              ' or spectral_model==3:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_band_alpha'])
        p_2=float(relevant_fermi_data.at[entry, 'flnc_band_beta'])
        E_break=float(relevant_fermi_data.at[entry, 'flnc_band_epeak'])
        low_HR=integrate.quad(lambda x: Band_function(x, E_break, p_1, p_2), \
                       low_e_range[0], low_e_range[1]) #Band
        high_HR=integrate.quad(lambda x: Band_function(x, E_break, p_1, p_2), \
                       high_e_range[0], high_e_range[1]) #Band
        HR=high_HR[0]/low_HR[0]
        
    elif spectral_model=='flnc_sbpl              ' or spectral_model==4:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_sbpl_indx1'])
        p_2=float(relevant_fermi_data.at[entry, 'flnc_sbpl_indx2'])
        E_break=float(relevant_fermi_data.at[entry, 'flnc_sbpl_brken'])
        smoothen=float(relevant_fermi_data.at[entry, 'flnc_sbpl_brksc'])
        low_HR=integrate.quad(lambda x: Smoothly_Broken_PL(x, E_break, p_1, p_2, smoothen), \
                       low_e_range[0], low_e_range[1])
        high_HR=integrate.quad(lambda x: Smoothly_Broken_PL(x, E_break, p_1, p_2, smoothen),\
                       high_e_range[0], high_e_range[1])
        HR=high_HR[0]/low_HR[0]
        
    else:
        if relevant_fermi_data.at[entry, 'flnc_plaw_index']:
            p_1=float(relevant_fermi_data.at[entry, 'flnc_plaw_index'])
            low_HR=integrate.quad(lambda x: power_law(x, p_1), low_e_range[0], \
                                  low_e_range[1])
            high_HR=integrate.quad(lambda x: power_law(x, p_1), \
                               high_e_range[0], high_e_range[1]) 
            #PL because that's the simplest one, unfortunately
            HR=high_HR[0]/low_HR[0]
        else:
            HR=np.NA
    return HR

#We had to approximate an error in these values as well, which we took to be related to the fluence

def swift_error_calculator(entry, fluence_table):
    if fluence_table.at[entry, ' 50_100kev_low '] != ' N/A ':
         HR_err_left = HR_err_func(float(fluence_table.at[entry, ' 50_100kev ']), \
                                float(fluence_table.at[entry, ' 25_50kev ']), \
                                float(fluence_table.at[entry,' 50_100kev_low ']), \
                                float(fluence_table.at[entry,' 25_50kev_low ']))
#         HR_err_left = (float(fluence_table.at[entry,' 50_100kev '])-\
#                         float(fluence_table.at[entry, ' 50_100kev_low ']))/\
#                             (float(fluence_table.at[entry, ' 25_50kev '])-\
#                                    float(fluence_table.at[entry,' 25_50kev_low ']))
    else:
        HR_err_left = 0
    if fluence_table.at[entry, ' 50_100kev_hi '] != ' N/A ':
         HR_err_right=HR_err_func(float(fluence_table.at[entry, ' 50_100kev ']),\
                                 float(fluence_table.at[entry, ' 25_50kev ']), \
                                 float(fluence_table.at[entry,' 50_100kev_hi ']), \
                             float(fluence_table.at[entry,' 25_50kev_hi ']))
#         HR_err_right = (float(fluence_table.at[entry, ' 50_100kev '])-\
#                         float(fluence_table.at[entry, ' 50_100kev_hi ']))/\
#                             (float(fluence_table.at[entry, ' 25_50kev '])-\
#                                    float(fluence_table.at[entry, ' 25_50kev_hi ']))
    else:
        HR_err_right = 0
    return [HR_err_left, HR_err_right]

#we also had to define all the spectral functions, so the program wouldn't get confused

def Compton_PL (variable, index, break_E):
    value=np.power(variable, index)*np.exp(-variable/break_E)
    return value

def Band_function (variable, break_E, index_1, index_2):
    index_1=abs(index_1)
    index_2=abs(index_2)
    value=np.where(variable < ((index_1-index_2)*break_E)/(index_1), \
        ((variable/100)**index_2)*np.exp(index_2-index_1)*(((index_1-index_2)*\
        break_E)/(100*(index_1)))**(index_1-index_2), ((variable/100)**index_1)*\
        np.exp(-((index_1)*variable)/break_E))
    return np.real(value)

def Smoothly_Broken_PL (variable, break_E, index_1, index_2, smoothen):
    value=(variable/break_E)**(index_1)*(0.5*(1+(variable/break_E)**(1/smoothen)))**\
    (-(index_1-index_2)/smoothen)
    return value

#I don't think we actually used this one though

def z_power_law(variable, index, redshift):
    value=(variable*(1+redshift))**(index)
    return value

#Swift has difficulty constraining a burst's peak energy, so we can only use two functions to estimate its prompt fluence
#we also have to read it from a separate table, as BAT doesn't store the PL and CPL flux and fluence data in the same table, but rather in four different ones
def swift_fluence_function(burst_index, column_name):
    if swift_fluences_reference.at[burst_index, ' Best-fit model']==' PL':
        fluence=swift_fluences_PL.at[burst_index, column_name]
    elif swift_fluences_reference.at[burst_index, ' Best-fit model']==' CPL':
        fluence=swift_fluences_CPL.at[burst_index, column_name]
    else:
        fluence=np.NaN
    return fluence

def swift_flux_function(burst_index, column_name):
    if swift_fluences_reference.at[burst_index, ' Best-fit model']==' PL':
        flux=swift_fluxes_PL.at[burst_index, column_name]
    elif swift_fluences_reference.at[burst_index, ' Best-fit model']==' CPL':
        flux=swift_fluxes_CPL.at[burst_index, column_name]
    else:
        flux=np.NaN
    return flux

#So we also had to figure out which model Swift-BAT fit to each burst
def swift_model_function(burst_index, column_name):
    if swift_fluences_reference.at[burst_index, ' Best-fit model']==' PL':
        model=swift_model_PL.at[burst_index, column_name]
    elif swift_fluences_reference.at[burst_index, ' Best-fit model']==' CPL':
        model=swift_model_CPL.at[burst_index, column_name]
    else:
        model=np.NaN
    return model

#whereas we can read of the peak energy of every Fermi-GBM burst that isn't fit by a power law from one table
def getting_the_GBM_E_peak(relevant_fermi_data, spectral_model, entry):
    if spectral_model=='flnc_plaw              ' or spectral_model==1:
        if relevant_fermi_data.at[entry, 'flnc_comp_epeak']:
            E_peak=float(relevant_fermi_data.at[entry, 'flnc_comp_epeak'])
        elif relevant_fermi_data.at[entry, 'flnc_band_epeak']:
            E_peak=float(relevant_fermi_data.at[entry, 'flnc_band_epeak'])
        else:
            E_peak=np.NaN
    elif spectral_model=='flnc_comp              ' or spectral_model==2:
        E_peak=float(relevant_fermi_data.at[entry, 'flnc_comp_epeak'])
    elif spectral_model=='flnc_band              ' or spectral_model==3:
        E_peak=float(relevant_fermi_data.at[entry, 'flnc_band_epeak'])
    elif spectral_model=='flnc_sbpl              ' or spectral_model==4:
        E_peak=float(relevant_fermi_data.at[entry, 'flnc_sbpl_brken'])
    else:
        E_peak=np.NaN
    return E_peak

#and they immediately come with an uncertianty that we don't have to derive ourselves
def getting_the_GBM_E_peak_err(relevant_fermi_data, spectral_model, entry):
    if spectral_model=='flnc_plaw              ' or spectral_model==1:
        if relevant_fermi_data.at[entry, 'flnc_comp_epeak']:
            E_peak_err=(abs(float(relevant_fermi_data.at[entry, 'flnc_comp_epeak_pos_err']))\
                +abs(float(relevant_fermi_data.at[entry, 'flnc_comp_epeak_neg_err'])))/2
        elif relevant_fermi_data.at[entry, 'flnc_band_epeak']:
            E_peak_err=(abs(float(relevant_fermi_data.at[entry, 'flnc_band_epeak_pos_err']))\
                +abs(float(relevant_fermi_data.at[entry, 'flnc_band_epeak_neg_err'])))/2
            #because np.mean SUDDENLY doesn't work on floats. >8(
        else:
            E_peak_err=np.NaN
    elif spectral_model=='flnc_comp              ' or spectral_model==2:
        E_peak_err=(abs(float(relevant_fermi_data.at[entry, 'flnc_comp_epeak_pos_err']))\
                +abs(float(relevant_fermi_data.at[entry, 'flnc_comp_epeak_neg_err'])))/2
    elif spectral_model=='flnc_band              ' or spectral_model==3:
        E_peak_err=(abs(float(relevant_fermi_data.at[entry, 'flnc_band_epeak_pos_err']))\
                +abs(float(relevant_fermi_data.at[entry, 'flnc_band_epeak_neg_err'])))/2
    elif spectral_model=='flnc_sbpl              ' or spectral_model==4:
        E_peak_err=(abs(float(relevant_fermi_data.at[entry, 'flnc_sbpl_brken_pos_err']))\
                +abs(float(relevant_fermi_data.at[entry, 'flnc_sbpl_brken_neg_err'])))/2
    else:
        E_peak_err=np.NaN
    return E_peak_err

#less than half of swift bursts are fit with a CPL, so this doesn't work as well for them
def getting_the_BAT_E_peak(entry):
    if swift_fluences_reference.at[entry, ' Best-fit model']==' CPL':
        if swift_model_CPL.at[entry, ' Epeak '] != ' N/A ':
            E_peak=float(swift_model_CPL.at[entry, ' Epeak '])
        else:
            E_peak= np.NaN
    elif swift_fluences_reference.at[entry, ' Best-fit model']==' PL':
        if swift_model_CPL.at[entry, ' Epeak ']:
            if swift_model_CPL.at[entry, ' Epeak '] != ' N/A ':
                E_peak=float(swift_model_CPL.at[entry, ' Epeak '])
            else:
                E_peak= np.NaN
        else:
            E_peak= np.NaN
    else:
        E_peak=np.NaN
    return E_peak

#we cna also use those functions from fermi to derive a good k-correction
def Fermi_k_func(relevant_fermi_data, spectral_model, entry, redshift):
#     pdb.set_trace()
    low_e_range=[10, 1000]
    high_e_range=[1/(1+redshift), 10000/(1+redshift)]
    if spectral_model=='flnc_plaw              ' or spectral_model==1:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_plaw_index'])
        denominator=integrate.quad(lambda x: power_law(x, p_1), low_e_range[0], \
                              low_e_range[1]) #PL
        numerator=integrate.quad(lambda x: power_law(x, p_1), high_e_range[0], \
                              high_e_range[1]) #PL
        k=numerator[0]/denominator[0]
        
    elif spectral_model=='flnc_comp              ' or spectral_model==2:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_comp_index'])
        E_break=float(relevant_fermi_data.at[entry, 'flnc_comp_epeak'])
        denominator=integrate.quad(lambda x:  Compton_PL(x, p_1, E_break), \
                       low_e_range[0], low_e_range[1]) #CPL
        numerator=integrate.quad(lambda x: Compton_PL(x, p_1, E_break), \
                       high_e_range[0], high_e_range[1]) #CPL
        k=numerator[0]/denominator[0]
        
    elif spectral_model=='flnc_band              ' or spectral_model==3:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_band_alpha'])
        p_2=float(relevant_fermi_data.at[entry, 'flnc_band_beta'])
        E_break=float(relevant_fermi_data.at[entry, 'flnc_band_epeak'])
        denominator=integrate.quad(lambda x: Band_function(x, E_break, p_1, p_2), \
                       low_e_range[0], low_e_range[1]) #Band
        numerator=integrate.quad(lambda x: Band_function(x, E_break, p_1, p_2), \
                       high_e_range[0], high_e_range[1]) #Band
        k=numerator[0]/denominator[0]
        
    elif spectral_model=='flnc_sbpl              ' or spectral_model==4:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_sbpl_indx1'])
        p_2=float(relevant_fermi_data.at[entry, 'flnc_sbpl_indx2'])
        E_break=float(relevant_fermi_data.at[entry, 'flnc_sbpl_brken'])
        smoothen=float(relevant_fermi_data.at[entry, 'flnc_sbpl_brksc'])
        denominator=integrate.quad(lambda x: Smoothly_Broken_PL(x, E_break, p_1, p_2, smoothen), \
                       low_e_range[0], low_e_range[1])
        numerator=integrate.quad(lambda x: Smoothly_Broken_PL(x, E_break, p_1, p_2, smoothen),\
                       high_e_range[0], high_e_range[1])
        k=numerator[0]/denominator[0]
        
    else:
        if relevant_fermi_data.at[entry, 'flnc_plaw_index']:
            p_1=float(relevant_fermi_data.at[entry, 'flnc_plaw_index'])
            denominator=integrate.quad(lambda x: power_law(x, p_1), low_e_range[0], \
                                  low_e_range[1])
            numerator=integrate.quad(lambda x: power_law(x, p_1), \
                               high_e_range[0], high_e_range[1]) 
            #PL because that's the simplest one, unfortunately
            k=numerator[0]/denominator[0]
        else:
            k=np.NA
    return k

#but here, XRT can fill in a lot of the gap from BAT
def XRT_k_func(relevant_xrt_data, entry, redshift):
#     pdb.set_trace()
    low_e_range=[0.2, 10]
    high_e_range=[0.2/(1+redshift), 10/(1+redshift)]
    if float(relevant_xrt_data.at[entry, 'plaw_index']):
        p_1=-abs(float(relevant_xrt_data.at[entry, 'plaw_index']))
        denominator=integrate.quad(lambda x: power_law(x, p_1), low_e_range[0], \
                              low_e_range[1]) #PL
        numerator=integrate.quad(lambda x: power_law(x, p_1), high_e_range[0], \
                              high_e_range[1]) #PL
        k=np.real(numerator[0]/denominator[0])
    else:
        k=np.NaN
    return k

#I don't think we need the kinetic energy--I think I jsut forgot to delete it
def kinetic_energy(nFn, beta, D_L, redshift, time, nu):
#     pdb.set_trace()
    time=((time/60)/60)/24
    nu=nu/(10**18)
    p=2*beta
    Y=1
    eps_e=0.1
    eps_b=0.01
    nFn_term=(nFn/(5.2*1e-14))**(4/(p+2))
    dist_term=((D_L*1e-28)**(8/(p+2)))
    rs_term=(1+redshift)**(-1)
    time_term=(time**((3*p-2)/(p+2)))*(1+Y)**(4/(p+2))
    if p != 2:
        fluxy_term=((6.73*((p-2)/(p-1))**(p-1))*((3.3*1e-6)**((p-2.3)/2)))**(-4/(p+2))
    else:
        #This is jsut a guess around the value given the other results, 
        #since p=2 is a discountinuity
        fluxy_term=((6.73*((p-1.9925)/(p-1))**(p-1))*((3.3*1e-6)**((p-2.3)/2)))**(-4/(p+2))
    carrier_dens_terms=((eps_e*1e1)**((4*(1-p))/(p+2)))*(eps_b*1e2)**((2-p)/(p+2))
    freq_term=(nu*1e-18)**((2*(p-2))/(p+2))
    E_k=nFn_term*dist_term*rs_term*time_term*fluxy_term*carrier_dens_terms*freq_term
    return E_k

#same thing with efficency
def efficiency(E_k, E_g):
    eta=E_g/(E_g+E_k)
    return eta

#again, we won't need this
def prompt_Eiso_error(d_lum, flu, relevant_fermi_data, spectral_model, entry, redshift):
    k=Fermi_k_func(sample_2_data, ghostie, fermi_placement, rs)
    S_err=(4*np.pi*((d_lum)**2)*1*k/(1+rs)*relevant_fermi_data[fermi_placement, \
                                                              'fluence_error'])**2
    low_e_range=[10, 1000]
    high_e_range=[1/(1+redshift), 10000/(1+redshift)]
    if spectral_model=='flnc_plaw              ' or spectral_model==1:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_plaw_index'])
        denominator=integrate.quad(lambda x: power_law(x, p_1), low_e_range[0], \
                              low_e_range[1]) #PL
        numerator=integrate.quad(lambda x: power_law(x, p_1), high_e_range[0], \
                              high_e_range[1]) #PL
        q=numerator[0]/denominator[0]
        
    elif spectral_model=='flnc_comp              ' or spectral_model==2:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_comp_index'])
        E_break=float(relevant_fermi_data.at[entry, 'flnc_comp_epeak'])
        denominator=integrate.quad(lambda x:  Compton_PL(x, p_1, E_break), \
                       low_e_range[0], low_e_range[1]) #CPL
        numerator=integrate.quad(lambda x: Compton_PL(x, p_1, E_break), \
                       high_e_range[0], high_e_range[1]) #CPL
        q=numerator[0]/denominator[0]
        
    elif spectral_model=='flnc_band              ' or spectral_model==3:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_band_alpha'])
        p_2=float(relevant_fermi_data.at[entry, 'flnc_band_beta'])
        E_break=float(relevant_fermi_data.at[entry, 'flnc_band_epeak'])
        denominator=integrate.quad(lambda x: Band_function(x, E_break, p_1, p_2), \
                       low_e_range[0], low_e_range[1]) #Band
        numerator=integrate.quad(lambda x: Band_function(x, E_break, p_1, p_2), \
                       high_e_range[0], high_e_range[1]) #Band
        q=numerator[0]/denominator[0]
        
    elif spectral_model=='flnc_sbpl              ' or spectral_model==4:
        p_1=float(relevant_fermi_data.at[entry, 'flnc_sbpl_indx1'])
        p_2=float(relevant_fermi_data.at[entry, 'flnc_sbpl_indx2'])
        E_break=float(relevant_fermi_data.at[entry, 'flnc_sbpl_brken'])
        smoothen=float(relevant_fermi_data.at[entry, 'flnc_sbpl_brksc'])
        denominator=integrate.quad(lambda x: Smoothly_Broken_PL(x, E_break, p_1, p_2, smoothen), \
                       low_e_range[0], low_e_range[1])
        numerator=integrate.quad(lambda x: Smoothly_Broken_PL(x, E_break, p_1, p_2, smoothen),\
                       high_e_range[0], high_e_range[1])
        q=numerator[0]/denominator[0]
        
    else:
        if relevant_fermi_data.at[entry, 'flnc_plaw_index']:
            p_1=float(relevant_fermi_data.at[entry, 'flnc_plaw_index'])
            denominator=integrate.quad(lambda x: power_law(x, p_1), low_e_range[0], \
                                  low_e_range[1])
            numerator=integrate.quad(lambda x: power_law(x, p_1), \
                               high_e_range[0], high_e_range[1]) 
            #PL because that's the simplest one, unfortunately
            k=numerator[0]/denominator[0]
        else:
            q=np.NA
    return q

#but this is where it gets interesting--we're taking an fluence integral under the entire afterglow segment, and normalizing it using the 11 and 24 hour fluxes.
#Is it particularly reliable? probably only for the bursts where the plateau and afterglow dominate over the steep decay.
def afterglow_fluence_finder(spectral_sample_data, ag_instrument_data, placement, \
                             front, back, error_front, error_back):
    row=spectral_sample_data['XRT row'][placement]
    model_selector=ag_instrument_data[' #breaks '][row]
    if model_selector>=0:
        power_1=-float(ag_instrument_data[' alpha_1 '][row])
        if model_selector>=1:
            break_1=int(float(ag_instrument_data[' break_1 '][row]))
            split_err=re.split(",", ag_instrument_data.at[placement, ' D_break_1 '].strip())
            if split_err and split_err[0] !='' and split_err[1] != '':
                split_err= list(map(float, split_err))
                err_break_1=np.mean([abs(split_err[0]), abs(split_err[1])])
            elif split_err[0] !='':
                err_break_1=float(split_err[0])
            elif split_err[1] != '':
                err_break_1=float(split_err[1])
            else:
                err_break_1=0
            power_2=-float(ag_instrument_data[' alpha_2 '][row])
            if model_selector>=2:
                break_2=int(float(ag_instrument_data[' break_2 '][row]))
                split_err=re.split(",", ag_instrument_data.at[placement, ' D_break_2 '].strip())
                if split_err and split_err[0] !='' and split_err[1] != '':
                    split_err= list(map(float, split_err))
                    err_break_2=np.mean([abs(split_err[0]), abs(split_err[1])])
                elif split_err[0] !='':
                    err_break_2=float(split_err[0])
                elif split_err[1] != '':
                    err_break_2=float(split_err[1])
                else:
                    err_break_2=0
                power_3=-float(ag_instrument_data[' alpha_3 '][row])
                if model_selector>=3:
                    break_3=int(float(ag_instrument_data[' break_3 '][row]))
                    split_err=re.split(",", ag_instrument_data.at[placement, ' D_break_3 '].strip())
                    if split_err and split_err[0] !='' and split_err[1] != '':
                        split_err= list(map(float, split_err))
                        err_break_3=np.mean([abs(split_err[0]), abs(split_err[1])])
                    elif split_err[0] !='':
                        err_break_3=float(split_err[0])
                    elif split_err[1] != '':
                        err_break_3=float(split_err[1])
                    else:
                        err_break_3=0
                    power_4=-float(ag_instrument_data[' alpha_4 '][row])
                    if model_selector>=4:
                        break_4=int(float(ag_instrument_data[' break_4 '][row]))
                        split_err=re.split(",", ag_instrument_data.at[placement, ' D_break_4 '].strip())
                        if split_err and split_err[0] !='' and split_err[1] != '':
                            split_err= list(map(float, split_err))
                            err_break_4=np.mean([abs(split_err[0]), abs(split_err[1])])
                        elif split_err[0] !='':
                            err_break_4=float(split_err[0])
                        elif split_err[1] != '':
                            err_break_4=float(split_err[1])
                        else:
                            err_break_4=0
                        power_5=-float(ag_instrument_data[' alpha_5 '][row])
                        if model_selector==5:
                            break_5=int(float(ag_instrument_data[' break_5 '][row]))
                            split_err=re.split(",", ag_instrument_data.at[placement, ' D_break_5 '].strip())
                            if split_err and split_err[0] !='' and split[1] != '':
                                split_err= list(map(float, split_err))
                                err_break_5=np.mean([abs(split_err[0]), abs(split_err[1])])
                            elif split_err[0] !='':
                                err_break_5=float(split_err[0])
                            elif split_err[1] != '':
                                err_break_5=float(split_err[1])
                            else:
                                err_break_5=0
                            power_6=-float(ag_instrument_data[' alpha_6 '][row]) 
                            initial_value=quintuply_broken_PL(1, break_1, break_2, break_3, \
                                        break_4, break_5, power_1, power_2, power_3, power_4, \
                                            power_5, power_6, 0, 11*60*60)
                            Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_11 '][row]),\
                             0.01)[0]
                            AG_fluence=integrate.quad(lambda t: quintuply_broken_PL(Q, \
                                break_1, break_2, break_3, break_4, break_5, power_1, \
                                power_2, power_3, power_4, power_5, power_6, 0, t), \
                                                      front, back)[0]
                            AG_fluence_err=np.sqrt((\
                                quintuply_broken_PL_derivative(Q, break_1, break_2, break_3,\
                                break_4, break_5, power_1, power_2, power_3, power_4, power_5, \
                                power_6, 0, front)*error_front)**2+(\
                                quintuply_broken_PL_derivative(Q, break_1, break_2, break_3,\
                                break_4, break_5, power_1, power_2, power_3, power_4, power_5, \
                                power_6, 0, back)*error_back)**2)
                        else:
                            initial_value=quadruply_broken_PL(1, break_1, break_2, break_3, break_4, \
                                                 power_1, power_2, power_3, power_4, power_5, \
                                                 0, 11*60*60)
                            Q=fsolve(lambda A: A*initial_value-float(\
                                        ag_instrument_data[' Flux_11 '][row]),0.01)[0]
                            AG_fluence=integrate.quad(lambda t: quadruply_broken_PL(Q, \
                                break_1, break_2, break_3, break_4, power_1, power_2, \
                                power_3, power_4, power_5, 0, t), front, back)[0]
                            AG_fluence_err=np.sqrt((quadruply_broken_PL_derivative(Q, \
                                break_1, break_2, break_3, break_4, power_1, power_2, \
                                power_3, power_4, power_5, 0, front)*error_front)**2+(\
                                quadruply_broken_PL_derivative(Q, \
                                break_1, break_2, break_3, break_4, power_1, power_2, \
                                power_3, power_4, power_5, 0, back)*error_back)**2)
                    else:
                        initial_value=triply_broken_PL(1, break_1, break_2, break_3, \
                                                 power_1, power_2, power_3, power_4, \
                                                 0, 11*60*60)
                        Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_11 '][row]),\
                         0.01)[0]
                        AG_fluence=integrate.quad(lambda t: triply_broken_PL(Q, \
                                break_1, break_2, break_3, power_1, power_2, \
                                power_3, power_4, 0, t), front, back)[0]
                        AG_fluence_err=np.sqrt((triply_broken_PL_derivative(Q, \
                                break_1, break_2, break_3, power_1, power_2, \
                                power_3, power_4, 0, front)*error_front)**2+(\
                                triply_broken_PL_derivative(Q, \
                                break_1, break_2, break_3, power_1, power_2, \
                                power_3, power_4, 0, back)*error_back)**2)
                else:
                    initial_value=doubly_broken_PL(1, break_1, break_2, \
                                                 power_1, power_2, power_3, \
                                                 0, 11*60*60)
                    Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_11 '][row]),\
                         0.01)[0]
                    AG_fluence=integrate.quad(lambda t: doubly_broken_PL(Q, \
                                break_1, break_2, power_1, power_2, \
                                power_3, 0, t), front, back)[0]
                    AG_fluence_err=np.sqrt((doubly_broken_PL_derivative(Q, \
                                break_1, break_2, power_1, power_2, \
                                power_3, 0, front)*error_front)**2+(\
                                doubly_broken_PL_derivative(Q, \
                                break_1, break_2, power_1, power_2, \
                                power_3, 0, back)*error_back)**2)
            else:
                initial_value=singly_broken_PL(1, break_1, power_1, power_2, 0, 11*60*60)
                Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_11 '][row]),\
                         0.01)[0]
                AG_fluence=integrate.quad(lambda t: singly_broken_PL(Q, \
                                break_1, power_1, power_2, 0, t), front, back)[0]
                AG_fluence_err=np.sqrt(\
                                (singly_broken_PL_derivative(1, break_1, power_1, power_2,\
                                0, front)*error_front**2\
                                 +(singly_broken_PL_derivative(1, break_1, power_1, \
                                power_2, 0, back)*error_back)**2))
        else:
            initial_value=power_law(power_1, 11*60*60)
            if ag_instrument_data[' Flux_11 '][row] != ' N/A ':
                Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_11 '][row]),\
                         0.01)[0]
            elif ag_instrument_data[' Flux_24 '][row] != " N/A ":
                Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_24 '][row]),\
                         0.01)[0]
            else:
                Q=np.NaN
            AG_fluence=integrate.quad(lambda t: Q*power_law(power_1, t), front, back)[0]
            AG_fluence_err=np.sqrt((PL_derivative(power_1, front)*error_front)**2\
                             +(PL_derivative(power_1, back)*error_back)**2)
    else:
        AG_fluence=np.NaN
    # print('The time starts at: {}s.'.format(front))
    # print('The time ends at: {}s.'.format(back))
    # print('There are {} breaks in the model.'.format(model_selector))
    # print('The normalization factor is {}.'.format(Q))
    # print('This corresponds to row {} in the XRT data.'.format(row))
    # print('That gives us a fluence of {} photons/cm^2'.format(AG_fluence))
    return AG_fluence

#and here's why we had to find all those errors earlier
def afterglow_fluence_err_finder(spectral_sample_data, ag_instrument_data, placement, \
                                front, back, error_front, error_back):
    row=spectral_sample_data['XRT row'][placement]
    model_selector=ag_instrument_data[' #breaks '][row]
    if model_selector>=0:
        power_1=-float(ag_instrument_data[' alpha_1 '][row])
        if model_selector>=1:
            break_1=int(float(ag_instrument_data[' break_1 '][row]))
            split_err=re.split(",", ag_instrument_data.at[placement, ' D_break_1 '].strip())
            if split_err and split_err[0] !='' and split_err[1] != '':
                split_err= list(map(float, split_err))
                err_break_1=np.mean([abs(split_err[0]), abs(split_err[1])])
            elif split_err[0] !='':
                err_break_1=float(split_err[0])
            elif split_err[1] != '':
                err_break_1=float(split_err[1])
            else:
                err_break_1=0
            power_2=-float(ag_instrument_data[' alpha_2 '][row])
            if model_selector>=2:
                break_2=int(float(ag_instrument_data[' break_2 '][row]))
                split_err=re.split(",", ag_instrument_data.at[placement, ' D_break_2 '].strip())
                if split_err and split_err[0] !='' and split_err[1] != '':
                    split_err= list(map(float, split_err))
                    err_break_2=np.mean([abs(split_err[0]), abs(split_err[1])])
                elif split_err[0] !='':
                    err_break_2=float(split_err[0])
                elif split_err[1] != '':
                    err_break_2=float(split_err[1])
                else:
                    err_break_2=0
                power_3=-float(ag_instrument_data[' alpha_3 '][row])
                if model_selector>=3:
                    break_3=int(float(ag_instrument_data[' break_3 '][row]))
                    split_err=re.split(",", ag_instrument_data.at[placement, ' D_break_3 '].strip())
                    if split_err and split_err[0] !='' and split_err[1] != '':
                        split_err= list(map(float, split_err))
                        err_break_3=np.mean([abs(split_err[0]), abs(split_err[1])])
                    elif split_err[0] !='':
                        err_break_3=float(split_err[0])
                    elif split_err[1] != '':
                        err_break_3=float(split_err[1])
                    else:
                        err_break_3=0
                    power_4=-float(ag_instrument_data[' alpha_4 '][row])
                    if model_selector>=4:
                        break_4=int(float(ag_instrument_data[' break_4 '][row]))
                        split_err=re.split(",", ag_instrument_data.at[placement, ' D_break_4 '].strip())
                        if split_err and split_err[0] !='' and split_err[1] != '':
                            split_err= list(map(float, split_err))
                            err_break_4=np.mean([abs(split_err[0]), abs(split_err[1])])
                        elif split_err[0] !='':
                            err_break_4=float(split_err[0])
                        elif split_err[1] != '':
                            err_break_4=float(split_err[1])
                        else:
                            err_break_4=0
                        power_5=-float(ag_instrument_data[' alpha_5 '][row])
                        if model_selector==5:
                            break_5=int(float(ag_instrument_data[' break_5 '][row]))
                            split_err=re.split(",", ag_instrument_data.at[placement, ' D_break_5 '].strip())
                            if split_err and split_err[0] !='' and split[1] != '':
                                split_err= list(map(float, split_err))
                                err_break_5=np.mean([abs(split_err[0]), abs(split_err[1])])
                            elif split_err[0] !='':
                                err_break_5=float(split_err[0])
                            elif split_err[1] != '':
                                err_break_5=float(split_err[1])
                            else:
                                err_break_5=0
                            power_6=-float(ag_instrument_data[' alpha_6 '][row]) 
                            initial_value=quintuply_broken_PL(1, break_1, break_2, break_3, \
                                        break_4, break_5, power_1, power_2, power_3, power_4, \
                                            power_5, power_6, 0, 11*60*60)
                            Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_11 '][row]),\
                             0.01)[0]
                            AG_fluence=integrate.quad(lambda t: quintuply_broken_PL(Q, \
                                break_1, break_2, break_3, break_4, break_5, power_1, \
                                power_2, power_3, power_4, power_5, power_6, 0, t), \
                                                      front, back)[0]
                            AG_fluence_err=np.sqrt((\
                                quintuply_broken_PL_derivative(Q, break_1, break_2, break_3,\
                                break_4, break_5, power_1, power_2, power_3, power_4, power_5, \
                                power_6, 0, front)*error_front)**2+(\
                                quintuply_broken_PL_derivative(Q, break_1, break_2, break_3,\
                                break_4, break_5, power_1, power_2, power_3, power_4, power_5, \
                                power_6, 0, back)*error_back)**2)
                        else:
                            initial_value=quadruply_broken_PL(1, break_1, break_2, break_3, break_4, \
                                                 power_1, power_2, power_3, power_4, power_5, \
                                                 0, 11*60*60)
                            Q=fsolve(lambda A: A*initial_value-float(\
                                        ag_instrument_data[' Flux_11 '][row]),0.01)[0]
                            AG_fluence=integrate.quad(lambda t: quadruply_broken_PL(Q, \
                                break_1, break_2, break_3, break_4, power_1, power_2, \
                                power_3, power_4, power_5, 0, t), front, back)[0]
                            AG_fluence_err=np.sqrt((quadruply_broken_PL_derivative(Q, \
                                break_1, break_2, break_3, break_4, power_1, power_2, \
                                power_3, power_4, power_5, 0, front)*error_front)**2+(\
                                quadruply_broken_PL_derivative(Q, \
                                break_1, break_2, break_3, break_4, power_1, power_2, \
                                power_3, power_4, power_5, 0, back)*error_back)**2)
                    else:
                        initial_value=triply_broken_PL(1, break_1, break_2, break_3, \
                                                 power_1, power_2, power_3, power_4, \
                                                 0, 11*60*60)
                        Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_11 '][row]),\
                         0.01)[0]
                        AG_fluence=integrate.quad(lambda t: triply_broken_PL(Q, \
                                break_1, break_2, break_3, power_1, power_2, \
                                power_3, power_4, 0, t), front, back)[0]
                        AG_fluence_err=np.sqrt((triply_broken_PL_derivative(Q, \
                                break_1, break_2, break_3, power_1, power_2, \
                                power_3, power_4, 0, front)*error_front)**2+(\
                                triply_broken_PL_derivative(Q, \
                                break_1, break_2, break_3, power_1, power_2, \
                                power_3, power_4, 0, back)*error_back)**2)
                else:
                    initial_value=doubly_broken_PL(1, break_1, break_2, \
                                                 power_1, power_2, power_3, \
                                                 0, 11*60*60)
                    Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_11 '][row]),\
                         0.01)[0]
                    AG_fluence=integrate.quad(lambda t: doubly_broken_PL(Q, \
                                break_1, break_2, power_1, power_2, \
                                power_3, 0, t), front, back)[0]
                    AG_fluence_err=np.sqrt((doubly_broken_PL_derivative(Q, \
                                break_1, break_2, power_1, power_2, \
                                power_3, 0, front)*error_front)**2+(\
                                doubly_broken_PL_derivative(Q, \
                                break_1, break_2, power_1, power_2, \
                                power_3, 0, back)*error_back)**2)
            else:
                initial_value=singly_broken_PL(1, break_1, power_1, power_2, 0, 11*60*60)
                Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_11 '][row]),\
                         0.01)[0]
                AG_fluence=integrate.quad(lambda t: singly_broken_PL(Q, \
                                break_1, power_1, power_2, 0, t), front, back)[0]
                AG_fluence_err=np.sqrt(\
                                (singly_broken_PL_derivative(1, break_1, power_1, power_2,\
                                0, front)*error_front**2\
                                 +(singly_broken_PL_derivative(1, break_1, power_1, \
                                power_2, 0, back)*error_back)**2))
        else:
            initial_value=power_law(power_1, 11*60*60)
            if ag_instrument_data[' Flux_11 '][row] != ' N/A ':
                Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_11 '][row]),\
                         0.01)[0]
            elif ag_instrument_data[' Flux_24 '][row] != " N/A ":
                Q=fsolve(lambda A: A*initial_value-float(ag_instrument_data[' Flux_24 '][row]),\
                         0.01)[0]
            else:
                Q=np.NaN
            AG_fluence=integrate.quad(lambda t: Q*power_law(power_1, t), front, back)[0]
            AG_fluence_err=np.sqrt((PL_derivative(power_1, front)*error_front)**2\
                             +(PL_derivative(power_1, back)*error_back)**2)
    else:
        AG_fluence=np.NaN
    return AG_fluence_err

# %%
##Data imports 1##
#this imports every.single.Swift.burst.ever seen. I cannot comment it out as it's used to figure out where and what dat ais available. You'd think there'd be an
#automatic way to do this, but... This isn't the real issue
GRB_list=["GRB230204B", "GRB230204A", "GRB230123A",
                                              "GRB230116D", "GRB221231A", "GRB221226B",\
                                              "GRB221216A", "GRB221202A", "GRB221201A",\
                                              "GRB221120A", "GRB221110A", "GRB221028A", \
                                              "GRB221027B", "GRB221024A", "GRB221016A", \
                                              "GRB221009A", "GRB221006A", "GRB220930A",\
                                              "GRB220921A", "GRB220907A", "GRB220826A",\
                                              "GRB220813A", "GRB220730A", "GRB220715B",\
                                              "GRB220714B", "GRB220711B", "GRB220710A",\
                                              "GRB220708B", "GRB220708A", "GRB220706A", \
                                              "GRB220701A", "GRB220627A", "GRB220623A", \
                                              "GRB220618A", "GRB220611A", "GRB220527A",\
                                              "GRB220521A", "GRB220518A", "GRB220511A", \
                                              "GRB220506A", "GRB220430A", "GRB220427A", \
                                              "GRB220412A", "GRB220408A", "GRB220403B", \
                                              "GRB220325A", "GRB220319A", "GRB220306B",\
                                              "GRB220305A", "GRB220302A", "GRB220219B",\
                                              "GRB220118A", "GRB220117B", "GRB220117A", \
                                              "GRB220107B", "GRB220107A", "GRB220101A",\
                                              "GRB211227A", "GRB211221A", "GRB211211A", \
                                              "GRB211207A", "GRB211129A", "GRB211107B",\
                                              "GRB211106A", "GRB211025A", "GRB211024B",\
                                              "GRB211023B", "GRB210930A", "GRB210928A",\
                                              "GRB210919A", "GRB210912A", "GRB210905A", \
                                              "GRB210901A", "GRB210827A", "GRB210824A",\
                                              "GRB210822A", "GRB210820A", "GRB210818A",\
                                              "GRB210807C", "GRB210807A", "GRB210802A", \
                                              "GRB210731A", "GRB210730A", "GRB210726A", \
                                              "GRB210725B", "GRB210725A", "GRB210724A", \
                                              "GRB210723A", "GRB210722A", "GRB210712A", \
                                              "GRB210708A", "GRB210706A", "GRB210704A", \
                                              "GRB210702A", "GRB210626A", "GRB210619B", \
                                              "GRB210618A", "GRB210610B", "GRB210610A",\
                                              "GRB210606A", "GRB210605B", "GRB210527A", \
                                              "GRB210517A", "GRB210515C", "GRB210514A", \
                                              "GRB210509A", "GRB210504A", "GRB210421A", \
                                              "GRB210420B", "GRB210419C", "GRB210419A", \
                                              "GRB210411C", "GRB210410A", "GRB210402A", \
                                              "GRB210323A", "GRB210321A", "GRB210318B", \
                                              "GRB210318A", "GRB210308A", "GRB210307A", \
                                              "GRB210306A", "GRB210305A", "GRB210226A", \
                                              "GRB210222B", "GRB210218A", "GRB210217A", \
                                              "GRB210212A", "GRB210211A", "GRB210210A", \
                                              "GRB210209A", "GRB210207B", "GRB210205A", \
                                              "GRB210204A", "GRB210112A", "GRB210104B", \
                                              "GRB210104A", "GRB210102B", "GRB201229A", \
                                              "GRB201223A", "GRB201221D", "GRB201221A", \
                                              "GRB201216C", "GRB201209A", "GRB201208A", \
                                              "GRB201203A", "GRB201128A", "GRB201116A", \
                                              "GRB201104B", "GRB201103B", "GRB201029A",\
                                              "GRB201027A", "GRB201026A", "GRB201024A", \
                                              "GRB201021C", "GRB201020B", "GRB201020A", \
                                              "GRB201017A", "GRB201015A", "GRB201014A", \
                                              "GRB201013A", "GRB201008A", "GRB201006A", \
                                              "GRB201001A", "GRB200925B", "GRB200922A",\
                                              "GRB200917A", "GRB200907B", "GRB200906A", \
                                              "GRB200901B", "GRB200901A", "GRB200829A", \
                                              "GRB200826A", "GRB200819A", "GRB200809B", \
                                              "GRB200806A", "GRB200803A", "GRB200729A", \
                                              "GRB200716C", "GRB200715A", "GRB200714E", \
                                              "GRB200713A", "GRB200711A", "GRB200630A", \
                                              "GRB200613A", "GRB200612A", "GRB200529A",\
                                              "GRB200528A", "GRB200524A", "GRB200522A", \
                                              "GRB200519A", "GRB200517A", "GRB200512A", \
                                              "GRB200509A", "GRB200425A", "GRB200424A", \
                                              "GRB200416A", "GRB200412B", "GRB200411A",\
                                              "GRB200410A", "GRB200409A", "GRB200324A", \
                                              "GRB200306C", "GRB200306A", "GRB200303A", \
                                              "GRB200228B", "GRB200227A", "GRB200224A", \
                                              "GRB200219C", "GRB200219A", "GRB200216B", \
                                              "GRB200215A", "GRB200205B", "GRB200205A", \
                                              "GRB200131A", "GRB200127A", "GRB200125A", \
                                              "GRB200122A", "GRB200109A", "GRB200107B", \
                                              "GRB191228A", "GRB191221B", "GRB191220A", \
                                              "GRB191123A", "GRB191122A", "GRB191106A", \
                                              "GRB191101A", "GRB191031D", "GRB191031C", \
                                              "GRB191029A", "GRB191024A", "GRB191019A", \
                                              "GRB191017B", "GRB191016A", "GRB191011A", \
                                              "GRB191004B", "GRB191004A", "GRB190926A", \
                                              "GRB190919B", "GRB190829A", "GRB190828B", \
                                              "GRB190824A", "GRB190823A", "GRB190821A", \
                                              "GRB190816A", "GRB190804B", "GRB190731A", \
                                              "GRB190719C", "GRB190718A", "GRB190708B", \
                                              "GRB190706B", "GRB190701A", "GRB190630C", \
                                              "GRB190630B", "GRB190627A", "GRB190613B", \
                                              "GRB190613A", "GRB190611A", "GRB190604B", \
                                              "GRB190531B", "GRB190530A", "GRB190519A", \
                                              "GRB190515B", "GRB190511A", "GRB190422A", \
                                              "GRB190404C", "GRB190324A", "GRB190320A", \
                                              "GRB190311A", "GRB190305A", "GRB190220A", \
                                              "GRB190219A", "GRB190211A", "GRB190204A", \
                                              "GRB190203A", "GRB190202A", "GRB190129B", \
                                              "GRB190123A", "GRB190114C", "GRB190114B", \
                                              "GRB190114A", "GRB190109B", "GRB190109A", \
                                              "GRB190106A", "GRB190103B", "GRB181228A", \
                                              "GRB181213A", "GRB181203A", "GRB181202A", \
                                              "GRB181201A", "GRB181126A", "GRB181123B", \
                                              "GRB181110A", "GRB181103A", "GRB181030A", \
                                              "GRB181023A", "GRB181022A", "GRB181020A", \
                                              "GRB181013A", "GRB181010A", "GRB181003A", \
                                              "GRB181002A", "GRB180930A", "GRB180925A", \
                                              "GRB180924A", "GRB180914B", "GRB180905A", \
                                              "GRB180904A", "GRB180828A", "GRB180823A", \
                                              "GRB180821A", "GRB180818B", "GRB180818A", \
                                              "GRB180812A", "GRB180809B", "GRB180809A", \
                                              "GRB180806A", "GRB180805B", "GRB180805A",\
                                              "GRB180728A", "GRB180727A", "GRB180721A", \
                                              "GRB180720C", "GRB180720B", "GRB180709A", \
                                              "GRB180706A", "GRB180705A", "GRB180704A", \
                                              "GRB180703A", "GRB180630A", "GRB180626A", \
                                              "GRB180624A", "GRB180623A", "GRB180620B", \
                                              "GRB180620A", "GRB180618A", "GRB180614A", \
                                              "GRB180613A", "GRB180602A", "GRB180514A", \
                                              "GRB180512A", "GRB180510B", "GRB180510A", \
                                              "GRB180504A", "GRB180425A", "GRB180418A", \
                                              "GRB180411A", "GRB180410A", "GRB180409A", \
                                              "GRB180404B", "GRB180404A", "GRB180402A", \
                                              "GRB180331B", "GRB180331A", "GRB180329B", \
                                              "GRB180325A", "GRB180324A", "GRB180316A", \
                                              "GRB180314B", "GRB180314A", "GRB180311A", \
                                              "GRB180305A", "GRB180224A", "GRB180222A", \
                                              "GRB180210A", "GRB180205A", "GRB180204A", \
                                              "GRB180115A", "GRB180111A", "GRB180103A", \
                                              "GRB180102A", "GRB171222A", "GRB171216A", \
                                              "GRB171212A", "GRB171211A", "GRB171209A", \
                                              "GRB171205A", "GRB171124A", "GRB171123A", \
                                              "GRB171120A", "GRB171115A", "GRB171102B", \
                                              "GRB171027A", "GRB171020A", "GRB171010B", \
                                              "GRB171010A", "GRB171007A", "GRB171004A", \
                                              "GRB171001A", "GRB170921A", "GRB170912B", \
                                              "GRB170912A", "GRB170906C", "GRB170906B", \
                                              "GRB170906A", "GRB170903A", "GRB170827A", \
                                              "GRB170822A", "GRB170813A", "GRB170810A", \
                                              "GRB170804A", "GRB170803A", "GRB170728B", \
                                              "GRB170728A", "GRB170714A", "GRB170711A", \
                                              "GRB170710A", "GRB170705A", "GRB170629A", \
                                              "GRB170626A", "GRB170607A", "GRB170604A", \
                                              "GRB170531B", "GRB170531A", "GRB170526A", \
                                              "GRB170524B", "GRB170524A", "GRB170519A", \
                                              "GRB170516A", "GRB170510A", "GRB170428A", \
                                              "GRB170419A", "GRB170405A", "GRB170331A", \
                                              "GRB170330A", "GRB170318B", "GRB170318A", \
                                              "GRB170317A", "GRB170311A", "GRB170306A", \
                                              "GRB170214A", "GRB170208B", "GRB170208A", \
                                              "GRB170206B", "GRB170205A", "GRB170202A", \
                                              "GRB170127B", "GRB170127A", "GRB170126A", \
                                              "GRB170115A", "GRB170113A", "GRB170111A", \
                                              "GRB161224A", "GRB161220A", "GRB161219B", \
                                              "GRB161217A", "GRB161214B", "GRB161202A", \
                                              "GRB161129A", "GRB161117B", "GRB161117A", \
                                              "GRB161113A", "GRB161108A", "GRB161105A", \
                                              "GRB161104A", "GRB161023A", "GRB161022A", \
                                              "GRB161017A", "GRB161015A", "GRB161014A", \
                                              "GRB161011A", "GRB161010A", "GRB161007A", \
                                              "GRB161004B", "GRB161004A", "GRB161001A", \
                                              "GRB160927A", "GRB160917A", "GRB160912A", \
                                              "GRB160910A", "GRB160905A", "GRB160826A", \
                                              "GRB160824A", "GRB160821B", "GRB160816A", \
                                              "GRB160815A", "GRB160804A", "GRB160801A", \
                                              "GRB160716A", "GRB160712A", "GRB160705B", \
                                              "GRB160703A", "GRB160630A", "GRB160629A", \
                                              "GRB160625B", "GRB160625A", "GRB160624A", \
                                              "GRB160623A", "GRB160611A", "GRB160607A", \
                                              "GRB160601A", "GRB160525B", "GRB160521B", \
                                              "GRB160509A", "GRB160506A", "GRB160504A", \
                                              "GRB160501A", "GRB160425A", "GRB160422A", \
                                              "GRB160417A", "GRB160412A", "GRB160411A", \
                                              "GRB160410A", "GRB160408A", "GRB160327A", \
                                              "GRB160325A", "GRB160321A", "GRB160314A", \
                                              "GRB160313A", "GRB160303A", "GRB160228A", \
                                              "GRB160227A", "GRB160225A", "GRB160223A", \
                                              "GRB160221A", "GRB160220B", "GRB160220A", \
                                              "GRB160216A", "GRB160203A", "GRB160131A", \
                                              "GRB160127A", "GRB160123A", "GRB160121A", \
                                              "GRB160119A", "GRB160117B", "GRB160117A", \
                                              "GRB160104A", "GRB160101A", "GRB151229A", \
                                              "GRB151228B", "GRB151215A", "GRB151212A", \
                                              "GRB151210A", "GRB151205A", "GRB151127A", \
                                              "GRB151120A", "GRB151118A", "GRB151114A", \
                                              "GRB151112A", "GRB151111A", "GRB151031A", \
                                              "GRB151029A", "GRB151027B", "GRB151027A", \
                                              "GRB151023A", "GRB151022A", "GRB151021A", \
                                              "GRB151006A", "GRB151004A", "GRB151001B", \
                                              "GRB151001A", "GRB150925A", "GRB150915A", \
                                              "GRB150911A", "GRB150910A", "GRB150907B", \
                                              "GRB150902A", "GRB150831B", "GRB150831A", \
                                              "GRB150821A", "GRB150819A", "GRB150818A", \
                                              "GRB150817A", "GRB150811A", "GRB150801B", \
                                              "GRB150728A", "GRB150727A", "GRB150724B", \
                                              "GRB150724A", "GRB150722A", "GRB150720A", \
                                              "GRB150716A", "GRB150711A", "GRB150710B", \
                                              "GRB150627A", "GRB150626B", "GRB150626A", \
                                              "GRB150622A", "GRB150616A", "GRB150615A", \
                                              "GRB150608A", "GRB150607A", "GRB150530A", \
                                              "GRB150527A", "GRB150523A", "GRB150518A", \
                                              "GRB150514A", "GRB150430A", "GRB150428B", \
                                              "GRB150428A", "GRB150424A", "GRB150423A", \
                                              "GRB150403A", "GRB150323C", "GRB150323B", \
                                              "GRB150323A", "GRB150318A", "GRB150317A", \
                                              "GRB150314A", "GRB150309A", "GRB150302A", \
                                              "GRB150301B", "GRB150301A", "GRB150222A", \
                                              "GRB150219A", "GRB150213B", "GRB150212A", \
                                              "GRB150211A", "GRB150206A", "GRB150204A", \
                                              "GRB150203A", "GRB150202B", "GRB150202A", \
                                              "GRB150201A", "GRB150120B", "GRB150120A", \
                                              "GRB150110B", "GRB150103A", "GRB150101B", \
                                              "GRB150101A", "GRB141225A", "GRB141221A", \
                                              "GRB141220A", "GRB141215A", "GRB141212B", \
                                              "GRB141212A", "GRB141130A", "GRB141121A", \
                                              "GRB141109B", "GRB141109A", "GRB141031B", \
                                              "GRB141031A", "GRB141028A", "GRB141026A", \
                                              "GRB141022A", "GRB141020A", "GRB141017A", \
                                              "GRB141015A", "GRB141005A", "GRB141004A", \
                                              "GRB140930B", "GRB140928A", "GRB140927A", \
                                              "GRB140919A", "GRB140916A", "GRB140909A", \
                                              "GRB140907A", "GRB140903A", "GRB140824A", \
                                              "GRB140818B", "GRB140818A", "GRB140817A", \
                                              "GRB140815A", "GRB140808A", "GRB140801A", \
                                              "GRB140730A", "GRB140719A", "GRB140716A", \
                                              "GRB140713A", "GRB140710A", "GRB140709B", \
                                              "GRB140709A", "GRB140706A", "GRB140703A", \
                                              "GRB140629A", "GRB140628A", "GRB140626A", \
                                              "GRB140623A", "GRB140622A", "GRB140620A", \
                                              "GRB140619B", "GRB140619A", "GRB140614B", \
                                              "GRB140614A", "GRB140611A", "GRB140610A", \
                                              "GRB140606B", "GRB140529A", "GRB140518A", \
                                              "GRB140516A", "GRB140515A", "GRB140512A", \
                                              "GRB140509A", "GRB140508A", "GRB140506A", \
                                              "GRB140502A", "GRB140430A", "GRB140428A", \
                                              "GRB140423A", "GRB140419A", "GRB140413A", \
                                              "GRB140412A", "GRB140408A", "GRB140331A", \
                                              "GRB140323A", "GRB140320C", "GRB140320B", \
                                              "GRB140320A", "GRB140318A", "GRB140311B", \
                                              "GRB140311A", "GRB140304A", "GRB140302A", \
                                              "GRB140301A", "GRB140226A", "GRB140215A", \
                                              "GRB140213A", "GRB140211A", "GRB140209A", \
                                              "GRB140206A", "GRB140129B", "GRB140129A", \
                                              "GRB140114A", "GRB140108A", "GRB140104B", \
                                              "GRB140103A", "GRB140102A", "GRB131231A", \
                                              "GRB131229A", "GRB131227A", "GRB131205A", \
                                              "GRB131202A", "GRB131128A", "GRB131127A", \
                                              "GRB131122A", "GRB131117A", "GRB131108A", \
                                              "GRB131105A", "GRB131103A", "GRB131031A", \
                                              "GRB131030A", "GRB131024B", "GRB131024A", \
                                              "GRB131018B", "GRB131018A", "GRB131014A", \
                                              "GRB131011A", "GRB131004A", "GRB131002B", \
                                              "GRB131002A", "GRB130929A", "GRB130925A", \
                                              "GRB130912A", "GRB130907A", "GRB130903A", \
                                              "GRB130831B", "GRB130831A", "GRB130822A", \
                                              "GRB130816B", "GRB130816A", "GRB130812A", \
                                              "GRB130807A", "GRB130806A", "GRB130803A", \
                                              "GRB130727A", "GRB130725B", "GRB130725A", \
                                              "GRB130722A", "GRB130716A", "GRB130702A", \
                                              "GRB130701A", "GRB130627B", "GRB130627A", \
                                              "GRB130625A", "GRB130623A", "GRB130615A", \
                                              "GRB130612A", "GRB130610A", "GRB130609B", \
                                              "GRB130609A", "GRB130608A", "GRB130606B", \
                                              "GRB130606A", "GRB130605A", "GRB130604A", \
                                              "GRB130603B", "GRB130603A", "GRB130529A", \
                                              "GRB130528A", "GRB130527A", "GRB130518A", \
                                              "GRB130515A", "GRB130514B", "GRB130514A", \
                                              "GRB130513A", "GRB130511A", "GRB130508A", \
                                              "GRB130505B", "GRB130505A", "GRB130504C", \
                                              "GRB130504A", "GRB130502B", "GRB130502A", \
                                              "GRB130427B", "GRB130427A", "GRB130420B", \
                                              "GRB130420A", "GRB130418A", "GRB130408A", \
                                              "GRB130327B", "GRB130327A", "GRB130315A", \
                                              "GRB130313A", "GRB130306A", "GRB130305A", \
                                              "GRB130211A", "GRB130206A", "GRB130131B", \
                                              "GRB130131A", "GRB130122A", "GRB130102A", \
                                              "GRB121229A", "GRB121226A", "GRB121217A", \
                                              "GRB121212A", "GRB121211A", "GRB121209A", \
                                              "GRB121202A", "GRB121201A", "GRB121128A", \
                                              "GRB121125A", "GRB121123A", "GRB121117A", \
                                              "GRB121108A", "GRB121102A", "GRB121031A", \
                                              "GRB121028A", "GRB121027A", "GRB121025A", \
                                              "GRB121024A", "GRB121017A", "GRB121011A", \
                                              "GRB121001A", "GRB120927A", "GRB120923A", \
                                              "GRB120922A", "GRB120911A", "GRB120909A", \
                                              "GRB120907A", "GRB120819A", "GRB120817A", \
                                              "GRB120816A", "GRB120815A", "GRB120811C", \
                                              "GRB120811A", "GRB120807A", "GRB120805A", \
                                              "GRB120804A", "GRB120803B", "GRB120803A", \
                                              "GRB120802A", "GRB120729A", "GRB120728A", \
                                              "GRB120724A", "GRB120722A", "GRB120716A", \
                                              "GRB120714A", "GRB120712A", "GRB120711B", \
                                              "GRB120711A", "GRB120703A", "GRB120701A", \
                                              "GRB120630A", "GRB120624B", "GRB120612A", \
                                              "GRB120521C", "GRB120521B", "GRB120521A", \
                                              "GRB120514A", "GRB120422A", "GRB120419A", \
                                              "GRB120404A", "GRB120403B", "GRB120401A", \
                                              "GRB120328A", "GRB120327A", "GRB120326A", \
                                              "GRB120324A", "GRB120320A", "GRB120312A", \
                                              "GRB120311B", "GRB120311A", "GRB120308A", \
                                              "GRB120305A", "GRB120302A", "GRB120224A", \
                                              "GRB120219A", "GRB120215A", "GRB120213A", \
                                              "GRB120212A", "GRB120211A", "GRB120202A", \
                                              "GRB120121A", "GRB120119A", "GRB120118B", \
                                              "GRB120116A", "GRB120106A", "GRB120102A", \
                                              "GRB111229A", "GRB111228A", "GRB111225A", \
                                              "GRB111222A", "GRB111215B", "GRB111215A", \
                                              "GRB111212A", "GRB111211A", "GRB111210A", \
                                              "GRB111209A", "GRB111208A", "GRB111205A", \
                                              "GRB111204A", "GRB111201A", "GRB111129A", \
                                              "GRB111123A", "GRB111121A", "GRB111117A", \
                                              "GRB111109A", "GRB111107A", "GRB111103B", \
                                              "GRB111029A", "GRB111022B", "GRB111022A", \
                                              "GRB111020A", "GRB111018A", "GRB111016A", \
                                              "GRB111008A", "GRB110928A", "GRB110921A", \
                                              "GRB110918A", "GRB110915B", "GRB110915A", \
                                              "GRB110903A", "GRB110820A", "GRB110818A", \
                                              "GRB110808A", "GRB110802A", "GRB110801A", \
                                              "GRB110731A", "GRB110726A", "GRB110719A", \
                                              "GRB110715A", "GRB110709B", "GRB110709A", \
                                              "GRB110708A", "GRB110625A", "GRB110610A", \
                                              "GRB110604A", "GRB110530A", "GRB110521A", \
                                              "GRB110520A", "GRB110503A", "GRB110428A", \
                                              "GRB110422A", "GRB110420A", "GRB110414A", \
                                              "GRB110411A", "GRB110407A", "GRB110402A", \
                                              "GRB110319B", "GRB110319A", "GRB110318B", \
                                              "GRB110315A", "GRB110312A", "GRB110305A", \
                                              "GRB110223B", "GRB110223A", "GRB110213B", \
                                              "GRB110213A", "GRB110210A", "GRB110208A", \
                                              "GRB110206A", "GRB110205A", "GRB110201A", \
                                              "GRB110128A", "GRB110119A", "GRB110112A", \
                                              "GRB110107A", "GRB110106B", "GRB110106A", \
                                              "GRB110102A", "GRB101225A", "GRB101224A", \
                                              "GRB101219B", "GRB101219A", "GRB101213A", \
                                              "GRB101204A", "GRB101201A", "GRB101117B", \
                                              "GRB101114A", "GRB101112A", "GRB101030A", \
                                              "GRB101024A", "GRB101023A", "GRB101017A", \
                                              "GRB101011A", "GRB101008A", "GRB100915A", \
                                              "GRB100909A", "GRB100906A", "GRB100905A", \
                                              "GRB100902A", "GRB100901A", "GRB100823A", \
                                              "GRB100816A", "GRB100814A", "GRB100807A", \
                                              "GRB100805A", "GRB100802A", "GRB100728B", \
                                              "GRB100728A", "GRB100727A", "GRB100725B", \
                                              "GRB100725A", "GRB100724A", "GRB100713A", \
                                              "GRB100704A", "GRB100702A", "GRB100628A", \
                                              "GRB100625A", "GRB100621A", "GRB100619A", \
                                              "GRB100615A", "GRB100614A", "GRB100606A", \
                                              "GRB100528A", "GRB100526B", "GRB100526A", \
                                              "GRB100522A", "GRB100518A", "GRB100514A", \
                                              "GRB100513A", "GRB100508A", "GRB100504A", \
                                              "GRB100425A", "GRB100424A", "GRB100420A", \
                                              "GRB100418A", "GRB100414A", "GRB100413A", \
                                              "GRB100331B", "GRB100316D", "GRB100316C", \
                                              "GRB100316B", "GRB100316A", "GRB100305A", \
                                              "GRB100302A", "GRB100219A", "GRB100213B", \
                                              "GRB100213A", "GRB100212A", "GRB100206A", \
                                              "GRB100205A", "GRB100117A", "GRB100115A", \
                                              "GRB100111A", "GRB100103A", "GRB091230", \
                                              "GRB091221", "GRB091208B", "GRB091202", \
                                              "GRB091130B", "GRB091127", "GRB091111", \
                                              "GRB091109B", "GRB091109A", "GRB091104", \
                                              "GRB091102", "GRB091029", "GRB091026", \
                                              "GRB091024", "GRB091020", "GRB091018", \
                                              "GRB091010", "GRB091003", "GRB090929B", \
                                              "GRB090927", "GRB090926B", "GRB090926A", \
                                              "GRB090915", "GRB090912", "GRB090904B", \
                                              "GRB090904A", "GRB090902B", "GRB090831C", \
                                              "GRB090827", "GRB090823", "GRB090817", \
                                              "GRB090814B", "GRB090814A", "GRB090813", \
                                              "GRB090812", "GRB090809", "GRB090807", \
                                              "GRB090728", "GRB090727", "GRB090726",\
                                              "GRB090720", "GRB090715B", "GRB090709A", \
                                              "GRB090702", "GRB090628", "GRB090625B", \
                                              "GRB090621B", "GRB090621A", "GRB090618", \
                                              "GRB090607", "GRB090531B", "GRB090531A", \
                                              "GRB090530", "GRB090529", "GRB090519", \
                                              "GRB090518", "GRB090516", "GRB090515", \
                                              "GRB090510", "GRB090429B", "GRB090429A", \
                                              "GRB090426", "GRB090424", "GRB090423", \
                                              "GRB090422", "GRB090419", "GRB090418A", \
                                              "GRB090417B", "GRB090407", "GRB090404", \
                                              "GRB090401B", "GRB090328A", "GRB090323", \
                                              "GRB090313", "GRB090309", "GRB090308", \
                                              "GRB090307", "GRB090306B", "GRB090205", \
                                              "GRB090201", "GRB090126A", "GRB090123", \
                                              "GRB090118", "GRB090117", "GRB090113", \
                                              "GRB090111", "GRB090107B", "GRB090102", \
                                              "GRB081230", "GRB081228", "GRB081226A", \
                                              "GRB081224", "GRB081222", "GRB081221", \
                                              "GRB081211", "GRB081210", "GRB081204", \
                                              "GRB081203B", "GRB081203A", "GRB081128", \
                                              "GRB081127", "GRB081126", "GRB081121", \
                                              "GRB081118", "GRB081109", "GRB081105", \
                                              "GRB081104", "GRB081102", "GRB081029", \
                                              "GRB081028", "GRB081025", "GRB081024A", \
                                              "GRB081016B", "GRB081016A", "GRB081012", \
                                              "GRB081011", "GRB081008", "GRB081007", \
                                              "GRB081003A", "GRB081001", "GRB080928", \
                                              "GRB080919", "GRB080916C", "GRB080916B", \
                                              "GRB080916A", "GRB080915A", "GRB080913", \
                                              "GRB080906", "GRB080905B", "GRB080905A", \
                                              "GRB080903", "GRB080828", "GRB080825B"]
# THIS IS. This line has to specify all the overlapping GBM-BAT GRBs, and it's not flexible because every single one had to be named in the standard way,
# which the Fermi-GBM data /doesn't currently include/. So I had to find all their names by hand. I would recommend future users find a better way.
GRB_search_list=['GRB220711B', 'GRB220618A', 'GRB220611A', 'GRB220518A', 'GRB220403B', \
                 'GRB220325A', 'GRB211211A', 'GRB211129A', 'GRB211024B', 'GRB210619B', \
                 'GRB210610A', 'GRB210514A', 'GRB210307A', 'GRB210226A', 'GRB210222B', \
                 'GRB201229A', 'GRB201027A', 'GRB201014A', 'GRB200901B', 'GRB200729A', \
                 'GRB200713A', 'GRB200630A', 'GRB200528A', 'GRB200522A', 'GRB200410A', \
                 'GRB200228B', 'GRB200205B', 'GRB191011A', 'GRB191004A', 'GRB190828B', \
                 'GRB190824A', 'GRB190821A', 'GRB190719C', 'GRB190718A', 'GRB190604B', \
                 'GRB190422A', 'GRB190311A', 'GRB190219A', 'GRB190114A', 'GRB181202A', \
                 'GRB181013A', 'GRB180925A', 'GRB180809B', 'GRB180805A', 'GRB180720C', \
                 'GRB180720B', 'GRB180706A', 'GRB180630A', 'GRB180626A', 'GRB180620A', \
                 'GRB180402A', 'GRB180331B', 'GRB171222A', 'GRB171209A', 'GRB171115A', \
                 'GRB170903A', 'GRB170822A', 'GRB170728B', 'GRB170714A', 'GRB170710A', \
                 'GRB170629A', 'GRB170626A', 'GRB170607A', 'GRB170604A', 'GRB170202A', \
                 'GRB170111A', 'GRB161129A', 'GRB161117A', 'GRB161113A', 'GRB161011A', \
                 'GRB160905A', 'GRB160824A', 'GRB160705B', 'GRB151229A', 'GRB151228B', \
                 'GRB151205A', 'GRB151112A', 'GRB151031A', 'GRB151021A', 'GRB150925A', \
                 'GRB150801B', 'GRB150626B', 'GRB150430A', 'GRB150428B', 'GRB150323C', \
                 'GRB150323A', 'GRB150213B', 'GRB150202A', 'GRB150101B', 'GRB141220A', \
                 'GRB141212A', 'GRB141026A', 'GRB140916A', 'GRB140817A', 'GRB140730A', \
                 'GRB140706A', 'GRB140703A', 'GRB140629A', 'GRB140518A', 'GRB140430A', \
                 'GRB140331A', 'GRB140301A', 'GRB140114A', 'GRB140103A', 'GRB131127A', \
                 'GRB131103A', 'GRB130831B', 'GRB130608A', 'GRB130504A', 'GRB130420A', \
                 'GRB130418A', 'GRB130211A', 'GRB121209A', 'GRB121125A', 'GRB121123A', \
                 'GRB121031A', 'GRB121027A', 'GRB120907A', 'GRB120811C', 'GRB120811A', \
                 'GRB120724A', 'GRB120612A', 'GRB120324A', 'GRB120116A', 'GRB111229A', \
                 'GRB111225A', 'GRB111129A', 'GRB110530A', 'GRB110319B', 'GRB110319A', \
                 'GRB110210A', 'GRB110106B', 'GRB110102A', 'GRB101219B', 'GRB101030A', \
                 'GRB101024A', 'GRB100915A', 'GRB100902A', 'GRB100727A', 'GRB100725B', \
                 'GRB100704A', 'GRB100621A', 'GRB100614A', 'GRB100514A', 'GRB100425A', \
                 'GRB100115A', 'GRB091130B', 'GRB090904B', 'GRB090904A', 'GRB090812', \
                 'GRB090621A', 'GRB090618', 'GRB090531A', 'GRB090429B', 'GRB090423', \
                 'GRB090418A', 'GRB090111', 'GRB081221', 'GRB081210', 'GRB081118', \
                 'GRB080915A', 'GRB080903', 'GRB080707']
#Yeah, this dictionary will make life 500x easier later--this definitely searches for known progenitors.
Known_Precursors = {}
Known_Precursors["Long Collapsars"] = ('GRB050416A   ', 'GRB081007   ', 'GRB091127   ', \
                                       'GRB050525A   ', 'GRB050824   ', 'GRB060218   ', \
                                       'GRB060729   ', 'GRB060904B   ', 'GRB070419A   ', \
                                       'GRB071025   ', 'GRB071112C   ', 'GRB080109   ', \
                                       'GRB080319B   ', 'GRB081007A   ', 'GRB090618   ', \
                                       'GRB091127   ', 'GRB100316D   ', 'GRB100418A   ', \
                                       'GRB101219B   ', 'GRB101225A   ', 'GRB111209A   ', \
                                       'GRB111211A   ', 'GRB111228A   ', 'GRB120422A   ', \
                                       'GRB120714B   ', 'GRB120729A   ', 'GRB130215A   ', \
                                       'GRB130427A   ', 'GRB130702A   ', 'GRB130831A   ', \
                                       'GRB140206A   ', 'GRB140606B   ', 'GRB150818A   '\
                                       'GRB161219B   ', 'GRB161228B   ', 'GRB171010A   ', \
                                       'GRB171205A   ', 'GRB180720B   ', 'GRB180728A   ', \
                                       'GRB190114C   ', 'GRB190829A   ', 'GRB221009A   ', \
                                       'GRB211023A   ', 'GRB200826A   ', 'GRB210210A   ')
Known_Precursors["Short Collapsars"] = ('GRB200826A   ')
Known_Precursors["Short Mergers"] = ('GRB130603B   ', 'GRB160821B   ', 'GRB200522A   ', \
                                     'GRB150101B   ', 'GRB160624A   ', 'GRB170817A   ', \
                                     'GRB070809   ')
Known_Precursors["Long Mergers"] = ('GRB211211A   ', 'GRB230307A   ', 'GRB120304B   ', \
                                    'GRB111005A   ', 'GRB060614   ')
Known_Precursors["Potentially Exotic"] = ("GRB210704A   ")
Known_Precursors["Galactic Detected"] = ("GRB050509B   ", "GRB050709   ", "GRB051210    ", \
                                         "GRB070714B   ", "GRB071227    ", "GRB080503    ", \
                                         "GRB080905A   ", "GRB090515    ", "GRB160303A   ")
###
Fermi_Precursors = {}
Fermi_Precursors["Long Collapsars"] = ('GRB050416461', 'GRB050525002', 'GRB050824966', \
                                       'GRB060218148', 'GRB060729800', 'GRB060904104', \
                                       'GRB070419447', 'GRB071025172', 'GRB071112772', \
                                       'GRB080319258', 'GRB081007224', 'GRB090618353', \
                                       'GRB091127976', 'GRB100316531', 'GRB100418882', \
                                       'GRB101219686', 'GRB101225776', 'GRB111209300', \
                                       'GRB111211928', 'GRB111228656', 'GRB120422300', \
                                       'GRB120714888', 'GRB120729455', 'GRB130215063', \
                                       'GRB130427324', 'GRB130702003', 'GRB130831544', \
                                       'GRB140206303', 'GRB140606133', 'GRB150818483', \
                                       'GRB161219783', 'GRB161228552' ,'GRB171010792', \
                                       'GRB171205306', 'GRB180720598', 'GRB180728728', \
                                       'GRB190114872', 'GRB190829830', 'GRB200826187', \
                                       'GRB210210083', 'GRB211023545', 'GRB221009553')
Fermi_Precursors["Short Collapsars"] = ('GRB200826187')
Fermi_Precursors["Short Mergers"] = ('GRB070809807', 'GRB130603659', 'GRB150101641', \
                                     'GRB160624477', 'GRB160821936', 'GRB170817528', \
                                     'GRB200522487')
Fermi_Precursors["Long Mergers"] = ('GRB060614530', 'GRB111005336', 'GRB120304248', \
                                    'GRB211211548', 'GRB230307655')
Fermi_Precursors["Potentially Exotic"] = ('GRB210704814')
Fermi_Precursors["Galactic Detected"] = ('GRB050509166', 'GRB050709942', 'GRB051210240', \
                                         'GRB070714207', 'GRB071227842', 'GRB080503518', \
                                         'GRB080905499', 'GRB090515198', 'GRB160303454')
###
Overlap_Precursors = {}
Overlap_Precursors["Long Collapsars"] = (050416461.0, 50525002.0, 050824966.0, \
                                       060218148.0, 060729800.0, 060904104.0, \
                                       070419447.0, 071025172.0, 071112772.0, \
                                       080319258.0, 081007224.0, 090618353.0, \
                                       091127976.0, 100316531.0, 100418882.0, \
                                       101219686.0, 101225776.0, 111209300.0, \
                                       111211928.0, 111228656.0, 120422300.0, \
                                       120714888.0, 120729455.0, 130215063.0, \
                                       130427324.0, 130702003.0, 130831544.0, \
                                       140206303.0, 140606133.0, 150818483.0, \
                                       161219783.0, 161228552.0, 171010792.0, \
                                       171205306.0, 180720598.0, 180728728.0, \
                                       190114872.0, 190829830.0, 200826187.0, \
                                       210210083.0, 211023545.0, 221009553.0)
Overlap_Precursors["Short Collapsars"] = (200826187.0)
Overlap_Precursors["Short Mergers"] = (070809807.0, 130603659.0, 150101641.0, \
                                     160624477.0, 160821936.0, 170817528.0, \
                                     200522487.0)
Overlap_Precursors["Long Mergers"] = (060614530.0, 111005336.0, 120304248.0, \
                                    211211548.0, 230307655.0)
Overlap_Precursors["Potentially Exotic"] = (210704814.0)
Overlap_Precursors["Galactic Detected"] = (050509166.0, 050709942.0, 051210240.0, \
                                         070714207.0, 071227842.0, 080503518.0, \
                                         080905499.0, 090515198.0, 160303454.0)

fermi_data=pd.read_csv('{}/22_May_2023_Fermi_Data.csv'.format(directory_path))
swift_data=pd.read_csv('{}/real swift bat data 3.csv'.format(directory_path))
swift_model_PL=pd.read_csv('{}/swift_models_3_PL.csv'.format(directory_path))
swift_model_CPL=pd.read_csv('{}/swift_models_3_PL.csv'.format(directory_path))
swift_fluxes_PL=pd.read_csv('{}/swift_fluxes_3_PL.csv'.format(directory_path))
swift_fluxes_CPL=pd.read_csv('{}/swift_fluxes_3_PL.csv'.format(directory_path))
swift_fluences_PL=pd.read_csv('{}/swift_fluences_3_PL.csv'.format(directory_path))
swift_fluences_CPL=pd.read_csv('{}/swift_fluences_3_PL.csv'.format(directory_path))
swift_fluences_reference=pd.read_csv('{}/swift_fluences_3_match.csv'.format(directory_path))
swift_redshifts=pd.read_csv('{}/swift_redshifts.csv'.format(directory_path))
second_redshifts=pd.read_csv("{}/changed_alternative_GCN_redshifts.csv".format(directory_path))
more_xrt_data=pd.read_csv('{}/XRT_most_detailed_data.csv'.format(directory_path))
weird_xrt_data=pd.DataFrame.copy(more_xrt_data)
for k in range(0, len(weird_xrt_data)):
    weird_xrt_data['GRB '][k]='0{}'.format(weird_xrt_data['GRB '][k])
average_decay_data=pd.read_csv('{}/Stage_2_decay_rates.csv'.format(directory_path))

# %%
linkin_park=[] #they wrote numb and that kinda sounds like num which is short for number...
#look, I'm tired of coming up with variable names
FOB=[] #this is going to be the GRB name because nothing matters, but genre-defining
MCR= [] #this can hold the start times, since these guys broke up first and stayed broke up
PatD=[] #and this the stop times because whatever, who cares what it's called
piloots_21=[] #had to think up another one...
Owl_City= [] # actually 2
AG_Goldstein_Data=np.zeros((len(GRB_list)+1,4)) 
#but I will force this name to make sense--it will not be something weird like post-9/11 punk
p=0
for i in range(0, len(GRB_list)):
    GRB_string=GRB_list[i]
    if GRB_string in GRB_search_list:
        good_WT_string='{}_WTCURVE.qdp'.format(GRB_string)
        good_PC_string='{}_PCCURVE.qdp'.format(GRB_string)
        bad_WT_string='{}_WTCURVE_incbad.qdp'.format(GRB_string)
        bad_PC_string='{}_PCCURVE_incbad.qdp'.format(GRB_string)
        if os.path.exists('{}/xrt_lcs/{}'.format(directory_path, good_WT_string)):
            data = ascii.read('{}/xrt_lcs/{}'.format(directory_path, good_WT_string), header_start=2)
            start_time_1, stop_time_1=friendly_loop(data)
            err_start_time_1=0.002
            err_stop_time_1=0.002
            if os.path.exists('{}/xrt_lcs/{}'.format(directory_path, good_PC_string)):
                data = ascii.read('{}/xrt_lcs/{}'.format(directory_path, good_PC_string), header_start=2)
                start_time_2, stop_time_2=friendly_loop(data)
                err_start_time_2=2.5
                err_stop_time_2=2.5
            elif os.path.exists('{}/xrt_lcs/{}'.format(directory_path, bad_PC_string)):
                data = ascii.read('{}/xrt_lcs/{}'.format(directory_path, bad_PC_string), header_start=2)
                start_time_2, stop_time_2=friendly_loop(data)
                err_start_time_2=2.5
                err_stop_time_2=2.5
            elif os.path.exists('{}/xrt_lcs/{}'.format(directory_path, bad_WT_string)):
                data = ascii.read('{}/xrt_lcs/{}'.format(directory_path, bad_WT_string), header_start=2)
                start_time_2, stop_time_2=friendly_loop(data)
                err_start_time_2=0.002
                err_stop_time_2=0.002
            else:
                start_time_2=start_time_1+10
                stop_time_2=stop_time_1-10
                err_start_time_2=2.5
                err_stop_time_2=2.5
            start_time=np.min([start_time_1, start_time_2])
            if np.where([start_time_1, start_time_2]==\
                        np.min([start_time_1, start_time_2]))[0]==0:
                err_start_time=err_start_time_1
            else:
                err_start_time=err_start_time_2
            stop_time=np.max([stop_time_1, stop_time_2])
            if np.where([stop_time_1, stop_time_2]==\
                        np.max([stop_time_1, stop_time_2]))[0]==0:
                err_stop_time=err_stop_time_1
            else:
                err_stop_time=err_stop_time_2
            linkin_park=np.append(linkin_park, GRB_search_list.index(GRB_string))
            FOB=np.append(FOB, GRB_string)
            MCR=np.append(MCR, round(start_time, 3))
            PatD=np.append(PatD, round(stop_time, 3))
            piloots_21=np.append(piloots_21, err_start_time)
            Owl_City=np.append(Owl_City, err_stop_time)
        elif os.path.exists('{}/xrt_lcs/{}'.format(directory_path, good_PC_string)):
            data = ascii.read('{}/xrt_lcs/{}'.format(directory_path, good_PC_string), header_start=2)
            start_time_1, stop_time_1=friendly_loop(data)
            err_start_time_1=2.5
            err_stop_time_1=2.5
            if os.path.exists('{}/xrt_lcs/{}'.format(directory_path,bad_WT_string)):
                data = ascii.read('{}/xrt_lcs/{}'.format(directory_path, bad_WT_string), header_start=2)
                start_time_2, stop_time_2=friendly_loop(data)
                err_start_time_2=0.002
                err_stop_time_2=0.002
            elif os.path.exists('{}/xrt_lcs/{}'.format(directory_path, bad_PC_string)):
                data = ascii.read('{}/xrt_lcs/{}'.format(directory_path, bad_PC_string), header_start=2)
                start_time_2, stop_time_2=friendly_loop(data)
                err_start_time_2=2.5
                err_stop_time_2=2.5
            else:
                start_time_2=start_time_1+10
                stop_time_2=stop_time_1-10
                err_start_time_2=2.5
                err_stop_time_2=2.5
            start_time=np.min([start_time_1, start_time_2])
            if np.where([start_time_1, start_time_2]==\
                        np.min([start_time_1, start_time_2]))[0]==0:
                err_start_time=err_start_time_1
            else:
                err_start_time=err_start_time_2
            stop_time=np.max([stop_time_1, stop_time_2])
            if np.where([stop_time_1, stop_time_2]==\
                        np.max([stop_time_1, stop_time_2]))[0]==0:
                err_stop_time=err_stop_time_1
            else:
                err_stop_time=err_stop_time_2
            linkin_park=np.append(linkin_park, GRB_search_list.index(GRB_string))
            FOB=np.append(FOB, GRB_string)
            MCR=np.append(MCR, round(start_time, 1))
            PatD=np.append(PatD, round(stop_time, 1))
            piloots_21=np.append(piloots_21, err_start_time)
            Owl_City=np.append(Owl_City, err_stop_time)
        elif os.path.exists('{}/xrt_lcs/{}'.format(directory_path, bad_WT_string)):
            data = ascii.read('{}/xrt_lcs/{}'.format(directory_path,bad_WT_string), header_start=2)
            start_time_1, stop_time_1=friendly_loop(data)
            err_start_time_1=0.002
            err_stop_time_1=0.002
            if os.path.exists('{}/xrt_lcs/{}'.format(directory_path, bad_PC_string)):
                data = ascii.read('{}/xrt_lcs/{}'.format(directory_path, bad_PC_string), header_start=2)
                start_time_2, stop_time_2=friendly_loop(data)
                err_start_time_2=2.5
                err_stop_time_2=2.5
            else:
                start_time_2=start_time_1+10
                stop_time_2=stop_time_1-10
                err_start_time_2=2.5
                err_stop_time_2=2.5
            start_time=np.min([start_time_1, start_time_2])
            if np.where([start_time_1, start_time_2]==\
                        np.min([start_time_1, start_time_2]))[0]==0:
                err_start_time=err_start_time_1
            else:
                err_start_time=err_start_time_2
            stop_time=np.max([stop_time_1, stop_time_2])
            if np.where([stop_time_1, stop_time_2]==\
                        np.max([stop_time_1, stop_time_2]))[0]==0:
                err_stop_time=err_stop_time_1
            else:
                err_stop_time=err_stop_time_2
            linkin_park=np.append(linkin_park, GRB_search_list.index(GRB_string))
            FOB=np.append(FOB, GRB_string)
            MCR=np.append(MCR, round(start_time, 3))
            PatD=np.append(PatD, round(stop_time, 3))
            piloots_21=np.append(piloots_21, err_start_time)
            Owl_City=np.append(Owl_City, err_stop_time)
        elif os.path.exists('{}/xrt_lcs/{}'.format(directory_path, bad_PC_string)):
            data = ascii.read('{}/xrt_lcs/{}'.format(directory_path, bad_PC_string), header_start=2)
            start_time, stop_time=friendly_loop(data)
            err_start_time=2.5
            err_stop_time=2.5
            linkin_park=np.append(linkin_park, GRB_search_list.index(GRB_string))
            FOB=np.append(FOB, GRB_string)
            MCR=np.append(MCR, round(start_time, 1))
            PatD=np.append(PatD, round(stop_time, 1))
            piloots_21=np.append(piloots_21, err_start_time)
            Owl_City=np.append(Owl_City, err_stop_time)
        else:
            print("{} files could not be found".format(GRB_string))
Linkin_Park=pd.DataFrame(linkin_park)
Fallout_Boi=pd.DataFrame(FOB) #I'm not typing out these dumb band names
MyChemRom=pd.DataFrame(MCR)
Panic=pd.DataFrame(PatD)
TwentyOnePilots=pd.DataFrame(piloots_21)
Adam_Young=pd.DataFrame(Owl_City)
Afterglow_of_Sample_3_Data=pd.DataFrame(data=pd.concat([Linkin_Park,Fallout_Boi, MyChemRom, \
                                    TwentyOnePilots, Panic, Adam_Young], axis=1).values,\
                                    columns=["Match_Num", "Name", "AG_Start_Time",\
                                    "AG_Start_Error", "AG_Stop_Time", "AG_Stop_Error"])

# %%
##Creating Sample 0##
tiny_removals=[]
singular_swift_fluences_list=np.zeros((1,8))
singular_swift_fluxes_list=np.zeros((1,4))
singular_swift_models_list=np.zeros((1,7))
tarjete=0
stage=0
for i in range(0,len(swift_fluences_reference)):
    #for each event
    placeholder=[[swift_data.at[i, 'GRBname '], swift_fluence_function(i, ' 25_50kev '), \
                swift_fluence_function(i, ' 25_50kev_low '), \
                swift_fluence_function(i, ' 25_50kev_hi '), \
                swift_fluence_function(i, ' 50_100kev '),\
                                          swift_fluence_function(i, ' 50_100kev_low '),\
                                          swift_fluence_function(i, ' 50_100kev_hi '), 
                                          swift_fluence_function(i, ' 15_350kev ')]]
    singular_swift_fluences_list=np.append(singular_swift_fluences_list, \
                placeholder, axis=0)
    other_placeholder=[[swift_data.at[i, 'GRBname '], swift_flux_function(i, ' 15_350kev '), \
                swift_fluence_function(i, ' 15_350kev_low '), \
                swift_fluence_function(i, ' 15_350kev_hi ')]]
    singular_swift_fluxes_list=np.append(singular_swift_fluxes_list, \
                other_placeholder, axis=0)
    third_placeholder=[[swift_data.at[i, 'GRBname '], swift_model_function(i, ' alpha '), \
                swift_model_function(i, ' alpha_low '), \
                        swift_model_function(i, ' alpha_hi '), \
                        swift_model_function(i, ' Epeak '), \
                        swift_model_function(i, ' Epeak_low '), \
                        swift_model_function(i, ' Epeak_hi ')]]
    singular_swift_models_list=np.append(singular_swift_models_list, \
                third_placeholder, axis=0)
singular_swift_fluences_list=singular_swift_fluences_list[1:]
singular_swift_fluxes_list=singular_swift_fluxes_list[1:]
singular_swift_models_list=singular_swift_models_list[1:]
swift_fluences=pd.DataFrame(singular_swift_fluences_list, columns=['GRBname ', ' 25_50kev ',\
                                                                  ' 25_50kev_low ', \
                                                                   ' 25_50kev_hi ', \
                                                                   ' 50_100kev ', \
                                                                   ' 50_100kev_low ', \
                                                                   ' 50_100kev_hi ', \
                                                                   ' 15_350kev '])
swift_fluxes=pd.DataFrame(singular_swift_fluxes_list, columns=['GRBname ', ' 15_350kev ',\
                                                                  ' 15_350kev_low ', \
                                                                   ' 15_350kev_hi '])
swift_models=pd.DataFrame(singular_swift_models_list, columns=['GRBname ', ' alpha ', \
                                                               ' alpha_low ', \
                                                               ' alpha_hi ', \
                                                              ' E_peak ', ' E_peak_low ', \
                                                              ' E_peak_hi '])
for i in range(0,len(swift_fluences_reference)):
    test_flux=swift_fluxes.at[i,' 15_350kev ']
    test_direction=swift_data['  RA_ground   '][i]
    test_time=swift_data['     T90      '][i]
    if test_direction==' N/A ' or test_direction=='N/A':
            tiny_removals.append(i)
    elif test_time==' N/A ' or test_time=='N/A':
            tiny_removals.append(i)
    elif test_flux==' N/A ':
            tiny_removals.append(i)
if (len(swift_fluxes)-1)<(len(swift_data)-1):
    for k in range(len(swift_fluxes)-1, len(swift_data)):
        tiny_removals.append(k)
sample_0_data = swift_data.drop(index=tiny_removals)
# print(len(sample_0_data))
# print(len(swift_fluences))
# print(len(swift_fluxes))
# print(len(swift_models))
# print(len(swift_data))

# %%
##Creating Samples 1 and 4##
removals=[]
xrt_matches=np.zeros((1,2))
xrt_stages=np.zeros((1,2))
xrt_avg_slope=np.zeros((1,2))
tarjete=0
stage=0
for i in range(0,len(swift_fluxes)):
    #for each event
    query=swift_data[' XRT_detection '][i]
    #query the XRT detection then add it to the first empty list
    #print(swift_fluxes.at[i,' 15_350kev '])
    test_flux=swift_fluxes.at[i,' 15_350kev ']
    test_direction=swift_data['  RA_ground   '][i]
    test_time=swift_data['     T90      '][i]
    #snap=swift_data['GRBname '][i].strip()[3:]
    #pdb.set_trace()
#     if swift_data['GRBname '][i]=='GRB080503    ':
#         pdb.set_trace()
    if len(np.where(swift_data['GRBname '][i].strip()[3:]==more_xrt_data['GRB '])[0])>0 \
    or len(np.where(swift_data['GRBname '][i][3:].strip()==weird_xrt_data['GRB '])[0])>0:
#         if swift_data['GRBname '][i].strip()[3:]=='211211A':
#             pdb.set_trace()
        if len(np.where(swift_data['GRBname '][i].strip()[3:]==more_xrt_data['GRB '])[0])>0:
            xrt_match=np.where(swift_data['GRBname '][i].strip()[3:]==\
                               more_xrt_data['GRB '])[0][0]
        elif len(np.where(swift_data['GRBname '][i].strip()[3:]==\
                          weird_xrt_data['GRB '])[0])>0:
            xrt_match=np.where(swift_data['GRBname '][i].strip()[3:]==
                               weird_xrt_data['GRB '])[0][0]
        #         xrt_trigger=xrt_data['Trigger Number'][xrt_match]
#         test_pindex=swift_pindex['GRB'][xrt_match]
#         xrt_flux=more_xrt_data[' Obs Flux_2 (pc) '][xrt_match]
        #they're off becuase I downloaded them at different times.
#         if xrt_flux==' N/A ' or xrt_flux=='N/A' or xrt_flux=='NaN':
#             removals.append(i)
        temporal_indices=np.zeros(6)
        for l in range(0, 6):
            flux_name=' alpha_{} '.format(l+1)
            if more_xrt_data[flux_name][xrt_match]==' N/A ' \
                or more_xrt_data[flux_name][xrt_match]=='N/A' \
                or more_xrt_data[flux_name][xrt_match]=='NaN':
                temporal_indices[l]=-100 #just throw it out of range
            else:
                temporal_indices[l]=more_xrt_data[flux_name][xrt_match]
        if len(np.where((temporal_indices > -0.75) & (temporal_indices < 0.75))[0])>0:
            tarjit=np.where((temporal_indices > -0.75) & (temporal_indices < 0.75))[0][0]
            flux_catch, xrt_flux, tarjete, state=defining_the_flux(tarjit, temporal_indices)
            slope_index=[]
            slope_sum=0
            for m in range(0, 6):
                if m==tarjit-1:
                    continue
                elif temporal_indices[m]==-100:
                    continue
                else:
                    slope_index=np.append(slope_index, temporal_indices[m])
            if len(slope_index)>0:
                slope_sum=np.average(slope_index)
            else:
                slope_sum=0
#         elif (np.where(abs(temporal_indices)==min(abs(temporal_indices)))[0][0]>0) and \
#         min(abs(temporal_indices))<100:
#             tarjit=np.where(abs(temporal_indices)==min(abs(temporal_indices)))[0][0]
#             flux_catch, xrt_flux, tarjete, state=defining_the_flux(tarjit, temporal_indices)
        else:
            xrt_flux=' N/A '
        if xrt_flux==' N/A ' or xrt_flux=='N/A' or xrt_flux=='NaN':
            removals.append(i)
#             #anyway, we don't wnat the ones without flux right now
# #         elif test_pindex==' N/A ' or test_pindex=='N/A' or test_pindex=='NaN':
# #             removals.append(i)
# #         elif query!='Yes' and query!=' Yes ':
# #             #and add it to the second one if it didn't detect or it's not sure
# #             removals.append(i)
        if test_flux==' N/A ':
            removals.append(i)
        elif test_direction==' N/A ' or test_direction=='N/A':
            removals.append(i)
        elif test_time==' N/A ' or test_time=='N/A' or test_time=='     N/A      ':
            removals.append(i)
#         elif xrt_trigger=='BAT/GUANO':
            #These ones probably won't show up in the the data, but I don't want them
            #because they aren't relevant to what I'm doing, being that they /probably/ 
            #aren't triggered by Fermi and the data ends up missing everything relevant
#             removals.append(i)
        else:
            xrt_matches=np.append(xrt_matches,[[int(i),int(xrt_match)]], axis=0)
#             xrt_stages=np.append(xrt_stages,[[int(state),int(tarjete)]], axis=0)
#             xrt_avg_slope=np.append(xrt_avg_slope,[[int(i),int(slope_sum)]], axis=0)
            continue
    else:
        removals.append(i)
if (len(swift_fluxes)-1)<(len(swift_data)-1):
    for k in range(len(swift_fluxes)-1, len(swift_data)):
        removals.append(k)
edited_swift_data = swift_data.drop(index=removals)
xrt_matches=xrt_matches[1:]
# xrt_stages=xrt_stages[1:]
# xrt_avg_slope=xrt_avg_slope[1:]
test_flux=[]
for j in edited_swift_data.index:
    test_flux.append(swift_fluxes.at[j,' 15_350kev '])
red_swift_removals=[]
r=0
redshift_matches=np.zeros((1,2))
redshift_xrt=np.zeros((1,2))
# redshift_stages=np.zeros((1,2))
# redshift_avg_slope=np.zeros((1,2))
for j in edited_swift_data.index:
    for i in range(0, len(swift_redshifts)):
        swift_string=edited_swift_data['GRBname '][j].strip()
        redshift_string=swift_redshifts['GRBname '][i].strip()
        redsweefer = swift_redshifts.at[i, ' z '] #I call a swiffer a sweefer and 
        #this error was going to be a serious issue unlike the stupid question marks
        if swift_string==redshift_string:
            redshift_matches=np.append(redshift_matches,[[int(i),int(j)]], axis=0)
            row_match=int(np.where(xrt_matches[:, 0]==j)[0][0])
            redshift_xrt=np.append(redshift_xrt,[[int(j),xrt_matches[row_match,1]]], \
                                       axis=0)
#             redshift_stages=np.append(redshift_stages,[[int(xrt_stages[row_match,0]),
#                                                         int(xrt_stages[row_match,1])]], \
#                                       axis=0)
#             redshift_avg_slope=np.append(redshift_avg_slope,[[int(i),\
#                                                     int(xrt_avg_slope[row_match, 1])]],\
#                                          axis=0)
            if 'or' in str(redsweefer):
#                 print(redsweefer)
#                 print(redshift_string)
                red_swift_removals.append(j)
                redshift_xrt=redshift_xrt[:len(redshift_xrt)-1,:]
#                 redshift_stages=redshift_stages[:len(redshift_stages)-1,:]
#                 redshift_avg_slope=redshift_avg_slope[:len(redshift_stages)-1,:]
            if '-' in str(redsweefer):
                #whatever jerk did this one is okay in my book, not like ALL THE OTHER AWFUL
                #NONNUMERIC CHARACTERS AHHHHHHHHH, jerks, all the rest of you
                redswiffer=re.split("-", redsweefer.strip())
                red_broom = list(map(float, redswiffer))
                swift_redshifts.at[i, ' z ']=np.mean([red_broom[0], red_broom[1]])
            break
        elif swift_string != redshift_string and i==len(swift_redshifts)-1:
            r=r+1
            red_swift_removals.append(j)
        else: 
            continue
other_redshift_matches=np.zeros((1,2))
for j in edited_swift_data.index:
    for i in range(0, len(second_redshifts)):
        swift_string_2=edited_swift_data['GRBname '][j].strip()
        redshift_string_2=second_redshifts['GRB'][i].strip()
        redmop = second_redshifts.at[i, 'zc']
        if swift_string_2==redshift_string_2:
            other_redshift_matches=np.append(other_redshift_matches,[[int(i),int(j)]], axis=0)
            if type(red_swift_removals)!=int:
                if len(np.where(red_swift_removals==j)[0])>0:
                    add_in=np.where(red_swift_removals==j)[0][0]
                    red_swift_removals=np.delete(red_swift_removals,add_in)
            elif type(red_swift_removals)==int and red_swift_removals==j:
                red_swift_removals=[]
            #I /think/ this will add the removed ones back in. Probably.
            row_match=int(np.where(xrt_matches[:, 0]==j)[0][0])
            redshift_xrt=np.append(redshift_xrt,[[int(j),xrt_matches[row_match,1]]],\
                                   axis=0)
#             redshift_stages=np.append(redshift_stages,[[int(xrt_stages[row_match,0]),\
#                                                         int(xrt_stages[row_match,1])]],\
#                                       axis=0)
            break
#         elif swift_string_2 != redshift_string_2 and i==len(second_redshifts)-1:
#             red_swift_removals=np.append(red_swift_removals, j)
#             break
        else: 
            continue
# weird_removals=list(set(red_swift_removals))
redshift_matches=redshift_matches[1:]
redshift_xrt=redshift_xrt[1:]
# redshift_stages=redshift_stages[1:]
# redshift_avg_slope=redshift_avg_slope[1:]
redshift_swift_data = edited_swift_data.drop(index=red_swift_removals)
# print(len(edited_swift_data))
# print(len(redshift_swift_data))
edited_swift_data = swift_data.drop(index=removals)
sample_1_data=edited_swift_data
#For some reason I won't pretend to understand, creating sample 4 HERE prevents errors in 
#creating sample 3. Oy vey.
sample_4_data=redshift_swift_data
# print(len(sample_1_data))
# print(len(sample_4_data))

# %%
##Creating Sample 2##
sample_2_data=fermi_data
# print(len(sample_2_data))

# %%
##Creating Sample 3##
##Data conversion##
edited_indices=list(np.where(edited_swift_data.index!=0))
t=Time(edited_swift_data.at[int(edited_swift_data.index[0]),\
                            '       Trig_time_UTC        '].strip(), \
       format='isot', scale='utc')
gorp=t.jd; #honestly at this point I gave up on naming variables
shoes=float(edited_swift_data.at[int(edited_swift_data.index[0]),'  RA_ground   '])
sorks=float(edited_swift_data.at[int(edited_swift_data.index[0]),'     T90      '])
gatorade=float(edited_swift_data.at[int(edited_swift_data.index[0]), ' Image_position_err '])
helmet=float(edited_swift_data.at[int(edited_swift_data.index[0]),'  DEC_ground   '])
nom=str(edited_swift_data.at[int(edited_swift_data.index[0]),'GRBname '].strip())
for i in list(edited_swift_data.index):
    if type(edited_swift_data.at[i,'       Trig_time_UTC        '])==str:
        t=Time(edited_swift_data.at[i,'       Trig_time_UTC        '].strip(), \
               format='isot', scale='utc')
        gorp=np.append(gorp, t.jd)
        shoes=np.append(shoes, float(edited_swift_data.at[i,'  RA_ground   ']))
        sorks=np.append(sorks, float(edited_swift_data.at[i,'     T90      ']))
        gatorade=np.append(gatorade, float(edited_swift_data.at[i,' Image_position_err ']))
        helmet=np.append(helmet, float(edited_swift_data.at[i,'  DEC_ground   ']))
        nom=np.append(nom, str(edited_swift_data.at[i,'GRBname '].strip()))
    else :
        gorp=np.append(gorp, 1)
        sorks=np.append(sorks, float(edited_swift_data.at[i,'     T90      ']))
        shoes=np.append(shoes, float(edited_swift_data.at[i,'  RA_ground   ']))
        gatorade=np.append(gatorade, float(edited_swift_data.at[i,' Image_position_err ']))
        helmet=np.append(helmet, float(edited_swift_data.at[i,'  DEC_ground   ']))
        nom=np.append(nom, str(edited_swift_data.at[i,'GRBname '].strip()))
t0=fermi_data.at[0,"trigger_time           "]
t0=re.split("/| |:", t0)
t0 = list(map(int, t0))
t0=datetime.datetime(t0[2], t0[0], t0[1], t0[3], t0[4], t0[5], int(t0[6]))
trail_mix = julian.to_jd(t0, fmt='jd') #honestly at this point I gave up on naming variables
soles=fermi_data.at[0,'ra        ']
soles=soles.split(' ')
soles=list(map(float, soles)) 
boots=((soles[0]+(soles[1]+soles[2]/60)/60)*15)
fabric=fermi_data.at[0,'dec      ']
fabric=fabric.split(' ')
fabric=list(map(float, fabric)) 
hat=(fabric[0]+(fabric[1]+fabric[2]/60)/60)
socks=fermi_data.at[0,'t90_error']
w0ter=fermi_data.at[0, 'error_radius']
for i in range(1,len(fermi_data)):
    ti83=fermi_data.at[i,"trigger_time           "]
    ti83=re.split("/| |:", ti83)
    ti83 = list(map(int, ti83)) 
    ti83=datetime.datetime(ti83[2], ti83[0], ti83[1], ti83[3], ti83[4], ti83[5], \
                           int(ti83[6]))
    trail_mix = np.append(trail_mix, julian.to_jd(ti83, fmt='jd'))
    soles=fermi_data.at[i,'ra        ']
    soles=soles.split(' ')
    soles=list(map(float, soles)) 
    soles=((soles[0]+(soles[1]+soles[2]/60)/60)*15)
    boots=np.append(boots, soles)
    fabric=fermi_data.at[i,'dec      ']
    fabric=fabric.split(' ')
    fabric=list(map(float, fabric)) 
    fabric=(fabric[0]+(fabric[1]+fabric[2]/60)/60)
    hat=np.append(hat, fabric)
    socks=np.append(socks, fermi_data.at[i,'t90_error'])
    w0ter=np.append(w0ter, fermi_data.at[i,'error_radius'])
##match data##
food = np.zeros((len(trail_mix),len(gorp)))
feet_covers = np.zeros((len(trail_mix),len(gorp)))
UV_protection = np.zeros((len(trail_mix),len(gorp)))
k=0;
for i in range(0, len(gorp)):
    for j in range(0, len(trail_mix)):
        food[j, i] = abs(trail_mix[j]-gorp[i]);
        feet_covers[j, i] = abs(boots[j]-shoes[i]);
        UV_protection[j, i] = abs(hat[j]-helmet[i]);
        if food[j,i]<0.0035:
            if w0ter[j]<0.1:
                if feet_covers[j,i]<3*gatorade[i]:
                    if UV_protection[j,i]<3*gatorade[i]:
                        k=k+1;
            else:
                if feet_covers[j,i]<3*w0ter[j]:
                    if UV_protection[j,i]<3*w0ter[j]:
                        k=k+1;
matches = np.zeros((k+1, 15))
potential_matches=np.zeros((k,3))
# print(k)
k=0;
err=0;
special_cases=[]
for i in range(0, len(gorp)):
    for j in range(0, len(trail_mix)):
        if food[j,i]<0.0035: #time
            if w0ter[j]<0.1:
                if feet_covers[j,i]<3*gatorade[i]: #RA
                    if UV_protection[j,i]<3*gatorade[i]: #dec
                        potential_matches[k,0]=food[j,i]
                        potential_matches[k,1]=w0ter[j]
                        potential_matches[k,2]=gatorade[i]
                        name=fermi_data.at[j,"name        "]
                        matches[k,0] = int(name[3]+name[4]+name[5]+name[6]+\
                                           name[7]+name[8]+name[9]+name[10]+name[11])
                        matches[k,1] = j
                        matches[k,2] = edited_swift_data.index[i]
                        if w0ter[j] != 0 and gatorade[i] != 0:
                            matches[k,3] = 1-NormalDist(mu=boots[j], sigma=w0ter[j]).overlap(\
                                        NormalDist(mu=shoes[i], sigma=gatorade[i]))*\
                                        NormalDist(mu=hat[j], sigma=w0ter[j]).overlap(\
                                        NormalDist(mu=helmet[i], sigma=gatorade[i])) #alpha
                            matches[k,4] = NormalDist(mu=boots[j], sigma=w0ter[j]).overlap(\
                                        NormalDist(mu=shoes[i], sigma=gatorade[i]))*\
                                        NormalDist(mu=hat[j], sigma=w0ter[j]).overlap(\
                                        NormalDist(mu=helmet[i], sigma=gatorade[i])) #beta
                        else:
                            matches[k,3]='NaN'
                            matches[k,4]='NaN'
                            err=err+1;
                        matches[k,5] = fermi_data.at[j,"t90     "]
                        matches[k,6] = fermi_data.at[j,"t90_error"]
                        matches[k,7] = fermi_data.at[j,'flux_256 ']
                        matches[k,8] = fermi_data.at[j,'fluence   ']
                        matches[k,9] = fermi_data.at[j,'fluence_error']
                        if swift_fluxes.at[edited_swift_data.index[i],\
                                                           ' 15_350kev '] != ' N/A ':
                            matches[k,10] = swift_fluxes.at[edited_swift_data.index[i],\
                                                               ' 15_350kev ']
                        else:
                            matches[k,10] = 0
                        if swift_fluences.at[edited_swift_data.index[i],\
                                                             ' 15_350kev '] != ' N/A ':
                            matches[k,11] =  swift_fluences.at[edited_swift_data.index[i],\
                                                             ' 15_350kev ']
                        else:
                            matches[k,11] = 0
                        if fermi_data.at[j,'flnc_best_fitting_model']==\
                        'flnc_plaw              ':
                            matches[k,12] = 1
                        elif fermi_data.at[j,'flnc_best_fitting_model']==\
                        'flnc_band              ':
                            matches[k,12] = 2
                        elif fermi_data.at[j,'flnc_best_fitting_model']==\
                        'flnc_comp              ':
                            matches[k,12] = 3
                        elif fermi_data.at[j,'flnc_best_fitting_model']==\
                        'flnc_sbpl              ':
                            matches[k,12] = 4
                        else:
                            matches[k,12] = 0
                        #pdb.set_trace()
                        xrt_xcall=np.where(xrt_matches[:,0]==\
                                                 edited_swift_data.index[i-1])[0][0]
                        matches[k,13]=int(xrt_matches[xrt_xcall, 1])
#                         matches[k,14] = int(xrt_stages[i-1,1])
#     #                     matches[k, 11]=xrt_data.at[xrt_call,\
#     #                                         'XRT Early Flux (0.3-10 keV) [10^-11 erg/cm^2/s]']
#     #                     print(i-1)
#     #                     print(xrt_call)
#                         if int(xrt_stages[i-1,0])==0:
#                             flux_recall=' Obs Flux_{} (pc) '.\
#                                 format(int(xrt_stages[i-1,1]))
#                             flux_change=' D_Obs Flux_{} (pc) '.\
#                                     format(int(xrt_stages[i-1,1]))
#                         elif int(xrt_stages[i-1,0])==1:
#                             flux_recall=' Obs Flux_{} (wt) '.\
#                                 format(int(xrt_stages[i-1,1]))
#                             flux_change=' D_Obs Flux_{} (wt) '.\
#                                     format(int(xrt_stages[i-1,1]))
#                         else:
#                             print("error, this burst ({}) didn't get an assigned flux.").\
#                             format(fermi_data.at[j,"name        "])
#     #                     print(fermi_data.at[j,"name        "])
#     #                     print(flux_recall)
#     #                     print(more_xrt_data[flux_recall][xrt_call])
#                         matches[k, 13]=more_xrt_data[flux_recall][xrt_call]
#                         if more_xrt_data[flux_change][xrt_call]==' , ':
#                             matches[k,14]=0.34*float(more_xrt_data[flux_recall][xrt_call])
#                             matches[k,15]=0.34*float(more_xrt_data[flux_recall][xrt_call])
#                             #this is about one standard deviation
#                         elif ',' in str(more_xrt_data[flux_change][xrt_call]):
#                             dual_input=re.split(",", str(more_xrt_data[flux_change][xrt_call])\
#                                                 .strip())
#                             headphones = list(map(float, dual_input))
#                             matches[k,14]=headphones[0]
#                             matches[k,15]=headphones[1]
#                             #I originally only left room for one value. That was a mistake.
#                         else:
#                             matches[k, 14]=more_xrt_data[flux_change][xrt_call]
#                             matches[k, 15]=more_xrt_data[flux_change][xrt_call]
#                         matches[k,16]=xrt_avg_slope[xrt_xcall, 1]
                        k=k+1
            else:
                if feet_covers[j,i]<3*w0ter[j]: #RA
                    if UV_protection[j,i]<3*w0ter[j]: #dec
                        special_cases=np.append(special_cases,fermi_data.at[j,"name        "])
            
                        potential_matches[k,0]=food[j,i]
                        potential_matches[k,1]=w0ter[j]
                        potential_matches[k,2]=gatorade[i]
                        name=fermi_data.at[j,"name        "]
                        matches[k,0] = int(name[3]+name[4]+name[5]+name[6]+\
                                           name[7]+name[8]+name[9]+name[10]+name[11])
                        matches[k,1] = j
                        matches[k,2] = edited_swift_data.index[i]
                        if w0ter[j] != 0 and gatorade[i] != 0:
                            matches[k,3] = 1-NormalDist(mu=boots[j], sigma=w0ter[j]).overlap(\
                                        NormalDist(mu=shoes[i], sigma=gatorade[i]))*\
                                        NormalDist(mu=hat[j], sigma=w0ter[j]).overlap(\
                                        NormalDist(mu=helmet[i], sigma=gatorade[i])) #alpha
                            matches[k,4] = NormalDist(mu=boots[j], sigma=w0ter[j]).overlap(\
                                        NormalDist(mu=shoes[i], sigma=gatorade[i]))*\
                                        NormalDist(mu=hat[j], sigma=w0ter[j]).overlap(\
                                        NormalDist(mu=helmet[i], sigma=gatorade[i])) #beta
                        else:
                            matches[k,3]='NaN'
                            matches[k,4]='NaN'
                            err=err+1;
                        matches[k,5] = fermi_data.at[j,"t90     "]
                        matches[k,6] = fermi_data.at[j,"t90_error"]
                        matches[k,7] = fermi_data.at[j,'flux_256 ']
                        matches[k,8] = fermi_data.at[j,'fluence   ']
                        matches[k,9] = fermi_data.at[j,'fluence_error']
                        if swift_fluxes.at[edited_swift_data.index[i],\
                                                           ' 15_350kev '] != ' N/A ':
                            matches[k,10] = swift_fluxes.at[edited_swift_data.index[i],\
                                                               ' 15_350kev ']
                        else:
                            matches[k,10] = 0
                        if swift_fluences.at[edited_swift_data.index[i],\
                                                             ' 15_350kev '] != ' N/A ':
                            matches[k,11] = swift_fluences.at[edited_swift_data.index[i],\
                                                             ' 15_350kev ']
                        else:
                            matches[k,11] = 0
                        if fermi_data.at[j,'flnc_best_fitting_model']==\
                        'flnc_plaw              ':
                            matches[k,12] = 1
                        elif fermi_data.at[j,'flnc_best_fitting_model']==\
                        'flnc_band              ':
                            matches[k,12] = 2
                        elif fermi_data.at[j,'flnc_best_fitting_model']==\
                        'flnc_comp              ':
                            matches[k,12] = 3
                        elif fermi_data.at[j,'flnc_best_fitting_model']==\
                        'flnc_sbpl              ':
                            matches[k,12] = 4
                        else:
                            matches[k,12] = 0
                        xrt_xcall=np.where(xrt_matches[:,0]==\
                                                 edited_swift_data.index[i-1])[0][0]
                        matches[k,13]=int(xrt_matches[xrt_xcall, 1])
#                         matches[k,14] = int(xrt_stages[i-1,1])
#                         if int(xrt_stages[i-1,0])==0:
#                             flux_recall=' Obs Flux_{} (pc) '.\
#                                 format(int(xrt_stages[i-1,1]))
#                             flux_change=' D_Obs Flux_{} (pc) '.\
#                                     format(int(xrt_stages[i-1,1]))
#                         elif int(xrt_stages[i-1,0])==1:
#                             flux_recall=' Obs Flux_{} (wt) '.\
#                                 format(int(xrt_stages[i-1,1]))
#                             flux_change=' D_Obs Flux_{} (wt) '.\
#                                     format(int(xrt_stages[i-1,1]))
#                         else:
#                             print("error, this burst ({}) didn't get an assigned flux.").\
#                             format(fermi_data.at[j,"name        "])
#                         matches[k, 13]=more_xrt_data[flux_recall][xrt_call]
#                         if more_xrt_data[flux_change][xrt_call]==' , ':
#                             matches[k,14]=0.34*float(more_xrt_data[flux_recall][xrt_call])
#                             matches[k,15]=0.34*float(more_xrt_data[flux_recall][xrt_call])
#                             #this is about one standard deviation
#                         elif ',' in str(more_xrt_data[flux_change][xrt_call]):
#                             dual_input=re.split(",", str(more_xrt_data[flux_change][xrt_call])\
#                                                 .strip())
#                             headphones = list(map(float, dual_input))
#                             matches[k,14]=headphones[0]
#                             matches[k,15]=headphones[1]
#                             #I originally only left room for one value. That was a mistake.
#                         else:
#                             matches[k, 14]=more_xrt_data[flux_change][xrt_call]
#                             matches[k, 15]=more_xrt_data[flux_change][xrt_call]
#                         matches[k,16]=xrt_avg_slope[xrt_xcall, 1]
                        k=k+1
regular_matches=pd.DataFrame(matches, columns=['name', 'Fermi row', 'Swift row',\
                                               'Type I error', \
                                       'Type II error', 't90', 't90_error', 'GBM flux', \
                                       'GBM fluence', 'Fluence error',\
                                               'BAT flux', 'BAT fluence',\
                                               'Spectral Model', 'XRT row', \
                                               'Plateau Stage'])#,\
#                                                '"Plateau" X-Ray Flux', \
#                                                "Pre-X-ray Flux Change", \
#                                                "Post-X-ray Flux Change"])
regular_matches=regular_matches[:-1]
sample_3_data=regular_matches
# print(special_cases)

# %%
# ##Sample 3 Diagnositics##
unknown_progenitor_HR=[]
unknown_progenitor_HR_err=np.zeros((1,2))
unknown_progenitor_t90=[]
unknown_progenitor_t90_err=[]
unknown_progenitor_Epeak_over_s=[]
unknown_progenitor_Epeak_over_s_err=[]
unknown_progenitor_Epeak=[]
unknown_progenitor_Epeak_err=[]
unknown_progenitor_flu=[]
unknown_progenitor_flu_err=[]
unknown_progenitor_ag_flu=[]
unknown_progenitor_ag_flu_err=[]
long_merger_HR=[]
long_merger_HR_err=np.zeros((1,2))
long_merger_t90=[]
long_merger_t90_err=[]
long_merger_Epeak_over_s=[]
long_merger_Epeak_over_s_err=[]
long_merger_Epeak=[]
long_merger_Epeak_err=[]
long_merger_flu=[]
long_merger_flu_err=[]
long_merger_ag_flu=[]
long_merger_ag_flu_err=[]
short_merger_HR=[]
short_merger_HR_err=np.zeros((1,2))
swift_short_merger_t90=[]
swift_short_merger_t90_err=[]
short_merger_t90=[]
short_merger_t90_err=[]
short_merger_Epeak_over_s=[]
short_merger_Epeak_over_s_err=[]
short_merger_Epeak=[]
short_merger_Epeak_err=[]
short_merger_flu=[]
short_merger_flu_err=[]
short_merger_ag_flu=[]
short_merger_ag_flu_err=[]
long_collapsar_HR=[]
long_collapsar_HR_err=np.zeros((1,2))
long_collapsar_t90=[]
long_collapsar_t90_err=[]
long_collapsar_Epeak_over_s=[]
long_collapsar_Epeak_over_s_err=[]
long_collapsar_Epeak=[]
long_collapsar_Epeak_err=[]
long_collapsar_flu=[]
long_collapsar_flu_err=[]
long_collapsar_ag_flu=[]
long_collapsar_ag_flu_err=[]
short_collapsar_HR=[]
short_collapsar_HR_err=np.zeros((1,2))
short_collapsar_t90=[]
short_collapsar_t90_err=[]
short_collapsar_Epeak_over_s=[]
short_collapsar_Epeak_over_s_err=[]
short_collapsar_Epeak=[]
short_collapsar_Epeak_err=[]
short_collapsar_flu=[]
short_collapsar_flu_err=[]
short_collapsar_ag_flu=[]
short_collapsar_ag_flu_err=[]
nbins=51
p=0
spectral_mask=np.isin(sample_3_data['Spectral Model'], 0.0)
spectral_sample_3_data=sample_3_data[~spectral_mask]
Refined_Goldstein_Data=np.zeros((len(spectral_sample_3_data),12))
lc_mask=np.isin(spectral_sample_3_data['name'], \
                Overlap_Precursors["Long Collapsars"])
# lc_sample_3_data=sample_3_data[lc_mask]
spectral_lc_sample_3_data=spectral_sample_3_data[lc_mask]
# print(lc_sample_3_data.columns)
for i in spectral_lc_sample_3_data.index.to_list():
    q=spectral_lc_sample_3_data.at[i, 'Fermi row']
    long_collapsar_HR=np.append(long_collapsar_HR, Fermi_HR_func(sample_2_data, \
                                sample_2_data.at[q, 'flnc_best_fitting_model'], q))
    long_collapsar_t90=np.append(long_collapsar_t90, float(sample_2_data.at[q, 't90     ']))
    long_collapsar_Epeak_over_s=np.append(long_collapsar_Epeak_over_s, \
                        (getting_the_GBM_E_peak(sample_2_data, \
                                    sample_2_data.at[q, 'flnc_best_fitting_model'], q))/\
                                       float(sample_2_data.at[q, 'fluence   ']))
    long_collapsar_Epeak_over_s_err=np.append(long_collapsar_Epeak_over_s_err, \
                        (getting_the_GBM_E_peak_err(sample_2_data, \
                                    sample_2_data.at[q, 'flnc_best_fitting_model'], q))/\
                                       float(sample_2_data.at[i, 'fluence   ']))
    long_collapsar_Epeak=np.append(long_collapsar_Epeak, getting_the_GBM_E_peak(sample_2_data, \
                                    sample_3_data.at[i, 'Spectral Model'], \
                                               sample_3_data.at[i, 'Fermi row']))
    long_collapsar_flu=np.append(long_collapsar_flu, float(sample_3_data.at[i, 'GBM fluence']))
    long_collapsar_Epeak_err=np.append(long_collapsar_Epeak, \
                                       getting_the_GBM_E_peak_err(sample_2_data, \
                                    sample_3_data.at[i, 'Spectral Model'], \
                                               sample_3_data.at[i, 'Fermi row']))
    long_collapsar_flu_err=np.append(long_collapsar_flu_err, \
                                     float(sample_3_data.at[i, 'Fluence error']))
    long_collapsar_ag_flu=np.append(long_collapsar_ag_flu, \
            afterglow_fluence_finder(spectral_lc_sample_3_data, more_xrt_data, i, \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Time"], \
            #heck if i remember why it's i-1 and I just wrote this 10 seconds ago. 
            #I discovered it by trial and error
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Error"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Error"]))
    long_collapsar_ag_flu_err=np.append(long_collapsar_ag_flu_err, \
             afterglow_fluence_err_finder(spectral_lc_sample_3_data, more_xrt_data, i, \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Error"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Error"]))
    long_collapsar_t90_err=np.append(long_collapsar_t90_err, float(sample_2_data.at[i,\
                                                                        't90_error']))
    Refined_Goldstein_Data[p,0]=sample_2_data.at[q, 'name        '][3:]
    Refined_Goldstein_Data[p,1]=long_collapsar_Epeak_over_s[-1]
    Refined_Goldstein_Data[p,2]=long_collapsar_Epeak_over_s_err[-1]
    Refined_Goldstein_Data[p,3]=long_collapsar_HR[-1]
    Refined_Goldstein_Data[p,4]=long_collapsar_t90[-1]
    Refined_Goldstein_Data[p,5]=long_collapsar_t90_err[-1]
    Refined_Goldstein_Data[p,6]=long_collapsar_Epeak[-1]
    Refined_Goldstein_Data[p,7]=long_collapsar_Epeak_err[-1]
    Refined_Goldstein_Data[p,8]=long_collapsar_flu[-1]
    Refined_Goldstein_Data[p,9]=long_collapsar_flu_err[-1]
    Refined_Goldstein_Data[p,10]=long_collapsar_ag_flu[-1]
    Refined_Goldstein_Data[p,11]=long_collapsar_ag_flu_err[-1]
    p=p+1
#     print('The Fermi name is {} and the BAT name is {}'.\
#         format(fermi_data.at[q, 'name        '], \
#         Afterglow_of_Sample_3_Data.loc[i-1, 'Name']))
# print(p)
long_collapsar_t90_err=np.where(spectral_lc_sample_3_data['t90_error'] != 'N/A', \
        spectral_lc_sample_3_data['t90_error'], 0)
if len(long_collapsar_t90_err)>0:
    for j in range(0, len(long_collapsar_t90_err)):
        long_collapsar_t90_err[j]=float(long_collapsar_t90_err[j])
        
sc_mask=np.isin(spectral_sample_3_data['name'], \
                Overlap_Precursors["Short Collapsars"])
spectral_sc_sample_3_data=spectral_sample_3_data[sc_mask]
for i in spectral_sc_sample_3_data.index.to_list():
    q=spectral_sample_3_data.at[i, 'Fermi row']
    short_collapsar_HR=np.append(short_collapsar_HR, Fermi_HR_func(sample_2_data, \
                                sample_2_data.at[q, 'flnc_best_fitting_model'], q))
    short_collapsar_t90=np.append(short_collapsar_t90, \
                                  float(sample_2_data.at[q, 't90     ']))
    short_collapsar_Epeak_over_s=np.append(short_collapsar_Epeak_over_s, \
                        (getting_the_GBM_E_peak(sample_2_data, \
                                sample_2_data.at[q, 'flnc_best_fitting_model'], q))/\
                                       float(sample_2_data.at[q, 'fluence   ']))
    short_collapsar_Epeak=np.append(short_collapsar_Epeak, getting_the_GBM_E_peak(sample_2_data, \
                                    sample_3_data.at[i, 'Spectral Model'], \
                                               sample_3_data.at[i, 'Fermi row']))
    short_collapsar_flu=np.append(short_collapsar_flu, float(sample_3_data.at[i, 'GBM fluence']))
    short_collapsar_Epeak_over_s_err=np.append(short_collapsar_Epeak_over_s, \
                        (getting_the_GBM_E_peak_err(sample_2_data, \
                                sample_2_data.at[q, 'flnc_best_fitting_model'], q))/\
                                       float(sample_2_data.at[q, 'fluence_error']))
    short_collapsar_Epeak_err=np.append(short_collapsar_Epeak, \
                                       getting_the_GBM_E_peak_err(sample_2_data, \
                                    sample_3_data.at[i, 'Spectral Model'], \
                                               sample_3_data.at[i, 'Fermi row']))
    short_collapsar_flu_err=np.append(short_collapsar_flu_err, \
                                     float(sample_3_data.at[i, 'Fluence error']))
    short_collapsar_t90_err=np.append(short_collapsar_t90_err, \
                                  float(sample_2_data.at[q, 't90_error']))
    short_collapsar_ag_flu=np.append(short_collapsar_ag_flu, \
            afterglow_fluence_finder(spectral_sc_sample_3_data, more_xrt_data, i, \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Time"], \
            #heck if i remember why it's i-1 and I just wrote this 10 seconds ago. 
            #I discovered it by trial and error
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Error"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Error"]))
    short_collapsar_ag_flu_err=np.append(short_collapsar_ag_flu_err, \
             afterglow_fluence_err_finder(spectral_sc_sample_3_data, more_xrt_data, i, \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Error"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Error"]))
    Refined_Goldstein_Data[p,0]=sample_2_data.at[q, 'name        '][3:]
    Refined_Goldstein_Data[p,1]=short_collapsar_Epeak_over_s[-1]
    Refined_Goldstein_Data[p,2]=short_collapsar_Epeak_over_s_err[-1]
    Refined_Goldstein_Data[p,3]=short_collapsar_HR[-1]
    Refined_Goldstein_Data[p,4]=short_collapsar_t90[-1]
    Refined_Goldstein_Data[p,5]=short_collapsar_t90_err[-1]
    Refined_Goldstein_Data[p,6]=short_collapsar_Epeak[-1]
    Refined_Goldstein_Data[p,7]=short_collapsar_Epeak_err[-1]
    Refined_Goldstein_Data[p,8]=short_collapsar_flu[-1]
    Refined_Goldstein_Data[p,9]=short_collapsar_flu_err[-1]
    Refined_Goldstein_Data[p,10]=short_collapsar_ag_flu[-1]
    Refined_Goldstein_Data[p,11]=short_collapsar_ag_flu_err[-1]
    p=p+1
# print(p)
short_collapsar_t90_err=np.where(spectral_sc_sample_3_data['t90_error'] != 'N/A', \
        spectral_sc_sample_3_data['t90_error'], 0)
if len(short_collapsar_t90_err)>0:
    for j in range(0, len(short_collapsar_t90_err)):
        short_collapsar_t90_err[j]=float(short_collapsar_t90_err[j])
        
lm_mask=np.isin(spectral_sample_3_data['name'], \
                Overlap_Precursors["Long Mergers"])
spectral_lm_sample_3_data=spectral_sample_3_data[lm_mask]
for i in spectral_lm_sample_3_data.index.to_list():
    q=spectral_sample_3_data.at[i, 'Fermi row']
    long_merger_HR=np.append(long_merger_HR, Fermi_HR_func(sample_2_data, \
                                sample_2_data.at[q, 'flnc_best_fitting_model'], q))
    long_merger_t90=np.append(long_merger_t90, float(sample_2_data.at[q, 't90     ']))
    long_merger_Epeak_over_s=np.append(long_merger_Epeak_over_s, \
                        (getting_the_GBM_E_peak(sample_2_data, \
                                sample_2_data.at[q, 'flnc_best_fitting_model'], q))/\
                                       float(sample_2_data.at[q, 'fluence   ']))
    long_merger_Epeak=np.append(long_merger_Epeak, getting_the_GBM_E_peak(sample_2_data, \
                                    sample_3_data.at[i, 'Spectral Model'], \
                                    sample_3_data.at[i, 'Fermi row']))
    long_merger_flu=np.append(long_merger_flu, float(sample_3_data.at[i, 'GBM fluence']))
    long_merger_Epeak_over_s_err=np.append(long_merger_Epeak_over_s_err, \
                        (getting_the_GBM_E_peak(sample_2_data, \
                                sample_2_data.at[q, 'flnc_best_fitting_model'], q))/\
                                       float(sample_2_data.at[q, 'fluence   ']))
    long_merger_Epeak_err=np.append(long_merger_Epeak, \
                                       getting_the_GBM_E_peak_err(sample_2_data, \
                                    sample_3_data.at[i, 'Spectral Model'], \
                                               sample_3_data.at[i, 'Fermi row']))
    long_merger_flu_err=np.append(long_merger_flu_err, \
                                     float(sample_3_data.at[i, 'Fluence error']))
    long_merger_ag_flu=np.append(long_merger_ag_flu, \
            afterglow_fluence_finder(spectral_lm_sample_3_data, more_xrt_data, i, \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Time"], \
            #heck if i remember why it's i-1 and I just wrote this 10 seconds ago. 
            #I discovered it by trial and error
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Error"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Error"]))
    long_merger_ag_flu_err=np.append(long_merger_ag_flu_err, \
             afterglow_fluence_err_finder(spectral_lm_sample_3_data, more_xrt_data, i, \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Error"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Error"]))
    long_merger_t90_err=np.append(long_merger_t90_err, float(sample_2_data.at[q, 't90_error']))
    Refined_Goldstein_Data[p,0]=sample_2_data.at[q, 'name        '][3:]
    Refined_Goldstein_Data[p,1]=long_merger_Epeak_over_s[-1]
    Refined_Goldstein_Data[p,2]=long_merger_Epeak_over_s_err[-1]
    Refined_Goldstein_Data[p,3]=long_merger_HR[-1]
    Refined_Goldstein_Data[p,4]=long_merger_t90[-1]
    Refined_Goldstein_Data[p,5]=long_merger_t90_err[-1]
    Refined_Goldstein_Data[p,6]=long_merger_Epeak[-1]
    Refined_Goldstein_Data[p,7]=long_merger_Epeak_err[-1]
    Refined_Goldstein_Data[p,8]=long_merger_flu[-1]
    Refined_Goldstein_Data[p,9]=long_merger_flu_err[-1]
    Refined_Goldstein_Data[p,10]=long_merger_ag_flu[-1]
    Refined_Goldstein_Data[p,11]=long_merger_ag_flu_err[-1]
    p=p+1
# print(p)
long_merger_t90_err=np.where(spectral_lm_sample_3_data['t90_error'] != 'N/A', \
        spectral_lm_sample_3_data['t90_error'], 0)
if len(long_merger_t90_err)>0:
    for j in range(0, len(long_merger_t90_err)):
        long_merger_t90_err[j]=float(long_merger_t90_err[j])

sm_mask=np.isin(spectral_sample_3_data['name'], \
                Overlap_Precursors["Short Mergers"])
spectral_sm_sample_3_data=spectral_sample_3_data[sm_mask]
for i in spectral_sm_sample_3_data.index.to_list():
    q=spectral_sample_3_data.at[i, 'Fermi row']
    short_merger_HR=np.append(short_merger_HR, Fermi_HR_func(sample_2_data, \
                                sample_2_data.at[q, 'flnc_best_fitting_model'], q))
    short_merger_t90=np.append(short_merger_t90, float(sample_2_data.at[q, 't90     ']))
    short_merger_Epeak_over_s=np.append(short_merger_Epeak_over_s, \
                       (getting_the_GBM_E_peak(sample_2_data, 
                            sample_2_data.at[q, 'flnc_best_fitting_model'], q))/\
                                       float(sample_2_data.at[q, 'fluence   ']))
    short_merger_Epeak=np.append(short_merger_Epeak, getting_the_GBM_E_peak(sample_2_data, \
                                    sample_3_data.at[i, 'Spectral Model'], \
                                               sample_3_data.at[i, 'Fermi row']))
    short_merger_flu=np.append(short_merger_flu, float(sample_3_data.at[i, 'GBM fluence']))
    short_merger_Epeak_over_s_err=np.append(short_merger_Epeak_over_s_err, \
                       (getting_the_GBM_E_peak_err(sample_2_data, 
                            sample_2_data.at[q, 'flnc_best_fitting_model'], q))/\
                                       float(sample_2_data.at[i, 'fluence   ']))
    short_merger_Epeak_err=np.append(short_merger_Epeak, \
                                       getting_the_GBM_E_peak_err(sample_2_data, \
                                    sample_3_data.at[i, 'Spectral Model'], \
                                               sample_3_data.at[i, 'Fermi row']))
    short_merger_flu_err=np.append(short_merger_flu_err, \
                                     float(sample_3_data.at[i, 'Fluence error']))
    short_merger_ag_flu=np.append(short_merger_ag_flu, \
            afterglow_fluence_finder(spectral_sm_sample_3_data, more_xrt_data, i, \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Time"], \
            #heck if i remember why it's i-1 and I just wrote this 10 seconds ago. 
            #I discovered it by trial and error
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Error"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Error"]))
    short_merger_ag_flu_err=np.append(short_merger_ag_flu_err, \
             afterglow_fluence_err_finder(spectral_sm_sample_3_data, more_xrt_data, i, \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Error"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Error"]))
    short_merger_t90_err=np.append(short_merger_t90_err, float(sample_2_data.at[q, \
                                                                        't90_error']))
    Refined_Goldstein_Data[p,0]=sample_2_data.at[q, 'name        '][3:]
    Refined_Goldstein_Data[p,1]=short_merger_Epeak_over_s[-1]
    Refined_Goldstein_Data[p,2]=short_merger_Epeak_over_s_err[-1]
    Refined_Goldstein_Data[p,3]=short_merger_HR[-1]
    Refined_Goldstein_Data[p,4]=short_merger_t90[-1]
    Refined_Goldstein_Data[p,5]=short_merger_t90_err[-1]
    Refined_Goldstein_Data[p,6]=short_merger_Epeak[-1]
    Refined_Goldstein_Data[p,7]=short_merger_Epeak_err[-1]
    Refined_Goldstein_Data[p,8]=short_merger_flu[-1]
    Refined_Goldstein_Data[p,9]=short_merger_flu_err[-1]
    Refined_Goldstein_Data[p,10]=short_merger_ag_flu[-1]
    Refined_Goldstein_Data[p,11]=short_merger_ag_flu_err[-1]
    p=p+1
# print(p)
short_merger_t90_err=np.where(spectral_sm_sample_3_data['t90_error'] != 'N/A', \
        spectral_sm_sample_3_data['t90_error'], 0)
if len(short_merger_t90_err)>0:
    for j in range(0, len(short_merger_t90_err)):
        short_merger_t90_err[j]=float(short_merger_t90_err[j])
        
spectral_unknown_progenitor_sample_3_data=spectral_sample_3_data[~(lm_mask|lc_mask|sm_mask|sc_mask)]
for i in spectral_unknown_progenitor_sample_3_data.index.to_list():
    q=spectral_sample_3_data.at[i, 'Fermi row']
    unknown_progenitor_HR=np.append(unknown_progenitor_HR, Fermi_HR_func(sample_2_data, \
                                sample_2_data.at[q, 'flnc_best_fitting_model'], q))
    unknown_progenitor_t90=np.append(unknown_progenitor_t90, \
                            float(sample_2_data.at[q, 't90     ']))
    unknown_progenitor_Epeak_over_s=np.append(unknown_progenitor_Epeak_over_s, \
                       (getting_the_GBM_E_peak(sample_2_data, sample_2_data.at[q, 'flnc_best_fitting_model'], i))/\
                                       float(sample_2_data.at[q, 'fluence   ']))
    unknown_progenitor_Epeak=np.append(unknown_progenitor_Epeak, getting_the_GBM_E_peak(sample_2_data, \
                                    sample_3_data.at[i, 'Spectral Model'], \
                                               sample_3_data.at[i, 'Fermi row']))
    unknown_progenitor_flu=np.append(unknown_progenitor_flu, float(sample_3_data.at[i, 'GBM fluence']))
    unknown_progenitor_Epeak_over_s_err=np.append(unknown_progenitor_Epeak_over_s_err, \
                       (getting_the_GBM_E_peak(sample_2_data, sample_2_data.at[q, 'flnc_best_fitting_model'], i))/\
                                       float(sample_2_data.at[q, 'fluence   ']))
    unknown_progenitor_Epeak_err=np.append(unknown_progenitor_Epeak, \
                                       getting_the_GBM_E_peak_err(sample_2_data, \
                                    sample_3_data.at[i, 'Spectral Model'], \
                                               sample_3_data.at[i, 'Fermi row']))
    unknown_progenitor_flu_err=np.append(unknown_progenitor_flu_err, \
                                     float(sample_3_data.at[i, 'Fluence error']))
    unknown_progenitor_ag_flu=np.append(unknown_progenitor_ag_flu, \
            afterglow_fluence_finder(spectral_unknown_progenitor_sample_3_data, \
            more_xrt_data, i, \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Time"], \
            #heck if i remember why it's i-1 and I just wrote this 10 seconds ago. 
            #I discovered it by trial and error
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Error"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Error"]))
    unknown_progenitor_ag_flu_err=np.append(unknown_progenitor_ag_flu_err, \
             afterglow_fluence_err_finder(spectral_unknown_progenitor_sample_3_data,\
            more_xrt_data, i, \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Time"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Start_Error"], \
            Afterglow_of_Sample_3_Data.at[i-1, "AG_Stop_Error"]))
    unknown_progenitor_t90_err=np.append(unknown_progenitor_t90_err, \
                            float(sample_2_data.at[q, 't90_error']))
    Refined_Goldstein_Data[p,0]=sample_2_data.at[q, 'name        '][3:]
    Refined_Goldstein_Data[p,1]=unknown_progenitor_Epeak_over_s[-1]
    Refined_Goldstein_Data[p,2]=unknown_progenitor_Epeak_over_s_err[-1]
    Refined_Goldstein_Data[p,3]=unknown_progenitor_HR[-1]
    Refined_Goldstein_Data[p,4]=unknown_progenitor_t90[-1]
    Refined_Goldstein_Data[p,5]=unknown_progenitor_t90_err[-1]
    Refined_Goldstein_Data[p,6]=unknown_progenitor_Epeak[-1]
    Refined_Goldstein_Data[p,7]=unknown_progenitor_Epeak_err[-1]
    Refined_Goldstein_Data[p,8]=unknown_progenitor_flu[-1]
    Refined_Goldstein_Data[p,9]=unknown_progenitor_flu_err[-1]
    Refined_Goldstein_Data[p,10]=unknown_progenitor_ag_flu[-1]
    Refined_Goldstein_Data[p,11]=unknown_progenitor_ag_flu_err[-1]
    p=p+1
# print(p)
unknown_progenitor_t90_err=np.where(spectral_unknown_progenitor_sample_3_data['t90_error'] != 'N/A', \
        spectral_unknown_progenitor_sample_3_data['t90_error'], 0)
if len(unknown_progenitor_t90_err)>0:
    for j in range(0, len(unknown_progenitor_t90_err)):
        unknown_progenitor_t90_err[j]=float(unknown_progenitor_t90_err[j])
# print("There are {} known long collapsars, {} known long \
# mergers, {} known short mergers, and {} known short collapsars in Sample 3.".\
# format(len(long_collapsar_t90), len(long_merger_t90), \
#        len(short_merger_t90), len(short_collapsar_t90)))
# print('There are {}/{} bursts in total'.format(np.count_nonzero(\
#                                                 ~np.isnan(long_collapsar_HR))\
#         +np.count_nonzero(~np.isnan(long_merger_HR))+\
#         np.count_nonzero(~np.isnan(short_merger_HR))+\
#         np.count_nonzero(~np.isnan(short_collapsar_HR))+\
#         np.count_nonzero(~np.isnan(unknown_progenitor_HR)), \
#                                               len(spectral_sample_3_data)))

Goldstein_Sample_3_Data=pd.DataFrame(Refined_Goldstein_Data, columns=["Name",\
                                        "Peak E. over Flue.", "E_P_Over_S_Err",\
                                        "Hardness Ratio", "t90", "t90_err", "Peak_E",\
                                        "Peak_E_Err", "Fluence", "Fluence_Err", \
                                        "AG_fluence", "AG_fluence_err"])
Goldstein_Sample_3_Data.to_excel("Goldstein_Sample_3_Data.xlsx")

# %%
# sample_3_data.columns


