#!/usr/bin/env python
# coding: utf-8

# All this should do is give you the predicted classification of a burst given its t90, peak energy, and fluence. It might even work too.

# In[ ]:


##Imports##
import math
import scipy
import gc
import numpy as np
#from svm import *
import scipy.stats as st
import re
from statistics import NormalDist
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pdb
import pandas as pd
from sklearn.datasets import load_svmlight_file
import matplotlib
import matplotlib.colors as mcolors
import sklearn
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
from matplotlib.collections import LineCollection
import matplotlib.cm as cm
import argparse
import pdb
from matplotlib.backends.backend_pdf import PdfPages
import os

cwd = os.getcwd()

#Known progenitors
Collapsars = [050416461.0, 50525002.0, 050824966.0, 060218148.0, 060729800.0, 060904104.0, 070419447.0, 071025172.0, 071112772.0, 080319258.0, 081007224.0, 090618353.0, 091127976.0,                100316531.0, 100418882.0, 101219686.0, 101225776.0, 111209300.0, 111211928.0, 111228656.0, 120422300.0, 120714888.0, 120729455.0, 130215063.0, 130427324.0, 130702003.0,                     130831544.0, 140206303.0, 140606133.0, 150818483.0, 161219783.0, 161228552.0, 171010792.0, 171205306.0, 180720598.0, 180728728.0, 190114872.0, 190829830.0, 200826187.0,                          210210083.0, 211023545.0, 221009553.0, 200826187.0, 150210935.0]
Mergers=[070809807.0, 130603659.0, 150101641.0, 160624477.0, 160821936.0, 170817528.0, 200522487.0, 060614530.0, 111005336.0, 120304248.0, 211211549.0, 230307655.0, 050509166.0,           050709942.0, 051210240.0, 070714207.0, 071227842.0, 080503518.0, 080905499.0, 090515198.0, 160303454.0]

parser = argparse.ArgumentParser()
parser.add_argument("-Ep", "--peakEnergy_keV", required=True)
parser.add_argument("-S", "--Fluence_erg", required=True)
parser.add_argument("-t90", "--Duration_s", required=True)
args = parser.parse_args()

Goldstein_data=pd.read_csv("Goldstein_Full_DataSet_w_Name.csv")
Goldstein_training=pd.read_csv("Goldstein_Training_Set.csv")
x1=Goldstein_training['Peak_E_over_Flu'].values
x2=Goldstein_training['t90'].values
X=np.vstack((x1,x2)).T
Y=Goldstein_training['Classification.'].values
Y[Y==2]=0
 
xx1=Goldstein_data['Peak_E_over_Flu'].values
xx2=Goldstein_data['t90'].values
Xtest=np.vstack((xx1,xx2)).T
 
model=sklearn.svm.SVC(kernel='rbf',probability=True,class_weight='balanced')
model.fit(X,Y)
p=model.predict_proba(Xtest)

X_single=np.log10(float(args.peakEnergy_keV)/float(args.Fluence_erg))
Y_single=np.log10(float(args.Duration_s))
Xtest_single=np.vstack((X_single, Y_single)).T
preds_single=model.predict_proba(Xtest_single)

All_Final_Data=pd.DataFrame(np.vstack((Goldstein_data['Name'].values, Goldstein_data['Peak_E_over_Flu'].values, Goldstein_data['t90'].values, p[:,0], p[:,1])).T,                             columns=['Name', 'Indicator', 't90', 'Pm', 'Pc'])
# print(All_Final_Data.sort_values(by=['Indicator'], ascending=False).to_string())
array_170817A=All_Final_Data[All_Final_Data['Name']==170817529.0]

#
X_SBOAT=[(np.log10(1000/3.1475e-03))]
Y_SBOAT=[np.log10(34.561)]
XSBOAT=np.vstack((X_SBOAT, Y_SBOAT)).T
preds_SBOAT=model.predict_proba(XSBOAT)
#
X_BOAT=[np.log10(1387/9.47e-02)]
Y_BOAT=[np.log10(289)]
XBOAT=np.vstack((X_SBOAT, Y_SBOAT)).T
preds_BOAT=model.predict_proba(XBOAT)
#
 
dfc=Goldstein_data[Goldstein_data['Classification.']=='collapsar']
dfm=Goldstein_data[Goldstein_data['Classification.']=='merger']
looong=Goldstein_data[Goldstein_data['t90']>np.log10(4.2)]
shrt=Goldstein_data[Goldstein_data['t90']<= np.log10(4.2)]
n_colors = 2
cmap = matplotlib.cm.get_cmap('cividis')

#
print("We predict a {:0.2f} chance that your burst was a merger and {:0.2f} that it was a collapsar.".format(preds_single[0][0], preds_single[0][1]))
if float(args.Fluence_erg)>1E-3:
    print("Your burst was unusually bright, and may falsely be classified as a collapsar as a result.")
##
with PdfPages('your_GRB_in_context.pdf', cwd) as export_pdf:
    pp=plt.scatter(10**xx2, 10**xx1,c=p[:,1],s=8,cmap='cividis', vmin=0, vmax=1)
    # plt.errorbar(y=10**Goldstein_data['Peak_E_over_Flu'], x=10**Goldstein_data['t90'], yerr=Goldstein_data['E_P_Over_S_Err'], xerr=Goldstein_data['t90_err'], \
                #  ls='none', ecolor=cm.cividis(p[:,1]))
    collapsars=plt.scatter(10**dfc['t90'],10**dfc['Peak_E_over_Flu'],c='black',label='known collapsar',s=25, marker="*")
    # plt.errorbar(x=10**dfc['t90'], y=10**dfc['Peak_E_over_Flu'], yerr=dfc['E_P_Over_S_Err'], xerr=dfc['t90_err'], \
                #  ls='none', c='black')
    mergers=plt.scatter(10**dfm['t90'],10**dfm['Peak_E_over_Flu'],c='orangered',label='known merger',s=25, marker="X")
    # plt.errorbar(x=10**dfm['t90'], y=10**dfm['Peak_E_over_Flu'], yerr=dfm['E_P_Over_S_Err'], xerr=dfm['t90_err'], \
                #  ls='none', c='orangered')
    first_GW_GRB=plt.scatter(10**array_170817A['t90'], 10**array_170817A['Indicator'], s=25, c='darkseagreen', label='GRB 170817A', marker="8")
    SBOAT_GRB=plt.scatter(10**Y_SBOAT[0], 10**X_SBOAT[0], s=25, c='peru', label='GRB 230307A', marker="s")
    BOAT_GRB=plt.scatter(10**Y_BOAT[0], 10**X_BOAT[0], s=25, c='olive', label='GRB 221009A', marker="P")
    your_burst=plt.scatter(10**Y_single, 10**X_single, s=25, c='darkorchid', label="Your Burst", marker="D")
    plt.legend(loc="lower left", prop={'size': 8}, handles=[collapsars, mergers, first_GW_GRB, SBOAT_GRB, BOAT_GRB, your_burst])
    plt.xlabel('Duration (s, GBM)')
    plt.ylabel('E$_p$/$S_{10 - 1000 keV}$ (keV cm$^{2}$ erg$^{-1}$, GBM)')
    plt.colorbar(pp,label='Collapsar Probability')
    plt.xscale('log')
    plt.yscale('log')
    export_pdf.savefig()
    plt.close()

gc.collect()
print('All Done.')