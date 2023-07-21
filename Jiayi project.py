# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


# Common imports
import sklearn
assert sklearn.__version__ >= "0.20"

import numpy as np
import os
import csv

# Using pandas to format datasets
import pandas as pd

# To make this notebook's output stable across runs
np.random.seed(42)

# Pretty figures
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
plt.style.use('seaborn-poster')
mpl.rc('axes', labelsize=14)
mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=12)

import math # For log and mean operations

#Paths
path1 = r"C:\Users\Maya\OneDrive - University of Calgary\Desktop\3featureinvestigation\GSE158854_Avon_Breast_Cancer_RSEM_FPKM.txt.gz" #PPBC study
#path2 = r'C:\Users\Maya\OneDrive - University of Calgary\Desktop\3featureinvestigation\GSM6181401_sample_5.txt.gz' #FSS study
path3 = r'C:\Users\Maya\OneDrive - University of Calgary\Desktop\3featureinvestigation\GSM7468914_SO_5602_a2_C1_Att_24hr_Readcount_by_Gene_Name.txt' #Cellular stress



### PPBC dataset
#Opening and cleaning up PPBC dataset by dropping zero values etc.
df_PPBC = pd.read_csv(path1, compression='gzip', header=0, sep='\t', quotechar='"')
df_PPBC = df_PPBC.set_index('gene_symbol')
df_PPBC = df_PPBC.loc[(df_PPBC!=0).all(axis=1)] # Drop all rows containing any zeroes
df_PPBC.head()
#Taking the difference between the pathological condition and control condition
df_PPBC['PPBC'] = (df_PPBC.filter(like='PPBC', axis=1)).mean(axis=1) # Taking the log and mean of all cancerous samples
df_PPBC['Nullip'] = (df_PPBC.filter(like='Nullip', axis=1)).mean(axis=1) # Taking the log and mean of all non-cancerous samples
df_PPBC['DEG_PPBC'] = df_PPBC['PPBC'] - df_PPBC['Nullip'] # Taking the difference between the cancerous and non-cancerous means to create the differential gene expression for each gene in the dataset
df_PPBC[['PPBC', 'Nullip', 'DEG_PPBC']].head()


### FSS dataset
#df_2 = pd.read_csv(path2, header=0, sep='\t', quotechar='"')
#df_2.head()


### Cell stress dataset
df_CS = pd.read_csv(path3, header=0, sep='\t', quotechar='"')
df_CS.drop(list(df_CS.filter(regex='Chr|Start|End|Strand|Length')), axis=1, inplace=True)
df_CS.rename(columns={'Reads': 'DEG_CS', 'Geneid': 'gene_symbol'}, inplace=True)
df_CS = df_CS.set_index('gene_symbol')
df_CS.head()

df_X = pd.merge(df_PPBC['DEG_PPBC'], df_CS['DEG_CS'], on='gene_symbol') # Merging dataframes
df_X = df_X.dropna(axis=0) # Drop all NaN values
df_X.drop(df_X[df_X['DEG_PPBC'] >= 50000].index, inplace = True)
df_X.drop(df_X[df_X['DEG_CS'] >= 50000].index, inplace = True)
df_X.head()


fig = plt.figure(figsize = (250, 250))
ax = plt.axes()
 
# Creating plot
x=df_X['DEG_PPBC']
y=df_X['DEG_CS']
my_cmap = plt.get_cmap('hsv')
ax.scatter(x, y, alpha = 0.8, c=(x+y), cmap=my_cmap)
plt.title("Postpartum breast cancer by gene expression")
ax.set_xlabel('PPBC', fontweight ='bold', labelpad=40)
ax.set_ylabel('Cell stress', fontweight ='bold', labelpad=40)
 
# show plot
plt.show()



###########################









# fig = plt.figure(figsize = (250, 250))
# ax = plt.axes()
 
# # Creating plot
# x=df_PPBC['DEG_PPBC']
# y=df_CS['DEG_CS']
# my_cmap = plt.get_cmap('hsv')
# ax.scatter3D(x, y, alpha = 0.8, c=(x+y), cmap=my_cmap)
# plt.title("Postpartum breast cancer by gene expression")
# ax.set_ylim([0, 5000])
# ax.set_xlim([0, 5000])
# ax.set_xlabel('S10_PPBC_FPKM', fontweight ='bold', labelpad=40)
# ax.set_ylabel('S16_PPBC_FPKM', fontweight ='bold', labelpad=40)
 
# # show plot
# plt.show()

# from matplotlib import cm
# from matplotlib.ticker import LinearLocator, FormatStrFormatter
# from mpl_toolkits.mplot3d import Axes3D


# fig = plt.figure(figsize = (250, 250))
# ax = plt.axes(projection ="3d")
 
# # Creating plot
# x=df_1['S10_PPBC_FPKM']
# y=df_1['S16_PPBC_FPKM']
# z=df_1['S2_PPBC_FPKM']
# my_cmap = plt.get_cmap('hsv')
# ax.scatter3D(x, y, z, alpha = 0.8, c=(x+y+z), cmap=my_cmap)
# plt.title("Postpartum breast cancer by gene expression")
# ax.set_ylim([0, 5000])
# ax.set_xlim([0, 5000])
# ax.set_zlim([0, 5000])
# ax.set_xlabel('S10_PPBC_FPKM', fontweight ='bold', labelpad=40)
# ax.set_ylabel('S16_PPBC_FPKM', fontweight ='bold', labelpad=40)
# ax.set_zlabel('S2_PPBC_FPKM', fontweight ='bold', labelpad=40)
 
# # show plot
# plt.show()


