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

path1 = r"C:\Users\Maya\OneDrive - University of Calgary\Desktop\GSE158854_Avon_Breast_Cancer_RSEM_FPKM.txt.gz"
path2 = 'GSM6181401_sample_5.txt.gz'
path3_1 = 'GSM7468915_SO_5602_a2_C1_Sus_24hr_Readcount_by_Gene_Name.txt'
path3_2 = 'GSM7468915_SO_5602_a2_C1_Sus_24hr_Readcount_by_Gene_Name.txt'

df_1 = pd.read_csv(path1, compression='gzip', header=0, sep='\t', quotechar='"')
df_1.head()


from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

me = plt.figure().gca(projection='3d')
me.scatter(df_1['S2_PPBC_FPKM'], df_1['S10_PPBC_FPKM'], df_1['S16_PPBC_FPKM'])
me.set_xlabel('S10_PPBC_FPKM')
me.set_ylabel('S16_PPBC_FPKM')
me.set_zlabel('S2_PPBC_FPKM')
plt.xlim([0, 10000])
plt.ylim([0, 105000])
plt.zlim([0, 5000])

plt.show()


