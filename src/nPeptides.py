#%% Import needed assets
import pandas as pd
from matplotlib import pyplot as plt

#%% Import data
df = pd.read_excel("../data/Data Day 1 and 2 and 3 OnlyPep NoDups w Blind.xlsx")
#%% Trim data
data = df.iloc[:,1:]
data = df.iloc[3:,1:]
group = df.iloc[1,1:]
day = df.iloc[0,1:]
blinds = df.iloc[2,1:]


#%% Define functions
def CountPeptidesPresent(column):
    pepN = 0
    for colN in range(column.shape[0]):
        if column[colN] > 0:
            pepN += 1
    return pepN

#%% Find number of peptides for each sample
