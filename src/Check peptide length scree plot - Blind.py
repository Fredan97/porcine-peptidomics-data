# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 15:19:29 2024

@author: fr1682fo
"""
import pandas as pd
import matplotlib.pyplot as plt
import math

#Check peptide length scree plot

#%% Import and trim data
df = pd.read_excel("../data/Data Day 1 and 2 and 3 OnlyPep NoDups w Blind.xlsx")
data = df.iloc[3:,:]
days = df.iloc[0:1,:]
groups = df.iloc[1:2,:]
blinds = df.iloc[2:3,:]
#%% Convert to log
for x in range(data.shape[0]):
    for y in range(data.shape[1]-1):
        if data.iloc[x][y+1] != 0:
            data.iloc[x,y+1] =  math.log(data.iloc[x][y+1],2)

#%%
def findLength(seq):
    length = 0
    for char in seq:
        if char not in ["(","+","1","5",".","9",")"]:
            length += 1
    return length

#%%
longest = 0
for sequence in df.iloc[2:,0]:
    if findLength(sequence) > longest:
        longest = findLength(sequence)


#%% Create color and linestyle dictionary
colorDict = {}
colorDict["S.a"] = "gold"
colorDict["P.a"] = "aqua"
colorDict["Ctrl"] = "black"
colorDict["Double"] = "lime"
colorDict["Acc Double"] = "red"

styleDict = {}
styleDict["No"] = "-"
styleDict["Yes"] = "--"


#%% Fix font sizes and styles
plt.rcParams['font.size']=7
plt.rcParams["font.family"] = "Arial"

#%% Create labels and handles for legend
grouphandles = [plt.Line2D([0,0],[0,0],color=colorDict[i], linestyle='-') for i in ["Ctrl","S.a","P.a"]]
grouplabels = ["Ctrl","S.a","P.a"]
blindhandles = [plt.Line2D([0,0],[0,0],color = "gray",linestyle = styleDict[i]) for i in ["No","Yes"]]
blindlabels = ["Known","Blind"]

#%% Sum and plot log2 of each length
fig, ax = plt.subplots()
x = list(range(1,longest+1))
for day in ["1"]:
    for blind in ["No","Yes"]:
        for group in ["S.a","P.a","Ctrl"]:
            y = [0]*len(x)
            for i in range(data.shape[0]):
                length = findLength(data.iloc[i,0])
                for j in range(1,data.shape[1]):
                    if (group == groups.iloc[0,j]) and (blinds.iloc[0,j]) == blind:
                        y[length-1] += data.iloc[i,j]
            sum_y = sum(y)
            if (sum_y != 0):
                sumstand_y = [float(k)/sum_y for k in y]
            else:
                sumstand_y = y
            ax.plot(x,sumstand_y,0.1,color = colorDict[group],linestyle = styleDict[blind])
plt.xlabel('Peptide length (AA)')
plt.ylabel('Relative intensity')  
plt.legend(grouphandles+blindhandles,grouplabels+blindlabels,loc = 'upper right')
plt.show()      


