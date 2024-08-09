# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 15:19:29 2024

@author: fr1682fo
"""
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import math

#Check peptide length scree plot

#%% Read data
df = pd.read_excel("../data/Data Day 1 and 2 and 3 OnlyPep NoDups.xlsx")
data = df.iloc[2:,:]
days = df.iloc[0:1,:]
groups = df.iloc[1:2,:]
#%% Convert to log
for x in range(data.shape[0]):
    for y in range(data.shape[1]-1):
        if data.iloc[x][y+1] != 0:
            data.iloc[[x],[y+1]] = math.log(data.iloc[x][y+1],2)

#%% Define needed functions
def findLength(seq):
    length = 0
    for char in seq:
        if char not in ["(","+","1","5",".","9",")"]:
            length += 1
    return length

#%% Find longest peptide in all samples
longest = 0
for sequence in df.iloc[2:,0]:
    if findLength(sequence) > longest:
        longest = findLength(sequence)


#%% Color and position dictionary
colorDict = {}
colorDict["S.a"] = "gold"
colorDict["P.a"] = "aqua"
colorDict["Ctrl"] = "black"
colorDict["Double"] = "lime"
colorDict["Acc Double"] = "red"


#%% Fix font sizes and styles
plt.rcParams['font.size']=7
plt.rcParams["font.family"] = "Arial"

#%% Create labels and handles for legend
linehandles = [plt.Line2D([0,0],[0,0],color=colorDict[i], linestyle='-') for i in ["S.a","P.a","Ctrl","Double","Acc Double"]]
linelabels = ["S.a","P.a","Ctrl","Double infection","Accidental double infection"]

    
             
#%% Sum log2 of each length
x = list(range(1,longest+1))
for day in ["1","2","3"]:
    fig, ax = plt.subplots()
    for group in ["S.a","P.a","Ctrl","Double","Acc Double"]:
        if not(day == "3" and group != "Double"):
            y = [0]*len(x)
            for i in range(data.shape[0]):
                length = findLength(data.iloc[i,0])
                for j in range(1,data.shape[1]):
                    if (group == groups.iloc[0,j]) and (str(days.iloc[0,j]) == day):
                        y[length-1] += data.iloc[i,j]
            sum_y = sum(y)
            max_y = max(y)
            if (max_y != 0):
                maxstand_y = [float(k)/max_y for k in y]
            else:
                maxstand_y = y
            if (sum_y != 0):
                sumstand_y = [float(k)/sum_y for k in y]
            else:
                sumstand_y = y
            ax.plot(x,sumstand_y,0.1,color = colorDict[group])
            plt.title(label = "Day "+day)
            plt.xlabel('Peptide length (AA)')
            plt.ylabel('Relative intensity')
    if (day == "3"):
        plt.legend(handles = [linehandles[3]],labels = [linelabels[3]],loc="upper right")
    else:
        plt.legend(linehandles,linelabels,loc = 'upper right')
    plt.show()
        



