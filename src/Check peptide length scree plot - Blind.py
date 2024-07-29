# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 15:19:29 2024

@author: fr1682fo
"""
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

#Check peptide length scree plot

#%%
df = pd.read_excel("Data Day 1 and 2 and 3 OnlyPep NoDups w Blind.xlsx")
data = df.iloc[3:,:]
days = df.iloc[0:1,:]
groups = df.iloc[1:2,:]
blinds = df.iloc[2:3,:]
#%%
for x in range(data.shape[0]):
    for y in range(data.shape[1]-1):
        if data.iloc[x][y+1] != 0:
            data.iloc[[x],[y+1]] = 1

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


#%% Color and position dictionary
colorDict = {}
colorDict["S.a"] = "gold"
colorDict["P.a"] = "aqua"
colorDict["Ctrl"] = "black"

groupPosDict = {}
groupPosDict["No"] = -1/2
groupPosDict["Yes"] = 1/2
#%% Plot day by day

#Count if present in any sample


width = 0.4

x = list(range(1,longest+1))
for group in ["S.a","P.a","Ctrl"]:
    fig, ax = plt.subplots()
    for blind in ["No","Yes"]:
        y = [0]*len(x)
        for i in range(data.shape[0]):
            length = findLength(data.iloc[i,0])
            n = 0
            for j in range(1,data.shape[1]):
                if (group == groups.iloc[0,j]) and (blinds.iloc[0,j] == blind) and 1 == days.iloc[0,j]:
                    n += data.iloc[i,j]
                    if n != 0:
                        y[length-1] += 1
                        break
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
        
        if blind == "Yes":
            col = "red"
        else:
            col = colorDict[group]
        ax.bar([x_pos+width*groupPosDict[blind] for x_pos in x],maxstand_y,width,color = col)
        plt.title(label = group)
        plt.xlabel('Peptide length (AA)')
        plt.ylabel('Relative abundance')
        plt.show()
        
        
        
    
             
#%% Count # of occurances
x = list(range(1,longest+1))
for day in ["1","2","3"]:
    fig, ax = plt.subplots()
    for group in ["S.a","P.a","Ctrl","Double","Acc Double"]:
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
        ax.bar([x_pos+width*groupPosDict[group] for x_pos in x],sumstand_y,width,color = colorDict[group])
        plt.title(label = "Day "+day)
        plt.xlabel('Peptide length (AA)')
        plt.ylabel('Relative abundance')
        plt.show()



