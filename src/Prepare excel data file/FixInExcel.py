# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 17:05:34 2023

@author: fr1682fo
"""


import pandas as pd
import os

def add_Variable(data_frame, final_data_frame, pos, sampleNum, columnName):
    #Calculate where the area column is located in data_frame
    colPos = 0
    for column in data_frame.columns:
        if column[0:4] != columnName:
            colPos += 1
        else:
            break
    
    #Add 1 extra to all areas
    for z in range(data_frame.shape[0]):
        data_frame.iloc[z,colPos] += 1
        
    
    

            
    k = 0
    for i in data_frame.get("Peptide"):
        if i[1] == "." and i[-2] == ".":
            x = i[2:-2]
        elif i[1] == ".":
            x = i[2:]
        elif i[-2] ==".":
            x = i[:-2]
        else:
            x = i
        isInList = False
        l = 0
        for y in final_data_frame.get("Peptide"):
            if x == y and data_frame.get("Protein Accession")[k] == final_data_frame.get("Protein")[l]:
                isInList = True
                break
            l = l+1
        if isInList == True:
            final_data_frame.loc[l, "Sample "+str(sampleNum)] = data_frame.loc[k][colPos]
        else:
            newRow = []
            newRow.append(data_frame.get("Protein Accession")[k])
            newRow.append(x)
            newRow.append(data_frame.loc[k,"Start"])
            newRow.append(data_frame.loc[k,"End"])
            for n in range(pos-1):
                newRow.append(0)
            newRow.append(data_frame.loc[k][colPos])
            for n in range(55-pos):
                newRow.append(0)
            final_data_frame.loc[len(final_data_frame)] = newRow
        k += 1
        


def findFolderName(folderNumber, directory):
    for names in os.listdir():
        start = names[:43]
        if start[-1] == "_":
            number = start[41]
        else:
            number = start[41:43]
        if folderNumber == number:
            return names
        



Data = pd.read_excel("Data.xlsx")
pos = 1
rang = list(range(46)) + list(range(53,57)) + list(range(65,70))
for x in rang:
    folderName = findFolderName(str(x+1),os.getcwd())
    location = (os.getcwd()+'\\'+folderName)
    fileName = location+"\\protein-peptides.csv"
    add_Variable(pd.read_csv(fileName),Data,pos,x+1, "Area")
    pos += 1
        
#%%
#Run the following row when you know the code above works as intended
#Data.to_excel("Data.xlsx", index = False)
          
    
