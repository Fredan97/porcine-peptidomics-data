import pandas as pd


#%%
data1 = pd.read_excel("Data Day 1.xlsx")
data2 = pd.read_excel("Data Day 2.xlsx")


#%%
finalData = pd.read_excel("Data Day 1 and 2.xlsx")


#%%
k = 0
for x in data1.get("Peptide"):
    isInDay2 = False
    l = 0
    newRow = []
    for y in data2.get("Peptide"):
        if x == y and data1.get("Protein")[k] == data2.get("Protein")[l]:
            isInDay2 = True
            
            for n in range(data1.shape[1]):
                newRow.append(data1.iloc[k][n])
            for n in range(4,data2.shape[1]):
                newRow.append(data2.iloc[l][n])
            finalData.loc[len(finalData)] = newRow
            k += 1
            break
        else:
            l += 1
    if not isInDay2:
        for n in range(data1.shape[1]):
            newRow.append(data1.iloc[k][n])
        for n in range(4,data2.shape[1]):
            newRow.append(0)
        finalData.loc[len(finalData)] = newRow
        k += 1
        


j = 0
for x in data2.get("Peptide"):
    isInFinalData = False
    k = 0
    newRow = []
    for y in finalData.get("Peptide"):
        if x == y and data2.get("Protein")[j] == finalData.get("Protein")[k]:
            isInFinalData = True
            break
        else:
            k += 1
    if not isInFinalData:
        for n in range(4):
            newRow.append(data2.iloc[j][n])
        for n in range(4,data1.shape[1]):
            newRow.append(0)
        for n in range(4,data2.shape[1]):
            newRow.append(data2.iloc[j][n])
        finalData.loc[len(finalData)] = newRow
    j += 1
    
    
#finalData.to_excel("Data Day 1 and 2.xlsx",index=False)