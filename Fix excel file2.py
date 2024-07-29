import pandas as pd


#%%
data3 = pd.read_excel("Data Day 3.xlsx")



#%%
finalData = pd.read_excel("Data Day 1 and 2 and 3.xlsx")


#%%
k = 0
for x in data3.get("Peptide"):
    emptyRow = [0]*109
    newRow = []
    l = 0
    exists = False
    for y in finalData.get("Peptide"):
        if exists == False:
            if x == y and data3.get("Protein")[k] == finalData.get("Protein")[l]:
                exists = True
                for i in range(4):
                    finalData.iloc[[l],[113+i]] = data3.iloc[k][4+i]
        l += 1
    if exists == False:
        for i in range(4):
            emptyRow.append(data3.iloc[k][4+i])
        for i in range(4):
            newRow.append(data3.iloc[k][i])
        for i in emptyRow:
            newRow.append(i)
        finalData.loc[len(finalData)] = newRow
    k += 1
            



    
#finalData.to_excel("Data Day 1 and 2 and 3.xlsx",index=False)