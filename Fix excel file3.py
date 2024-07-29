import pandas as pd


#%%
dataBlind = pd.read_excel("Data_rerun_day1_OnlyPep_Nodups.xlsx")



#%%
finalData = pd.read_excel("Data Day 1 and 2 and 3 OnlyPep NoDups w Blind.xlsx")


#%%
k = 0
for x in dataBlind.get("Peptide"):
    emptyRow = [0]*113
    newRow = []
    l = 0
    exists = False
    for y in finalData.get("Variable"):
        if exists == False:
            if x == y:
                exists = True
                for i in range(12):
                    finalData.iloc[[l],[114+i]] = dataBlind.iloc[k][3+i]
        l += 1
    if exists == False:
        for i in range(12):
            emptyRow.append(dataBlind.iloc[k][3+i])
        for i in range(1):
            newRow.append(dataBlind.iloc[k][i])
        for i in emptyRow:
            newRow.append(i)
        finalData.loc[len(finalData)] = newRow
    k += 1
            



    
#finalData.to_excel("Data Day 1 and 2 and 3.xlsx",index=False)