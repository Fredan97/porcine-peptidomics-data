#%% Import needed assets
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib_venn

#%% Import data
df = pd.read_excel("../data/Data Day 1 and 2 and 3 OnlyPep NoDups w Blind.xlsx")
#%% Trim excess data
data = df.iloc[3:,:]
groups = df.iloc[1,:]
blinds = df.iloc[2,:]
#%% Define used functions
def findRangeOfGroup(groups,blinds,blind,group):
    rangeOfGroup = []
    for i in range(len(groups)):
        if groups[i] == group and blinds[i] == blind:
            rangeOfGroup.append(i)
    return rangeOfGroup

def lookForExistingPeptide(rowN,columnRange,data):
    for colN in columnRange:
        if data.iloc[rowN][colN] > 0:
            return True
    return False

def makeBooleanDataframe(booleanDf, data,saRange,paRange,ctrlRange,saBlindRange,paBlindRange,ctrlBlindRange):
    for rowN in range(data.shape[0]):
        newRow = []
        newRow.append(data.iloc[rowN][0])
        inSa = lookForExistingPeptide(rowN, saRange, data)
        inPa = lookForExistingPeptide(rowN, paRange, data)
        inCtrl = lookForExistingPeptide(rowN, ctrlRange, data)
        inSaBlind = lookForExistingPeptide(rowN, saBlindRange, data)
        inPaBlind = lookForExistingPeptide(rowN, paBlindRange, data)
        inCtrlBlind = lookForExistingPeptide(rowN, ctrlBlindRange, data)
        for x in [inSa,inPa,inCtrl,inSaBlind,inPaBlind,inCtrlBlind]: newRow.append(x)
        booleanDf.loc[len(booleanDf)] = newRow
        
def makeSets(booleanDf):
    setSa = set()
    setPa = set()
    setCtrl = set()
    setSaBlind = set()
    setPaBlind = set()
    setCtrlBlind = set()
    setList = [setSa, setPa, setCtrl, setSaBlind, setPaBlind, setCtrlBlind]
    for rowN in range(booleanDf.shape[0]):
        for colN in range(len(setList)):
            if booleanDf.iloc[rowN][colN+1] == True:
                setList[colN].add(booleanDf.iloc[rowN][0])
    return setList
#%% Convert to boolean dataframes

#Create empty dataframes
booleanDf = pd.DataFrame(columns = ["Peptide","S.a","P.a","Ctrl","S.a-blind","P.a-blind","Ctrl-blind"])

#Define range of groups in dataframes
saRange = findRangeOfGroup(groups,blinds,"No","S.a")
paRange = findRangeOfGroup(groups,blinds,"No","P.a")
ctrlRange = findRangeOfGroup(groups,blinds,"No","Ctrl")
saBlindRange = findRangeOfGroup(groups,blinds,"Yes","S.a")
paBlindRange = findRangeOfGroup(groups,blinds,"Yes","P.a")
ctrlBlindRange = findRangeOfGroup(groups,blinds,"Yes","Ctrl")


#Go through datasets and convert to boolean
makeBooleanDataframe(booleanDf,data,saRange,paRange,ctrlRange,saBlindRange,paBlindRange,ctrlBlindRange)


#%% Make sets from boolean dataframes
[setSa, setPa, setCtrl, setSaBlind, setPaBlind, setCtrlBlind] = makeSets(booleanDf)


#%% Fix font sizes and styles
plt.rcParams['font.size']=7
plt.rcParams["font.family"] = "Arial"



#%% Make venn diagrams
plt.figure(figsize=(3,3))
matplotlib_venn.venn3(subsets = (setSa, setPa, setCtrl),set_labels = ("S.a","P.a","Ctrl"))
plt.figure(figsize=(3,3))
matplotlib_venn.venn3((setSaBlind,setPaBlind,setCtrlBlind),set_labels = ("S.a-blind","P.a-blind","Ctrl-blind"))
plt.figure(figsize=(3,3))
matplotlib_venn.venn2((setSa,setSaBlind),set_labels = ("S.a","S.a-blind"))
plt.figure(figsize=(3,3))
matplotlib_venn.venn2((setPa,setPaBlind),set_labels = ("P.a","P.a-blind"))
plt.figure(figsize=(3,3))
matplotlib_venn.venn2((setCtrl,setCtrlBlind),set_labels = ("Ctrl","Ctrl-blind"))