#%% Import needed assets
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib_venn

#%% Import data
df = pd.read_excel("../data/Data Day 1.xlsx")
dfblind = pd.read_excel("../data/Data_rerun_day1_OnlyPep_Nodups.xlsx")
#%% Trim excess data
data = df.iloc[:,1:51]
data = data.drop(columns = ["Start", "End"])
datablind = dfblind.drop(columns = ["Start", "End"])

#%% Define used functions
def lookForExistingPeptide(rowN, columnRange,data):
    for colN in columnRange:
        if data.iloc[rowN][colN] > 0:
            return True
    return False

def makeBooleanDataframe(booleanDf, data,saRange,paRange,ctrlRange):
    for rowN in range(data.shape[0]):
        newRow = []
        newRow.append(data.iloc[rowN][0])
        inSa = lookForExistingPeptide(rowN, saRange, data)
        inPa = lookForExistingPeptide(rowN, paRange, data)
        inCtrl = lookForExistingPeptide(rowN, ctrlRange, data)
        for x in [inSa,inPa,inCtrl]: newRow.append(x)
        booleanDf.loc[len(booleanDf)] = newRow
        
def makeSets(booleanDf):
    setSa = set()
    setPa = set()
    setCtrl = set()
    setList = [setSa, setPa, setCtrl]
    for rowN in range(booleanDf.shape[0]):
        for colN in [1,2,3]:
            if booleanDf.iloc[rowN][colN] == True:
                setList[colN-1].add(booleanDf.iloc[rowN][0])
    return setList
#%% Convert to boolean dataframes

#Create empty dataframes
booleanDf = pd.DataFrame(columns = ["Peptide","S.a","P.a","Ctrl"])
booleanDfBlind = pd.DataFrame(columns = ["Peptide","S.a-blind","P.a-blind","Ctrl-blind"])

#Define range of groups in dataframes
saRange = range(1,18)
paRange = range(18,35)
ctrlRange = range(35,48)
saBlindRange = [9,10,11,12]
paBlindRange = [1,2,4,8]
ctrlBlindRange = [3,5,6,7]


#Go through datasets and convert to boolean
makeBooleanDataframe(booleanDf,data,saRange,paRange,ctrlRange)
makeBooleanDataframe(booleanDfBlind,datablind,saBlindRange,paBlindRange,ctrlBlindRange)


#%% Make sets from boolean dataframes
[setSa, setPa, setCtrl] = makeSets(booleanDf)
[setSaBlind, setPaBlind, setCtrlBlind] = makeSets(booleanDfBlind)


#%% Fix font sizes and styles
plt.rcParams['pdf.fonttype']=42
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