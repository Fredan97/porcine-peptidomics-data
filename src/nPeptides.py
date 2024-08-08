#%% Import needed assets
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import random as rd

#%% Import data
df = pd.read_excel("../data/Data Day 1 and 2 and 3 OnlyPep NoDups.xlsx")
#%% Trim data
data = df.iloc[2:,1:]
group = df.iloc[1,1:]
day = df.iloc[0,1:]


#%% Define functions
def CountPeptidesPresent(column):
    pepN = 0
    for colN in range(column.shape[0]):
        if (column.iloc[colN] > 0):
            pepN += 1
    return pepN

#%% Find number of peptides for each sample
pepN = []
for colN in range(data.shape[1]):
    pepN.append(CountPeptidesPresent(data.iloc[:,colN]))

#%%Prepare dictionaries for plotting properties
posMap = group.map({"Ctrl": 0, "S.a":1,"P.a":2,"Double":3,"Acc Double":4}).tolist()

colorMap = group.map({"S.a": 0, "P.a": 1, "Ctrl": 2, "Double": 3, "Acc Double": 5}).tolist()

colors = []
for n in colorMap:
    colors.append(sns.color_palette()[n])


markers = day.map({1: "o", 2: "s", 3: "X"}).tolist()

#%% Fix font sizes and styles
plt.rcParams['pdf.fonttype']=42
plt.rcParams["font.family"] = "Arial"

#%% Plot the results
plt.figure(figsize = (5,4.5))
for i in range(len(pepN)):
    plt.scatter(
        posMap[i]*0.6+day[i]*0.2+rd.uniform(0, 0.1),
        pepN[i],
        color = colors[i],
        marker = markers[i],
        edgecolor='w',
        linewidths = 0.6)


unique_days = [1, 2, 3]
day_handles = [plt.Line2D([0,0],[0,0],color='gray', marker={1: "o", 2: "s", 3: "X"}[day], linestyle='') for day in unique_days]
day_labels = [f'Day {day}' for day in unique_days]
group_labels = ["Ctrl","S.a","P.a","Double infection","Accidental double infection"]
group_handles = [plt.Line2D([0,0],[0,0],color=sns.color_palette()[i], marker='o', linestyle='') for i in [2,0,1,3,5]]
plt.legend(day_handles + group_handles, day_labels + group_labels, loc='upper right')
plt.tick_params(labelbottom = False, bottom = False)
plt.ylabel("Observed peptides")
plt.xlim([0,4])

plt.show()