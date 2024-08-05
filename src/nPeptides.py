#%% Import needed assets
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

#%% Import data
df = pd.read_excel("../data/Data Day 1 and 2 and 3 OnlyPep NoDups w Blind.xlsx")
#%% Trim data
data = df.iloc[:,1:]
data = df.iloc[3:,1:]
group = df.iloc[1,1:]
day = df.iloc[0,1:]
blinds = df.iloc[2,1:]


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
colorMap = group.map({"S.a": 0, "P.a": 1, "Ctrl": 2, "Double": 3, "Acc Double": 5}).tolist()
colors = []
for n in colorMap:
    colors.append(sns.color_palette()[n])

edgeColors = blinds.map({"No": "w", "Yes": "black"}).tolist()

markers = day.map({1: "o", 2: "s", 3: "X"}).tolist()

#%% Plot the results
plt.figure(figsize = (10,8))
for i in range(len(pepN)):
    plt.scatter(
        i,
        pepN[i],
        color = colors[i],
        edgecolors = edgeColors[i],
        marker = markers[i])


unique_days = [1, 2, 3]
day_handles = [plt.Line2D([0,0],[0,0],color='gray', marker={1: "o", 2: "s", 3: "X"}[day], linestyle='') for day in unique_days]
day_labels = [f'Day {day}' for day in unique_days]
group_labels = list({"S.a": 0, "P.a": 1, "Ctrl": 2, "Double": 3, "Acc Double": 5}.keys())
group_handles = [plt.Line2D([0,0],[0,0],color=sns.color_palette()[i], marker='o', linestyle='') for i in [0,1,2,3,5]]
blind_labels = list({"Blinded"})
blind_handles = [plt.Line2D([0,0],[0,0],color='gray',marker = 'o',markeredgecolor = 'black', linestyle ='')]
plt.legend(day_handles + group_handles + blind_handles, day_labels + group_labels + blind_labels, loc='upper right')

plt.show()