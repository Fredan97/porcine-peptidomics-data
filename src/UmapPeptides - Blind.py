from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import pandas as pd
import umap
import math

#%%
reducer = umap.UMAP(random_state=3)
#%%
df = pd.read_excel("../data/Data Day 1 and 2 and 3 OnlyPep NoDups w Blind.xlsx")
data = df.iloc[3:,1:]
group = df.iloc[1,1:]
day = df.iloc[0,1:]
blinds = df.iloc[2,1:]

#%%Change data into log

for x in range(data.shape[0]):
    for y in range(data.shape[1]):
        if data.iloc[x][y] != 0:
            data.iloc[x,y] = math.log(data.iloc[x][y],2)

#%% Transpose data
transpdata=data.transpose()


#%% Perform dimensionality reduction
scaled_data = StandardScaler().fit_transform(transpdata)
embedding = reducer.fit_transform(scaled_data)

#%% Fix font sizes and styles
plt.rcParams['font.size']=7
plt.rcParams["font.family"] = "Arial"

#%%
# Prepare colors and markers
colors = group.map({"S.a": "gold", "P.a": "aqua", "Ctrl": "black", "Double": "lime", "Acc Double": "red"}).tolist()
edgeColors = blinds.map({"No": "w", "Yes": "red"}).tolist()
markers = day.map({1: "o", 2: "s", 3: "X"}).tolist()

# Create the scatter plot
plt.figure(figsize=(5,5))
unique_days = [1, 2, 3]
marker_map = {1: "o", 2: "s", 3: "X"}

# Plot each day's data
for unique_day in unique_days:
    idx = [i for i, d in enumerate(day) if d == unique_day]
    plt.scatter(
        embedding[idx, 0],
        embedding[idx, 1],
        c=[colors[i] for i in idx],
        label=f'Day {unique_day}',
        marker=marker_map[unique_day],
        edgecolor=[edgeColors[i] for i in idx],
        s=30
    )

handlecolors = ["black","gold","aqua","lime","red"]
# Create a legend for days and markers
day_handles = [plt.Line2D([0,0],[0,0],color='gray', marker=marker_map[day], linestyle='') for day in unique_days]
day_labels = [f'Day {day}' for day in unique_days]
group_labels = list(["Contorl","S.a", "P.a", "Double infection", "Accidental double infection"])
group_handles = [plt.Line2D([0,0],[0,0],color=handlecolors[i], marker='o', linestyle='',markersize=5) for i in [0,1,2,3,4]]
blind_labels = list({"Blinded"})
blind_handles = [plt.Line2D([0,0],[0,0],color='gray',marker = 'o',markeredgecolor = 'red', linestyle ='',markersize=5)]
plt.legend(day_handles + group_handles + blind_handles, day_labels + group_labels + blind_labels, loc='lower left')

# Add title and labels
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')

plt.show()

