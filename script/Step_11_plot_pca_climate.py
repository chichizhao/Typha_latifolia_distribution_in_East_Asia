#!/bin/bash
# function: plot PCA for 3 populations with climate data west  east and middle
# usage: python plot_pca_climate.py
# Author: chichi
# DATA
# data set 
# west
west1=[5.0625,11.025,23.6588,1363.15,28.8,-17.8,46.6,19.1333,-10.3333,20.8667,-12.35,245,39,9,38.5787,91,31,90,39]
west2=[7,12.85,27.5161,1257.62,30.4,-16.3,46.7,19.9,-7.03333,21.7833,-8.95,228,32,10,36.0555,89,33,77,36]
west3=[6.04583,11.4583,23.9714,1379.24,29.2,-18.6,47.8,19.95,-9.55,21.8167,-11.8167,227,34,7,44.1082,93,26,82,27]
west4=[4.01667,13.1667,29.0015,1201.79,26.2,-19.2,45.4,16.0833,-9.8,17.9,-11.5667,314,43,14,33.4372,116,50,103,50]
west5=[6.57917,12.7917,27.9294,1202.94,27.5,-18.3,45.8,18.3833,-9.73333,19.7167,-9.73333,232,36,9,45.3339,97,29,92,29]
west6=[8.07917,11.275,21.8085,1521,31.8,-19.9,51.7,23.2167,-8.75,24.7167,-12.55,159,22,6,41.7985,63,20,56,20]
west7=[8.07917,11.275,21.8085,1521,31.8,-19.9,51.7,23.2167,-8.75,24.7167,-12.55,159,22,6,41.7985,63,20,56,20]
west8 =[4.320833,12.925,27.73605,1218.092,26.2,-20.4,46.6,16.15,-12,17.76667,-12,273,43,8,50.08178,118,28,111,28]
# East
East7=[9.3375,9.458333,23.58687,1127.973,27.7,-12.4,40.1,22.5,-5.316667,22.5,-5.316667,930,255,11,101.6927,572,39,572,39]
East8=[5.6875,11.70833,26.07647,1241.594,26.2,-18.7,44.9,20.06667,-10.58333,20.06667,-10.58333,555,130,3,94.21339,325,13,325,13]
East9=[1.43333,12.5167,23.179,1527.43,25.5,-28.5,54,19.0333,-18.6333,19.0333,-18.6333,621,160,5,105.082,402,18,402,18]
East10=[3.225,11.28333,22.47676,1449.451,26.6,-23.6,50.2,19.96667,-15.8,19.96667,-15.8,568,125,5,89.34017,330,22,330,22]
East11=[0.754167,11.90833,21.41787,1640.22,26.6,-29,55.6,19.63333,-16.3,19.63333,-20.71667,543,130,4,98.63185,336,18,336,19]
East12=[2.7875,13.30833,25.64226,1465.51,27,-24.9,51.9,19.96667,-16.15,19.96667,-16.15,446,144,2,122.0104,321,7,321,7]
East13=[1.7,12.5167,22.6751,1597.67,26.9,-28.3,55.2,20.1333,-15.0667,20.1333,-19.3167,473,145,2,115.977,325,10,325,10]
East14=[1.55417,10.8417,20.3408,1567.94,26.7,-26.6,53.3,19.6,-19.1167,19.6,-19.1167,632,140,8,85.9582,351,30,351,30]
East15=[2.483333,11.53333,22.05226,1520.437,27.1,-25.2,52.3,20.1,-17.45,20.1,-17.45,594,140,4,94.55128,351,18,351,18]
East16=[-3.71667,14.7833,25.621,1634.83,23.2,-34.5,57.7,15.15,-21.3,15.15,-25.0167,449,129,3,114.708,313,11,313,12]
East17=[2.4875,11.60833,22.11111,1528.533,27.2,-25.3,52.5,20.16667,-17.56667,20.16667,-17.56667,589,140,5,95.98936,354,19,354,19]
East18=[-4.55,14.86667,25.99067,1626.907,23.3,-33.9,57.2,14.81667,-24.93333,14.81667,-24.93333,518,133,5,102.5903,334,17,334,17]
East19=[-4.15,15.6,26.53061,1639.818,23.2,-35.6,58.8,14.8,-21.75,14.8,-25.46667,467,133,4,112.8864,322,13,322,14]



# Middle
Middle1=[4.9375,10.4917,33.413,727.384,19.5,-11.9,31.4,12.85,-4.51667,13.5167,-4.51667,682,112,3,77.9401,319,15,315,15]
Middle2=[15.0833,10.0333,30.2209,868.537,31.3,-1.9,33.2,24.45,4.06667,25.65,4.06667,809,150,9,68.5778,382,33,350,33]
Middle3=[13.2375,10.2417,29.1785,923.278,30.4,-4.7,35.1,23.0333,1.56667,24.5667,1.56667,631,111,6,69.5901,301,22,261,22]

# Baikal Lake
Baikal1 = [-0.354167,11.95833,23.72685,1438.815,24.3,-26.1,50.4,16.56667,-15.85,16.56667,-18.98333,410,101,6,92.74624,249,22,249,29]
Baikal2 = [-1.375,10.96667,22.6584,1379.278,21.9,-26.5,48.4,14.88333,-15.93333,14.88333,-19.26667,352,89,4,99.29326,224,15,224,18]
Baikal3 = [-0.15,11,24.22907,1274.188,21.8,-23.6,45.4,14.93333,-14.33333,14.93333,-16.58333,407,101,6,93.71016,248,21,248,25]


# PCA analysis
# Step 1: Import the necessary libraries
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
import plotly.express as px

# Step 2: Combine all the lists into a single data frame
data = np.array([west1, west2, west3, west4, west5, west6, west7,west8, East7, East8, East9, East10, East11, East12, East13, East14, East15, East16, East17, East18, East19, Middle1, Middle2, Middle3, Baikal1, Baikal2, Baikal3])
df = pd.DataFrame(data)
# add the group name to the data frame
df['group'] = ['WA','WA', 'WA', 'WC', 'WC', 'WB','WB','WC', 'EA', 'EA', 'EB', 'EB', 'EB', 'EB', 'EB', 'EB', 'EC', 'EC', 'EC', 'EC', 'EC', 'MA', 'MA', 'MA', 'BA', 'BA', 'BA']
#print(df)          
# add the column name to the data frame
df.columns = ['bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19','group']
#print(df)
# pca analysis
features = ['bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19']   
#print(df)
# write the data to a csv file
df.to_csv('PCA_East_Asia_climate_TOP_3_8pops.csv', index=False)
# standardize the data
scaler = StandardScaler()
df_scaled = scaler.fit_transform(df[features])
#print(df_scaled)
# apply PCA
pca = PCA()
# fit PCA with the data
components = pca.fit_transform(df[features])
#print(components)
# Get the feature importance
importance = np.abs(pca.components_[0])


# Sort the feature labels based on the importance
feature_names = np.array(features)
sorted_idx = np.argsort(importance)
sorted_feature_names = feature_names[sorted_idx]

# Create labels with PC number, variance and original variable name
labels = {
    str(i): f"PC {i+1} ({var:.2f}%) - {sorted_feature_names[i]}"
    for i, var in enumerate(pca.explained_variance_ratio_ * 100)
}
print(labels)
#print(components)
# here we want plot the pca result with matplotlib
import matplotlib.pyplot as plt
# plot the pca result with matplotlib
#colors = {'west':'red', 'East':'blue', 'Middle':'green', 'Baikal':'gray'}
colors = { 'WA':'#F3740B', 'WB':'#EFA38A', 'WC':'#D2352C','BA':'#808080','EC':'#ADD9ED', 'EB':'#30BBCE','EA':'#4A7DB4','MA':'#68AC56'}
# color_set=['#30BBCE','#F3740B','#D2352C','#4A7DB4','#68AC56','#808080','#EFA38A','#ADD9ED']

# create the lagend plot for the pca result
fig, ax = plt.subplots()
fig.set_size_inches(2.5,2.5)
ax.set_xlim(0.2,7.6)
ax.set_ylim(8.2,15.6)
# set the  font in arail bold
plt.rcParams['font.sans-serif'] = ['Arial']
# plt.rcParams['font.weight'] = 'bold'
# aff the ticks of the subplots
ax.set_xticks([])
ax.set_yticks([])
# off the margin of the subplots
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
# add the pca result
for i in range(3):
    for j in range(3):
        if i < j: 
        # in the origial plot, we have 19*19 subplots, so we need to change the position of the subplots
            standard_i = max(components.__abs__()[:,i])
            #print(standard_i)
            standard_j = max(components.__abs__()[:,j])
            print(standard_j)
            # according to the standard_i and standard_j, we can change the position of the subplots into this plot
            new_i = (14-4*i)+ 1.5*components[:,i]/standard_i
            new_j = (2+4*(j-i-1))+ 1.5*components[:,j]/standard_j
            print(new_i)
            print(new_j)
            ax.scatter(new_j, new_i, c=df['group'].apply(lambda x: colors[x]), s=5)
            # add box
            box = plt.Rectangle((0.4+4*(j-i-1),12.4-4*i),3.2,3.2,fc='none',ec='k',lw=1)
            ax.add_artist(box)
            # add the 0 line
            line = plt.Line2D([0.4+4*(j-i-1),3.6+4*(j-i-1)],[14-4*i,14-4*i],color='k',linestyle='--',linewidth=0.5,alpha=0.4)
            ax.add_artist(line)
            line = plt.Line2D([2+4*(j-i-1),2+4*(j-i-1)],[12.4-4*i,15.6-4*i],color='k',linestyle='--',linewidth=0.5,alpha=0.4)
            ax.add_artist(line)
            # add the label of the pca
            ax.text(2+4*(j-i-1),12.2-4*i-0.1,labels[str(j)],fontsize=6,horizontalalignment='center',verticalalignment='center')
            ax.text(0.2+4*(j-i-1)-0.1,14-4*i,labels[str(i)],fontsize=6,horizontalalignment='center',verticalalignment='center',rotation=90)

#add the west, east and middle legend
#box = plt.Rectangle((4.4,8.4),3.2,3.2,fc='none',ec='k',lw=1)
#ax.add_artist(box)
for i in range(8):
    if i%2 == 0:
        circle = plt.Circle((5,8.7+(i//2)*0.8), 0.1, color=list(colors.values())[i])
        ax.add_artist(circle)
        ax.text(5.4,8.68+(i//2)*0.8,list(colors.keys())[i],fontsize=6,verticalalignment='center')
    else:
        circle = plt.Circle((6.5,8.7+(i//2)*0.8), 0.1, color=list(colors.values())[i])
        ax.add_artist(circle)
        ax.text(6.9,8.68+(i//2)*0.8,list(colors.keys())[i],fontsize=6,verticalalignment='center')
 

# save the figure
plt.savefig('PCA_East_Asia_climate_TOP_3_8pops.pdf', dpi=600, bbox_inches='tight')



