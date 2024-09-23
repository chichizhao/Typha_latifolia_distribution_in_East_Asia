#!/bin/python
# -*- coding: utf-8 -*-
# Author: chichi
# function: plot the demography with momi estimated parameters
# usage: python plot_domography.py

import matplotlib.pyplot as plt
from math import log10 as lg


fig, ax = plt.subplots(1, 1, figsize=(16, 12))
# set the font style in arial
plt.rcParams["font.family"] = "arial"

# set the y axis scale
ax.set_ylim([-50000, 6.5e5])

# set the y axis in exponential format
region = ["West China", "Baikal", "East China", "Middle China"]
color_set = ['#EB100C','#73786A','#4BC2C9','#6AA84F']

# add the line of 598432
line = plt.Line2D([10, 75], [598432, 598432], color="#9D9788", linestyle="--", linewidth=1, alpha=0.5)
ax.add_line(line)

# add the line of 394032
line = plt.Line2D([8, 75], [394032, 394032], color="#C78770", linestyle="--", linewidth=1, alpha=0.5)
ax.add_line(line)

# add the line of 331514
line = plt.Line2D([10, 75], [331514, 331514], color="#C49C81", linestyle="--", linewidth=1, alpha=0.5)
ax.add_line(line)

# add the line of 291830
line = plt.Line2D([20, 75], [291830, 291830], color="#B59A7D", linestyle="--", linewidth=1, alpha=0.5)
ax.add_line(line)

# add the line of 83689
line = plt.Line2D([60, 75], [83689, 83689], color="#4EB2DD", linestyle="--", linewidth=1, alpha=0.5)
ax.add_line(line)

# add the line of 68976
line = plt.Line2D([34, 75], [68976, 68976], color="#B88570", linestyle="--", linewidth=1, alpha=0.5)
ax.add_line(line)

# add the line of 25073
line = plt.Line2D([58, 75], [25073, 25073], color="#69B0D9", linestyle="--", linewidth=1, alpha=0.5)
ax.add_line(line)


# add the region name
box = plt.Rectangle((0, -40000), lg(97265)+5+lg(109422)+5+lg(126084), 10000, fc="#EB100C", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((0, -40000), lg(97265)+5+lg(109422)+5+lg(126084), 60000, fc="#EB100C", edgecolor="none", alpha=0.1)
ax.add_patch(box)
ax.text((lg(97265)+5+lg(109422)+5+lg(126084))/2, -55000, "Northeast China", ha='center', va='center', fontsize=14, color='black')

box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5, -40000), lg(264083), 10000, fc="#73786A", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5, -40000), lg(264083),  60000, fc="#73786A", edgecolor="none", alpha=0.1)
ax.add_patch(box)
ax.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)/2, -55000, "Baikal Lake", ha='center', va='center', fontsize=14, color='black')

box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5, -40000), lg(25729)+5+lg(224581)+lg(2440)+5, 10000, fc="#4BC2C9", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5, -40000), lg(25729)+5+lg(224581)+lg(2440)+5, 60000, fc="#4BC2C9", edgecolor="none", alpha=0.1)
ax.add_patch(box)
ax.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)/2, -55000, "Northeast China", ha='center', va='center', fontsize=14, color='black')

box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5+lg(2440)+5, -40000), lg(22567), 10000, fc="#6AA84F", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5+lg(2440)+5, -40000), lg(22567), 60000, fc="#6AA84F", edgecolor="none", alpha=0.1)
ax.add_patch(box)

ax.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5+lg(2440)+5+lg(22567)/2, -55000, "Central China", ha='center', va='center', fontsize=14, color='black')

# add the line of 598432

color_set=['#30BBCE','#F3740B','#D2352C','#4A7DB4','#68AC56','#808080','#EFA38A','#ADD9ED']

# set the x axis scale
ax.set_xlim([0, 98])
box = plt.Rectangle((0, 0),lg(97265), 394032, fc="#D2352C", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5, 0),lg(109422), 331514, fc="#EFA38A", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5, 0),lg(126084), 291830, fc="#F3740B", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5, 0),lg(264083), 68976, fc="#808080", edgecolor="none") 
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5, 0),lg(25729), 68976, fc="#ADD9ED", edgecolor="none")
ax.add_patch(box)

box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5, 0),lg(224581), 25073, fc="#30BBCE", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5, 0),lg(2440), 25073, fc="#6BA9EC", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5+lg(2440)+5, 0),lg(22567), 83689, fc="#68AC56", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle(((lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2)/2,394032),lg(630200), 598432-394032, fc="#C78770", edgecolor="none")
ax.add_patch(box)
# add the link between the boxes
box = plt.Rectangle((0, 394032),lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(1e8)/4-lg(10016663)/2)/2+lg(71278)/2, 5000, fc="#CCCCCC", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2,331514),lg(71278), 394032-331514, fc="#C49C81", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5, 329700),lg(109422)+5+lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4+lg(10016663)/2, 5000, fc="#CCCCCC", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2, 291830),lg(10016663), 331514-291830, fc="#B59A7D", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5, 291830),lg(126084)+5+lg(264083)+2.5+lg(10e8)/2, 5000, fc="#CCCCCC", edgecolor="none")
ax.add_patch(box)

box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+2.5-lg(10e8)/2, 68976),lg(10e8), 291830-68976, fc="#97ADB7", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5, 63976),lg(25729)+lg(264083)+5, 5000, fc="#CCCCCC", edgecolor="none")
ax.add_patch(box)

box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+2.5-lg(673200)/2, 25073),lg(673200), 83689-25073, fc="#4EB2DD", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5, 24073),lg(2440)+5+lg(224581), 5000, fc="#CCCCCC", edgecolor="none") 
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+2.5+lg(673200)/2 +(5+lg(2440)+5-2.5-lg(673200)/2 )/2-lg(139479300)/2, 83689),lg(139479300), 598432-83689, fc="#56B0B0", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+2.5-lg(673200)/2, 80689),5+lg(2440)+5-2.5+lg(673200)/2+lg(22567), 5000, fc="#CCCCCC", edgecolor="none")
ax.add_patch(box)

box = plt.Rectangle((((lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2)/2+((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+2.5+lg(673200)/2 +(5+lg(2440)+5-2.5-lg(673200)/2 )/2-lg(139479300)/2)-(lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2)/2+lg(630200)/2))/2,598432),lg(633490), 68976, fc="#9D9788", edgecolor="none")
ax.add_patch(box)
box = plt.Rectangle(((lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2)/2, 598432),(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+2.5+lg(673200)/2 +(5+lg(2440)+5-2.5-lg(673200)/2 )/2-lg(139479300)/2)-(lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2)/2+lg(139479300), 5000, fc="#CCCCCC", edgecolor="none")
ax.add_patch(box)
#box = plt.Rectangle((97265,382622),1504812.51, 689760-382622, fc="#CCCCCC", edgecolor="none")
#ax.add_patch(box)

# add the imgration arrow
arrow = plt.Arrow(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+4.95, 13216, -4.8, 0, width=15000, color="black")
ax.add_patch(arrow)
#add the text
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+1, 13216+5000, "18.4%", fontsize=10)
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729), 13216-15000, "13.2 kya", fontsize=10)
arrow = plt.Arrow(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5+lg(2440)+4.95, 23533, -4.8, 0, width=15000, color="black")

ax.add_patch(arrow)
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5+lg(2440)+1, 23533+5000, "33.0%", fontsize=10)
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5+lg(2440), 23533-15000, "23.5 kya", fontsize=10)
# off the axis
ax.axis("off")

# add the group name
plt.text(0+1, -23000, "WC", fontsize=12)
plt.text(lg(97265)+5+1, -23000, "WB", fontsize=12)
plt.text(lg(97265)+5+1+lg(109422)+5, -23000, "WA", fontsize=12)
plt.text(lg(97265)+5+2+lg(109422)+5+lg(126084)+5, -23000, "BA", fontsize=12)
plt.text(lg(97265)+5+2+lg(109422)+5+lg(126084)+5+lg(264083)+5, -23000, "EC", fontsize=12)
plt.text(lg(97265)+5+2+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5, -23000, "EB", fontsize=12)
plt.text(lg(97265)+5+1+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5, -23000, "EA", fontsize=12)
plt.text(lg(97265)+5+1.5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5+lg(2440)+5, -23000, "MA", fontsize=12)



# add the number of the population
plt.text(lg(97265)/2, 394032/2, "97,625", fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)/2, 331514/2, "109,422", fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)/2, 291830/2, "126,084", fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)/2, 68976/2, "264,083", fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)/2, 68976/2, "25,729", fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)/2, 25073/2, "224,581", fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5+lg(2440)/2, 25073/2, "2,440", fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+5+lg(2440)+5+lg(22567)/2, 83689/2, "22,567", fontsize=10, horizontalalignment="center")

plt.text(((lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2)/2+((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+2.5+lg(673200)/2 +(5+lg(2440)+5-2.5-lg(673200)/2 )/2-lg(139479300)/2)-(lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2)/2+lg(630200)/2))/2+lg(633490)/2,598432+68976/2,"244,141",fontsize=10, horizontalalignment="center")
plt.text(((lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2)/2+((lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+2.5+lg(673200)/2 +(5+lg(2440)+5-2.5-lg(673200)/2 )/2-lg(139479300)/2)-(lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2)/2+lg(630200)/2))/2+lg(633490)/2,583432,"598 kya",fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+2.5+lg(673200)/2 +(5+lg(2440)+5-2.5-lg(673200)/2 )/2,(598432-83689)/2+83689,"1,278,659",fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+2.5+lg(673200)/2 +(5+lg(2440)+5-2.5-lg(673200)/2 )/2,68689,"84 kya",fontsize=10, horizontalalignment="center")

plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+2.5,(83689-25073)/2+25073,"49,491", fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+5+lg(25729)+5+lg(224581)+2.5,14073,"25 kya", fontsize=10, horizontalalignment="center")

plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+2.5,(291830-68976)/2+68976,"100,000,000", fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+5+lg(264083)+2.5,48976,"69 kya", fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4,(331514-291830)/2+291830,"10,016,663",fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+5+lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4,276830,"292 kya",fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2,(394032-331514)/2+331514, "71,278",fontsize=10, horizontalalignment="center")
plt.text(lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2,316514, "332 kya",fontsize=10, horizontalalignment="center")

plt.text((lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2)/2+lg(630200)/2,(580103-394032)/2+394032,"144,569",fontsize=10, horizontalalignment="center")
plt.text((lg(97265)+5+lg(109422)+2.5+(lg(126084)+2.5+lg(264083)/2+1.25-lg(10e8)/4-lg(10016663)/2)/2-lg(71278)/2)/2+lg(630200)/2,379032,"394 kya",fontsize=10, horizontalalignment="center")

# add the history temperature data
import pandas as pd
import numpy as np

data = pd.read_csv('/home/chichi/data/china/china3/momi/EDC_dD_temp_estim.tab', sep='\t', header=0)
# print(data.head())
# remove the rows with missing values
data = data.dropna()
print(data.head())
age = data['Age model [ka]']*1000
T = data['delta T [°C]']*1
box = plt.Rectangle((75, 0), 20, 650000, fc="none", edgecolor="black", linewidth=1)
ax.add_patch(box)
line = plt.Line2D( 82-T, age, color="red", linestyle="-", linewidth=0.5, alpha=0.5)
ax.add_line(line)
ax.text(85, -50000, "  Antarctic Temperature (°C) \n relative to present", fontsize=14, color="black", ha="center", va="center")
box = plt.Rectangle((77, -8000), 0.1, 8000, color="black")
ax.add_patch(box)
ax.text(77, -20000, "5", fontsize=14, color="black", ha="center", va="center")
box = plt.Rectangle((82, -8000), 0.1, 8000, color="black")
ax.add_patch(box)
ax.text(82, -20000, "0", fontsize=14, color="black", ha="center", va="center")
box = plt.Rectangle((87, -8000), 0.1, 8000, color="black")
ax.add_patch(box)
ax.text(87, -20000, "-5", fontsize=14, color="black", ha="center", va="center")
box = plt.Rectangle((92, -8000), 0.1, 8000, color="black")
ax.add_patch(box)
ax.text(92, -20000, "-10", fontsize=14, color="black", ha="center", va="center")

ax.add_patch(box)
for i in range(0, 7):
    box = plt.Rectangle((95, i*100000), 0.8, 1000, color="black")
    ax.add_patch(box)
    ax.text(96.5, i*100000, str(i*100), fontsize=14, color="black", ha="center", va="center",rotation=90)
ax.text(98.5, 300000, "Times before Prensent (kya) (Jouzel et al., 2007)", fontsize=14, color="black", ha="center", va="center",rotation=90)

# the 0 celsius line
line = plt.Line2D([82, 82], [0, 700000], color="black", linestyle="--", linewidth=0.5, alpha=0.2)
ax.add_line(line)

# add the legend line of 598 kya

arrow = plt.Arrow(75, 598432, 12, 0, width=10000, color="#A29D8E")
ax.add_patch(arrow)
ax.text(77, 598432+5000, "598 kya", fontsize=12, color="#A29D8E")
# add the line of y = 394 kya
arrow = plt.Arrow(75, 394032, 7.5, 0, width=10000, color="#C14336")
ax.add_patch(arrow)
ax.text(75.5, 394032-12000, "394 kya", fontsize=12, color="#C14336")
# add the line of y = 332 kya
arrow = plt.Arrow(75, 331514, 3.2, 0, width=10000, color="#E4A68F")
ax.add_patch(arrow)
ax.text(75.5, 331514+8000, "332 kya", fontsize=12, color="#E4A68F")
# add the line of y = 291 kya
arrow = plt.Arrow(75, 291830, 10.5, 0, width=10000, color="#E47B31")
ax.add_patch(arrow)
ax.text(77, 291830+5000, "291 kya", fontsize=12, color="#E47B31")
# add the line of y = 84 kya
#line = plt.Line2D([75, 82+2.4], [83689, 83689], color="#6DAEAF", linestyle="-", linewidth=2)
arrow = plt.Arrow(75, 83689, 8.5, 0, width=10000, color="#6DAEAF")
ax.add_patch(arrow)
ax.text(77, 83689+5000, "84 kya", fontsize=12, color="#6DAEAF")
# add the line of y = 69 kya
arrow = plt.Arrow(75, 68976, 12, 0, width=10000, color="#B88570")
ax.add_patch(arrow)
ax.text(77, 68976-15000, "69 kya", fontsize=12, color="#B88570")

#add the line of y = 25 kya
arrow = plt.Arrow(75, 25073, 15, 0, width=10000, color="#69B0D9")
ax.add_patch(arrow)
ax.text(77, 25073+5000, "25 kya", fontsize=12, color="#69B0D9")

plt.savefig("momi_domography3.pdf", dpi=600, bbox_inches="tight")
