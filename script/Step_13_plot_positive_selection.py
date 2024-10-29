#!/bin/py
# -*- coding: utf-8 -*-
# function: plot the RAISD results with manhattan plot
# usage: python plot_raisd.py -i raisd_result -o output_prefix
# author: chichi

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl



chr_len ={'Chr01': 24370473,'Chr02': 18857030,'Chr03': 18714887,'Chr04': 18467694,'Chr05': 17512124,'Chr06': 14553481,'Chr07': 14503471,'Chr08': 13751786,'Chr09': 11865621,'Chr10': 11839088,'Chr11': 11614203,'Chr12': 10714389,'Chr13': 10242974,'Chr14': 9831573,'Chr15': 9461071,'Contig1': 452046,'Contig2': 203240}
chr_start_pos= {'Chr01': 0,'Chr02': 24370473,'Chr03': 43227503,'Chr04': 61942390,'Chr05': 80410084,'Chr06': 97922208,'Chr07': 112775689,'Chr08': 127279160,'Chr09': 141030946,'Chr10': 152896567,'Chr11': 164735655,'Chr12': 176349858,'Chr13': 187064247,'Chr14': 197307221,'Chr15': 207138794,'Contig1': 216599865,'Contig2': 261846911}


#def plot_raisd(input_file, output_prefix):  
data ={}
with open("/home/chichi/data/china/china3/raisd/RAiSD_Report.MA_clean") as f:
    for line in f:
        line = line.strip()

        if line.startswith("// Ch"):
            chr = line.split( )[1]
            print(chr)
            data[line.split( )[1]] = []
            print(chr)
        elif line.startswith("// Contig"):
            break
        else:
            data[chr].append(line.split( ))
print(data.keys())
"""
# here need the 99% snp cutoff line
snp_values = []
for chr in data:
    for i in range(len(data[chr])):
        snp_values.append(float(data[chr][i][6]))
snp_values.sort()
print(snp_values[int(len(snp_values)*0.999)])
# 99.9% snp cutoff line 23.45
print(snp_values[int(len(snp_values)*0.995)])
# 99.5% snp cutoff line 13.41
print(snp_values[int(len(snp_values)*0.99)])
# 99% snp cutoff line 10.48
print(snp_values[int(len(snp_values)*0.95)])
# 95% snp cutoff line 4.838
"""




# plot the manhattan plot for var of each chr
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
ax.set_xlim(-3500000, 220599865)
ax.set_ylim(-5, 880)
box = plt.Rectangle((0, -3), 216599865, 1, fc='#68AC56', ec='#68AC56', lw=0.5)
ax.add_patch(box)
num = 0
ax.axis('off')
# here we need mark the significant snp with their name 
#for chr in data:
gff_file = "/home/chichi/data/china/china3/gwas/pasa2.longest.filter.gff3.cds"
gff_data = pd.read_csv(gff_file, sep="\t", header=None)
print(gff_data.head())

# read the gene annotation file for the genome
annotation = "/home/chichi/data/china/china3/gwas/func_merg.xls"
annotation_data = pd.read_csv(annotation, sep='\t',index_col= False)
annotation_dict = {}
for i in range(len(annotation_data)):
    if str(annotation_data.loc[i, 'PFAM']).startswith('PF'):
        annotation_dict[annotation_data.loc[i, '#ID']] = [annotation_data.loc[i, 'PFAM'].split('(')[0], annotation_data.loc[i, 'PFAM'], annotation_data.loc[i, 'GO'], annotation_data.loc[i, 'KEGG'], annotation_data.loc[i, 'NR'], annotation_data.loc[i, 'Swissprot']]
    else:
        annotation_dict[annotation_data.loc[i, '#ID']] = [annotation_data.loc[i, 'PFAM'], annotation_data.loc[i, 'PFAM'], annotation_data.loc[i, 'GO'], annotation_data.loc[i, 'KEGG'], annotation_data.loc[i, 'NR'], annotation_data.loc[i, 'Swissprot']]
# print(annotation_dict.keys())

MA_sig_gene = []

for chr in data:
#for chr in {'Chr01'}:
    chr_x = []
    chr_y = []
    chr_x_sig = []
    chr_y_sig = []
    chr_x_sig_marker = []
    chr_y_sig_marker = []

    gff_data_chr = gff_data[gff_data[0] == chr]
    print(gff_data_chr.head())
    for i in range(len(data[chr])):
        chr_x.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
        if float(data[chr][i][6]) > 10.48:
            chr_x_sig.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
            chr_y_sig.append(float(data[chr][i][6]))
            # to determine the significant snp in which gene or not
            if float(data[chr][i][6]) > 10.48:
                for index, row in gff_data_chr.iterrows():
                    if int(row[3]) <= int(data[chr][i][0]) <= int(row[4]):
                        #print(row[3], row[4], data[chr][i][0])
                        pos = int(data[chr][i][0])
                        id = gff_data_chr.loc[index, 8].split('Parent=')[1]

                        if str(id) in annotation_dict:
                            print(annotation_dict[id])
                            MA_sig_gene.append([chr, pos, str(id), annotation_dict[id][0], annotation_dict[id][1], annotation_dict[id][2], annotation_dict[id][3], annotation_dict[id][4], annotation_dict[id][5]])
                        
                        chr_x_sig_marker.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
                        chr_y_sig_marker.append(float(data[chr][i][6]))
                        break
        
        chr_y.append(float(data[chr][i][6]))


    num += 1# add the chromosome length scale
    if num % 2 == 0:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#68AC56", alpha=1, rasterized=True)
    else:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#68AC56", alpha=1, rasterized=True)
    #if chr == "Chr02":
    #    break

# add the significant line 10.04
#ax.scatter(x, y, s=1, label=chr, color = "grey", alpha=0.2, rasterized=True)
#ax.scatter(x_sig, y_sig, s=1, color = "grey", alpha=0.5, rasterized=True)
#ax.scatter(x_sig_marker, y_sig_marker, s=3, color = "#68AC56", alpha=1, rasterized=True)       
#ax.axhline(y=10.48, color='b', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [10.48, 10.48], color='b', linestyle='--', linewidth=0.5, alpha=0.5)
#ax.axhline(y=23.45, color='r', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [23.45, 23.45], color='r', linestyle='--', linewidth=0.5, alpha=0.5)


for i in range(1, 16):
    if i < 10:
        chr = 'Chr0' + str(i)
    else:
        chr = 'Chr' + str(i)
    ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -20, str(i), ha='center', va='center', fontsize=10, color='black')
    # add the chromosome region line 
    #ax.axvline(x=chr_start_pos[chr], color='black', linestyle='-', alpha=0.2)
    if i % 2 == 1:
        # ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -10, str(i), ha='center', va='center', fontsize=7, color='black')
        box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 870, fc='grey', ec='none', lw=0.5, alpha=0.2)
        ax.add_patch(box)
    #box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 100, fc='grey', ec='none', lw=0.5, alpha=0.2)
#ax.axvline(x=chr_start_pos['Contig1'], color='black', linestyle='-')
#ax.axvline(x=chr_start_pos["Chr01"], color='black', linestyle='-')
box = plt.Rectangle((-3500000, 0), 1, 90, fc='black', ec='black', lw=1)
ax.add_patch(box)
#ax.axhline(y=14, color='black', linestyle='--', alpha=0.5)

# add the y axis label
#ax.text(-1800000, 10.5, "99%", ha='right', va='center', fontsize=10, color='blue')
#ax.text(-1800000, 30, "99.9%", ha='right', va='center', fontsize=10, color='red')
ax.text(-12500000, 400, "-log10(µ)", ha='center', va='center', fontsize=7, color='black', rotation= 90)
ax.text(-1000000, -20, "Chr", ha='center', va='center', fontsize=10, color='black')
ax.text(-3700000, 0, "0", ha='right', va='center', fontsize=7, color='black')
#ax.text(-200000, 20, "20", ha='right', va='center', fontsize=7)
ax.text(-3700000, 40, "40", ha='right', va='center', fontsize=7)
#ax.text(-200000, 60, "60", ha='right', va='center', fontsize=10)
ax.text(-3700000, 80, "80", ha='right', va='center', fontsize=7)
#ax.text(-200000, 100, "100", ha='right', va='center', fontsize=10)

box = plt.Rectangle((216599865, 0), 4000000, 100, fc='#68AC56', ec='none', lw=0.5, alpha=0.7)
ax.add_patch(box)
ax.text(216599865 + 2000000, 50, "MA", ha='center', va='center', fontsize=7, color='black', rotation= 90)

# write the gene annotation in to file
# sort the gene by the significant snp
MA_sig_gene = sorted(MA_sig_gene, key=lambda x: (x[0], int(x[1])))
with open("MA_sig_gene.xls", 'w') as f:
    f.write("chr\tpos\tgene_id\tPFAM_family\tPFAM\tGO\tKEGG\tNR\tSwissprot\tp_wald\n")
    for item in MA_sig_gene:
        #f.write("%s\n" % str(item).replace('[','').replace(']','').replace('\'',''))
        #f.write("%s\n" % '\t'.join(map(str, item)))
        f.write("%s\n" % '\t'.join(map(str, item)))


data ={}
with open("/home/chichi/data/china/china3/raisd/RAiSD_Report.EA_clean") as f:
    for line in f:
        line = line.strip()
        if line.startswith("// Ch"):
            chr = line.split( )[1]
            print(chr)
            data[line.split( )[1]] = []
            print(chr)
        elif line.startswith("// Contig"):
            break
        else:
            data[chr].append(line.split( ))
print(data.keys())
"""
# here need the 99% snp cutoff line
snp_values = []
for chr in data:
    for i in range(len(data[chr])):
        snp_values.append(float(data[chr][i][6]))
snp_values.sort()
print(snp_values[int(len(snp_values)*0.999)])
# 99.9% snp cutoff line 17.71
print(snp_values[int(len(snp_values)*0.995)])
# 99.5% snp cutoff line 9.233
print(snp_values[int(len(snp_values)*0.99)])
# 99% snp cutoff line 7.099
print(snp_values[int(len(snp_values)*0.95)])
# 95% snp cutoff line 3.668
"""
box = plt.Rectangle((0, 100), 216599865, 10, fc='white', ec='none', lw=0.5)
ax.add_patch(box)
box = plt.Rectangle((0, 107), 216599865, 1, fc='#4A7DB4', ec='#4A7DB4', lw=0.5)
ax.add_patch(box)


EA_sig_gene = []

for chr in data:
#for chr in {'Chr15'}:
    chr_x = []
    chr_y = []
    chr_x_sig = []
    chr_y_sig = []
    chr_x_sig_marker = []
    chr_y_sig_marker = []

    gff_data_chr = gff_data[gff_data[0] == chr]
    print(gff_data_chr.head())
    for i in range(len(data[chr])):
        
        chr_x.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
        if float(data[chr][i][6]) > 7.099:
        
            chr_x_sig.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
            chr_y_sig.append(float(data[chr][i][6])/80*100 + 100+10)
            # to determine the significant snp in which gene or not
            if float(data[chr][i][6]) > 7.099:
                for index, row in gff_data_chr.iterrows():
                    if int(row[3]) <= int(data[chr][i][0]) <= int(row[4]):
                        #print(row[3], row[4], data[chr][i][0]) 
                        pos = int(data[chr][i][0])
                        id = gff_data_chr.loc[index, 8].split('Parent=')[1]

                        if str(id) in annotation_dict:
                            print(annotation_dict[id])
                            EA_sig_gene.append([chr, pos, str(id), annotation_dict[id][0], annotation_dict[id][1], annotation_dict[id][2], annotation_dict[id][3], annotation_dict[id][4], annotation_dict[id][5]])
                                               
                        chr_x_sig_marker.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
                        chr_y_sig_marker.append(float(data[chr][i][6])/80*100 + 100+10)
                        break
        chr_y.append(float(data[chr][i][6])/80*100 + 100+10)


    num += 1# add the chromosome length scale
    #if chr == "Chr02":
    #    break
    if num % 2 == 1:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#4A7DB4", alpha=1, rasterized=True)
    else:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#4A7DB4", alpha=1, rasterized=True)
        
#ax.axhline(y=7.099/80*100+100+10, color='b', linestyle='--', linewidth=0.5, alpha=0.5)
#ax.axhline(y=17.71/80*100 + 100+10, color='r', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [7.099/80*100+100+10, 7.099/80*100+100+10], color='b', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [17.71/80*100 + 100+10, 17.71/80*100 + 100+10], color='r', linestyle='--', linewidth=0.5, alpha=0.5)
# print the number of significant snp
# print(len(positive_select_pos))
box = plt.Rectangle((-3500000, 110), 1, 100, fc='black', ec='black', lw=1)
ax.add_patch(box)


# add the y axis label
#ax.text(-11500000, 150, "-log10(µ)", ha='center', va='center', fontsize=7, color='black', rotation= 90)
#ax.text(0, -8, "chr", ha='center', va='center', fontsize=14, color='black')
ax.text(-3700000, 110, "0", ha='right', va='center', fontsize=7, color='black')
#ax.text(-200000, 120, "20", ha='right', va='center', fontsize=10)
ax.text(-3700000, 110+30/80*100, "30", ha='right', va='center', fontsize=7)
#ax.text(-200000, 160, "60", ha='right', va='center', fontsize=10)
ax.text(-3700000, 110 +60/80*100, "60", ha='right', va='center', fontsize=7)

box = plt.Rectangle((216599865, 110), 4000000, 100, fc='#4A7DB4', ec = 'none', lw=0.5, alpha=0.7)
ax.add_patch(box)
ax.text(216599865 + 2000000, 160, "EA", ha='center', va='center', fontsize=7, color='black', rotation= 90)

# write the gene annotation in to file
# sort the gene by the significant snp
EA_sig_gene = sorted(EA_sig_gene, key=lambda x: (x[0], int(x[1])))
with open("EA_sig_gene.xls", 'w') as f:
    f.write("chr\tpos\tgene_id\tPFAM_family\tPFAM\tGO\tKEGG\tNR\tSwissprot\tp_wald\n")
    for item in EA_sig_gene:
        #f.write("%s\n" % str(item).replace('[','').replace(']','').replace('\'',''))
        #f.write("%s\n" % '\t'.join(map(str, item)))
        f.write("%s\n" % '\t'.join(map(str, item)))
        

      

data ={}
with open("/home/chichi/data/china/china3/raisd/RAiSD_Report.EB_clean") as f:
    for line in f:
        line = line.strip()

        if line.startswith("// Ch"):
            chr = line.split( )[1]
            print(chr)
            data[line.split( )[1]] = []
            print(chr)
        elif line.startswith("// Contig"):
            break
        else:
            data[chr].append(line.split( ))
print(data.keys())
"""
# here need the 99% snp cutoff line
snp_values = []
for chr in data:
    for i in range(len(data[chr])):
        snp_values.append(float(data[chr][i][6]))
snp_values.sort()
print(snp_values[int(len(snp_values)*0.999)])
# 99.9% snp cutoff line 23.79
print(snp_values[int(len(snp_values)*0.995)])
# 99.5% snp cutoff line 13.68
print(snp_values[int(len(snp_values)*0.99)])
# 99% snp cutoff line 9.641
print(snp_values[int(len(snp_values)*0.95)])
# 95% snp cutoff line 4.319
"""

box = plt.Rectangle((0, 210), 216599865, 10, fc='white', ec='none', lw=0.5)
ax.add_patch(box)

box = plt.Rectangle((0, 217), 216599865, 1, fc='#30BBCE', ec='#30BBCE', lw=0.5)
ax.add_patch(box)

EB_sig_gene = []

for chr in data:
#for chr in {'Chr01'}:
    chr_x = []
    chr_y = []
    chr_x_sig = []
    chr_y_sig = []
    chr_x_sig_marker = []
    chr_y_sig_marker = []
    
    gff_data_chr = gff_data[gff_data[0] == chr]
    print(gff_data_chr.head())
    for i in range(len(data[chr])):
        chr_x.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
        if float(data[chr][i][6]) > 9.641:
            chr_x_sig.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
            chr_y_sig.append(float(data[chr][i][6])/80*100 + 200+20)
            # to determine the significant snp in which gene or not
            if float(data[chr][i][6]) > 9.641:
                for index, row in gff_data_chr.iterrows():
                    if int(row[3]) <= int(data[chr][i][0]) <= int(row[4]):
                        #print(row[3], row[4], data[chr][i][0])
                        pos = int(data[chr][i][0])
                        id = gff_data_chr.loc[index, 8].split('Parent=')[1]

                        if str(id) in annotation_dict:
                            print(annotation_dict[id])
                            EB_sig_gene.append([chr, pos, str(id), annotation_dict[id][0], annotation_dict[id][1], annotation_dict[id][2], annotation_dict[id][3], annotation_dict[id][4], annotation_dict[id][5]])
                                               
                        chr_x_sig_marker.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
                        chr_y_sig_marker.append(float(data[chr][i][6])/80*100 + 200+20)
                        break

        chr_y.append(float(data[chr][i][6])/80*100 + 200+20)
    #print(chr, positive_select_pos)


    num += 1# add the chromosome length scale
    #if chr == "Chr02":
    #    break
    if num % 2 == 1:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#30BBCE", alpha=1, rasterized=True)
    else:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#30BBCE", alpha=1, rasterized=True)
        

# add the significant line 10.04

#ax.axhline(y=229.641, color='b', linestyle='--', linewidth=0.5, alpha=0.5)
#ax.axhline(y=243.79, color='r', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [9.641/80*100+200+20, 9.641/80*100+200+20], color='b', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [23.79/80*100 + 200+20, 23.79/80*100 + 200+20], color='r', linestyle='--', linewidth=0.5, alpha=0.5)


# print the number of significant snp
# print(len(positive_select_pos))
box = plt.Rectangle((-3500000, 220), 1, 100, fc='black', ec='black', lw=1)
ax.add_patch(box)

#ax.axhline(y=14, color='black', linestyle='--', alpha=0.5)

# add the y axis label
#ax.text(-1800000, 10, "99%", ha='right', va='center', fontsize=10, color='blue')
#ax.text(-1800000, 30, "99.9%", ha='right', va='center', fontsize=10, color='red')
#ax.text(-11500000, 230, "-log10(µ)", ha='center', va='center', fontsize=7, color='black', rotation= 90)
#ax.text(0, -8, "chr", ha='center', va='center', fontsize=14, color='black')
ax.text(-3700000, 220, "0", ha='right', va='center', fontsize=7, color='black')
ax.text(-3700000, 220+30/80*100, "30", ha='right', va='center', fontsize=7)
#ax.text(-200000, 40, "40", ha='right', va='center', fontsize=10)
ax.text(-3700000, 220+60/80*100, "60", ha='right', va='center', fontsize=7)
#ax.text(-200000, 80, "80", ha='right', va='center', fontsize=10)
box = plt.Rectangle((216599865, 220), 4000000, 100, fc='#30BBCE', ec='none', lw=0.5, alpha=0.7)
ax.add_patch(box)
ax.text(216599865 + 2000000, 270, "EB", ha='center', va='center', fontsize=7, color='black', rotation= 90)

# write the gene annotation in to file
# sort the gene by the significant snp
EB_sig_gene = sorted(EB_sig_gene, key=lambda x: (x[0], int(x[1])))
with open("EB_sig_gene.xls", 'w') as f:
    f.write("chr\tpos\tgene_id\tPFAM_family\tPFAM\tGO\tKEGG\tNR\tSwissprot\tp_wald\n")
    for item in EB_sig_gene:
        #f.write("%s\n" % str(item).replace('[','').replace(']','').replace('\'',''))
        #f.write("%s\n" % '\t'.join(map(str, item)))
        f.write("%s\n" % '\t'.join(map(str, item)))
        
        



data ={}
with open("/home/chichi/data/china/china3/raisd/RAiSD_Report.EC_clean") as f:
    for line in f:
        line = line.strip()
        if line.startswith("// Ch"):
            chr = line.split( )[1]
            print(chr)
            data[line.split( )[1]] = []
            print(chr)
        elif line.startswith("// Contig"):
            break
        else:
            data[chr].append(line.split( ))
print(data.keys())
"""
# here need the 99% snp cutoff line
snp_values = []
for chr in data:
    for i in range(len(data[chr])):
        snp_values.append(float(data[chr][i][6]))
snp_values.sort()
print(snp_values[int(len(snp_values)*0.999)])
# 99.9% snp cutoff line 8.807
print(snp_values[int(len(snp_values)*0.995)])
# 99.5% snp cutoff line 3.72
print(snp_values[int(len(snp_values)*0.99)])
# 99% snp cutoff line 2.326
print(snp_values[int(len(snp_values)*0.95)])
# 95% snp cutoff line 0.3973    
"""
box = plt.Rectangle((0, 320), 216599865, 10, fc='white', ec='none', lw=0.5)
ax.add_patch(box)

box = plt.Rectangle((0, 327), 216599865, 1, fc='#ADD9ED', ec='#ADD9ED', lw=0.5)
ax.add_patch(box)


EC_sig_gene = []

for chr in data:
#for chr in {'Chr01'}:    

    chr_x = []
    chr_y = []
    chr_x_sig = []
    chr_y_sig = []
    chr_x_sig_marker = []
    chr_y_sig_marker = []

    gff_data_chr = gff_data[gff_data[0] == chr]
    print(gff_data_chr.head())
    for i in range(len(data[chr])):
        chr_x.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
        if float(data[chr][i][6]) > 8.807:
            chr_x_sig.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
            chr_y_sig.append(float(data[chr][i][6])/130*100 + 300+30)
            # to determine the significant snp in which gene or not
            if float(data[chr][i][6]) > 8.807:
                for index, row in gff_data_chr.iterrows():
                    if int(row[3]) <= int(data[chr][i][0]) <= int(row[4]):
                        #print(row[3], row[4], data[chr][i][0])
                        pos = int(data[chr][i][0])
                        id = gff_data_chr.loc[index, 8].split('Parent=')[1]
                        
                        if str(id) in annotation_dict:
                            print(annotation_dict[id])
                            EC_sig_gene.append([chr, pos, str(id), annotation_dict[id][0], annotation_dict[id][1], annotation_dict[id][2], annotation_dict[id][3], annotation_dict[id][4], annotation_dict[id][5]])
                                               
                        chr_x_sig_marker.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
                        chr_y_sig_marker.append(float(data[chr][i][6])/130*100 + 300+30)
                        break
        chr_y.append(float(data[chr][i][6])/130*100 + 300+30)
    #print(chr, positive_select_pos)


    num += 1# add the chromosome length scale
    #if chr == "Chr02":
    #    break
    if num % 2 == 1:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#ADD9ED", alpha=1, rasterized=True)
    else:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#ADD9ED", alpha=1, rasterized=True)

#ax.axhline(y=3.72, color='b', linestyle='--', linewidth=0.5, alpha=0.5)
#ax.axhline(y=8.807/130*100+300+30, color='r', linestyle='--', linewidth=0.5, alpha=0.5)
#ax.plot([-3500000, 216599865], [3.72, 3.72], color='b', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [8.807/130*100+300+30, 8.807/130*100+300+30], color='r', linestyle='--', linewidth=0.5, alpha=0.5)

# print the number of significant snp
# print(len(positive_select_pos))
box = plt.Rectangle((-3500000, 330), 1, 100, fc='black', ec='black', lw=1)
ax.add_patch(box)

# write the gene annotation in to file
# sort the gene by the significant snp
EC_sig_gene = sorted(EC_sig_gene, key=lambda x: (x[0], int(x[1])))
with open("EC_sig_gene.xls", 'w') as f: 
    f.write("chr\tpos\tgene_id\tPFAM_family\tPFAM\tGO\tKEGG\tNR\tSwissprot\tp_wald\n")
    for item in EC_sig_gene:
        #f.write("%s\n" % str(item).replace('[','').replace(']','').replace('\'',''))
        #f.write("%s\n" % '\t'.join(map(str, item)))
        f.write("%s\n" % '\t'.join(map(str, item)))

# add the y axis label
#ax.text(-1800000, 10.5, "99.5%", ha='right', va='center', fontsize=10, color='blue')
#ax.text(-1800000, 10, "99.9%", ha='right', va='center', fontsize=10, color='red')
#ax.text(-11500000, 80, "-log10(µ)", ha='center', va='center', fontsize=12, color='black', rotation= 90)
#ax.text(0, -13, "chr", ha='center', va='center', fontsize=14, color='black')
ax.text(-3700000, 330, "0", ha='right', va='center', fontsize=7, color='black')
#ax.text(-200000, 20, "20", ha='right', va='center', fontsize=10)
ax.text(-3700000, 50/130*100+300+30, "50", ha='right', va='center', fontsize=7)
#ax.text(-200000, 60, "60", ha='right', va='center', fontsize=10)
#ax.text(-200000, 80, "80", ha='right', va='center', fontsize=10)
ax.text(-3700000, 100/130*100+300+30, "100", ha='right', va='center', fontsize=7)
#ax.text(-200000, 120, "120", ha='right', va='center', fontsize=10)
#ax.text(-200000, 140, "140", ha='right', va='center', fontsize=10)

box = plt.Rectangle((216599865, 330), 4000000, 100, fc='#ADD9ED', ec='none', lw=0.5, alpha=0.7)
ax.add_patch(box)
ax.text(216599865 + 2000000, 380, "EC", ha='center', va='center', fontsize=7, color='black', rotation= 90)


data ={}
with open("/home/chichi/data/china/china3/raisd/RAiSD_Report.BA_clean") as f:
    for line in f:
        line = line.strip()

        if line.startswith("// Ch"):
            chr = line.split( )[1]
            print(chr)
            data[line.split( )[1]] = []
            print(chr)
        elif line.startswith("// Contig"):
            break
        else:
            data[chr].append(line.split( ))
print(data.keys())
"""
# here need the 99% snp cutoff line
snp_values = []
for chr in data:
    for i in range(len(data[chr])):
        snp_values.append(float(data[chr][i][6]))
snp_values.sort()
print(snp_values[int(len(snp_values)*0.999)])
# 99.9% snp cutoff line 26.47
print(snp_values[int(len(snp_values)*0.995)])
# 99.5% snp cutoff line 13.16
print(snp_values[int(len(snp_values)*0.99)])
# 99% snp cutoff line 8.225
print(snp_values[int(len(snp_values)*0.95)])
# 95% snp cutoff line 1.333
"""
box = plt.Rectangle((0, 430), 216599865, 10, fc='white', ec='none', lw=0.5)
ax.add_patch(box)
box = plt.Rectangle((0, 437), 216599865, 1, fc='black', ec='black', lw=0.5)
ax.add_patch(box)


BA_sig_gene = []

for chr in data:
#for chr in {'Chr15'}:
    chr_x = []
    chr_y = []
    chr_x_sig = []
    chr_y_sig = []
    chr_x_sig_marker = []
    chr_y_sig_marker = []

    gff_data_chr = gff_data[gff_data[0] == chr]
    print(gff_data_chr.head())
    for i in range(len(data[chr])):
        chr_x.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
        if float(data[chr][i][6]) > 13.16:
            chr_x_sig.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
            chr_y_sig.append(float(data[chr][i][6])/140*100 + 440)
            # to determine the significant snp in which gene or not
            if float(data[chr][i][6]) > 13.16:
                for index, row in gff_data_chr.iterrows():
                    if int(row[3]) <= int(data[chr][i][0]) <= int(row[4]):
                        #print(row[3], row[4], data[chr][i][0])
                        pos = int(data[chr][i][0])
                        id = gff_data_chr.loc[index, 8].split('Parent=')[1]

                        if str(id) in annotation_dict:
                            print(annotation_dict[id])
                            BA_sig_gene.append([chr, pos, str(id), annotation_dict[id][0], annotation_dict[id][1], annotation_dict[id][2], annotation_dict[id][3], annotation_dict[id][4], annotation_dict[id][5]])
                        
                        
                        chr_x_sig_marker.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
                        chr_y_sig_marker.append(float(data[chr][i][6])/140*100 + 440)
                        break
        chr_y.append(float(data[chr][i][6])/140*100 + 440)
    #print(chr, positive_select_pos)


    num += 1# add the chromosome length scale
    #if chr == "Chr02":
    #    break
    if num % 2 == 1:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "red", alpha=1, rasterized=True)
    else:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "red", alpha=1, rasterized=True)

#ax.axhline(y=13.16/140*100 + 440, color='b', linestyle='--', linewidth=0.5, alpha=0.5)
#ax.axhline(y=26.47/140*100 + 440, color='r', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [13.16/140*100 + 440, 13.16/140*100 + 440], color='b', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [26.47/140*100 + 440, 26.47/140*100 + 440], color='r', linestyle='--', linewidth=0.5, alpha=0.5)
# print the number of significant snp
# print(len(positive_select_pos))
box = plt.Rectangle((-3500000, 440), 1, 100, fc='black', ec='black', lw=1)
ax.add_patch(box)

# write the gene annotation in to file
# sort the gene by the significant snp
BA_sig_gene = sorted(BA_sig_gene, key=lambda x: (x[0], int(x[1])))
with open("BA_sig_gene.xls", 'w') as f:
    f.write("chr\tpos\tgene_id\tPFAM_family\tPFAM\tGO\tKEGG\tNR\tSwissprot\tp_wald\n")
    for item in BA_sig_gene:
        #f.write("%s\n" % str(item).replace('[','').replace(']','').replace('\'',''))
        #f.write("%s\n" % '\t'.join(map(str, item)))
        f.write("%s\n" % '\t'.join(map(str, item)))


#ax.axhline(y=14, color='black', linestyle='--', alpha=0.5)

# add the y axis label
#ax.text(-1800000, 10.5, "99.5%", ha='right', va='center', fontsize=10, color='blue')
#ax.text(-1800000, 30, "99.9%", ha='right', va='center', fontsize=10, color='red')
#ax.text(-11500000, 80, "-log10(µ)", ha='center', va='center', fontsize=14, color='black', rotation= 90)
#ax.text(0, -14, "chr", ha='center', va='center', fontsize=14, color='black')
ax.text(-3700000, 440, "0", ha='right', va='center', fontsize=7, color='black')
#ax.text(-200000, 20, "20", ha='right', va='center', fontsize=10)
#ax.text(-200000, 40, "40", ha='right', va='center', fontsize=10)
ax.text(-3700000, 60/140*100+440, "60", ha='right', va='center', fontsize=7)
#ax.text(-200000, 80, "80", ha='right', va='center', fontsize=10)
#ax.text(-200000, 100, "100", ha='right', va='center', fontsize=10)
ax.text(-3700000, 120/140*100+440, "120", ha='right', va='center', fontsize=7)
#ax.text(-200000, 140, "140", ha='right', va='center', fontsize=10)

box = plt.Rectangle((216599865, 440), 4000000, 100, fc='grey', ec='none', lw=0.5, alpha=0.7)
ax.add_patch(box)
ax.text(216599865 + 2000000, 490, "BA", ha='center', va='center', fontsize=7, color='black', rotation= 90)

data ={}
with open("/home/chichi/data/china/china3/raisd/RAiSD_Report.WA_clean") as f:
    for line in f:
        line = line.strip()

        if line.startswith("// Ch"):
            chr = line.split( )[1]
            print(chr)
            data[line.split( )[1]] = []
            print(chr)
        elif line.startswith("// Contig"):
            break
        else:
            data[chr].append(line.split( ))
print(data.keys())
"""
# here need the 99% snp cutoff line
snp_values = []
for chr in data:
    for i in range(len(data[chr])):
        snp_values.append(float(data[chr][i][6]))
snp_values.sort()
print(snp_values[int(len(snp_values)*0.999)])
# 99.9% snp cutoff line 79.59
print(snp_values[int(len(snp_values)*0.995)])
# 99.5% snp cutoff line 33.09
print(snp_values[int(len(snp_values)*0.99)])
# 99% snp cutoff line 20.77
print(snp_values[int(len(snp_values)*0.95)])
# 95% snp cutoff line 4.217
"""
box = plt.Rectangle((0, 540), 216599865, 10, fc='white', ec='none', lw=0.5)
ax.add_patch(box)
box = plt.Rectangle((0, 547), 216599865, 1, fc='#F3740B', ec='#F3740B', lw=0.5)
ax.add_patch(box)

WA_sig_gene = []

for chr in data:

#for chr in {'Chr15'}:
    chr_x = []
    chr_y = []
    chr_x_sig = []
    chr_y_sig = []
    chr_x_sig_marker = []
    chr_y_sig_marker = []
    
    gff_data_chr = gff_data[gff_data[0] == chr]
    print(gff_data_chr.head())
    for i in range(len(data[chr])):
        chr_x.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
        if float(data[chr][i][6]) > 33.09:
            chr_x_sig.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
            chr_y_sig.append(float(data[chr][i][6])/500*100 + 550)
            # to determine the significant snp in which gene or not
            if float(data[chr][i][6]) > 33.09:
                for index, row in gff_data_chr.iterrows():
                    if int(row[3]) <= int(data[chr][i][0]) <= int(row[4]):
                        #print(row[3], row[4], data[chr][i][0])
                        pos = int(data[chr][i][0])
                        id = gff_data_chr.loc[index, 8].split('Parent=')[1]
                        
                        if str(id) in annotation_dict:
                            print(annotation_dict[id])
                            WA_sig_gene.append([chr, pos, str(id), annotation_dict[id][0], annotation_dict[id][1], annotation_dict[id][2], annotation_dict[id][3], annotation_dict[id][4], annotation_dict[id][5]])
                        chr_x_sig_marker.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
                        chr_y_sig_marker.append(float(data[chr][i][6])/500*100 + 550)  
                        break
        chr_y.append(float(data[chr][i][6])/500*100 + 550)
    #print(chr, positive_select_pos)


    num += 1# add the chromosome length scale
    if num % 2 == 1:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#F3740B", alpha=1, rasterized=True)
    else:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#F3740B", alpha=1, rasterized=True)
    #if chr == "Chr02":
    #    break

#ax.axhline(y=33.09/500*100 + 550, color='b', linestyle='--', linewidth=0.5, alpha=0.5)
#ax.axhline(y=79.59/500*100 + 550, color='r', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [33.09/500*100 + 550, 33.09/500*100 + 550], color='b', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [79.59/500*100 + 550, 79.59/500*100 + 550], color='r', linestyle='--', linewidth=0.5, alpha=0.5)
# print the number of significant snp
# print(len(positive_select_pos))
box = plt.Rectangle((-3500000, 550), 1, 100, fc='black', ec='black', lw=1)
ax.add_patch(box)

# write the gene annotation in to file
# sort the gene by the significant snp
WA_sig_gene = sorted(WA_sig_gene, key=lambda x: (x[0], int(x[1])))
with open("WA_sig_gene.xls", 'w') as f:
    f.write("chr\tpos\tgene_id\tPFAM_family\tPFAM\tGO\tKEGG\tNR\tSwissprot\tp_wald\n")
    for item in WA_sig_gene:
        #f.write("%s\n" % str(item).replace('[','').replace(']','').replace('\'',''))
        #f.write("%s\n" % '\t'.join(map(str, item)))
        f.write("%s\n" % '\t'.join(map(str, item)))


#ax.axhline(y=14, color='black', linestyle='--', alpha=0.5)

# add the y axis label
#ax.text(-1800000, 50, "99.5%", ha='right', va='center', fontsize=10, color='blue')
#ax.text(-1800000, 150, "99.9%", ha='right', va='center', fontsize=10, color='red')
#ax.text(-11500000, 300, "-log10(µ)", ha='center', va='center', fontsize=12, color='black', rotation= 90)
#ax.text(0, -50, "chr", ha='center', va='center', fontsize=14, color='black')
ax.text(-3700000, 0, "0", ha='right', va='center', fontsize=7, color='black')
#ax.text(-200000, 100, "100", ha='right', va='center', fontsize=10)
ax.text(-3700000, 200/500*100 + 550, "200", ha='right', va='center', fontsize=7)
#ax.text(-200000, 300, "300", ha='right', va='center', fontsize=10)
ax.text(-3700000, 400/500*100 + 550, "400", ha='right', va='center', fontsize=7)
#ax.text(-200000, 500, "500", ha='right', va='center', fontsize=10)

box = plt.Rectangle((216599865, 550), 4000000, 100, fc='#F3740B', ec='none', lw=0.5, alpha=0.7)
ax.add_patch(box)
ax.text(216599865 + 2000000, 600, "WA", ha='center', va='center', fontsize=7, color='black', rotation= 90)

    


data ={}
with open("/home/chichi/data/china/china3/raisd/RAiSD_Report.WB_clean")  as f:
    for line in f:
        line = line.strip()

        if line.startswith("// Ch"):
            chr = line.split( )[1]
            print(chr)
            data[line.split( )[1]] = []
            print(chr)
        elif line.startswith("// Contig"):
            break
        else:
            data[chr].append(line.split( ))
print(data.keys())
"""
# here need the 99% snp cutoff line
snp_values = []
for chr in data:
    for i in range(len(data[chr])):
        snp_values.append(float(data[chr][i][6]))
snp_values.sort()
print(snp_values[int(len(snp_values)*0.999)])
# 99.9% snp cutoff line 168.7
print(snp_values[int(len(snp_values)*0.995)])
# 99.5% snp cutoff line 61.96
print(snp_values[int(len(snp_values)*0.99)])
# 99% snp cutoff line 36.44
print(snp_values[int(len(snp_values)*0.95)])
# 95% snp cutoff line 6.869
"""
box = plt.Rectangle((0, 650), 216599865, 10, fc='white', ec='none', lw=0.5)
ax.add_patch(box)
box = plt.Rectangle((0, 657), 216599865, 1, fc='#EFA38A', ec='#EFA38A', lw=0.5)
ax.add_patch(box)


WB_sig_gene = []


for chr in data:
#for chr in {'Chr15'}:
    chr_x = []
    chr_y = []
    chr_x_sig = []
    chr_y_sig = []
    chr_x_sig_marker = []
    chr_y_sig_marker = []
    
    gff_data_chr = gff_data[gff_data[0] == chr]
    print(gff_data_chr.head())
    for i in range(len(data[chr])):
        chr_x.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
        if float(data[chr][i][6]) > 61.96:
            chr_x_sig.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
            chr_y_sig.append(float(data[chr][i][6])/600*100 + 660)
            # to determine the significant snp in which gene or not
            if float(data[chr][i][6]) > 61.96:
                for index, row in gff_data_chr.iterrows():
                    if int(row[3]) <= int(data[chr][i][0]) <= int(row[4]):
                        #print(row[3], row[4], data[chr][i][0])
                        pos = int(data[chr][i][0])
                        id = gff_data_chr.loc[index, 8].split('Parent=')[1]
                        
                        if str(id) in annotation_dict:
                            print(annotation_dict[id])
                            WB_sig_gene.append([chr, pos, str(id), annotation_dict[id][0], annotation_dict[id][1], annotation_dict[id][2], annotation_dict[id][3], annotation_dict[id][4], annotation_dict[id][5]])
                                        
                        chr_x_sig_marker.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
                        chr_y_sig_marker.append(float(data[chr][i][6])/600*100 + 660)
                        break
        chr_y.append(float(data[chr][i][6])/600*100 + 660)  

    num += 1# add the chromosome length scale
    #if chr == "Chr02":
    #    break
    if num % 2 == 1:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#EFA38A", alpha=1, rasterized=True)
    else:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#EFA38A", alpha=1, rasterized=True)
        
#ax.axhline(y=61.96/600*100 + 660, color='b', linestyle='--', linewidth=0.5, alpha=0.5)
#ax.axhline(y=168.7/600*100 + 660, color='r', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [61.96/600*100 + 660, 61.96/600*100 + 660], color='b', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [168.7/600*100 + 660, 168.7/600*100 + 660], color='r', linestyle='--', linewidth=0.5, alpha=0.5)

# print the number of significant snp
# print(len(positive_select_pos))

# write the gene annotation in to file
# sort the gene by the significant snp
WB_sig_gene = sorted(WB_sig_gene, key=lambda x: (x[0], int(x[1])))
with open("WB_sig_gene.xls", 'w') as f:
    f.write("chr\tpos\tgene_id\tPFAM_family\tPFAM\tGO\tKEGG\tNR\tSwissprot\tp_wald\n")
    for item in WB_sig_gene:
        #f.write("%s\n" % str(item).replace('[','').replace(']','').replace('\'',''))
        #f.write("%s\n" % '\t'.join(map(str, item)))
        f.write("%s\n" % '\t'.join(map(str, item)))
        

box = plt.Rectangle((-3500000, 660), 1, 100, fc='black', ec='black', lw=1)
ax.add_patch(box)
# add the y axis label
#ax.text(-1800000, 50, "99.5%", ha='right', va='center', fontsize=10, color='blue')
#ax.text(-1800000, 150, "99.9%", ha='right', va='center', fontsize=10, color='red')
#ax.text(-11500000, 400, "-log10(µ)", ha='center', va='center', fontsize=12, color='black', rotation= 90)
#ax.text(0, -60, "chr", ha='center', va='center', fontsize=14, color='black')
ax.text(-3700000, 0/600*100 + 660, "0", ha='right', va='center', fontsize=7, color='black')
#ax.text(-200000, 100, "100", ha='right', va='center', fontsize=10)
ax.text(-3700000, 250/600*100 + 660, "250", ha='right', va='center', fontsize=7)
#ax.text(-200000, 300, "300", ha='right', va='center', fontsize=10)
#ax.text(-200000, 400, "400", ha='right', va='center', fontsize=10)
ax.text(-3700000, 500/600*100 + 660, "500", ha='right', va='center', fontsize=7)

#ax.text(-200000, 600, "600", ha='right', va='center', fontsize=10)
box = plt.Rectangle((216599865, 660), 4000000, 100, fc='#EFA38A', ec='none', lw=0.5, alpha=0.7)
ax.add_patch(box)
ax.text(216599865 + 2000000, 710, "WB", ha='center', va='center', fontsize=7, color='black', rotation= 90)


  
data ={}
with open("/home/chichi/data/china/china3/raisd/RAiSD_Report.WC_clean") as f:
    for line in f:
        line = line.strip()
        if line.startswith("// Ch"):
            chr = line.split( )[1]
            print(chr)
            data[line.split( )[1]] = []
            print(chr)
        elif line.startswith("// Contig"):
            break
        else:
            data[chr].append(line.split( ))
print(data.keys())
"""
# here need the 99% snp cutoff line
snp_values = []
for chr in data:
    for i in range(len(data[chr])):
        snp_values.append(float(data[chr][i][6]))
snp_values.sort()
print(snp_values[int(len(snp_values)*0.999)])
# 99.9% snp cutoff line 56.98
print(snp_values[int(len(snp_values)*0.995)])
# 99.5% snp cutoff line 21.37
print(snp_values[int(len(snp_values)*0.99)])
# 99% snp cutoff line 14.52
print(snp_values[int(len(snp_values)*0.95)])
# 95% snp cutoff line 2.735
"""
box = plt.Rectangle((0, 760), 216599865, 10, fc='white', ec='none', lw=0.5)
ax.add_patch(box)
box = plt.Rectangle((0, 767), 216599865, 1, fc='#D2352C', ec='#D2352C', lw=0.5)
ax.add_patch(box)

WC_sig_gene = []

for chr in data:
#for chr in {'Chr15'}:
    chr_x = []
    chr_y = []
    chr_x_sig = []
    chr_y_sig = []
    chr_x_sig_marker = []
    chr_y_sig_marker = []
    gff_data_chr = gff_data[gff_data[0] == chr]
    print(gff_data_chr.head())
    for i in range(len(data[chr])):
        chr_x.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
        if float(data[chr][i][6]) > 21.37:
            chr_x_sig.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
            chr_y_sig.append(float(data[chr][i][6])/420*100 + 770)
            # to determine the significant snp in which gene or not
            if float(data[chr][i][6]) > 21.37:
                for index, row in gff_data_chr.iterrows():
                    if int(row[3]) <= int(data[chr][i][0]) <= int(row[4]):
                        #print(row[3], row[4], data[chr][i][0])
                        pos = int(data[chr][i][0])
                        id = gff_data_chr.loc[index, 8].split('Parent=')[1]
                        
                        if str(id) in annotation_dict:
                            print(annotation_dict[id])
                            WC_sig_gene.append([chr, pos, str(id), annotation_dict[id][0], annotation_dict[id][1], annotation_dict[id][2], annotation_dict[id][3], annotation_dict[id][4], annotation_dict[id][5]])
                                            
                        chr_x_sig_marker.append(int(data[chr][i][0])+int(chr_start_pos[chr]))
                        chr_y_sig_marker.append(float(data[chr][i][6])/420*100 + 770)
                        break
        chr_y.append(float(data[chr][i][6])/420*100 + 770)

    num += 1# add the chromosome length scale
    #if chr == "Chr02":
    #    break
    if num % 2 == 1:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='grey', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#D2352C", alpha=1, rasterized=True)
    else:
        ax.scatter(chr_x, chr_y, s=1, label=chr, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig, chr_y_sig, s=1, facecolors='black', edgecolors='none', alpha=0.5, rasterized=True)
        ax.scatter(chr_x_sig_marker, chr_y_sig_marker, s=1.5, color = "#D2352C", alpha=1, rasterized=True)

#ax.axhline(y=21.37/420*100 + 770, color='b', linestyle='--', linewidth=0.5, alpha=0.5)
#ax.axhline(y=56.98/420*100 + 770, color='r', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [21.37/420*100 + 770, 21.37/420*100 + 770], color='b', linestyle='--', linewidth=0.5, alpha=0.5)
ax.plot([-3500000, 216599865], [56.98/420*100 + 770, 56.98/420*100 + 770], color='r', linestyle='--', linewidth=0.5, alpha=0.5)

# write the gene annotation in to file
# sort the gene by the significant snp
WC_sig_gene = sorted(WC_sig_gene, key=lambda x: (x[0], int(x[1])))
with open("WC_sig_gene.xls", 'w') as f:
    f.write("chr\tpos\tgene_id\tPFAM_family\tPFAM\tGO\tKEGG\tNR\tSwissprot\tp_wald\n")
    for item in WC_sig_gene:
        #f.write("%s\n" % str(item).replace('[','').replace(']','').replace('\'',''))
        #f.write("%s\n" % '\t'.join(map(str, item)))
        f.write("%s\n" % '\t'.join(map(str, item)))
        

box = plt.Rectangle((-3500000, 767), 1, 100, fc='black', ec='black', lw=1)
ax.add_patch(box)

# add the y axis label
#ax.text(-1800000, 50, "99.5%", ha='right', va='center', fontsize=10, color='blue')
#ax.text(-1800000, 150, "99.9%", ha='right', va='center', fontsize=10, color='red')
#ax.text(-11500000, 300, "-log10(µ)", ha='center', va='center', fontsize=12, color='black', rotation= 90)
#ax.text(0, -42, "chr", ha='center', va='center', fontsize=14, color='black')
ax.text(-3700000, 0/420*100 + 770, "0", ha='right', va='center', fontsize=7, color='black')
ax.text(-3700000, 200/420*100 + 770, "200", ha='right', va='center', fontsize=7)
ax.text(-3700000, 400/420*100 + 770, "400", ha='right', va='center', fontsize=7)

box = plt.Rectangle((216599865, 770), 4000000, 100, fc='#D2352C', ec='none', lw=0.5, alpha=0.7)
ax.add_patch(box)
ax.text(216599865 + 2000000, 820, "WC", ha='center', va='center', fontsize=7, color='black', rotation= 90)

       




plt.savefig("raisd_all2.pdf", dpi=1200, bbox_inches='tight')
