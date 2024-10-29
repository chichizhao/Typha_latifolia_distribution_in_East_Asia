#!/bin/python
# function: plot the manhanttan plot of the gwas data and annotate the genes
# usage: python plot_sig_and_annotegene.py
#

#!/bin/python
# function: plot the manhanttan plot of the gwas data
# usage: python plot_snp.py

# step 1: read the gwas data with p value
# file: sig_SNPs_bio_5_mod_sub_env_data4.part5_merge3_filter_variants.csv
# read the data with pandas
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

chr_len ={'Chr01': 24370473,'Chr02': 18857030,'Chr03': 18714887,'Chr04': 18467694,'Chr05': 17512124,'Chr06': 14553481,'Chr07': 14503471,'Chr08': 13751786,'Chr09': 11865621,'Chr10': 11839088,'Chr11': 11614203,'Chr12': 10714389,'Chr13': 10242974,'Chr14': 9831573,'Chr15': 9461071,'Contig1': 452046,'Contig2': 203240}
chr_start_pos= {'Chr01': 0,'Chr02': 24370473,'Chr03': 43227503,'Chr04': 61942390,'Chr05': 80410084,'Chr06': 97922208,'Chr07': 112775689,'Chr08': 127279160,'Chr09': 141030946,'Chr10': 152896567,'Chr11': 164735655,'Chr12': 176349858,'Chr13': 187064247,'Chr14': 197307221,'Chr15': 207138794,'Contig1': 216599865,'Contig2': 261846911}


file = '/home/chichi/data/china/china3/gwas/bio_all/Linear_Mixed_Model/bio_8/bio_8_20240108_161859/bio_8_mod_sub_env_data4.part2_merge3_filter_variants.assoc2.txt'
data = pd.read_csv(file, sep='\t',index_col= False)
print(data.head())

# step 2: plot the manhattan plot
fig = plt.figure(figsize=(14, 2.5))
ax = fig.add_subplot(111)
# off the axis
ax.axis('off')
ax.set_xlim(-4000000, 223599865)
ax.set_ylim(-0.2, 20)
num = 0
#print(len(data))
# remove the data with p_wald < 1e-4
# Create a new column 'chr_str' with formatted chromosome numbers
data['chr_str'] = 'Chr' + data['chr'].apply(lambda x: f'{x:02}')
print('new column chr_str created')
# Create a new column 'color' based on the 'p_wald' values
data['color'] = np.where(data['p_wald'] >= 7.46e-08, 'grey', 'blue')

# Create a new column 'x' for the x-coordinates of the scatter plot
data['x'] = data.apply(lambda row: chr_start_pos[row['chr_str']] + row['ps'], axis=1)

print('new column x created')
# Create a new column 'y' for the y-coordinates of the scatter plot
data['y'] = -np.log10(data['p_wald'])
print('new column y created')

# Plot all rows
#for _, row in data.iterrows():
#    ax.scatter(row['x'], row['y'], color=row['color'], s=2,alpha=0.5,rasterized=True)

# plot the data with with the chromosomes 
for i in range(1, 16):
    if i < 10:
        chr = 'Chr0' + str(i)
    else:
        chr = 'Chr' + str(i)
    chr_data_x = data[data['chr_str'] == chr]['x']
    chr_data_y = data[data['chr_str'] == chr]['y']
    chr_data_color = data[data['chr_str'] == chr]['color']
    if i%2 == 0:
        # replace the grey color with black
        chr_data_color = chr_data_color.replace('grey', 'black')
    ax.scatter(chr_data_x, chr_data_y, color=chr_data_color, s=2,alpha=0.5,rasterized=True)
#for _, row in data.iloc[:5000].iterrows():
#    ax.scatter(row['x'], row['y'], color=row['color'], s=1)
print('plot the manhattan plot')
# annotate the genes 
## read the gene annotation file
gff ="/home/chichi/data/china/china3/gwas/pasa2.longest.filter.gff3.cds"
gff_data = pd.read_csv(gff, sep='\t',index_col= False)
#print(gff_data.head())

## read the gene annotation file 
annotation = "/home/chichi/data/china/china3/gwas/func_merg.xls"
annotation_data = pd.read_csv(annotation, sep='\t',index_col= False)
annotation_dict = {}
for i in range(len(annotation_data)):
    if str(annotation_data.loc[i, 'PFAM']).startswith('PF'):
        annotation_dict[annotation_data.loc[i, '#ID']] = [annotation_data.loc[i, 'PFAM'].split('(')[0], annotation_data.loc[i, 'PFAM'], annotation_data.loc[i, 'GO'], annotation_data.loc[i, 'KEGG'], annotation_data.loc[i, 'NR'], annotation_data.loc[i, 'Swissprot']]

file = '/home/chichi/data/china/china3/gwas/bio_all/Linear_Mixed_Model/bio_8/bio_8_20240108_161859/best_p-values/p_wald_bio_8_mod_sub_env_data4.part2_merge3_filter_variants_top0.001.csv'
data = pd.read_csv(file, sep=',',index_col= False)

print('gene annotation file read')
gene_list_info = []


# here we need check the snp in CDS region or not
for i in range(len(data)):
    chr = data.loc[i, 'chr']
    pos = data.loc[i, 'ps']
    if chr < 10:
        chr = 'Chr0' + str(chr)
    else:
        chr = 'Chr' + str(chr)
    #print(chr)
    # check the snp in the CDS region or not
    gff_chr = gff_data[gff_data['Contig1'] == chr]
    #print(gff_chr)
    if -np.log10(data.loc[i, 'p_wald']) >= 7.46:
        for j in range(len(gff_chr)):
            if pos >= gff_chr.iloc[j, 3] and pos <= gff_chr.iloc[j, 4]:
                #ax.text(chr_start_pos[chr] + pos, -np.log10(data.loc[i, 'p_wald']), gff_chr.iloc[j, 8], ha='center', va='center', fontsize=8, color='black')
                id = gff_chr.iloc[j, 8].split('Parent=')[1]
                if id in annotation_dict:
                    # ax.text(chr_start_pos[chr] + pos, -np.log10(data.loc[i, 'p_wald'])+1, annotation_dict[id], ha='center', va='center', fontsize=12, color='black')
                    # add the gene information to the gene_list_info
                    gene_list_info.append([chr, pos, str(id), annotation_dict[id][0], annotation_dict[id][1], annotation_dict[id][2], annotation_dict[id][3], annotation_dict[id][4], annotation_dict[id][5], -np.log10(data.loc[i, 'p_wald'])])
                    # highlight the gene point
                    ax.scatter(chr_start_pos[chr] + pos, -np.log10(data.loc[i, 'p_wald']), color='red', s=5,rasterized=True)
                break
#sort the gene_list_info based on the chromosome and position
gene_list_info = sorted(gene_list_info, key=lambda x: (x[0], x[1]))            
with open('bio8_top0.01_annotated_genes.xls', 'w') as f:  
    #f.write("chr\tpos\tgene_id\tPFAM_family\tPFAM\tGO\tKEGG\tNR\tSwissprot\n")
    #f.write("chr,pos,gene_id,PFAM_family,PFAM,GO,KEGG,NR,Swissprot\n")
    f.write("chr\tpos\tgene_id\tPFAM_family\tPFAM\tGO\tKEGG\tNR\tSwissprot\tp_wald\n")
    for item in gene_list_info:
        item = item[:-1]
        #f.write("%s\n" % str(item).replace('[','').replace(']','').replace('\'',''))
        #f.write("%s\n" % '\t'.join(map(str, item)))
        f.write("%s\n" % '\t'.join(map(str, item)))
print( gene_list_info)

# remove the snp on the same gene
new_gene_list_info = []
for i in range(len(gene_list_info)):
    if i == 0:
        new_gene_list_info.append(gene_list_info[i])
    else:
        if gene_list_info[i][2] != gene_list_info[i-1][2]:
            new_gene_list_info.append(gene_list_info[i])

# add the gene code in the manhattan plot
for i in range(len(new_gene_list_info)):
    add_text = new_gene_list_info[i][2].replace('evm.model.','').replace('.','-')
    ax.text(chr_start_pos[new_gene_list_info[i][0]] + new_gene_list_info[i][1], new_gene_list_info[i][9]+1, add_text, ha='center', va='center', fontsize=12, color='black')

# add the line at the significant threshold of 7.46e-08
#ax.axhline(y=-np.log10(7.46e-08), color='b', linestyle='--', alpha=0.5)
plt.plot([0, 216599865], [-np.log10(7.46e-08), -np.log10(7.46e-08)], color='blue', linestyle='--', alpha=0.5)
# add the mark of the 0
#ax.axhline(y=0, color='black', linestyle='-')
plt.plot([-4000000, 223599865], [-0.1, -0.1], color='black', linestyle='-')
box = plt.Rectangle((-4000000, 0), 100000, 20, fc='none', ec='black', lw=1)
ax.add_patch(box)


# add the mark of the chromosome
for i in range(1, 16):
    if i < 10:
        chr = 'Chr0' + str(i)
    else:
        chr = 'Chr' + str(i)
    ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -1.4, str(i), ha='center', va='center', fontsize=14, color='black')
    # add the chromosome region line 
    #if i%2 == 1:
    #    box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 20, fc='grey', ec='none', lw=1, alpha=0.5)
    #    ax.add_patch(box)
#ax.axvline(x=chr_start_pos['Contig1'], color='black', linestyle='-')
#ax.axvline(x=chr_start_pos["Chr01"], color='black', linestyle='-')
# add the y axis label
#ax.text(-0.1, 7.46, '-log10(7.46e-08)', ha='center', va='center', fontsize=12, color='black')
ax.text(-13000000, 10, '-log10(p_wald)', ha='center', va='center', fontsize=14, color='black', rotation=90) 
ax.text(-4700000, 0, '0', ha='right', va='center', fontsize=14, color='black')
ax.text(-4700000, 5, '5', ha='right', va='center', fontsize=14, color='black')
#ax.text(-3000000, 7.46, '7.46', ha='center', va='center', fontsize=12, color='blue')
ax.text(-4700000, 10, '10', ha='right', va='center', fontsize=14, color='black')
ax.text(-4700000, 15, '15', ha='right', va='center', fontsize=14, color='black')
ax.text(-4700000, 20, '20', ha='right', va='center', fontsize=14, color='black')
ax.text(0, -1.4, 'Chr', ha='center', va='center', fontsize=14, color='black')
#box = plt.Rectangle((0, 0), 216599865, 20, fc='none', ec='black', lw=3)
#ax.add_patch(box)
box = plt.Rectangle((218599865, 0), 5000000, 20, fc='grey', ec='black', lw=1, alpha=0.1)
ax.add_patch(box)
# add the bio2 label
ax.text(216599865+4500000, 17, 'bio8', ha='center', va='center', fontsize=14, color='black',rotation=90)
# add the x-axis and y-axis label
# ax.set_xlabel('Position')
# ax.set_ylabel('-log10(p-value)')
# ax.set_title('Manhattan plot of GWAS data')
plt.savefig('manhattan_bio8.pdf', dpi=300, bbox_inches='tight')