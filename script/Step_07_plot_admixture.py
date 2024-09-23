#!/bin/python
# -*- coding: utf-8 -*-
# function: plot admixture results of k=4, k=5, k=6 in name
# Author: CHICHI
# Usage: python plot_admixture.py
#### basic description ####


import argparse
import matplotlib.pyplot as plt
rename_dict = {'TKS':'WC3', 'Xu4273':'WC2', 'Xu4319':'WC1', 'HTB':'WB2', 'Xu4436':'WB1','Xu4243':'WA3', 'Xu4267':'WA2', 'Xu4217':'WA1','T37':'BA1','T28':'BA2','T34':'BA3','Xu4624':'EC5','Xu7165':'EC4','Xu269':'EC3','Xu513':'EC2','291':'EC1','Xu254':'EB6','Xu627':'EB5','Xu552':'EB4','Xu407':'EB3','QTH':'EB2','Xu320':'EB1','TM':'EA2','Xu3748':'EA1', 'Xu6531':'MA3','Xu6547':'MA2','jzg':'MA1'}


def get_sample_k6():
    sample = []
    with open("popmap2.txt", 'r') as f:
        for line in f:
            sample.append(line.strip().split()[0])
    return sample

def get_q_k6():
    q = []
    with open("merge3_filter_variants.6_7.Q", 'r') as f:
        for line in f:
            q.append(line.strip().split())
    return q

def get_order_k6(q,sample):
    # here we use the strategy the that based the groups, than the order of the first column, second column and so on
    samples_num = len(q)
    data = {}
    site = ['TKS', 'Xu4273', 'Xu4319', 'HTB', 'Xu4436', 'Xu4243', 'Xu4267', 'Xu4217','T37','T28','T34','Xu4624','Xu7165','Xu269','Xu513','291','QTH','Xu254','Xu320','Xu627','Xu552','Xu407','TM','Xu3748', 'Xu6531','Xu6547','jzg']
    for i in site:
        data[i] = []
    for i in range(samples_num):
        for j in range(len(site)):
            if sample[i].startswith(site[j]):
                data[site[j]].append([sample[i], i , q[i]])
    #print(data)
    #print(data)
    # add we need sort the data in each group
    

    for i in data:
        if i in ['TKS', 'Xu4273', 'Xu4319', 'HTB', 'Xu4436', 'Xu4243', 'Xu4267', 'Xu4217']:
            #sprint(data[i][1][2][1])
            data[i] = sorted(data[i], key=lambda x: x[2][1], reverse= False)
        if i in ['T37']:
            data[i] = sorted(data[i], key=lambda x: x[2][4], reverse= True)
            
        if i in ['291','Xu513']:
            data[i] = sorted(data[i], key=lambda x: x[2][5], reverse=True)
        if i in ['Xu407']:
            data[i] = sorted(data[i], key=lambda x: x[2][0], reverse=True)
        if i in ['QTH','Xu320','TM']:
            data[i] = sorted(data[i], key=lambda x: x[2][0], reverse=True)
        if i in ['Xu3748']:
            data[i] = sorted(data[i], key=lambda x: x[2][0], reverse=False)
    print(data)
    
    # then we need to get the order of the samples
    return data


def get_sample_k5():
    sample = []
    with open("popmap2.txt", 'r') as f:
        for line in f:
            sample.append(line.strip().split()[0])
    return sample

def get_q_k5():
    q = []
    with open("merge3_filter_variants.5_3.Q", 'r') as f:
        for line in f:
            q.append(line.strip().split())
    return q

def get_order_k5(q,sample):
    # here we use the strategy the that based the groups, than the order of the first column, second column and so on
    samples_num = len(q)
    data = {}
    site = ['TKS', 'Xu4273', 'Xu4319', 'HTB', 'Xu4436','Xu4243', 'Xu4267', 'Xu4217','T37','T28','T34','Xu4624','Xu7165','Xu269','Xu513','291','Xu254','Xu627','Xu552','Xu407','QTH','Xu320','TM','Xu3748', 'Xu6531','Xu6547','jzg']
    for i in site:
        data[i] = []
    for i in range(samples_num):
        for j in range(len(site)):
            if sample[i].startswith(site[j]):
                data[site[j]].append([sample[i], i , q[i]])
    #print(data)
    #print(data)
    # add we need sort the data in each group
    

    for i in data:
        if i in ['TKS', 'Xu4273', 'Xu4319', 'HTB', 'Xu4436', 'Xu4243', 'Xu4267', 'Xu4217']:
            #sprint(data[i][1][2][1])
            data[i] = sorted(data[i], key=lambda x: x[2][4], reverse= False)
        if i in ['T37']:
            data[i] = sorted(data[i], key=lambda x: x[2][4], reverse= True)
            
        if i in ['291','Xu513']:
            data[i] = sorted(data[i], key=lambda x: x[2][2], reverse=True)
        if i in ['Xu407']:
            data[i] = sorted(data[i], key=lambda x: x[2][0], reverse=True)
        if i in ['QTH','Xu320','TM']:
            data[i] = sorted(data[i], key=lambda x: x[2][0], reverse=True)
        if i in ['Xu3748']:
            data[i] = sorted(data[i], key=lambda x: x[2][0], reverse=True)
    print(data)
    # then we need to get the order of the samples
    return data

def get_sample_k4():
    sample = []
    with open("popmap2.txt", 'r') as f:
        for line in f:
            sample.append(line.strip().split()[0])
    return sample

def get_q_k4():
    q = []
    with open("merge3_filter_variants.4_7.Q", 'r') as f:
        for line in f:
            q.append(line.strip().split())
    return q

def get_order_k4(q,sample):
    # here we use the strategy the that based the groups, than the order of the first column, second column and so on
    samples_num = len(q)
    data = {}
    site = ['TKS', 'Xu4273', 'Xu4319', 'HTB', 'Xu4436','Xu4243', 'Xu4267', 'Xu4217','T37','T28','T34','Xu4624','Xu7165','Xu269','Xu513','291','Xu254','Xu627','Xu552','Xu407','QTH','Xu320','TM','Xu3748', 'Xu6531','Xu6547','jzg']
    for i in site:
        data[i] = []
    for i in range(samples_num):
        for j in range(len(site)):
            if sample[i].startswith(site[j]):
                data[site[j]].append([sample[i], i , q[i]])
    #print(data)
    #print(data)
    # add we need sort the data in each group
    
    for i in data:
        if i in ['TKS', 'Xu4273', 'Xu4319', 'HTB', 'Xu4436', 'Xu4243', 'Xu4267', 'Xu4217']:
            #sprint(data[i][1][2][1])
            data[i] = sorted(data[i], key=lambda x: x[2][0], reverse= False)
        if i in ['T37']:
            data[i] = sorted(data[i], key=lambda x: x[2][0], reverse= False)
            
        if i in ['291','Xu513']:
            data[i] = sorted(data[i], key=lambda x: x[2][0], reverse=True)
        if i in ['Xu407']:
            data[i] = sorted(data[i], key=lambda x: x[2][0], reverse=True)
        if i in ['QTH','Xu320','TM']:
            data[i] = sorted(data[i], key=lambda x: x[2][2], reverse=False)
        if i in ['Xu3748']:
            data[i] = sorted(data[i], key=lambda x: x[2][0], reverse=True)
    print(data)
    # then we need to get the order of the samples
    return data


def plot_admixture_k4_k5_k6(data_k4, data_k5, data_k6):
    fig, ax = plt.subplots(figsize=(45, 8))
    plt.axis('off')

    region = ["West China", "Baikal", "East China", "Middle China"]
    color_set = ['#EB100C','#73786A','#4BC2C9','#6AA84F']
    region_dict = {"West China":["TKS", "Xu4273", "Xu4319", "HTB", "Xu4436","Xu4243", "Xu4267", "Xu4217"], "Baikal":["T37", "T28", "T34"], "East China":["Xu4624", "Xu7165", "Xu269", "Xu513", "291","Xu254", "Xu627", "Xu552", "Xu407", "QTH", "Xu320", "TM", "Xu3748"], "Middle China":["Xu6531", "Xu6547", "jzg"]}
    start = 0
    postion = 0
    postion2 = 0
    for i in range(len(region)):
        for j in range(len(region_dict[region[i]])):
            for k in range(len(data_k6[region_dict[region[i]][j]])):
                postion += 1
            postion += 0.5
        box = plt.Rectangle((start+postion2, 0), postion-postion2-0.5, 4.7, fc=color_set[i], alpha=0.1)
        ax.add_patch(box)
        box = plt.Rectangle((start+postion2-0.05, 0-0.01), postion-postion2-0.5+0.1, 4.72, fc="none", ec=color_set[i], lw=1)
        ax.add_patch(box)
        box = plt.Rectangle((start+postion2-0.05, 4.6), postion-postion2-0.5+0.1, 0.1, fc=color_set[i], alpha=1)
        ax.add_patch(box)
        ax.text(start+postion2+(postion-postion2)/2, 4.8, region[i], ha='center', va='center', fontsize=15, color='black')
        postion2 = postion    

    # we need add the group name
    # add the color bar of the box plot
    group_name = ['WC', 'WB', 'WA', 'BA', 'EC', 'EB', 'EA', 'MA']
    color_set = ['#D2352C','#EFA38A','#F3740B','#808080','#ADD9ED','#30BBCE','#6BA9EC','#68AC56']
    group_name_dict = {"WC":["TKS", "Xu4273", "Xu4319"], "WB":["HTB", "Xu4436"], "WA":["Xu4243", "Xu4267", "Xu4217"], "BA":["T37", "T28", "T34"], "EC":["Xu4624", "Xu7165", "Xu269", "Xu513", "291"], "EB":["Xu254", "Xu627", "Xu552", "Xu407", "QTH", "Xu320"], "EA":["TM", "Xu3748"], "MA":["Xu6531", "Xu6547", "jzg"]}
    start = 0
    postion = 0
    postion2 = 0
    for i in range(len(group_name)):
        for j in range(len(group_name_dict[group_name[i]])):
            for k in range(len(data_k6[group_name_dict[group_name[i]][j]])):
                postion += 1
            postion += 0.5
        box = plt.Rectangle((start+postion2, 0), postion-postion2-0.5, 4.35, fc=color_set[i], alpha=0.2)
        ax.add_patch(box)
        box = plt.Rectangle((start+postion2-0.05, 0-0.01), postion-postion2-0.5+0.1, 4.35, fc="none", ec=color_set[i], lw=0.5)
        ax.add_patch(box)
        box = plt.Rectangle((start+postion2-0.05, 4.3), postion-postion2-0.5+0.1, 0.05, fc=color_set[i], alpha=1)
        ax.add_patch(box)
        ax.text(start+postion2+(postion-postion2)/2, 4.47, group_name[i], ha='center', va='center', fontsize=15, color='black')
        postion2 = postion

    # add the group name




    color_set = ['#30BBCE','#F3740B','#D2352C','#4A7DB4','#68AC56','#808080']
    color_set2 = ['#D2352C','#34587F','#4A7DB4','#68AC56','#808080','#FF7F7F']
    gropus = [['HTB', 'Xu4436'], ['TKS', 'Xu4273', 'Xu4319'], ['Xu4243', 'Xu4267', 'Xu4217']]
    site = ['TKS', 'Xu4273', 'Xu4319', 'HTB', 'Xu4436','Xu4243', 'Xu4267', 'Xu4217','T37','T28','T34','Xu4624','Xu7165','Xu269','Xu513','291','Xu254','Xu627','Xu552','Xu407','QTH','Xu320','TM','Xu3748', 'Xu6531','Xu6547','jzg']
    # here we plot the data 

    # off the axis
    postion = 0
#    for i in range(len(gropus)):
#        postion += 0.2
    for j in range(len(site)):
        postion += 0.5
        print(j)
        # here the j is the site name
        site_sample_number = len(data_k6[site[j]])
        print(site_sample_number)
        for k in range(site_sample_number):
            box_data_k6 = data_k6[site[j]][k][2]
            start = 3
            for m in range(6):
                ax.bar(postion, float(box_data_k6[m]), bottom=start, color=color_set[m], edgecolor='white', width=1)
                start += float(box_data_k6[m])
            postion += 1
        
        # add the group name
        box = plt.Rectangle((postion-site_sample_number-0.5, 4.02), site_sample_number, 0.05, fc="black", alpha=0.5)
        ax.add_patch(box)
        if site[j] == 'Xu254':
            print(site[j])
            #ax.text(postion-site_sample_number/2-0.5, 4.15, site[j], ha='center', va='center', fontsize=12, color='black')
            new_name = rename_dict[site[j]]
            ax.text(postion-site_sample_number/2-0.5, 4.15, new_name, ha='center', va='center', fontsize=12, color='black')
        elif site[j] == '291':
            #ax.text(postion-site_sample_number/2-0.5, 4.15, "Xu291", ha='center', va='center', fontsize=12, color='black')
            new_name = rename_dict[site[j]]
            ax.text(postion-site_sample_number/2-0.5, 4.15, new_name, ha='center', va='center', fontsize=12, color='black')
        else:
            #ax.text(postion-site_sample_number/2-0.5, 4.15, site[j], ha='center', va='center', fontsize=13, color='black')
            new_name = rename_dict[site[j]]
            ax.text(postion-site_sample_number/2-0.5, 4.15, new_name, ha='center', va='center', fontsize=12, color='black')
    
    # add the legend 
    # k = 6 
    ax.text(-1, 3.5, "K=6", ha='center', va='center', fontsize=18, color='black',rotation=90)



    color_set = ['#4A7DB4','#68AC56','#808080','#D2352C','#F3740B']
    color_set2 = ['#D2352C','#34587F','#4A7DB4','#68AC56','#808080','#FF7F7F']
    gropus = [['HTB', 'Xu4436'], ['TKS', 'Xu4273', 'Xu4319'], ['Xu4243', 'Xu4267', 'Xu4217']]
    site = ['TKS', 'Xu4273', 'Xu4319', 'HTB', 'Xu4436','Xu4243', 'Xu4267', 'Xu4217','T37','T28','T34','Xu4624','Xu7165','Xu269','Xu513','291','Xu254','Xu627','Xu552','Xu407','QTH','Xu320','TM','Xu3748', 'Xu6531','Xu6547','jzg']
    # here we plot the data 
    #fig, ax = plt.subplots(figsize=(45, 2))
    #plt.axis('off')
    # off the axis
    postion = 0
#    for i in range(len(gropus)):
#        postion += 0.2
    for j in range(len(site)):
        postion += 0.5
        print(j)
        # here the j is the site name
        site_sample_number = len(data_k5[site[j]])
        print(site_sample_number)
        for k in range(site_sample_number):
            box_data_k5 = data_k5[site[j]][k][2]
            start = 1.5
            for m in range(5):
                ax.bar(postion, float(box_data_k5[m]), bottom=start, color=color_set[m], edgecolor='white', width=1)
                start += float(box_data_k5[m])
            postion += 1
        
        # add the group name
        box = plt.Rectangle((postion-site_sample_number-0.5, 2.52), site_sample_number, 0.05, fc="black", alpha=0.5)
        ax.add_patch(box)
        # add the site name
        if site[j] == 'Xu254':
            print(site[j])
            #ax.text(postion-site_sample_number/2-0.5, 2.65, site[j], ha='center', va='center', fontsize=12, color='black')
            new_name = rename_dict[site[j]]
            ax.text(postion-site_sample_number/2-0.5, 2.65, new_name, ha='center', va='center', fontsize=12, color='black')
        elif site[j] == '291':
            #ax.text(postion-site_sample_number/2-0.5, 2.65, "Xu291", ha='center', va='center', fontsize=12, color='black')
            new_name = rename_dict[site[j]]
            ax.text(postion-site_sample_number/2-0.5, 2.65, new_name, ha='center', va='center', fontsize=12, color='black')
        else:
            #ax.text(postion-site_sample_number/2-0.5, 2.65, site[j], ha='center', va='center', fontsize=13, color='black')
            new_name = rename_dict[site[j]]
            ax.text(postion-site_sample_number/2-0.5, 2.65, new_name, ha='center', va='center', fontsize=12, color='black')
    #add the legend
    ax.text(-1, 2.0, "K=5", ha='center', va='center', fontsize=18, color='black',rotation=90)


    color_set = ['#808080','#D2352C','#68AC56','#4A7DB4']
    color_set2 = ['#D2352C','#34587F','#4A7DB4','#68AC56']
    gropus = [['HTB', 'Xu4436'], ['TKS', 'Xu4273', 'Xu4319'], ['Xu4243', 'Xu4267', 'Xu4217']]
    site = ['TKS', 'Xu4273', 'Xu4319', 'HTB', 'Xu4436','Xu4243', 'Xu4267', 'Xu4217','T37','T28','T34','Xu4624','Xu7165','Xu269','Xu513','291','Xu254','Xu627','Xu552','Xu407','QTH','Xu320','TM','Xu3748', 'Xu6531','Xu6547','jzg']
    # here we plot the data 
    # off the axis
    postion = 0
#    for i in range(len(gropus)):
#        postion += 0.2
    for j in range(len(site)):
        postion += 0.5
        print(j)
        # here the j is the site name
        site_sample_number = len(data_k4[site[j]])
        print(site_sample_number)
        for k in range(site_sample_number):
            box_data_k4 = data_k4[site[j]][k][2]
            start = 0
            for m in range(4):
                ax.bar(postion, float(box_data_k4[m]), bottom=start, color=color_set[m], edgecolor='white', width=1)
                start += float(box_data_k4[m])
            postion += 1
        
        # add the group name
        box = plt.Rectangle((postion-site_sample_number-0.5, 1.02), site_sample_number, 0.05, fc="black", alpha=0.5)
        ax.add_patch(box)
        # add the site name
        if site[j] == 'Xu254':
            print(site[j])
            #ax.text(postion-site_sample_number/2-0.5, 1.15, site[j], ha='center', va='center', fontsize=12, color='black')
            new_name = rename_dict[site[j]]
            ax.text(postion-site_sample_number/2-0.5, 1.15, new_name, ha='center', va='center', fontsize=12, color='black')
        elif site[j] == '291':
            #ax.text(postion-site_sample_number/2-0.5, 1.15, "Xu291", ha='center', va='center', fontsize=12, color='black')
            new_name = rename_dict[site[j]]
            ax.text(postion-site_sample_number/2-0.5, 1.15, new_name, ha='center', va='center', fontsize=12, color='black')
        else:
            #ax.text(postion-site_sample_number/2-0.5, 1.15, site[j], ha='center', va='center', fontsize=13, color='black')
            new_name = rename_dict[site[j]]
            ax.text(postion-site_sample_number/2-0.5, 1.15, new_name, ha='center', va='center', fontsize=12, color='black')
    # add the legend
    ax.text(-1, 0.5, "K=4", ha='center', va='center', fontsize=18, color='black',rotation=90)


    plt.savefig("k4_k5_k6.svg", dpi=600, bbox_inches='tight')
 
 
def main():
    samples = get_sample_k6()
    q = get_q_k6()
    data_k6 = get_order_k6(q, samples)
    
    samples_k5 = get_sample_k5()
    q_k5 = get_q_k5()
    data_k5 = get_order_k5(q_k5, samples_k5)

    samples_k4 = get_sample_k4()
    q_k4 = get_q_k4()
    data_k4 = get_order_k4(q_k4, samples_k4)

    #plot_admixture(data, args.output)
    plot_admixture_k4_k5_k6(data_k4, data_k5, data_k6)


if __name__ == '__main__':
    main()   