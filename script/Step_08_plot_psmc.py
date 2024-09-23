#!/bin/python
# -*- coding: utf-8 -*-
# Path: script/plot_psmc.py
# Author: chichi
# function: plot psmc result with chichi's understanding
# reference: https://github.com/shengxinzhuan/easy_psmc_plot
# usage python plot_psmc.py -namelist popmap.txt -psmcdir psmc_result -out psmc_plot

import sys
import os
import argparse
import click
import matplotlib.pyplot as plt
import numpy as np

# we need set some parameters 
mutation_rate = 5.51e-9 # we estimate from the whole genome data
generation_time = 1 # we set it as 1 year for one generation
size = 100 # we set the bin size as 100 in default in psmc analysis
def get_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-namelist', '--namelist', type=str, help='the namelist of the samples')
    parser.add_argument('-psmcdir', '--psmcdir', type=str, help='the directory of the psmc result')
    parser.add_argument('-out', '--out', type=str, help='the output directory of the psmc plot')
    return parser



# here is the idea of the psmc plot
# we want to plot three regions of the psmc result in one figure 
# and mark the three regions with three different colors 
# and for the samples, from the same site, we may want to get the average plot of this site
# while the above may not that true, it is just a initial idea now.
def read_name_and_group(filename):
    name_group_dict = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            name = line.split('\t')[0]
            region = line.split('\t')[1]
            print(name, region)
            name_group_dict[name] = [region]
    return name_group_dict



# the following function is get the psmc result from the psmcfa file
def psmc_fun(filename, size, mutation, generation_time):
    # Read the raw file
    psmc_rlt = open(filename, "r")
    result = psmc_rlt.read()
    psmc_rlt.close()
    # Getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    #print(last_block)
    # this is the last block of psmc result. for example, the last block of psmc result is like this: the rpund of 50 
    last_block = last_block.split('\n')
    #print(last_block)
    # convert the last block to a list
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2] =="RS":
            #print(line)
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))
    #print(time_windows)
    #print(estimated_lambdas)        
    # based on the result of the time_windows and estimated_lambdas, get the information of the last round of psmc
    '''
    rs period   time    estimated_lambda  
    RS	0	0.000000	1407962.503178	0.000456	0.000000	0.000000
    '''
    # Getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the estimated parameters
    #print(result)
    result = result[-1].split('\n')[0]
    #print(result)
    # comments by chichi 20230-09-14, this result use the last round of psmc
    # it is the last pa line, so it is the last round of psmc
    result = result.split(' ')
    #print(result)
    theta = float(result[1])
    #print(theta)
    N0 = theta/(4*mutation)/size
    # here the size is the bin size, here it is 100 in default
    # mutation is the mutation rate we difined
    # theta is pa line in the last block of psmc result, and for the NO, the theta is the first element of the pa line 
    # Scalling times and sizes
    # the method one we use to estimate the time and size, we know the mutation rate.
    times = [generation_time * 2 * N0 * i for i in time_windows]
    # set the scale of the last round of psmc, here the generation time is 1.
    sizes = [N0 * i for i in estimated_lambdas]
    # we get the real population size that we estimate.
    #print(times)
    #print(sizes)
    # Remove the false positive result
    raw_dict = {}
    false_result = sizes[-1]
    #print(false_result)
    for i in range(len(sizes)):
        raw_dict[times[i]] = sizes[i]
    times = []
    sizes = []
    for k, v in raw_dict.items():
        if str(v) != str(false_result):
            times.append(k)
            sizes.append(v)
        else:
            break
    #print(times)
    #print(sizes)
    return(times, sizes)

def plot_result(result_dict, out):
    fig = plt.figure(figsize=(16, 10))
    plt.rcParams['font.sans-serif'] = ['Arial']

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Years before present', fontsize=14)
    plt.ylabel('Effective population size', fontsize=14)
    plt.xlim(500, 3.5e6)
    plt.ylim(100, 3e8)
    box = plt.Rectangle((500,100), 3.5e6-500, 3e8-100, facecolor='none', edgecolor='black', linewidth=2)
    plt.gca().add_patch(box)
    
    # add the climate date of last  800 kya
    import pandas as pd
    import numpy as np

    data = pd.read_csv('/home/chichi/data/china/china3/momi/EDC_dD_temp_estim.tab', sep='\t', header=0)
    print(data.head())

    # remove the rows with missing values
    data = data.dropna()
    print(data.head())
    age = data['Age model [ka]']*1000
    T = (data['delta T [°C]']/10) +7.5 
    # turn in to exponential scale in 10
    T = 10**T
    line = plt.plot(age, T, color='black', linewidth=0.5,alpha=0.5)
    ax = plt.gca().add_line(line[0])
    
    box = plt.Rectangle((800000, 2e6), 2000, 2e8-2e6, facecolor='none', edgecolor='black', linewidth=1)
    plt.gca().add_patch(box)
    for i in range(1,5):
        martT = 10**((7.5+(i*5-15)/10))
        Rectangle = plt.Rectangle((800000, martT), 50000, 1e4, facecolor='none', edgecolor='black', linewidth=1)
        plt.gca().add_patch(Rectangle)
        plt.text(1000000, martT, str(i*5-15), fontsize=14, va='center', ha='center')
    
    plt.text(1.3e6, 10**(7.2), 'Antarctic Temperature (°C) \n relative to present', fontsize=14, va='center', ha='center', rotation=90)
    # add the line of 0 degree
    line = plt.plot([500,8e5], [10**7.5, 10**7.5], color='black', linewidth=0.5,alpha=0.5, linestyle='--')
    plt.gca().add_line(line[0])
    
    # here we need set the color of the regions
    colors = { 'WA':'#F3740B', 'WB':'#EFA38A', 'WC':'#D2352C','BA':'#808080','EC':'#ADD9ED', 'EB':'#30BBCE','EA':'#4A7DB4','MA':'#68AC56'}

    WA = {}
    WB = {}
    WC = {}
    BA = {}
    EC = {}
    EB = {}
    EA = {}
    MA = {}
    WA_times = []
    WB_times = []
    WC_times = []
    BA_times = []
    EC_times = []
    EB_times = []
    EA_times = []
    MA_times = []
    
    
    EAST = {}
    WEST = {}
    MIDDLE = {}
    NORTH = {}
    EAST_times = []
    WEST_times = []
    MIDDLE_times = []
    NORTH_times = []

    for name in result_dict.keys():
        region = result_dict[name][0]
        times = result_dict[name][1]
        sizes = result_dict[name][2]
        # here we plot the psmc result with the times and sizes, and it is steps plot
        #print(times)
        #print(sizes)
        if region == 'WA':
            WA[name] = [times, sizes]
            WA_times = WA_times + times
        elif region == 'WB':
            WB[name] = [times, sizes]
            WB_times = WB_times + times
        elif region == 'WC':
            WC[name] = [times, sizes]
            WC_times = WC_times + times
        elif region == 'BA':
            BA[name] = [times, sizes]
            BA_times = BA_times + times
        elif region == 'EC':
            EC[name] = [times, sizes]
            EC_times = EC_times + times
        elif region == 'EB':
            EB[name] = [times, sizes]
            EB_times = EB_times + times
        elif region == 'EA':
            EA[name] = [times, sizes]
            EA_times = EA_times + times
        elif region == 'MA':
            MA[name] = [times, sizes]
            MA_times = MA_times + times
        plt.step(times, sizes, color=colors[region], alpha=0.5, label=region,linewidth=0.2)
    # sort the WA times
    WA_times = sorted(WA_times)
    # so we have the time points of the WA region
    # so we can calculate the each time region of the average size
    WA_average_size = calculate_average_zise(WA, WA_times)
    #print(average_size)
    WA_times, WA_average_size = smooth_carve(WA_times, WA_average_size)
    # add the average steps plot of the WA region
    plt.step(WA_times, WA_average_size, color=colors['WA'], label='WA', linewidth=5, alpha=0.7)    
    
    # sort the WB times
    WB_times = sorted(WB_times)
    WB_average_size = calculate_average_zise(WB, WB_times)
    WB_times, WB_average_size = smooth_carve(WB_times, WB_average_size)
    plt.step(WB_times, WB_average_size, color=colors['WB'], label='WB', linewidth=5, alpha=0.7)
    
    # sort the WC times
    WC_times = sorted(WC_times)
    WC_average_size = calculate_average_zise(WC, WC_times)
    WC_times, WC_average_size = smooth_carve(WC_times, WC_average_size)
    plt.step(WC_times, WC_average_size, color=colors['WC'], label='WC', linewidth=5, alpha=0.7)
    
    # sort the BA times
    BA_times = sorted(BA_times)
    BA_average_size = calculate_average_zise(BA, BA_times)
    BA_times, BA_average_size = smooth_carve(BA_times, BA_average_size)
    plt.step(BA_times, BA_average_size, color=colors['BA'], label='BA', linewidth=5, alpha=0.7)
    
    # sort the EC times
    EC_times = sorted(EC_times)
    EC_average_size = calculate_average_zise(EC, EC_times)
    EC_times, EC_average_size = smooth_carve(EC_times, EC_average_size)
    plt.step(EC_times, EC_average_size, color=colors['EC'], label='EC', linewidth=5, alpha=0.7)
    
    # sort the EB times
    EB_times = sorted(EB_times)
    EB_average_size = calculate_average_zise(EB, EB_times)
    EB_times, EB_average_size = smooth_carve(EB_times, EB_average_size)
    plt.step(EB_times, EB_average_size, color=colors['EB'], label='EB', linewidth=5, alpha=0.7)
    
    # sort the EA times
    EA_times = sorted(EA_times)
    EA_average_size = calculate_average_zise(EA, EA_times)
    EA_times, EA_average_size = smooth_carve(EA_times, EA_average_size)
    plt.step(EA_times, EA_average_size, color=colors['EA'], label='EA', linewidth=5, alpha=0.7)
    
    # sort the MA times
    MA_times = sorted(MA_times)
    MA_average_size = calculate_average_zise(MA, MA_times)
    MA_times, MA_average_size = smooth_carve(MA_times, MA_average_size)
    plt.step(MA_times, MA_average_size, color=colors['MA'], label='MA', linewidth=5, alpha=0.7)
    


    # set x axis and y axis font size
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # change the axis mark inside the plot
    plt.tick_params(axis='x', which='major', direction='in', length=8, width=1.5)
    plt.tick_params(axis='y', which='major', direction='in', length=8, width=1.5)
    # off the minor tick
    plt.tick_params(axis='x', which='minor', bottom=False, top=False)
    plt.tick_params(axis='y', which='minor', left=False, right=False)
    # set axis line width
    plt.rcParams['axes.linewidth'] = 1.5

    # add the lengend
    box = plt.Rectangle((3.3e5, 4e3), 3.2e5, 1.5e3, facecolor=colors['WA'], alpha=0.9)
    plt.text(7e5, 4e3, 'WA', fontsize=14)
    plt.gca().add_patch(box)
    box = plt.Rectangle((3.3e5, 2e3), 3.2e5, 0.6e3, facecolor=colors['WB'], alpha=0.9)
    plt.text(7e5, 2e3, 'WB', fontsize=14)
    plt.gca().add_patch(box)
    box = plt.Rectangle((3.3e5, 1e3), 3.2e5, 3e2, facecolor=colors['WC'], alpha=0.9)
    plt.text(7e5, 1e3, 'WC', fontsize=14)
    plt.gca().add_patch(box)
    box = plt.Rectangle((3.3e5, 5e2), 3.2e5, 1.5e2, facecolor=colors['BA'], alpha=0.9)
    plt.text(7e5, 5e2, 'BA', fontsize=14)
    plt.gca().add_patch(box)
    
    box = plt.Rectangle((1e6, 4e3), 9e5, 1.5e3, facecolor=colors['EC'], alpha=0.9)
    plt.text(2e6, 4e3, 'EC', fontsize=14)
    plt.gca().add_patch(box)
    box = plt.Rectangle((1e6, 2e3), 9e5, 0.6e3, facecolor=colors['EB'], alpha=0.9)
    plt.text(2e6, 2e3, 'EB', fontsize=14)
    plt.gca().add_patch(box)
    box = plt.Rectangle((1e6, 1e3), 9e5, 3e2, facecolor=colors['EA'], alpha=0.9)
    plt.text(2e6, 1e3, 'EA', fontsize=14)
    plt.gca().add_patch(box)
    box = plt.Rectangle((1e6, 5e2), 9e5, 1.5e2, facecolor=colors['MA'], alpha=0.9)
    plt.text(2e6, 5e2, 'MA', fontsize=14)
    plt.gca().add_patch(box)
    

    
######################################################
    # add three important time region# 
    # 6000 years ago : Mid-Holocene
    # add the region of the mid-holocene from 5000 to 7000 years ago
    #box = plt.Rectangle((5e3, 1e2), 2e3, 3e6-1e2, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.2) 
    #plt.gca().add_patch(box)
    arrow = plt.arrow(6e3, 1e2, 0, 1.2e7-1e2, head_width=1e3, head_length=8e6, fc='grey', ec='grey',width=500,alpha=0.1)
    plt.gca().add_patch(arrow)
    plt.text(6e3, 3e6, '~6 kya', fontsize=14, ha='left')
    plt.text(6e3, 1.5e6, 'Mid-Holocene', fontsize=14, ha='left')

    # Last glacial maximum (LGM; ~21,000 years BP)
    # add the region of the last glacial maximum from 20000 to 22000 years ago
    #box = plt.Rectangle((1.8e4, 1e2), 4e3, 3e6-1e2, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.2)
    #plt.gca().add_patch(box)
    arrow = plt.arrow(2.1e4, 1e2, 0, 2e6-1e2, head_width=4e3, head_length=5e5, fc='grey', ec='grey',width=2000,alpha=0.1)
    plt.gca().add_patch(arrow)
    plt.text(2.1e4, 1e6, '~21 kya', fontsize=14, ha='left')
    plt.text(2.1e4, 6e5, 'LGM', fontsize=14, ha='left')

    # 130000 years ago : Last inter-glacial (LIG; ~120,000 - 140,000 years BP)
    # add the region of the last inter-glacial from 120000 to 140000 years ago
    #box = plt.Rectangle((1.2e5, 1e2), 2e4, 3e6-1e2, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.2)
    #plt.gca().add_patch(box)
    arrow = plt.arrow(1.3e5, 1e2, 0, 3e6-1e2, head_width=2e4, head_length=1e6, fc='grey', ec='grey',width=9000,alpha=0.1)
    plt.gca().add_patch(arrow)
    plt.text(1.3e5, 2e6, '~130 kya', fontsize=14, ha='left')
    plt.text(1.3e5, 1.3e6, 'LIG', fontsize=14, ha='left')
    
    # add the year of some mark region for development
    # the division year  2 kya
    arrow = plt.arrow(2.4e3, 1e2, 0, 2.4e3-1e2, head_width=2e2, head_length=1e3, fc='#4A7DB4', ec='#4A7DB4',width=50,alpha=0.5)
    plt.gca().add_patch(arrow)
    plt.text(2.6e3, 1.5e2, '2.4 kya', fontsize=14, va='center')
    
    # the year of 1kya

    arrow = plt.arrow(1e3, 1e2, 0, 2.2e3-1e2, head_width=1e2, head_length=5e2, fc='#F3740B', ec='#F3740B',width=30,alpha=0.5)
    plt.gca().add_patch(arrow)
    plt.text(1e3, 1.5e2, '1 kya', fontsize=14, va='center')
    
    # add the cited date
    
    plt.text(4e4, 6e7, 'Jouzel et al., 2007', fontsize=14, ha='center')
    


    

    # 400 000 years ago : mark region
    #box = plt.Rectangle((4e5, 1e2), 1e5, 3e6-1e2, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.2)
    #plt.gca().add_patch(box)
######################################################


    plt.savefig("psmc_8pop.pdf", dpi=600, bbox_inches='tight')
    plt.close()

def calculate_average_zise(dict, times):
    # base one each some samples in the dict has it period of the psmc result.
    # the data has such character:
    # 1. it goes with step plot with corresponding times and sizes
    # 2. each samples's time line points are not the same
    # the solutions is gather all the time points of the samples, sort them.
    # from the very first time point, we calculate the average size of the this first time period.
    # and keep the thoughts going on. for the end of this study line, it may not that all the samples reach the end of the time line.
    # so we depends on how many sample go there and calculate the average size of this period to represent the average size of this period, also this population.
    # we get the sorted time points line: times
    # we get the dict of the samples: dict
    summary_times = times
    average_size = [0] * len(summary_times)
    average_sample_number = [0]* len(summary_times)
    #print(len(summary_times))
    #print(len(average_size))
    for sample in dict.keys():
        sample_times = dict[sample][0]
        sample_sizes = dict[sample][1]
        j = 0
        k = 0
        for i in range(len(summary_times)):
            #print(i)
            #print(j)
            if summary_times[i] == sample_times[j]:
                for t in range(k, i+1):
                    average_size[t] = average_size[t] + sample_sizes[j]
                    average_sample_number[t] = average_sample_number[t] + 1
                k = i
                j = j + 1
                if j == len(sample_times):
                    break
            elif summary_times[i] > sample_times[j]:
                break
    #print(average_size)
    #print(average_sample_number)
    average_size = [average_size[i]/average_sample_number[i] for i in range(len(average_size))]
    return average_size
def smooth_carve(time, average_zise):
    # here we want to smooth the carve of the psmc result
    # here is my solution: gather 4 time points,in to one and give the new time point the average size of the 4 time points
    step = 5
    smooth_time = []
    smooth_size = []
    for i in range(0, len(time), step):
        smooth_time.append(time[i])
        smooth_size.append(sum(average_zise[i:i+step])/step)
    return smooth_time, smooth_size

if __name__ == '__main__':
    args = get_argparser().parse_args()
    namelist = args.namelist
    psmcdir = args.psmcdir
    out = args.out
    name_group_dict = read_name_and_group(namelist)
    result_dict = {}
    for name in name_group_dict.keys():
        region = name_group_dict[name][0]
        psmcfile = os.path.join(psmcdir, name + '.psmc')
        print(psmcfile)
        times, sizes = psmc_fun(psmcfile, size, mutation_rate, generation_time)
        result_dict[name] = [region, times, sizes]
    #print(result_dict)
    plot_result(result_dict, out)

