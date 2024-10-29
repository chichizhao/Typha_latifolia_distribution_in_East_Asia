#!/bin/python
# function: rona calculate use python script 
# usage: python rona.py



import os
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy import stats
import statsmodels.api as sm




def rona_pred(gen, present, future, significance, weighted):
    # gen, present, future are expected to be pandas DataFrame
    # significance, weighted are boolean values
    gen = pd.DataFrame(gen)
    present = pd.DataFrame(present)
    future = pd.DataFrame(future)

    rona_output = pd.DataFrame(np.nan, index=gen.index, columns=gen.columns)
    allfreq_output = pd.DataFrame(np.nan, index=gen.index, columns=gen.columns)
    reg_output = pd.DataFrame(np.nan, index=gen.index, columns=["statistic", "pvalue", "Rsquared", "adjRsquared", "sigma"])
    print(rona_output)

    for j in gen.columns:
        for i in gen.index:
            X = sm.add_constant(present)
            y = gen.loc[i]
            #print(X)
            #print(y)
            #print(i)
            #print(j)
            model = sm.OLS(y, X).fit()
            #print(model.params)
            rona_output.loc[i, j] = abs((model.params[1] * future.loc[j] + model.params[0]) - gen.loc[i, j])[0]
            #print(rona_output)
            allfreq_output.loc[i, j] = (model.params[1] * future.loc[j] + model.params[0])[0]
            reg_output.loc[i] = [model.fvalue, model.pvalues[1], model.rsquared, model.rsquared_adj, model.mse_resid]
    #print(reg_output)
    #print(rona_output)#
    #print(allfreq_output)
    SEM = rona_output.apply(lambda x: x.sem(), axis=1)
    #print(SEM)
    if significance == "TRUE":
        reg_output = reg_output[reg_output['pvalue'] < 0.05]
        rona_output = rona_output[reg_output.index]
        allfreq_output = allfreq_output[reg_output.index]

    avg_rona_output = pd.DataFrame(rona_output.mean(axis=1), columns=["Avg_unweighted_RONA"])
    #print(avg_rona_output)

    if weighted == "TRUE":
        avg_rona_output["Avg_weighted_RONA"] = rona_output.apply(lambda x: np.average(x, weights=reg_output.loc[x.index, "Rsquared"]), axis=1)

    rona_output.columns = ["rona_" + str(col) for col in rona_output.columns]
    allfreq_output.columns = ["pred_" + str(col) for col in allfreq_output.columns]
    final_output = [reg_output, rona_output, avg_rona_output, allfreq_output, SEM]
    return final_output



# here we need read the ref_feq current matrix and set the name of the matrix  of each population is the same as the name of the population
ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio1.matrix", index_col=0, sep=" ")
print(ref_freq)
# read the current environment matrix
current_env = pd.read_csv("/home/chichi/data/china/china3/future/rona/current_env.txt", index_col=0, sep="\t")
#print(current_env)
# we need rearrage the current_env matrix row's name order to be the same as the ref_freq matrix column's name order
order = ref_freq.columns
#print(order)
# change the order of the current_env matrix
current_env = current_env.reindex(order)
print(current_env)

# read the future env matrix 
ssp126_2081_2100 = pd.read_csv("/home/chichi/data/china/china3/future/rona/future_env_ssp126_2081_2100.txt", index_col=0, sep="\t")
#print(ssp126_2081_2100)
# change the order of the future_env matrix
ssp126_2081_2100 = ssp126_2081_2100.reindex(order)
#print(ssp126_2081_2100)

ssp126_2061_2080 = pd.read_csv("/home/chichi/data/china/china3/future/rona/future_env_ssp126_2061_2080.txt", index_col=0, sep="\t")
#print(ssp126_2061_2080)
# change the order of the future_env matrix
ssp126_2061_2080 = ssp126_2061_2080.reindex(order)
print(ssp126_2061_2080)

ssp370_2081_2100 = pd.read_csv("/home/chichi/data/china/china3/future/rona/future_env_ssp370_2081_2100.txt", index_col=0, sep="\t")
#print(ssp370_2081-2100)
# change the order of the future_env matrix
ssp370_2081_2100 = ssp370_2081_2100.reindex(order)
#print(ssp370_2081_2100)

ssp370_2061_2080 = pd.read_csv("/home/chichi/data/china/china3/future/rona/future_env_ssp370_2061_2080.txt", index_col=0, sep="\t")
#print(ssp370_2061_2080)
# change the order of the future_env matrix
ssp370_2061_2080 = ssp370_2061_2080.reindex(order)
#print(ssp370_2061_2080)


# here we need do this process in batch for all the bios and all the future env
future_env = [ssp126_2081_2100, ssp126_2061_2080, ssp370_2081_2100, ssp370_2061_2080]
future_env_name = ["ssp126_2081_2100", "ssp126_2061_2080", "ssp370_2081_2100", "ssp370_2061_2080"]


# for bio1
for i in range(4):
    future_env_bio1 = future_env[i].bio_1
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_1, future_env_bio1, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio1.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio1.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio1.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio1.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio1.txt", sep="\t")



# for bio2
ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio2.matrix", index_col=0, sep=" ")
for i in range(4):
    future_env_bio2 = future_env[i].bio_2
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_2, future_env_bio2, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio2.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio2.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio2.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio2.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio2.txt", sep="\t")

# for bio3
ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio3.matrix", index_col=0, sep=" ")
for i in range(4):
    future_env_bio3 = future_env[i].bio_3
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_3, future_env_bio3, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio3.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio3.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio3.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio3.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio3.txt", sep="\t")

# for bio4
ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio4.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio4 = future_env[i].bio_4
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_4, future_env_bio4, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio4.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio4.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio4.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio4.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio4.txt", sep="\t")

# for bio6
ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio6.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio6 = future_env[i].bio_6
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_6, future_env_bio6, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio6.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio6.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio6.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio6.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio6.txt", sep="\t")

# for bio7
ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio7.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio7 = future_env[i].bio_7
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_7, future_env_bio7, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio7.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio7.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio7.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio7.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio7.txt", sep="\t")

# for bio8

ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio8.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio8 = future_env[i].bio_8
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_8, future_env_bio8, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio8.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio8.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio8.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio8.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio8.txt", sep="\t")

# for bio9

ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio9.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio9 = future_env[i].bio_9
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_9, future_env_bio9, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio9.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio9.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio9.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio9.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio9.txt", sep="\t")


# for bio11

ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio11.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio11 = future_env[i].bio_11
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_11, future_env_bio11, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio11.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio11.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio11.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio11.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio11.txt", sep="\t")



# for bio12

ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio12.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio12 = future_env[i].bio_12
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_12, future_env_bio12, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio12.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio12.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio12.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio12.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio12.txt", sep="\t")

# for bio13
ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio13.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio13 = future_env[i].bio_13
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_13, future_env_bio13, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio13.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio13.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio13.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio13.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio13.txt", sep="\t")


# for bio14

ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio14.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio14 = future_env[i].bio_14
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_14, future_env_bio14, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio14.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio14.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio14.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio14.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio14.txt", sep="\t")

# for bio15

ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio15.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio15 = future_env[i].bio_15
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_15, future_env_bio15, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio15.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio15.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio15.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio15.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio15.txt", sep="\t")


# for bio16

ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio16.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio16 = future_env[i].bio_16
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_16, future_env_bio16, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio16.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio16.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio16.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio16.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio16.txt", sep="\t")


# for bio17

ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio17.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio17 = future_env[i].bio_17
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_17, future_env_bio17, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio17.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio17.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio17.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio17.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio17.txt", sep="\t")


# for bio18

ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio18.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio18 = future_env[i].bio_18
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_18, future_env_bio18, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio18.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio18.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio18.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio18.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio18.txt", sep="\t")
  
# for bio19

ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio19.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio19 = future_env[i].bio_19
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_19, future_env_bio19, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio19.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio19.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio19.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio19.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio19.txt", sep="\t")
 
# for bio5
ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio5.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio5 = future_env[i].bio_5
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_5, future_env_bio5, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio5.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio5.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio5.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio5.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio5.txt", sep="\t")

# for bio10

ref_freq = pd.read_csv("/home/chichi/data/china/china3/future/rona/bios/afs_bio10.matrix", index_col=0, sep=" ")

for i in range(4):
    future_env_bio10 = future_env[i].bio_10
    reg_output, rona_output, avg_rona_output, allfreq_output, SEM = rona_pred(ref_freq, current_env.bio_10, future_env_bio10, "FALSE", "FALSE")

    # save the results of the analysis to the file
    reg_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "reg_output_" + future_env_name[i] + "_bio10.txt", sep="\t")
    rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "rona_output_" + future_env_name[i] + "_bio10.txt", sep="\t")
    avg_rona_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "avg_rona_output_" + future_env_name[i] + "_bio10.txt", sep="\t")
    allfreq_output.to_csv("/home/chichi/data/china/china3/future/rona/" + "allfreq_output_" + future_env_name[i] + "_bio10.txt", sep="\t")
    SEM.to_csv("/home/chichi/data/china/china3/future/rona/" + "SEM_" + future_env_name[i] + "_bio10.txt", sep="\t")