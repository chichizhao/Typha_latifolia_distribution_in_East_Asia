#!/bin/python
# -*- coding: utf-8 -*-
# function: with momi2 package to infer the divide time and population size
# Author: CHICHI
import momi
import logging
logging.basicConfig(level=logging.INFO,filename="momi.log")
sfs = momi.Sfs.load("path/momi/")
model = momi.DemographicModel(N_e=1e5, gen_time=1, muts_per_gen=5.51e-7)
model.set_data(sfs,2160000)
# INFO:momi.demo_model:{it: 10, KLDivergence: Autograd ArrayBox with value 3.280359771938302, N_WC: 124159.49499916375, N_WB: 122988.51267479421, N_WA: 142815.32950939925, N_BA: 300000.0000000004, N_EC: 20139.169080407835, N_EB: 202909.96327607005, N_EA: 5120.638660783585, N_MA: 9999.99999999999, t_split_BA_EC: 67624.14265097171, N_BAA: 98029356.75602588, t_split_WA_BA: 291218.1655915705, N_WAA: 10500109.105950443, t_split_WB_WA: 331591.9753108786, N_WBA: 1648907.6671376142, t_split_WC_WB: 391898.20493432484, N_WCA: 278525.3296718002, t_split_EA_EB: 22841.697318267867, N_EAA: 1998777.1977695106, t_split_MA_EA: 83268.00840613396, N_MAA: 3765772.2333919, t_split_WC_MA: 599562.0179841466, N_A: 759913.9681148246, p_EB_EC: 0.2, p_MA_EA: 0.5, t_pulse_EB_EC: 12656.070419839421, t_pulse_MA_EA: 23607.474885973716}
# add the evolution leaves
model.add_leaf("WC", N="N_WC")
model.add_leaf("WB", N="N_WB")
model.add_leaf("WA", N="N_WA")
model.add_leaf("BA", N="N_BA")
model.add_leaf("EC", N="N_EC")
model.add_leaf("EB", N="N_EB")
model.add_leaf("EA", N="N_EA")
model.add_leaf("MA", N="N_MA")
# add the evolution relationship
model.move_lineages("EC", "BA", t="t_split_BA_EC",N="N_BAA")
model.move_lineages("BA", "WA", t="t_split_WA_BA",N="N_WAA")
model.move_lineages("WA", "WB", t="t_split_WB_WA",N="N_WBA")
model.move_lineages("WB", "WC", t="t_split_WC_WB",N="N_WCA")
model.move_lineages("EB", "EA", t="t_split_EA_EB",N="N_EAA")
model.move_lineages("EA", "MA", t="t_split_MA_EA",N="N_MAA")
model.move_lineages("MA", "WC", t="t_split_WC_MA",N="N_A")
# add migration
model.move_lineages("EC", "EB", t="t_pulse_EB_EC", p="p_EB_EC")
model.move_lineages("EA", "MA", t="t_pulse_MA_EA", p="p_MA_EA")

# Add the intitial guess for the parameters
model.add_size_param("N_WC", lower=7e4, upper=1.5e5)
model.add_size_param("N_WB", lower=1e5, upper=1.5e5)
model.add_size_param("N_WA", lower=1e5, upper=1.5e5)
model.add_size_param("N_BA", lower=1.5e5, upper=3e5)
model.add_size_param("N_EC", lower=2e4, upper=5e4)
model.add_size_param("N_EB", lower=2e5, upper=3e5)
model.add_size_param("N_EA", lower=1e3, upper=2e4)
model.add_size_param("N_MA", lower=1e4, upper=5e4)
model.add_time_param("t_split_BA_EC", lower=5e4, upper=8e4)
model.add_size_param("N_BAA", 1e8)
model.add_time_param("t_split_WA_BA", lower=2.5e5, upper=3e5)
model.add_size_param("N_WAA", 1e7)
model.add_time_param("t_split_WB_WA", lower=3.3e5, upper=3.5e5)
model.add_size_param("N_WBA", 1e4)
model.add_time_param("t_split_WC_WB", lower=3.6e5, upper=4e5)
model.add_size_param("N_WCA", 1e4)
model.add_time_param("t_split_EA_EB", lower=2e4, upper=5e4)
model.add_size_param("N_EAA", 1e4)
model.add_time_param("t_split_MA_EA", lower=8e4, upper=1.6e5)
model.add_size_param("N_MAA", 1e6)
model.add_time_param("t_split_WC_MA", lower=5e5, upper=6.5e5)
model.add_size_param("N_A", 1e5)
model.add_pulse_param("p_EB_EC", lower=0.15, upper=0.20)
model.add_pulse_param("p_MA_EA", lower=0.2, upper=0.5)
model.add_time_param("t_pulse_EB_EC", lower=1e4, upper=3.5e4)
model.add_time_param("t_pulse_MA_EA", lower=1e4, upper=2.5e4)


results = []
model_copy = model.copy()
model_copy.set_params(model.get_params(),randomize=True)
results.append(model_copy.optimize(method="L-BFGS-B"),options={'maxiter': 2,"ftol":1e-7})

