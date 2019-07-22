#!/bin/bash

# first run, two chains per run
# PCAdapt detected between two and three populations for mullus so we run the model for both options
BayeScan2.1_macos64bits mul_bayesc_input_2pop.txt -o mul_bayesc_2pop_run1_ch1 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100
BayeScan2.1_macos64bits mul_bayesc_input_2pop.txt -o mul_bayesc_2pop_run1_ch2 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100

BayeScan2.1_macos64bits mul_bayesc_input_3pop.txt -o mul_bayesc_3pop_run1_ch1 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100
BayeScan2.1_macos64bits mul_bayesc_input_3pop.txt -o mul_bayesc_3pop_run1_ch2 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100

BayeScan2.1_macos64bits dip_bayesc_input_2pop.txt -o dip_bayesc_2pop_run1_ch1 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100
BayeScan2.1_macos64bits dip_bayesc_input_2pop.txt -o dip_bayesc_2pop_run1_ch2 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100

# second run, two chains per run
# here we try with a higher prior odds
BayeScan2.1_macos64bits mul_bayesc_input_2pop.txt -o mul_bayesc_2pop_run2_ch1 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10000
BayeScan2.1_macos64bits mul_bayesc_input_2pop.txt -o mul_bayesc_2pop_run2_ch2 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10000

BayeScan2.1_macos64bits mul_bayesc_input_3pop.txt -o mul_bayesc_2pop_run2_ch1 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10000
BayeScan2.1_macos64bits mul_bayesc_input_3pop.txt -o mul_bayesc_2pop_run2_ch2 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10000

BayeScan2.1_macos64bits dip_bayesc_input_2pop.txt -o dip_bayesc_2pop_run2_ch1 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10000
BayeScan2.1_macos64bits dip_bayesc_input_2pop.txt -o dip_bayesc_2pop_run2_ch2 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10000


# third run
# try again with prior odds = 100, and a greater number of iterations

BayeScan2.1_macos64bits mul_bayesc_input_2pop.txt -o mul_bayesc_2pop_run3_ch1 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100
BayeScan2.1_macos64bits mul_bayesc_input_2pop.txt -o mul_bayesc_2pop_run3_ch2 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100

BayeScan2.1_macos64bits mul_bayesc_input_3pop.txt -o mul_bayesc_3pop_run3_ch1 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100
BayeScan2.1_macos64bits mul_bayesc_input_3pop.txt -o mul_bayesc_3pop_run3_ch2 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100

BayeScan2.1_macos64bits dip_bayesc_input_2pop.txt -o dip_bayesc_2pop_run3_ch1 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100
BayeScan2.1_macos64bits dip_bayesc_input_2pop.txt -o dip_bayesc_2pop_run3_ch2 -n 10000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100