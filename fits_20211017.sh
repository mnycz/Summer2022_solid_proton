#!/bin/bash
#Script to carry out fits of sin^2th for 10/15/2021 updated unpolarized ep/eD simulations with unfolding
#NOTE: first argument to bash script must be "./eAll.Linux.x86_64.exe" which is the name of the fitting code executable

##### change the file location to where the bins file is located
$1 /w/eic-scshelf2104/users/gsevans/thirdWeekSULIs22/files/Small_Grid_Fixed_backup.csv
#$1 /w/eic-scshelf2104/users/gsevans/thirdWeekSULIs22/files/Small_Grid_Fixed.csv
#$1 /w/eic-scshelf2104/users/gsevans/thirdWeekSULIs22/files/rate_Q2x_acc_finebin_50uA_eAll_commited4fe_rod_6mm_11_LH2_100files_1e6.csv

#$1  /w/eic-scshelf2104/users/gsevans/5thWeekSULIs22/files/Small_Grid_Fixed.csv
wait
