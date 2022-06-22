#!/bin/bash
#Script to carry out fits of sin^2th for 10/15/2021 updated unpolarized ep/eD simulations with unfolding
#NOTE: first argument to bash script must be "./eAll.Linux.x86_64.exe" which is the name of the fitting code executable


$1 ep 36.8 0 true /w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_ep_allq2_w4PDF_5_100.dat
wait
$1 ep 44.8 1 true /w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_ep_allq2_w4PDF_10_100.dat
wait
$1 ep 100.0 2 true /w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_ep_allq2_w4PDF_10_275.dat
wait
$1 ep 15.4 3 true /w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_ep_allq2_w4PDF_18_275.dat
wait
$1 ep 100.0 4 true /w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_ep_allq2_w4PDF_18_275.dat
wait

#$1 ed 36.8 true /w/halla-scifs17exp/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_eD_allq2_w4PDF_5_41.dat
#wait

#$1 ed 36.8 0 true /w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_eD_allq2_w4PDF_5_100.dat  
#wait 
#$1 ed 44.8 1 true /w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_eD_allq2_w4PDF_10_100.dat
#wait
#$1 ed 100.0 2 true /w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_eD_allq2_w4PDF_10_137.dat
#wait
#$1 ed 15.4 3 true /w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_eD_allq2_w4PDF_18_137.dat
#wait
#$1 ed 100.0 4 true /w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_eD_allq2_w4PDF_18_137.dat
#wait
#$1 ed 10.0 5 true /w/halla-scshelf2102/triton/mnycz/EIC/code_Alex/CT18NLO/New/asymmetry_unfolded4_eD_allq2_w4PDF_18_137.dat
#wait

