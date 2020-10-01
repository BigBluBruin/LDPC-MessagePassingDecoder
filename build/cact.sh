#! /bin/sh



for cc in minsum
do 
for ii in 2.20

do
    cat Result_n_1296_k_648_check_$cc\_Para_$ii\0000_ind\_*.txt > result_n_2193_k_1161_SNR_$ii\_Rate_0.50_check_$cc\.txt
done
done
