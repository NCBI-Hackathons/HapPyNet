head -5 ~/SRR_Acc_List.txt | parallel --verbose "python read_sra.py"
