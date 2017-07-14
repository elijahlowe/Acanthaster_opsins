
#awk command looks at column 3 (gene id) and take only the first occurrence 
awk '!a[$3]++' gbr_hmmer_seq_table.txt | grep -f generatic_tf_ids /dev/stdin | less