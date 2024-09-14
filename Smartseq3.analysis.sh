
# step1
cd ~/YJ/Smartseq3/L2cell_EML
nohup ~/miniconda3/envs/R4/bin/Rscript ~/YJ/Smartseq3/Smartseq3.merge_demultiplexed_fastq.edited.R \
--dir /home/jcbu/YJ/Smartseq3/L2cell_EML \
--pigz ~/miniconda3/bin/pigz --threads 20 &



# step2: run zUMIs
nohup bash ~/software/zUMIs/zUMIs.sh  \
-y ~/YJ/Smartseq3/L2cell_EML_output/Smartseq3.yaml &
