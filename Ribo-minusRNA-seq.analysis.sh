####################################################################################################
#########################   RNA-seq analysis pipeline  without spike-in     ########################
####################################################################################################

###file_path
input_path=~/rawdata/
file_path=~/RNAseq_batch/
raw_dir=${file_path}/00_rawdata
logs_dir=${file_path}/logs
fastqc_dir=${file_path}/01_rawdata_QC
trimmedFastq_dir=${file_path}/02_trimmed_data
trimmedFastq_log_dir=${file_path}/logs/trimmed_data_logs
trimmeddata_fastqc_dir=${file_path}/03_trimmed_data_QC
rRNAdepletion_dir=${file_path}/04_rRNA_depletion
align_exp_dir=${file_path}/05_STAR_mm10
alignexp_log_dir=${file_path}/logs/align_mm10
bw_dir=${file_path}/06_bw
RSEM_dir=${file_path}/07_RSEM
sampleinfo=${file_path}/sampleinfo.txt


#step 0
####check files
mkdir -p $logs_dir
for id in `ls ${input_path}/md5.txt`
do
path1=$(dirname $id)
echo $path1
cd $path1
md5sum -c $id
done > ${logs_dir}/checkfiles.log


#step 1.1
####change file name#####
mkdir -p ${raw_dir}

cat $sampleinfo| while read id;
do
arr=($id)
sample1=${arr[0]}
sample2=${arr[1]}
fq1=$(ls ${input_path}/Rawdata/*/*_R1.fq.gz | grep "${sample1}_R1.fq.gz")
fq2=$(ls ${input_path}/Rawdata/*/*_R2.fq.gz | grep "${sample1}_R2.fq.gz")
ln -s $fq1 ${raw_dir}/${sample2}_R1.fastq.gz
ln -s $fq2 ${raw_dir}/${sample2}_R2.fastq.gz
done

#step 1.2
####fastqc of raw data ####
mkdir -p ${fastqc_dir}

in_path=${raw_dir}
out_path=${fastqc_dir}

nohup_number=0
for file in `ls ${in_path}/*R1.fastq.gz`
do
    ID=$(basename $file)
    fq1=${file}
    fq2=${file/R1/R2}
    if [ ! -s ${out_path}/"${ID/.fastq.gz/_fastqc.zip}" ]
    then
    echo "Generating file: $path/00-1_rawFastqc/"${ID}_R1_fastqc.zip"..."
    fastqc $fq1 -t 5  -o ${out_path}/  &
    fastqc $fq2 -t 5  -o ${out_path}/  &
    nohup_number=`echo $nohup_number+2 | bc`
    fi

    if [[ $nohup_number -eq 28 ]]
    then
        wait
        echo "waiting..."
        nohup_number=0
    fi
done

wait

#step 1.3
### merge reports of fastqc
multiqc ${fastqc_dir}/ -n rawdata_multiqc -o ${fastqc_dir}/


#step 2.1
### Trimming adapters  (trimmomatic)
mkdir -p ${trimmedFastq_dir}
mkdir -p ${trimmedFastq_log_dir}

nohup_number=0
for fq1 in `ls ${raw_dir}/*R1.fastq.gz`
do
sample="$(basename ${fq1/_R1.fastq.gz/})"
fq2=${fq1/_R1.fastq.gz/_R2.fastq.gz}
    if [ ! -s ${trimmedFastq_dir}/${sample}.R1.paired.fq.gz ]
    then
     java -jar ~/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 10 -phred33 \
     $fq1 $fq2 \
     ${trimmedFastq_dir}/${sample}.R1.paired.fq.gz \
     ${trimmedFastq_dir}/${sample}.R1.unpaired.fq.gz \
     ${trimmedFastq_dir}/${sample}.R2.paired.fq.gz \
     ${trimmedFastq_dir}/${sample}.R2.unpaired.fq.gz \
     ILLUMINACLIP:/home/jcbu/software/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    > ${trimmedFastq_log_dir}/${sample}.log 2>${trimmedFastq_log_dir}/${sample}.err &

     nohup_number=`echo $nohup_number+1 | bc`
    fi


    if [[ $nohup_number -eq 28 ]]
    then
        wait
        echo "waiting..."
        nohup_number=0
    fi
done

wait


#step 2.2
#####qc for trimmed data####
mkdir -p ${trimmeddata_fastqc_dir}

nohup_number=0
for file in `ls ${trimmedFastq_dir}/*.R1.paired.fq.gz`
do
    ID=$(basename $file)
    fq1=${file}
    fq2=${file/R1/R2}
    if [ ! -s ${trimmeddata_fastqc_dir}/"${ID/.paired.fq.gz/_fastqc.html}" ]
    then
    echo "Generating file: ${trimmeddata_fastqc_dir}/"${ID/.paired.fq.gz/_fastqc.html}"..."
    fastqc $fq1 -t 5  -o ${trimmeddata_fastqc_dir}/  &
    fastqc $fq2 -t 5  -o ${trimmeddata_fastqc_dir}/  &

     nohup_number=`echo $nohup_number+2 | bc`
    fi

    if [[ $nohup_number -eq 28 ]]
    then
        wait
        echo "waiting..."
        nohup_number=0
    fi
done



wait

### merge reports of fastqc
multiqc ${trimmeddata_fastqc_dir}/ -n trimmeddata_multiqc -o ${trimmeddata_fastqc_dir}/


#step 2.3
#####rRNA depletion of trimmed data####
mkdir -p ${rRNAdepletion_dir}

for fq1 in `ls ${trimmedFastq_dir}/*.R1.paired.fq.gz`
do
fq2=${fq1/R1.paired.fq.gz/R2.paired.fq.gz}
sample="$(basename ${fq1/.R1.paired.fq.gz/})"
if [ ! -s ${rRNAdepletion_dir}/${sample}.rmrRNA.R1.fq.gz ]
then
  ~/software/bbmap/bbduk.sh  \
  in1=$fq1 \
  in2=$fq2 \
  ref=${mouse_rRNA_ref}  \
  out1=${rRNAdepletion_dir}/${sample}.rmrRNA.R1.fq.gz \
  out2=${rRNAdepletion_dir}/${sample}.rmrRNA.R2.fq.gz \
  outm1=${rRNAdepletion_dir}/${sample}.rRNA.R1.fq.gz \
  outm2=${rRNAdepletion_dir}/${sample}.rRNA.R2.fq.gz \
  2> ${rRNAdepletion_dir}/${sample}.log

fi
done

wait

#step 3.1
###Aligning to experimental genome#####
echo -e "\n***************************\nalign of experimental genome begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${align_exp_dir}
mkdir -p ${alignexp_log_dir}

for fq1 in `ls ${rRNAdepletion_dir}/*.rmrRNA.R1.fq.gz`
do
fq2=${fq1/.rmrRNA.R1.fq.gz/.rmrRNA.R2.fq.gz}
sample="$(basename ${fq1/.rmrRNA.R1.fq.gz/})"
mkdir -p ${align_exp_dir}/${sample}
if [ ! -s ${align_exp_dir}/${sample}/${sample}Aligned.toTranscriptome.out.bam ]
then
    STAR  \
    --runThreadN 20 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --genomeDir    ${STAR_index}   \
    --readFilesIn $fq1 $fq2  \
    --readFilesCommand  zcat \
    --outFileNamePrefix ${align_exp_dir}/${sample}/${sample} \
    --outStd BAM_SortedByCoordinate \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --quantMode TranscriptomeSAM GeneCounts \
    1> ${align_exp_dir}/${sample}/${sample}.srt.bam \
    2> ${align_exp_dir}/${sample}/${sample}_STARalign.err

    samtools flagstat ${align_exp_dir}/${sample}/${sample}.srt.bam > ${alignexp_log_dir}/${sample}.stat
    samtools view -h -q 255 ${align_exp_dir}/${sample}/${sample}.srt.bam | samtools view -bh \
    1>${align_exp_dir}/${sample}/${sample}.srt.unique.bam \
    2>${alignexp_log_dir}/${sample}.GetUniqueReads.log \
    samtools index  ${align_exp_dir}/${sample}/${sample}.srt.unique.bam

fi
done

wait

#step 4.1
### Making RPKM-normalized bigWig files
echo -e "\n***************************\ntracks need to be done....\n***************************"

mkdir -p ${bw_dir}

for bam_file in `ls ${align_exp_dir}`
do
  sample="$(basename ${bam_file})"
    if [ ! -s "${bw_dir}/${sample}.bw" ]
    then
        bamCoverage -b ${align_exp_dir}/${sample}/${sample}.srt.unique.bam \
        --binSize 5 \
        --normalizeUsing RPKM \
        -p 20 \
        -o ${bw_dir}/${sample}.bw 2>${bw_dir}/${sample}.log
    fi
done

wait

##step 5.1
## read count and FPKM calculation using RSEM
mkdir -p ${RSEM_dir}

cd ${RSEM_dir}
for bam_file in `ls ${align_exp_dir}`
do
  sample="$(basename ${bam_file})"
    if [ ! -s "${RSEM_dir}/${sample}.genes.results" ]
    then
        rsem-calculate-expression --alignments --paired-end -p 10 --no-bam-output  \
        ${align_exp_dir}/${sample}/${sample}Aligned.toTranscriptome.out.bam \
        ${RSEM_index} \
        ${sample}
    fi
done
