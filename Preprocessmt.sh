#QC
parallel -j 3 --xapply 'kneaddata -i {1} -i {2} -o kneaddata_out -v -db /software_users/liuxl/human_reference/human_transcriptome  -t 20 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output --trimmomatic /software_users/liuxl/Trimmomatic-0.36 --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50"' ::: *_R1.fq.gz ::: *_R2.fq.gz
 
#View original reads quality
cd kneaddata_out/
kneaddata_read_count_table --input kneaddata_out --output ./kneaddata_read_counts.txt
 
#level out
mkdir pairend/
mv *R1_kneaddata_paired_2.fastq pairend/
mv *R1_kneaddata_paired_1.fastq pairend/
mkdir log
mv *.log log
mv kneaddata_read_counts.txt log
rm *
cd pairend/
mkdir level_out
for i in `ls *_R1_kneaddata_paired_1.fastq |sed 's/_R1_kneaddata_paired_1.fastq//'`;do head -n 70511928 ${i}_R1_kneaddata_paired_1.fastq > level_out/${i}_level_1.fastq; done
for i in `ls *_R1_kneaddata_paired_1.fastq |sed 's/_R1_kneaddata_paired_1.fastq//'`;do head -n 70511928 ${i}_R1_kneaddata_paired_2.fastq > level_out/${i}_level_2.fastq; done
 
#Merge pairend data
cd level_out
mkdir cat_reads
for i in `ls *_1.fastq |sed 's/_1.fastq//'`;do cat ${i}_1.fastq ${i}_2.fastq | awk '{if(NR%4==1) print "@"NR;else print $0}' >cat_reads/${i}.fastq;done
 
#humann2
parallel -j 3 'humann2 --threads 15 --input {} --output humann2_out' ::: cat_reads/*.fastq
cd humann2_out
 
#Composition annotation and functional analysis
mkdir metaphlan 
cp *_humann2_temp/*_metaphlan_bowtie2.txt metaphlan/
rm -rf *_humann2_temp
for i in `ls *_pathabundance.tsv |sed 's/_pathabundance.tsv//'`;do metaphlan2.py --input_type bowtie2out metaphlan/${i}_metaphlan_bowtie2.txt > metaphlan/profiled_${i};done
merge_metaphlan_tables.py metaphlan/profi* >speciesmt.txt
for i in `ls *_pathabundance.tsv |sed 's/_pathabundance.tsv//'`;do humann2_renorm_table --input ${i} _pathabundance.tsv --output ${i}_pathabundance_relab.tsv --units relab;done
for i in `ls *_pathabundance.tsv |sed 's/_pathabundance.tsv//'`;do humann2_renorm_table --input ${i} _genefamilies.tsv --output ${i}_genefamilies_relab.tsv --units relab;done
humann2_join_tables --input /data/liuxl/tet_raw_data/RS19IMWJ10/kneaddata_out/pairend/level_out/humann2_out --output genefamiliesmt.tsv --file_name genefamilies_relab
humann2_join_tables --input /data/liuxl/tet_raw_data/RS19IMWJ10/kneaddata_out/pairend/level_out/humann2_out --output pathcoveragemt.tsv --file_name pathcoverage
humann2_join_tables --input /data/liuxl/tet_raw_data/RS19IMWJ10/kneaddata_out/pairend/level_out/humann2_out --output pathaboundancemt.tsv --file_name pathabundance_relab
 
#Compressing 
cd ..
gunzip cat_reads/*.fastq
gunzip *.fastq
gunzip ../*.fastq
