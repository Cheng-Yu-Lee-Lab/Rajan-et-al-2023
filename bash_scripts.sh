##Transcription factor analysis 
#Loop through directories containing .fastq files for all CUT&RUN data to generate bam files
for d in */* 
	do	
		files=($d/*)		
		fastqc ${files[0]}
		fastqc ${files[1]}
		name=${d#*/}
		cutadapt -j 32 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${d}/R1.trimmed.fastq -p ${d}/R2.trimmed.fastq ${files[0]} ${files[1]} --json=${d}/cutadapt_report.json
		bowtie2 -x BDGP6 -1 ${d}/R1.trimmed.fastq -2 ${d}/R2.trimmed.fastq --local --very-sensitive --no-mixed --no-discordant --no-unal --dovetail -I 10 -X 1000 -p 32 -S ${d}/aligned.sam > ${d}/bowtie_err.txt 2> ${d}/${name}_bowtie 
		samtools view -bS ${d}/aligned.sam > ${d}/aligned.bam -@ 32
		samtools sort ${d}/aligned.bam -o ${d}/aligned.sorted.bam -@ 32
		samtools index ${d}/aligned.sorted.bam
		
		#Generate initial bigwig file for visulaization 
		bamCoverage -b ${d}/aligned.sorted.bam -o ${d}/track.small.bw --binSize 5 -p max --normalizeUsing RPKM --ignoreDuplicates --maxFragmentLength 120
	done;

#Replicates were confirmed to be similar by correlation analysis and by inspection of tracks in genome, for visualization purposes, replicates were merged and merged bigwigs were generated
samtools merge mergedBams/$genotype.bam --threads 16 $rep1/aligned.sorted.bam $rep2/aligned.sorted.bam
samtools sort -o sortedBams/$genotype.bam --threads 16 mergedBams/$genotype.bam 
bamCoverage -b sortedBams/$sample.bam -o bigwigs/$sample.small.bw --binSize 5 -p max --normalizeUsing RPKM --ignoreDuplicates --maxFragmentLength 120

#z-score was calculated for all bigwigs in bigwigs/ using script in r_scripts.Rmd
multiBigwigSummary bins -b bigwigs/fruC-myc.small.bw bigwigs/fruCOM.small.bw -o figures/correlation_fru.npz -p max

#Generate correlation plot in FigS3.B
plotCorrelation -in figures/correlation_fru.npz -c pearson -p scatterplot --skipZeros --removeOutliers --log1p --xRange 0 9 --yRange 0 9 -o figures/correlation_fru_pearson.png

#For Transcription factor peak files, run modified CUT&RUNTools2 which is compatible with dm6 genome, keep peak file from peakcalling/macs2.narrow/{name}.narrowPeak
#In CUT&RUNTools2 json file, spike_in=FALSE, frag120=TRUE
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-fruCMYC_1.json 4362-AR-3_AGGAACAC-GACATTCC_S141
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-fruCMYC_2.json 5127-AR-27_GGTATAGG-CTGGAGTA_S65
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-fruCOM_1.json 6044-AR-1_TCGGATTC-TAGTTGCG_S96
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-fruCOM_2.json 6044-AR-2_CTGTACCA-AGTCTGTG_S97
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-SuH_1.json 6044-AR-10_GTCCTAAG-TGGTAGCT_S105
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-SuH_2.json 6044-AR-11_GCGTTAGA-GTATGCTG_S106
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-GAF_1.json 6044-AR-8_GCAATTCC-TGTTCGAG_S103
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-GAF_2.json 6044-AR-9_CTCAGAAG-AACTTGCC_S104
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-SuZ12_1.json 6675-AR-1_CGAATACG-GTCGGTAA_S94
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-SuZ12_2.json 6675-AR-2_GTCCTTGA-TCAGACGA_S95
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-caf1_1.json 6675-AR-5_GTCGATTG-ACGGTCTT_S98
./CUT-RUNTools-2.0-master/run_bulkModule.sh bulk-config-caf1_2.json 6675-AR-6_ATAACGCC-TGAACCTG_S99

#From H3K9me3 Cut&Run data, align using code above, call peaks with macs2
macs2 callpeak -t sortedBams/H3K9me3.bam -c 4362-iGG.bam -g 142573017 -n data/H3K9me3_peaks
awk '($5) > 100' data/H3K9me3_peaks.narrowPeak | bedtools merge -i stdin -d 10000 | awk '($3-$2) >= 3000' > data/heterochromatin_domains.bed

#Merge CUT&RUNTools2 replicates into one file and remove high signal and heterochromatin domains 
#Peak files were annotated with Homer (note, depending on algorithm, chromosomes need to have chr added or removed at start of chromosome names) 

bedtools sort -i data/fruC-myc_peaks/fruC-myc_cut_tools_merged.narrowPeak | awk '($5) > 100' | bedtools intersect -a stdin -b regions/dm6.blacklist.noCHR.bed -v | bedtools intersect -a stdin -b regions/dm6.TA_repeat.noCHR.bed -v | bedtools intersect -a stdin -b data/heterochromatin_domains.bed -v | bedtools merge -i stdin | cat | awk -F'\t' -v OFS='\t' '{ $1 = "chr" $1 }1' > data/fruC-myc_peaks/fruC-myc-PEAKS.bed
annotatePeaks.pl data/fruC-myc_peaks/fruC-myc-PEAKS.bed dm6 > data/fruC-myc_peaks/fruC-myc-PEAKS.txt

bedtools sort -i data/fruCOM_peaks/fruCOM_cut_tools_merged.narrowPeak | awk '($5) > 100' | bedtools intersect -a stdin -b regions/dm6.blacklist.noCHR.bed -v | bedtools intersect -a stdin -b regions/dm6.TA_repeat.noCHR.bed -v | bedtools intersect -a stdin -b data/heterochromatin_domains.bed -v | bedtools merge -i stdin | cat | awk -F'\t' -v OFS='\t' '{ $1 = "chr" $1 }1' > data/fruCOM_peaks/fruCOM-PEAKS.bed
annotatePeaks.pl data/fruCOM_peaks/fruCOM-PEAKS.bed dm6 > data/fruCOM_peaks/fruCOM-PEAKS.txt

bedtools sort -i data/SuH_peaks/SuH_cut_tools_merged.narrowPeak | bedtools intersect -a stdin -b regions/dm6.blacklist.noCHR.bed -v | bedtools intersect -a stdin -b regions/dm6.TA_repeat.noCHR.bed -v | bedtools intersect -a stdin -b data/heterochromatin_domains.bed -v | bedtools merge -i stdin | cat | awk -F'\t' -v OFS='\t' '{ $1 = "chr" $1 }1' > data/SuH_peaks/SuH-PEAKS.bed
annotatePeaks.pl data/SuH_peaks/SuH-PEAKS.bed dm6 > data/SuH_peaks/SuH-PEAKS.txt

bedtools sort -i data/GAF_peaks/GAF_cut_tools_merged.narrowPeak | bedtools intersect -a stdin -b regions/dm6.blacklist.noCHR.bed -v | bedtools intersect -a stdin -b regions/dm6.TA_repeat.noCHR.bed -v | bedtools intersect -a stdin -b data/heterochromatin_domains.bed -v | bedtools merge -i stdin | cat | awk -F'\t' -v OFS='\t' '{ $1 = "chr" $1 }1' > data/GAF_peaks/GAF-PEAKS.bed
annotatePeaks.pl data/GAF_peaks/GAF-PEAKS.bed dm6 > data/GAF_peaks/GAF-PEAKS.txt

#Random FruC-myc peak dataset generated for control analysis
bedtools shuffle -i data/fruC-myc_peaks/fruC-myc-PEAKS.noCHR.bed -excl regions/dm6.blacklist_heterochromatin.noCHR.bed -g regions/dm6.genome | cat | awk -F'\t' -v OFS='\t' '{ $1 = "chr" $1 }1' > data/random_peaks/randomPEAKS.bed
annotatePeaks.pl data/random_peaks/randomPEAKS.bed dm6 > data/random_peaks/randomPEAKS.txt

#Random Su(H) peak dataset generated for control analysis
bedtools shuffle -i data/SuH_peaks/SuH-PEAKS.noCHR.bed -excl regions/dm6.blacklist_heterochromatin.noCHR.bed -g regions/dm6.genome | cat | awk -F'\t' -v OFS='\t' '{ $1 = "chr" $1 }1' > data/random_peaks/randomPEAKS.SuH.bed
annotatePeaks.pl data/random_peaks/randomPEAKS.SuH.bed dm6 > data/random_peaks/randomPEAKS.SuH.txt


#Genomics distributions and peak overlaps were visualized from Homer annotations in Microsoft Excel (Fig3.A, D, E, F, H) 

#Bedfiles were seperated into promoters or regulatory based on Homer annotations 
#Heatmaps were generated using deeptools using code below: 

#Fig3B
computeMatrix reference-point \
 -S bigwigs/z-score/fruC-myc.small.bw  \
 -R data/ATAC_peaks/ATAC-promoter.bed data/ATAC_peaks/ATAC-regulatory.bed \
 --referencePoint center \
 -a 500 -b 500 \
 -out figures/ATAC.tab.gz \
 -p max
 
plotHeatmap \
 -m figures/ATAC.tab.gz \
 -out figures/ATAC \
 --legendLocation none \
 --xAxisLabel "Distance" \
 --heatmapHeight 10 \

#Fig3G
computeMatrix reference-point \
 -S  bigwigs/z-score/fruC-myc.small.bw  \
 -R data/SuH_peaks/SuH_promoter.bed data/SuH_peaks/SuH_regulatory.bed \
 --referencePoint center \
 -a 500 -b 500 \
 -out figures/SuH.tab.gz \
 -p max
 
plotHeatmap \
 -m figures/SuH.tab.gz \
 -out figures/SuH \
 --legendLocation none \
 --xAxisLabel "Distance" \
 --heatmapHeight 10 \

#Fig3I
computeMatrix reference-point \
 -S  bigwigs/z-score/fruC-myc.small.bw  \
 -R data/GAF_peaks/GAF_promoter.bed data/GAF_peaks/GAF_regulatory.bed \
 --referencePoint center \
 -a 500 -b 500 \
 -out figures/GAF.tab.gz \
 -p max
 
plotHeatmap \
 -m figures/GAF.tab.gz \
 -out figures/GAF \
 --legendLocation none \
 --xAxisLabel "Distance" \
 --heatmapHeight 10 \

#FigS3.E
#Motif analyis: Files were prepared for iCisTarget 
bedtools slop -b 200 -i data/fruCOM_peaks/fruCOM-peaks.bed -g regions/dm6.genome > data/motifs/iCisTarget/all-slopped.bed
#Bedfile was provided to https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/ 
#Motif analysis in HOMER:
findMotifsGenome.pl  data/motifs/iCisTarget/all-slopped.bed dm6 data/fruCOM_peaks/HOMER -size given

#Motif preparation for XSTREME. fasta file provided to https://meme-suite.org/meme/tools/xstreme
bedtools getfasta -fi regions/dm6.fa -bed data/motifs/iCisTarget/all-slopped.bed > data/XSTREME/fruC.fa

#Heatmaps for Fig6.B, C. fru_nb.bed was gotten from taking fruC-myc peaks which were annotated as being associated with a gene that is expressed uniquely in NBs (from scRNA-seq data, see scRNA-seq.ipynb for gene list generation)  
#fru_all.bed is equivalent to data/fruC-myc_peaks/fruC-myc-PEAKS.bed, but contains just the summits from the narrowPeak file 
computeMatrix reference-point \
 -S bigwigs/z-score/Caf1.small.bw bigwigs/z-score/SuZ12.small.bw  \
 -R data/fru_all.bed \
 --referencePoint center \
 -a 2000 -b 2000 \
 -out figures/final/PRC2.tab.gz \
 -p max
 
plotHeatmap \
 -m figures/final/PRC2.tab.gz \
 -out figures/final/PRC2_fruC \
 --legendLocation none \
 --xAxisLabel "Distance" \
 --heatmapHeight 10 \
 --yMin 0 \
 --yMax 3 \
 --zMin 0 \
 --zMax 5

computeMatrix reference-point \
 -S bigwigs/z-score/Caf1.small.bw bigwigs/z-score/SuZ12.small.bw  \
 -R data/fru_nb.bed \
 --referencePoint center \
 -a 2000 -b 2000 \
 -out figures/final/PRC2_nb.tab.gz \
 -p max
 
plotHeatmap \
 -m figures/final/PRC2_nb.tab.gz \
 -out figures/final/PRC2_nb_fruC \
 --legendLocation none \
 --xAxisLabel "Distance" \
 --heatmapHeight 10 \
 --yMin 0 \
 --yMax 3 \
 --zMin 0 \
 --zMax 5

#----------------------------------------------------------------------------------------------------------------------------------------------------------#
#For Histone data: 
##Transcription factor analysis 
#Loop through directories containing .fastq files for all CUT&RUN data to generate bam files and run peakcalling with gopeaks, and generate input files for diffReps
for d in */* 
	do	
		files=($d/*)		
		fastqc ${files[0]}
		fastqc ${files[1]}
		name=${d#*/}
		cutadapt -j 32 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${d}/R1.trimmed.fastq -p ${d}/R2.trimmed.fastq ${files[0]} ${files[1]} --json=${d}/cutadapt_report.json
		bowtie2 -x BDGP6 -1 ${d}/R1.trimmed.fastq -2 ${d}/R2.trimmed.fastq --local --very-sensitive --no-mixed --no-discordant --no-unal --dovetail -I 10 -X 1000 -p 32 -S ${d}/aligned.sam > ${d}/bowtie_err.txt 2> ${d}/${name}_bowtie 
		samtools view -bS ${d}/aligned.sam > ${d}/aligned.bam -@ 32
		samtools sort ${d}/aligned.bam -o ${d}/aligned.sorted.bam -@ 32
		samtools index ${d}/aligned.sorted.bam
		
		gopeaks -b ${d}/aligned.sorted.bam -c bam/4362-iGG.bam -o regions/gopeaks/broad/$name --chromsize regions/dm6.genome --broad 
        gopeaks -b ${d}/aligned.sorted.bam -c bam/4362-iGG.bam -o regions/gopeaks/narrow/$name --chromsize regions/dm6.genome 
		
		bedtools bamtobed -i ${d}/aligned.sorted.bam > data/diffReps_inputs/$name.bed
	done;
	
#Generate diffReps reports
diffReps.pl --treatment data/diffReps_inputs/H3K27ac-fruDF_rep1_5127.bed data/diffReps_inputs/H3K27ac-fruDF_rep2_5127.bed --control data/diffReps_inputs/H3K27ac_rep1_5127.bed data/diffReps_inputs/H3K27ac_rep2_5127.bed --pval 0.005 --report diffReps/Reports/H3K27ac_fruDelvCntrl --chrlen regions/dm6.genome --nproc 16 --window 500
diffReps.pl --treatment data/diffReps_inputs/H3K27me3-fruDF_rep1_5127.bed data/diffReps_inputs/H3K27me3-fruDF_rep2_5127.bed --control data/diffReps_inputs/H3K27me3_rep1_5127.bed data/diffReps_inputs/H3K27me3_rep2_5127.bed --pval 0.005 --report diffReps/Reports/H3K27me3_fruDelvCntrl --chrlen regions/dm6.genome --nproc 16 --window 500
diffReps.pl --treatment data/diffReps_inputs/H3K4me3-fruDF_rep1_6044.bed data/diffReps_inputs/H3K4me3-fruDF_rep2_6044.bed --control data/diffReps_inputs/H3K4me3_rep1_6044.bed data/diffReps_inputs/H3K4me3_rep2_6044.bed --pval 0.005 --report diffReps/Reports/H3K4me3_fruDelvCntrl --chrlen regions/dm6.genome --nproc 16 --window 500
#Visulazation of volcano plots generated in python file (python-scripts.ipynb) 

#Merge gopeaks reps for each histone mark
cat regions/gopeaks/broad/H3K27ac_rep1_5127_peaks.bed regions/gopeaks/broad/H3K27ac_rep2_5127_peaks.bed regions/gopeaks/broad/H3K27ac-fruDF_rep1_5127_peaks.bed regions/gopeaks/broad/H3K27ac-fruDF_rep2_5127_peaks.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > regions/gopeaks/merged/H3K27ac_broad.bed
cat regions/gopeaks/broad/H3K27me3_rep1_5127_peaks.bed regions/gopeaks/broad/H3K27me3_rep2_5127_peaks.bed regions/gopeaks/broad/H3K27me3-fruDF_rep1_5127_peaks.bed regions/gopeaks/broad/H3K27me3-fruDF_rep2_5127_peaks.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > regions/gopeaks/merged/H3K27me3_broad.bed
cat regions/gopeaks/narrow/H3K4me3_rep1_6044_peaks.bed regions/gopeaks/narrow/H3K4me3_rep2_6044_peaks.bed regions/gopeaks/narrow/H3K4me3-fruDF_rep1_6044_peaks.bed regions/gopeaks/narrow/H3K4me3-fruDF_rep2_6044_peaks.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > regions/gopeaks/merged/H3K4me3_narrow.bed

#Remove blacklisted regions, and prepare .saf file for TMM normalized bigwigs
for f in regions/gopeaks/merged/*.bed
	do
		temp=${f#*/}
		sample=$(basename "$f" .bed)
		echo $sample
		bedtools intersect -a $f -b regions/dm6.blacklist.noCHR.bed -v | bedtools intersect -a stdin -b regions/dm6.TA_repeat.noCHR.bed -v > regions/gopeaks/blacklisted/${sample}.blacklisted.bed
		awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' regions/gopeaks/merged/${sample}.bed > regions/gopeaks/saf/${sample}.saf
	done

#Calculate scalefactors in R (see R-scripts.rmd) using edgeR

#Make individual bigwig replicates
mkdir bigwigs/TMM_norm
bamCoverage --bam bam/H3K4me3_rep1_6044.bam -o bigwigs/TMM_norm/H3K4me3_rep1_6044.TMM.bw --binSize 5 --scaleFactor 0.13041953 --ignoreDuplicates --numberOfProcessors 16
bamCoverage --bam bam/H3K4me3_rep2_6044.bam -o bigwigs/TMM_norm/H3K4me3_rep2_6044.TMM.bw --binSize 5 --scaleFactor 0.12401646 --ignoreDuplicates --numberOfProcessors 16
bamCoverage --bam bam/H3K4me3-fruDF_rep1_6044.bam -o bigwigs/TMM_norm/H3K4me3-fruDF_rep1_6044.TMM.bw --binSize 5 --scaleFactor 0.06902642 --ignoreDuplicates --numberOfProcessors 16
bamCoverage --bam bam/H3K4me3-fruDF_rep2_6044.bam -o bigwigs/TMM_norm/H3K4me3-fruDF_rep2_6044.TMM.bw --binSize 5 --scaleFactor 0.08052365 --ignoreDuplicates --numberOfProcessors 16
bamCoverage --bam bam/H3K27ac_rep1_5127.bam -o bigwigs/TMM_norm/H3K27ac_rep1_5127.TMM.bw --binSize 5 --scaleFactor 0.11900730 --ignoreDuplicates --numberOfProcessors 16
bamCoverage --bam bam/H3K27ac_rep2_5127.bam -o bigwigs/TMM_norm/H3K27ac_rep2_5127.TMM.bw --binSize 5 --scaleFactor 0.02546946 --ignoreDuplicates --numberOfProcessors 16
bamCoverage --bam bam/H3K27ac-fruDF_rep1_5127.bam -o bigwigs/TMM_norm/H3K27ac-fruDF_rep1_5127.TMM.bw --binSize 5 --scaleFactor 0.07163300 --ignoreDuplicates --numberOfProcessors 16
bamCoverage --bam bam/H3K27ac-fruDF_rep2_5127.bam -o bigwigs/TMM_norm/H3K27ac-fruDF_rep2_5127.TMM.bw --binSize 5 --scaleFactor 0.11726911 --ignoreDuplicates --numberOfProcessors 16
bamCoverage --bam bam/H3K27me3_rep1_5127.bam -o bigwigs/TMM_norm/H3K27me3_rep1_5127.TMM.bw --binSize 5 --scaleFactor 0.12723577 --ignoreDuplicates --numberOfProcessors 16
bamCoverage --bam bam/H3K27me3_rep2_5127.bam -o bigwigs/TMM_norm/H3K27me3_rep2_5127.TMM.bw --binSize 5 --scaleFactor 0.07776931 --ignoreDuplicates --numberOfProcessors 16
bamCoverage --bam bam/H3K27me3-fruDF_rep1_5127.bam -o bigwigs/TMM_norm/H3K27me3-fruDF_rep1_5127.TMM.bw --binSize 5 --scaleFactor 0.34230829 --ignoreDuplicates --numberOfProcessors 16
bamCoverage --bam bam/H3K27me3-fruDF_rep2_5127.bam -o bigwigs/TMM_norm/H3K27me3-fruDF_rep2_5127.TMM.bw --binSize 5 --scaleFactor 0.16682195 --ignoreDuplicates --numberOfProcessors 16

#Merge bigwig replicatesmkdir bigwigs/TMM_norm/merged
bigwigCompare -b1 bigwigs/TMM_norm/H3K4me3_rep1_6044.TMM.bw -b2 bigwigs/TMM_norm/H3K4me3_rep2_6044.TMM.bw --operation mean --binSize 10 --blackListFileName regions/dm6.blacklist.noCHR.bed -p max -o bigwigs/TMM_norm/merged/H3K4me3.bw
bigwigCompare -b1 bigwigs/TMM_norm/H3K4me3-fruDF_rep1_6044.TMM.bw -b2 bigwigs/TMM_norm/H3K4me3-fruDF_rep2_6044.TMM.bw --operation mean --binSize 10 --blackListFileName regions/dm6.blacklist.noCHR.bed -p max -o bigwigs/TMM_norm/merged/H3K4me3-fruDF.bw

bigwigCompare -b1 bigwigs/TMM_norm/H3K27ac_rep1_5127.TMM.bw -b2 bigwigs/TMM_norm/H3K27ac_rep2_5127.TMM.bw --operation mean --binSize 10 --blackListFileName regions/dm6.blacklist.noCHR.bed -p max/2 -o bigwigs/TMM_norm/merged/H3K27ac.bw
bigwigCompare -b1 bigwigs/TMM_norm/H3K27ac-fruDF_rep1_5127.TMM.bw -b2 bigwigs/TMM_norm/H3K27ac-fruDF_rep2_5127.TMM.bw --operation mean --binSize 10 --blackListFileName regions/dm6.blacklist.noCHR.bed -p max -o bigwigs/TMM_norm/merged/H3K27ac-fruDF.bw

bigwigCompare -b1 bigwigs/TMM_norm/H3K27me3_rep1_5127.TMM.bw -b2 bigwigs/TMM_norm/H3K27me3_rep2_5127.TMM.bw --operation mean --binSize 10 --blackListFileName regions/dm6.blacklist.noCHR.bed -p max -o bigwigs/TMM_norm/merged/H3K27me3.bw
bigwigCompare -b1 bigwigs/TMM_norm/H3K27me3-fruDF_rep1_5127.TMM.bw -b2 bigwigs/TMM_norm/H3K27me3-fruDF_rep2_5127.TMM.bw --operation mean --binSize 10 --blackListFileName regions/dm6.blacklist.noCHR.bed -p max -o bigwigs/TMM_norm/merged/H3K27me3-fruDF.bw

#Make z-score bigwig files from merged bigwigs using R-script 

#Heatmaps for Fig5.E, 5SH. 
computeMatrix reference-point \
 -S  bigwigs/TMM_norm/merged/H3K27me3.bw bigwigs/TMM_norm/merged/H3K27me3-fruDF.bw  \
 -R data/fru_all.bed \
 --referencePoint center \
 -a 2000 -b 2000 \
 -out figures/final/fruALL_K27me3.tab.gz \
 -p max
     
plotHeatmap \
 -m figures/final/fruALL_K27me3.tab.gz \
 -out figures/final/fruALL_K27me3 \
 --legendLocation none \
 --xAxisLabel "Distance" \
 --heatmapHeight 10  \
 --yMin 0 \
 --yMax 6.5 \
 --zMin 0 \
 --zMax 16

computeMatrix reference-point \
 -S  bigwigs/TMM_norm/merged/H3K27me3.bw bigwigs/TMM_norm/merged/H3K27me3-fruDF.bw  \
 -R data/fru_nb.bed \
 --referencePoint center \
 -a 2000 -b 2000 \
 -out figures/final/fruNB_K27me3.tab.gz \
 -p max
     
plotHeatmap \
 -m figures/final/fruNB_K27me3.tab.gz \
 -out figures/final/fruNB_K27me3 \
 --legendLocation none \
 --xAxisLabel "Distance" \
 --heatmapHeight 10  \
 --yMin 0 \
 --yMax 6.5 \
 --zMin 0 \
 --zMax 16

#Make 500bP binned genome
genome_size_file=regions/dm6.genome
bin_size=500

output_file="regions/binned_${bin_size}.txt"

while IFS=$'\t' read -r chromosome size; do
    echo "Processing chromosome $chromosome"
    num_bins=$(( $size / $bin_size ))
    echo "There are $num_bins bins since the chromosome size is $size and with bin size $bin_size"
    for ((i=1; i<=$num_bins; i++)); do
        bin_start=$(( ($i - 1) * $bin_size + 1 ))
        bin_end=$(( $i * $bin_size ))
        echo -e "${chromosome}\t${bin_start}\t${bin_end}\t${chromosome}_${i}" >> ${output_file}
    done
done < $genome_size_file

#Subset bam files for equal number of reads

#Get Fru-Unbound Bins and Fru-Bound bins
bedtools intersect -a ${output_file} -b data/fruC-myc_peaks/fruC-myc-PEAKS.noCHR.bed > "regions/binned_${bin_size}_fru_bound.txt"
bedtools intersect -a ${output_file} -b data/fruC-myc_peaks/fruC-myc-PEAKS.noCHR.bed -v > "regions/binned_${bin_size}_fru_unbound.txt"

bamstrings="bam/H3K4me3_rep1_6044.bam bam/H3K4me3_rep2_6044.bam bam/H3K4me3-fruDF_rep1_6044.bam bam/H3K4me3-fruDF_rep2_6044.bam bam/H3K27ac_rep1_5127.bam bam/H3K27ac_rep2_5127.bam bam/H3K27ac-fruDF_rep1_5127.bam bam/H3K27ac-fruDF_rep2_5127.bam bam/H3K27me3_rep1_5127.bam bam/H3K27me3_rep2_5127.bam bam/H3K27me3-fruDF_rep1_5127.bam bam/H3K27me3-fruDF_rep2_5127.bam "
mkdir sub10Mbams

for bam in $bamstrings;
    do
    echo $bam
    output=$(basename "$bam")
    echo "sub10Mbams/${output}"

    samtools view -H $bam > header.sam

    samtools view $bam > reads.sam
    shuf -n 10000000 reads.sam > reads_10M.sam
    cat header.sam reads_10M.sam | samtools view -bS | samtools sort -o "sub10Mbams/${output}"
    samtools index "sub10Mbams/${output}"

    rm reads.sam
    rm reads_10M.sam
    rm header.sam
done;

#Get reads / bin -> Visualize data in python file (python-scripts.ipynb), Figure 5F

mkdir hists

bedtools multicov -bams sub10Mbams/*
-bed regions/binned_${bin_size}_fru_unbound.txt > "hists/multi_cov_subsample_fru_unbound.txt"

bedtools multicov -bams sub10Mbams/*
-bed "regions/binned_${bin_size}_fru_bound.txt" > "hists/multi_cov_subsample_fru_bound.txt"

#Call Polycomb domains using parameters from: https://www.pnas.org/doi/10.1073/pnas.1716299115#supplementary-materials

mkdir regions/sicer

sicer -t bam/H3K27me3_rep1_5127.bam -s dm6 -w 500 -f 0 -egf 0.7 -g 2000 -cpu 5 -o regions/sicer
sicer -t bam/H3K27me3_rep2_5127.bam -s dm6 -w 500 -f 0 -egf 0.7 -g 2000 -cpu 5 -o regions/sicer

awk '($4) > 500' regions/sicer/H3K27me3_rep1_5127-W500-G2000.bed | awk '($3-$2) > 3000' | bedtools merge -i stdin -d 10000 > regions/sicer/H3K27me3_rep1_5127-W500-G2000.top.bed
awk '($4) > 500' regions/sicer/H3K27me3_rep2_5127-W500-G2000.bed | awk '($3-$2) > 3000' | bedtools merge -i stdin -d 10000 > regions/sicer/H3K27me3_rep2_5127-W500-G2000.top.bed
bedtools intersect -a regions/sicer/H3K27me3_rep1_5127-W500-G2000.top.bed -b regions/sicer/H3K27me3_rep2_5127-W500-G2000.top.bed -wa > regions/sicer/polycomb.bed

#Get regions either bound by fru-peaks in non Pc domains or Pc domains 
bedtools intersect -a data/fruC-myc_peaks/fruC-myc-PEAKS.noCHR.bed -b regions/sicer/polycomb.bed -v > regions/sicer/fru_no_poly.bed

#Convert to SAF files, calculate reads/basepair in R (see R-scripts.rmd) using featureCounts, Figure 5G
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' regions/sicer/polycomb.bed > regions/sicer/polycomb.saf
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' regions/sicer/fru_no_poly.bed > regions/sicer/fru_no_poly.saf