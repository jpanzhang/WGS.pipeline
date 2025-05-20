
# WGS.pipeline
The main steps, software, and code/parameters used from off-machine data to the VCF file are described below. 

#### (1) Quality control of raw reads
#### Software: fastp v0.23.4
Code: fastp_install_dir/fastp -w 8 -i sample.R1.fq -I sample.R2.fq -o sample.R1.clean.fq -O sample.R2.clean.fq 

#### (2) Mapping reads and sorting the BAM file
#### Software: Sentieon v202308
Code: Sentieon_install_dir/bin/sentieon bwa mem -R "@RG\tID:rg_sample\tSM:sample\tPL:$PL" -t 40 -K 10000000 sheep.genome. fa sample.R1.clean.fq sample.R2.clean.fq | util sort -t 40 --sam2bam -o sample.sorted.bam -i -

#### (3) Mark Duplicates
#### Software: Sentieon v202308
Code: Sentieon_install_dir/bin/sentieon driver -t 40 -i sample-sorted.bam --algo LocusCollector --fun score_info sample.score.txt 
Sentieon_install_dir/bin/sentieon driver -t 40 -i sample-sorted.bam --algo Dedup --score_info sample.score.txt --metrics sample.dedup_metrics.txt sample.deduped.bam 

#### (4) Sequencing accuracy assessment and alignment statistics (optional)
#### Software: Sentieon v202308, PanDepth v2.25, and Samtools v1.13
#### Code for calculating reads quality: 
Sentieon_install_dir/bin/sentieon driver -r sheep.genome.fa -t 40 -i sample-sorted.bam \
--algo MeanQualityByCycle sample.mq_metrics.txt \
--algo QualDistribution sample.qd_metrics.txt \
--algo GCBias --summary sample.gc_summary.txt sample.gc_metrics.txt \
--algo AlignmentStat sample.aln_metrics.txt \
--algo BaseDistributionByCycle sample.bd_metrics.txt \
--algo QualityYield sample.qy_metrics.txt \
--algo InsertSizeMetricAlgo sample.is_metrics.txt) 
#### Code for calculating genome coverage: 
PanDepth_install_dir/pandepth -i sample.deduped.bam -t 8 -o genome.coverage.stat.txt
#### Code for calculating mapping rate: 
Samtools_install_dir/Samtools flagstat -@40 -i sample.deduped.bam > sample.mapped.stat.txt

#### (5) Variants calling for each sample
#### Software: Sentieon v202308
Code: Sentieon_install_dir/bin/sentieon driver -r sheep.genome.fa -t 40 -i sample.deduped.bam --algo Haplotyper --emit_conf=30 --call_conf=30 --emit_mode gvcf sample.gvcf.vcf.gz 

#### (6) Variant joint calling for all samples
#### Software: Sentieon v202308
Code: GVCF_list=gvcf.gz.list.txt #paths to all gVCF files
GVCF_inputs=$(awk '{info=info" -v "$0} END {print info}' $GVCF_list) 
Sentieon_install_dir/bin/sentieon driver -t 40 -r sheep.genome.fa --algo GVCFtyper $GVCF_inputs joint-calling.vcf 

#### (7) Variants extraction and hard filtration
#### Software: GATK v4.1.8.1
#### Code for SNPs: 
GATK_install_dir/gatk SelectVariants -R sheep.fa -V joint-calling.vcf -O filtered.SNP.vcf \
-select '((QD>=2.0 && MQ>=40.0 && FS<=60.0 && SOR<=3.0) && (QD>=2.0 && MQ<40.0 && FS<=60.0 && SOR<=3.0 || MQRankSum>=-12.5 && ReadPosRankSum>=-8.0) && vc.isSNP())'
#### Code for InDels: 
GATK_install_dir/gatk SelectVariants -R sheep.fa -V joint-calling.vcf -O filtered.InDel.vcf \
-select '((QD>=2.0 && FS<=200.0 && SOR<=10.0) && (QD>=2.0 && FS<=200.0 && SOR<=10.0 || ReadPosRankSum>=-20.0) && vc.isIndel())'

#### (8) Further filtration
#### Software: VCFtools v0.1.16
#### Code for SNPs: 
Vcftools_install_dir/vcftools --vcf filtered.SNP.vcf --max-alleles 2 --min-alleles 2 --min-meanDP 5 --max-missing 0.9 --remove-filtered-all --recode --recode-INFO-all --out SNP.vcf
#### Code for InDels: 
Vcftools_install_dir/vcftools --vcf filtered.InDel.vcf --min-meanDP 5 --max-missing 0.9 --remove-filtered-all --recode --recode-INFO-all --out InDel.vcf

