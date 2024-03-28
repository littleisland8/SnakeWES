rule bwaTumor:
	input:
		reads=["SnakeWES/data/fastq/{sample}.tumor.R1.tr.fastq.gz", "SnakeWES/data/fastq/{sample}.tumor.R2.tr.fastq.gz"],
		idx=multiext(config['genome'], ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
	output:
		temp("SnakeWES/alignments/{sample}.tumor.srt.bam")
	log:
		"SnakeWES/logs/{sample}.bwaTumor.log"
	message:
		"Mapping reads with bwa-mem2 on tumor {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.bwaTumor.txt"
	params:
		extra=r"-R '@RG\tID:{sample}\tSM:{sample}_tumor'",
		sort="samtools",  
		sort_order="coordinate"
	threads: 1
	wrapper:
		"v3.3.3/bio/bwa-mem2/mem"

rule bwaControl:    
	input: 
		reads=["SnakeWES/data/fastq/{sample}.control.R1.tr.fastq.gz", "SnakeWES/data/fastq/{sample}.control.R2.tr.fastq.gz"],
		idx=multiext(config['genome'], ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
	output: 
		temp("SnakeWES/alignments/{sample}.control.srt.bam")
	threads: 1
	message:
		"Mapping reads with bwa-mem2 on control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.bwaControl.txt"
	params:
		extra=r"-R '@RG\tID:{sample}\tSM:{sample}_control'",
		sort="samtools",  
		sort_order="coordinate"
	log:
		"SnakeWES/logs/{sample}.bwaControl.log"
	wrapper:
		"v3.3.3/bio/bwa-mem2/mem"

rule bwaGermline:    
	input: 
		reads=["SnakeWES/data/fastq/{sample}.germline.R1.tr.fastq.gz", "SnakeWES/data/fastq/{sample}.germline.R2.tr.fastq.gz"],
		idx=multiext(config['genome'], ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
	output: 
		temp("SnakeWES/alignments/{sample}.germline.srt.bam")
	threads: 1
	message:
		"Mapping reads with bwa-mem2 on germline {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.bwaGermline.txt"
	params:
		extra=r"-R '@RG\tID:{sample}\tSM:{sample}_germline'",
		sort="samtools",  
		sort_order="coordinate"
	log:
		"SnakeWES/logs/{sample}.bwaGermline.log"
	wrapper:
		"v3.3.3/bio/bwa-mem2/mem"

rule markDuplicatesTumor:
	input:
		bams="SnakeWES/alignments/{sample}.tumor.srt.bam"
	output:
		bam=temp("SnakeWES/alignments/{sample}.tumor.dd.bam"),
		bai=temp("SnakeWES/alignments/{sample}.tumor.dd.bai"),
		metrics="SnakeWES/data/{sample}.tumor.metrics.txt"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.markDuplicatesTumor.log"
	message:
		"Mark and remove PCR-optical duplicates"
	benchmark:
		"benchmarks/{sample}.markDuplicatesTumor.txt"
	params:
		extra="--REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT"
	resources:
		mem_mb=4096    
	log:
		"SnakeWES/logs/{sample}.markDuplicatesTumor.log"
	wrapper:
		 "v3.3.3/bio/picard/markduplicates"

rule markDuplicatesControl:
	input:
		bams="SnakeWES/alignments/{sample}.control.srt.bam"
	output:
		bam=temp("SnakeWES/alignments/{sample}.control.dd.bam"),
		bai=temp("SnakeWES/alignments/{sample}.control.dd.bai"),
		metrics="SnakeWES/data/{sample}.control.metrics.txt"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.markDuplicatesControl.log"
	message:
		"Mark and remove PCR-optical duplicates"
	benchmark:
		"benchmarks/{sample}.markDuplicatesControl.txt"
	params:
		extra="--REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT"
	resources:
		mem_mb=4096    
	log:
		"SnakeWES/logs/{sample}.markDuplicatesControl.log"
	wrapper:
		 "v3.3.3/bio/picard/markduplicates"

rule markDuplicatesGermline:
	input:
		bams="SnakeWES/alignments/{sample}.germline.srt.bam"
	output:
		bam=temp("SnakeWES/alignments/{sample}.germline.dd.bam"),
		bai=temp("SnakeWES/alignments/{sample}.germline.dd.bai"),
		metrics="SnakeWES/data/{sample}.germline.metrics.txt"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.markDuplicatesGermline.log"
	message:
		"Mark and remove PCR-optical duplicates"
	benchmark:
		"benchmarks/{sample}.markDuplicatesGermline.txt"
	params:
		extra="--REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT"
	resources:
		mem_mb=4096    
	log:
		"SnakeWES/logs/{sample}.markDuplicatesGermline.log"
	wrapper:
		 "v3.3.3/bio/picard/markduplicates"

rule baseRecalibratorTumor:
	input:
		bam="SnakeWES/alignments/{sample}.tumor.dd.bam",
		bai="SnakeWES/alignments/{sample}.tumor.dd.bai",
		ref=config['genome'],
		dict="SnakeWES/resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		known=["SnakeWES/resources/dbSNP.b156.filt.vcf.gz", 
		"SnakeWES/resources/Homo_sapiens_assembly38.known_indels.vcf.gz",
		config["gnomAD"], "SnakeWES/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"],  # optional known sites - single or a list
	output:
		recal_table="SnakeWES/data/{sample}.tumor.recal.table",
	threads: 1
	benchmark:
		"SnakeWES/benchmarks/{sample}.baseRecalibratorTumor.txt"
	log:
		"SnakeWES/logs/{sample}.baseRecalibratorTumor.log"
	message:
		"Perform base recalibration model"
	params:
		extra="-L " + config['intervals']
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/baserecalibrator"

rule baseRecalibratorControl:
	input:
		bam="SnakeWES/alignments/{sample}.control.dd.bam",
		bai="SnakeWES/alignments/{sample}.control.dd.bai",
		ref=config['genome'],
		dict="SnakeWES/resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		known=["SnakeWES/resources/dbSNP.b156.filt.vcf.gz", 
		"SnakeWES/resources/Homo_sapiens_assembly38.known_indels.vcf.gz",
		config["gnomAD"], "SnakeWES/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"],  # optional known sites - single or a list
	output:
		recal_table="SnakeWES/data/{sample}.control.recal.table"
	threads: 1
	benchmark:
		"benchmarks/{sample}.baseRecalibratorControl.txt"
	log:
		"SnakeWES/logs/{sample}.baseRecalibratorControl.log"
	message:
		"Perform base recalibration model"
	params:
		extra="-L " + config['intervals']
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/baserecalibrator"

rule baseRecalibratorGermline:
	input:
		bam="SnakeWES/alignments/{sample}.germline.dd.bam",
		bai="SnakeWES/alignments/{sample}.germline.dd.bai",
		ref=config['genome'],
		dict="SnakeWES/resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		known=["SnakeWES/resources/dbSNP.b156.filt.vcf.gz", 
		"SnakeWES/resources/Homo_sapiens_assembly38.known_indels.vcf.gz",
		config["gnomAD"], "SnakeWES/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"],  # optional known sites - single or a list
	output:
		recal_table="SnakeWES/data/{sample}.germline.recal.table"
	threads: 1
	benchmark:
		"SnakeWES/benchmarks/{sample}.baseRecalibratorGermline.txt"
	log:
		"SnakeWES/logs/{sample}.baseRecalibratorGermline.log"
	message:
		"Perform base recalibration model"
	params:
		extra="-L " + config['intervals']		
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/baserecalibrator"

rule ApplyBQSRTumor:
	input:
		bam="SnakeWES/alignments/{sample}.tumor.dd.bam",
		bai="SnakeWES/alignments/{sample}.tumor.dd.bai",	
		ref=config['genome'],
		dict="SnakeWES/resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		recal_table="SnakeWES/data/{sample}.tumor.recal.table",
	output:
		bam="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		bai="SnakeWES/alignments/{sample}.tumor.dd.rec.bai"
	threads: 1
	benchmark:
		"benchmarks/{sample}.ApplyBQSRTumor.txt"
	log:
		"SnakeWES/logs/recal/{sample}.ApplyBQSRTumor.log",
	params:
		extra="-L " + config['intervals']
		#embed_ref=True,  # embed the reference in cram output
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/applybqsr"

rule ApplyBQSRControl:
	input:
		bam="SnakeWES/alignments/{sample}.control.dd.bam",
		bai="SnakeWES/alignments/{sample}.control.dd.bai",   
		ref=config['genome'],
		dict="SnakeWES/resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		recal_table="SnakeWES/data/{sample}.control.recal.table",
	output:
		bam="SnakeWES/alignments/{sample}.control.dd.rec.bam",
		bai="SnakeWES/alignments/{sample}.control.dd.rec.bai"
	threads: 1
	benchmark:
		"SnakeWES/benchmarks/{sample}.ApplyBQSRControl.txt"
	log:
		"SnakeWES/logs/recal/{sample}.ApplyBQSRControl.log",
	params:
		extra="-L " + config['intervals']
		#embed_ref=True,  # embed the reference in cram output
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/applybqsr"

rule ApplyBQSRGermline:
	input:
		bam="SnakeWES/alignments/{sample}.germline.dd.bam",
		bai="SnakeWES/alignments/{sample}.germline.dd.bai",   
		ref=config['genome'],
		dict="resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		recal_table="SnakeWES/data/{sample}.germline.recal.table",
	output:
		bam="SnakeWES/alignments/{sample}.germline.dd.rec.bam",
		bai="SnakeWES/alignments/{sample}.germline.dd.rec.bai"
	threads: 1
	benchmark:
		"SnakeWES/benchmarks/{sample}.ApplyBQSRGermline.txt"
	log:
		"SnakeWES/logs/recal/{sample}.ApplyBQSRGermline.log",
	params:
		extra="-L " + config['intervals']
		#embed_ref=True,  # embed the reference in cram output
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/applybqsr"
