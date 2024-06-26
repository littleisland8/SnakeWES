rule fastpTumor: ## aggiungere flag per i threads
	input:
		sample=["SnakeWES/data/fastq/{sample}.tumor.R1.fastq.gz", "SnakeWES/data/fastq/{sample}.tumor.R2.fastq.gz"]
	output:
		trimmed=["SnakeWES/data/fastq/{sample}.tumor.R1.tr.fastq.gz", "SnakeWES/data/fastq/{sample}.tumor.R2.tr.fastq.gz"],
		json="SnakeWES/data/{sample}.tumor.json",
		failed="SnakeWES/data/{sample}.tumor.failedreads.txt",
		html="SnakeWES/data/{sample}.tumor.html",
		unpaired1="SnakeWES/data/{sample}.tumor.u1.fq.gz",
		unpaired2="SnakeWES/data/{sample}.tumor.u2.fq.gz"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.fastpTumor.log"
	message:
		"Trimming with fastp"	
	params:
		#adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
		extra="--length_required 30 --detect_adapter_for_pe --disable_quality_filtering"
	benchmark:
		"SnakeWES/benchmarks/{sample}.fastpTumor.txt"
	wrapper:
		"v3.3.3/bio/fastp"

rule fastpControl: ## aggiungere flag per i threads
	input:
		sample=["SnakeWES/data/fastq/{sample}.control.R1.fastq.gz", "SnakeWES/data/fastq/{sample}.control.R2.fastq.gz"]
	output:
		trimmed=["SnakeWES/data/fastq/{sample}.control.R1.tr.fastq.gz", "SnakeWES/data/fastq/{sample}.control.R2.tr.fastq.gz"],
		json="SnakeWES/data/{sample}.control.json",
		failed="SnakeWES/data/{sample}.control.failedreads.txt",
		html="SnakeWES/data/{sample}.control.html",
		unpaired1="SnakeWES/data/{sample}.control.u1.fq.gz",
		unpaired2="SnakeWES/data/{sample}.control.u2.fq.gz"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.fastpControl.log"
	message:
		"Trimming with fastp"	
	params:
		#adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
		extra="--length_required 30 --detect_adapter_for_pe --disable_quality_filtering"
	benchmark:
		"SnakeWES/benchmarks/{sample}.fastpControl.txt"
	wrapper:
		"v3.3.3/bio/fastp"

rule fastpGermline: ## aggiungere flag per i threads
	input:
		sample=["SnakeWES/data/fastq/{sample}.germline.R1.fastq.gz", "SnakeWES/data/fastq/{sample}.germline.R2.fastq.gz"]
	output:
		trimmed=["SnakeWES/data/fastq/{sample}.germline.R1.tr.fastq.gz", "SnakeWES/data/fastq/{sample}.germline.R2.tr.fastq.gz"],
		json="SnakeWES/data/{sample}.germline.json",
		failed="SnakeWES/data/{sample}.germline.failedreads.txt",
		html="SnakeWES/data/{sample}.germline.html",
		unpaired1="SnakeWES/data/{sample}.germline.u1.fq.gz",
		unpaired2="SnakeWES/data/{sample}.germline.u2.fq.gz"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.fastpGermline.log"
	message:
		"Trimming with fastp"	
	params:
		#adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
		extra="--length_required 30 --detect_adapter_for_pe --disable_quality_filtering"
	benchmark:
		"SnakeWES/benchmarks/{sample}.fastpGermline.txt"
	wrapper:
		"v3.3.3/bio/fastp"