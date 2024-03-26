rule fastqcUntrimmedTumor:
	input:
		"SnakeWES/data/fastq/{sample}.tumor.{strand}.fastq.gz"
	output:
		html="SnakeWES/qc/{sample}_tumor_{strand}_fastqc.html",
		zip="SnakeWES/qc/{sample}_tumor_{strand}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra = "--quiet"
	log:
		"SnakeWES/logs/{sample}_{strand}.fastqcUntrimmedTumor.log"
	threads: 1
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"

rule fastqcUntrimmedControl:
	input:
		"SnakeWES/data/fastq/{sample}.control.{strand}.fastq.gz"
	output:
		html="SnakeWES/qc/{sample}_control_{strand}_fastqc.html",
		zip="SnakeWES/qc/{sample}_control_{strand}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra = "--quiet"
	log:
		"SnakeWES/logs/{sample}_{strand}.fastqcUntrimmedControl.log"
	threads: 1
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"

rule fastqcUntrimmedGermline:
	input:
		"SnakeWES/data/fastq/{sample}.germline.{strand}.fastq.gz"
	output:
		html="SnakeWES/qc/{sample}_germline_{strand}_fastqc.html",
		zip="SnakeWES/qc/{sample}_germline_{strand}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra = "--quiet"
	log:
		"SnakeWES/logs/{sample}_{strand}.fastqcUntrimmedGermline.log"
	threads: 1
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"

rule fastqcTrimmedTumor:
	input:
		"SnakeWES/data/fastq/{sample}.tumor.{strand}.tr.fastq.gz"
	output:
		html="SnakeWES/qc/{sample}_tumor_{strand}_tr_fastqc.html",
		zip="SnakeWES/qc/{sample}_tumor_{strand}_tr_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra="--quiet"
	log:
		"SnakeWES/logs/{sample}_{strand}.fastqcUntrimmedTumor.log"
	threads: 1
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"

rule fastqcTrimmedControl:
	input:
		"SnakeWES/data/fastq/{sample}.control.{strand}.tr.fastq.gz"
	output:
		html="SnakeWES/qc/{sample}_control_{strand}_tr_fastqc.html",
		zip="SnakeWES/qc/{sample}_control_{strand}_tr_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra="--quiet"
	log:
		"SnakeWES/logs/{sample}_{strand}.fastqcTrimmedControl.log"
	threads: 1
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"

rule fastqcTrimmedGermline:
	input:
		"SnakeWES/data/fastq/{sample}.germline.{strand}.tr.fastq.gz"
	output:
		html="SnakeWES/qc/{sample}_germline_{strand}_tr_fastqc.html",
		zip="SnakeWES/qc/{sample}_germline_{strand}_tr_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra="--quiet"
	log:
		"SnakeWES/logs/{sample}_{strand}.fastqcTrimmedGermline.log"
	threads: 1
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"