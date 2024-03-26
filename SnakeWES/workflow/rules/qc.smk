rule fastqc:
    input:
        "SnakeWES/data/fastq/{sample}.R1.tr.fastq.gz"
    output:
        html="SnakeWES/qc/fastqc/{sample}.R1.tr.html",
        zip="SnakeWES/qc/fastqc/{sample}.R1.tr.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "SnakeWES/logs/fastqc/{sample}.log"
    threads: 1
    resources:
        mem_mb=4096 
    wrapper:
        "v3.3.3/bio/fastqc"


rule samtools_stats:
    input:
        bam="SnakeWES/data/bam/{sample}.dd.rec.bam",
        bed=config['interval']
    output:
        "SnakeWES/qc/samtools_stats/{sample}.txt"
    log:
        "SnakeWES/logs/stats/{sample}.samtools.log",
    wrapper:
        "v3.3.3/bio/samtools/stats"
