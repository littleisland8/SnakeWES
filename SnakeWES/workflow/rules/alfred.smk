rule AlfredTumorStats:
	input:
		bam="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		bai="SnakeWES/alignments/{sample}.tumor.dd.rec.bai"
	output:
		"SnakeWES/qc/{sample}.tumor.qc.tsv.gz"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.tumor.alfred.log"
	conda:
		"../envs/alfred.yaml"
	params:
		ref=config["genome"],
		interval=config["intervals"]	
	shell:
		"alfred qc -r {params.ref} -b {params.interval} -o {output} {input.bam} 2>{log}"

rule AlfredControlStats:
	input:
		bam="SnakeWES/alignments/{sample}.control.dd.rec.bam",
		bai="SnakeWES/alignments/{sample}.control.dd.rec.bai"
	output:
		"SnakeWES/qc/{sample}.control.qc.tsv.gz"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.control.alfred.log"
	conda:
		"../envs/alfred.yaml"
	params:
		ref=config["genome"],
		interval=config["intervals"]	
	shell:
		"alfred qc -r {params.ref} -b {params.interval} -o {output} {input.bam} 2>{log}"

rule AlfredGermlineStats:
	input:
		bam="SnakeWES/alignments/{sample}.germline.dd.rec.bam",
		bai="SnakeWES/alignments/{sample}.germline.dd.rec.bai"
	output:
		"SnakeWES/qc/{sample}.germline.qc.tsv.gz"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.germline.alfred.log"
	conda:
		"../envs/alfred.yaml"
	params:
		ref=config["genome"],
		interval=config["intervals"]	
	shell:
		"alfred qc -r {params.ref} -b {params.interval} -o {output} {input.bam} 2>{log}"

rule PlotAlfredTumorStats:
	input:
		"SnakeWES/qc/{sample}.tumor.qc.tsv.gz"
	output:
		"SnakeWES/qc/{sample}.tumor.qc.tsv.gz.pdf"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.tumor.alfred.plot.log"
	conda:
		"../envs/alfred.yaml"
	params:
		script="SnakeWES/workflow/scripts/stats.R"
	shell:
		"Rscript {params.script} {input} 2>{log}"

rule PlotAlfredControlStats:
	input:
		"SnakeWES/qc/{sample}.control.qc.tsv.gz"
	output:
		"SnakeWES/qc/{sample}.control.qc.tsv.gz.pdf"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.control.alfred.plot.log"
	conda:
		"../envs/alfred.yaml"
	params:
		script="SnakeWES/workflow/scripts/stats.R"
	shell:
		"Rscript {params.script} {input} 2>{log}"

rule PlotAlfredGermlineStats:
	input:
		"SnakeWES/qc/{sample}.germline.qc.tsv.gz"
	output:
		"SnakeWES/qc/{sample}.germline.qc.tsv.gz.pdf"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.germline.alfred.plot.log"
	conda:
		"../envs/alfred.yaml"
	params:
		script="SnakeWES/workflow/scripts/stats.R"
	shell:
		"Rscript {params.script} {input} 2>{log}"