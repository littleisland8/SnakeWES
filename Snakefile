configfile: "SnakeWES/config/conf.yaml"

include: "SnakeWES/workflow/rules/fastqc.smk"
include: "SnakeWES/workflow/rules/trim.smk"
include: "SnakeWES/workflow/rules/depth.smk"
include: "SnakeWES/workflow/rules/align.smk"
include: "SnakeWES/workflow/rules/alfred.smk"
include: "SnakeWES/workflow/rules/vep.smk"

if config["run_mode"] == "PoN":

	include: "SnakeWES/workflow/rules/PoN.smk"
	include: "SnakeWES/workflow/rules/pureCN.smk"

	rule PoN:
		input:
			expand(f"qc/{{sample}}_tumor_{{strand}}_fastqc.html", sample=config["samples"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}_control_{{strand}}_fastqc.html", sample=config["controls"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}_tumor_{{strand}}_tr_fastqc.html", sample=config["samples"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}_control_{{strand}}_tr_fastqc.html", sample=config["controls"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}.tumor.qc.tsv.gz.pdf", sample=config["samples"].values()),
			expand(f"qc/{{sample}}.control.qc.tsv.gz.pdf", sample=config["controls"].values()),
			expand(f"results/{{sample}}_tumor/{{sample}}.{{caller}}.filtOnCtrs.vep.tsv", sample=config["samples"].values(), caller=config["callers"].values()),
			expand(f"results/multisample.{{caller}}.vep.tsv", caller=config["callers"].values()),


	rule PoNOneSample:
		input:
			expand(f"qc/{{sample}}_tumor_{{strand}}_fastqc.html", sample=config["samples"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}_control_{{strand}}_fastqc.html", sample=config["controls"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}_tumor_{{strand}}_tr_fastqc.html", sample=config["samples"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}_control_{{strand}}_tr_fastqc.html", sample=config["controls"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}.tumor.qc.tsv.gz.pdf", sample=config["samples"].values()),
			expand(f"qc/{{sample}}.control.qc.tsv.gz.pdf", sample=config["controls"].values()),
			expand(f"results/{{sample}}_tumor/{{sample}}.{{caller}}.filtOnCtrs.vep.tsv", sample=config["samples"].values(), caller=config["callers"].values()),

	rule PureCN:
		input:
			expand(f"alignments/{{sample}}.tumor.coveragePureCN.loess.txt", sample=config["samples"].values()),
			expand(f"alignments/{{sample}}.control.coveragePureCN.loess.txt", sample=config["controls"].values()),
			expand(f"results/{{sample}}_tumor/{{sample}}.{{caller}}.filtOnCtr.vep.tmp01.tsv", sample=config["samples"].values(), caller=config["callers"].values()),
			expand(f"results/multisample.{{caller}}.vep.tmp01.tsv", caller=config["callers"].values()) 

elif config["run_mode"] == "Nocontrol":

	include: "SnakeWES/workflow/rules/nocontrol.smk"
	
	rule Nocontrols:
		input:
			expand(f"qc/{{sample}}_tumor_{{strand}}_fastqc.html", sample=config["samples"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}_tumor_{{strand}}_tr_fastqc.html", sample=config["samples"].values(), strand=config["strand"].values()),
			expand(f"alignments/{{sample}}.tumor.dd.rec.bam", sample=config["samples"].values()),
			expand(f"qc/{{sample}}.tumor.qc.tsv.gz.pdf", sample=config["samples"].values()),
			expand(f"results/multisample.{{caller}}.vcf.gz", caller=config["callers"].values()),
			expand(f"results/{{sample}}_tumor/{{sample}}.{{caller}}.filt.vep.vcf.gz", sample=config["samples"].values(), caller=config["callers"].values()),
			expand(f"results/multisample.{{caller}}.vep.vcf.gz", caller=config["callers"].values()),
			expand(f"results/{{sample}}_tumor/{{sample}}.{{caller}}.filt.vep.tsv", sample=config["samples"].values(), caller=config["callers"].values()),
			expand(f"results/multisample.{{caller}}.vep.tsv", caller=config["callers"].values()),


	rule NocontrolsOneSample:
		input:
			expand(f"qc/{{sample}}_tumor_{{strand}}_fastqc.html", sample=config["samples"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}_tumor_{{strand}}_tr_fastqc.html", sample=config["samples"].values(), strand=config["strand"].values()),
			expand(f"alignments/{{sample}}.tumor.dd.rec.bam", sample=config["samples"].values()),
			expand(f"qc/{{sample}}.tumor.qc.tsv.gz.pdf", sample=config["samples"].values()),
			expand(f"results/{{sample}}_tumor/{{sample}}.{{caller}}.vcf.gz", sample=config["samples"].values(), caller=config["callers"].values())

elif config["run_mode"] == "Paired":

	include: "SnakeWES/workflow/rules/paired.smk"

	rule Paired:
		input:
			expand(f"qc/{{sample}}_{{type}}_{{strand}}_fastqc.html", sample=config["samples"].values(), type=config["paired"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}_{{type}}_{{strand}}_tr_fastqc.html", sample=config["samples"].values(), type=config["paired"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}.{{type}}.qc.tsv.gz.pdf", sample=config["samples"].values(),type=config["paired"].values()),
			expand(f"results/{{sample}}.varscan.paired.norm.vcf.gz", sample=config["samples"].values()),
			expand(f"results/{{sample}}.mutect2.paired.filtered.norm.vcf.gz", sample=config["samples"].values()),
			f"results/multisample.mutect2.paired.tmp01.tsv",
			f"results/multisample.varscan.paired.tmp01.tsv"

	rule PairedOneSample:
		input:
			expand(f"qc/{{sample}}_{{type}}_{{strand}}_fastqc.html", sample=config["samples"].values(), type=config["paired"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}_{{type}}_{{strand}}_tr_fastqc.html", sample=config["samples"].values(), type=config["paired"].values(), strand=config["strand"].values()),
			expand(f"qc/{{sample}}.{{type}}.qc.tsv.gz.pdf", sample=config["samples"].values(),type=config["paired"].values()),
			expand(f"results/{{sample}}.varscan.paired.snp.vcf.gz.tbi",sample=config["samples"].values()),
			expand(f"results/{{sample}}.varscan.paired.indel.vcf.gz.tbi",sample=config["samples"].values()),
			expand(f"results/{{sample}}.varscan.paired.vcf.gz", sample=config["samples"].values()),
			expand(f"results/{{sample}}.mutect2.paired.filtered.vep.vcf.gz.tbi", sample=config["samples"].values()),
			expand(f"results/{{sample}}.varscan.paired.vep.vcf.gz.tbi", sample=config["samples"].values()),


else:

	print("[Error] invalid choice: choose one from 'Nocontrol', 'PoN', 'Paired'")



#rule fastqc:
#	input:
#		expand(f"qc/{{sample}}_tumor_{{strand}}_fastqc.html", sample=config["samples"].values(), strand=config["strand"].values()),
#		expand(f"qc/{{sample}}_control_{{strand}}_fastqc.html", sample=config["controls"].values(), strand=config["strand"].values()),
#		expand(f"qc/{{sample}}_tumor_{{strand}}_tr_fastqc.html", sample=config["samples"].values(), strand=config["strand"].values()),
#		expand(f"qc/{{sample}}_control_{{strand}}_tr_fastqc.html", sample=config["controls"].values(), strand=config["strand"].values()),

#rule align:
#	input:
#		expand(f"alignments/{{sample}}.tumor.dd.rec.bam", sample=config["samples"].values()),
#		expand(f"alignments/{{sample}}.control.dd.rec.bam", sample=config["controls"].values())

#rule bamqc:
#	input:
#		expand(f"qc/{{sample}}.tumor.qc.tsv.gz.pdf", sample=config["samples"].values()),
#		expand(f"qc/{{sample}}.control.qc.tsv.gz.pdf", sample=config["controls"].values())

	
#rule callers:
#	input:
#		expand(f"results/{{sample}}_tumor/{{sample}}.{{caller}}.filtOnCtr.vcf.gz", sample=config["samples"].values(), caller=config["callers"].values()),
#		expand(f"results/{{sample}}_control/{{sample}}.{{caller}}.filt.vcf.gz", sample=config["controls"].values(), caller=config["callers"].values()),
#		expand(f"results/multisample.{{caller}}.vcf.gz", caller=config["callers"].values())

#rule pureCN:
#	input:
#		expand(f"alignments/{{sample}}.tumor.coveragePureCN.loess.txt", sample=config["samples"].values()),
#		expand(f"alignments/{{sample}}.control.coveragePureCN.loess.txt", sample=config["controls"].values())

#rule vep:
#	input:
#		expand(f"results/{{sample}}_tumor/{{sample}}.{{caller}}.filtOnCtr.vep.tmp01.tsv", sample=config["samples"].values(), caller=config["callers"].values()),
#		expand(f"results/multisample.{{caller}}.vep.tmp01.tsv", caller=config["callers"].values()) 
		#expand(f"results/{{sample}}_tumor/{{sample}}.{{caller}}.filtOnCtr.vep.vcf.gz", sample=config["samples"].values(), caller=config["callers"].values()),
