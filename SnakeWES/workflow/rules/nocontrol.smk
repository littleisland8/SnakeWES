########################################################################### Mutect2 ###########################################################################

rule Mutect2TumorNocontrol: 
	input:
		fasta=config["genome"],
		map="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		intervals=config["intervals"],
		pon=config["gatk_pon"],
		germline=config["gnomAD"]
	output:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.mutect2.vcf.gz",
		bam="SnakeWES/alignments/{sample}.tumor.mutect2.bam",        
		f1r2="SnakeWES/data/{sample}.tumor.f1r2.tar.gz"
	message:
		"Mutect2 calling with {wildcards.sample}"
	threads: 1
	resources:
		mem_mb=4096
	log:
		"SnakeWES/logs/{sample}.Mutect2TumorNocontrol.log",
	params:
		extra="--max-reads-per-alignment-start 0" 
	wrapper:
		"v3.3.3/bio/gatk/mutect"

rule OrientationModelTumorNocontrol:
	input:
		f1r2="SnakeWES/data/{sample}.tumor.f1r2.tar.gz",
	output:
		"SnakeWES/data/{sample}.tumor.artifacts_prior.tar.gz",
	resources:
		mem_mb=4096,
	log:
		"SnakeWES/logs/{sample}.OrientationModelTumorNocontrol.log",
	wrapper:
		"v3.3.3/bio/gatk/learnreadorientationmodel"

rule GetPileupSummariesTumorNocontrol:
	input:
		bam="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		intervals=config["intervals"],
		variants=config["gnomAD"]
	output:
		"SnakeWES/data/{sample}.tumor.pileupTable"
	threads: 1
	resources:
		mem_mb=4096,
	params:
		extra="",
	log:
		"SnakeWES/logs/{sample}.GetPileupSummariesTumorNocontrol.log",
	wrapper:
		"v3.3.3/bio/gatk/getpileupsummaries"

rule CalculateContaminationTumorNocontrol:
	input:
		"SnakeWES/data/{sample}.tumor.pileupTable"
	output:
		contamination="SnakeWES/data/{sample}.tumor.contaminationTable",
		segmentation="SnakeWES/data/{sample}.tumor.tumorSeg.txt"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.CalculateContaminationTumorNocontrol.log",
	conda:
		"../envs/gatk4.yaml"
	params:
		mem_mb="-Xmx4G"
		#extra="--tumor-segmentation SnakeWES/data/{wildcards.sample}.tumseg.txt"
	shell:
		"gatk --java-options {params.mem_mb} CalculateContamination -I {input} -O {output.contamination} --tumor-segmentation {output.segmentation} > {log} 2>&1"


rule FilterMutectCallsTumorNocontrol:
	input:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.mutect2.vcf.gz",
		ref=config["genome"],
		bam="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		intervals=config["intervals"],
		contamination="SnakeWES/data/{sample}.tumor.contaminationTable", # from gatk CalculateContamination
		segmentation="SnakeWES/data/{sample}.tumor.tumorSeg.txt", # from gatk CalculateContamination
		f1r2="SnakeWES/data/{sample}.tumor.artifacts_prior.tar.gz" # from gatk LearnReadOrientationBias
	output:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.mutect2.FiltMut.vcf.gz"
	threads: 10
	log:
		"SnakeWES/logs/{sample}.FilterMutectCallsTumorNocontrol.log",
	params:
		#extra="--tumor-segmentation SnakeWES/data/{wildcard.sample}.tumseg.txt",  # optional arguments, see GATK docs
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/gatk/filtermutectcalls"

rule CleanFilterMutectTumorOutputNocontrol:
	input:
		"SnakeWES/results/{sample}_tumor/{sample}.mutect2.FiltMut.vcf.gz"
	output:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.mutect2.filt.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.mutect2.filt.vcf.gz.tbi"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.CleanFilterMutectTumorOutputNocontrol.log"
	params:	
		excl=config["chr_to_exclude"],
		depth=config['filtering_tumors']['min_depth'],
		vaf=config['filtering_tumors']['vaf'],
		alt=config['filtering_tumors']['alt_depth'],
		ref=config["genome"],
	conda:
		"../envs/bcftools.yaml"
	shell:	
		"""
		bcftools view -Ov -i'FILTER == "PASS" || FILTER == "germline" & POPAF > 2' {input} |
		grep -v -f {params.excl} | 
		bcftools norm -m - -f {params.ref} | 
		bcftools view -i "FORMAT/DP[0] >= {params.depth} & FORMAT/AD[0:1] >= {params.alt} & FORMAT/AF >= {params.vaf}" -Oz -o {output.vcf} &&
		tabix -p vcf {output.vcf}
		"""

rule MergeMutect2TumorOutputNocontrol:
	input:
		vcf=expand(f"SnakeWES/results/{{sample}}_tumor/{{sample}}.mutect2.filt.vep.vcf.gz", sample=config["samples"].values()),
		tbi=expand(f"SnakeWES/results/{{sample}}_tumor/{{sample}}.mutect2.filt.vep.vcf.gz.tbi", sample=config["samples"].values())
	output:
		vcf="SnakeWES/results/multisample.mutect2.vcf.gz",
		tbi="SnakeWES/results/multisample.mutect2.vcf.gz.tbi"
	threads: 1
	log:
		"SnakeWES/logs/multisample.MergeMutect2TumorOutputNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} && "
		"tabix -p vcf {output.vcf}"


########################################################################### Bcftools ###########################################################################

rule BcftoolsCallOnTumorsNocontrol:  
	input:
		bam="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		bai="SnakeWES/alignments/{sample}.tumor.dd.rec.bai"
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.bcftools.vcf.gz"
	threads: 1
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	conda:
		"../envs/bcftools.yaml"
	log:
		"SnakeWES/logs/{sample}.BcftoolsCallOnTumorsNocontrol.log"
	shell:
		"bcftools mpileup -Ou -d 10000 -R {params.intervals} -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -f {params.ref} {input.bam} | "
		"bcftools call -mv | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2>{log}"

rule BcftoolsIndexTumorsNocontrol:
	input:
		"SnakeWES/results/{sample}_tumor/{sample}.bcftools.vcf.gz"
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.bcftools.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/tabix.yaml"
	log:
		"SnakeWES/logs/{sample}.BcftoolsIndexTumorsNocontrol.log"
	shell:
		"tabix -p vcf {input}"

rule BcftoolsAddAfFieldOnTumorsNocontrol:  
	input:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.bcftools.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.bcftools.vcf.gz.tbi"
	output:
		vcf=temp("SnakeWES/results/{sample}_tumor/{sample}.bcftools.addaf.vcf.gz")
#		tsv=temp("SnakeWES/results/{sample}_tumor/{sample}.annot.tsv.gz"),
#		tsv_tbi=temp("SnakeWES/results/{sample}_tumor/{sample}.annot.tsv.gz.tbi")
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"SnakeWES/logs/{sample}.BcftoolsAddAfFieldOnTumorsNoncontrol.log"
	params:
		tsv="SnakeWES/results/{sample}_tumor/{sample}.bcftools.annot.tsv.gz",
		tsv_tbi="SnakeWES/results/{sample}_tumor/{sample}.bcftools.annot.tsv.gz.tbi"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%AD{{0}}\t%AD{{1}}]' {input.vcf} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {params.tsv} 2>{log} && 
		tabix -b2 -e2 {params.tsv} 2>>{log} && 
		bcftools annotate -a {params.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF {input.vcf} -Oz -o {output.vcf} 2>>{log} && 
		rm {params.tsv} 2>>{log} && 
		rm {params.tsv_tbi} 2>>{log}
		"""

rule BcftoolsFilterTumorsNocontrol: 
	input:
		"SnakeWES/results/{sample}_tumor/{sample}.bcftools.addaf.vcf.gz"
	output:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.bcftools.filt.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.bcftools.filt.vcf.gz.tbi"
	threads: 1
	params:
		config["chr_to_exclude"],
		excl=config["chr_to_exclude"],
		depth=config['filtering_tumors']['min_depth'],
		vaf=config['filtering_tumors']['vaf'],
		alt=config['filtering_tumors']['alt_depth']
	conda:
		"../envs/bcftools.yaml"
	log:
		"SnakeWES/logs/{sample}.BcftoolsFilterTumorsNocontrol.log"
	shell:
		"bcftools view -i 'FORMAT/DP >= {params.depth} & FORMAT/AF >= {params.vaf} & FORMAT/AD[0:1] >= {params.alt}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2>{log} && "
		"tabix -p vcf {output.vcf} 2>>{log}"

rule MergeBcftoolsTumorOutputNocontrol:
	input:
		vcf=expand(f"SnakeWES/results/{{sample}}_tumor/{{sample}}.bcftools.filt.vcf.gz", sample=config["samples"].values()),
		tbi=expand(f"SnakeWES/results/{{sample}}_tumor/{{sample}}.bcftools.filt.vcf.gz.tbi", sample=config["samples"].values())
	output:
		vcf="SnakeWES/results/multisample.bcftools.vcf.gz",
		tbi="SnakeWES/results/multisample.bcftools.vcf.gz.tbi"
	threads: 1
	log:
		"SnakeWES/logs/multisample.MergeBcftoolsTumorOutputNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} && "
		"tabix -p vcf {output.vcf} 2>>{log}"


########################################################################### Varscan ###########################################################################

rule VarscanCallSnvTumorsNocontrol: 
	input:
		bam="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		bai="SnakeWES/alignments/{sample}.tumor.dd.rec.bai"
	output:
		vcf_snv=temp("SnakeWES/results/{sample}_tumor/{sample}.varscan.snv.vcf.gz"),
		tbi_snv=temp("SnakeWES/results/{sample}_tumor/{sample}.varscan.snv.vcf.gz.tbi")
	threads: 1
	params:
		intervals=config["intervals"],
		ref=config["genome"]
	conda:
		"../envs/varscan.yaml"
	log:
		"SnakeWES/logs/{sample}.VarscanCallSnvTumorsNocontrol.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} |"
		"varscan mpileup2snp --strand-filter 1 --output-vcf --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 |"
		"bgzip -c > {output.vcf_snv} 2>{log} &&"
		"tabix -p vcf {output.vcf_snv} 2>>{log}"

rule VarscanCallIndelTumorsNocontrol: 
	input:
		bam="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		bai="SnakeWES/alignments/{sample}.tumor.dd.rec.bai"
	output:
		vcf_indel=temp("SnakeWES/results/{sample}_tumor/{sample}.varscan.indel.vcf.gz"),
		tbi_indel=temp("SnakeWES/results/{sample}_tumor/{sample}.varscan.indel.vcf.gz.tbi")
	threads: 1
	params:
		intervals=config["intervals"],
		ref=config["genome"]
	conda:
		"../envs/varscan.yaml"
	log:
		"SnakeWES/logs/{sample}.VarscanCallIndelTumorsNocontrol.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} | "
		"varscan mpileup2indel --strand-filter 1 --output-vcf --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 | "
		"bgzip -c > {output.vcf_indel} 2>{log} && "
		"tabix -p vcf {output.vcf_indel} 2>>{log}"

rule VarscanConcatTumorNocontrol: 
	input:
		vcf_snv="SnakeWES/results/{sample}_tumor/{sample}.varscan.snv.vcf.gz",
		vcf_indel="SnakeWES/results/{sample}_tumor/{sample}.varscan.indel.vcf.gz",
		tbi_snv="SnakeWES/results/{sample}_tumor/{sample}.varscan.snv.vcf.gz.tbi",
		tbi_indel="SnakeWES/results/{sample}_tumor/{sample}.varscan.indel.vcf.gz.tbi"
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.varscan.concat.tmp.vcf.gz"
	params: 
		config["genome"]
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"SnakeWES/logs/{sample}.VarscanConcatTumorNocontrol.log"
	shell: 
		"bcftools concat -a -Ov {input.vcf_snv} {input.vcf_indel} | "
		"bcftools norm -Oz -m - -f {params} -Oz -o {output} 2>{log}"

rule VarscanModAFTumorNoncontrol:
	input:
		"SnakeWES/results/{sample}_tumor/{sample}.varscan.concat.tmp.vcf.gz"
	output:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.varscan.concat.vcf.gz",
		tsv="SnakeWES/results/{sample}_tumor/{sample}.varscan.annot.tsv.gz",
		tsv_tbi="SnakeWES/results/{sample}_tumor/{sample}.varscan.annot.tsv.gz.tbi"
	log:
		"SnakeWES/logs/{sample}.VarscanModAFTumorNoncontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%FREQ]' {input} | awk -F'\t' 'BEGIN {{OFS="\t"}} {{ $5 = $5 / 100; print }}'| bgzip -c > {output.tsv} 2>{log} && 
		tabix -b2 -e2 {output.tsv} 2>>{log} &&
		bcftools annotate -x FORMAT/FREQ -a {output.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF  {input} -Oz -o {output.vcf} 2>>{log}
		"""

rule VarscanReheaderTumorNontrol:
	input:
		"SnakeWES/results/{sample}_tumor/{sample}.varscan.concat.vcf.gz"
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.varscan.vcf.gz"
	threads: 1
	params:
		txt="SnakeWES/results/{sample}_tumor/{sample}.rh.txt"
	conda:
		"../envs/bcftools.yaml"
	log:
		"SnakeWES/logs/{sample}.VarscanReheaderTumorNontrol.log"
	shell: 
		"echo {wildcards.sample} > {params.txt} 2>{log} && "
		"bcftools reheader -s {params.txt} -o {output} {input} 2>>{log} && "
		"rm {params.txt} 2>>{log}"

rule FilterVarscanTumorsNocontrol: 
	input:
		"SnakeWES/results/{sample}_tumor/{sample}.varscan.vcf.gz"
	output:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.varscan.filt.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.varscan.filt.vcf.gz.tbi"
	threads: 1
	params:
		excl=config["chr_to_exclude"],
		alt=config["filtering_tumors"]["alt_depth"],
		min_depth=config["filtering_tumors"]["min_depth"],
		vaf=config["filtering_tumors"]["vaf"]
	conda:
		"../envs/bcftools.yaml"
	log:
		"SnakeWES/logs/{sample}.FilterVarscanTumorsNocontrol.log"
	shell:
		"bcftools view -Ov -i'FORMAT/DP >= {params.min_depth} & FORMAT/AD >= {params.alt} & FORMAT/AF >= {params.vaf}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2>{log} && "
		"tabix -f -p vcf {output.vcf} 2>>{log}"

rule MergeVarscanTumorsNocontrol: ## difference between single curly and doucle curly brackets ##
	input:
		vcf=expand(f"SnakeWES/results/{{sample}}_tumor/{{sample}}.varscan.filt.vcf.gz", sample=config['samples'].values()),
		tbi=expand(f"SnakeWES/results/{{sample}}_tumor/{{sample}}.varscan.filt.vcf.gz.tbi", sample=config['samples'].values())
	output:
		vcf="SnakeWES/results/multisample.varscan.vcf.gz",
		tbi="SnakeWES/results/multisample.varscan.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"SnakeWES/logs/multisample.varscan.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} &&"
		"tabix -p vcf {output.vcf} 2>>{log}"

########################################################################### Freebayes ###########################################################################

rule FreebayesCallNormTumorNocontrol: 
	input:
		bam="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		bai="SnakeWES/alignments/{sample}.tumor.dd.rec.bai"
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.freebayes.vcf.gz"
	threads: 1
	conda:
		"../envs/freebayes.yaml"
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	log:
		"SnakeWES/logs/{sample}.freebayes_call_and_norm_on_samples.log"
	shell:
		"freebayes -f {params.ref} -F 0.01 -C 2 -t <(zcat {params.intervals}) --pooled-continuous {input.bam} | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2>{log}"

rule FreebayesIndexTumorsNocontrol:
	input:
		"SnakeWES/results/{sample}_tumor/{sample}.freebayes.vcf.gz"
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.freebayes.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/tabix.yaml"
	log:
		"SnakeWES/logs/{sample}.FreebayesIndexTumorsNocontrol.log"
	shell:
		"tabix -p vcf {input}"

rule FreebayesAddAFTumorNocontrol:
	input:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.freebayes.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.freebayes.vcf.gz.tbi"
	output:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.freebayes.annotAF.vcf.gz",
		#tsv="SnakeWES/results/{sample}_tumor/{sample}.freebayes.annot.tsv.gz",
		#tsv_tbi="SnakeWES/results/{sample}_tumor/{sample}.freebayes.annot.tsv.gz.tbi"
	log:
		"SnakeWES/logs/{sample}.FreebayesAddAFTumorNocontrol.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	params:
		tsv="SnakeWES/results/{sample}_tumor/{sample}.freebayes.annot.tsv.gz",
		tsv_tbi="SnakeWES/results/{sample}_tumor/{sample}.freebayes.annot.tsv.gz.tbi"		
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%RO\t%AO\n' {input.vcf} 2>{log} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {params.tsv} 2>>{log} && 
		tabix -b2 -e2 {params.tsv} 2>>{log} && 
		bcftools annotate -x INFO/AF -a {params.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF -Oz -o {output.vcf} {input.vcf} 2>>{log} &&
		rm {params.tsv} 
		rm {params.tsv_tbi}
		"""

rule FreebayesFilterTumorNocontrol: 
	input:
		"SnakeWES/results/{sample}_tumor/{sample}.freebayes.annotAF.vcf.gz"
	output:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.freebayes.filt.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.freebayes.filt.vcf.gz.tbi"
	log:
		"SnakeWES/logs/{sample}.FreebayesFilterTumorNocontrol.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	params:
		excl=config["chr_to_exclude"],
		depth=config['filtering_tumors']['min_depth'],
		vaf=config['filtering_tumors']['vaf'],
		alt=config['filtering_tumors']['alt_depth']
	shell:
		"bcftools view -i 'QUAL > 1 & INFO/DP >= {params.depth} & FORMAT/AF >= {params.vaf} & INFO/AC >= {params.alt}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf} 2>>{log}"

rule MergeFreebayesTumorOutputNocontrol:
	input:
		vcf=expand(f"SnakeWES/results/{{sample}}_tumor/{{sample}}.freebayes.filt.vcf.gz", sample=config["samples"].values()),
		tbi=expand(f"SnakeWES/results/{{sample}}_tumor/{{sample}}.freebayes.filt.vcf.gz.tbi", sample=config["samples"].values())
	output:
		vcf="SnakeWES/results/multisample.freebayes.vcf.gz",
		tbi="SnakeWES/results/multisample.freebayes.vcf.gz.tbi"
	threads: 1
	log:
		"SnakeWES/logs/multisample.MergeFreebayesTumorOutputNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} && "
		"tabix -p vcf {output.vcf} 2>>{log}"


########################################################################### Strelka 2 ###########################################################################

rule Strelka2Configuration:
	input:
		"SnakeWES/alignments/{sample}.tumor.dd.rec.bam"
	output:
		"SnakeWES/results/{sample}_tumor/strelka2/runWorkflow.py"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.Strelka2Configuration.log"
	params:
		target=config["Strelka2_intervals"],
		ref=config["genome"],
		strelka2=config["Strelka2"]
	shell:
		"{params.strelka2} --bam {input} --referenceFasta {params.ref} --runDir SnakeWES/results/{wildcards.sample}_tumor/strelka2/ --callRegions {params.target} 2>{log}"


rule RunStrelka2:
	input:
		"SnakeWES/results/{sample}_tumor/strelka2/runWorkflow.py"
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.strelka2.vcf.gz"
	threads: 20
	log:
		"SnakeWES/logs/{sample}.RunStrelka2.log"
	params:
		tmp="SnakeWES/results/{sample}_tumor/strelka2/results/variants/variants.vcf.gz",
		ref=config["genome"]
	shell:
		"{input} -m local -j {threads} && bcftools norm -O z -m - -f {params.ref} -o {output} {params.tmp}"

rule Strelka2IndexTumorsNocontrol:
	input:
		"SnakeWES/results/{sample}_tumor/{sample}.strelka2.vcf.gz"
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.strelka2.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/tabix.yaml"
	log:
		"SnakeWES/logs/{sample}.FreebayesIndexTumorsNocontrol.log"
	shell:
		"tabix -p vcf {input}"

rule Strelka2AddAfFieldOnTumorsNocontrol:  
	input:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.strelka2.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.strelka2.vcf.gz.tbi"
	output:
		vcf=temp("SnakeWES/results/{sample}_tumor/{sample}.strelka2.addaf.vcf.gz")
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"SnakeWES/logs/{sample}.Strelka2AddAfFieldOnTumorsNocontrol.log"
	params:
		tsv="SnakeWES/results/{sample}_tumor/{sample}.strelka2.annot.tsv.gz",
		tsv_tbi="SnakeWES/results/{sample}_tumor/{sample}.strelka2.annot.tsv.gz.tbi"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%AD{{0}}\t%AD{{1}}]' {input.vcf} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {params.tsv} 2>{log} && 
		tabix -b2 -e2 {params.tsv} 2>>{log} && 
		bcftools annotate -a {params.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF {input.vcf} -Oz -o {output.vcf} 2>>{log} && 
		rm {params.tsv} 2>>{log} && 
		rm {params.tsv_tbi} 2>>{log}
		"""

rule Strelka2FilterTumorNocontrol: 
	input:
		"SnakeWES/results/{sample}_tumor/{sample}.strelka2.addaf.vcf.gz"
	output:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.strelka2.filt.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.strelka2.filt.vcf.gz.tbi"
	log:
		"SnakeWES/logs/{sample}.Strelka2FilterTumorNocontrol.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	params:
		excl=config["chr_to_exclude"],
		depth=config['filtering_tumors']['min_depth'],
		vaf=config['filtering_tumors']['vaf'],
		alt=config['filtering_tumors']['alt_depth']
	shell:
		"bcftools view -i 'FORMAT/DP >= {params.depth} & FORMAT/AF >= {params.vaf} & FORMAT/AD[0:1] >= {params.alt}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf} 2>>{log}"

########################################################################### Make Consensus ###########################################################################

rule MakeConsensus:
	input:
		vcf=expand(f"SnakeWES/results/{{sample}}_tumor/{{sample}}.{{caller}}.vcf.gz", sample=config["samples"].values(), caller=config["callers"].values()),
		tbi=expand(f"SnakeWES/results/{{sample}}_tumor/{{sample}}.{{caller}}.vcf.gz.tbi", sample=config["samples"].values(), caller=config["callers"].values())
	output:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.consensus.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.consensus.vcf.gz.tbi"
	threads:1
	log:
		"SnakeWES/logs/{sample}.MakeConsensus.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		callers=list(config["callers"].values())[0],
		files= ["SnakeWES/results/{sample}_tumor/{sample}." + x + ".vcf.gz" for x in list(config["callers"].values())]
	shell:
		"bcftools isec -w 1 -n~111 -O z -o {output.vcf} {params.files} && tabix -p vcf {output.vcf}"


########################################################################### VEP ###########################################################################

rule VepSingleSampleNocontrol:
	input:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.{caller}.filt.vcf.gz",
		cache=config["vepcache"]
	output:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.{caller}.filt.vep.vcf.gz",  
		tbi="SnakeWES/results/{sample}_tumor/{sample}.{caller}.filt.vep.vcf.gz.tbi",
		html="SnakeWES/results/{sample}_tumor/{sample}.{caller}.filt.vep.html"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.Vep{caller}SingleSampleNocontrol.log"
	conda:
		"../envs/vep.yaml"
	params:
		ref=config["genome"],
		clinvar=config["clinvar"],
		dbNSFP=config["dbNSFP"],
	shell:
		"vep -i {input.vcf} -o {output.vcf} --fork {threads} --compress_output bgzip --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --dir_cache {input.cache} --cache_version 110 --merged --fasta {params.ref} --format vcf --symbol --no_intergenic --merged --cache --pick --pick_allele --vcf --plugin dbNSFP,{params.dbNSFP},SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,FATHMM_converted_rankscore,FATHMM_pred --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN 2>{log} && tabix {output.vcf} 2>>{log}"

rule VepMultisampleNocontrol:
	input:
		vcf="SnakeWES/results/multisample.{caller}.vcf.gz",
		cache=config["vepcache"]
	output:
		vcf="SnakeWES/results/multisample.{caller}.vep.vcf.gz",
		tbi="SnakeWES/results/multisample.{caller}.vep.vcf.gz.tbi",
		html="SnakeWES/results/multisample.{caller}.vep.html"
	threads:1
	log:
		"SnakeWES/logs/Vep{caller}MultisamplePoN.log"
	conda:
		"../envs/vep.yaml"
	params:
		ref=config["genome"],
		clinvar=config["clinvar"],
		dbNSFP=config["dbNSFP"],
	shell:
		"vep -i {input.vcf} -o {output.vcf} --fork {threads} --compress_output bgzip --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --dir_cache {input.cache} --cache_version 110 --merged --fasta {params.ref} --format vcf --symbol --no_intergenic --merged --cache --pick --pick_allele --vcf --plugin dbNSFP,{params.dbNSFP},SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,FATHMM_converted_rankscore,FATHMM_pred --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN 2>{log} && tabix {output.vcf} 2>>{log}"

########################################################################### SPLIT VEP ###########################################################################

#MUTECT2

rule SplitVepMutect2SingleSampleNocontrol:
	input:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.mutect2.filt.vep.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.mutect2.filt.vep.vcf.gz.tbi",
		html="SnakeWES/results/{sample}_tumor/{sample}.mutect2.filt.vep.html",
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.mutect2.filt.vep.tmp01.tsv"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.SplitVepMutect2SingleSampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ[\t%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepMutect2MultisampleNocontrol:
	input:
		vcf="SnakeWES/results/multisample.mutect2.vep.vcf.gz",
		tbi="SnakeWES/results/multisample.mutect2.vep.vcf.gz.tbi",
		html="SnakeWES/results/multisample.mutect2.vep.html"
	output:
		"SnakeWES/results/multisample.mutect2.vep.tmp01.tsv"
	threads: 1
	log:
		"SnakeWES/logs/SplitVepMutect2MultisampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ[\t%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# FREEBAYES

rule SplitVepFreebayesSingleSampleNocontrol:
	input:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.freebayes.filt.vep.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.freebayes.filt.vep.vcf.gz.tbi",
		html="SnakeWES/results/{sample}_tumor/{sample}.freebayes.filt.vep.html",
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.freebayes.filt.vep.tmp01.tsv"
	threads:1
	log:
		"SnakeWES/logs/{sample}.SplitVepFreebayesSingleSampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ[\t%RO\t%AO\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepFreebayesMultisampleNocontrol:
	input:
		vcf="SnakeWES/results/multisample.freebayes.vep.vcf.gz",
		tbi="SnakeWES/results/multisample.freebayes.vep.vcf.gz.tbi",
		html="SnakeWES/results/multisample.freebayes.vep.html"
	output:
		"SnakeWES/results/multisample.freebayes.vep.tmp01.tsv"
	threads: 1 
	log:
		"SnakeWES/logs/SplitVepFreebayesMultisampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ[\t%RO\t%AO\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# BCFTOOLS 

rule SplitVepBcftoolsSingleSampleNocontrol:
	input:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.bcftools.filt.vep.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.bcftools.filt.vep.vcf.gz.tbi",
		html="SnakeWES/results/{sample}_tumor/{sample}.bcftools.filt.vep.html",
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.bcftools.filt.vep.tmp01.tsv"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.SplitVepBcftoolsSingleSampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ[\t%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepBcftoolsMultisampleNocontrol:
	input:
		vcf="SnakeWES/results/multisample.bcftools.vep.vcf.gz",
		tbi="SnakeWES/results/multisample.bcftools.vep.vcf.gz.tbi",
		html="SnakeWES/results/multisample.bcftools.vep.html"
	output:
		"SnakeWES/results/multisample.bcftools.vep.tmp01.tsv"
	threads: 1
	log:
		"SnakeWES/logs/SplitVepBcftoolsMultisampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ[\t%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# VARSCAN 

rule SplitVepVarscanSingleSampleNocontrol:
	input:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.varscan.filt.vep.vcf.gz",
		tbi="SnakeWES/results/{sample}_tumor/{sample}.varscan.filt.vep.vcf.gz.tbi",
		html="SnakeWES/results/{sample}_tumor/{sample}.varscan.filt.vep.html",
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.varscan.filt.vep.tmp01.tsv"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.SplitVepVarscanSingleSampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ[\t%RD\t%AD\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepVarscanMultisampleNocontrol:
	input:
		vcf="SnakeWES/results/multisample.varscan.vep.vcf.gz",
		tbi="SnakeWES/results/multisample.varscan.vep.vcf.gz.tbi",
		html="SnakeWES/results/multisample.varscan.vep.html"
	output:
		"SnakeWES/results/multisample.varscan.vep.tmp01.tsv"
	threads: 1
	log:
		"SnakeWES/logs/SplitVepVarscanMultisampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ[\t%RD\t%AD\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""
###################################### ADD SAMPLES NAMES ############################################################


rule AddSampleNamesMultisampleNocontrol:
	input:
		vcf="SnakeWES/results/multisample.{caller}.vep.vcf.gz",
		tsv="SnakeWES/results/multisample.{caller}.vep.tmp01.tsv"
	output:
		"SnakeWES/results/multisample.{caller}.vep.tmp02.tsv"
	threads: 1
	params: 
		txt="SnakeWES/results/{caller}_sample_list.txt", 
		script="SnakeWES/workflow/scripts/add_sample_names_wes.py"
	log:
		"SnakeWES/logs/{caller}.AddSampleNamesMultisamplePoN.log"
	shell:
		"bcftools query -l {input.vcf} | sed 's/_tumor//g' > {params.txt} 2>{log} && "
		"python {params.script} {input.tsv} {output} {params.txt} 2>>{log} && "
		"rm {params.txt} 2>>{log}"


##################################### Fill empty field with NA #####################################################


rule FillEmptyFieldwithNAoneSampleNocontrol:
	input:
		"SnakeWES/results/{sample}_tumor/{sample}.{caller}.filt.vep.tmp01.tsv"
	output:
		"SnakeWES/results/{sample}_tumor/{sample}.{caller}.filt.vep.tmp02.tsv"
	threads: 1 
	log:
		"SnakeWES/logs/{sample}.{caller}.FillEmptyFieldwithNAoneSamplePoN.log"
	shell:
		"""
		awk 'BEGIN{{FS=OFS="\t"}} {{for(i=1;i<=NF;i++) if($i==".") $i="NA"}}1' {input} > {output}
		"""

rule FillEmptyFieldwithNAmultisampleNocontrol:
	input:
		"SnakeWES/results/multisample.{caller}.vep.tmp02.tsv"
	output:
		"SnakeWES/results/multisample.{caller}.vep.tmp03.tsv"
	threads: 1 
	log:
		"SnakeWES/logs/{caller}.FillEmptyFieldwithNAoneSamplePoN.log"
	shell:
		"""
		awk 'BEGIN{{FS=OFS="\t"}} {{for(i=1;i<=NF;i++) if($i==".") $i="NA"}}1' {input} > {output}
		"""
####################################   FORMATTING HEADER   ############################################################

rule FormattingHeaderSingleSample:
	input:
		vcf="SnakeWES/results/{sample}_tumor/{sample}.varscan.filt.vep.vcf.gz",
		fi_tsv="SnakeWES/results/{sample}_tumor/{sample}.{caller}.filt.vep.tmp02.tsv"
	output:
		fo_tsv="SnakeWES/results/{sample}_tumor/{sample}.{caller}.filt.vep.tsv"
	params:
		"SnakeWES/workflow/scripts/formatting_header.R"
	threads: 1 
	log:
		"SnakeWES/logs/{sample}.{caller}.FormattingHeaderSingleSample.log"
	conda:
		"../envs/formatHeader.yaml"
	shell:
		"Rscript {params} {input.vcf} {input.fi_tsv} {output.fo_tsv}"

rule FormattingHeaderMultisample:
	input:
		vcf="SnakeWES/results/multisample.{caller}.vep.vcf.gz",
		fi_tsv="SnakeWES/results/multisample.{caller}.vep.tmp03.tsv"
	output:
		fo_tsv="SnakeWES/results/multisample.{caller}.vep.tsv"
	params:
		"SnakeWES/workflow/scripts/formatting_header.R"
	threads: 1 
	log:
		"SnakeWES/logs/{caller}.FormattingHeaderMultisample.log"
	conda:
		"../envs/formatHeader.yaml"
	shell:
		"Rscript {params} {input.vcf} {input.fi_tsv} {output.fo_tsv}"

