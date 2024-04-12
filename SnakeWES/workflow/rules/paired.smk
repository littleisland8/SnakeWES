rule GetPileupSummariesTumor:
	input:
		bam="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		intervals=config["intervals"],
		variants=config["gnomAD"]
	output:
		"SnakeWES/data/{sample}.tumor.pileup.table"
	threads: 1
	resources:
		mem_mb=4096,
	params:
		extra="",
	log:
		"SnakeWES/logs/{sample}.GetPileupSummariesTumor.log",
	wrapper:
		"v3.3.3/bio/gatk/getpileupsummaries"

rule CalculateContaminationTumor:
	input:
		"SnakeWES/data/{sample}.tumor.pileup.table"
	output:
		contamination="SnakeWES/data/{sample}.tumor.contamination.table",
		segmentation="SnakeWES/data/{sample}.tumor.tumorseg.txt"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.CalculateContaminationTumor.log",
	conda:
		"../envs/gatk4.yaml"
	params:
		#java_opts="-XX:ParallelGCThreads=" + str(config["threads"]),
		mem_mb="-Xmx4G"
		#extra="--tumor-segmentation data/{wildcards.sample}.tumseg.txt"
	shell:
		"gatk --java-options {params.mem_mb} CalculateContamination -I {input} -O {output.contamination} --tumor-segmentation {output.segmentation} > {log} 2>&1"

rule GetPileupSummariesGermline:
	input:
		bam="SnakeWES/alignments/{sample}.germline.dd.rec.bam",
		intervals=config["intervals"],
		variants=config["gnomAD"]
	output:
		"SnakeWES/data/{sample}.germline.pileup.table"
	threads: 1
	resources:
		mem_mb=4096,
	params:
		extra="",
	log:
		"SnakeWES/logs/{sample}.GetPileupSummariesGermlinePaired.log",
	wrapper:
		"v3.3.3/bio/gatk/getpileupsummaries"

rule CalculateContaminationPaired:
	input:
		tumor="SnakeWES/data/{sample}.tumor.pileup.table",
		germline="SnakeWES/data/{sample}.germline.pileup.table"
	output:
		contamination="SnakeWES/data/{sample}.paired.contamination.table",
		segmentation="SnakeWES/data/{sample}.paired.tumorseg.txt"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.CalculateContaminationPaired.log",
	conda:
		"../envs/gatk4.yaml"
	params:
		#java_opts="-XX:ParallelGCThreads=" + str(config["threads"]),
		mem_mb="-Xmx4G"
	shell:
		"gatk --java-options {params.mem_mb} CalculateContamination -I {input.tumor} -matched {input.germline} -O {output.contamination} --tumor-segmentation {output.segmentation} 2>&1> {log}"

rule Mutect2Paired:
	input:
		bamT="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		baiT="SnakeWES/alignments/{sample}.tumor.dd.rec.bai",
		bamC="SnakeWES/alignments/{sample}.germline.dd.rec.bam",
		baiC="SnakeWES/alignments/{sample}.germline.dd.rec.bai",
	output:
		vcf="SnakeWES/results/{sample}.mutect2.paired.vcf.gz",
	threads: 10
	log:
		"SnakeWES/logs/{sample}.mutect.paired.log"
	conda:
		"../envs/gatk4.yaml"
	params:
		ref=config["genome"],
		interval=config["intervals"],
		gnomAD=config["gnomAD"],
		pon=config["gatk_pon"]
	shell:
		"""
		gatk --java-options "-Xmx4G -XX:ParallelGCThreads={threads}" Mutect2 -R {params.ref} -I {input.bamT} -I {input.bamC} -pon {params.pon} --max-reads-per-alignment-start 0 --max-mnp-distance 0 -normal {wildcards.sample}_germline -L {params.interval} --native-pair-hmm-threads {threads} --af-of-alleles-not-in-resource 0.0000025 --germline-resource {params.gnomAD} -O {output.vcf} 2>{log}
		"""

rule FilterMutectCallsPaired:
	input:
		vcf="SnakeWES/results/{sample}.mutect2.paired.vcf.gz",
		ref=config["genome"],
		#bam="data/bam/{sample}.dd.rec.bam",
		intervals=config["intervals"],
		contamination="SnakeWES/data/{sample}.paired.contamination.table", # from gatk CalculateContamination
		segmentation="SnakeWES/data/{sample}.paired.tumorseg.txt", # from gatk CalculateContamination
		#f1r2="data/{sample}.artifacts_prior.tar.gz" # from gatk LearnReadOrientationBias
	output:
		vcf="SnakeWES/results/{sample}.mutect2.paired.filtered.vcf.gz"
	log:
		"SnakeWES/logs/{sample}.FilterMutectCallsPaired.log",
	params:
		#extra="--tumor-segmentation data/{wildcard.sample}.tumseg.txt",  # optional arguments, see GATK docs
		#java_opts="-XX:ParallelGCThreads=" + str(config["threads"])  # optional
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/gatk/filtermutectcalls"

rule NormMutect2Paired:
	input:
		"SnakeWES/results/{sample}.mutect2.paired.filtered.vcf.gz"
	output:
		vcf="SnakeWES/results/{sample}.mutect2.paired.norm.vcf.gz",
		tbi="SnakeWES/results/{sample}.mutect2.paired.norm.vcf.gz.tbi"
	log:
		"SnakeWES/logs/{sample}.NormMutect2Paired.log",
	conda:
		"../envs/bcftools.yaml"
	params:
		genome=config["genome"]
	shell:
		"bcftools norm -m - -f {params.genome} -O z -o {output.vcf} - && tabix -p vcf {output.vcf} 2>{log}"

#######################################################################################   VARSCAN   #######################################################################################

rule GermlineMpileup:
    input:
        # single or list of bam files
        bam="SnakeWES/alignments/{sample}.germline.dd.rec.bam",
        reference_genome=config["genome"],
    output:
        "SnakeWES/data/mpileup/{sample}.germline.mpileup.gz",
    log:
        "SnakeWES/logs/{sample}.GermlineMpileup.log",
    params:
        extra="-d 10000",  # optional
    wrapper:
        "v3.3.3/bio/samtools/mpileup"

rule TumorMpileup:
    input:
        # single or list of bam files
        bam="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
        reference_genome=config["genome"],
    output:
        "SnakeWES/data/mpileup/{sample}.tumor.mpileup.gz",
    log:
        "SnakeWES/logs/{sample}.TumorMpileup.log",
    params:
        extra="-d 10000",  # optional
    wrapper:
        "v3.3.3/bio/samtools/mpileup"

rule BgzipGermlineMpileup:
	input:
		"SnakeWES/data/mpileup/{sample}.germline.mpileup.gz"
	output:
		"SnakeWES/data/mpileup/{sample}.germline.mpileup"
	log:
		"logs/{sample}.BgzipGermlineMpileup.log"
	threads:1
	shell:
		"bgzip -d {input}"

rule BgzipTumorMpileup:
	input:
		"SnakeWES/data/mpileup/{sample}.tumor.mpileup.gz"
	output:
		"SnakeWES/data/mpileup/{sample}.tumor.mpileup"
	log:
		"SnakeWES/logs/{sample}.BgzipTumorMpileup.log"
	threads:1
	shell:
		"bgzip -d {input}"

rule VarscanSomatic:
	input:
		normal_pileup = "SnakeWES/data/mpileup/{sample}.germline.mpileup",
		tumor_pileup = "SnakeWES/data/mpileup/{sample}.tumor.mpileup"
	output:
		snp=temp("SnakeWES/results/{sample}.varscan.paired.snp.vcf"),
		indel=temp("SnakeWES/results/{sample}.varscan.paired.indel.vcf")
	log:
		"SnakeWES/logs/{sample}.VarscanSomatic.log"
	conda:
		"../envs/varscan.yaml"
	message:
		"Calling somatic variants {wildcards.sample}"
	threads:1
	shell:
		"varscan somatic {input.normal_pileup} {input.tumor_pileup} --output-vcf --min-avg-qual 15 --p-value 0.05 --min-var-freq 0.03  --output-snp {output.snp} --output-indel {output.indel} 2>{log}"

rule VarscanBgzipIndex:
	input:
		snp="SnakeWES/results/{sample}.varscan.paired.snp.vcf",
		indel="SnakeWES/results/{sample}.varscan.paired.indel.vcf"
	output:
		snpbg=temp("SnakeWES/results/{sample}.varscan.paired.snp.vcf.gz"),
		indelbg=temp("SnakeWES/results/{sample}.varscan.paired.indel.vcf.gz"),
		snpix=temp("SnakeWES/results/{sample}.varscan.paired.snp.vcf.gz.tbi"),
		indelix=temp("SnakeWES/results/{sample}.varscan.paired.indel.vcf.gz.tbi")
	log:
		"SnakeWES/logs/{sample}.VarscanBgzipIndex.log"
	threads:1
	shell:
		"bgzip {input.snp} && tabix {output.snpbg} 2>{log} && bgzip {input.indel} && tabix {output.indelbg} 2>{log}"


rule VarscanConcatTumorPaired: 
	input:
		vcf_snv="SnakeWES/results/{sample}.varscan.paired.snp.vcf.gz",
		vcf_indel="SnakeWES/results/{sample}.varscan.paired.indel.vcf.gz",
		tbi_snv="SnakeWES/results/{sample}.varscan.paired.snp.vcf.gz.tbi",
		tbi_indel="SnakeWES/results/{sample}.varscan.paired.indel.vcf.gz.tbi"
	output:
		"SnakeWES/results/{sample}.varscan.paired.vcf.gz"
	params: 
		config["genome"]
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"SnakeWES/logs/{sample}.VarscanConcatTumorPaired.log"
	shell: 
		"bcftools concat -a -Ov -o {output} {input.vcf_snv} {input.vcf_indel} 2>{log}"

rule NormVarscanPaired:
	input:
		"SnakeWES/results/{sample}.varscan.paired.vcf.gz"
	output:
		vcf="SnakeWES/results/{sample}.varscan.paired.norm.vcf.gz",
		tbi="SnakeWES/results/{sample}.varscan.paired.norm.vcf.gz.tbi"
	log:
		"SnakeWES/logs/{sample}.NormVarscanPaired.log",
	conda:
		"../envs/bcftools.yaml"
	params:
		genome=config["genome"],
		txt="SnakeWES/results/{sample}.varscan.rh.txt"
	shell:
		"echo {wildcards.sample}_germline > {params.txt} && echo {wildcards.sample}_tumor >> {params.txt} && bcftools view {input} |bcftools reheader -s {params.txt} |bcftools norm -m - -f {params.genome} -O z -o {output.vcf} - && tabix -p vcf {output.vcf} && rm {params.txt} 2>{log}"



#######################################################################################   STRELKA2   #######################################################################################

rule Strelka2ConfigurationPaired:
	input:
		bamT="SnakeWES/alignments/{sample}.tumor.dd.rec.bam",
		baiT="SnakeWES/alignments/{sample}.tumor.dd.rec.bam.bai",
		bamC="SnakeWES/alignments/{sample}.germline.dd.rec.bam",
		baiC="SnakeWES/alignments/{sample}.germline.dd.rec.bam.bai"
	output:
		"SnakeWES/results/{sample}_strelka2Paired/runWorkflow.py"
	threads: 1
	log:
		"SnakeWES/logs/{sample}.Strelka2ConfigurationPaired.log"
	params:
		target=config["intervals"],
		ref=config["genome"],
		strelka2=config["Strelka2_somatic"]
	shell:
		"{params.strelka2} --normalBam {input.bamC} --tumorBam {input.bamT} --referenceFasta {params.ref} --runDir SnakeWES/results/{wildcards.sample}_strelka2Paired/ --callRegions {params.target} 2>{log}"


rule RunStrelka2Paired:
	input:
		"SnakeWES/results/{sample}_strelka2Paired/runWorkflow.py"
	output:
		snps="SnakeWES/results/{sample}_strelka2Paired/results/variants/somatic.snvs.vcf.gz",
		indels="SnakeWES/results/{sample}_strelka2Paired/results/variants/somatic.indels.vcf.gz",
		snps_tbi="SnakeWES/results/{sample}_strelka2Paired/results/variants/somatic.snvs.vcf.gz.tbi",
		indels_tbi="SnakeWES/results/{sample}_strelka2Paired/results/variants/somatic.indels.vcf.gz.tbi"
	threads: 20
	log:
		"SnakeWES/logs/{sample}.RunStrelka2Paired.log"
	shell:
		"{input} -m local -j {threads}"

rule Strelka2ConcatTumorPaired: 
	input:
		vcf_snv="SnakeWES/results/{sample}_strelka2Paired/results/variants/somatic.snvs.vcf.gz",
		vcf_indel="SnakeWES/results/{sample}_strelka2Paired/results/variants/somatic.indels.vcf.gz",
		tbi_snv="SnakeWES/results/{sample}_strelka2Paired/results/variants/somatic.snvs.vcf.gz.tbi",
		tbi_indel="SnakeWES/results/{sample}_strelka2Paired/results/variants/somatic.indels.vcf.gz.tbi"
	output:
		"SnakeWES/results/{sample}.strelka2.paired.vcf.gz"
	params: 
		config["genome"]
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"SnakeWES/logs/{sample}.Strelka2ConcatTumorPaired.log"
	shell: 
		"bcftools concat -a -Oz {input.vcf_snv} {input.vcf_indel} 2>{log}"

rule Strelka2IndexTumorsPaired:
	input:
		"SnakeWES/results/{sample}.strelka2.paired.vcf.gz"
	output:
		"SnakeWES/results/{sample}.strelka2.paired.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/tabix.yaml"
	log:
		"SnakeWES/logs/{sample}.Strelka2IndexTumorsPaired.log"
	shell:
		"tabix -p vcf {input}"

rule Strelka2AddAfFieldOnTumorsPaired:  
	input:
		vcf="SnakeWES/results/{sample}.strelka2.paired.vcf.gz",
		tbi="SnakeWES/results/{sample}.strelka2.paired.vcf.gz.tbi"
	output:
		vcf=temp("SnakeWES/results/{sample}.strelka2.paired.addaf.vcf.gz")
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"SnakeWES/logs/{sample}.Strelka2AddAfFieldOnTumorsPaired.log"
	params:
		tsv="SnakeWES/results/{sample}.strelka2.paired.annot.tsv.gz",
		tsv_tbi="SnakeWES/results/{sample}.strelka2.paired.annot.tsv.gz.tbi"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%AD{{0}}\t%AD{{1}}]' {input.vcf} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {params.tsv} 2>{log} && 
		tabix -b2 -e2 {params.tsv} 2>>{log} && 
		bcftools annotate -a {params.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF {input.vcf} -Oz -o {output.vcf} 2>>{log} && 
		rm {params.tsv} 2>>{log} && 
		rm {params.tsv_tbi} 2>>{log}
		"""

rule NormStrelka2Paired:
	input:
		"SnakeWES/results/{sample}.strelka2.paired.addaf.vcf.gz"
	output:
		vcf="SnakeWES/results/{sample}.strelka2.paired.norm.vcf.gz",
		tbi="SnakeWES/results/{sample}.strelka2.paired.norm.vcf.gz.tbi"
	log:
		"SnakeWES/logs/{sample}.NormStrelka2Paired.log",
	conda:
		"../envs/bcftools.yaml"
	params:
		genome=config["genome"],
		txt="SnakeWES/results/{sample}.strelka2.rh.txt"
	shell:
		"echo {wildcards.sample}_germline > {params.txt} && echo {wildcards.sample}_tumor >> {params.txt} && bcftools view {input} |bcftools reheader -s {params.txt} |bcftools norm -m - -f {params.genome} -O z -o {output.vcf} - && tabix -p vcf {output.vcf} && rm {params.txt} 2>{log}"


#######################################################################################   VEP   #######################################################################################
rule VepMutectPaired:
	input:
		"results/{sample}.mutect2.paired.filtered.norm.vcf.gz"
	output:
		calls="results/{sample}.mutect2.paired.filtered.vep.vcf.gz",
		tbi="results/{sample}.mutect2.paired.filtered.vep.vcf.gz.tbi",
		html="results/{sample}.mutect2.paired.filtered.vep.html"
	threads: 10
	log:
		"logs/{sample}.VepMutectPaired.log"
	conda:
		"../envs/vep.yaml"
	params:
		ref=config["genome"],
		clinvar=config["clinvar"],
		dbNSFP=config["dbNSFP"],
		cache=config["vepcache"]
	shell:
		"vep -i {input} -o {output.calls} --fork {threads} --compress_output bgzip --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --dir_cache {params.cache} --cache_version 110 --merged --fasta {params.ref} --format vcf --symbol --no_intergenic --merged --cache --pick --pick_allele --vcf --plugin dbNSFP,{params.dbNSFP},SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,FATHMM_converted_rankscore,FATHMM_pred --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN 2>{log} && tabix {output.calls} 2>>{log}"

rule VepVarscanPaired:
	input:
		"results/{sample}.varscan.paired.norm.vcf.gz"
	output:
		calls="results/{sample}.varscan.paired.vep.vcf.gz",
		tbi="results/{sample}.varscan.paired.vep.vcf.gz.tbi",
		html="results/{sample}.varscan.paired.vep.html"
	threads: 10
	log:
		"logs/{sample}.VepVarscanPaired.log"
	conda:
		"../envs/vep.yaml"
	params:
		ref=config["genome"],
		clinvar=config["clinvar"],
		dbNSFP=config["dbNSFP"],
		cache=config["vepcache"]
	shell:
		"vep -i {input} -o {output.calls} --fork {threads} --compress_output bgzip --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --dir_cache {params.cache} --cache_version 110 --merged --fasta {params.ref} --format vcf --symbol --no_intergenic --merged --cache --pick --pick_allele --vcf --plugin dbNSFP,{params.dbNSFP},SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,FATHMM_converted_rankscore,FATHMM_pred --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN 2>{log} && tabix {output.calls} 2>>{log}"

#######################################################################################   Merge   #######################################################################################

rule MultisamplePairedMutect2:
	input:
		vcf=expand(f"results/{{sample}}.mutect2.paired.filtered.vep.vcf.gz", sample=config["samples"].values()),
		tbi=expand(f"results/{{sample}}.mutect2.paired.filtered.vep.vcf.gz.tbi", sample=config["samples"].values())
	output:
		temp("results/multisample.mutect2.paired.vep.tmp.vcf.gz")
	threads: 1
	log:
		"logs/MultisamplePairedMutect2.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -O z -o {output} {input.vcf} 2>{log}"

rule MultisamplePairedVarscan:
	input:
		vcf=expand(f"results/{{sample}}.varscan.paired.vep.vcf.gz", sample=config["samples"].values()),
		tbi=expand(f"results/{{sample}}.varscan.paired.vep.vcf.gz.tbi", sample=config["samples"].values())
	output:
		temp("results/multisample.varscan.paired.vep.tmp.vcf.gz")
	threads: 1
	log:
		"logs/MultisamplePairedVarscan.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -O z -o {output} {input.vcf} 2>{log}"

rule FormatMultisamplePairedMutect2:
	input:
		"results/multisample.mutect2.paired.vep.tmp.vcf.gz"
	output:
		"results/multisample.mutect2.paired.vep.vcf.gz"
	threads: 1
	log:
		"logs/FormatMultisamplePairedVarscan.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		tumor="results/tumor.sample.txt"
	shell:
		"""
		bcftools query -l {input} |grep "_tumor" > {params.tumor} 2>{log} && bcftools view -S {params.tumor} -O z -o {output} {input} 2>>{log}
		"""

rule FormatMultisamplePairedVarscan:
	input:
		"results/multisample.varscan.paired.vep.tmp.vcf.gz"
	output:
		"results/multisample.varscan.paired.vep.vcf.gz"
	threads: 1
	log:
		"logs/FormatMultisamplePairedVarscan.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		tumor="results/tumor.sample.txt"
	shell:
		"""
		bcftools query -l {input} |grep "tumor" > {params.tumor} 2>{log} && bcftools view -S {params.tumor} -O z -o {output} {input} 2>>{log}
		"""

rule ParseAnnotationVepMutect2:
	input:
		"results/multisample.mutect2.paired.vep.vcf.gz"
	output:
		"results/multisample.mutect2.paired.tmp01.tsv"
	threads: 1
	log:
		"logs/ParseAnnotationVepMutect2.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		""
	shell:
		"""
		bcftools +split-vep {input} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ[\t%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n"  -d -A tab > {output} 2>{log}
		"""

rule ParseAnnotationVepVarScan:
	input:
		"results/multisample.varscan.paired.vep.vcf.gz"
	output:
		f"results/multisample.varscan.paired.tmp01.tsv"
	threads: 1
	log:
		"logs/ParseAnnotationVepVarScan.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		""
	shell:
		"""
		bcftools +split-vep {input} -f  "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ[\t%RD\t%AD\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""


#######################################################################################   One Sample Paired Analysis   #######################################################################################

rule ParseAnnotationVepMutect2SingleSample:
	input:
		"results/{sample}.mutect2.paired.filtered.vep.vcf.gz"
	output:
		"results/{sample}.mutect2.paired.filtered.vep.tsv"
	threads: 1
	log:
		"logs/{sample}.ParseAnnotationVepMutect2SingleSample.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		sample="{sample}_tumor",
		header="resources/header.vep.txt"
	shell:
		"""
		cat {params.header} <(bcftools view -s {params.sample} {input} |bcftools +split-vep -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ[\t%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n"  -d -A tab |tr "," "\t") > {output} 2>{log}
		"""

rule ParseAnnotationVepVarScanSingleSample:
	input:
		"results/{sample}.varscan.paired.vep.vcf.gz"
	output:
		"results/{sample}.varscan.paired.vep.tsv"
	threads: 1
	log:
		"logs/{sample}.ParseAnnotationVepVarScanSingleSample.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		sample="{sample}_tumor",
		header="resources/header.vep.txt"
	shell:
		"""
		cat {params.header} <(bcftools view -s {params.sample} {input} |bcftools +split-vep -f  "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ[\t%RD\t%AD\t%AF][\t%GT]\n" -d -A tab) > {output} 2>{log}
		"""
