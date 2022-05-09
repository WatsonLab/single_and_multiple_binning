import glob
import os

def get_run_accessions(file):
	accns = list()

	meta = open(file, mode="r")

	for row in meta:
		arr = row.strip('\n').split('\t')
		accns.append(arr[0])

	return(accns)


SETONE = get_run_accessions("samples.tsv")
SETTWO = SETONE


print(SETONE)

localrules: combine_coverage, filter, mash_build

rule all:
	input: expand("metabat_single/{sample}/", sample=SETONE), expand("metabat_combined/{sample}/", sample=SETONE), expand("checkm/single/{sample}/checkm.tsv", sample=SETONE), expand("checkm/combined/{sample}/checkm.tsv", sample=SETONE), expand("filtered/single/{sample}", sample=SETONE), expand("filtered/combined/{sample}", sample=SETONE), "checkm_collated/single_checkm.tsv"

rule download_fastq:
	output: 
		R1="fastq/{id}_1.fastq.gz",
		R2="fastq/{id}_2.fastq.gz"
	params:
		ID="{id}"
	conda: "envs/fastqdump.yaml"
	shell:
		'''
		fastq-dump --outdir fastq --split-spot --split-files --gzip {params.ID}
		'''

rule trim_paired:
	input:
		R1="fastq/{id}_1.fastq.gz",
		R2="fastq/{id}_2.fastq.gz"
	output:
		R1="trimmed/{id}_1_val_1.fq.gz",
		R2="trimmed/{id}_2_val_2.fq.gz"
	conda: "envs/trim.yaml"
	threads: 2
	shell:
		'''
		trim_galore --cores 2 --output_dir trimmed --illumina --gzip --paired {input.R1} {input.R2} && rm {input.R1} {input.R2}
		'''

rule megahit:
	input:
		R1="trimmed/{id}_1_val_1.fq.gz",
		R2="trimmed/{id}_2_val_2.fq.gz"
	params:
		di="megahit/{id}"
	output: 
		di="megahit/{id}/",
		fa="megahit/{id}/final.contigs.fa"
	conda: "envs/megahit.yaml"
	threads: 8
	shell: "mkdir -p {params.di} && megahit --continue --k-list 27,47,67,87 --kmin-1pass -m 0.95 --min-contig-len 1000 -t {threads} -1 {input.R1} -2 {input.R2} -o {params.di} && rm {params.di}/intermediate_contigs"


rule bwa_index:
	input: "megahit/{id}/final.contigs.fa"
	output: 
		pac="bwa_indices/{id}.fa.pac",
	params:
		idx="bwa_indices/{id}.fa"
	conda: "envs/bwa.yaml"
	threads: 8
	shell:
		'''
		bwa index -p {params.idx} {input}
		'''

rule bwa_mem:
	input:
		R1="trimmed/{id}_1_val_1.fq.gz",
		R2="trimmed/{id}_2_val_2.fq.gz",
		pac="bwa_indices/{id2}.fa.pac"
	output: 
		cov="coverage/{id}.{id2}.txt" 
	params:
		idx="bwa_indices/{id2}.fa",
		bam="bam/{id}.{id2}.bam",
		fla="bam/{id}.{id2}.bam.flagstat"
	conda: "envs/bwa.yaml"
	threads: 8
	shell: 
		'''
		mkdir -p bam
		bwa mem -t 8 {params.idx} {input.R1} {input.R2} | samtools sort -@8 -m 500M -o {params.bam} -
		samtools index {params.bam}

		samtools flagstat {params.bam} > {params.fla}
	
		jgi_summarize_bam_contig_depths --outputDepth {output.cov} {params.bam} && rm {params.bam} {params.bam}.bai
		'''

rule combine_coverage:
	input: expand("coverage/{sample}.{id}.txt", sample=SETONE, allow_missing=True)
	output: "combined_coverage/{id}.txt"
	shell:
		'''
		perl scripts/combine.pl {input} > {output}	
		'''

rule metabat_single:
	input:
		cov="coverage/{id}.{id}.txt",
		asm="megahit/{id}/final.contigs.fa"
	output:
		unb="metabat_single/{id}/{id}.unbinned.fa",
		dir=directory("metabat_single/{id}/")
	params:
		out="metabat_single/{id}/{id}"
	threads: 8
	conda: "envs/metabat2.yaml"
	shell:
		'''
		metabat2 --minContig 1500 -t {threads} -i {input.asm} -a {input.cov} --unbinned -o {params.out}
		'''


rule metabat_combined:
	input:
		cov="combined_coverage/{id}.txt",
		asm="megahit/{id}/final.contigs.fa"
	output:
		unb="metabat_combined/{id}/{id}.unbinned.fa",
		dir=directory("metabat_combined/{id}/")
	params:
		out="metabat_combined/{id}/{id}"
	threads: 8
	conda: "envs/metabat2.yaml"
	shell:
		'''
		metabat2 --minContig 1500 -t {threads} -i {input.asm} -a {input.cov} --unbinned -o {params.out}
		'''


rule checkm_data:
	output: directory("checkm_data")
	conda: "envs/checkm.yaml"
	shell:
		'''
		mkdir -p {output}
		cd {output}
		wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
		gunzip < checkm_data_2015_01_16.tar.gz | tar xvf -
		cd ..
		echo checkm_data | checkm data setRoot checkm_data
		'''


rule checkm_single:
	input: 
		gin="metabat_single/{id}/",
		dat=ancient("checkm_data")
	output: 
		dir=directory("checkm/single/{id}"),
		tsv="checkm/single/{id}/checkm.tsv"
	threads: 8
	conda: "envs/checkm.yaml"
	shell:
		'''
		checkm lineage_wf --tab_table -f {output.tsv} -t {threads} -x fa {input.gin} {output.dir}
		'''
	 

rule checkm_combined:
	input: 
		gin="metabat_combined/{id}/",
		dat=ancient("checkm_data")
	output: 
		dir=directory("checkm/combined/{id}"),
		tsv="checkm/combined/{id}/checkm.tsv"
	threads: 8
	conda: "envs/checkm.yaml"
	shell:
		'''
		checkm lineage_wf --tab_table -f {output.tsv} -t {threads} -x fa {input.gin} {output.dir}
		'''

rule filter:
	input:
		sinc="checkm/single/{id}/checkm.tsv",
		comc="checkm/combined/{id}/checkm.tsv",
		sinb="metabat_single/{id}",
		comb="metabat_combined/{id}"
	output:
		sin=directory("filtered/single/{id}"),
		com=directory("filtered/combined/{id}")
	shell:
		'''
		perl scripts/filter_genomes.pl {input.sinc} {input.sinb} {output.sin}
		perl scripts/filter_genomes.pl {input.comc} {input.comb} {output.com}
		'''


rule collate_checkm:
	input: 
		hed=expand("checkm/single/{sample}/checkm.tsv", sample=SETONE[1]),
		all=expand("checkm/single/{sample}", sample=SETONE)
	output:
		sin="checkm_collated/single_checkm.tsv",
		com="checkm_collated/combined_checkm.tsv"
	shell:
		'''
		mkdir -p checkm_collated
		touch {output.sin}
		touch {output.com}
		
		cat {input.hed} | head -1 >> {output.sin}
		cat checkm/single/*/checkm.tsv | grep -v Bin >> {output.sin}

		cat {input.hed} | head -1 >> {output.com}
		cat checkm/combined/*/checkm.tsv | grep -v Bin >> {output.com}
		'''
	
