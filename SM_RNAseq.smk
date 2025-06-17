import re, os, subprocess
import pandas as pd

# command used to create processing env
# mamba create -n processing -y -c conda-forge -c bioconda deeptools pyCRAC salmon star subread fastx_toolkit samtools

#paths
path = "00_raw/"
name_elem = '_1.fastq.gz'  #write file ending here

#references
STAR_INDEX = "/home/tomasz.turowski/seq_references/hg41/hg41_STAR_index/"
GTF = "/home/tomasz.turowski/seq_references/hg41/hg41_annotation_gencode_tRNA_rRNA.gtf"

#parsing file names and preparatory jobs
longName = [n.replace(name_elem,"") for n in os.listdir(path) if n.endswith(name_elem)]
SAMPLES = ["_".join(n.split("_")[:3]) for n in longName]

df_names = pd.DataFrame({
	'longName' : longName,
	'name'	  : SAMPLES}).set_index("name")

print(df_names)
df_names.to_csv('names.tab', sep="\t")

#SnakeMake pipeline

########## OUTPUTS ##########

rule all:
	input:
		expand("02_alignment/{sample}.bam",sample=SAMPLES),
		expand("02_alignment/{sample}.bam.bai",sample=SAMPLES),
		"cleaninig.done",
		"03_FetaureCounts/featureCounts_multimappers.list",
		"03_FetaureCounts/featureCounts_uniq.list",
		expand("04_BigWig/{sample}_raw_plus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_raw_minus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_CPM_plus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_CPM_minus.bw",sample=SAMPLES)


########## PREPROCESSING ##########

########## ALIGNMENT ##########

rule load_genome:
	input:
		index_check = STAR_INDEX+"SAindex",
	params:
		index_dir = STAR_INDEX
	output:
		touch('loading.done')
	conda:
		"envs/processing.yml"
	shell:
		"STAR --genomeLoad LoadAndExit --genomeDir {params.index_dir}"

rule align:
	input:
		read1 = "00_raw/{sample}_1.fastq.gz",
		read2 = "00_raw/{sample}_2.fastq.gz",
		idx = 'loading.done'
	params:
		index_dir = STAR_INDEX,
		prefix = "02_alignment/{sample}_STAR_"
	output:
		bam = "02_alignment/{sample}_STAR_Aligned.out.bam"
	conda:
		"envs/processing.yml"
	shell:
		"STAR --outFileNamePrefix {params.prefix} --readFilesCommand zcat --genomeDir {params.index_dir} --genomeLoad LoadAndKeep --outSAMtype BAM Unsorted --readFilesIn {input.read1} {input.read2}"

########## POSTPROCESSING ##########

rule sort:
	input:
		bam = "02_alignment/{sample}_STAR_Aligned.out.bam"
	output:
		bam = "02_alignment/{sample}.bam",
		bai = "02_alignment/{sample}.bam.bai"
	conda:
		"envs/processing.yml"
	shell:
		"""
		samtools sort {input.bam} > {output.bam}
		samtools index {output.bam}
		"""

rule clean:
	input:
		expand("02_alignment/{sample}.bam.bai",sample=SAMPLES)
	output:
		touch("cleaninig.done")
	conda:
		"envs/processing.yml"
	shell:
		"""
		rm -r logs/
		rm Aligned.out.sam Log.final.out Log.out Log.progress.out SJ.out.tab
		"""

rule featureCounts:
	input:
		bam = expand("02_alignment/{sample}.bam",sample=SAMPLES) #use list of files
	output:
		multi = "03_FetaureCounts/featureCounts_multimappers.list",
		uniq = "03_FetaureCounts/featureCounts_uniq.list"
	params:
		gtf=GTF
	conda:
		"envs/processing.yml"
	shell:
		"""
		featureCounts -M -s 1 -a {params.gtf} -o {output.multi} {input.bam}
		featureCounts -s 1 -a {params.gtf} -o {output.uniq} {input.bam}
		"""

rule BigWigs_CPM:
	input:
		bam = "02_alignment/{sample}.bam",
		bai = "02_alignment/{sample}.bam.bai"
	output:
		bwP = "04_BigWig/{sample}_CPM_plus.bw",
		bwM = "04_BigWig/{sample}_CPM_minus.bw"
	conda:
		"envs/processing.yml"
	shell:
		"""
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwP} --filterRNAstrand reverse --normalizeUsing CPM --binSize 1
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwM} --filterRNAstrand forward --normalizeUsing CPM --binSize 1
		"""

rule BigWigs_raw:
	input:
		bam = "02_alignment/{sample}.bam",
		bai = "02_alignment/{sample}.bam.bai"
	output:
		bwP = "04_BigWig/{sample}_raw_plus.bw",
		bwM = "04_BigWig/{sample}_raw_minus.bw"
	conda:
		"envs/processing.yml"
	shell:
		"""
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwP} --filterRNAstrand reverse --binSize 1
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwM} --filterRNAstrand forward --binSize 1
		"""