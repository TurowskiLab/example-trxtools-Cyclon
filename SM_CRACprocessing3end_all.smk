import re, os, subprocess
import pandas as pd

# command used to create processing env
# mamba create -n processing -y -c conda-forge -c bioconda deeptools pyCRAC salmon star subread fastx_toolkit samtools

#paths
path = "00_raw/"
name_elem = '.fastq.gz'  #write file ending here

#references
STAR_INDEX = "/home/tomasz.turowski/seq_references/hg41/hg41_STAR_index/"
GTF = "/home/tomasz.turowski/seq_references/hg41/hg41_annotation_gencode_tRNA_rRNA.gtf"

#parsing file names and preparatory jobs
longName = [n.strip(name_elem) for n in os.listdir(path) if n.endswith(name_elem)]
adaptors = [n.split("_")[0] for n in longName]
barcodes = [n.split("_")[1] for n in longName]
SAMPLES = ["_".join(n.split("_")[2:]) for n in longName]

def parseBarcode(bc):
    toGrep = "".join(re.findall(r'([ATCG])', bc))
    [firstN,postN]=[len(i) for i in bc.split(toGrep)]
    return firstN, toGrep, postN

def readLengths(n):
	f = path+n+name_elem
	command = "zcat "+f+" | head -40 | scripts/fastqReadsLength.awk | cut -f1"
	return int(subprocess.run([command],shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

df_names = pd.DataFrame({
	'adaptor' : adaptors,
	'barcode' : barcodes,
	'barcode_pattern': [re.sub(r'([ATCG])', 'X', bc) for bc in barcodes],
	'bcLen'	  : [len(bc)+1 for bc in barcodes], #+1 for fastx_trimmer
	'longName' : longName,
	'readLength': [readLengths(n) for n in longName],
	'name'	  : SAMPLES}).set_index("name")

print(df_names)
df_names.to_csv('names.tab', sep="\t")

d1_name = df_names['longName'].to_dict()
d2_bcLen = df_names['bcLen'].to_dict()
d3_as = df_names['adaptor'].to_dict()
d4_readLen = df_names['readLength'].to_dict()
d5_bcPattern = df_names['barcode_pattern'].to_dict()

#SnakeMake pipeline

########## OUTPUTS ##########

rule all: #TODO eventually produce all files, esp. BWs
	input:
		expand("01_preprocessing/01a_{sample}_noumi.fastq.gz", sample=SAMPLES),
		expand("01_preprocessing/01d_{sample}_flexbar.fastq.gz",sample=SAMPLES),
		# expand("01_preprocessing/01d_{sample}_cutadapt.fastq.gz",sample=SAMPLES),
		# expand("01_preprocessing/01e_{sample}_3end.fasta.gz", sample=SAMPLES),
		# expand("02_alignment/{sample}.bam",sample=SAMPLES),
		# expand("02_alignment/{sample}.bam.bai",sample=SAMPLES),
		expand("02b_alignment_all/{sample}_all.bam", sample=SAMPLES),
		expand("02b_alignment_all/{sample}_all.bam.bai", sample=SAMPLES),
		expand("02b_alignment_all/{sample}_all_dedup.bam", sample=SAMPLES),
		expand("02b_alignment_all/{sample}_all_dedup.bam.bai", sample=SAMPLES),
		"03_FetaureCounts/featureCounts_multimappers.list",
		"03_FetaureCounts/featureCounts_uniq.list",
		expand("04_BigWig/{sample}_all_dedup_raw_plus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_all_dedup_raw_minus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_all_dedup_CPM_plus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_all_dedup_CPM_minus.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}.sam",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_5end_fwd.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_5end_rev.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_3end_fwd.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_3end_rev.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_3end_polyA_fwd.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_3end_polyA_rev.bw",sample=SAMPLES)


########## PREPROCESSING ##########

def bcFile(wildcards):
	sample_name = wildcards.sample
	return path+d1_name[sample_name]+name_elem

def bcPattern(wildcards):
	sample_name = wildcards.sample
	return d5_bcPattern[sample_name]

# rule QC:
# 	input:
# 		bcFile
# 	params:
# 		"01_preprocessing/01a_{sample}_QC"
# 	output:
# 		"01_preprocessing/01a_{sample}_QC.fastq.gz"
# 	conda:
# 		"envs/flexbar.yml"
# 	shell:
# 		"flexbar -r {input} -t {params} -q TAIL -n 4 -qf i1.8 -qt 20 -z GZ"

rule umitools_extract:
	input:
		bcFile
	params:
		bc = bcPattern
	output:
		"01_preprocessing/01a_{sample}_noumi.fastq.gz"
	conda:
		"envs/umitools_fix.yml" #hotfix from commit d380e8d as per issue #688
	shell:
		"umi_tools extract --stdin={input} --stdout={output} --bc-pattern={params.bc}"

# rule umitools_extract:
# 	input:
# 		bcFile
# 	params:
# 		bc = bcPattern
# 	output:
# 		"01_preprocessing/01a_{sample}_noumi.fastq.gz"
# 	conda:
# 		"envs/umitools.yml"
# 	shell:
# 		"gunzip -c {input} | umi_tools extract --stdout={output} --bc-pattern={params.bc} | tr -d '\r' | gzip > {output}" 

# rule collapsing:
# 	input:
# 		"01_preprocessing/01a_{sample}_noumi.fastq.gz"
# 		# bcFile
# 	output:
# 		"01_preprocessing/01b_{sample}_comp.fasta.gz"
# 	conda:
# 		"envs/processing.yml"
# 	shell:
# 		"gunzip -c {input} | fastx_collapser | gzip > {output}"

#functions for debarcoding
def bcLen(wildcards):
	sample_name = wildcards.sample
	return d2_bcLen[sample_name]

rule debarcoding:
	input:
		"01_preprocessing/01a_{sample}_noumi.fastq.gz" #changed from 01b to 01a
		# bcFile
	params:
		bcLen = bcLen
	output:
		"01_preprocessing/01c_{sample}_deBC.fastq.gz"
	conda:
		"envs/processing.yml"
	shell:
		"gunzip -c {input} | fastx_trimmer -f {params.bcLen} | gzip > {output}"
		# "gunzip -c {input} | tr -d '\r' | fastx_trimmer -f {params.bcLen} | gzip > {output}"

#functions for adaptor sequence
def adSeq(wildcards):
	sample_name = wildcards.sample
	return d3_as[sample_name]

rule flexbar_3end_trimming:
	input:
		"01_preprocessing/01c_{sample}_deBC.fastq.gz"
		# bcFile
	params:
		adSeq = adSeq,
		n = "01_preprocessing/01d_{sample}_flexbar"
	output:
		"01_preprocessing/01d_{sample}_flexbar.fastq.gz"
	conda:
		"envs/flexbar.yml" 
	shell:
		"flexbar -r {input} -t {params.n} -as {params.adSeq} -ao 4 -u 3 -m 7 -n 4 -bt RIGHT -z GZ"

rule cutadapt:
	input:
		"01_preprocessing/01c_{sample}_deBC.fastq.gz"
	params:
		adSeq = adSeq
	output:
		"01_preprocessing/01d_{sample}_cutadapt.fastq.gz"
	conda:
		"envs/cutadapt.yml"
	shell:
		"cutadapt -a {params.adSeq} -O 4 -m 7 -j 4 -o {output} {input}"

def maxLen(wildcards):
	sample_name = wildcards.sample
	readLen = d4_readLen[sample_name]
	bcLen = d2_bcLen[sample_name]
	return readLen-bcLen-4

rule length_filtering:
	input:
		"01_preprocessing/01d_{sample}_flexbar.fasta.gz"
	params:
		maxLen = maxLen
	output:
		"01_preprocessing/01e_{sample}_3end.fasta.gz"
	shell:
		"zcat {input} | ./scripts/lenFilterFastaMax.awk -v var={params.maxLen} | gzip > {output}"

########## ALIGNMENT ##########

rule align_all:
	input:
		reads_all = "01_preprocessing/01d_{sample}_flexbar.fastq.gz",
	params:
		index_dir = STAR_INDEX,
		prefix_all = "02b_alignment_all/{sample}_all_STAR_"
	output:
		bam_all = "02b_alignment_all/{sample}_all_STAR_Aligned.out.bam"
	conda:
		"envs/processing.yml"
	shell:
		"STAR --outFileNamePrefix {params.prefix_all} --readFilesCommand zcat --genomeDir {params.index_dir} --genomeLoad LoadAndRemove --outSAMtype BAM Unsorted --readFilesIn {input.reads_all}"

# rule align:
# 	input:
# 		reads = "01_preprocessing/01e_{sample}_3end.fasta.gz",
# 	params:
# 		index_dir = STAR_INDEX,
# 		prefix = "02_alignment/{sample}_STAR_"
# 	output:
# 		bam = "02_alignment/{sample}_STAR_Aligned.out.bam"
# 	conda:
# 		"envs/processing.yml"
# 	shell:
# 		"STAR --outFileNamePrefix {params.prefix} --readFilesCommand zcat --genomeDir {params.index_dir} --genomeLoad LoadAndRemove --outSAMtype BAM Unsorted --readFilesIn {input.reads}"

########## POSTPROCESSING ##########

rule sort_all:
	input:
		bam_all = "02b_alignment_all/{sample}_all_STAR_Aligned.out.bam"
	output:
		bam_all = "02b_alignment_all/{sample}_all.bam",
		bai_all = "02b_alignment_all/{sample}_all.bam.bai"
	conda:
		"envs/processing.yml"
	shell:
		"""
		samtools sort {input.bam_all} > {output.bam_all}
		samtools index {output.bam_all}
		"""

# rule sort:
# 	input:
# 		bam = "02_alignment/{sample}_STAR_Aligned.out.bam"
# 	output:
# 		bam = "02_alignment/{sample}.bam",
# 		bai = "02_alignment/{sample}.bam.bai"
# 	conda:
# 		"envs/processing.yml"
# 	shell:
# 		"""
# 		samtools sort {input.bam} > {output.bam}
# 		samtools index {output.bam}
# 		"""

rule umitools_dedup_all:
	input:
		bam = "02b_alignment_all/{sample}_all.bam",
		bai = "02b_alignment_all/{sample}_all.bam.bai"
	output:
		bam = temp("02b_alignment_all/{sample}_all_dedup_unsorted.bam")
	conda:
		"envs/umitools_fix.yml"
	shell:
		"""
		umi_tools dedup -I {input.bam} -S {output.bam}
		"""

# rule umitools_dedup:
# 	input:
# 		bam = "02_alignment/{sample}.bam",
# 		bai = "02b_alignment_all/{sample}.bam.bai"
# 	output:
# 		bam = "02_alignment/{sample}_dedup.bam"
# 	conda:
# 		"envs/umitools.yml"
# 	shell:
# 		"""
# 		umi_tools dedup -I {input.bam} -S {output.bam} --umi-separator=":"
# """

rule sort_dedup:
	input:
		bam_all = "02b_alignment_all/{sample}_all_dedup_unsorted.bam"
	output:
		bam_all = "02b_alignment_all/{sample}_all_dedup.bam",
		bai_all = "02b_alignment_all/{sample}_all_dedup.bam.bai"
	conda:
		"envs/processing.yml"
	shell:
		"""
		samtools sort {input.bam_all} > {output.bam_all}
		samtools index {output.bam_all}
		"""

rule featureCounts:
	input:
		bam = expand("02b_alignment_all/{sample}_all_dedup.bam",sample=SAMPLES) #use list of files
	output:
		multi = "03_FetaureCounts/featureCounts_multimappers.list",
		uniq = "03_FetaureCounts/featureCounts_uniq.list"
	params:
		gtf=GTF
	conda:
		"envs/processing.yml"
	shell:
		# Use -M for multimappers, -s 1 for strandedness, -O for overlapping reads
		# -a for annotation file, -o for output file
		"""
		featureCounts -M -s 1 -O -a {params.gtf} -o {output.multi} {input.bam}
		featureCounts -s 1 -O -a {params.gtf} -o {output.uniq} {input.bam}
		"""

rule BigWigs_CPM:
	input:
		bam = "02b_alignment_all/{sample}_all_dedup.bam",
		bai = "02b_alignment_all/{sample}_all_dedup.bam.bai"
	output:
		bwP = "04_BigWig/{sample}_all_dedup_CPM_plus.bw",
		bwM = "04_BigWig/{sample}_all_dedup_CPM_minus.bw"
	conda:
		"envs/processing.yml"
	shell:
		"""
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwP} --filterRNAstrand reverse --normalizeUsing CPM --binSize 1
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwM} --filterRNAstrand forward --normalizeUsing CPM --binSize 1
		"""

rule BigWigs_raw:
	input:
		bam = "02b_alignment_all/{sample}_all_dedup.bam",
		bai = "02b_alignment_all/{sample}_all_dedup.bam.bai"
	output:
		bwP = "04_BigWig/{sample}_all_dedup_raw_plus.bw",
		bwM = "04_BigWig/{sample}_all_dedup_raw_minus.bw"
	conda:
		"envs/processing.yml"
	shell:
		"""
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwP} --filterRNAstrand reverse --binSize 1
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwM} --filterRNAstrand forward --binSize 1
		"""

rule bam2sam:
	input:
		bam = "02b_alignment_all/{sample}.bam",
		bai = "02b_alignment_allt/{sample}.bam.bai"
	output:
		sam = "04_BigWig/{sample}.sam"
	conda:
		"envs/processing.yml"
	shell:
		"""
		samtools view -h {input.bam} > {output.sam}
		"""

# rule trxtools_5end:
# 	input:
# 		"04_BigWig/{sample}.sam"
# 	output:
# 		"04_BigWig/{sample}_PROFILE_5end_fwd.bw",
# 		"04_BigWig/{sample}_PROFILE_5end_rev.bw"
# 	conda:
# 		"envs/processing.yml"
# 	shell:
# 		"""
# 		SAM2profilesGenomic.py -u 5end -f {input}
# 		"""

# rule trxtools_3end:
# 	input:
# 		"04_BigWig/{sample}.sam"
# 	output:
# 		"04_BigWig/{sample}_PROFILE_3end_fwd.bw",
# 		"04_BigWig/{sample}_PROFILE_3end_rev.bw",
# 		"04_BigWig/{sample}_PROFILE_3end_polyA_fwd.bw",
# 		"04_BigWig/{sample}_PROFILE_3end_polyA_rev.bw"
# 	conda:
# 		"envs/processing.yml"
# 	shell:
# 		"""
# 		SAM2profilesGenomic.py -u 3end -n -f {input}
# 		"""