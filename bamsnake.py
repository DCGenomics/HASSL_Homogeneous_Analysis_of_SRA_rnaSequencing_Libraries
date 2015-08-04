

REF = 'GRCh38.p4'
FASTAREF='/home/ubuntu/russ/ncbi/GCF_000001405.30_GRCh38.p4_genomic.fna'
HISATREF = "/home/ubuntu/refs/hisat_index/GRCh38.p4"

# DATASETS = "SRR1295542".split() 
THREADS = 10 

SAMPLES = 'SRR2089122'.split()

rule all: 
	input: expand("{sample}.GRCh38.p4.hisat.sam", sample=SAMPLES), expand("{sample}.hisat.novel.splicesites.txt", sample=SAMPLES), expand("{sample}.transferred.log", sample=SAMPLES), expand("{sample}.transferred", sample=SAMPLES), expand("{sample}.transferred.splices", sample=SAMPLES)


# SRA -> PILEUP -> RAW COUNTS OFF NCBI GFF3 

rule transfer_logs_s3: 
	output: touch("{sample}.transferred.log")
	input:  "{sample}.hisat.log"
	message: "transferring {input}'s logs to S3"
	shell: "s3cmd put {sample}.hisat.1.log s3://ncbi-hackathon-aug/rnamapping/"

rule transfer_bam_s3: 
	output: touch("{sample}.transferred")
	input: "{sample}.GRCh38.p4.hisat.sorted.bam"
	message: "transferring {input} to S3"
	shell: "s3cmd put {input} s3://ncbi-hackathon-aug/rnamapping/"

rule transfer_splices_s3:
	output: touch("{sample}.transferred.splices")
	input: "{sample}.hisat.novel.splicesites.txt"
	message: "transferring hisat splice sites {input} to s3"
	shell: "s3cmd put {input} s3://ncbi-hackathon-aug/rnamapping/"	


rule sort_bam:
	output: "{sample}.GRCh38.p4.hisat.sorted.bam"
	input: "{sample}.GRCh38.p4.hisat.bam"
	message: "sorting {intput} to {output}"
	shell: "samtools sort {input} {output}.sorted"

rule sam_to_bam:
	output: "{sample}.GRCh38.p4.hisat.bam"
	input: "{sample}.GRCh38.p4.hisat.sam"
	message: "converting sam to bam: {input} to {output}"
	shell: "samtools view -bS {input} > {output}"


rule hisat_alignment_two:
	output: "{sample}.GRCh38.p4.hisat.sam", "{sample}.hisat.log"
	input: "{sample}.hisat.novel.splicesites.txt"
	message: "running second pass alignment"
	shell: "hisat -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {THREADS} --sra-acc {sample} --mm -t -S {sample}.GRCh38.p4.hisat.sam --novel-splicesite-infile {sample}.hisat.novel.splicesite.txt 2>> {sample}.hisat.log"


rule hisat_alignment_one: 
	output: "{sample}.hisat.novel.splicesites.txt"
	message: "hisat aligning reads from {sample} to GRCh38.p4 with {THREADS} to produce splicesites"
	shell: "hisat -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x {HISATREF} -p {THREADS} --sra-acc {wildcards.sample} --mm -t -S {wildcards.sample}.GRCh38.p4.firstpass.sam --novel-splicesite-outfile {output} 2> {wildcards.sample}.hisat.log"




