configfile: "config.yaml"

rule all:
	input:
		expand("data/SRR214{acc}_1.fastq", acc=config["accs"]),
		expand("data/SRR214{acc}_2.fastq", acc=config["accs"])

rule get_data:
	output:
		expand("data/SRR214{acc}.sra", acc=config["accs"])
	shell:
		"wget -r -np -nd -k -P data/ ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP062/SRP062005/"

rule extract_fastq:
	input:
		expand("data/SRR214{acc}.sra", acc=config["accs"])
	output:
		expand("data/SRR214{acc}_1.fastq", acc=config["accs"]),
		expand("data/SRR214{acc}_2.fastq", acc=config["accs"])
	shell:
		"fastq-dump --split-files {input} --outdir data/"
