#rule clean_get_data:
#	output:
#		expand("data/SRR214{acc}.sra", acc=config["accs"])
#	shell:
#		"wget -r -np -nd -k -P data/ ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP062/SRP062005/"

#rule clean_extract_fastq:
#	input:
#		expand("data/SRR214{acc}.sra", acc=config["accs"])
#	output:
#		expand("data/SRR214{acc}_1.fastq", acc=config["accs"]),
#		expand("data/SRR214{acc}_2.fastq", acc=config["accs"])
#	shell:
#		"fastq-dump --split-files {input} --outdir data/"

rule clean_names_file:
	input:
		expand("data/SRR214{acc}_1.fastq", acc=config["accs"]),
                expand("data/SRR214{acc}_2.fastq", acc=config["accs"])
	output:
		"data/human.files"
	shell:
		"Rscript code/human_names.R"

rule clean_make_contigs:
	input:
		expand("data/SRR214{acc}_1.fastq", acc=config["accs"]),
                expand("data/SRR214{acc}_2.fastq", acc=config["accs"]),
		files="data/human.files"
	output:
		expand("data/human.{prefix}.contigs.{suffix}", prefix=("trim","scrap"), suffix=("fasta", "qual")),
		expand("data/human.contigs.{suffix2}", suffix2=("report", "groups"))
	shell:
		"mothur '#set.dir(input=data, output=data); make.contigs(file={input.files}, processors=12)'"

rule clean_screen_seqs_contigs:
	input:
		fasta="data/human.trim.contigs.fasta",
		groups="data/human.contigs.groups"
	output:
		"data/human.trim.contigs.good.fasta",
		"data/human.trim.contigs.bad.accnos",
		"data/human.contigs.good.groups"
	shell:
		"mothur '#set.dir(input=data, output=data); screen.seqs(fasta={input.fasta}, group={input.groups}, maxambig=0, maxlength=275, maxhomop=8);'"

rule clean_unique_seqs_contigs:
	input:
		fasta="data/human.trim.contigs.good.fasta"
	output:
		"data/human.trim.contigs.good.names",
		"data/human.trim.contigs.good.unique.fasta"
	shell:
		"mothur '#set.dir(input=data, output=data); unique.seqs(fasta={input.fasta})'"

rule clean_count_seqs:
	input:
		names="data/human.trim.contigs.good.names",
		groups="data/human.contigs.good.groups"
	output:
		"data/human.trim.contigs.good.count_table"
	shell:
		"mothur '#set.dir(input=data, output=data); count.seqs(name={input.names}, group={input.groups});'"

rule clean_align_seqs:
	input:
		fasta="data/human.trim.contigs.good.unique.fasta",
		reference="data/references/silva.v4.align"
	output:
		"data/human.trim.contigs.good.unique.align",
		"data/human.trim.contigs.good.unique.align.report",
		"data/human.trim.contigs.good.unique.flip.accnos"
	shell:
		"mothur '#set.dir(input=data, output=data); align.seqs(fasta={input.fasta}, reference={input.reference}, processors=2);'"

rule clean_screen_seqs_silva:
	input:
		count="data/human.trim.contigs.good.count_table",
		fasta="data/human.trim.contigs.good.unique.align"
	output:
		"data/human.trim.contigs.good.unique.good.align",
		"data/human.trim.contigs.good.unique.bad.accnos",
		"data/human.trim.contigs.good.good.count_table"
	shell:
		"mothur '#set.dir(input=data, output=data); screen.seqs(fasta={input.fasta}, count={input.count}, start=5, end=860);'"

rule clean_filter_seqs:
	input:
		fasta="data/human.trim.contigs.good.unique.good.align"
	output:
		"data/human.filter",
		"data/human.trim.contigs.good.unique.good.filter.fasta"
	shell:
		"mothur '#set.dir(input=data, output=data); filter.seqs(fasta={input.fasta}, vertical=T, trump=.);'"
rule clean_unique_seqs_silva:
	input:
		count="data/human.trim.contigs.good.good.count_table",
		fasta="data/human.trim.contigs.good.unique.good.filter.fasta"
	output:
		"data/human.trim.contigs.good.unique.good.filter.count_table",
		"data/human.trim.contigs.good.unique.good.filter.unique.fasta"
	shell:
		"mothur '#set.dir(input=data, output=data); unique.seqs(fasta={input.fasta}, count={input.count});'"

rule clean_precluster:
	input:
		count="data/human.trim.contigs.good.unique.good.filter.count_table",
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.fasta"
	output:
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.count_table"
	shell:
		"mothur '#set.dir(input=data, output=data); pre.cluster(fasta={input.fasta}, count={input.count}, diffs=2);'"

rule clean_chimera_uchime:
	input:
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		count="data/human.trim.contigs.good.unique.good.filter.unique.precluster.count_table"
	output:
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.chimeras",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos"
	shell:
		"mothur '#set.dir(input=data, output=data); chimera.uchime(fasta={input.fasta}, count={input.count}, dereplicate=T);'"

rule clean_remove_seqs:
	input:
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		accnos="data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos"
	output:
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta"
	shell:
		"mothur '#set.dir(input=data, output=data); remove.seqs(fasta={input.fasta}, accnos={input.accnos});'"

rule clean_classify_seqs:
	input:
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta",
		count="data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table",
		trainset_fasta="data/references/trainset14_032015.pds.fasta",
		trainset_tax="data/references/trainset14_032015.pds.tax"
	output:
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary"
	shell:
		"mothur '#set.dir(input=data, output=data); classify.seqs(fasta={input.fasta}, count={input.count}, reference={input.trainset_fasta}, taxonomy={input.trainset_tax}, cutoff=80);'"

rule clean_remove_lineage:
	input:
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta",
		count="data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table",
		taxonomy="data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy"
	output:
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table"
	shell:
		"mothur '#set.dir(input=data, output=data); remove.lineage(fasta={input.fasta}, count={input.count}, taxonomy={input.taxonomy}, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);'"

rule clean_rename_final_files:
	input:
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta",
		count="data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table",
		taxonomy="data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy"
	output:
		"data/human.fasta",
		"data/human.count_table",
		"data/human.taxonomy"
	shell:
		"mv {input.fasta} data/human.fasta; "
		"mv {input.count} data/human.count_table; "
		"mv {input.taxonomy} data/human.taxonomy"
"""
rule clean_dists:
	input:
		fasta="data/human.fasta"
	output:
		"data/human.dist"
	shell:
		"mothur '#set.dir(input=data, output=data); dist.seqs(fasta={input.fasta}, cutoff=0.03)'"
"""
rule clean_cluster:
	input:
		fasta="data/human.fasta",
		count="data/human.count_table",
		taxonomy="data/human.taxonomy"
	output:
		"data/human.dist",
		"data/human.opti_mcc.list",
		"data/human.opti_mcc.sensspec"
	shell:
		"mothur '#set.dir(input=data, output=data); cluster.split(fasta={input.fasta}, count={input.count}, taxonomy={input.taxonomy}, splitmethod=classify, taxlevel=5, cutoff=.03);'"

rule clean_make_shared:
	input:
		list="data/human.opti_mcc.list",
		count="data/human.count_table"
	output:
		"data/human.opti_mcc.shared"
	shell:
		"mothur '#set.dir(input=data, output=data); make.shared(list=human.opti_mcc.list, count=human.count_table, label=0.3);'"

rule clean_classify_otu:
	input:
		list="data/human.opti_mcc.list",
		count="data/human.count_table",
		taxonomy="data/human.taxonomy"
	output:
		"data/human.opti_mcc.0.03.cons.taxonomy",
		"data/human.opti_mcc.0.03.cons.tax.summary"
	shell:
		"mothur '#set.dir(input=data, output=data); classify.otu(list={input.list}, count={input.count}, taxonomy={input.taxonomy}, label=0.03);'"

rule clean_sub_sample:
	input:
		shared="data/human.opti_mcc.shared"
	output:
		"data/human.opti_mcc.0.03.subsample.shared"
	shell:
		"mothur '#set.dir(input=data, output=data); sub.sample(shared={input.shared}, size=10000)'"

rule clean_filter_shared:
	input:
		shared="data/human.opti_mcc.0.03.subsample.shared"
	output:
		"data/human.opti_mcc.0.03.subsample.0.03.filter.shared"
	shell:
		"mothur '#set.dir(input=data, output=data); filter.shared(shared={input.shared}, minpercentsamples=5, makerare=F, minpercent=0.0001)'"
