#rule dirty_get_data:
#	output:
#		expand("data/SRR214{acc}.sra", acc=config["accs"])
#	shell:
#		"wget -r -np -nd -k -P data/ ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP062/SRP062005/"

#rule dirty_extract_fastq:
#	input:
#		expand("data/SRR214{acc}.sra", acc=config["accs"])
#	output:
#		expand("data/SRR214{acc}_1.fastq", acc=config["accs"]),
#		expand("data/SRR214{acc}_2.fastq", acc=config["accs"])
#	shell:
#		"fastq-dump --split-files {input} --outdir data/"

"""
rule dirty_chimera_uchime:
	input:
		fasta="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		count="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.count_table"
	output:
		"data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table",
		"data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.chimeras",
		"data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos"
	shell:
		"mothur '#set.dir(input=data, output=data); chimera.uchime(fasta={input.fasta}, count={input.count}, dereplicate=T);'"

rule dirty_remove_seqs:
	input:
		fasta="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		accnos="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos"
	output:
		"data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta"
	shell:
		"mothur '#set.dir(input=data, output=data); remove.seqs(fasta={input.fasta}, accnos={input.accnos});'"
"""
rule dirty_copy:
	input:
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.count_table"
	output:
		"data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		"data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.count_table"
	shell:
		"cp data/human.trim.contigs.good.unique.good.filter.unique.precluster.fasta data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.fasta; "
		"cp data/human.trim.contigs.good.unique.good.filter.unique.precluster.count_table data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.count_table"

rule dirty_classify_seqs:
	input:
		fasta="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		count="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.count_table",
		trainset_fasta="data/references/trainset14_032015.pds.fasta",
		trainset_tax="data/references/trainset14_032015.pds.tax"
	output:
		"data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.pds.wang.taxonomy",
		"data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.pds.wang.tax.summary"
	shell:
		"mothur '#set.dir(input=data, output=data); classify.seqs(fasta={input.fasta}, count={input.count}, reference={input.trainset_fasta}, taxonomy={input.trainset_tax}, cutoff=80);'"

rule dirty_remove_lineage:
	input:
		fasta="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		count="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.count_table",
		taxonomy="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.pds.wang.taxonomy"
	output:
		"data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.pds.wang.pick.taxonomy",
		"data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta",
		"data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table"
	shell:
		"mothur '#set.dir(input=data, output=data); remove.lineage(fasta={input.fasta}, count={input.count}, taxonomy={input.taxonomy}, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);'"

rule dirty_rename_final_files:
	input:
		fasta="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta",
		count="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table",
		taxonomy="data/dirty.trim.contigs.good.unique.good.filter.unique.precluster.pds.wang.pick.taxonomy"
	output:
		"data/dirty.fasta",
		"data/dirty.count_table",
		"data/dirty.taxonomy"
	shell:
		"mv {input.fasta} data/dirty.fasta; "
		"mv {input.count} data/dirty.count_table; "
		"mv {input.taxonomy} data/dirty.taxonomy"
"""
rule dirty_dists:
	input:
		fasta="data/dirty.fasta"
	output:
		"data/dirty.dist"
	shell:
		"mothur '#set.dir(input=data, output=data); dist.seqs(fasta={input.fasta}, cutoff=0.03)'"
"""
rule dirty_cluster:
	input:
		fasta="data/dirty.fasta",
		count="data/dirty.count_table",
		taxonomy="data/dirty.taxonomy"
	output:
		"data/dirty.dist",
		"data/dirty.opti_mcc.list",
		"data/dirty.opti_mcc.sensspec"
	shell:
		"mothur '#set.dir(input=data, output=data); cluster.split(fasta={input.fasta}, count={input.count}, taxonomy={input.taxonomy}, splitmethod=classify, taxlevel=5, cutoff=.03);'"

rule dirty_make_shared:
	input:
		list="data/dirty.opti_mcc.list",
		count="data/dirty.count_table"
	output:
		"data/dirty.opti_mcc.shared"
	shell:
		"mothur '#set.dir(input=data, output=data); make.shared(list=dirty.opti_mcc.list, count=dirty.count_table, label=0.3);'"

rule dirty_classify_otu:
	input:
		list="data/dirty.opti_mcc.list",
		count="data/dirty.count_table",
		taxonomy="data/dirty.taxonomy"
	output:
		"data/dirty.opti_mcc.0.03.cons.taxonomy",
		"data/dirty.opti_mcc.0.03.cons.tax.summary"
	shell:
		"mothur '#set.dir(input=data, output=data); classify.otu(list={input.list}, count={input.count}, taxonomy={input.taxonomy}, label=0.03);'"

rule dirty_sub_sample:
	input:
		shared="data/dirty.opti_mcc.shared"
	output:
		"data/dirty.opti_mcc.0.03.subsample.shared"
	shell:
		"mothur '#set.dir(input=data, output=data); sub.sample(shared={input.shared}, size=10000)'"

rule dirty_filter_shared:
	input:
		shared="data/dirty.opti_mcc.0.03.subsample.shared"
	output:
		"data/dirty.opti_mcc.0.03.subsample.0.03.filter.shared"
	shell:
		"mothur '#set.dir(input=data, output=data); filter.shared(shared={input.shared}, minpercentsamples=5, makerare=F, minpercent=0.0001)'"
