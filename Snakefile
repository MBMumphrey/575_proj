configfile: "config.yaml"

include: 'code/process-clean.smk'
include: 'code/process-dirty.smk'

rule all:
	input:
		"data/human.opti_mcc.0.03.subsample.0.03.filter.shared",
		"data/human.opti_mcc.0.03.cons.taxonomy",
		"data/dirty.opti_mcc.0.03.subsample.0.03.filter.shared",
		"data/dirty.opti_mcc.0.03.cons.taxonomy"
