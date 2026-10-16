### Canonical branch: competitive mapping against the hybrid transcriptome + salmon quantification

CANONICAL_TARGETS = [RESULTS+"counts/tx/quant.tsv"] + \
	expand(RESULTS+"bam/tx/{sample}.txt", sample = SAMPLES)

### Canonical workflow perform competetive hybrid mapping on a hybrid reference generated per run.
### If you stick to the model it might make sense to keep the reference stable. 
rule Map:
	input:
		reads = lambda wildcards: config["samples"][wildcards.sample],
		ref = RESULTS + "refs/Hybrd_transcriptome." + REFNAME + ".fa"
	output:
		bam = RESULTS+"bam/tx/{sample}.bam"
	threads:
		config["threads"]
	params:
		ref = RESULTS + "refs/Hybrd_transcriptome." + REFNAME + ".fa"
	conda:
		"../envs/mapping.yaml"

	shell:"""
		minimap2 -t {threads} -ax map-ont -Y {params.ref} {input.reads} | \
		samtools view -bS - | \
		samtools sort -@ {threads} -o {output.bam};
	"""
### Simple mapping statistics gathering. We'd like to have good map rates.
rule flagstat:
	input:
		bam = RESULTS+"bam/tx/{sample}.bam"
	output:
		stat = RESULTS+"bam/tx/{sample}.txt"
	conda:
		"../envs/mapping.yaml"
	
	shell:"""
		samtools flagstat {input.bam} -@ {threads} > {output.stat};
	"""
### Generation quantifcation by salmon for seperate mapping
rule salmonQuant:
	input:
		bam = RESULTS+"bam/tx/{sample}.bam",
		trp = RESULTS + "refs/Hybrd_transcriptome." + REFNAME + ".fa"
	output:
		trp_quant = RESULTS+"counts/{sample}/quant.sf",
		trp_dir = directory(RESULTS+"counts/{sample}/")
	threads:
		config["threads"]
	conda:
		"../envs/salmon.yaml"
	shell:
		"salmon quant --noErrorModel -p {threads} -l U -t {input.trp} -a {input.bam} -o {output.trp_dir}"
### Merge tpms
rule salmonMerge:
	input:
		dirs = expand(RESULTS+"counts/{sample}/", sample=SAMPLES, species = SPECIES)
	output:
		tx = RESULTS+"counts/tx/quant.tsv",
	conda:
		"../envs/salmon.yaml"
	shell:"""
		salmon quantmerge --quants {input.dirs} --output {output.tx};
	"""

# rule salmonMergeGene:
# 	input:
# 		dirs = expand(RESULTS+"counts/{sample}_{species}/", sample=SAMPLES, species = SPECIES)
# 	output:
# 		mFs = RESULTS+"counts/gene/quant.tsv"
# 	conda:
# 		"salmon"
# 	shell:"""
# 		salmon quantmerge --genes --quants {input.dirs} --output {output}
# 	"""
