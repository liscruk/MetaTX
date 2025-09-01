### Perform spliced alignment against hybrid reference genome
rule minimap2splice:
	input:
		reads = lambda wildcards: config["samples"][wildcards.sample],
		ref = RESULTS + "refs/Hybrid_genome." + REFNAME + ".fa"
	output:
		RESULTS+"bam/gn/{sample}.bam"
	threads:
		config["threads"]
	conda:
		"minimap"
	shell:"""

		minimap2 -t {threads} -ax splice -s 80 -G 200k -secondary=no  {input.index} {input.reads} \
		| samtools sort -@ 20 | samtools view -hbS -@ 20 > {output};
		samtools index {output}

	"""

### Perform mapping statistics for genomic mapping
rule flagstat:
	input:
		bam = RESULTS+"bam/gn/{sample}.bam"
	output:
		stat = RESULTS+"bam/gn/{sample}.txt"
	conda:
		"MetaTx"
	shell:"""
		samtools flagstat {input.bam} -@ {threads} > {output.stat};

	"""

rule runFlair:
	input:
		bam = RESULTS+"bam/gn/{sample}.bam",
		hybg = RESULTS + "refs/Hybrid_genome." + REFNAME + ".fa",
		hybgtf = RESULTS + "refs/Hybrid_annotation." + REFNAME + ".gtf",


	output:
		nx = RESULTS+"flair/{sample}",
		dummyOut = RESULTS+"flair/{sample}.isoforms.bed"

	conda:
		"flair"
	threads:
		config["threads"]
	shell:"""
		flair transcriptome -b {input.bam} \
		-g {input.hybg} \
		-t {threads} \
		-f {hybgtf} \
		--stringent \
		-o {output.nx}
		
	"""
