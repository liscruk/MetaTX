### Novel branch: spliced genomic mapping against the hybrid genome + FLAIR isoform discovery
### 1. Map each sample to the hybrid genome (spliced)
### 2. Pool all samples and build one isoform set with `flair transcriptome` (FLAIR >= 3.0)
### 3. Quantify every sample against that common isoform set with `flair quantify`

HYB_GENOME = RESULTS + "refs/Hybrid_genome." + REFNAME + ".fa"
HYB_GTF = RESULTS + "refs/Hybrid_annotation." + REFNAME + ".gtf"
FLAIR_PREFIX = RESULTS + "flair/" + REFNAME

NOVEL_TARGETS = expand(RESULTS+"bam/gn/{sample}.txt", sample = SAMPLES) + [
	FLAIR_PREFIX + ".isoforms.gtf",
	FLAIR_PREFIX + ".quant.counts.tsv",
	FLAIR_PREFIX + ".quant.tpm.tsv"
]

### Spliced alignment against the hybrid genome (parameters as recommended by FLAIR)
rule minimap2splice:
	input:
		reads = lambda wildcards: config["samples"][wildcards.sample],
		ref = HYB_GENOME
	output:
		bam = RESULTS+"bam/gn/{sample}.bam",
		bai = RESULTS+"bam/gn/{sample}.bam.bai"
	threads:
		config["threads"]
	conda:
		"../envs/mapping.yaml"
	shell:"""
		minimap2 -t {threads} -ax splice -s 80 -G 200k --secondary=no {input.ref} {input.reads} | \
		samtools sort -@ {threads} -o {output.bam} -
		samtools index {output.bam}
	"""

### Mapping statistics for genomic mapping
rule flagstatGenomic:
	input:
		bam = RESULTS+"bam/gn/{sample}.bam"
	output:
		stat = RESULTS+"bam/gn/{sample}.txt"
	conda:
		"../envs/mapping.yaml"
	shell:"""
		samtools flagstat {input.bam} > {output.stat}
	"""

### Pool all samples so isoforms are discovered once across the whole experiment
rule mergeGenomic:
	input:
		bams = expand(RESULTS+"bam/gn/{sample}.bam", sample = SAMPLES)
	output:
		bam = RESULTS+"flair/merged.bam",
		bai = RESULTS+"flair/merged.bam.bai"
	threads:
		config["threads"]
	conda:
		"../envs/mapping.yaml"
	shell:"""
		samtools merge -f -@ {threads} -o {output.bam} {input.bams}
		samtools index {output.bam}
	"""

### Isoform discovery (correct + collapse in one step)
rule flairTranscriptome:
	input:
		bam = RESULTS+"flair/merged.bam",
		bai = RESULTS+"flair/merged.bam.bai",
		genome = HYB_GENOME,
		gtf = HYB_GTF
	output:
		bed = FLAIR_PREFIX + ".isoforms.bed",
		gtf = FLAIR_PREFIX + ".isoforms.gtf",
		fa = FLAIR_PREFIX + ".isoforms.fa"
	params:
		prefix = FLAIR_PREFIX
	threads:
		config["threads"]
	conda:
		"../envs/flair.yaml"
	shell:"""
		flair transcriptome -b {input.bam} \
		-g {input.genome} \
		-f {input.gtf} \
		-t {threads} \
		--stringent \
		-o {params.prefix}
	"""

### Reads manifest for flair quantify: sample, condition, batch, reads
### Condition/batch are placeholders; --sample_id_only keeps the sample names as column headers.
rule flairManifest:
	output:
		manifest = RESULTS + "flair/reads_manifest.tsv"
	run:
		with open(output.manifest, "w") as fh:
			for s, fq in SAMPLES.items():
				fh.write(f"{s}\tcond\tbatch1\t{fq}\n")

### Quantify each sample against the common isoform set
rule flairQuantify:
	input:
		manifest = RESULTS + "flair/reads_manifest.tsv",
		fa = FLAIR_PREFIX + ".isoforms.fa",
		bed = FLAIR_PREFIX + ".isoforms.bed"
	output:
		counts = FLAIR_PREFIX + ".quant.counts.tsv",
		tpm = FLAIR_PREFIX + ".quant.tpm.tsv"
	params:
		prefix = FLAIR_PREFIX + ".quant",
		tmp = RESULTS + "flair/tmp"
	threads:
		config["threads"]
	conda:
		"../envs/flair.yaml"
	shell:"""
		flair quantify -r {input.manifest} \
		-i {input.fa} \
		--isoform_bed {input.bed} \
		--sample_id_only \
		--tpm \
		--temp_dir {params.tmp} \
		-t {threads} \
		-o {params.prefix}
	"""
