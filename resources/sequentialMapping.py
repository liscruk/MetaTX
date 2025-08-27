rule MapSequential:
	input:
		reads = lambda wildcards: config["samples"][wildcards.sample]
	output:
		# bam = RESULTS+"bam/genomic/{sample}.{species}.bam"
		hostBam = RESULTS + "bam/genomic/host/{sample}."+SPECIES[0]+".bam",
		xenoBam = RESULTS + "bam/genomic/xeno/{sample}."+SPECIES[1]+".bam",
		unmapped = RESULTS + "bam/genomic/unmapped/{sample}."+SPECIES[0]+".fq"

	threads:
		config["threads"]
	params:
		# ref = lambda wildcards: refHuman if wildcards.species == SPECIES[0] else refXeno
		hostref = refHuman,
		xenoref = refXeno

	shell:"""

		minimap2 -t {threads} -ax splice -uf -Y  "{params.hostref}" "{input.reads}" | \
		tee >(samtools view -b -q 30 - | samtools sort -@ {threads} -o "{output.hostBam}") | \
		samtools fastq -f 4 - > "{output.unmapped}"

		minimap2 -t {threads} -ax splice -Y "{params.xenoref}" "{output.unmapped}" | \
		samtools view -bS - | samtools sort -@ {threads} -o "{output.xenoBam}"
	"""
