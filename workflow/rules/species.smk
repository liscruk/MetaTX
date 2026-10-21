### Species composition: reads per species and the ratio between them, per mapping space.
###   tx = hybrid transcriptome BAMs (canonical branch), also reports the salmon EM split
###   gn = hybrid genome BAMs (novel branch)
### Each read is assigned by its primary alignment; reads with primary/secondary hits in both
### species are additionally counted as cross-species ambiguous (needs secondaries in the BAM:
### canonical maps with -N 100, the novel branch uses --secondary=no, so ambiguity is 0 there).

SPACES = (["tx"] if MODE in ("canonical", "both") else []) + (["gn"] if MODE in ("novel", "both") else [])
SPECIES_SCRIPTS = os.path.join(workflow.basedir, "scripts")

SPECIES_TARGETS = expand(RESULTS + "species/{space}/species_summary.tsv", space = SPACES) + \
	expand(RESULTS + "species/{space}/species_composition.{ext}", space = SPACES, ext = ["pdf", "png"])

wildcard_constraints:
	space = "tx|gn"

### Reference name -> species tables
rule txSpeciesMap:
	input:
		hx = txHuman,
		xx = txXeno
	output:
		RESULTS + "refs/tx2species." + REFNAME + ".tsv"
	params:
		hs = SPECIES[0],
		xs = SPECIES[1]
	shell:"""
		dc() {{ case "$1" in *.gz) gzip -dc "$1" ;; *) cat "$1" ;; esac; }}
		ids() {{ awk -v s="$1" '/^>/{{print substr($1, 2) "\\t" s}}'; }}
		{{ dc {input.hx} | ids {params.hs}; dc {input.xx} | ids {params.xs}; }} > {output}
	"""

rule gnSpeciesMap:
	input:
		RESULTS + "refs/Hybrid_genome." + REFNAME + ".fa.fai"
	output:
		RESULTS + "refs/gn2species." + REFNAME + ".tsv"
	params:
		hs = SPECIES[0],
		xs = SPECIES[1]
	shell:"""
		awk -F'\\t' -v a={params.hs} -v b={params.xs} \
		'index($1, a "_") == 1 {{print $1 "\\t" a; next}} index($1, b "_") == 1 {{print $1 "\\t" b}}' \
		{input} > {output}
	"""

### Per-sample read counts per species.
### tx BAMs are already grouped by read name (minimap2 order); gn BAMs are coordinate-sorted -> collate.
rule speciesCounts:
	input:
		bam = RESULTS + "bam/{space}/{sample}.bam",
		smap = lambda w: RESULTS + "refs/" + w.space + "2species." + REFNAME + ".tsv"
	output:
		RESULTS + "species/{space}/{sample}.counts.tsv"
	params:
		script = os.path.join(SPECIES_SCRIPTS, "species_counts.py"),
		tmp = RESULTS + "species/{space}/collate_{sample}",
		species = " ".join(SPECIES[:2])
	threads:
		config["threads"]
	conda:
		"../envs/species.yaml"
	shell:"""
		if [ "{wildcards.space}" = "gn" ]; then
			samtools collate -O -u -@ {threads} -T {params.tmp} {input.bam}
		else
			samtools view -u {input.bam}
		fi | samtools view -F 0x800 - | \
		python {params.script} --map {input.smap} --species {params.species} \
			--sample {wildcards.sample} --space {wildcards.space} --out {output}
	"""

### Merge samples, compute ratios, plot
rule speciesSummary:
	input:
		counts = expand(RESULTS + "species/{{space}}/{sample}.counts.tsv", sample = SAMPLES),
		salmon = lambda w: expand(RESULTS + "counts/{sample}/quant.sf", sample = SAMPLES) if w.space == "tx" else [],
		txmap = lambda w: RESULTS + "refs/tx2species." + REFNAME + ".tsv" if w.space == "tx" else []
	output:
		tsv = RESULTS + "species/{space}/species_summary.tsv",
		pdf = RESULTS + "species/{space}/species_composition.pdf",
		png = RESULTS + "species/{space}/species_composition.png"
	params:
		script = os.path.join(SPECIES_SCRIPTS, "species_summary.py"),
		species = " ".join(SPECIES[:2]),
		salmon = lambda w, input: ("--salmon " + " ".join(input.salmon) + " --txmap " + input.txmap) if w.space == "tx" else "",
		title = lambda w: "Species composition - " + ("transcriptome (tx)" if w.space == "tx" else "genome (gn)"),
		prefix = RESULTS + "species/{space}/species_composition"
	conda:
		"../envs/species.yaml"
	shell:"""
		python {params.script} --counts {input.counts} --species {params.species} {params.salmon} \
			--title "{params.title}" --tsv {output.tsv} --plot {params.prefix}
	"""
