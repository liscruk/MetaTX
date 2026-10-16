# np_metatx

Snakemake pipeline for meta-/xeno-transcriptomics on Oxford Nanopore long reads.
Reads from mixed samples (e.g. human host + mouse xenograft, or human + *Aspergillus fumigatus*)
are mapped **competitively** against a hybrid host + xeno reference and quantified per transcript,
so species assignment follows from the transcript of origin.

## Modes

| Mode | What it does |
|------|--------------|
| `Canonical` (default) | Build hybrid transcriptome → `minimap2 -ax map-ont` → `samtools flagstat` → `salmon quant` (alignment mode) → `salmon quantmerge` |
| `Novel` | Build hybrid genome/annotation → spliced genomic mapping (`minimap2 -ax splice`) → FLAIR isoform discovery (*work in progress*) |
| `Both` | Both of the above |

## Layout

```
config/config.default.yml      template config (copy and adapt; project configs are git-ignored)
workflow/snakefile             entry point: config, mode switch, rule all, hybrid reference (makeRef)
workflow/rules/canonical.smk   canonical branch: transcriptome mapping + salmon quantification
workflow/rules/novel.smk       novel branch: spliced genomic mapping + FLAIR (work in progress)
workflow/envs/*.yaml           conda environments (mapping, salmon, flair)
resources/sequentialMapping.py experimental sequential host→xeno mapping (not wired in)
```

## Requirements

- Snakemake ≥ 7 with conda/mamba
- Host and xeno references: genome FASTA, GTF (not GFF3) and transcriptome FASTA.
  Transcript IDs must be distinct between species (e.g. `ENST…` vs `ENSMUST…`).

## Usage

```bash
cp config/config.default.yml config/config_myproject.yml
# edit paths, samples, species and mode

snakemake -s workflow/snakefile \
  --configfile config/config_myproject.yml \
  --use-conda --cores 10
```

## Outputs (in `results`)

```
refs/Hybrid_genome.<hs.xx>.fa         concatenated genomes
refs/Hybrid_annotation.<hs.xx>.gtf    concatenated annotations
refs/Hybrd_transcriptome.<hs.xx>.fa   concatenated transcriptomes
bam/tx/<sample>.bam / .txt            transcriptome alignments + flagstat
counts/<sample>/quant.sf              per-sample salmon quantification
counts/tx/quant.tsv                   merged transcript table (all samples)
```

## License

See `LICENSE`.
