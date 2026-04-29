# assembly_qc

`assembly_qc` is a Snakemake workflow for assembly quality control. The workflow starts from haplotype-resolved assembly FASTA files, performs contamination/adaptor/rDNA cleaning, then generates assembly statistics, reference-alignment summaries, QV/completeness metrics, BUSCO/compleasm summaries, and human-specific visualization outputs.

## Main outputs

For each sample and available haplotype, the workflow can generate:

- cleaned FASTA files after short-contig filtering, FCS-GX/FCS-adaptor screening, optional mitochondrial extraction, rDNA splitting, and final FASTA renaming
- scaffold/contig assembly statistics, including N50-style summaries and telomere annotations
- assembly-to-reference PAF, SAF, chromosome coverage summaries, and chain files
- Merqury QV and completeness outputs from Illumina reads
- optional trio Merqury hap-mer plots when parental IDs are provided
- compleasm/BUSCO gene completeness summaries
- contig length plots, ideograms, ploidy plots, and assembly-eval config files
- human-only ModDotPlot summaries for acrocentric regions
- optional SMaHT DSA metadata tables when `SMAHT_METADATA: True`

## Repository layout

```text
.
├── Snakefile
├── config/
│   ├── config_asm_qc.yaml
│   └── manifest_asm_qc.tab
├── rules/
│   ├── merqury.smk
│   ├── saffire.smk
│   ├── compleasm.smk
│   ├── fasta_stats.smk
│   ├── moddotplot.smk
│   └── plots.smk
├── scripts/
├── runcluster
└── runlocal
```

## Configuration

The default config is:

```yaml
MANIFEST: config/manifest_asm_qc.tab

REF:
  CHM13:
    PATH: /net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta
    CHROMS: /net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/genome.txt
    CYTO: /net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/anno/cyto.bed
  GRCh38:
    PATH: /net/eichler/vol28/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa
    CHROMS: /net/eichler/vol28/eee_shared/assemblies/hg38/no_alt/genome.txt
    CYTO: /net/eichler/vol28/eee_shared/assemblies/hg38/no_alt/anno/cyto.bed

DB_PATH: /net/eichler/vol28/eee_shared/buscodb/
LINEAGE: primates
MODE: busco

TAXID: 9606
INCLUDE_MITO: True
HPRC_NAMING: True
```

### Key config fields

| Field | Description |
| --- | --- |
| `MANIFEST` | Tab-delimited manifest describing assemblies, read FOFNs, and optional trio information. |
| `REF` | Reference genomes to use for minimap2/Saffire outputs, ideograms, ploidy plots, and chain files. Each reference requires `PATH` and `CHROMS`; `CYTO` is optional. |
| `DB_PATH` | BUSCO/compleasm database path. |
| `LINEAGE` | compleasm lineage, such as `primates`. |
| `MODE` | compleasm mode, such as `busco`. |
| `TAXID` | Expected organism taxid for FCS-GX screening. Human-specific downstream outputs are enabled for `9606`. |
| `INCLUDE_MITO` | Whether mitochondrial FASTA extraction should be included when supported by the workflow settings. |
| `HPRC_NAMING` | Whether final FASTA names should follow HPRC-style naming logic. |

## Manifest

The manifest is tab-delimited. Current expected columns are:

```text
SAMPLE  NAME  H1  H2  UNASSIGNED  ILLUMINA  FOFN  TRIO  MO_ID  FA_ID
```

Example:

```text
SAMPLE_verkko  SAMPLE  test_data/asm/verkko.hap1.fasta  test_data/asm/verkko.hap2.fasta  test_data/asm/verkko.unassigned.fasta  SAMPLE  test_data/fofn/reads.fofn  NO  NA  NA
```

### Manifest columns

| Column | Description |
| --- | --- |
| `SAMPLE` | Workflow sample ID. This is used in `results/{SAMPLE}`. |
| `NAME` | Display or output naming value used by the workflow. |
| `H1` | Haplotype 1 assembly FASTA. Required for a sample to be included in the main workflow. |
| `H2` | Haplotype 2 assembly FASTA. Optional. |
| `UNASSIGNED` | Unassigned assembly FASTA. Optional. |
| `ILLUMINA` | Illumina read sample ID used to locate/read the FOFN for Merqury. Multiple assemblies can share the same Illumina ID. |
| `FOFN` | File of file names for Illumina FASTQ inputs used by Merqury/meryl. One FASTQ path per line. |
| `TRIO` | `YES` or `NO`; enables trio Merqury outputs when parental IDs are available. |
| `MO_ID` | Maternal Illumina sample ID for trio Merqury. Use `NA` if not available. |
| `FA_ID` | Paternal Illumina sample ID for trio Merqury. Use `NA` if not available. |

## Running the workflow

Run a dry-run first:

```bash
./runcluster 30 -np
```

Run on the cluster:

```bash
./runcluster 30
```

Run locally:

```bash
./runlocal 30
```

`runcluster` and `runlocal` both resolve the repository path automatically, load `miniconda/4.12.0`, create `log/`, and pass the remaining arguments to Snakemake. `runcluster` submits jobs through DRMAA/SGE with `mfree`, `h_rt`, `threads`, and `-V -cwd -j y` settings. Both wrappers use Singularity/Apptainer and bind `/net`.

To use a different config file:

```bash
./runcluster 30 --configfile path/to/config.yaml
```

## Partial workflow targets

You can run selected output groups instead of the full workflow:

```bash
./runcluster 30 get_cleaned_fasta_only
./runcluster 30 get_stats_only
./runcluster 30 get_saf_only
./runcluster 30 get_qv_only
./runcluster 30 get_busco_only
./runcluster 30 get_plots_only
./runcluster 30 get_moddotplots_only
```

| Target | Outputs |
| --- | --- |
| `get_cleaned_fasta_only` | Final cleaned FASTA files under `results/{sample}/contamination_screening/outputs/final_fasta/`. |
| `get_stats_only` | Per-haplotype FASTA statistics under `results/{sample}/stats/outputs/summary/`. |
| `get_saf_only` | Saffire SAF files, chromosome coverage summaries, trimmed PAFs, and human chain files. |
| `get_qv_only` | Merqury QV and completeness outputs. Trio plots are included when `TRIO=YES`. |
| `get_busco_only` | compleasm/BUSCO summary TSV. |
| `get_plots_only` | Contig length plots, ideograms, ploidy plots, and assembly-eval config files. |
| `get_moddotplots_only` | Human acrocentric ModDotPlot summary and contig statistics. |

## Workflow overview

The full `all` target creates:

```text
results/{sample}/complete_flag/all_done
```

for every sample with a valid `H1` assembly in the manifest.

Internally, the workflow includes these rule modules:

1. `fcs_gx.smk`  
   Links/indexes raw FASTA input, removes very short contigs, runs FCS-GX/adaptor screening, detects mitochondrial and rDNA sequence, splits rDNA, and writes final cleaned FASTA outputs.

2. `merqury.smk`  
   Builds meryl databases from Illumina FASTQ FOFNs, computes Merqury QV/completeness, and optionally generates trio hap-mer outputs.

3. `saffire.smk`  
   Prepares reference FASTA resources, aligns cleaned assemblies to each configured reference with minimap2, trims/filters PAFs, generates SAF/chromosome coverage summaries, and creates chain files for human assemblies.

4. `compleasm.smk`  
   Runs compleasm and summarizes BUSCO-style completeness metrics.

5. `fasta_stats.smk`  
   Splits scaffolds into contigs, calculates telomere annotations, computes per-haplotype and full-genome assembly statistics, and summarizes stats tables.

6. `moddotplot.smk`  
   For human assemblies, identifies acrocentric contigs from CHM13 alignments and generates ModDotPlot-related summaries/plots.

7. `plots.smk`  
   Generates contig length plots, ideogram plots, ploidy plots, and assembly-eval configuration files.

## Important output locations

```text
results/{sample}/
├── complete_flag/
│   └── all_done
├── contamination_screening/
│   └── outputs/
│       ├── final_fasta/{sample}_{hap}.fasta
│       ├── contig_fasta/{sample}_{hap}.fasta
│       ├── mito_fasta/{hap}-mt.fasta
│       └── rdna_fasta/{hap}-rdna.fasta
├── stats/
│   └── outputs/
│       └── summary/{hap}.summary.stats
├── saffire/
│   └── outputs/
│       ├── trimmed_pafs/{ref}/{hap}.minimap2.trimmed.paf
│       ├── safs/{ref}/{hap}.minimap2.saf
│       └── chrom_cov/{ref}/{hap}.minimap2.chrom_cov.tsv
├── chain_files/
│   └── outputs/
│       ├── {sample}_{hap}_To_{ref}.chain
│       └── {sample}_{hap}_To_{ref}.invert.chain
├── merqury/
│   └── outputs/
│       ├── {sample}.qv
│       └── {sample}_{hap}.qv
├── compleasm/
│   └── outputs/
│       └── summary/{sample}.summary.tsv
├── plots/
│   └── outputs/
│       ├── contig_length/{hap}.scaffold.scatter_logged.png
│       ├── ideo/{ref}/pdf/{sample}.minimap2.ideoplot.pdf
│       ├── ideo/{ref}/pdf/{sample}.minimap2.ideoplot_wide.pdf
│       └── ploidy/{ref}/pdf/{sample}.minimap2.ploidy.pdf
└── moddotplot/
    └── outputs/
        ├── summary/{hap}.generated_acros.tsv
        └── contig_stats/{hap}.CHM13_lifted_contigs.tsv
```

Shared resources are written under:

```text
resources/
├── reference/{ref}/
│   ├── genome.fa
│   ├── genome.fa.fai
│   └── genome_index.txt
├── meryl/{ILLUMINA}/
└── acro_target_beds/
```

## Notes

- Human-specific outputs, including ModDotPlot and chain files, depend on `TAXID: 9606`.
- `H1` is required for a sample to be selected by the default workflow.
- `H2` and `UNASSIGNED` are optional; the workflow runs only for haplotypes present in the manifest.
- `FOFN` is used for Merqury/meryl input and should list Illumina FASTQ files one per line.
- The workflow uses Singularity/Apptainer containers through Snakemake’s `--use-singularity` option.

