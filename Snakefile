import pandas as pd
import os
import sys

# define your configuration via python
# or define a yaml e.g.

configfile: "config/config_asm_qc.yaml"
MANIFEST = config.get('MANIFEST', 'config/manifest_asm_qc.tab')

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
IDEO_PLOT_SCRIPT = "/net/eichler/vol28/7200/software/pipelines/compteam_tools/ideo_plot.py"
N50_SCRIPT = "/net/eichler/vol28/7200/software/pipelines/compteam_tools/n50"
PLOIDY_PLOT_SCRIPT = f"{SNAKEMAKE_DIR}/scripts/ploidy.R"
TAXID = config.get("TAXID", "9606")
INCLUDE_MITO = bool(config.get("INCLUDE_MITO", False))
HPRC_NAMING = bool(config.get("HPRC_NAMING", False))
SMAHT_METADATA = bool(config.get("SMAHT_METADATA", False))
REF_DICT = config["REF"]

## manifest ==================

def get_asm_manifest_df(manifest_df):
    add_haps = {"H2": "hap2", "UNASSIGNED": "un"}
    rows = []
    for idx, row in manifest_df.iterrows():
        if "H1" in manifest_df.columns and pd.notna(row["H1"]) and str(row["H1"]).strip():
            rows.append({"SAMPLE": row["SAMPLE"], "HAP": "hap1", "FASTA": row["H1"]})
        for col in add_haps:
            if col in manifest_df.columns and pd.notna(row[col]) and str(row[col]).strip():
                rows.append({"SAMPLE": row["SAMPLE"], "HAP": add_haps[col], "FASTA": row[col]})
    return pd.DataFrame(rows)

full_manifest_df = pd.read_csv(
    MANIFEST, header=0, sep="\t", comment="#",
    na_values=["", "NA", "na", "N/A"]
)


conv_manifest_df = get_asm_manifest_df(full_manifest_df)
sample_list = sorted(full_manifest_df["SAMPLE"].astype(str).unique())
full_manifest_df.set_index("SAMPLE", inplace=True)
samples_with_asm = [
    s for s in full_manifest_df.index
    if pd.notna(full_manifest_df.at[s, "H1"]) and str(full_manifest_df.at[s, "H1"]).strip()
]

groups = conv_manifest_df[["SAMPLE","HAP"]].drop_duplicates().copy()
conv_manifest_df.set_index(["SAMPLE","HAP"], inplace=True)


wildcard_constraints:
    sample = "|".join(sample_list),
    hap = "|".join(["hap1","hap2","un"]),


## rule functions =========

def get_sample_fcs_final_outputs(wildcards):
    sample = wildcards.sample
    sample_sub = groups[groups["SAMPLE"] == sample]
    outputs = [
        f"results/{sample}/contamination_screening/outputs/final_fasta/{sample}_{row.HAP}.fasta"
        for idx, row in sample_sub.iterrows()
    ]
    if INCLUDE_MITO and (int(TAXID) == 9696):
        outputs += [
            f"results/{sample}/contamination_screening/outputs/mito_fasta/{row.HAP}-mt.fasta"
            for idx, row in sample_sub.iterrows()
        ]
    return outputs

def get_sample_merqury_final_outputs(wildcards):
    sample = wildcards.sample
    trio_result = find_trios(wildcards)
    if len(trio_result) > 0:
        return [f"results/{sample}/merqury/outputs/{sample}.qv"]+trio_result
    else:
        return [f"results/{sample}/merqury/outputs/{sample}.qv"]


def get_sample_saffire_final_outputs(wildcards):
    sample = wildcards.sample
    sample_sub = groups[groups["SAMPLE"] == sample]
    final_outputs = [
        f"results/{sample}/saffire/outputs/safs/{ref}/{row.HAP}.minimap2.saf"
        for idx, row in sample_sub.iterrows()
        for ref in REF_DICT
    ] + [
        f"results/{sample}/saffire/outputs/chrom_cov/{ref}/{row.HAP}.minimap2.chrom_cov.tsv"
        for idx, row in sample_sub.iterrows()
        for ref in REF_DICT
    ] + [
        f"results/{sample}/saffire/outputs/chrom_cov/{ref}/{row.HAP}.minimap2.contigs_chrom_cov.tsv"
        for idx, row in sample_sub.iterrows()
        for ref in REF_DICT
    ]
    chain_files = [
        f"results/{sample}/chain_files/outputs/{sample}_{row.HAP}_To_{ref}.chain"
        for idx, row in sample_sub.iterrows()
        for ref in REF_DICT
    ]
    if int(TAXID) == 9606:
        final_outputs += chain_files
    return final_outputs

def get_sample_compleasm_final_outputs(wildcards):
    sample = wildcards.sample
    return f"results/{sample}/compleasm/outputs/summary/{sample}.summary.tsv"


def get_sample_fasta_stats_outputs(wildcards):
    sample = wildcards.sample
    sample_sub = groups[groups["SAMPLE"] == sample]
    final_outputs = [
        f"results/{sample}/stats/outputs/summary/{row.HAP}.summary.stats"
        for idx, row in sample_sub.iterrows()
    ]
    return final_outputs

def get_sample_moddotplot_outputs(wildcards):
    if not int(TAXID) == 9606:
        final_outputs = []
    else:
        sample = wildcards.sample
        sample_sub = groups[groups["SAMPLE"] == sample]
        final_outputs = [
            f"results/{sample}/moddotplot/outputs/summary/{row.HAP}.generated_acros.tsv"
            for idx, row in sample_sub.iterrows()
        ] + [
            f"results/{sample}/moddotplot/outputs/contig_stats/{row.HAP}.CHM13_lifted_contigs.tsv"
            for idx, row in sample_sub.iterrows()
        ]
    return final_outputs

def get_sample_plots_outputs(wildcards):
    sample = wildcards.sample
    sample_sub = groups[groups["SAMPLE"] == sample]
    final_outputs = [
        f"results/{sample}/plots/outputs/contig_length/{row.HAP}.scaffold.scatter_logged.png"
        for idx, row in sample_sub.iterrows()
    ] + [
        f"results/{sample}/plots/outputs/ideo/{ref}/pdf/{sample}.minimap2.ideoplot.pdf"
        for idx, row in sample_sub.iterrows()
        for ref in REF_DICT
    ] + [
        f"results/{sample}/plots/outputs/ideo/{ref}/pdf/{sample}.minimap2.ideoplot_wide.pdf"
        for idx, row in sample_sub.iterrows()
        for ref in REF_DICT
    ] + [
        f"results/{sample}/plots/outputs/ploidy/CHM13/pdf/{sample}.minimap2.ploidy.pdf"
        for idx, row in sample_sub.iterrows()
    ] + [
        f"results/{sample}/assembly_eval_config/output/config_file/{sample}.config.yaml"
        for idx, row in sample_sub.iterrows()
    ]
    return final_outputs

def get_sample_smaht_dsa_metadata(wildcards):
    sample = wildcards.sample
    sample_sub = groups[groups["SAMPLE"] == sample]
    if SMAHT_METADATA:
        final_outputs = [
            f"results/{sample}/smaht_dsa_metadata/outputs/{sample}.{row.HAP}.DSA_ExternalQualityMetric.tsv"
            for idx, row in sample_sub.iterrows()
        ]
    else:
        final_outputs = []
    return final_outputs

def get_all_outputs(which_one):
    outputs = []
    for sample in samples_with_asm:
        wildcards = type("WC", (), {"sample": sample})()
        if which_one == "cleaned_fasta":
            outputs.extend(get_sample_fcs_final_outputs(wildcards))
        elif which_one == "saf":
            outputs.extend(get_sample_saffire_final_outputs(wildcards))
        elif which_one == "plots":
            outputs.extend(get_sample_plots_outputs(wildcards))
        elif which_one == "moddot_plots":
            outputs.extend(get_sample_moddotplot_outputs(wildcards))
        elif which_one == "qv":
            outputs.extend(get_sample_merqury_final_outputs(wildcards))
        elif which_one == "smaht_dsa_metadata":
            outputs.extend(get_sample_smaht_dsa_metadata(wildcards))
        elif which_one == "stats":
            outputs.extend(get_sample_fasta_stats_outputs(wildcards))
    return outputs

localrules: all, gather_outputs_per_sample



rule all:
    input:
        expand("results/{sample}/complete_flag/all_done",
            sample=samples_with_asm
        )

rule get_cleaned_fasta_only:
    input:
        lambda wildcards: get_all_outputs(which_one = "cleaned_fasta")

rule get_qv_only:
    input:
        lambda wildcards: get_all_outputs(which_one = "qv")

rule get_saf_only:
    input:
        lambda wildcards: get_all_outputs(which_one = "saf")

rule get_plots_only:
    input:
        lambda wildcards: get_all_outputs(which_one = "plots")

rule get_moddotplots_only:
    input:
        lambda wildcards: get_all_outputs(which_one = "moddot_plots")

rule get_metadata_only:
    input:
        lambda wildcards: get_all_outputs(which_one = "smaht_dsa_metadata")

rule get_stats_only:
    input:
        lambda wildcards: get_all_outputs(which_one = "stats")
        
rule get_busco_only:
    input:
        expand("results/{sample}/compleasm/outputs/summary/{sample}.summary.tsv",
            sample=samples_with_asm
        )


rule gather_outputs_per_sample:
    input:
        get_sample_fcs_final_outputs,
        get_sample_merqury_final_outputs,
        get_sample_saffire_final_outputs,
        get_sample_compleasm_final_outputs,
        get_sample_fasta_stats_outputs,
        get_sample_moddotplot_outputs,
        get_sample_plots_outputs,
        get_sample_smaht_dsa_metadata
    output:
        flag = touch("results/{sample}/complete_flag/all_done")


##===include MUST BE HERE.
# include: "rules/fix_sex_chr.smk"
include: "rules/fcs_gx.smk"
include: "rules/merqury.smk"
include: "rules/saffire.smk"
include: "rules/compleasm.smk"
include: "rules/fasta_stats.smk"
include: "rules/moddotplot.smk"
include: "rules/plots.smk"


##=========================

