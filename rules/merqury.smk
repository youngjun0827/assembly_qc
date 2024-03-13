import glob
import pandas as pd
import os

if not os.path.exists("log"):
    os.makedirs("log")


SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)


configfile: "config/config_asm_qc.yaml"
MANIFEST = config.get('MANIFEST', 'config/manifest_asm_qc.tab')


manifest_df = pd.read_csv(
    MANIFEST, header=0, index_col="SAMPLE", sep="\t", comment="#"
).fillna("NA")



wildcard_constraints:
    sample="|".join(manifest_df["ILLUMINA"].unique()),
    asm="|".join(manifest_df.index),
    read="\d+",



def find_fastq(wildcards):
    fofn_df = pd.read_csv(
        list(manifest_df[manifest_df["ILLUMINA"] == wildcards.sample]["FOFN"].unique())[
            0
        ],
        sep="\t",
        header=None,
        names=["fastq"],
    )
    return fofn_df.at[int(wildcards.read), "fastq"]


def agg_reads(wildcards):
    fofn_df = pd.read_csv(
        list(manifest_df[manifest_df["ILLUMINA"] == wildcards.sample]["FOFN"].unique())[
            0
        ],
        sep="\t",
        header=None,
        names=["fastq"],
    )
    return directory(
        expand(
            rules.run_meryl.output.meryl, sample=wildcards.sample, read=fofn_df.index
        )
    )


def find_cleaned_hap_one(wildcards):
    return f"QC_results/contamination_screening/results/{wildcards.asm}_hap1/fasta/{wildcards.asm}_hap1.fasta"


def find_cleaned_hap_two(wildcards):
    return f"QC_results/contamination_screening/results/{wildcards.asm}_hap2/fasta/{wildcards.asm}_hap2.fasta"


def find_meryl(wildcards):
    return directory(
        rules.meryl_combine.output.meryl.format(
            sample=manifest_df.at[wildcards.asm, "ILLUMINA"],
        )
    )


def find_trios(wildcards):
    return expand(
        "QC_results/merqury/results/{asm}/trio/{asm}_trio.spectra-asm.st.png",
        asm=manifest_df[manifest_df["TRIO"] == "YES"].index,
    )


def find_mat_meryl(wildcards):
    mother = manifest_df.loc[
        (manifest_df["ILLUMINA"] == wildcards.sample) & (manifest_df["MO_ID"] != "NA")
    ].iloc[0]["MO_ID"]
    return directory(expand(rules.meryl_combine.output.meryl, sample=mother))


def find_pat_meryl(wildcards):
    father = manifest_df.loc[
        (manifest_df["ILLUMINA"] == wildcards.sample) & (manifest_df["FA_ID"] != "NA")
    ].iloc[0]["FA_ID"]
    return directory(expand(rules.meryl_combine.output.meryl, sample=father))


def find_hapmer_hist(wildcards):
    return expand(
        rules.hapmers.output.inherited_hist,
        sample=manifest_df.at[wildcards.asm, "ILLUMINA"],
    )


def find_hapmers(wildcards):
    sample = manifest_df.at[wildcards.asm, "ILLUMINA"]
    mother = manifest_df.at[wildcards.asm, "MO_ID"]
    father = manifest_df.at[wildcards.asm, "FA_ID"]
    return " ".join(
        [
            os.path.abspath(f"QC_results/merqury/meryl/{sample}/{sample}_all.meryl"),
            os.path.abspath(
                f"QC_results/merqury/meryl/{sample}/hapmers/{sample}_all_and_{mother}_all.only.meryl"
            ),
            os.path.abspath(
                f"QC_results/merqury/meryl/{sample}/hapmers/{sample}_all_and_{father}_all.only.meryl"
            ),
        ]
    )


localrules:
    all,
    agg_hapmers,


rule all:
    input:
        expand(
            "QC_results/merqury/results/{asm}/{asm}.qv",
            asm=manifest_df.loc[manifest_df["H1"] != "NA"].index,
        ),
        find_trios,


rule non_trio_only:
    input:
        expand(
            "QC_results/merqury/results/{asm}/{asm}.qv",
            asm=manifest_df.loc[manifest_df["H1"] != "NA"].index,
        ),


rule trio_only:
    input:
        find_trios,


rule run_meryl:
    input:
        fastq=find_fastq,
    output:
        meryl=temp(directory("QC_results/merqury/meryl/{sample}/{read}.meryl")),
    resources:
        mem=lambda wildcards, attempt: attempt * 120,
        hrs=48,
    threads: 1
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "meryl/1.0",
    shell:
        """
        meryl k=21 count memory={resources.mem} {input.fastq} output {output.meryl}
        """


rule meryl_combine:
    input:
        meryl=agg_reads,
    output:
        meryl=directory("QC_results/merqury/meryl/{sample}/{sample}_all.meryl"),
    resources:
        mem=8,
        hrs=48,
    threads: 1
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "meryl/1.0",
    shell:
        """
        if [[ $( echo {input.meryl} | wc -w ) == 1 ]]; then
            cp -rl {input.meryl} {output.meryl}
        else
            meryl union-sum {input.meryl} output {output.meryl}
        fi
        """


rule merqury_script:
    input:
        meryl=find_meryl,
        hap_one=find_cleaned_hap_one,
    output:
        run_script="QC_results/merqury/results/{asm}/{asm}_run.sh",
    params:
        hap_two=find_cleaned_hap_two,
    resources:
        mem=4,
        hrs=2,
    threads: 1
    run:
        meryl_abs = os.path.abspath(input.meryl)
        hap_one_abs = os.path.abspath(input.hap_one)
        if params.hap_two.upper() == "NA":
            asm_all = [hap_one_abs]
        else:
            asm_all = [hap_one_abs, os.path.abspath(params.hap_two)]
        asm_all = " ".join(asm_all)
        with open(output.run_script, "w") as outfile:
            outfile.write("#!/usr/bin/env bash \n")
            outfile.write(f"merqury.sh {meryl_abs} {asm_all} {wildcards.asm}\n")
        os.chmod(output.run_script, 0o755)


rule merqury_run:
    input:
        run_script=rules.merqury_script.output.run_script,
    output:
        png="QC_results/merqury/results/{asm}/{asm}.qv",
    resources:
        mem=4,
        hrs=96,
    threads: 16
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "pcre2/10.35",
        "hdf5/1.10.1",
        "R/4.0.0",
        "meryl/1.0",
        "samtools/1.12",
        "bedtools/2.28.0",
        "igvtools/2.8.2",
        "merqury/1.1",
    shell:
        """
        pushd QC_results/merqury/results/{wildcards.asm}/; ./$( basename {input.run_script} ); popd
        """


rule hapmers:
    input:
        mat_meryl=find_mat_meryl,
        pat_meryl=find_pat_meryl,
        asm_meryl=rules.meryl_combine.output.meryl,
    output:
        inherited_hist="QC_results/merqury/meryl/{sample}/hapmers/inherited_hapmers.hist",
    resources:
        mem=4,
        hrs=96,
    threads: 8
    shell:
        """
        pushd $( dirname {output.inherited_hist} )
        $( dirname $( which merqury.sh ) )/trio/hapmers.sh ../../../{input.mat_meryl} ../../../{input.pat_meryl} ../../../{input.asm_meryl}
        popd
        """


rule agg_hapmers:
    input:
        find_hapmer_hist,
    output:
        hflag=touch(temp("QC_results/merqury/results/{asm}/trio/.hapmers_done")),


rule merqury_trio_script:
    input:
        hapmer_flag=rules.agg_hapmers.output.hflag,
        asm_meryl=find_meryl,
        hap_one=find_cleaned_hap_one,
        hap_two=find_cleaned_hap_two,
    output:
        run_script="QC_results/merqury/results/{asm}/trio/{asm}_run.sh",
    params:
        meryl_hapmers=find_hapmers,
    resources:
        mem=4,
        hrs=2,
    threads: 1
    run:
        meryl_all = params.meryl_hapmers
        asm_all = " ".join(
            [os.path.abspath(input.hap_one), os.path.abspath(input.hap_two)]
        )
        out_asm_name = "_".join([wildcards.asm, "trio"])
        with open(output.run_script, "w") as outfile:
            outfile.write("#!/usr/bin/env bash \n")
            outfile.write(f"merqury.sh {meryl_all} {asm_all} {out_asm_name}\n")
        os.chmod(output.run_script, 0o755)


rule merqury_trio_run:
    input:
        run_script=rules.merqury_trio_script.output.run_script,
    output:
        png="QC_results/merqury/results/{asm}/trio/{asm}_trio.spectra-asm.st.png",
    resources:
        mem=4,
        hrs=96,
    threads: 16
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "pcre2/10.35",
        "hdf5/1.10.1",
        "R/4.0.0",
        "meryl/1.0",
        "samtools/1.12",
        "bedtools/2.28.0",
        "igvtools/2.8.2",
        "merqury/1.1",
    shell:
        """
        pushd QC_results/merqury/results/{wildcards.asm}/trio; ./$( basename {input.run_script} ); popd
        """
