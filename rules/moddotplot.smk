"""
Adopted from acro_parm.smk created by Jiadong Lin
"""

import os
import pandas as pd
import numpy as np
import re
import glob
import pysam

human_acro_bed = [
    ["chr13","0","22508596"],
    ["chr14","0","17708411"],
    ["chr15","0","22694466"],
    ["chr21","0","16306378"],
    ["chr22","0","20711065"]
]

bed_df = pd.DataFrame(human_acro_bed, columns = ["chr", "start", "end"], dtype=str)
bed_df["NAME"] = bed_df["chr"] + "_" + (bed_df["start"].astype(int) + 1).astype(str) + "_" + bed_df["end"]
bed_df.set_index("NAME", drop=True, inplace=True)

acros = ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']

try:
    REF = config["REF"]["CHM13"]["PATH"]
except KeyError:
    REF = "/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta"

# rustybam/0.1.16


def find_tigs_flag(wildcards):
    return [ f"results/{wildcards.sample}/moddotplot/flags/tag_contigs.{wildcards.hap}.{region}.done" for region in bed_df.index ]

def find_pq_tigs(wildcards):
    return [ f"results/{wildcards.sample}/moddotplot/work/find_tigs/beds/{wildcards.hap}/{region}_pq_contig.bed" for region in bed_df.index ]

def find_p_tigs(wildcards): 
    return [ f"results/{wildcards.sample}/moddotplot/work/find_tigs/beds/{wildcards.hap}/{region}_p_contig.bed" for region in bed_df.index ]

def get_all_flag(wildcards):
    return [ f"results/{wildcards.sample}/moddotplot/flags/pq_selfplot.{wildcards.hap}.{region}.done" for region in bed_df.index]

rule summarize_moddot_results:
    input:
        get_all_flag
    output:
        tsv = "results/{sample}/moddotplot/outputs/summary/{hap}.generated_acros.tsv"
    resources:
        mem=4,
        hrs=1,
    threads: 1
    run:
        header = ["Haplotype"]+acros
        called = [[wildcards.sample]]
        plot_dir=f"results/{wildcards.sample}/moddotplot/outputs/pq_selfplots/{wildcards.hap}"
        for chrom in acros:
            region_name = bed_df[bed_df["chr"] == chrom].index[0]
            pdf_files = glob.glob(f"{plot_dir}/{region_name}_*_FULL.pdf")
            if len(pdf_files) > 0:
                contig_names = [os.path.basename(pdf).replace(f"{region_name}_","").replace("_FULL.pdf","") for pdf in pdf_files]
                called[0].append(",".join(contig_names))
            else:
                called[0].append("NA")
        result_df = pd.DataFrame(called, columns = header)
        result_df.to_csv(output.tsv, sep="\t", index=False)
                

checkpoint subset_target_region:
    output:
        bed="resources/acro_target_beds/{region}.bed",
    resources:
        mem=4,
        hrs=2,
    threads: 1
    run:
        out_df = bed_df[bed_df["chr"] == str(wildcards.region).split("_")[0]]
        out_df[["chr", "start", "end"]].to_csv(output.bed, sep="\t", header=False, index=False)

rule liftover:
    input:
        bed=rules.subset_target_region.output.bed,
        paf="results/{sample}/saffire/outputs/trimmed_pafs/CHM13/{hap}.minimap2.trimmed.paf"
    output:
        paf= "results/{sample}/moddotplot/work/liftover/CHM13/pafs/{hap}/{region}.paf",
        flag = touch("results/{sample}/moddotplot/flags/liftover.{hap}.{region}.done"),
    resources:
        mem=lambda wildcards, attempt: 48 * attempt,
        hrs=2,
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    shell: """
        set -euo pipefail    
        rustybam liftover --bed {input.bed} {input.paf} > {output.paf}
        """

rule trim_paf_moddot:
    input:
        paf= rules.liftover.output.paf
    output:
        paf = "results/{sample}/moddotplot/work/liftover/CHM13/trimmed_pafs/{hap}/{region}.paf",
        flag = "results/{sample}/moddotplot/flags/trim_paf_moddot.{hap}.{region}.done"
    threads: 8
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    resources:
        mem = 12,
        hrs = 48
    shell: """
        set -euo pipefail    
        rustybam trim-paf {input.paf} > {output.paf}
        touch {output.flag}
        """

rule paf_stats:
    input:
        paf=rules.trim_paf_moddot.output.paf
    output:
        stats = "results/{sample}/moddotplot/work/liftover/CHM13/paf_stats/{hap}/{region}.stats",
    resources:
        mem=12,
        hrs=2,
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    shell: """
        set -euo pipefail    
        tmp_output="{output.stats}.tmp"
        rustybam stats --paf --qbed {input.paf} > $tmp_output
        mv $tmp_output {output.stats}
        """

rule tag_contigs:
    input:
        paf = rules.trim_paf_moddot.output.paf,
        paf_stats = rules.paf_stats.output.stats
    output:
        flag = "results/{sample}/moddotplot/flags/tag_contigs.{hap}.{region}.done"
    params:
        pq = "results/{sample}/moddotplot/work/find_tigs/beds/{hap}/{region}_pq_contig.bed",
        p = "results/{sample}/moddotplot/work/find_tigs/beds/{hap}/{region}_p_contig.bed",
        aln = "results/{sample}/moddotplot/work/find_tigs/pafs/{hap}/{region}_contig_aln.paf",
    threads: 1,
    resources:
        mem=10,
        hrs=24,
        disk_free=1
    script:
        f"{SNAKEMAKE_DIR}/scripts/tag_contigs.py"

rule get_pq_fa:
    input:
        flag = rules.tag_contigs.output.flag,
        hap = rules.rename_fasta.output.final_fasta
    output: 
        flag = "results/{sample}/moddotplot/flags/get_pq_fa.{hap}.{region}.done",
    params:
        fa = "results/{sample}/moddotplot/work/get_pq_tigs/fasta/{hap}/{region}.pq_contig.fa",
        bed = "results/{sample}/moddotplot/work/find_tigs/beds/{hap}/{region}_pq_contig.bed",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    run:
        os.makedirs(f"results/{wildcards.sample}/moddotplot/work/get_pq_tigs/fasta/{wildcards.hap}", exist_ok=True)

        input_bed = params.bed
        input_hap = pysam.FastaFile(input.hap)
        output_fa = params.fa
        region_name = wildcards.region
        output_fa_data = []
        if (os.path.isfile(input_bed)) and (os.path.getsize(input_bed) > 0):
            with open(input_bed) as finp_bed:
                finp_bed_lines = finp_bed.read().strip().split("\n")
            for bed_line in finp_bed_lines:
                token = bed_line.split("\t")
                contig_name, start, end = token[:3]
                start, end = int(start)+1, int(end)
                output_fa_seq_name = f">{region_name}_{contig_name}"
                subseq = input_hap.fetch(contig_name, start, end)
                output_fa_data += [output_fa_seq_name, subseq]
        with open(output_fa,"w") as fout:
            fout.write("\n".join(output_fa_data))
        if len(output_fa_data) > 0:
            shell(f"samtools faidx {output_fa}")
        with open(str(output.flag), "w") as f_flag:
            print(f"{wildcards.sample}:{wildcards.hap}:{wildcards.region} Done.", file=f_flag)

rule pq_selfplot:
    input:
        tag_contig_flag = rules.tag_contigs.output.flag,
        get_pq_fa_flag = rules.get_pq_fa.output.flag,
    output:
        flag = "results/{sample}/moddotplot/flags/pq_selfplot.{hap}.{region}.done"
    params:
        fa = "results/{sample}/moddotplot/work/get_pq_tigs/fasta/{hap}/{region}.pq_contig.fa",
        bed = "results/{sample}/moddotplot/work/find_tigs/beds/{hap}/{region}_pq_contig.bed"
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1,
    singularity:
        "docker://eichlerlab/moddotplot:0.9.0"
    shell: """
        set -euo pipefail    
        outdir=$(dirname {output.flag} | sed "s/flags/outputs\/pq_selfplots\/{wildcards.hap}/g")
        if [[ ! -d $outdir ]];then
            mkdir -p $outdir
        fi
        if [[ -s {params.bed} ]]; then
            moddotplot static -f {params.fa} --no-hist --no-bed -o $outdir
            echo "{wildcards.sample}:{wildcards.hap}:{wildcards.region}" > {output.flag}
        else
            echo "{wildcards.sample}:{wildcards.hap}:{wildcards.region}" > {output.flag}
        fi
        """

rule group_tigs:
    input:
        flag = find_tigs_flag,
        paf = "results/{sample}/saffire/outputs/trimmed_pafs/CHM13/{hap}.minimap2.trimmed.paf",
    output:
        tigs = "results/{sample}/moddotplot/outputs/contig_stats/{hap}.CHM13_lifted_contigs.tsv",
        tab = "results/{sample}/moddotplot/outputs/contig_stats/{hap}.CHM13_lifted_contigs.tab"
    params:
        p_tigs = find_p_tigs,
        pq_tigs = find_pq_tigs,
    threads: 1
    resources:
        mem=10,
        hrs=24,
        disk_free=1

    run:
        sample = wildcards.sample
        input_pq_tigs = params.pq_tigs
        input_p_tigs = params.p_tigs
        input_paf = input.paf
         
        acros = ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']
        fout = open(output.tigs, 'w')
        print('Sample\tArm\tTig\tTig_start\tTig_end\tSize\tChrom',file=fout)
        haps = open(output.tab, 'w')
        print('SAMPLE\tHAP\tBED', file=haps)
        
        qry_chrom_tracker = {}
        qrylen_dict = {}
        with open(input_paf) as f:
            lines = f.read().strip().split("\n")
        for line in lines:
            query, query_len, query_start, query_end, strand, target, target_len, target_start, target_end, num_matches, alignment_len = line.strip().split("\t")[:11]
            qrylen_dict[query] = query_len
            if int(alignment_len) < 200000:
                continue
            if query not in qry_chrom_tracker:
                qry_chrom_tracker[query] = []
            if not target in qry_chrom_tracker[query]:
                qry_chrom_tracker[query].append(target)
        for pq_tig in input_pq_tigs:
            if not os.path.getsize(pq_tig)>0:
                continue
            region = "_".join(os.path.basename(pq_tig).split("_")[:3])
            with open(pq_tig) as finp_pq_tig:
                pq_tig_lines = finp_pq_tig.read().strip().split("\n")
            for pq_tig_line in pq_tig_lines:
                entries = pq_tig_line.strip().split('\t')
                qlen = qrylen_dict[entries[0]]
                print(f"{sample}\tpq\t{entries[0]}\t{entries[1]}\t{entries[2]}\t{qlen}\t{region}", file=fout)
                print(f"{sample}\t{wildcards.hap}\t{pq_tig}", file=haps)
        for p_tig in input_p_tigs:
            if not os.path.getsize(p_tig)>0:
                continue
            region = "_".join(os.path.basename(p_tig).split("_")[:3])
            with open(p_tig) as finp_p_tig:
                p_tig_lines = finp_p_tig.read().strip().split("\n")
            for p_tig_line in p_tig_lines:
                entries = p_tig_line.strip().split('\t')
                qlen = qrylen_dict[entries[0]]
                not_acro = False
                try:
                    for ele in qry_chrom_tracker[entries[0]]:
                        if ele not in acros:
                            not_acro = True
                            break
                    if not not_acro:
                        print(f"{sample}\tp\t{entries[0]}\t{entries[1]}\t{entries[2]}\t{qlen}\t{region}", file=fout)
                except KeyError:
                    continue

        fout.close()
        haps.close()
