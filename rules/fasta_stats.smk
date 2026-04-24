import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio import SeqIO
import pysam
import glob

SNAKEMAKE_ROOT_DIR = os.path.dirname(workflow.snakefile).replace("/rules","")

# NOTE
# The whole-genome–level summary was added because of the SMAHT DSA metadata requirements.
# As the number of required outputs increased, the rule ended up growing unintentionally, it created some redundant steps.
# You might consider calculating all the statistics directly from the original scaffold-level FASTA instead of splitting it.
#-----------------------------------------


def find_fasta_set_for_full_genome(which_one):
    def inner(wildcards):
        
        try:
            h2_val = full_manifest_df.at[wildcards.sample, "H2"]
        except KeyError:
            h2_val = None 

        has_h2 = isinstance(h2_val, str) and h2_val.strip() != "" and os.path.isfile(h2_val)
        if which_one == "scaffold":
            fasta_type = "final"
        elif which_one == "contig":
            fasta_type = "contig"
        haps = ["hap1", "hap2"] if has_h2 else ["hap1"]

        outputs = [
            f"results/{wildcards.sample}/contamination_screening/outputs/{fasta_type}_fasta/{wildcards.sample}_{hap}.fasta"
            for hap in haps
        ]
        if INCLUDE_MITO and (int(TAXID) == 9606): # try to find and add mt fasta for only human as set INCLUDE_MITO as true.
            mito_fasta_candidates = glob.glob(f"results/{wildcards.sample}/contamination_screening/outputs/mito_fasta/*mt.fasta")
            outputs += [ mito_fasta_candidate for mito_fasta_candidate in mito_fasta_candidates if (os.path.getsize(mito_fasta_candidate) > 0) ]

        return outputs[0] if len(outputs) == 1 else outputs
    return inner


def find_hap2_and_un_for_nucfreq(wildcards):
    
    try:
        h2_val = full_manifest_df.at[wildcards.sample, "H2"]
    except KeyError:
        h2_val = None 

    try:
        un_val = full_manifest_df.at[wildcards.sample, "UNASSIGNED"]
    except KeyError:
        un_val = None 


    has_h2 = isinstance(h2_val, str) and h2_val.strip() != "" and os.path.isfile(h2_val)
    has_un = isinstance(un_val, str) and un_val.strip() != "" and os.path.isfile(un_val)

    haps = []
    if has_h2:
        haps.append("hap2")
    if has_un:
        haps.append("un")

    outputs = [
        f"results/{wildcards.sample}/contamination_screening/outputs/final_fasta/{wildcards.sample}_{hap}.fasta"
        for hap in haps
    ]
    return outputs[0] if len(outputs) == 1 else outputs


def find_all_hap_summary(wildcards):
    try:
        h2_val = full_manifest_df.at[wildcards.sample, "H2"]
    except KeyError:
        h2_val = None 

    has_h2 = isinstance(h2_val, str) and h2_val.strip() != "" and os.path.isfile(h2_val)
    haps = ["hap1", "hap2"] if has_h2 else ["hap1"]
    outputs = [
        f"results/{wildcards.sample}/stats/outputs/summary_by_hap/{hap}.summary.stats"
        for hap in haps
    ]
    return outputs[0] if len(outputs) == 1 else outputs

def get_stats_only_inputs():
    outputs = []
    for sample in samples_with_asm:
        result = find_all_hap_summary(type("WC", (), {"sample": sample})())
        if isinstance(result, str):
            outputs.append(result)
        else:
            outputs.extend(result)
    return outputs

rule split_scaffolds:
    input:
        scaffold_fasta=rules.rename_fasta.output.final_fasta,
    output:
        contig_fasta = "results/{sample}/contamination_screening/outputs/contig_fasta/{sample}_{hap}.fasta"
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: attempt * 8,
        hrs=12,
    run:

        contig_fasta = output.contig_fasta
        contig_fasta_records = list()
        gap_pattern = re.compile(r"N+")
        fasta_records = list(SeqIO.parse(input.scaffold_fasta, "fasta"))

        for record in fasta_records:

            seq_name = record.id
            raw_sequence = str(record.seq)
            gaps = list(gap_pattern.finditer(raw_sequence))

            non_gap_coords = []
            prev_gap_end = 0

            if len(gaps) == 0:
                contig_fasta_records.append(SeqRecord(Seq(raw_sequence), id=seq_name, description=""))
                continue

            for gap in gaps:
                g_start, g_end = gap.start(), gap.end() # g_start; 0-based
                g_start += 1 # converted 1-based
                if g_start == 1:
                    prev_gap_end = g_end
                    continue
                non_gap_coords.append((prev_gap_end+1, g_start-1)) # 1-based coords for non-gap segment
                prev_gap_end = g_end    
            if prev_gap_end < len(raw_sequence):
                non_gap_coords.append((prev_gap_end+1, len(raw_sequence)))

            for non_gap_coord in non_gap_coords:
                non_gap_start, non_gap_end = non_gap_coord
                contig_name = f"{seq_name}_sub_{non_gap_start}-{non_gap_end}"
                contig_sequence = raw_sequence[non_gap_start-1:non_gap_end]
                contig_fasta_records.append(SeqRecord(Seq(contig_sequence), id=contig_name, description=""))

        with open(contig_fasta, "w") as fout:
            fasta_writer = FastaWriter(fout, wrap=None)
            fasta_writer.write_file(contig_fasta_records)
        pysam.faidx(contig_fasta)


rule combine_full_genome:
    input:
        all_fasta = find_fasta_set_for_full_genome(which_one="scaffold")
    output:
        full_genome = "results/{sample}/stats/work/full_genome/{sample}.fasta",
        full_genome_fai = "results/{sample}/stats/work/full_genome/{sample}.fasta.fai",
    threads: 1,
    resources:
        hrs=1,
        mem=8,
    run:
        if isinstance(input.all_fasta, str):
            input_all_fasta = [input.all_fasta]
        else:
            input_all_fasta = input.all_fasta

        records = []
        for fasta in input_all_fasta:
            records.extend(list(SeqIO.parse(fasta, "fasta")))
        output_full_genome = output.full_genome
        outdir = os.path.dirname(output_full_genome)
        os.makedirs(outdir, exist_ok=True)
        with open(output_full_genome, "w") as fout:
            writer = FastaWriter(fout, wrap=None)
            writer.write_file(records)
        pysam.faidx(output_full_genome)

rule combine_hap2_and_unassigned_fasta:
    input:
        hap2_and_un_fasta = find_hap2_and_un_for_nucfreq
    output:
        flag = "results/{sample}/contamination_screening/flags/combine_hap2_and_unassigned_fasta.done"
    params:
        hap_and_un_fasta = "results/{sample}/contamination_screening/outputs/final_fasta/{sample}_hap2_with_un.fasta"
    threads: 1,
    resources:
        hrs=1,
        mem=8,
    run:
        output_merged_fasta = params.hap_and_un_fasta
        if isinstance(input.hap2_and_un_fasta, str):
            input_hap2_and_un_fasta = [input.hap2_and_un_fasta]
        else:
            input_hap2_and_un_fasta = input.hap2_and_un_fasta
        out_flag = output.flag
        outdir = os.path.dirname(out_flag)
        os.makedirs(outdir, exist_ok=True)        
        if len(input_hap2_and_un_fasta) == 2:
            records = []
            for fasta in input_hap2_and_un_fasta:
                records.extend(list(SeqIO.parse(fasta, "fasta")))
            with open(output_merged_fasta, "w") as fout_fasta:
                writer = FastaWriter(fout_fasta, wrap=None)
                writer.write_file(records)
            pysam.faidx(output_merged_fasta)
        with open(out_flag,"w") as fout:
            print (f"{wildcards.sample}.done", file=fout)

rule get_assebmly_eval_config:        
    input:
        rules.combine_hap2_and_unassigned_fasta.output.flag
    output:
        assembly_eval_config = "results/{sample}/assembly_eval_config/output/config_file/{sample}.config.yaml"
    threads: 1,
    resources:
        hrs=1,
        mem=4,
    run:
        final_fasta_files = glob.glob(f"results/{wildcards.sample}/contamination_screening/outputs/final_fasta/{wildcards.sample}_*.fasta")
        fasta_dict = dict()
        for fasta in final_fasta_files:
            fasta_hap = os.path.basename(fasta).replace(f"{wildcards.sample}_","").split(".")[0]
            fasta_path = os.path.abspath(fasta)
            fasta_dict[fasta_hap] = fasta_path
        hap1_fasta = fasta_dict["hap1"]
        if "hap2_with_un" in fasta_dict:
            hap2_fasta = fasta_dict["hap2_with_un"]
        else:
            if "hap2" in fasta_dict:
                hap2_fasta = fasta_dict["hap2"]
            else:
                hap2_fasta = hap1_fasta
        with open(output.assembly_eval_config, "w") as fout:
            print (f"""
{wildcards.sample}:
  fofns:
    HiFi: HIFI_FOFN
  type_map:
    HiFi: winnowmap
  nuc_opts: -y 100
  asm_h1: {hap1_fasta}
  asm_h2: {hap2_fasta}
  repeat_mask: False
  species: Human
                """, file=fout
            )            

rule combine_full_genome_contigs:
    input:
        all_fasta = find_fasta_set_for_full_genome(which_one="contig")
    output:
        full_genome = "results/{sample}/stats/work/full_genome/contig_fasta/{sample}.fasta",
        full_genome_fai = "results/{sample}/stats/work/full_genome/contig_fasta/{sample}.fasta.fai",
    threads: 1,
    resources:
        hrs=1,
        mem=8,
    run:
        if isinstance(input.all_fasta, str):
            input_all_fasta = [input.all_fasta]
        else:
            input_all_fasta = input.all_fasta

        records = []
        for fasta in input_all_fasta:
            records.extend(list(SeqIO.parse(fasta, "fasta")))
        output_full_genome = output.full_genome
        outdir = os.path.dirname(output_full_genome)
        os.makedirs(outdir, exist_ok=True)
        with open(output_full_genome, "w") as fout:
            writer = FastaWriter(fout, wrap=None)
            writer.write_file(records)
        pysam.faidx(output_full_genome)

rule get_telo_stats:
    input:
        scaffold_fasta = rules.rename_fasta.output.final_fasta,
        contig_fasta = rules.split_scaffolds.output.contig_fasta
    output:
        scaffold_telo_tbl = "results/{sample}/stats/work/telo/{hap}.scaffold.telo.tbl",
        contig_telo_tbl = "results/{sample}/stats/work/telo/{hap}.contig.telo.tbl",
        flag = touch("results/{sample}/stats/flags/get_telo_stats.{hap}.done")
    singularity:
        "docker://eichlerlab/binf-basics:0.1",
    threads: 1,
    resources:
        hrs=1,
        mem=6,
    shell: """
        set -euo pipefail
        echo -e "seq_name\tstart\tend\tseq_length" > {output.scaffold_telo_tbl} && seqtk telo {input.scaffold_fasta} >> {output.scaffold_telo_tbl}
        echo -e "seq_name\tstart\tend\tseq_length" > {output.contig_telo_tbl} && seqtk telo {input.contig_fasta} >> {output.contig_telo_tbl}
        """

rule get_full_genome_telo_stats:
    input:
        full_genome_fasta = rules.combine_full_genome.output.full_genome,
        full_genome_contig_fasta = rules.combine_full_genome_contigs.output.full_genome
    output:
        telo_tbl = "results/{sample}/stats/work/telo/full_genome.telo.tbl",
        contig_telo_tbl = "results/{sample}/stats/work/telo/full_genome_contigs.telo.tbl",
        flag = touch("results/{sample}/stats/flags/get_full_genome_telo_stats.done")
    singularity:
        "docker://eichlerlab/binf-basics:0.1",
    threads: 1,
    resources:
        hrs=1,
        mem=8,
    shell: """
        set -euo pipefail
        echo -e "seq_name\tstart\tend\tseq_length" > {output.telo_tbl} && seqtk telo {input.full_genome_fasta} >> {output.telo_tbl}
        echo -e "seq_name\tstart\tend\tseq_length" > {output.contig_telo_tbl} && seqtk telo {input.full_genome_contig_fasta} >> {output.contig_telo_tbl}
        """    

rule get_contig_stats:
    input:
        fasta = rules.split_scaffolds.output.contig_fasta,
        telo_tbl = rules.get_telo_stats.output.contig_telo_tbl
    output:
        stats = "results/{sample}/stats/work/fasta_stats/{hap}.contig.stats",
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=4,
    script:
        f"{SNAKEMAKE_ROOT_DIR}/scripts/fasta_stats.py"

rule get_scaffold_stats:
    input:
        fasta = rules.rename_fasta.output.final_fasta,
        telo_tbl = rules.get_telo_stats.output.scaffold_telo_tbl,
    output:
        stats = "results/{sample}/stats/work/fasta_stats/{hap}.scaffold.stats",
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=4,
    script:
        f"{SNAKEMAKE_ROOT_DIR}/scripts/fasta_stats.py"

rule get_full_genome_stats:
    input:
        fasta = rules.combine_full_genome.output.full_genome,
        telo_tbl = rules.get_full_genome_telo_stats.output.telo_tbl
    output:
        stats = "results/{sample}/stats/work/fasta_stats/full_genome.stats"
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=4,
    script:
        f"{SNAKEMAKE_ROOT_DIR}/scripts/fasta_stats.py"

rule get_full_genome_contig_stats:
    input:
        fasta = rules.combine_full_genome_contigs.output.full_genome,
        telo_tbl = rules.get_full_genome_telo_stats.output.contig_telo_tbl
    output:
        stats = "results/{sample}/stats/work/fasta_stats/full_genome.contig.stats"
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=4,
    script:
        f"{SNAKEMAKE_ROOT_DIR}/scripts/fasta_stats.py"

rule summary_hap_stats:
    input:
        scaffold_stats = rules.get_scaffold_stats.output.stats,
        contig_stats = rules.get_contig_stats.output.stats,
        covered_chrom_tsv = "results/{sample}/saffire/outputs/chrom_cov/CHM13/{hap}.minimap2.chrom_cov.tsv",
        contigs_covered_chrom_tsv = "results/{sample}/saffire/outputs/chrom_cov/CHM13/{hap}.minimap2.contigs_chrom_cov.tsv",
        sample_qv = rules.merqury_run.output.qv,
        busco_result = rules.summarize_compleasm_results.output.summary,
    output:
        hap_summary = "results/{sample}/stats/outputs/summary/{hap}.summary.stats",
    params:
        un_qv = "results/{sample}/merqury/outputs/{sample}_un.qv"
    threads: 1,
    resources:
        hrs=1,
        mem=4,
    run:
        scaffold_stats = input.scaffold_stats
        contig_stats = input.contig_stats
        covered_chrom_tsv = input.covered_chrom_tsv
        contigs_covered_chrom_tsv = input.contigs_covered_chrom_tsv
        df_scaffold = pd.read_csv(scaffold_stats, sep="\t").add_suffix("_scaffold")
        df_scaffold = df_scaffold.rename(columns = {"sample_scaffold":"sample", "haplotype_scaffold":"haplotype"})
        
        df_contig = pd.read_csv(contig_stats, sep="\t").add_suffix("_contig")
        df_contig = df_contig.rename(columns = {"sample_contig":"sample", "haplotype_contig":"haplotype"})

        df = pd.merge(df_scaffold, df_contig, on=["sample","haplotype"], how="outer")

        assigned_chrom_df = pd.read_csv(contigs_covered_chrom_tsv, sep="\t", usecols=["contig", "best_target_chrom"])
        assigned_chrom_df = assigned_chrom_df[assigned_chrom_df["best_target_chrom"] != "undefined"].copy()
        count_df = (
            assigned_chrom_df.groupby("best_target_chrom")["contig"]
            .nunique()
            .reset_index(name="n_contigs")
        )
        count_df["n_gaps"] = count_df["n_contigs"] - 1


        if wildcards.hap == "un":
            with open(params.un_qv) as f_un_qv:
                token = f_un_qv.read().strip().split("\n")[-1].split("\t")
                qv = float(token[3])
        else: # hap1 or hap2
            qv_df = pd.read_csv(input.sample_qv, sep="\t", header=None, names=["haplotype","error_count", "total_count", "qv","error_rate"]).drop_duplicates().set_index("haplotype")
            qv_df = qv_df[pd.to_numeric(qv_df["qv"], errors="coerce").notna()]
            qv = float(qv_df.loc[f"{wildcards.sample}_{wildcards.hap}","qv"])
            
        df["quality_value"] = qv

        covered_chrom_df = pd.read_csv(input.covered_chrom_tsv, sep="\t").set_index("chrom")
        df["chrX_ratio"] = covered_chrom_df.loc["chrX", "covered_pct"]
        df["chrY_ratio"] = covered_chrom_df.loc["chrY", "covered_pct"]

        df["gaps_between_scaffolds"] = count_df["n_gaps"].sum()



        busco_df = pd.read_csv(input.busco_result, sep="\t").set_index(["ASSEMBLY","HAPLOTYPE"])
        busco_df["percent_single_copy"] = (busco_df["S"]/busco_df["N"]*100).round(4)
        busco_df["percent_multi_copy"] = (busco_df["D"]/busco_df["N"]*100).round(4)
        busco_df["percent_fragmented_copy"] = (busco_df["F"]/busco_df["N"]*100).round(4)
        busco_df["percent_missing_copy"] = (busco_df["M"]/busco_df["N"]*100).round(4)
        
        df["percent_single_copy"] = busco_df.loc[(wildcards.sample, wildcards.hap), "percent_single_copy"]
        df["percent_multi_copy"] = busco_df.loc[(wildcards.sample, wildcards.hap), "percent_multi_copy"]
        df["percent_fragmented_copy"] = busco_df.loc[(wildcards.sample, wildcards.hap), "percent_fragmented_copy"]
        df["percent_missing_copy"] = busco_df.loc[(wildcards.sample, wildcards.hap), "percent_missing_copy"]
        df["number_of_contigs"] = df["total_num_of_contigs_scaffold"] - df["num_of_contigs_with_N_scaffold"]
        df["number_of_near_T2T_scaffolds"] = df["num_of_near_t2t_contigs_scaffold"] - df["num_of_near_t2t_contigs_contig"]

        df = df.rename(columns = {
                "num_of_contigs_with_N_scaffold":"number_of_scaffolds",
                "num_of_near_t2t_contigs_contig":"number_of_near_T2T_contigs",
                "l50_contig":"contig_l50",
                "l50_scaffold":"scaffold_l50",
                "n50_contig":"contig_n50",
                "n50_scaffold":"scaffold_n50",
                "aun_contig":"contig_aun",
                "aun_scaffold":"scaffold_aun",
                "bases_contig":"total_ungapped_length",
                "bases_scaffold":"genome_size",
                "gc_content_contig":"gc_content",
                "largest_size_contig":"largest_contig_size",
            }
        )
        selected_columns = ["chrX_ratio","chrY_ratio","genome_size","quality_value",\
        "number_of_contigs", "number_of_scaffolds","gaps_between_scaffolds",\
        "number_of_near_T2T_contigs", "number_of_near_T2T_scaffolds", "contig_l50","contig_n50","scaffold_l50", "scaffold_n50","total_ungapped_length", "largest_contig_size","gc_content",\
        "contig_aun", "scaffold_aun", \
        "percent_single_copy","percent_multi_copy","percent_fragmented_copy","percent_missing_copy"
        ]
        df = df[selected_columns]
        df.to_csv(output.hap_summary, sep="\t", index=False)

rule smaht_hap_metadata:
    input:
        scaffold_stats = rules.get_scaffold_stats.output.stats,
        contig_stats = rules.get_contig_stats.output.stats,
        contigs_covered_chrom_tsv = "results/{sample}/saffire/outputs/chrom_cov/CHM13/{hap}.minimap2.contigs_chrom_cov.tsv",
        sample_qv = rules.merqury_run.output.qv,
        busco_result = rules.summarize_compleasm_results.output.summary,
    output:
        hap_metadata = "results/{sample}/smaht_dsa_metadata/outputs/{sample}.{hap}.DSA_ExternalQualityMetric.tsv",
    params:
        un_qv = "results/{sample}/merqury/outputs/{sample}_un.qv"
    threads: 1,
    resources:
        hrs=1,
        mem=4,
    run:
        scaffold_stats = input.scaffold_stats
        contig_stats = input.contig_stats
        contigs_covered_chrom_tsv = input.contigs_covered_chrom_tsv
        df_scaffold = pd.read_csv(scaffold_stats, sep="\t").add_suffix("_scaffold")
        df_scaffold = df_scaffold.rename(columns = {"sample_scaffold":"sample", "haplotype_scaffold":"haplotype"})
        
        df_contig = pd.read_csv(contig_stats, sep="\t").add_suffix("_contig")
        df_contig = df_contig.rename(columns = {"sample_contig":"sample", "haplotype_contig":"haplotype"})

        df = pd.merge(df_scaffold, df_contig, on=["sample","haplotype"], how="outer")
        if wildcards.hap == "un":
            with open(params.un_qv) as f_un_qv:
                token = f_un_qv.read().strip().split("\n")[-1].split("\t")
                qv = float(token[3])
        else: # hap1 or hap2
            qv_df = pd.read_csv(input.sample_qv, sep="\t", header=None, names=["haplotype","error_count", "total_count", "qv","error_rate"]).drop_duplicates().set_index("haplotype")
            qv_df = qv_df[pd.to_numeric(qv_df["qv"], errors="coerce").notna()]
            qv = float(qv_df.loc[f"{wildcards.sample}_{wildcards.hap}","qv"])
            
        df["quality_value"] = qv
        
        assigned_chrom_df = pd.read_csv(contigs_covered_chrom_tsv, sep="\t", usecols=["contig", "best_target_chrom"])
        assigned_chrom_df = assigned_chrom_df[assigned_chrom_df["best_target_chrom"] != "undefined"].copy()
        count_df = (
            assigned_chrom_df.groupby("best_target_chrom")["contig"]
            .nunique()
            .reset_index(name="n_contigs")
        )
        count_df["n_gaps"] = count_df["n_contigs"] - 1
        df["gaps_between_scaffolds"] = count_df["n_gaps"].sum()

        busco_df = pd.read_csv(input.busco_result, sep="\t").set_index(["ASSEMBLY","HAPLOTYPE"])
        busco_df["percent_single_copy"] = (busco_df["S"]/busco_df["N"]*100).round(4)
        busco_df["percent_multi_copy"] = (busco_df["D"]/busco_df["N"]*100).round(4)
        busco_df["percent_fragmented_copy"] = (busco_df["F"]/busco_df["N"]*100).round(4)
        busco_df["percent_missing_copy"] = (busco_df["M"]/busco_df["N"]*100).round(4)
        
        df["percent_single_copy"] = busco_df.loc[(wildcards.sample, wildcards.hap), "percent_single_copy"]
        df["percent_multi_copy"] = busco_df.loc[(wildcards.sample, wildcards.hap), "percent_multi_copy"]
        df["percent_fragmented"] = busco_df.loc[(wildcards.sample, wildcards.hap), "percent_fragmented_copy"]
        df["percent_missing"] = busco_df.loc[(wildcards.sample, wildcards.hap), "percent_missing_copy"]
        df["number_of_contigs"] = df["total_num_of_contigs_scaffold"] - df["num_of_contigs_with_N_scaffold"]
        df["submitted_id"] = f"UWSC_EXTERNAL-QUALITY-METRIC_{wildcards.sample.upper()}-{wildcards.hap.upper()}"
        df["category"] = "Donor Specific Assembly"
        df["number_of_chromosomes"] = assigned_chrom_df["best_target_chrom"].nunique()
        
        df = df.rename(columns = {
                "num_of_contigs_with_N_scaffold":"number_of_scaffolds",
                "l50_contig":"contig_l50",
                "l50_scaffold":"scaffold_l50",
                "n50_contig":"contig_n50",
                "n50_scaffold":"scaffold_n50",
                "bases_contig":"total_ungapped_length",
                "gc_content_contig":"gc_content",
                "largest_size_contig":"largest_contig_size",
            }
        )
        selected_columns = [
            "submitted_id","category",\
            "contig_l50","contig_n50","gaps_between_scaffolds","gc_content","largest_contig_size","number_of_contigs","number_of_scaffolds",\
            "percent_single_copy","percent_multi_copy","percent_fragmented","percent_missing",\
            "quality_value",\
            "scaffold_l50","scaffold_n50","total_ungapped_length","number_of_chromosomes"
            ]
        df = df[selected_columns]
        df.to_csv(output.hap_metadata, sep="\t", index=False)
