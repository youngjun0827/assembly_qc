from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import pysam
import numpy as np
import re
import gzip
from collections import defaultdict

FCS_SOURCE_DIR = "/net/eichler/vol28/7200/software/pipelines/foreign_contamination_screen"
MITO_DB = config.get("mito", f"{FCS_SOURCE_DIR}/db/mito_ebv.fa")
RDNA_DB = config.get("rdna", f"{FCS_SOURCE_DIR}/db/rdna.fa")

manifest_df = conv_manifest_df.copy()

#-----------------------------------------


def get_raw_fasta(wildcards):
    return manifest_df.loc[(wildcards.sample, wildcards.hap), "FASTA"]

def find_gx_report(wildcards):
    (IDS,) = glob_wildcards("results/{sample}/contamination_screening/work/fcs_gx/outputs/{hap}/{gx_name}.fcs_gx_report.txt".format(
            sample=wildcards.sample,
            hap=wildcards.hap, 
            gx_name="{gx_name}"
        )
    )

    return expand("results/{sample}/contamination_screening/work/fcs_gx/outputs/{hap}/{gx_name}.fcs_gx_report.txt",
        sample=wildcards.sample,
        hap=wildcards.hap,
        gx_name=IDS,
    )


wildcard_constraints:
    sample = "|".join(sample_list),
    sub = "|".join(["gx", "adaptor"]),


localrules:
    run_fcs_gx,


rule index_fasta:
    input:
        fasta = get_raw_fasta
    output:
        link_fasta = "results/{sample}/contamination_screening/work/temp/raw_fasta/{hap}.fasta",
        fai = "results/{sample}/contamination_screening/work/temp/raw_fasta/{hap}.fasta.fai"
    threads: 1
    resources:
        mem = 8,
        hrs = 4,
    shell: """
        ln -s {input.fasta} {output.link_fasta}
        samtools faidx {output.link_fasta}
        """


rule remove_short_contigs:
    input:
        fasta = rules.index_fasta.output.link_fasta,
        fai = rules.index_fasta.output.fai
    output:
        filtered_fasta = "results/{sample}/contamination_screening/work/fastq_cleaning/short_contig_flitering/cleaned_fasta/{hap}.filtered.fasta",
        filtered_fai = "results/{sample}/contamination_screening/work/fastq_cleaning/short_contig_flitering/cleaned_fasta/{hap}.filtered.fasta.fai"
    threads: 1
    resources:
        mem = 24,
        hrs = 2,
    run:
        filtered_fasta = output.filtered_fasta
        below_ten_fasta = filtered_fasta.replace(".filtered.fasta",".below_10bp.fasta")

        try:
            raw_fasta_records = list(SeqIO.parse(input.fasta, "fasta"))
        except UnicodeDecodeError:
            with gzip.open(input.fasta, "rt") as input_fasta:
                raw_fasta_records = list(SeqIO.parse(input_fasta, "fasta"))
        filtered_fasta_records = list()
        below_ten_fasta_records = list()

        for raw_fasta_record in raw_fasta_records:
            contig_name = str(raw_fasta_record.id)
            contig_seq = str(raw_fasta_record.seq)
            if len(contig_seq) < 10:
                below_ten_fasta_records.append(SeqRecord(Seq(contig_seq), id=contig_name, description=""))
            else:
                filtered_fasta_records.append(SeqRecord(Seq(contig_seq), id=contig_name, description=""))

        with open(filtered_fasta, "w") as fout_f:
            fasta_writer = FastaWriter(fout_f, wrap=None)
            fasta_writer.write_file(filtered_fasta_records)
        pysam.faidx(filtered_fasta) # indexing
        if len(below_ten_fasta_records) > 0:
            with open(below_ten_fasta, "w") as fout_b:
                fasta_writer = FastaWriter(fout_b, wrap=None)
                fasta_writer.write_file(below_ten_fasta_records)
            pysam.faidx(below_ten_fasta) # indexing


rule run_fcs_gx:
    input:
        filt_fasta = rules.remove_short_contigs.output.filtered_fasta
    output:
        dir = directory("results/{sample}/contamination_screening/work/fcs_gx/outputs/{hap}"),
        flag = "results/{sample}/contamination_screening/flags/run_fcs_gx.{hap}.done",
    threads: 1
    resources:
        mem = 16,
        hrs = 12,
        load = 25,
    params:
        taxid=TAXID,
        GXDB_LOC="/data/scratch/GXDB/gxdb/",
        fcs_img=f"{FCS_SOURCE_DIR}/images/fcs-gx.sif",
        fcs_script=f"{FCS_SOURCE_DIR}/fcs.py",
    shell: """
        set -euo pipefail  
        python3 {params.fcs_script} --image {params.fcs_img} screen genome --fasta {input.filt_fasta} --out-dir {output.dir} --gx-db {params.GXDB_LOC}  --tax-id {params.taxid}
        touch {output.flag}
        """


rule run_fcs_adapter:
    input:
        filt_fasta = rules.remove_short_contigs.output.filtered_fasta
    output: # "results/{sample}/contamination_screening/work/fcs_gx/{hap}.{gx_name}.fcs_gx_report.txt"
        report_txt = "results/{sample}/contamination_screening/work/fcs_adaptor/outputs/{hap}/fcs_adaptor_report.txt",
    threads: 1
    singularity:
        f"{FCS_SOURCE_DIR}/images/fcs-adaptor.sif"
    resources:
        mem = 16,
        hrs = 12,
        load = 25,
    params:
        taxid=TAXID,
        fcs_adaptor_img = f"{FCS_SOURCE_DIR}/fcsadaptor/fcs-adaptor.sif",
        fcs_adaptor = f"{FCS_SOURCE_DIR}/fcsadaptor/run_fcsadaptor.sh",
    shell: """
        /app/fcs/bin/av_screen_x -o $( dirname {output.report_txt} ) --euk {input.filt_fasta}
        """


rule blast_mito:
    input:
        filt_fasta = rules.remove_short_contigs.output.filtered_fasta,
        mito_db = MITO_DB,
    output:
        blast_out = "results/{sample}/contamination_screening/work/blast/outputs/{hap}-mito.txt",
        flag = "results/{sample}/contamination_screening/flags/blast_mito.{hap}.done",
    threads: 8
    resources:
        mem = 4,
        hrs = 12,
    singularity: "docker://eichlerlab/ncbi-tk:0.1"
    shell: """
        set -euo pipefail    
        blastn -num_threads {threads} -query {input.filt_fasta} -db {input.mito_db} -outfmt 6 -evalue 1e-30 > {output.blast_out}
        touch {output.flag}
        """


rule filter_mito:
    input:
        blast_out = rules.blast_mito.output.blast_out,
        filt_fai = rules.remove_short_contigs.output.filtered_fai,
    output:
        mito_bed = "results/{sample}/contamination_screening/work/blast/beds/{hap}.mito.bed",
        flag = "results/{sample}/contamination_screening/flags/fliter_mito.{hap}.done",
    threads: 1
    resources:
        mem = 4,
        hrs = 12,
    singularity: "docker://eichlerlab/binf-basics:0.1"
    shell: """
        set -euo pipefail    
        bedtools coverage -a <( awk -vOFS="\\t" '{{print $1,"0",$2}}' {input.filt_fai} ) -b <( cut -f 1,7,8 {input.blast_out} ) | awk -vOFS="\\t" '{{ if ($NF >= 0.5) print $1,$2,$3}}' > {output.mito_bed}
        touch {output.flag}
        """


rule extract_mito:
    input:
        blast_out = rules.blast_mito.output.blast_out,
        filt_fasta = rules.remove_short_contigs.output.filtered_fasta,
    output:
        mito_fasta = "results/{sample}/contamination_screening/outputs/mito_fasta/{hap}-mt.fasta",
        flag = "results/{sample}/contamination_screening/flags/extract_mito.{hap}.done"
    params:
        mt_id = "NC_012920.1",
        mt_len = 16569,
        min_pid = 99.0,
        min_len = 16000
    threads: 1
    resources:
        mem = 4,
        hrs = 4,
    run:
        blast_header = ["qseqid","sseqid","pid","length","mismatch","gapopen","qstart_blast","qend_blast","sstart","send","evalue","bitscore"]
        blast_df = pd.read_csv(input.blast_out, sep="\t", header=None, names=blast_header)
        full_mt = (
            (blast_df["sseqid"] == params.mt_id) & 
            (
                ((blast_df["sstart"]==1) & (blast_df["send"]==params.mt_len))|
                ((blast_df["sstart"]==params.mt_len) & (blast_df["send"]==1))
            )
        )

        mito_df = blast_df[full_mt & (blast_df["length"] >= params.min_len) & (blast_df["pid"] >= params.min_pid)].copy()

        output_mito_fasta = output.mito_fasta
        if mito_df.empty:
            open(output_mito_fasta, "w").close()
            open(output.flag, "w").close()
            return

        mito_df.loc[:, "strand"] = mito_df.apply(lambda row: "+" if row["sstart"] < row["send"] else "-", axis=1) # bug fixed. q{start|end}_blast -> s{start|end}
        mito_df["qstart"] = mito_df[["qstart_blast","qend_blast"]].min(axis=1)
        mito_df["qend"] = mito_df[["qstart_blast","qend_blast"]].max(axis=1)
        mito_df["qstart_bed"] = mito_df["qstart"] - 1
        mito_df = mito_df.sort_values(
            by=["qseqid","bitscore","pid","length","mismatch","gapopen"],
            ascending=[True, False, False, False, True, True]
        )
        best = mito_df.drop_duplicates(subset=["qseqid"], keep="first").copy()

        fasta_dict = SeqIO.to_dict(SeqIO.parse(input.filt_fasta, "fasta"))

        out_records = []
        check_dup_seq = dict()
        num_id = 1
        num_rep_id = 0
        for idx, row in best.iterrows():
            contig = row["qseqid"]

            start = int(row["qstart_bed"])
            end = int(row["qend"])
            subrec = fasta_dict[contig][start:end]
            subrec_seq = str(subrec.seq).upper()
            if row["strand"] == "-":
                subrec = subrec.reverse_complement(id=True, name=True, description=True)
            new_id = f"{contig}-mt{num_id}"
            if HPRC_NAMING:
                name = full_manifest_df.loc[wildcards.sample, "NAME"]
                new_id = f"{name}#mt#{contig}-mt{num_id}"

            subrec.id = new_id
            subrec.name = new_id
            subrec.description = ""

            if not subrec_seq in check_dup_seq:
                num_rep_id += 1
                check_dup_seq[subrec_seq] = []
                check_dup_seq[subrec_seq].append(new_id)
                out_records.append(subrec)
                num_id += 1
            else:
                check_dup_seq[subrec_seq].append(new_id)

        with open(output_mito_fasta, "w") as fout:
            fasta_writer = FastaWriter(fout, wrap=None)
            fasta_writer.write_file(out_records)
        pysam.faidx(output_mito_fasta)
        open(output.flag, "w").close()


rule trim_bed:
    input:
        filt_fai = rules.remove_short_contigs.output.filtered_fai,
        gx_report = find_gx_report,
        flag = rules.run_fcs_gx.output.flag,
        adapt_report = rules.run_fcs_adapter.output.report_txt,
        mito_bed = rules.filter_mito.output.mito_bed,
    output:
        trim_file = "results/{sample}/contamination_screening/work/fastq_cleaning/gx_adapt_cleaning/beds/{hap}.trim.bed",
    threads: 1
    resources:
        mem = 16,
        hrs = 12,
    run:
        fai_df = pd.read_csv(
            input.filt_fai,
            sep="\t",
            header=None,
            names=["#seq_id", "length", "offset", "linebases", "linewidth"],
        )
        order = fai_df["#seq_id"].tolist()
        order_map = {seq: i for i, seq in enumerate(order)}

        gx_reports = input.gx_report
        out_df = pd.DataFrame()
        mito_df = pd.read_csv(
            input.mito_bed,
            sep="\t",
            header=None,
            names=["#seq_id", "start_pos", "end_pos"],
        )
        mito_df["reason"] = "mito_ebv_rdna"
        out_df = pd.concat([out_df, mito_df])
        df_gx = pd.concat(
            [
                pd.read_csv(
                    x,
                    sep="\t",
                    header=None,
                    comment="#",
                    names=[
                        "#seq_id",
                        "start_pos",
                        "end_pos",
                        "seq_len",
                        "action",
                        "div",
                        "agg_cont_cov",
                        "top_tax_name",
                    ],
                )
                for x in gx_reports
            ]
        )
        df_gx["start_pos"] = df_gx["start_pos"] - 1
        df_gx["reason"] = "foreign_contam"
        out_df = pd.concat(
            [out_df, df_gx[["#seq_id", "start_pos", "end_pos", "reason"]]]
        )
        df_adapt = pd.read_csv(input.adapt_report, sep="\t")
        df_adapt_trim = df_adapt.loc[df_adapt["action"] == "ACTION_TRIM"].copy()
        df_adapt_exc = df_adapt.loc[df_adapt["action"] == "ACTION_EXCLUDE"].copy()
        if len(df_adapt_trim) > 0:
            df_adapt_trim["#seq_id"] = df_adapt_trim["#accession"]
            df_adapt_trim["range"] = df_adapt_trim["range"].str.split(",")
            df_adapt_trim = df_adapt_trim.explode("range")
            df_adapt_trim["start_pos"] = (
                df_adapt_trim["range"]
                .str.split("..", expand=True, regex=False)[0]
                .astype(int)
                - 1
            )
            df_adapt_trim["end_pos"] = (
                df_adapt_trim["range"]
                .str.split("..", expand=True, regex=False)[1]
                .astype(int)
            )
            df_adapt_trim["reason"] = "adapter"
            out_df = pd.concat(
                [out_df, df_adapt_trim[["#seq_id", "start_pos", "end_pos", "reason"]]]
            )
        if len(df_adapt_exc) > 0:
            df_adapt_exc["#seq_id"] = df_adapt_exc["#accession"]
            df_adapt_exc["start_pos"] = 0
            df_adapt_exc["end_pos"] = df_adapt_exc["length"]
            df_adapt_exc["reason"] = "adapter"
            out_df = pd.concat(
                [out_df, df_adapt_exc[["#seq_id", "start_pos", "end_pos", "reason"]]]
            )

        out_df["sort_key"] = out_df["#seq_id"].map(order_map)
        out_df = (
            out_df.sort_values(["sort_key","start_pos"]).drop(columns="sort_key")
        )

        out_df[["#seq_id", "start_pos", "end_pos", "reason"]].to_csv(
            output.trim_file, sep="\t", header=False, index=False
        )


rule coerce_bed:
    input:
        trim_file = rules.trim_bed.output.trim_file,
        filt_fai = rules.remove_short_contigs.output.filtered_fai,
    output:
        regions_file = "results/{sample}/contamination_screening/work/fastq_cleaning/gx_adapt_cleaning/beds/{hap}.regions.out",
        flag = "results/{sample}/contamination_screening/flags/coerce_bed.{hap}.done"
    threads: 1
    resources:
        mem = 4,
        hrs = 12,
    singularity: "docker://eichlerlab/binf-basics:0.1"
    # removed sort to keep the order of contig names
    shell: """
        set -euo pipefail    
        bedtools complement -i {input.trim_file} -g {input.filt_fai}| awk '{{if ( $3 - $2 > 1000 ) print $1":"$2+1"-"$3}}' > {output.regions_file}
        touch {output.flag}
        """


rule trim_sequence:
    input:
        regions_file=rules.coerce_bed.output.regions_file,
        filt_fasta=rules.remove_short_contigs.output.filtered_fasta,
    output:
        cleaned_fasta = "results/{sample}/contamination_screening/work/fastq_cleaning/gx_adapt_cleaning/cleaned_fasta/{hap}.gx_adapt_cleaned.fasta",
        cleaned_index="results/{sample}/contamination_screening/work/fastq_cleaning/gx_adapt_cleaning/cleaned_fasta/{hap}.gx_adapt_cleaned.fasta.fai", ## contamination & adaptor cleaned
    threads: 1
    singularity: "docker://eichlerlab/binf-basics:0.1"
    resources:
        mem=4,
        hrs=12,
    ### removed sed 's/:/#/g'
    shell: """
        samtools faidx -r {input.regions_file} {input.filt_fasta} > {output.cleaned_fasta}
        samtools faidx {output.cleaned_fasta}
        """


rule blast_rdna:
    input:
        filt_fasta=rules.trim_sequence.output.cleaned_fasta,
        rdna_db=RDNA_DB,
    output:
        blast_out="results/{sample}/contamination_screening/work/blast/outputs/{hap}-rdna.txt",
        flag = "results/{sample}/contamination_screening/flags/blast_rdna.{hap}.done"
    threads: 8
    resources:
        mem=4,
        hrs=12,
    singularity: "docker://eichlerlab/ncbi-tk:0.1"
    shell: """
        set -euo pipefail
        blastn -num_threads {threads} -query {input.filt_fasta} -db {input.rdna_db} -outfmt 6 -evalue 1e-30 > {output.blast_out}
        touch {output.flag}
        """


rule filter_rdna:
    input:
        blast_out = rules.blast_rdna.output.blast_out,
        fai = rules.trim_sequence.output.cleaned_index,
    output:
        rdna_ctg = "results/{sample}/contamination_screening/work/fastq_cleaning/rdna_cleaning/beds/{hap}-rdna.bed",
        other_ctg = "results/{sample}/contamination_screening/work/fastq_cleaning/rdna_cleaning/beds/{hap}-rdna_clean.bed",
        flag = "results/{sample}/contamination_screening/flags/filter_rdna.{hap}.done"
    threads: 1
    resources:
        mem = 4,
        hrs = 12,
    singularity: "docker://eichlerlab/binf-basics:0.1"
    shell: """
        set -euo pipefail
        bedtools coverage -a <( awk -vOFS="\\t" '{{print $1,"0",$2}}' {input.fai} ) -b <( cut -f 1,7,8 {input.blast_out} ) | awk -vOFS="\\t" '{{ if ($NF >= 0.95) print $1}}' > {output.rdna_ctg}
        if [[ $( cat {output.rdna_ctg} | wc -l ) == 0 ]]; then
            cp {input.fai} {output.other_ctg}
        else
            cat {output.rdna_ctg} | tr "\\n" "|" | sed 's/|$//' | xargs -i grep -vE {{}} {input.fai} |  cut -f 1 > {output.other_ctg}
        fi 
        touch {output.flag}
        """


rule split_rdna:
    input:
        fasta = rules.trim_sequence.output.cleaned_fasta,
        fai = rules.trim_sequence.output.cleaned_index,
        rdna = rules.filter_rdna.output.rdna_ctg,
        others = rules.filter_rdna.output.other_ctg,
    output:
        cleaned_fasta = "results/{sample}/contamination_screening/work/rdna_cleaning/cleaned_fasta/{hap}.gx_adapt_rdna_cleaned.fasta", # adapt & r-DNA-cleaned
        cleaned_fai = "results/{sample}/contamination_screening/work/rdna_cleaning/cleaned_fasta/{hap}.gx_adapt_rdna_cleaned.fasta.fai",
        rdna_fasta = "results/{sample}/contamination_screening/outputs/rdna_fasta/{hap}-rdna.fasta",
        flag = "results/{sample}/contamination_screening/flags/split_rdna.{hap}.done"
    threads: 1
    resources:
        mem = 4,
        hrs = 12,
    singularity: "docker://eichlerlab/binf-basics:0.1"
    shell: """
        set -euo pipefail
        if [[ $( wc -l {input.rdna} | awk '{{print $1}}' ) == 0 ]]; then
            cp {input.fasta} {output.cleaned_fasta}
            touch {output.rdna_fasta}
        else
            samtools faidx -r {input.rdna} {input.fasta} > {output.rdna_fasta}
            samtools faidx -r {input.others} {input.fasta} > {output.cleaned_fasta}
        fi 
        samtools faidx {output.cleaned_fasta}
        touch {output.flag}
        """

rule rename_fasta:
    input:
        fasta = rules.remove_short_contigs.output.filtered_fasta,
        cleaned_fasta = rules.split_rdna.output.cleaned_fasta,
        mito_fasta = rules.extract_mito.output.mito_fasta,
    output:
        final_fasta = "results/{sample}/contamination_screening/outputs/final_fasta/{sample}_{hap}.fasta",
        final_fai = "results/{sample}/contamination_screening/outputs/final_fasta/{sample}_{hap}.fasta.fai"
    threads: 1
    resources:
        mem = 8,
        hrs = 12,
    run:
        final_fasta = output.final_fasta

        split_counter = defaultdict(int)
        original_fasta = pysam.FastaFile(input.fasta)
        cleaned_fasta_records = list(SeqIO.parse(input.cleaned_fasta, "fasta"))
        final_fasta_records = list()

        for record in cleaned_fasta_records:
            seq_name = str(record.id)
            #original_seq_name = seq_name.split(":")[0]
            original_seq_name = re.sub(r":[^:]+?$", "", seq_name)
            cleaned_sequence = str(record.seq)
            raw_sequence = original_fasta.fetch(original_seq_name)
            if cleaned_sequence == raw_sequence:
                cleaned_seq_name = original_seq_name
            else:
                split_counter[original_seq_name] += 1
                cleaned_seq_name = f"{original_seq_name}-split{split_counter[original_seq_name]}"
                # cleaned_seq_name = seq_name.replace(":","_trim_")
            if HPRC_NAMING:
                name = full_manifest_df.loc[wildcards.sample, "NAME"]
                hap_name = str(wildcards.hap).replace("hap","")
                cleaned_seq_name = f"{name}#{hap_name}#{cleaned_seq_name}"
            final_fasta_records.append(SeqRecord(Seq(cleaned_sequence), id=cleaned_seq_name, description=""))
        with open(final_fasta, "w") as fout:
            fasta_writer = FastaWriter(fout, wrap=None)
            fasta_writer.write_file(final_fasta_records)

        pysam.faidx(final_fasta) # indexing
