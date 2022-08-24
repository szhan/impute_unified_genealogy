import click
from datetime import datetime
import gzip
import sys

import numpy as np

import tskit
import tsinfer

sys.path.append("python/")
import read_vcf
import make_compatible_sample_data as make_compatible
import util


@click.command()
@click.option(
    "--in_trees", "-i1",
    type=click.Path(exists=True),
    required=True,
    help="Input tree sequence consisted of reference genomes (binary format)"
)
@click.option(
    "--in_bcf", "-i2",
    type=click.Path(exists=True),
    required=True,
    help="Input BCF containing missing sites to impute (binary format)"
)
@click.option(
    "--out_dir", "-o",
    type=click.Path(exists=True, dir_okay=True),
    required=True,
    help="Output directory"
)
@click.option(
    "--out_prefix", "-p",
    type=str,
    required=True,
    help="Prefix of output files"
)
@click.option(
    "--remove_leaves",
    is_flag=True,
    help="Remove leaves when making an ancestors tree sequence from the input tree sequence"
)
@click.option(
    "--verbose",
    is_flag=True,
    help="Print out site information after each processing step.",
)
def run_pipeline(in_trees, in_bcf, out_dir, out_prefix, remove_leaves, verbose):
    start_dt = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    out_samples_file = out_dir + "/" + out_prefix + ".samples"
    out_bcf_file = out_dir + "/" + out_prefix + ".bcf"
    out_csc_file = out_dir + "/" + out_prefix + ".csv"

    ### Create an ancestor ts from the reference genomes
    # Multiallelic sites are automatically removed when generating an ancestor ts.
    # Sites which are biallelic in the full sample set but monoallelic in the ref. sample set are removed.
    # So, only biallelic sites are retained in the ancestor ts.
    print(f"INFO: Loading tree sequence from {in_trees}")
    ts_ref = tskit.load(file=in_trees)
    ts_ref_seq_len = ts_ref.sequence_length

    print(f"INFO: Generating ancestors tree sequence")
    ts_anc = tsinfer.eval_util.make_ancestors_ts(ts=ts_ref, remove_leaves=remove_leaves)

    if verbose:
        print(
            f"TS anc has {ts_anc.num_samples} sample genomes ({ts_anc.sequence_length} bp)"
        )
        print(f"TS anc has {ts_anc.num_sites} sites and {ts_anc.num_trees} trees")
        print("TS anc")
        util.count_sites_by_type(ts_anc)

    ### Create a SampleData object holding the query genomes
    print(f"INFO: Loading VCF from {in_bcf}")
    sd_query_raw = read_vcf.create_sample_data_from_vcf(
        in_file=in_bcf, out_file=out_samples_file, max_right_position=ts_ref_seq_len - 1)

    if verbose:
        print(
            f"SD query has {sd_query_raw.num_samples} sample genomes ({sd_query_raw.sequence_length} bp)"
        )
        print(f"SD query has {sd_query_raw.num_sites} sites")
        print("SD query")
        util.count_sites_by_type(sd_query_raw)

    print("INFO: Processing VCF")
    sd_query_pre = make_compatible.make_compatible_sample_data(sample_data=sd_query_raw, ancestors_ts=ts_anc)

    ### Impute the query genomes
    print(f"INFO: Imputing query genomes")
    ts_imputed = tsinfer.match_samples(sample_data=sd_query_pre, ancestors_ts=ts_anc)

if __name__ == "__main__":
    run_pipeline()
