import click
from datetime import datetime
import gzip
import json
import sys

import numpy as np

import msprime
import tskit
import tsinfer
from tsinfer import make_ancestors_ts

sys.path.append("python/")
import read_vcf
import make_compatible_sample_data as make_compatible
#import measures


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

    print(f"INFO: Generating ancestors tree sequence")
    ts_anc = tsinfer.eval_util.make_ancestors_ts(ts=ts_ref, remove_leaves=remove_leaves)

    if verbose:
        print(
            f"TS anc has {ts_anc.num_samples} sample genomes ({ts_anc.sequence_length} bp)"
        )
        print(f"TS anc has {ts_anc.num_sites} sites and {ts_anc.num_trees} trees")
        print("TS anc")
        #util.count_sites_by_type(ts_anc)

    ### Create a SampleData object holding the query genomes
    print(f"INFO: Loading VCF from {in_bcf}")
    sd_query_raw = read_vcf.create_sample_data_from_vcf(in_file=in_bcf, out_file=out_samples_file)

    if verbose:
        print(
            f"SD query has {sd_query_raw.num_samples} sample genomes ({sd_query_raw.sequence_length} bp)"
        )
        print(f"SD query has {sd_query_raw.num_sites} sites")
        print("SD query")
        #util.count_sites_by_type(sd_query_raw)

    #assert util.check_site_positions_ts_issubset_sd(ts_anc, sd_query)

    print("INFO: Processing VCF")
    sd_query_pre = make_compatible.make_compatible_sample_data(sample_data=sd_query_raw, ancestors_ts=ts_anc)

    ### Impute the query genomes
    print(f"INFO: Imputing query genomes")
    ts_imputed = tsinfer.match_samples(sample_data=sd_query_pre, ancestors_ts=ts_anc)

    # ### Evaluate imputation performance
    # ts_ref_site_positions = [s.position for s in ts_ref.sites()]
    # sd_query_site_positions = [s.position for s in sd_query_pre.sites()]
    # ts_imputed_site_positions = [s.position for s in ts_imputed.sites()]

    # assert len(ts_ref_site_positions) == len(sd_query_site_positions)
    # assert len(ts_ref_site_positions) == len(ts_imputed_site_positions)

    # assert set(ts_ref_site_positions) == set(sd_query_site_positions)
    # assert set(ts_ref_site_positions) == set(ts_imputed_site_positions)

    # results = None
    # for v_ref, v_query_pre, v_query_post in zip(
    #     ts_ref.variants(),  # Reference genomes from which to get the minor allele and MAF
    #     sd_query_pre.variants(),  # Query genomes before imputation
    #     ts_imputed.variants(),  # Query genomes with masked sites imputed
    # ):
    #     if v_query_imputed.site.position in masked_site_positions:
    #         # CHECK that ancestral states are identical.
    #         assert (
    #             v_ref.alleles[0] == sd_query.sites_alleles[v_query_true.site.id][0]
    #         )
    #         assert (
    #             v_ref.alleles[0]
    #             == sd_query_masked.sites_alleles[v_query_masked.site.id][0]
    #         )
    #         assert v_ref.alleles[0] == v_query_imputed.alleles[0]

    #         # TODO:
    #         #   Why doesn't `v.num_alleles` always reflect the number of genotypes
    #         #   after simplifying?
    #         if len(set(v_ref.genotypes)) == 1:
    #             # Monoallelic sites in `ts_ref` are not imputed
    #             # TODO: Revisit
    #             continue

    #         assert v_ref.num_alleles == 2
    #         assert set(v_query_masked.genotypes) == set([-1])
    #         assert not np.any(v_query_imputed.genotypes == -1)

    #         # Note: A minor allele in `ts_ref` may be a major allele in `sd_query`
    #         freqs_ref = v_ref.frequencies()
    #         af_0 = freqs_ref[v_ref.alleles[0]]
    #         af_1 = freqs_ref[v_ref.alleles[1]]

    #         # Get MAF from `ts_ref`
    #         # Definition of a minor allele: < 0.50
    #         if af_1 < af_0:
    #             minor_allele_index = 1
    #             maf = af_1
    #         else:
    #             minor_allele_index = 0
    #             maf = af_0

    #         # Assess imputation performance
    #         total_concordance = measures.compute_concordance(
    #             genotypes_true=v_query_true.genotypes,
    #             genotypes_imputed=v_query_imputed.genotypes,
    #         )
    #         iqs = measures.compute_iqs(
    #             genotypes_true=v_query_true.genotypes,
    #             genotypes_imputed=v_query_imputed.genotypes,
    #         )

    #         # line.shape = (1, 4)
    #         line = np.array(
    #             [
    #                 [v_ref.site.position, maf, total_concordance, iqs],
    #             ]
    #         )
    #         if results is None:
    #             results = line
    #         else:
    #             results = np.append(results, line, axis=0)

    # end_dt = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    # ### Write results
    # header_text = (
    #     "\n".join(
    #         [
    #             "#" + "start_timestamp" + "=" + f"{start_dt}",
    #             "#" + "end_timestamp" + "=" + f"{end_dt}",
    #             "#" + "tskit" + "=" + f"{tskit.__version__}",
    #             "#" + "tsinfer" + "=" + f"{tsinfer.__version__}",
    #             "#" + "size_ref" + "=" + f"{ts_ref.num_samples}",
    #             "#" + "size_query" + "=" + f"{sd_query.num_samples}",
    #         ]
    #     )
    #     + "\n"
    # )

    # header_text += ",".join(["position", "maf", "total_concordance", "iqs"])

    # np.savetxt(
    #     out_csv,
    #     results,
    #     fmt="%.10f",
    #     delimiter=",",
    #     newline="\n",
    #     comments="",
    #     header=header_text,
    # )


if __name__ == "__main__":
    run_pipeline()
