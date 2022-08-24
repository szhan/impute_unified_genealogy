import tskit
import tsinfer
import numpy as np
import datetime

import sys
sys.path.append("python/")
import measures


def function(ts_ref, sd_query, ts_imputed):
    """
    Evaluate imputation performance

    :param ts_ref:
    :param sd_query:
    :param ts_imputed:
    :type ts_ref: tskit.TreeSequence
    :type sd_query: tsinfer.SampleData
    :type ts_imputed: tskit.TreeSequence
    :return: None
    :rtype: None
    """
    results = None
    for v_ref, v_query_pre, v_query_post in zip(
        ts_ref.variants(),  # Reference genomes from which to get the minor allele and MAF
        sd_query.variants(),  # Query genomes before imputation
        ts_imputed.variants(),  # Query genomes with masked sites imputed
    ):
        if v_query_imputed.site.position in masked_site_positions:
            # CHECK that ancestral states are identical.
            assert (
                v_ref.alleles[0] == sd_query.sites_alleles[v_query_true.site.id][0]
            )
            assert (
                v_ref.alleles[0]
                == sd_query_masked.sites_alleles[v_query_masked.site.id][0]
            )
            assert v_ref.alleles[0] == v_query_imputed.alleles[0]

            # TODO:
            #   Why doesn't `v.num_alleles` always reflect the number of genotypes
            #   after simplifying?
            if len(set(v_ref.genotypes)) == 1:
                # Monoallelic sites in `ts_ref` are not imputed
                # TODO: Revisit
                continue

            assert v_ref.num_alleles == 2
            assert set(v_query_masked.genotypes) == set([-1])
            assert not np.any(v_query_imputed.genotypes == -1)

            # Note: A minor allele in `ts_ref` may be a major allele in `sd_query`
            freqs_ref = v_ref.frequencies()
            af_0 = freqs_ref[v_ref.alleles[0]]
            af_1 = freqs_ref[v_ref.alleles[1]]

            # Get MAF from `ts_ref`
            # Definition of a minor allele: < 0.50
            if af_1 < af_0:
                minor_allele_index = 1
                maf = af_1
            else:
                minor_allele_index = 0
                maf = af_0

            # Assess imputation performance
            total_concordance = measures.compute_concordance(
                genotypes_true=v_query_true.genotypes,
                genotypes_imputed=v_query_imputed.genotypes,
            )
            iqs = measures.compute_iqs(
                genotypes_true=v_query_true.genotypes,
                genotypes_imputed=v_query_imputed.genotypes,
            )

            # line.shape = (1, 4)
            line = np.array(
                [
                    [v_ref.site.position, maf, total_concordance, iqs],
                ]
            )
            if results is None:
                results = line
            else:
                results = np.append(results, line, axis=0)

    end_dt = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    ### Write results
    header_text = (
        "\n".join(
            [
                "#" + "start_timestamp" + "=" + f"{start_dt}",
                "#" + "end_timestamp" + "=" + f"{end_dt}",
                "#" + "tskit" + "=" + f"{tskit.__version__}",
                "#" + "tsinfer" + "=" + f"{tsinfer.__version__}",
                "#" + "size_ref" + "=" + f"{ts_ref.num_samples}",
                "#" + "size_query" + "=" + f"{sd_query.num_samples}",
            ]
        )
        + "\n"
    )

    header_text += ",".join(["position", "maf", "total_concordance", "iqs"])

    np.savetxt(
        out_csv,
        results,
        fmt="%.10f",
        delimiter=",",
        newline="\n",
        comments="",
        header=header_text,
    )