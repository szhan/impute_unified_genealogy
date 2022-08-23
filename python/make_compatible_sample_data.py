import tskit
import numpy as np


def make_compatible_sample_data(sample_data, ancestors_ts):
    """
    Make an editable copy of a `sample_data` object, and edit it so that:
    (1) the derived alleles in `sample_data` not in `ancestors_ts` are marked as MISSING;
    (2) the allele list in `new_sample_data` corresponds to the allele list in `ancestors_ts`.

    N.B. Two `SampleData` attributes `sites_alleles` and `sites_genotypes`,
    which are not explained in the tsinfer API doc, are used to facilitate the editing.

    :param SampleData sample_data:
    :param TreeSequence ancestors_ts:
    :return SampleData:
    """
    new_sample_data = sample_data.copy()

    # Iterate through the sites in `ancestors_ts` using one generator,
    # while iterating through the sites in `sample_data` using another generator,
    # letting the latter generator catch up.
    sd_variants = sample_data.variants()
    sd_v = next(sd_variants)
    for ts_site in ancestors_ts.sites():
        while sd_v.site.position != ts_site.position:
            # Sites in `samples_data` but not in `ancestors_ts` are not imputed.
            # Also, leave them as is in the `sample_data`, but keep track of them.
            sd_v = next(sd_variants)

        sd_site_id = sd_v.site.id  # Site id in `sample_data`

        # CHECK that all the sites in `ancestors_ts` are biallelic.
        assert len(ts_site.alleles) == 2

        # Get the derived allele in `ancestors_ts` in nucleotide space
        ts_ancestral_allele = ts_site.ancestral_state
        ts_derived_allele = ts_site.alleles - {ts_ancestral_allele}
        assert len(ts_derived_allele) == 1  # CHECK
        ts_derived_allele = tuple(ts_derived_allele)[0]

        # CHECK that the ancestral allele should be the same
        # in both `ancestors_ts` and `sample_data`.
        assert ts_ancestral_allele == sd_v.alleles[0]

        if ts_derived_allele not in sd_v.alleles:
            # Case 1:
            # If the derived alleles in the `sample_data` are not in `ancestors_ts`,
            # then mark them as missing.
            #
            # The site in `sample_data` may be mono-, bi-, or multiallelic.
            #
            # We cannot determine whether the extra derived alleles in `sample_data`
            # are derived from 0 or 1 in `ancestors_ts` anyway.
            new_sample_data.sites_genotypes[sd_site_id] = np.where(
                sd_v.genotypes != 0,  # Keep if ancestral
                tskit.MISSING_DATA,  # Otherwise, flag as missing
                0,
            )
            print(
                f"Site {sd_site_id} has no matching derived alleles in the query samples."
            )
            # Update allele list
            new_sample_data.sites_alleles[sd_site_id] = [ts_ancestral_allele]
        else:
            # The allele lists in `ancestors_ts` and `sample_data` may be different.
            ts_derived_allele_index = sd_v.alleles.index(ts_derived_allele)

            if ts_derived_allele_index == 1:
                # Case 2:
                # Both the ancestral and derived alleles correspond exactly.
                if len(sd_v.alleles) == 2:
                    continue
                # Case 3:
                # The derived allele in `ancestors_ts` is indexed as 1 in `sample_data`,
                # so mark alleles >= 2 as missing.
                new_sample_data.sites_genotypes[sd_site_id] = np.where(
                    sd_v.genotypes > 1,  # 0 and 1 should be kept "as is"
                    tskit.MISSING_DATA,  # Otherwise, flag as missing
                    sd_v.genotypes,
                )
                print(
                    f"Site {sd_site_id} has extra derived allele(s) in the query samples (set as missing)."
                )
            else:
                # Case 4:
                #   The derived allele in `ancestors_ts` is NOT indexed as 1 in `sample_data`,
                #   so the alleles in `sample_data` needs to be reordered,
                #   such that the 1-indexed allele is also indexed as 1 in `ancestors_ts`.
                new_sample_data.sites_genotypes[sd_site_id] = np.where(
                    sd_v.genotypes == 0,
                    0,  # Leave ancestral allele "as is"
                    np.where(
                        sd_v.genotypes == ts_derived_allele_index,
                        1,  # Change it to 1 so that it corresponds to `ancestors_ts`
                        tskit.MISSING_DATA,  # Otherwise, mark as missing
                    ),
                )
                print(
                    f"Site {sd_site_id} has the target derived allele at a different index."
                )
            # Update allele list
            new_sample_data.sites_alleles[sd_site_id] = [
                ts_ancestral_allele,
                ts_derived_allele,
            ]

    new_sample_data.finalise()

    return new_sample_data
