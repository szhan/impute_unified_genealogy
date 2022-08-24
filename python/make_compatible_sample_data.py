import tskit
import tsinfer
import numpy as np
from tqdm import tqdm


def make_compatible_sample_data(sample_data, ancestors_ts):
    """
    Make an editable copy of a `sample_data` object, and edit it so that:
    (1) the derived alleles in `sample_data` not in `ancestors_ts` are marked as MISSING;
    (2) the allele list in `new_sample_data` corresponds to the allele list in `ancestors_ts`.
    (3) sites in `ancestors_ts` but not in `sample_data` are added to `new_sample_data` with all the genotypes MISSING.

    All the sites in `sample_data` and `ancestors_ts` must be biallelic.

    N.B. Two `SampleData` attributes `sites_alleles` and `sites_genotypes`,
    which are not explained in the tsinfer API doc, are used to facilitate the editing.

    :param SampleData sample_data:
    :param TreeSequence ancestors_ts:
    :return: Edited copy of sample_data.
    :rtype SampleData:
    """
    ts_site_pos = ancestors_ts.sites_position
    sd_site_pos = sample_data.sites_position[:]
    all_site_pos = sorted(set(ts_site_pos).union(set(sd_site_pos)))

    with tsinfer.SampleData(sequence_length=ancestors_ts.sequence_length) as new_sd:
        # Add sites
        for pos in tqdm(all_site_pos):
            if pos in ts_site_pos and pos not in sd_site_pos:
                # Case 1:
                # Site found in `ancestors_ts` but not `sample_data`
                # Add the site to `new_sample_data` with all genotypes MISSING.
                ts_site = ancestors_ts.site(position=pos)
                assert len(ts_site.alleles) == 2,\
                    f"Non-biallelic site {ts_site.alleles}"
                ts_ancestral_state = ts_site.ancestral_state
                ts_derived_state = list(ts_site.alleles - {ts_ancestral_state})[0]
                new_sd.add_site(
                    position=pos,
                    genotypes=np.full(sample_data.num_samples, tskit.MISSING_DATA),
                    alleles=[ts_ancestral_state, ts_derived_state]
                )
            elif pos in ts_site_pos and pos in sd_site_pos:
                # Case 2:
                # Site found in both `ancestors_ts` and `sample_data`
                # Align the allele lists and genotypes if unaligned.
                # Add the site to `new_sample_data` with (aligned) genotypes from `sample_data`.
                ts_site = ancestors_ts.site(position=pos)
                sd_site_id = sd_site_pos.tolist().index(pos)
                sd_site_alleles = sample_data.sites_alleles[sd_site_id]
                assert len(ts_site.alleles) == 2,\
                    f"Non-biallelic site {ts_site.alleles}"
                assert len(set(sd_site_alleles) - {None}) == 2,\
                    f"Non-biallelic site {sd_site_alleles}"
                ts_ancestral_state = ts_site.ancestral_state
                ts_derived_state = list(ts_site.alleles - {ts_ancestral_state})[0]
                if list(ts_site.alleles) == sd_site_alleles:
                    # Case 2a:
                    # Both alleles are in `ancestors_ts` and `sample_data`.
                    # Already aligned, so no need to realign.
                    new_sd.add_site(
                        position=pos,
                        genotypes=sample_data.sites_genotypes[sd_site_id],
                        alleles=[ts_ancestral_state, ts_derived_state]
                    )
                elif set(ts_site.alleles) == set(sd_site_alleles):
                    # Case 2b:
                    # Both alleles are in `ancestors_ts` and `sample_data`.
                    # Align them by flipping the alleles in `sample_data`.
                    sd_site_gt = sample_data.sites_genotypes[sd_site_id]
                    new_gt = np.where(
                        sd_site_gt == tskit.MISSING_DATA,
                        tskit.MISSING_DATA,
                        np.where(sd_site_gt == 0, 1, 0) # Flip
                    )
                    new_sd.add_site(
                        position=pos,
                        genotypes=new_gt,
                        alleles=[ts_ancestral_state, ts_derived_state]
                    )
                else:
                    # Case 2c:
                    # The allele(s) present in `sample_data` but absent in `ancestor_ts`
                    # is always incorrectly imputed.
                    # It is best to ignore these sites when assess imputation performance.
                    new_sd.add_site(
                        position=pos,
                        genotypes=np.full(sample_data.num_samples, tskit.MISSING_DATA),
                        alleles=[ts_ancestral_state, ts_derived_state]
                    )
            elif pos not in ts_site_pos and pos in sd_site_pos:
                # Case 3:
                # Site found in `sample_data` but not `ancestors_ts`
                # Add the site to `new_sample_data` with the original genotypes from `sample_data`.
                sd_site_id = sd_site_pos.tolist().index(pos)
                sd_site_alleles = sample_data.sites_alleles[sd_site_id]
                assert len(set(sd_site_alleles) - {None}) == 2,\
                    f"Non-biallelic site {sd_site_alleles}"
                new_sd.add_site(
                    position=pos,
                    genotypes=sample_data.sites_genotypes[sd_site_id],
                    alleles=sample_data.sites_alleles[sd_site_id]
                )
            else:
                raise ValueError(f"Position {pos} must be in the tree sequence and/or sample data.")

        # Add individuals
        for ind in sample_data.individuals():
            new_sd.add_individual(
                ploidy=len(ind.samples),
                population=ind.population,
                metadata=ind.metadata
            )

    return new_sd
