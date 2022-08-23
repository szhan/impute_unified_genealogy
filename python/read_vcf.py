import click
import numpy as np
import tskit
import tsinfer
import cyvcf2


def add_haploid_individuals(vcf, sample_data):
    """
    Add individuals named in cyvcf2.VCF object to a SampleData object,
    including indiividual names in the metadata field.

    This only processes variable sites from haploid individuals.

    TODO: Parse anno file.

    :param cyvcf2.VCF vcf:
    :param tsinfer.SampleData sample_data:
    """
    for name in vcf.samples:
        sample_data.add_individual(
            ploidy = 1,
            metadata = {"name": name}
        )


def add_haploid_sites(vcf, sample_data, max_right_position, verbose):
    """
    Read the sites in the VCF and add them to the SampleData object,
    reordering the alleles to put the ancestral allele first, if it is available.

    This only processes variable sites from haploid individuals.

    :param cyvcf2.VCF vcf:
    :param tsinfer.SampleData sample_data:
    :param float max_right_position: Rightmost 0-based coordinate position
    :param bool verbose: If True, print information about duplicate sites. (default=False)
    """
    pos = 0
    for v in vcf:
        # Check for duplicate positions. Keep the first encountered.
        if pos == v.POS:
            if verbose:
                print(f"WARN: Duplicate positions for variant at position {pos}.")
            continue
        else:
            pos = v.POS

        if max_right_position != None:
            if pos > max_right_position:
                continue

        ancestral = v.INFO.get("AA", v.REF)

        old_alleles = [v.REF] + v.ALT
        new_alleles = [ancestral] + list(set(old_alleles) - {ancestral})

        allele_index = {
            old_index: new_alleles.index(allele)
            for old_index, allele in enumerate(old_alleles)
        }

        if v.num_unknown > 0:
            allele_index[-1] = tskit.MISSING_DATA
            new_alleles += [None]

        assert np.all(v.ploidy == 1)
        new_genotypes = [allele_index[g[0]] for g in v.genotypes]

        sample_data.add_site(pos, genotypes=new_genotypes, alleles=new_alleles)


def create_sample_data_from_vcf(in_file, out_file, max_right_position=None, verbose=False):
    """
    Convert a VCF file into a SampleData object,
    and write it to a samples file for input to tsinfer.
    """
    vcf = cyvcf2.VCF(in_file, strict_gt=True)

    with tsinfer.SampleData(
        path=out_file
    ) as sample_data:
        add_haploid_individuals(vcf, sample_data)
        add_haploid_sites(vcf, sample_data, max_right_position, verbose)

    return(sample_data)