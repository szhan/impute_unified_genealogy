import numpy as np
import tskit
import tsinfer
import cyvcf2


def get_chromosome_length(vcf):
    assert len(vcf.seqlens) == 1
    return vcf.seqlens[0]


def add_individuals(vcf, samples):
    for name in vcf.samples:
        samples.add_individual(
            ploidy = 1,
            metadata = {"name": name}
        )


def add_haploid_sites(vcf, samples):
    """
    Read the sites in the VCF and add them to the SampleData object,
    reordering the alleles to put the ancestral allele first, if it is available.
    """
    pos = 0
    for v in vcf:
        # Check for duplicate positions. Keep the first encountered.
        if pos == v.POS:
            print(f"WARN: Duplicate positions for variant at position {pos}.")
            continue
        else:
            pos = v.POS

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

        genotypes = [
            allele_index[old_index]
            for row in v.genotypes
            for old_index in row[0:ploidy_level] # ???
        ]

        samples.add_haploid_site(pos, genotypes=genotypes, alleles=new_alleles)


def create_sample_data_from_vcf(in_file):
    vcf = cyvcf2.VCF(in_file, strict_gt=True)

    with tsinfer.SampleData(
        sequence_length = get_chromosome_length(vcf)
    ) as samples:
        add_individuals(vcf, samples)
        add_haploid_sites(vcf, samples)

    return(samples)