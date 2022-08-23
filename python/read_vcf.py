import click
import numpy as np
import tskit
import tsinfer
import cyvcf2


def add_haploid_individuals(vcf, sample_data):
    """
    TODO: Parse anno file.
    """
    for name in vcf.samples:
        sample_data.add_individual(
            ploidy = 1,
            metadata = {"name": name}
        )


def add_haploid_sites(vcf, sample_data):
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

        assert np.all(v.ploidy == 1)
        new_genotypes = [allele_index[g[0]] for g in v.genotypes]

        sample_data.add_site(pos, genotypes=new_genotypes, alleles=new_alleles)


@click.command()
@click.option('--in_file', '-i', required=True, help="Input VCF/BCF file")
@click.option('--out_file', '-o', required=True, help="Output samples file")
def create_sample_data_from_vcf(in_file, out_file):
    """
    Convert a VCF/BCF file into a SampleData object,
    and write it to a samples file for input to tsinfer.
    """
    vcf = cyvcf2.VCF(in_file, strict_gt=True)

    with tsinfer.SampleData(
        path=out_file
    ) as sample_data:
        add_haploid_individuals(vcf, sample_data)
        add_haploid_sites(vcf, sample_data)

    return(sample_data)


if __name__ == "__main__":
    create_sample_data_from_vcf()