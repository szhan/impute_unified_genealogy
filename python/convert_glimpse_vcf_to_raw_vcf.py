import datetime
import click
import numpy as np
import cyvcf2


@click.command()
@click.option('--in_file', '-i', required=True, help="Input BCF file")
@click.option('--out_file', '-o', required=True, help="Output compressed BCF file")
@click.option('--verbose', default=False, help="Print extra information about the variants")
def convert_glimpse_vcf_to_raw_vcf(in_file, out_file, verbose):
    """
    Obtain VCF/BCF with missing genotypes from a GLIMPSE imputed VCF/BCF file.
    Since sequencing coverage in ancient WGS data is often too low for phasing,
    the genotypes in the rebuilt VCF/BCF are all represented as haploid.

    The following numbers are denoted as:
    0 = REF
    1 = ALT
    -1 = UNKNOWN

    Per-sample genotypes are assigned UNKNOWN if at least one condition is met:
    1) Non-SNP;
    2) Total depth is 0;
    3) Total depth is -1; and
    4) Ref. depth is equal to alt. depth (ambiguous read support).

    :param str in_file: Input BCF file
    :param str out_file: Output BCF file
    :param bool verbose: If True, print extra info (default = False)
    """
    REF_GENOTYPE = (0, False,)
    ALT_GENOTYPE = (1, False,)
    UNK_GENOTYPE = (-1, False,)

    num_sites = 0  # Total number of variable sites
    num_snps = 0
    num_indels = 0
    num_svs = 0

    # The strict_gt setting affects v.gt_types but not v.genotypes.
    # default:
    #   0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
    # strict_gt=True:
    #   0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
    vcf = cyvcf2.VCF(in_file, strict_gt=True)
    DT = datetime.datetime.now()
    RAW_HEADER = " ".join([
        "##convert_glimpse_vcf_to_raw_vcf" + ";",
        f"Date={DT.day}/{DT.month}/{DT.year} - {DT.hour}:{DT.minute}:{DT.second}"
    ])
    vcf.add_to_header(RAW_HEADER)

    w = cyvcf2.Writer(out_file, tmpl=vcf, mode="wb")
    w.write_header()

    pos = 0
    for v in vcf:
        # Check for duplicate positions. Keep the first encountered.
        if pos == v.POS:
            if verbose:
                print(f"WARN: Duplicate positions for variant at position {pos}.")
            continue
        else:
            pos = v.POS

        num_sites += 1
        num_snps = num_snps + 1 if v.is_snp else num_snps
        num_indels = num_indels + 1 if v.is_indel else num_indels
        num_svs = num_svs + 1 if v.is_sv else num_svs

        REF_DP = v.gt_ref_depths
        ALT_DP = v.gt_alt_depths
        TOT_DP = v.gt_depths

        num_samples = v.num_called + v.num_unknown

        # Clear INFO fields and select FORMAT fields
        # Site-level data
        del v.INFO["RAF"]
        del v.INFO["AF"]
        del v.INFO["INFO"]
        # Sample-level data
        v.set_format("GP", np.full(num_samples, '.', dtype=np.bytes_))
        v.set_format("HS", np.full(num_samples, '.', dtype=np.bytes_))
        v.set_format("PL", np.full(num_samples, '.', dtype=np.bytes_))

        # If it is not a SNP, then assign UNKNOWN.
        # TODO: Is this appropriate?
        if not v.is_snp:
            if verbose:
                print(" ".join(str(x) for x in [
                      "Non-SNP", ":", v.CHROM, v.POS, v.num_unknown, v.num_het]))
            v.genotypes = [UNK_GENOTYPE] * num_samples
            v.genotypes = v.genotypes
            w.write_record(v)
            continue

        # If the depths of all the samples in a site are 0 or -1, then assign UNKNOWN.
        # TODO: Ask Ali why total depth can be either 0 or -1,
        #       i.e. both ref. depth and alt. depth are either 0 or -1.
        if np.all(TOT_DP == 0) or np.all(TOT_DP) == -1:
            v.genotypes = [UNK_GENOTYPE] * num_samples
            v.genotypes = v.genotypes
            w.write_record(v)
            continue

        # If DP favors REF, then assign REF.
        # If DP favors ALT, then assign ALT.
        # Otherwise, assign UNKNOWN.
        is_refdp_gt_altdp = REF_DP > ALT_DP
        is_refdp_lt_altdp = REF_DP < ALT_DP
        is_refdp_et_altdp = REF_DP == ALT_DP

        is_totdp_zero = TOT_DP == 0
        is_totdp_negative_one = TOT_DP == -1
        is_unknown = is_refdp_et_altdp | is_totdp_zero | is_totdp_negative_one

        for i in np.arange(num_samples):
            if is_refdp_gt_altdp[i]:
                tmp_gt = REF_GENOTYPE
            elif is_refdp_lt_altdp[i]:
                tmp_gt = ALT_GENOTYPE
            else:
                tmp_gt = UNK_GENOTYPE
            v.genotypes[i] = tmp_gt

        v.genotypes = v.genotypes
        w.write_record(v)

    w.close()
    vcf.close()


if __name__ == "__main__":
    convert_glimpse_vcf_to_raw_vcf()
