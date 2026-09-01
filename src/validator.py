

import sys
import pandas as pd
import gzip


def vcf_line_to_hgvs_g(vcf_line: str) -> str:
    """
    Parses a single VCF line and outputs the HGVS genomic (g.) syntax.
    Expects standard VCF columns: CHROM, POS, ID, REF, ALT...
    """
    # Skip header lines
    if vcf_line.startswith('#'):
        return "Header line skipped."
        
    # Split VCF columns by tab
    columns = vcf_line.strip().split('\t')
    if len(columns) < 5:
        return "Error: Invalid VCF line format (fewer than 5 columns)."
        
    chrom, pos, _, ref, alt = columns[:5]
    
    # Basic data validation for a standard SNP
    if len(ref) != 1 or len(alt) != 1 or ref == alt:
        return f"Error: Line at position {pos} is not a valid single nucleotide substitution (SNP)."
    
    # Construct HGVS g. syntax: {accession}:g.{position}{ref}>{alt}
    hgvs_g = f"{chrom}:{pos}_{ref}>{alt}"
    return hgvs_g

# Example Usage:
if __name__ == "__main__":
    table=sys.argv[1]
    vcf=sys.argv[2]
    hgvs_gs = []
    with gzip.open(vcf, 'rt') as vcf_file:
        for line in vcf_file:
            hgvs_g = vcf_line_to_hgvs_g(line)
            if not line.startswith('#'):
                if ">" in hgvs_g:
                    hgvs_gs.append(hgvs_g)
    # features
    features_df = pd.read_csv(table)
    t_hgvs_gs = features_df['Variant'].tolist()
    if set(t_hgvs_gs).difference(set(hgvs_gs)):
        print("Variants in feature table but not in VCF:")
        print(set(t_hgvs_gs).difference(set(hgvs_gs)))
    if set(hgvs_gs).difference(set(t_hgvs_gs)):
        print("Variants in VCF but not in feature table:")
        print(set(hgvs_gs).difference(set(t_hgvs_gs)))
