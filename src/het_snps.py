import argparse
import pandas as pd

from count_hets import HetSite


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--het-sites',
                        help='Heterozygous marker positions file',
                        required=True)
    parser.add_argument('--tumor-bam',
                        help='Tumor BAM or CRAM file',
                        required=True)
    parser.add_argument('--normal-bam',
                        help='Normal BAM or CRAM file',
                        required=True)
    parser.add_argument('--out',
                        help='Output TSV path',
                        required=True)
    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    args = get_args()
    hs = HetSite(args['het_sites'],
                 args['tumor_bam'],
                 args['normal_bam'])
    try:
        output_df = pd.DataFrame(hs.output)
        output_df = output_df[['contig', 'position', 'ref_allele', 'alt_allele',
                               't_ref_count', 't_alt_count', 'n_ref_count', 'n_alt_count']]
        output_df.to_csv(args['out'], sep='\t', index=False)
    finally:
        hs.close()


if __name__ == '__main__':
    main()
