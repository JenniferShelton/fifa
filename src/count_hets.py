
import pandas as pd
import pysam
import sys
import argparse


class HetSite():
    def __init__(self, het_sites_file, tumor_bam, normal_bam):
        self.output = {'contig': [],
                       'position': [],
                       'ref_allele': [],
                       'alt_allele': [],
                       't_ref_count': [],
                       't_alt_count': [],
                       'n_ref_count': [],
                       'n_alt_count': []
                       }
        self.tumor_bam = self.open_bam(tumor_bam)
        self.normal_bam = self.open_bam(normal_bam)
        self.het_sites = self.load_sites(het_sites_file)
        for row in self.het_sites.itertuples(index=False):
            self.load_het_site_counts(row.chrom,
                                      int(row.position),
                                      row.ref,
                                      row.alt)


    def open_bam(self, bam):
        samfile = pysam.AlignmentFile(bam, "rb")
        return samfile

    def fetch_pileups(self, chrom, pos, samfile):
        pileup = samfile.pileup(chrom, int(pos)-1, int(pos),
                                max_depth=10000,
                                truncate=True)
        return pileup

    def load_sites(self, het_sites_file):
        '''Load marker sites and normalize to columns: chrom, position, ref, alt.'''
        def normalize_columns(columns):
            return [str(c).strip().lower() for c in columns]

        df = pd.read_csv(het_sites_file, sep='\t', comment='@', dtype=str)
        df.columns = normalize_columns(df.columns)

        if {'contig', 'position', 'ref_allele', 'alt_allele'}.issubset(df.columns):
            sites = df.rename(columns={'contig': 'chrom',
                                       'ref_allele': 'ref',
                                       'alt_allele': 'alt'})
        elif {'chrom', 'position', 'ref', 'alt'}.issubset(df.columns):
            sites = df
        elif {'chrom', 'start', 'ref_alt'}.issubset(df.columns):
            sites = df
            sites['position'] = sites['start']
            ref_alt_split = sites['ref_alt'].astype(str).str.split('/', n=1, expand=True)
            sites['ref'] = ref_alt_split[0]
            sites['alt'] = ref_alt_split[1]
        else:
            # Fallback for headerless marker tables.
            raw = pd.read_csv(het_sites_file, sep='\t', header=None, comment='@', dtype=str)
            if raw.shape[1] >= 5:
                raw = raw.iloc[:, :5]
                raw.columns = ['chrom', 'start', 'end', 'strand', 'ref_alt']
                raw['position'] = raw['start']
                ref_alt_split = raw['ref_alt'].astype(str).str.split('/', n=1, expand=True)
                raw['ref'] = ref_alt_split[0]
                raw['alt'] = ref_alt_split[1]
                sites = raw
            elif raw.shape[1] >= 4:
                raw = raw.iloc[:, :4]
                raw.columns = ['chrom', 'position', 'ref', 'alt']
                sites = raw
            else:
                raise ValueError(
                    'Unrecognized markers format. Expected either contig/position/ref_allele/alt_allele, '
                    'chrom/position/ref/alt, or a 5-column chrom/start/end/strand/ref_alt table.'
                )

        sites = sites[['chrom', 'position', 'ref', 'alt']].copy()
        sites = sites.dropna()
        sites['position'] = sites['position'].astype(int)
        sites['ref'] = sites['ref'].str.upper()
        sites['alt'] = sites['alt'].str.upper()
        # Keep SNPs only.
        sites = sites[(sites['ref'].str.len() == 1) & (sites['alt'].str.len() == 1)]
        return sites

    def load_het_site_counts(self, contig, position, ref_allele, alt_allele):
        t_ref_count, t_alt_count = self.read_pileup_return_count(self.tumor_bam,
                                                                  contig,
                                                                  position,
                                                                  ref_allele,
                                                                  alt_allele,
                                                                  MIN_MQ=10,
                                                                  MIN_BQ=10)
        n_ref_count, n_alt_count = self.read_pileup_return_count(self.normal_bam,
                                                                  contig,
                                                                  position,
                                                                  ref_allele,
                                                                  alt_allele,
                                                                  MIN_MQ=10,
                                                                  MIN_BQ=10)
        self.output['contig'].append(contig)
        self.output['position'].append(position)
        self.output['ref_allele'].append(ref_allele)
        self.output['alt_allele'].append(alt_allele)
        self.output['t_ref_count'].append(t_ref_count)
        self.output['t_alt_count'].append(t_alt_count)
        self.output['n_ref_count'].append(n_ref_count)
        self.output['n_alt_count'].append(n_alt_count)
    
    def read_pileup_return_count(self, samfile, chrom, pos, ref,
                                alt, MIN_MQ=10, MIN_BQ=10,
                                testing=False):
        '''
        self.t_ref_count = row['t_ref_count']
        self.t_alt_count = row['t_alt_count']
        self.n_ref_count = row['n_ref_count']
        self.n_alt_count = row['n_alt_count']
        '''
        ref_reads = set()
        alt_reads = set()
        other_reads = set()
        pileup = self.fetch_pileups(chrom, pos, samfile)
        ref_len = len(ref)
        alt_len = len(alt)
        for pileupcolumn in pileup:
            if pileupcolumn.pos != pos - 1:
                continue
            for pileupread in pileupcolumn.pileups:
                pos_in_read = pileupread.query_position
                if pos_in_read is None:
                    continue

                aln = pileupread.alignment
                if testing:
                    if aln.is_duplicate:
                        print('is_duplicate')
                        sys.exit(0)
                    if aln.is_qcfail:
                        print('is_qcfail')
                        sys.exit(0)

                if aln.mapping_quality < MIN_MQ \
                        or aln.is_supplementary \
                        or aln.is_secondary:
                    continue

                base_qualities = aln.query_qualities
                if base_qualities is None or base_qualities[pos_in_read] < MIN_BQ:
                    continue

                seq = aln.query_sequence
                if seq is None:
                    continue

                read_name = aln.query_name
                # Compare only small local slices at the current read position.
                if seq[pos_in_read:pos_in_read + ref_len] == ref:
                    ref_reads.add(read_name)
                elif seq[pos_in_read:pos_in_read + alt_len] == alt:
                    alt_reads.add(read_name)
                else:
                    other_reads.add(read_name)
            break
        # report no ref is observed
        if len(ref_reads) == 0 and len(other_reads) > 0:
            print('WARNING: No reads supporting reference allele {0} at {1}:{2}.'.format(ref, chrom, pos))
            print(alt_reads, other_reads)
        # check sets to make sure reads don't show up in multiple sets
        # supporting multiple calls
        ref_reads_set = ref_reads - alt_reads - other_reads
        alt_reads_set = alt_reads - ref_reads - other_reads
        # tally set in ref and alt, non-ref/alt, all reads
        ref_count = len(ref_reads_set)
        alt_count = len(alt_reads_set)
        return ref_count, alt_count

    def close(self):
        self.tumor_bam.close()
        self.normal_bam.close()
        

def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-id',
                        help='Sample ID',
                        required=False
                       )
    parser.add_argument('--tumor-bam',
                        help='Tumor BAM or CRAM file',
                        required=True
                       )
    parser.add_argument('--normal-bam',
                        help='Normal BAM or CRAM file',
                        required=True
                       )
    parser.add_argument('--het-site-calls-out',
                        help='''The het sites should be in a TSV file, with the columns: contig, position, ref_allele, alt_allele, t_ref_count, 
                        t_alt_count, n_ref_count, n_alt_count. Also the columns `total_reads` which contains 
                        the number of reads before filtering, and `map_Q0_reads`. If present (in the VCF), these are also
                        used to filter het-sites.''',
                        required=True
                       )
    parser.add_argument('--het-sites',
                        help='Heterozygous marker positions file',
                        required=True
                       )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    '''
        Returns counts for a same of heterozygous sites in tumor and normal BAM files. 
        The het sites should be in a TSV file, with the columns: contig, position, 
        ref_allele, alt_allele, t_ref_count, t_alt_count, n_ref_count, n_alt_count
    '''
    args = get_args()
    hs = HetSite(args['het_sites'],
                 args['tumor_bam'],
                 args['normal_bam'])
    try:
        output_df = pd.DataFrame(hs.output)
        output_df = output_df[['contig', 'position', 'ref_allele', 'alt_allele',
                               't_ref_count', 't_alt_count', 'n_ref_count', 'n_alt_count']]
        output_df.to_csv(args['het_site_calls_out'],
                         sep='\t', index=False)
    finally:
        hs.close()


if __name__ == '__main__':
    main()