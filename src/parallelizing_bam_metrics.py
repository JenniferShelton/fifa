#!/usr/bin/env python

##########################################################################
#    USAGE: Extract features from BAM file for a set of variants in a VCF file,
#    with most-up-to-date paralellization scheme
##########################################################################

################################################################################
### MODULES ####################################################################

import pysam 
import pysamstats
import numpy as np # type: ignore
import pandas as pd # type: ignore
import csv
import os, sys
import re
import queue as queue_module
from scipy.stats import entropy # type: ignore
from Bio.SeqUtils import gc_fraction # type: ignore
from concurrent.futures import ProcessPoolExecutor
import traceback
import subprocess
from multiprocessing import Queue, Lock, Process, Pool, Manager
import time
import logging
import warnings
from metrics_dictionary import MetricsDictionary
from metrics_dictionary import safe_median
from metrics_dictionary import to_pyrimidine_context
import unicodedata

################################################################### /MODULES ###
################################################################################

global logger

def is_read_filtered(read):
    """
    Check if read passes
    our immitation of the Mutect2 Filters
    """
    
    mutect_filtered = False

    # Check SAM flags
    if read.is_unmapped or read.is_secondary or read.is_duplicate:
        mutect_filtered = True
    elif not read.is_proper_pair or read.mate_is_unmapped:
        mutect_filtered = True

    # Check read properties
    elif read.query_length <= 1 or read.mapq <= 10 or read.mapq == 255 or read.has_tag('SA') or not read.cigar:
        mutect_filtered = True

    # Check CIGAR string
    elif not re.fullmatch(r'(?=.*[MD=X])([0-9]+[MIDSHP=X])+', read.cigarstring) or \
         re.search(r'(D|I)[0-9]*(?=(D|I))', read.cigarstring) or \
         re.match(r'^[0-9]*[SH]*[0-9]*D', read.cigarstring) or re.match(r'D[0-9]*[SH]*$', read.cigarstring):
        mutect_filtered = True

    # Check alignment conditions
    elif read.query_alignment_start < 0 or read.query_alignment_end < read.query_alignment_start or \
         read.infer_read_length() <= 0:
        mutect_filtered = True

    return mutect_filtered

def process_cigar_tupples(read, reference_pos):
    """
    Process cigar string to get distance to effective 5' and 3' end
    and position of the variant in the read (both with and without
    counting clipped bases). 
    Iterates through cigar string operators, while keeping track of the corresponding
    position in the reference sequence. 

    Args:
        read (pysam.AllignedSegment)
        reference_pos (int): position of the variant on the reference sequence

    Returns:
        int: position of the variant in the read, excluding hard clipped bases
        int: position of the variant in the read, including hard clipped bases
        int: distance from variant to read effective 5' end
        int: distance from variant to read effective 3' end
        int: length of clipped bases (soft and hard clipped)
        
        Returns None if read is missing cigarstring or if the variant is a deletion
    """    
    if not(read.cigarstring):
        return (None, None, None, None)
    
    ref_pos_in_read = read.reference_start
    real_position_in_read = 0 # position of the variant in the read, including hard clipped bases
    distance_to_5prime = distance_to_3prime = 0 
    
    clipped_length = 0
    
    for operation, length in read.cigartuples:
        if operation == 0: # Match or mismatch
            if ref_pos_in_read <= reference_pos:  
                if reference_pos < ref_pos_in_read + length:
                    real_position_in_read += reference_pos - ref_pos_in_read
                    distance_to_5prime += reference_pos - ref_pos_in_read
                    distance_to_3prime += length - (reference_pos - ref_pos_in_read)
                else: 
                    real_position_in_read += length
                    distance_to_5prime += length
            else:
                distance_to_3prime += length 
            ref_pos_in_read += length

        elif operation == 1:  # Insertion
            if ref_pos_in_read <= reference_pos:
                distance_to_5prime += length
                real_position_in_read += length
            else:
                distance_to_3prime += length 
        
        elif operation == 2:  # Deletion
            if ref_pos_in_read <= reference_pos < ref_pos_in_read + length:
                return (None, None, None, None) 
            
            ref_pos_in_read += length
        
        elif operation == 4:  # Soft clipping
            clipped_length += length
            if ref_pos_in_read <= reference_pos: 
                real_position_in_read += length

        elif operation == 5:  # Hard clipping
            clipped_length += length
            if ref_pos_in_read <= reference_pos: 
                real_position_in_read += length
    
    return real_position_in_read, distance_to_5prime, distance_to_3prime, clipped_length

def get_mismatch_and_insertion_positions(read):
    """
    Counts total number of insertions plus mismatches, 
    and extracts base qualities at those positions with 
    pysam.AllignedSegment.query_qualities

    Returns:
        int: number of mismatches in the read
        list: base qualities of bases where a mismatch occurs. 
    """    
    mismatch_positions = []
    current_position = 0
    
    # Parse the MD tag for mismatches
    md_tag = read.get_tag('MD')
    i = 0
    while i < len(md_tag):
        if md_tag[i].isdigit():
            num = ''
            while i < len(md_tag) and md_tag[i].isdigit():
                num += md_tag[i]
                i += 1
            current_position += int(num)
        elif md_tag[i] == '^':
            while i < len(md_tag) and not md_tag[i].isdigit():
                i += 1
        else:
            mismatch_positions.append(current_position)
            current_position += 1
            i += 1
            
    mismatch_base_quals = [read.query_qualities[pos] for pos in mismatch_positions] if mismatch_positions else None
    
    # Parse the CIGAR string to count insertions in the total mismatches count
    cigar = read.cigarstring
    num_mismatches = len(mismatch_positions)
    i = 0
    while i < len(cigar):
        num = ''
        while cigar[i].isdigit() and i < len(cigar):
            num += cigar[i]
            i += 1
        if cigar[i] == 'I':
            for x in range(1, int(num)+1):
                num_mismatches += 1
            i += 1
        else:
            i += 1
    return num_mismatches, mismatch_base_quals


################################################################################
###GENERATE WINDOW-BASED METRICS################################################

def get_coverage_in_window(bamfile, fastafile, chrom, left, right):
    """
    Extract coverage using BAM file for a specified window in the chromosome.
    Window is denoted by left and right endpoints, and cannot exceed the range of the 
    contig in the reference sequence.

    If left and right endpoints are outside of contig range, function select instead
    accordingly to the leftmost and rightmost ends of the contig, narrowing the window range. 

    Args:
        left (int): leftmost end of the window
        right (int): rightmost end of the window

    Returns:
        int: median coverage in the specified window
        float: standard deviation of coverege in thewindow
    """    
    window_coverage = pysamstats.load_pileup(
        type="coverage_ext", alignmentfile=bamfile, fafile=fastafile, 
        chrom=chrom, start=left, end=right, truncate=False
    )
    
    # Ensure that left and right endpoints meet contig ranges
    indeces = np.where((window_coverage.pos >= left) & (window_coverage.pos <= right))[0]
    index_left, index_right = (indeces[0], indeces[-1]) if indeces.size > 0 else (0, 0)
    
    window_data = window_coverage[index_left:index_right].reads_all
    
    if window_data.size <= 1:
        return None, None
    
    return np.median(window_data), np.var(window_data)
    
def get_read_fractions(bamfile, chrom, left, right):
    """
    Extracts counts of reads that are duplicated,
    poor mapq, or failing other features within a
    specified window around the variant.
    Leftmost and rightmost positions denote the ends of 
    the basepair window. 

    Args:
        left (int): leftmost end of the window
        right (int): rightmost end of the window

    Returns:
        int: median fragment lengths in the window
        int: fraction of duplicated reads
        int: franction of reads that are improperly paired
        int: median mapq of reads in the window
        int: fraction of reads that are filtered by mutect2
    """    

    coverage = pysamstats.load_coverage(bamfile, chrom=chrom, start=left, end=right, truncate=True)
    num_total_reads = coverage['reads_all'].sum()
    ## In theory according to the documentation you should be able to get the num_mapq0, etc from 
    ## pysamstats, but for some reason I can't find those attributes in the recarray

    if num_total_reads == 0:
        return (0,0,0,0,0,0)
    
    num_improper_paired = num_total_reads - coverage['reads_pp'].sum()
    frag_lenths, mapqs = [], []
    num_mapq0 = num_duplicates = num_filtered_mutect = 0
    
    for read in bamfile.fetch(chrom, left, right):        
        if read.query_alignment_end > bamfile.get_reference_length(chrom) or is_read_filtered(read):
            num_filtered_mutect += 1
        else:
            frag_lenths.append(abs(read.template_length))
            mapqs.append(read.mapping_quality)
            
        num_mapq0 += read.mapping_quality == 0
        num_improper_paired += not read.is_proper_pair
        num_duplicates += read.is_duplicate
    
    return (
        safe_median(frag_lenths),
        num_duplicates / num_total_reads,
        num_mapq0 / num_total_reads,
        num_improper_paired / num_total_reads,
        safe_median(mapqs),
        num_filtered_mutect / num_total_reads
    )

################################################################################
###############################################GENERATE WINDOW-BASED METRICS####

#################################################################################
### SPLIT VARIANTS FOR PROCESSING ###############################################

def iter_pileup_columns(bamfile, chrom, pos):
    pileup_kwargs = {
        "contig": chrom,
        "start": pos - 1,
        "stop": pos,
        "min_base_quality": 10,
        "min_mapping_quality": 10,
        "multiple_iterators": False,
    }

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="multiple_iterators not implemented for CRAM",
            category=UserWarning,
        )
        try:
            return bamfile.pileup(**pileup_kwargs)
        except TypeError:
            pileup_kwargs.pop("multiple_iterators", None)
            return bamfile.pileup(**pileup_kwargs)


def process_variant(queue, sample, cohort, bam_path, ref_seq, iolock, final_dictionary):
    bamfile = pysam.AlignmentFile(bam_path, "rb", reference_filename=ref_seq)
    fastafile = pysam.FastaFile(ref_seq)

    while True:
        rec, sample_data, vcf_index = queue.get()
        
        if rec is None:
            bamfile.close()
            fastafile.close()
            break
        chrom = rec['CHROM']
        pos = rec['POS']
        ref = rec['REF']
        alt = rec['ALT']

        variant_id = '{0}:{1}_{2}>{3}'.format(chrom, pos, ref, alt)

        metrics = MetricsDictionary(cohort=cohort, sample=sample, index=vcf_index) 
        try :
            if 'Label' in rec.keys():
                metrics.set_metric('Label', rec['Label'])

            ad_values = sample_data.get('AD')
            af_values = sample_data.get('AF')
            dp_value = sample_data.get('DP')

            ad_is_missing = ad_values is None or (isinstance(ad_values, (list, tuple)) and len(ad_values) == 0)
            af_is_missing = af_values is None or (isinstance(af_values, (list, tuple)) and len(af_values) == 0)
            dp_is_missing = dp_value is None or (isinstance(dp_value, (list, tuple)) and len(dp_value) == 0)

            if not dp_is_missing:
                if isinstance(dp_value, (list, tuple)):
                    tumor_depth = sum(dp_value)
                else:
                    tumor_depth = int(dp_value)
                metrics.set_metric('tumor_depth', tumor_depth)
            else:
                metrics.set_metric('tumor_depth', 0)

            if not af_is_missing:
                if isinstance(af_values, (list, tuple)):
                    tumor_vaf = af_values[0] if af_values else 0
                else:
                    tumor_vaf = af_values
                metrics.set_metric('tumor_VAF', tumor_vaf)
            else:
                metrics.set_metric('tumor_VAF', 0)

            for pileupcolumn in iter_pileup_columns(bamfile, chrom, pos):
                if pileupcolumn.pos == pos - 1:
                    for pileupread in pileupcolumn.pileups:
                        metrics.increment_metric('num_total_reads')
                        query_index = pileupread.query_position
                        read = pileupread.alignment
                        
                        ## This is probably inefficient... there must be ways to do this with pysam without
                        ## the need for extra methods 
                        read_index, distance_5prime, distance_3prime, clipped_length = \
                            process_cigar_tupples(read, pos)

                        if query_index is None \
                                or query_index < 0 \
                                or read_index is None:
                            # Pretty sure that this is a problem when a read spans an indel
                            continue

                        base = read.query_sequence[query_index]
                        is_ref = base == ref
                        is_var = base == alt

                        
                        if is_read_filtered(read):
                            metrics.increment_metric('tumor_reads_filtered')

                        if is_ref or is_var: 
                            prefix='tumor_ref' if is_ref else 'tumor_var'
                            metrics.increment_metric('tumor_ref_count' if is_ref else 'tumor_var_count')
                            metrics.increment_metric(f'{prefix}_num_minus_strand' if read.is_reverse else f'{prefix}_num_plus_strand')
                            
                            metrics.add_metric(f'{prefix}_base_qualities', read.query_qualities[query_index])
                            metrics.add_metric(f'{prefix}_read_frag_length', read.infer_read_length())
                            metrics.add_metric(f'{prefix}_avg_pos_as_fraction', (read_index / (read.infer_read_length() / 2)))
                            metrics.add_metric(f'{prefix}_distances_to_5p_end', distance_5prime)
                            metrics.add_metric(f'{prefix}_distances_to_3p_end', distance_3prime)
                            metrics.add_metric(f'{prefix}_clipped_length', clipped_length)

                            if read.is_paired:
                                ## this is because we used to track the mapq of unpaired reads seperatly 
                                metrics.add_metric(f'{prefix}_mapping_quality', read.mapping_quality)

                            if read.has_tag('MD'):
                                num_mismatches, mismatch_base_quals = get_mismatch_and_insertion_positions(read)
                                metrics.add_metric('avg_num_mismatches', num_mismatches / read.infer_read_length())
                                if mismatch_base_quals:
                                    metrics.add_metric('avg_sum_mismatch_base_quals',(sum(mismatch_base_quals)))   
                        else:
                            metrics.increment_metric('tumor_other_bases_count')                   

            pileup_depth = metrics.get_metric('num_total_reads')
            if dp_is_missing and pileup_depth:
                metrics.set_metric('tumor_depth', pileup_depth)

            if af_is_missing:
                tumor_var_count = metrics.get_metric('tumor_var_count')
                tumor_depth_for_vaf = metrics.get_metric('tumor_depth') or pileup_depth
                if tumor_depth_for_vaf:
                    metrics.set_metric('tumor_VAF', tumor_var_count / tumor_depth_for_vaf)
                else:
                    metrics.set_metric('tumor_VAF', 0)
            metrics.aggregate_base_metrics(ref, alt) 

            ## Get Window-Based Metrics
            left = max(pos - 500, 0)
            right = min(pos + 500, bamfile.get_reference_length(chrom))

            left_window = [max(0, left - 500), left]
            right_window = [right, min(right + 500, bamfile.get_reference_length(chrom))]
            
            alignment_seq = fastafile.fetch(chrom, left, right)

            metrics.set_metric('window_gc_cont', gc_fraction(alignment_seq))
            metrics.set_metric('window_seq_entropy', entropy(np.unique(list(alignment_seq), return_counts=True)[1] / len(alignment_seq), base=2))
            
            median_cov, cov_variance = get_coverage_in_window(bamfile, fastafile, chrom, left, right)
            if median_cov is None or median_cov == 0:
                metrics.set_metric('window_min_cov_ratio', None)
                metrics.set_metric('window_max_cov_ratio', None)
                continue
            metrics.set_metric('window_median_cov', median_cov)
            metrics.set_metric('window_cov_variance', cov_variance)
            
            left_window_median_cov = get_coverage_in_window(bamfile, fastafile, chrom, *left_window)[0]
            right_window_median_cov = get_coverage_in_window(bamfile, fastafile, chrom, *right_window)[0]

            left_window_median_cov = left_window_median_cov if left_window_median_cov is not None else 0
            right_window_median_cov = right_window_median_cov if right_window_median_cov is not None else 0

            metrics.update_coverage_ratios(left=left_window_median_cov, right=right_window_median_cov)

            fractions = get_read_fractions(bamfile, chrom, left, right)
            metrics.set_metric('window_median_frag_len', fractions[0])
            metrics.set_metric('window_dup_frac', fractions[1])
            metrics.set_metric('window_multi_frac', fractions[2])
            metrics.set_metric('window_improper_frac', fractions[3])
            metrics.set_metric('window_median_mapq', fractions[4])
            metrics.set_metric('window_read_filter_frac', fractions[5])

            ## Extract flanking bases and format using universal convention:
            ## 5' flank(s) [Ref>Alt] 3' flank(s), reverse complemented to the pyrimidine convention
            sequence = fastafile.fetch(chrom, pos - 3, pos + 2)
            trinucleotide_context = f'{sequence[1]}[{ref}>{alt}]{sequence[3]}'
            pentanucleotide_context = f'{sequence[0:2]}[{ref}>{alt}]{sequence[3:5]}'
            metrics.set_metric('trinucleotide_context', to_pyrimidine_context(trinucleotide_context))
            metrics.set_metric('pentanucleotide_context', to_pyrimidine_context(pentanucleotide_context))
        
        except Exception as e:
            '''
            if the commands in the block fail the entire process will be terminated, and the error will be logged.
            '''
            logger.error(traceback.format_exc())
            logger.error(f"process_variant (very large) error trap: An error occurred: {e}")
            # metrics.set_metric('tumor_depth', 0)
            # metrics.set_metric('tumor_VAF', 0)
            # metrics.set_metric('window_gc_cont', 0)
            # metrics.set_metric('window_seq_entropy', 0)
            # metrics.set_metric('window_median_cov', 0)
            # metrics.set_metric('window_cov_variance', 0)
            # metrics.set_metric('window_min_cov_ratio', None)
            # metrics.set_metric('window_max_cov_ratio', None)
            # metrics.set_metric('trinucleotide_context', "")
            # metrics.set_metric('pentanucleotide_context', "")
            raise SystemExit(1)

        iolock.acquire()
        final_dictionary[variant_id] = metrics.get_all_metrics()
        iolock.release()

def read_vcf(sample, label, vcf_path, queue, num_threads, processes=None):
    def fail_if_worker_failed():
        if not processes:
            return
        failed_processes = [P for P in processes if P.exitcode not in (None, 0)]
        if failed_processes:
            failed_process = failed_processes[0]
            exit_code = failed_process.exitcode if failed_process.exitcode is not None else 1
            raise SystemExit(exit_code)

    vcffile = pysam.VariantFile(vcf_path) 
    for index, rec in enumerate(vcffile.fetch()):
        fail_if_worker_failed()
        if rec.ref not in ['A', 'C', 'T', 'G'] or rec.alts[0] not in ['A', 'C', 'T', 'G']:
            continue

        sample_data = dict(rec.samples[sample])

        rec_dict = {'CHROM' : str(rec.chrom),
               'POS' : int(rec.pos),
               'REF' : str(rec.ref),
               'ALT' : str(rec.alts[0])}
        
        if label[0] and label[1]:
            rec_dict['Label'] = 1 if rec.info.get(label[0]) == label[1] else 0

        queued = False
        while not queued:
            fail_if_worker_failed()
            try:
                queue.put((rec_dict, sample_data, index), timeout=0.2)
                queued = True
            except queue_module.Full:
                continue

    for i in range(int(num_threads)):
        queued = False
        while not queued:
            fail_if_worker_failed()
            try:
                queue.put((None, None, None), timeout=0.2)
                queued = True
            except queue_module.Full:
                continue

    vcffile.close()

def get_mobster_tail_scores(sample, vcf_path, out_path, mobster_scores, mobster_fit_rds=None):
    outfile = os.path.join(os.path.dirname(out_path), f'{sample}_mobster.csv')
    fitfile = os.path.join(os.path.dirname(out_path), f'{sample}_mobster_fit.rds')
    generated_fitfile = False

    script_path = os.path.abspath(__file__)
    directory_name = os.path.dirname(script_path)

    if mobster_fit_rds is not None:
        fitfile = mobster_fit_rds
        if not os.path.exists(fitfile):
            logger.error("Provided MOBSTER fit file does not exist: %s", fitfile)
            raise SystemExit(1)
    else:
        fit_process = subprocess.run(
            [directory_name + '/run_mobster_fit.R', sample, vcf_path, fitfile],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
            )

        logger.info(fit_process.stdout.decode('unicode_escape'))
        logger.error(fit_process.stderr.decode('unicode_escape'))
        if fit_process.returncode != 0:
            logger.error(
                "MOBSTER fit subprocess failed with exit code %s for sample %s",
                fit_process.returncode,
                sample,
            )
            raise SystemExit(fit_process.returncode)
        generated_fitfile = True

    sample_data_process = subprocess.run(
        [directory_name + '/run_mobster_sample_data.R', sample, vcf_path, fitfile, outfile],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
        )

    logger.info(sample_data_process.stdout.decode('unicode_escape'))
    logger.error(sample_data_process.stderr.decode('unicode_escape'))
    if sample_data_process.returncode != 0:
        logger.error(
            "MOBSTER sample-data subprocess failed with exit code %s for sample %s",
            sample_data_process.returncode,
            sample,
        )
        raise SystemExit(sample_data_process.returncode)
       
    with open(outfile, newline='') as mfile:
        logger.info(f"Finished MOBSTER calculations for {sample}")
        reader = csv.reader(mfile, delimiter=',')
        header = next(reader)
        for row in reader: 
            sample, chrom, pos, REF, ALT, Tail = row
            if REF not in ['A', 'C', 'T', 'G'] or ALT not in ['A', 'C', 'T', 'G']:
                continue
            variant_id = '{0}:{1}_{2}>{3}'.format(chrom, pos, REF, ALT)
            mobster_scores[variant_id] = {'Tail': Tail}
    mfile.close()
    if generated_fitfile and os.path.exists(fitfile):
        os.remove(fitfile)
    os.remove(outfile)

def extract_all_features(bam_path, vcf_path, ref_seq, sample, cohort, label, num_threads, output_file, skip_mobster=False, mobster_fit_rds=None):
    is_compressed = vcf_path.endswith('.gz')
    has_index = os.path.exists(vcf_path + '.csi') or os.path.exists(vcf_path + '.tbi')
    if not is_compressed:
        logger.error("VCF File is not compressed. Please bgzip the VCF file before running feature extraction.")
        exit(1)
    if not has_index:
        logger.error("VCF File is not indexed. Please index the VCF file before running feature extraction.")
        exit(1) 

    start = time.time()
    
    queue = Queue(maxsize=500) 

    iolock = Lock()
    to_df = []
    num_vars = 0

    with Manager() as manager:
        mobster_scores = manager.dict() if not skip_mobster else None
        all_features = manager.dict()

        pool = [Process(target=process_variant, args=(queue, sample, cohort, bam_path, ref_seq, iolock, all_features)) for i in range(int(num_threads))]
        if not skip_mobster:
            mobster_process = Process(target=get_mobster_tail_scores, args=(sample, vcf_path, output_file, mobster_scores, mobster_fit_rds), name="mobster_tail_scores")
            pool.insert(0, mobster_process)
        for P in pool:
            P.start()

        failed_process = None
        failed_exit_code = None
        try:
            read_vcf(sample, label, vcf_path, queue, num_threads, processes=pool)
        except SystemExit as e:
            failed_processes = [P for P in pool if P.exitcode not in (None, 0)]
            failed_process = failed_processes[0] if failed_processes else None
            failed_exit_code = failed_process.exitcode if failed_process is not None else (e.code if isinstance(e.code, int) else 1)
            logger.error(
                "Detected worker failure during VCF enqueue; terminating remaining workers. failed_exit_code: %s",
                failed_exit_code
            )
            for P in pool:
                if P.is_alive():
                    P.terminate()

        if failed_exit_code is None:
            while True:
                alive_processes = [P for P in pool if P.is_alive()]
                failed_processes = [P for P in pool if P.exitcode not in (None, 0)]

                if failed_processes:
                    failed_process = failed_processes[0]
                    logger.error(
                        "Worker %s failed with exit code %s; terminating remaining workers",
                        failed_process.name,
                        failed_process.exitcode,
                    )
                    for P in pool:
                        if P.is_alive():
                            P.terminate()
                    break

                if not alive_processes:
                    break

                time.sleep(0.2)

        for P in pool:
            P.join(timeout=5)
            if P.is_alive():
                logger.error("Worker %s did not terminate cleanly; force killing", P.name)
                try:
                    P.kill()
                except AttributeError:
                    P.terminate()
                P.join(timeout=5)

        if failed_process is not None:
            exit_code = failed_process.exitcode if failed_process.exitcode is not None else 1
            logger.error(
                "Feature extraction exiting with worker failure code %s from %s",
                exit_code,
                failed_process.name,
            )
            raise SystemExit(exit_code)
        if failed_exit_code is not None:
            logger.error("Feature extraction exiting with worker failure code %s", failed_exit_code)
            raise SystemExit(failed_exit_code)
        
        if skip_mobster:
            result = {variant: all_features.get(variant, {}) for variant in all_features.keys()}
        else:
            ## Takes care of cases when MOBSTER doesn't run succesfully
            result = {variant: {**all_features.get(variant, {}),**(mobster_scores.get(variant, {'Tail': 1}))}
                for variant in all_features.keys()}

        num_vars = len(all_features)
        to_df = [{'Variant': variant, **metric} for variant, metric in result.items()]
    
    pd.DataFrame(to_df).sort_values("vcf_index").drop("vcf_index", axis=1).to_csv(output_file, index=False)

    end = time.time() - start
    logger.info(f"Finished all metrics for {num_vars} vars in {sample} in {end}")
    logger.info(f"Features stored: \n{output_file}")

    ## Make a bunch of runs with different numbers of CPUs and see how it scales
    ## ask Jen about nice seff command
    ## Add num threads to the options 

#################################################################################
################################################ SPLIT VARIANTS for PROCESSING ##

#################################################################################
### PROCESS INPUT FILES #########################################################
    
def process_sample(sample, cohort, vcf_path, bam_path, ref_seq, output_file, label, num_threads, skip_mobster=False, mobster_fit_rds=None): 
    # try:
    if os.path.isfile(vcf_path) and os.path.isfile(bam_path) and os.path.isfile(ref_seq):
        logger.info(f'Processing BAM file for sample: {sample}')
        extract_all_features(bam_path, vcf_path, ref_seq, sample, cohort, label, num_threads, output_file, skip_mobster=skip_mobster, mobster_fit_rds=mobster_fit_rds)
            
    else:
        logger.error("There is an issue with one of your input files.")
        if not (os.path.isfile(vcf_path)):
            logger.error(f" The error is with your VCF File path: {vcf_path}")
        if not (os.path.isfile(bam_path)):
            logger.error(f" The error is with your BAM File path: {bam_path}")
        if not (os.path.isfile(ref_seq)):
            logger.error(f" The error is with your reference seq: {ref_seq}")
        raise FileNotFoundError("One or more input files are missing. Please verify that the paths are correct.")
    # except Exception as e:
    #     logger.error(f"process_sample error trap: Error processing {sample}: {e}")
    #     traceback.print_exc()

def process_bam_file(outpath, label, num_threads, sample, vcf_file, 
bam_file, ref_seq, cohort=None, skip_mobster=False, mobster_fit_rds=None):
    global logger
    logger = logging.getLogger(__name__)
    
    try:
        if os.path.isfile(outpath):
            output_file = outpath 
            
            if not (outpath.endswith('_extracted_features.csv')):
                logger.info("WARNING: If output path is a file, it must end with '_extracted_features.csv' due to dependencies "
                "in other parts of the code.")
                base_path, _ = os.path.splitext(outpath)
                output_file = base_path + "_extracted_features.csv"
                logger.info("New Output file : " + output_file)
        else:
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            output_file=os.path.join(outpath, f"{sample}_extracted_features.csv")
    except Exception as e:
        logger.error(f"process_bam_file error trap:An error occurred: {e}")
    process_sample(sample, cohort, vcf_file, bam_file, ref_seq, output_file, label, num_threads, skip_mobster=skip_mobster, mobster_fit_rds=mobster_fit_rds)
    # except Exception as e:
    #     logger.error(f"process_bam_file error trap:An error occurred: {e}")

#################################################################################
######################################################### PROCESS INPUT FILES ###