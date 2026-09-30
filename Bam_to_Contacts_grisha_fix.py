#!/usr/bin/env python
import pysam
import argparse
import sys

def process_bam(filename):
    return pysam.AlignmentFile(filename, "rb").fetch(until_eof=True)

def extract_nm_tag(read):
    return str(read.get_tag('NM')) if read and read.has_tag('NM') else '*'

def extract_other_tags(read, other_tags):
    if not read:
        return '*'
    tags = []
    # Use try/except on get_tag to avoid two C-calls (has_tag + get_tag) per tag.
    # get_tag will raise KeyError if tag is absent (pysam behavior).
    for tag in other_tags:
        try:
            value = read.get_tag(tag)
        except (KeyError, AttributeError):
            continue
        tags.append(f'{tag}:"{value}"')
    return ';'.join(tags) if tags else '*'

def format_alignment(read):
    # Return minimal typed values; convert to str() only when writing output
    if read.is_unmapped:
        return ['*'] * 7  # chr, start, end, strand, cigar, NM, mapq
    # Cache properties and minimize C->Python crossings
    ref = read.reference_name
    start = read.reference_start + 1  # keep int (1-based)
    end = read.reference_end          # int (1-based)
    strand = '-' if read.is_reverse else '+'
    cigar = read.cigarstring
    # Try to get NM once; this avoids an extra has_tag() call per read.
    try:
        nm = read.get_tag('NM')
    except (KeyError, AttributeError):
        nm = '*'
    mapq = read.mapping_quality
    return [ref, start, end, strand, cigar, nm, mapq]

def format_secondary_alignment(read):
    strand = '-' if read.is_reverse else '+'
    nm = read.get_tag('NM') if read.has_tag('NM') else '*'
    # keep start as int here; caller may stringify
    return f"({read.reference_name},{read.reference_start},{strand}{read.cigarstring},{nm})"

def get_mapping_status(alignments, secondary_alignments):
    # alignments expected to be a list where format_alignment returns ['*']*7 for unmapped reads
    if alignments and alignments[0] == '*':
        return 'N'
    return 'M' if secondary_alignments else 'U'

def get_base_name(query_name):
    if not query_name:
        return None
    i = query_name.rfind('.')
    return query_name[i+1:] if i != -1 else query_name

def process_reads(r1_iter, r2_iter, other_tags, rna_mode='STAR', dna_mode='STAR'):
    # 1. SUPPLEMENTARY READS FILTRATION
    r1_iter_filtered = (r for r in r1_iter if not r.is_supplementary)  # Generator for RNA
    r2_iter_filtered = (r for r in r2_iter if not r.is_supplementary)  # Generator for DNA
    # Local bindings to speed up hot loop (reduce global lookups)
    fmt_secondary = format_secondary_alignment
    get_base = get_base_name
    # Convert other_tags to a set for fast membership tests and precompute
    # the sets of wanted keys for primary-read tag scans. This avoids
    # rebuilding the set on every read and allows early exit when all
    # required tags are found.
    set_other_tags = set(other_tags) if other_tags else set()
    wanted_primary_bwa = set_other_tags.union({'NM', 'XA'})
    wanted_primary_nonbwa = set_other_tags.union({'NM'})
    # 2. INITIALIZING THE FIRST READS
    #r1, r2 = next(r1_iter, None), next(r2_iter, None)
    r1, r2 = next(r1_iter_filtered, None), next(r2_iter_filtered, None)

    while r1 is not None or r2 is not None:
        query_name = r1.query_name if r1 else r2.query_name
        base_name = get_base(query_name)
        # r?_data: [alignments, secondary_alignments_list_or_None, other_tags_str]
        r1_data = [None, None, '*']
        r2_data = [None, None, '*']

        # Process R1 reads (RNA)
        while r1:
            # Avoid repeated get_base_name() calls: compare suffix using
            # endswith('.' + base_name) or full equality. endswith is
            # implemented in C and is faster than repeated rfind calls.
            r1_qname = r1.query_name
            if not (r1_qname == base_name or r1_qname.endswith('.' + base_name)):
                break
            # BWA: use XA tag for secondary alignments and treat primary alignments normally
            if rna_mode == 'BWA':
                if r1.is_unmapped:
                    r1_data[0] = ['*'] * 7
                    r1_data[2] = '*'
                else:
                    # single-pass scan of tags (avoid allocating a dict)
                    nm = '*'
                    xa = None
                    ot = []
                    if hasattr(r1, 'get_tags'):
                        remaining = len(wanted_primary_bwa)
                        for k, v in r1.get_tags():
                            if k == 'NM':
                                nm = v
                                if 'NM' in wanted_primary_bwa:
                                    remaining -= 1
                            elif k == 'XA':
                                xa = v
                                if 'XA' in wanted_primary_bwa:
                                    remaining -= 1
                            elif k in set_other_tags:
                                ot.append(f'{k}:"{v}"')
                                remaining -= 1
                            if remaining <= 0:
                                break
                    ref = r1.reference_name
                    start = r1.reference_start + 1
                    end = r1.reference_end
                    strand = '-' if r1.is_reverse else '+'
                    cigar = r1.cigarstring
                    mapq = r1.mapping_quality
                    r1_data[0] = [ref, start, end, strand, cigar, nm, mapq]
                    r1_data[2] = ';'.join(ot) if ot else '*'
                    if xa:
                        if r1_data[1] is None:
                            r1_data[1] = []
                        r1_data[1].append(xa)
            else:
                # non-BWA (STAR/HISAT): secondary reads are marked by is_secondary,
                # primary reads should be formatted like format_alignment and other_tags collected
                if r1.is_secondary:
                    if r1_data[1] is None:
                        r1_data[1] = []
                    r1_data[1].append(fmt_secondary(r1))
                else:
                    if r1.is_unmapped:
                        r1_data[0] = ['*'] * 7
                        r1_data[2] = '*'
                    else:
                        # single-pass tag scan for non-BWA primary reads
                        nm = '*'
                        ot = []
                        if hasattr(r1, 'get_tags'):
                            remaining = len(wanted_primary_nonbwa)
                            for k, v in r1.get_tags():
                                if k == 'NM':
                                    nm = v
                                    if 'NM' in wanted_primary_nonbwa:
                                        remaining -= 1
                                elif k in set_other_tags:
                                    ot.append(f'{k}:"{v}"')
                                    remaining -= 1
                                if remaining <= 0:
                                    break
                        ref = r1.reference_name
                        start = r1.reference_start + 1
                        end = r1.reference_end
                        strand = '-' if r1.is_reverse else '+'
                        cigar = r1.cigarstring
                        mapq = r1.mapping_quality
                        r1_data[0] = [ref, start, end, strand, cigar, nm, mapq]
                        r1_data[2] = ';'.join(ot) if ot else '*'
            r1 = next(r1_iter_filtered, None)
 
        # Process R2 reads (DNA)
        while r2:
            r2_qname = r2.query_name
            if not (r2_qname == base_name or r2_qname.endswith('.' + base_name)):
                break
            # BWA branch
            if dna_mode == 'BWA':
                if r2.is_unmapped:
                    r2_data[0] = ['*'] * 7
                    r2_data[2] = '*'
                else:
                    # single-pass scan of tags for BWA
                    nm = '*'
                    xa = None
                    ot = []
                    if hasattr(r2, 'get_tags'):
                        remaining = len(wanted_primary_bwa)
                        for k, v in r2.get_tags():
                            if k == 'NM':
                                nm = v
                                if 'NM' in wanted_primary_bwa:
                                    remaining -= 1
                            elif k == 'XA':
                                xa = v
                                if 'XA' in wanted_primary_bwa:
                                    remaining -= 1
                            elif k in set_other_tags:
                                ot.append(f'{k}:"{v}"')
                                remaining -= 1
                            if remaining <= 0:
                                break
                    ref = r2.reference_name
                    start = r2.reference_start + 1
                    end = r2.reference_end
                    strand = '-' if r2.is_reverse else '+'
                    cigar = r2.cigarstring
                    mapq = r2.mapping_quality
                    r2_data[0] = [ref, start, end, strand, cigar, nm, mapq]
                    r2_data[2] = ';'.join(ot) if ot else '*'
                    if xa:
                        if r2_data[1] is None:
                            r2_data[1] = []
                        r2_data[1].append(xa)
            else:
                # non-BWA (STAR/HISAT): secondary reads are handled via format_secondary_alignment
                if r2.is_secondary:
                    if r2_data[1] is None:
                        r2_data[1] = []
                    r2_data[1].append(fmt_secondary(r2))
                else:
                    if r2.is_unmapped:
                        r2_data[0] = ['*'] * 7
                        r2_data[2] = '*'
                    else:
                        # single-pass tag scan for non-BWA primary reads
                        nm = '*'
                        ot = []
                        if hasattr(r2, 'get_tags'):
                            remaining = len(wanted_primary_nonbwa)
                            for k, v in r2.get_tags():
                                if k == 'NM':
                                    nm = v
                                    if 'NM' in wanted_primary_nonbwa:
                                        remaining -= 1
                                elif k in set_other_tags:
                                    ot.append(f'{k}:"{v}"')
                                    remaining -= 1
                                if remaining <= 0:
                                    break
                        ref = r2.reference_name
                        start = r2.reference_start + 1
                        end = r2.reference_end
                        strand = '-' if r2.is_reverse else '+'
                        cigar = r2.cigarstring
                        mapq = r2.mapping_quality
                        r2_data[0] = [ref, start, end, strand, cigar, nm, mapq]
                        r2_data[2] = ';'.join(ot) if ot else '*'
            r2 = next(r2_iter_filtered, None)

        yield query_name, r1_data, r2_data

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process BAM files and output alignment information.')
    parser.add_argument('-r1', '--rna_bam', required=True, help='Input RNA BAM file (R1)')
    parser.add_argument('-r2', '--dna_bam', help='Input DNA BAM file (R2)')
    parser.add_argument('-mr', '--rna_mode', choices=['BWA', 'STAR', 'HISAT'], required=True, 
                        help='RNA alignment mode: BWA, STAR, or HISAT')
    parser.add_argument('-md', '--dna_mode', choices=['BWA', 'STAR', 'HISAT'], 
                        help='DNA alignment mode: BWA, STAR, or HISAT (defaults to RNA mode if not specified)')
    parser.add_argument('-e', '--exp_type', choices=['OTA_PE', 'OTA_SE', 'ATA', 'RNA_SEQ_PE', 'RNA_SEQ_SE'], required=True,
                        help='Experiment type: OTA_PE, OTA_SE, ATA, RNA_SEQ_PE, or RNA_SEQ_SE')
    parser.add_argument('-t', '--other_tags', nargs='+', default=["NH"],
                        help='List of additional SAM tags to extract (default: NH)')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Output file prefix')
    return parser.parse_args()

def write_header(file, exp_type):
    if exp_type == 'ATA':
        header = ["read_id", "ATA_pairtype", 
                "rna_chr", "rna_start", "rna_end", "rna_strand", "rna_cigar", "rna_NM", "rna_mapq",
                "dna_chr", "dna_start", "dna_end", "dna_strand", "dna_cigar", "dna_NM", "dna_mapq",
                "rna_secondary_alignments", "dna_secondary_alignments",
                "rna_other_tags", "dna_other_tags"]
    elif exp_type == 'OTA_PE':
        header = ["read_id", "OTA_PE_pairtype", 
                "dna1_chr", "dna1_start", "dna1_end", "dna1_strand", "dna1_cigar", "dna1_NM", "dna1_mapq",
                "dna2_chr", "dna2_start", "dna2_end", "dna2_strand", "dna2_cigar", "dna2_NM", "dna2_mapq",
                "dna1_secondary_alignments", "dna2_secondary_alignments",
                "dna1_other_tags", "dna2_other_tags"]
    elif exp_type == 'OTA_SE':          
        header = ["read_id", "OTA_SE_pairtype", 
                "dna1_chr", "dna1_start", "dna1_end", "dna1_strand", "dna1_cigar", "dna1_NM", "dna1_mapq",
                "dna2_chr", "dna2_start", "dna2_end", "dna2_strand", "dna2_cigar", "dna2_NM", "dna2_mapq",
                "dna1_secondary_alignments", "dna2_secondary_alignments",
                "dna1_other_tags", "dna2_other_tags"]
    elif exp_type == 'RNA_SEQ_PE':          
        header = ["read_id", "RNAseq_PE_pairtype", 
                "rna1_chr", "rna1_start", "rna1_end", "rna1_strand", "rna1_cigar", "rna1_NM", "rna1_mapq",
                "rna2_chr", "rna2_start", "rna2_end", "rna2_strand", "rna2_cigar", "rna2_NM", "rna2_mapq",
                "rna1_secondary_alignments", "rna2_secondary_alignments",
                "rna1_other_tags", "rna2_other_tags"]
    elif exp_type == 'RNA_SEQ_SE':          
        header = ["read_id", "RNAseq_SE_pairtype", 
                "rna1_chr", "rna1_start", "rna1_end", "rna1_strand", "rna1_cigar", "rna1_NM", "rna1_mapq",
                "rna2_chr", "rna2_start", "rna2_end", "rna2_strand", "rna2_cigar", "rna2_NM", "rna2_mapq",
                "rna1_secondary_alignments", "rna2_secondary_alignments",
                "rna1_other_tags", "rna2_other_tags"]
              
    file.write('\t'.join(header) + '\n')

def are_bams_queryname_sorted(file1, file2):
    """Checks the sorting of two BAM files by queryname"""
    for filename in [file1, file2]:
        try:
            with pysam.AlignmentFile(filename, "rb") as bam:
                header = bam.header
                if not ('HD' in header and 'SO' in header['HD'] and header['HD']['SO'] == 'queryname'):
                    return False
        except:
            return False
    return True


def count_unique_read_ids(filename):
    """Counts unique read IDs in a BAM file. Assumes the file is sorted by queryname."""
    unique_count = 0
    prev_read_id = None
    try:
        with pysam.AlignmentFile(filename, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                current_read_id = read.query_name
                if prev_read_id is None or current_read_id != prev_read_id:
                    unique_count += 1
                    prev_read_id = current_read_id
        return unique_count
    except Exception as e:
        return -1  # Return -1 if an error occurs.
    

def main():
    args = parse_arguments()
    
    # Set default DNA mode to RNA mode if not specified
    if not args.dna_mode:
        args.dna_mode = args.rna_mode
    
    try:
        unique_file = f"{args.prefix}_Unique_RNA.tab.rc"
        other_file = f"{args.prefix}_Other.tab.rc"

        # For RNA-seq SE / OTA_SE we do NOT want to process the same BAM twice.
        # Keep processing_exp_type as requested but avoid setting args.dna_bam = args.rna_bam
        # (that would cause redundant work). For SE experiments we'll pass an empty
        # iterator as r2 to `process_reads` below.
        if args.exp_type in ['RNA_SEQ_SE','OTA_SE']:
            processing_exp_type = args.exp_type
        elif args.exp_type == 'RNA_SEQ_PE':
            if not args.dna_bam:
                # args.dna_bam = args.rna_bam  # Use the same file for PE RNA-seq
                raise ValueError("BAM file (-r2) is required for PE RNA-seq experiments")
            processing_exp_type = 'RNA_SEQ_PE'
        else:
            processing_exp_type = args.exp_type
            if not args.dna_bam:
                raise ValueError("DNA BAM file (-r2) is required for RNA-DNA experiments")
        
        if processing_exp_type in ['ATA', 'OTA_PE', 'RNA_SEQ_PE']:
            bam_files_sorted = are_bams_queryname_sorted(args.rna_bam, args.dna_bam)
            if not bam_files_sorted:
                raise ValueError("Both BAM files must be sorted by query name for paired-end experiments.")
            else:
                count1 = count_unique_read_ids(args.rna_bam) #samtools?
                count2 = count_unique_read_ids(args.dna_bam)
                if (count1 != count2) or count1 == -1 or count2 == -1:
                    raise ValueError("Two BAM files have different numbers of unique read identifiers, but should have the same number.")
                
        # use a larger buffering size for output files to reduce syscalls
        # increase internal buffer to match our TARGET_BYTES for fewer syscalls
        with open(unique_file, 'w', buffering=262144) as f_unique, open(other_file, 'w', buffering=262144) as f_other:
            write_header(f_unique, processing_exp_type)
            write_header(f_other, processing_exp_type)

            # Prepare iterators: for SE experiments we avoid opening the DNA BAM
            # and pass an empty iterator for r2 to prevent redundant work.
            r1_iter = process_bam(args.rna_bam)
            if processing_exp_type in ['RNA_SEQ_SE', 'OTA_SE']:
                r2_iter = iter(())
            else:
                r2_iter = process_bam(args.dna_bam)

            # batched writing to reduce Python->syscall overhead
            # flush after a fixed number of lines (deterministic)
            FLUSH_LINES = 1000
            buf_unique = []
            buf_other = []

            for query_name, r1_data, r2_data in process_reads(
                r1_iter,
                r2_iter,
                args.other_tags,
                args.rna_mode,
                args.dna_mode
            ):
                # compute r1 status normally (r?_data: [alignments, secondary_list_or_None, other_tags_str])
                r1_status = get_mapping_status(r1_data[0], r1_data[1])
                # for single-end experiments (SE) we explicitly set r2_status to ''
                if processing_exp_type in ['RNA_SEQ_SE', 'OTA_SE']:
                    r2_status = ''  # No r2 data for SE experiments, so status is empty string
                else:
                    r2_status = get_mapping_status(r2_data[0], r2_data[1])
                pairtype = f"{r1_status}{r2_status}"

                r1_alignments = r1_data[0] if r1_data[0] else ['*'] * 7
                r2_alignments = r2_data[0] if r2_data[0] else ['*'] * 7

                # convert ints -> str only once here
                r1_alignments_str = [str(x) for x in r1_alignments]
                r2_alignments_str = [str(x) for x in r2_alignments]

                r1_secondary = ';'.join(r1_data[1]) if r1_data[1] else '*'
                r2_secondary = ';'.join(r2_data[1]) if r2_data[1] else '*'

                r1_other = r1_data[2] if r1_data[2] else '*'
                r2_other = r2_data[2] if r2_data[2] else '*'

                line = [
                    query_name, pairtype,
                    *r1_alignments_str,
                    *r2_alignments_str,
                    r1_secondary, r2_secondary,
                    r1_other, r2_other
                ]

                output_line = '\t'.join(line) + '\n'

                s = output_line
                if pairtype in ['UU', 'UM', 'U']:
                    buf_unique.append(s)
                    if len(buf_unique) >= FLUSH_LINES:
                        f_unique.write(''.join(buf_unique))
                        buf_unique.clear()
                else:
                    buf_other.append(s)
                    if len(buf_other) >= FLUSH_LINES:
                        f_other.write(''.join(buf_other))
                        buf_other.clear()

            # flush remaining buffers
            if buf_unique:
                f_unique.write(''.join(buf_unique))
            if buf_other:
                f_other.write(''.join(buf_other))

    except IOError as e:
        print(f"Error processing files: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    # Optional line-by-line profiling with line_profiler:
    # Run as: python Bam_to_Contacts_grisha_fix.py --profile [other args...]
    if '--profile' in sys.argv:
        try:
            from line_profiler import LineProfiler
        except ImportError:
            print("line_profiler is not installed. Install with: pip install line_profiler", file=sys.stderr)
            sys.exit(1)
        # --profile requires --profile-output <file>
        if '--profile-output' not in sys.argv:
            print("When using --profile you must provide --profile-output <file>", file=sys.stderr)
            sys.exit(1)
        # extract profile output path
        po_idx = sys.argv.index('--profile-output')
        try:
            profile_out = sys.argv[po_idx + 1]
        except Exception:
            print("Missing value for --profile-output <file>", file=sys.stderr)
            sys.exit(1)
        # remove the flags/values so argparse in main() won't see them
        for token in ['--profile', '--profile-output', profile_out]:
            if token in sys.argv:
                sys.argv.remove(token)

        lp = LineProfiler()
        # register functions to profile
        try:
            lp.add_function(process_bam)
            lp.add_function(process_reads)
            lp.add_function(format_alignment)
            lp.add_function(format_secondary_alignment)
            lp.add_function(extract_other_tags)
            lp.add_function(extract_nm_tag)
            lp.add_function(get_mapping_status)
        except Exception:
            # If any name isn't available for some reason, ignore and proceed
            pass

        # run main under the profiler and write results to file
        lp.runcall(main)
        try:
            # write profile output with buffering to reduce I/O overhead
            with open(profile_out, 'w', buffering=65536) as f:
                lp.print_stats(stream=f)
        except Exception as e:
            print(f"Failed to write profile output to {profile_out}: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        main()