#!/usr/bin/env python

import pysam
import optparse

#Options
parser=optparse.OptionParser()
parser.add_option('-i', '--inbam', default = "",help='input bam file', type='str')
parser.add_option('-o', '--outfile', default= '', help='output file with breakpoints', type='str')
parser.add_option('-s', '--segments', default= 4, help='number of segments that constitute the correct positioning', type='int')

options, args=parser.parse_args()

# Open BAM file
bamfile = pysam.AlignmentFile(options.inbam, "rb")

def ref_span_from_cigar(cigarstring):
    ref_len = 0
    tmp = ''
    for c in cigarstring:
        if c.isdigit():
            tmp += c
        else:
            if c in 'MDN=X':
                ref_len += int(tmp)
            tmp = ''
    return ref_len

def leading_clip(cigarstring):
    tmp = ''
    for c in cigarstring:
        if c.isdigit():
            tmp += c
        else:
            if c in ('S', 'H'): #if c == 'S':
                return int(tmp)
            else:
                return 0
            tmp = ''
    return 0
"""
def merge_circular_segments(segments, chrom_sizes):
    merged = []
    i = 0

    while i < len(segments):
        cur = segments[i]
        lead_s, chrom, start, end, strand = cur

        if (
            i < len(segments) - 1
            and chrom in chrom_sizes
        ):
            next_seg = segments[i + 1]
            _, chrom2, start2, end2, strand2 = next_seg

            # check circular merge condition
            if (
                chrom == chrom2
                and strand == strand2
                and end == chrom_sizes[chrom]
                and start2 == 1
            ):
                # merge them
                merged.append((lead_s, chrom, start, end2, strand))
                i += 2
                continue

        merged.append(cur)
        i += 1

    return merged
 
chrom_sizes = {
    "pUC19": 2686
}
"""

discardedcounter = 0
totalcounter = 0

with open(options.outfile, 'w') as ouf:
	for read in bamfile:
	    # only process primary alignments to avoid repeats
	    if read.is_unmapped or read.is_secondary or read.is_supplementary:
	        continue
	
	    segments = []
	    totalcounter +=1
	
	    # --- Primary ---
	    ref_start = read.reference_start + 1
	    ref_end   = read.reference_end
	    strand    = '-' if read.is_reverse else '+'
	    cigar     = read.cigarstring
	    lead_s    = leading_clip(cigar)
	
	    segments.append((lead_s, read.reference_name, ref_start, ref_end, strand))
	
	    # --- SA tag ---
	    if read.has_tag('SA'):
	        sa_entries = read.get_tag('SA').strip().split(';')
	
	        for entry in sa_entries:
	            if not entry:
	                continue
	
	            fields = entry.split(',')
	            if len(fields) != 6:
	                continue
	
	            chrom, pos, strand_sa, cigar_sa, mapq, nm = fields
	            pos = int(pos)
	
	            ref_len = ref_span_from_cigar(cigar_sa)
	            ref_end_sa = pos + ref_len - 1
	            lead_s_sa = leading_clip(cigar_sa)
	
	            segments.append((lead_s_sa, chrom, pos, ref_end_sa, strand_sa))
	
	    # keep only reads with exactly 4 segments total
	    # it allows to only deal with the ones with obvious insertion of splint
	    if len(segments) != options.segments:
	        discardedcounter +=1
	        continue
	
	    # sort by read order (left clip)
	    segments.sort(key=lambda x: x[0])
	    #segments = merge_circular_segments(segments, chrom_sizes)
	
	    # output
	    out = [f"{chrom}:{start}-{end}" for _, chrom, start, end, strand in segments] #out = [f"{chrom}:{start}-{end}:{strand}" for _, chrom, start, end, strand in segments]
	    ouf.write(read.query_name+ '\t'+ '\t'.join(out) + '\n')
	
	
	print("Total processed reads: ", str(totalcounter))
	print("Discarded reads with incorrect number of segments: ", str(discardedcounter))
