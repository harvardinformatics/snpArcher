# Predict, from a cSRA run's SEQUENCE table, exactly what download_sra's
# aligned-run extraction will emit:
#
#   sam-dump -u | samtools collate | samtools fastq -F 0x900 -1 -2 -0 -s
#
# Input is one line per spot, as written by
#
#   vdb-dump -T SEQUENCE -C READ_LEN,PRIMARY_ALIGNMENT_ID,ALIGNMENT_COUNT -f tab
#
# i.e. three tab-separated columns, each holding one comma-separated value per
# read in the spot. Output is key=value lines.
#
# sam-dump writes a read as a primary or unaligned record when it is non-empty
# and either has a primary alignment (the aligned pass) or sits in a spot where
# at least one read has ALIGNMENT_COUNT 0 (the unaligned pass visits only
# those spots: sra-tools 3.2.1, tools/external/sra-pileup/sam-unaligned.c,
# read_seq_row and dump_seq_prim_row_sam). Two kinds of stored read therefore
# never reach the FASTQ:
#
#   - zero-length reads. bam-load stores a mate-less BAM record as a two-read
#     spot whose second read is an empty technical read. The surviving read is
#     written without the paired flag, so samtools sends it to -0 (unpaired).
#
#   - reads whose only alignments are secondary (PRIMARY_ALIGNMENT_ID 0 with
#     ALIGNMENT_COUNT > 0) when the mate is aligned too. Both passes skip them;
#     they appear only as 0x100 copies, which -F 0x900 rightly drops, and the
#     primary mate is left as a singleton (-s).
#
# Spots holding anything but two reads are counted and not modelled; the rule
# refuses them rather than guess.

BEGIN { FS = "\t" }

{
    spots++
    nr = split($1, len, ", ")
    if (split($2, pid, ", ") != nr || split($3, ac, ", ") != nr) {
        malformed++
        next
    }
    if (nr != 2) {
        unmodelled++
        next
    }

    visited = (ac[1] == 0 || ac[2] == 0)
    nonempty = 0
    emitted = 0
    for (i = 1; i <= 2; i++) {
        if (pid[i] != 0) primary++
        if (len[i] == 0) { zero_length++; continue }
        nonempty++
        if (pid[i] != 0 || visited) emitted++
        else lost_secondary_only++
    }

    if (emitted == 2) pairs++
    else if (emitted == 1 && nonempty == 2) singletons++
    else if (emitted == 1) unpaired++
}

END {
    print "spots=" spots + 0
    print "spots_unmodelled=" unmodelled + 0
    print "lines_malformed=" malformed + 0
    print "reads_primary=" primary + 0
    print "reads_zero_length=" zero_length + 0
    print "reads_lost_secondary_only=" lost_secondary_only + 0
    print "pairs=" pairs + 0
    print "singletons=" singletons + 0
    print "unpaired=" unpaired + 0
}
