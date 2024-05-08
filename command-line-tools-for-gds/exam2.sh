# Problems 1-5
# use athal_wu_0_A.bam

# Problem 1. How many alignments does the set contain?
NUM_ALIGNMENTS=$(samtools view gencommand_proj2_data/athal_wu_0_A.bam | wc -l)
echo "Number of alignments: $NUM_ALIGNMENTS"

# Problem 2. How many alignments show the read's mate unmapped?
NUM_UNMAPPED_ALIGNMENTS=$(samtools view gencommand_proj2_data/athal_wu_0_A.bam | cut -f7 | grep "*" | wc -l)
echo "Number of alignments: $NUM_UNMAPPED_ALIGNMENTS"

# Problem 3. How many alignments contain a deletion?
NUM_WITH_DELETION=$(samtools view gencommand_proj2_data/athal_wu_0_A.bam | cut -f6 | grep "D" | wc -l)
echo "Number of alignments with a deletion: $NUM_WITH_DELETION"

# Problem 4. How many alignments show the read’s mate mapped to the same chromosome?
NUM_MATE_SAME_CHR=$(samtools view gencommand_proj2_data/athal_wu_0_A.bam | cut -f7 | grep "=" | wc -l)
echo "Number of alignments showing the read's mate mapped to the same chromosome: $NUM_MATE_SAME_CHR"

# Problem 5. How many alignments are spliced?
NUM_SPLICED=$(samtools view gencommand_proj2_data/athal_wu_0_A.bam | cut -f6 | grep "N" | wc -l)
echo "Number of spliced alignments: $NUM_SPLICED"


# Problems 6-10
# Extract only the alignments in the range “Chr3:11,777,000-11,794,000”, corresponding to a locus of interest. For this alignment set:

# sort and index
samtools sort gencommand_proj2_data/athal_wu_0_A.bam gencommand_proj2_data/athal_wu_0_A.sorted
samtools index gencommand_proj2_data/athal_wu_0_A.sorted.bam

# Problem 6. How many alignments does the set contain?
NUM_ALIGNMENTS_2=$(samtools view gencommand_proj2_data/athal_wu_0_A.sorted.bam "Chr3:11,777,000-11,794,000" | wc -l)
echo "Number of alignments: $NUM_ALIGNMENTS_2"

# Problem 7. How many alignments show the read's mate unmapped?
NUM_UNMAPPED_ALIGNMENTS_2=$(samtools view gencommand_proj2_data/athal_wu_0_A.sorted.bam "Chr3:11,777,000-11,794,000" | cut -f7 | grep "*" | wc -l)
echo "Number of alignments: $NUM_UNMAPPED_ALIGNMENTS_2"

# Problem 8. How many alignments contain a deletion?
NUM_WITH_DELETION_2=$(samtools view gencommand_proj2_data/athal_wu_0_A.sorted.bam "Chr3:11,777,000-11,794,000" | cut -f6 | grep "D" | wc -l)
echo "Number of alignments with a deletion: $NUM_WITH_DELETION_2"

# Problem 9. How many alignments show the read’s mate mapped to the same chromosome?
NUM_MATE_SAME_CHR_2=$(samtools view gencommand_proj2_data/athal_wu_0_A.sorted.bam "Chr3:11,777,000-11,794,000" | cut -f7 | grep "=" | wc -l)
echo "Number of alignments showing the read's mate mapped to the same chromosome: $NUM_MATE_SAME_CHR_2"

# Problem 10. How many alignments are spliced?
NUM_SPLICED_2=$(samtools view gencommand_proj2_data/athal_wu_0_A.sorted.bam "Chr3:11,777,000-11,794,000" | cut -f6 | grep "N" | wc -l)
echo "Number of spliced alignments: $NUM_SPLICED_2"


# Problems 11-15
# Determine general information about the alignment process from the original BAM file.

# Problem 11. How many sequences are in the genome file?
echo "Number of sequences in the genome file: $(samtools view -H gencommand_proj2_data/athal_wu_0_A.bam | grep -i "@sq" | wc -l)"

# Problem 12. What is the length of the first sequence in the genome file?
echo "Length of the first sequence in the genome file: $(samtools view -H gencommand_proj2_data/athal_wu_0_A.bam | grep -i "@sq" | head -n 1 | cut -f3 | grep -o -E '[0-9]+'))"

# Problem 13. What alignment tool was used?
echo "Alignment tool used: $(samtools view -H gencommand_proj2_data/athal_wu_0_A.bam | grep -i "@PG" | head -n 1 | cut -f2)"

# Problem 14. What is the read identifier (name) for the first alignment?
echo "Name of the first alignemnt: $(samtools view -H gencommand_proj2_data/athal_wu_0_A.bam | head -n 1 | cut -f1)"

# Problem 15. What is the start position of this read’s mate on the genome? Give this as ‘chrom:pos’ if the read was mapped, or ‘*” if unmapped.
samtools view gencommand_proj2_data/athal_wu_0_A.bam | head -n 1


# Problems 16-20
# Using BEDtools, examine how many of the alignments at point 2 in the specified range overlap exons at the locus of interest.
# Use the BEDtools ‘-wo’ option to only report non-zero overlaps. The list of exons is given in the included ‘athal_wu_0_A_annot.gtf’ GTF file.

# Problem 16. How many overlaps (each overlap is reported on one line) are reported?
echo "Number of overlaps reported: $(bedtools intersect -wo -a gencommand_proj2_data/athal_wu_0_A_annot.gtf -b gencommand_proj2_data/athal_wu_0_A.bam | wc -l)"

# Problem 17. How many of these are 10 bases or longer?
echo "Number that are 10 or more bases long: $(bedtools intersect -wo -a gencommand_proj2_data/athal_wu_0_A_annot.gtf -b gencommand_proj2_data/athal_wu_0_A.bam | cut -f16 | awk '$1 > 9' | wc -l)"

# Problem 18. How many alignments overlap the annotations?
echo "Number of alignments that overlap annotations: $(bedtools intersect -wo -a gencommand_proj2_data/athal_wu_0_A_annot.gtf -b gencommand_proj2_data/athal_wu_0_A.bam | wc -l)"

# Problem 19. Conversely, how many exons have reads mapped to them?
echo "Number of exons that have reads mapped to them: $(bedtools intersect -wo -a gencommand_proj2_data/athal_wu_0_A_annot.gtf -b gencommand_proj2_data/athal_wu_0_A.bam  |  cut -f4 | sort -u | wc -l)"

# Problem 20. Question 20 If you were to convert the transcript annotations in the file “athal_wu_0_A_annot.gtf” into BED format, how many BED records would be generated?
echo "Number of BED records: $(bedtools intersect -wo -a gencommand_proj2_data/athal_wu_0_A_annot.gtf -b gencommand_proj2_data/athal_wu_0_A.bam | cut -f9 | cut -d " " -f4 | sort -u | wc -l)"
