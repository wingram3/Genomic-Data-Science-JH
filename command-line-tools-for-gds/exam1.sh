# Problem 1. How many chromosomes are there in the genome?
NUM_CHROMOSOMES=$(grep ">" gencommand_proj1_data/apple.genome | wc -l)
echo "Number of chromosomes: $NUM_CHROMOSOMES"

# Problem 2. How many genes?
NUM_GENES=$(cut -f1 gencommand_proj1_data/apple.genes | uniq | wc -l)
echo "Number of genes: $NUM_GENES"

# Problem 3. How many transcript variants?
NUM_TV=$(cut -f2 gencommand_proj1_data/apple.genes | uniq | wc -l)
echo "Number of transcript variants: $NUM_TV"

# Problem 4. How many genes have a single splice variant?
NUM_SSV=$(cut -f1 gencommand_proj1_data/apple.genes | uniq -c | grep " 1 " | wc -l)
echo "Number of genes with a single splice variant: $NUM_SSV"

# Problem 5. How many genes have two or more splice variants?
NUM_MSV=$(cut -f1 gencommand_proj1_data/apple.genes | uniq -c | grep -v " 1 " | wc -l)
echo "Number of genes with two or more splice variants: $NUM_MSV"

# Problem 6. How many genes are there on the ~@~X+~@~Y strand?
NUM_GPS=$(sort -k 4 -k 1 gencommand_proj1_data/apple.genes | grep "+" | wc -l)
echo "Number of genes on the + strand: $NUM_GPS"

# Problem 7. How many genes are on the '-' strand?
NUM_GMS=$(sort -k 4 gencommand_proj1_data/apple.genes | grep -v "+" | wc -l)
echo "Number of genes on the - strand: $NUM_GMS"

# Problem 8. How many genes are there on chromosome chr1?
NUM_CHR1=$(sort -k 4 -k 1 gencommand_proj1_data/apple.genes | grep "chr1" | wc -l)
echo "Number of genes on chromosome 1: $NUM_CHR1"

# Problem 9. How many genes are there on chromosome chr2?
NUM_CHR2=$(sort -k 4 -k 1 gencommand_proj1_data/apple.genes | grep "chr2" | wc -l)
echo "Number of genes on chromosome 2: $NUM_CHR2"

# Problem 10. How many genes are there on chromosome chr3?
NUM_CHR3=$(sort -k 4 -k 1 gencommand_proj1_data/apple.genes | grep "chr3" | wc -l)
echo "Number of genes on chromosome 3: $NUM_CHR3"

# Problem 11. How many transcripts are there on chr1"
NUM_TR_CHR1=$(sort -k 3 -k 5n gencommand_proj1_data/apple.genes | grep "chr1" | wc -l)
echo "Number of transcripts on chromosome 1: $NUM_TR_CHR1"

# Problem 12. How many transcripts are there on chr1"
NUM_TR_CHR2=$(sort -k 3 -k 5n gencommand_proj1_data/apple.genes | grep "chr2" | wc -l)
echo "Number of transcripts on chromosome 2: $NUM_TR_CHR2"

# Problem 13. How many transcripts are there on chr1"
NUM_TR_CHR3=$(sort -k 3 -k 5n gencommand_proj1_data/apple.genes | grep "chr3" | wc -l)
echo "Number of transcripts on chromosome 3: $NUM_TR_CHR3"

# Problem 14. How many genes are in common between condition A and condition B?
cut -f1 gencommand_proj1_data/apple.conditionA | sort -u > /tmp/apple.conditionA.sorted.unique
cut -f1 gencommand_proj1_data/apple.conditionB | sort -u > /tmp/apple.conditionB.sorted.unique
NUM_COMMON_AB=$(comm -1 -2 /tmp/apple.conditionA.sorted.unique /tmp/apple.conditionB.sorted.unique | wc -l)
echo "Number of genes in common between conditions A and B: $NUM_COMMON_AB"

# Problem 15. How many genes are specific to condition A?
cut -f1 gencommand_proj1_data/apple.conditionA | sort -u > /tmp/apple.conditionA.sorted.unique
cut -f1 gencommand_proj1_data/apple.conditionB | sort -u > /tmp/apple.conditionB.sorted.unique
NUM_SPEC_A=$(comm -2 -3 /tmp/apple.conditionA.sorted.unique /tmp/apple.conditionB.sorted.unique | wc -l)
echo "Number of genes specific to condition A: $NUM_SPEC_A"

# Problem 16. How many genes are specific to condition B?
cut -f1 gencommand_proj1_data/apple.conditionA | sort -u > /tmp/apple.conditionA.sorted.unique
cut -f1 gencommand_proj1_data/apple.conditionB | sort -u > /tmp/apple.conditionB.sorted.unique
NUM_SPEC_B=$(comm -1 -3 /tmp/apple.conditionA.sorted.unique /tmp/apple.conditionB.sorted.unique | wc -l)
echo "Number of genes specific to condition B: $NUM_SPEC_B"

# Problem 17. How many genes are common to all three conditions?
cut -f1 gencommand_proj1_data/apple.conditionA | sort -u > /tmp/apple.conditionA.sorted.unique
cut -f1 gencommand_proj1_data/apple.conditionB | sort -u > /tmp/apple.conditionB.sorted.unique
comm -1 -2 /tmp/apple.conditionA.sorted.unique /tmp/apple.conditionB.sorted.unique > /tmp/apple.condition_common_to_A_B.sorted.unique
cut -f1 gencommand_proj1_data/apple.conditionC | sort -u > /tmp/apple.conditionC.sorted.unique
NUM_COMM_ABC=$(comm -1 -2 /tmp/apple.condition_common_to_A_B.sorted.unique /tmp/apple.conditionC.sorted.unique | wc -l)
echo "Number of genes common to A, B, and C: $NUM_COMM_ABC"
