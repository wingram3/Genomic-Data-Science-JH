from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

sequence = "TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG"

result_handle = NCBIWWW.qblast('blastn', 'nt', sequence)
blast_record = NCBIXML.read(result_handle)

for alignment in blast_record.alignments[:5]:  # Top 5 results
    for hsp in alignment.hsps:
        print(f"****Alignment****")
        print(f"sequence: {alignment.title}")
        print(f"length: {alignment.length}")
        print(f"e value: {hsp.expect}")
        print(f"score: {hsp.score}")
        print(hsp.query)
        print(hsp.match)
        print(hsp.sbjct)

