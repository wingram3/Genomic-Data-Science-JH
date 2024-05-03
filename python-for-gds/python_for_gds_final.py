from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class FastaParser:
    """
    object to perform various operations on dna sequences from a FASTA file

    Attributes:
        records: a list to contain Bio.SeqRecord iterators containing info on dna sequences from 
        the input FASTA file 
        """

    def __init__(self):
        self.records = []


    def parse_fasta(self, filename: str) -> List[SeqRecord]:
        """
        Parses a fasta file

        Args:
            filename: the input fasta file with dna sequences to be parsed

        Returns:
            A list containing Bio.SeqRecord iterator, one for each dna sequence in the FASTA file
        """
        with open(filename) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                 self.records.append(record)
        return self.records


    def count_num_records(self) -> int:
        """
        Counts the number of records in list of dna records.

        Returns:
            An int equal to the number of records in the records list   
        """
        return len(self.records)


    def get_seq_lengths(self) -> Dict[str, int]:
        """
        Gets the length of each dna sequence from a dict of records, returns a dict with lengths

        Returns:
            A dict mapping dna sequence names to their respective lengths, 
            or how many nucelotide bases there are in the sequence
        """
        return {record.id: len(record.seq) for record in self.records}



    def get_longest_sequences(self) -> Dict[str, int]:
        """
        Gets the longest sequence/sequences from a dict of dna sequences lengths

        Returns:
            A dict mapping the name of the longest dna sequence in the input FASTA file to its 
            respective length. If there are multiple sequences of maximum length, the dict will
            include both of them.
        """
        lengths = self.get_seq_lengths()
        max_length = max(lengths.values())
        return {record.id: len(record.seq) for record in self.records if len(record.seq) == max_length}


    def get_shortest_sequences(self) -> Dict[str, int]:
        """
        Gets the shortest sequence/sequences from a dict of dna sequences lengths

        Returns:
            A dict mapping the ID of the shortest dna sequence in the input FASTA file to its 
            respective length. If there are multiple sequences of minimum length, the dict will
            include both of them.
        """
        lengths = self.get_seq_lengths()
        min_length = min(lengths.values())
        return {record.id: len(record.seq) for record in self.records if len(record.seq) == min_length}


    def _get_forward_reading_frames(self) -> Dict[str, Dict[str, List[str]]]:
        """
        Gets each of the three forward reading frames for all dna sequences in the input FASTA file

        Returns:
            A dict of dictsmapping record IDs to a list of their three respective forward 
            reading frames
        """
        forward_reading_frames = {}
        for record in self.records:
            forward_reading_frames[record.id] = {}
            frame1 = [str(record.seq[i : i+3]) for i in range(0, len(record.seq), 3)]
            frame2 = [str(record.seq[i : i+3]) for i in range(1, len(record.seq), 3)]
            frame3 = [str(record.seq[i : i+3]) for i in range(2, len(record.seq), 3)]
            forward_reading_frames[record.id][1] = frame1
            forward_reading_frames[record.id][2] = frame2
            forward_reading_frames[record.id][3] = frame3
        return forward_reading_frames

