from .globals import GLOBALS
import numpy as np
import psutil
import time


class Sequence:
    def __init__(self, fasta: str = None, sequence: str = None, name: str = None, seq_type: str = 'DNA'):
        if fasta:
            fasta = fasta.split('\n', 1) 
            self.sequence = fasta[1].replace('\n', '').upper()
            self.name = fasta[0].split(' ')[0].replace('>', '')
            self.length = len(self.sequence)
            self.type = seq_type
        else:
            if sequence is None:
                sequence = ''
            sequence = sequence.replace('\n', '')
            self.sequence = sequence.upper()
            self.name = name
            self.length = len(sequence)
            self.type = seq_type
            
        
    def __str__(self):
        sequence = self.sequence
        if self.length > GLOBALS.CHAR_MAX_PRINT:
            sequence = sequence[:GLOBALS.CHAR_MAX_PRINT] + f'... omitting {self.length - GLOBALS.CHAR_MAX_PRINT} characters'
        return 'ID: ' + self.name + '; Length: ' + str(self.length) + '\n' + sequence
    
    def __getitem__(self, location: slice):            
        return Sequence(sequence = self.sequence[location], name = str(location) + ' from ' + self.name)
    
    def __len__(self):
        return self.length
    
    def __add__(self, second_sequence):
        return Sequence(sequence = self.sequence + second_sequence.sequence, name = self.name)
        
    def _matching_matrix(self, pattern):
        # Create distance matrix
        matrix = np.zeros((len(pattern) + 1, self.length + 1), dtype = int)
        
        # Initialize first column
        for x in range(len(pattern) + 1):
            matrix[x, 0] = x

        # Initialize first row
        for y in range(1,len(pattern)+1):
            matrix[0, y] = 0
            
        # Fill rest of the matrix
        for y in range(1, len(pattern) + 1):
            for x in range(1, len(self.sequence) + 1):
                horizontal_distance = matrix[y, x-1] + 1
                vertical_distance = matrix[y-1, x] + 1
                if self.sequence[x-1] == pattern[y-1]:
                    diagonal_distance = matrix[y-1, x-1]
                else:
                    diagonal_distance = matrix[y-1, x-1] + 1
                matrix[y, x] = min(horizontal_distance, vertical_distance, diagonal_distance)
                    
        return matrix
    
    def reverse(self):
        return Sequence(sequence = self.sequence[::-1], name = self.name + ' reverse', seq_type = self.type)
    
    def find(self, pattern):
        if type(pattern) == str:
            return self.sequence.find(pattern)
        elif type(pattern) == Sequence:
            return self.sequence.find(pattern.sequence)
        
    def exact_match(self, pattern):
        """Returns the indices of all exact matches"""
        matches = []
        sequence = self.sequence
        while True:
            match = sequence.find(pattern)
            if match >= 0:
                matches.append(match)
                sequence = sequence[:match] + sequence[match + len(pattern):]
                print(match)
            else:
                break
        
        return matches
    
    def approximate_match(self, pattern, max_wait_for_memory: int = float('inf')):
        """
        Queries this sequence with another pattern of nucleotides using an efficient approximate matching algorithm.

        Parameters
        ----------
        pattern : str
            The pattern which will be compared to this sequence.

        Returns
        -------
        A list containing the location(s) within the sequence which most closely match the pattern
        and
        The number of edits between this sequence and the pattern
        """
        est_memory_needed = len(pattern) * len(self) * 1.2
        timeout_counter = 0
        while psutil.virtual_memory()[1] < est_memory_needed:
            time.sleep(1)
            timeout_counter += 1
            if timeout_counter > max_wait_for_memory:
                raise TimeoutError('Approximate matching timed out waiting for memory to become available')
        if type(pattern) == Sequence:
            pattern = pattern.sequence
        matrix = self._matching_matrix(pattern)
        min_distance = min(matrix[-1])
        occurrences = [[idx] for idx in range(len(matrix[-1])) if matrix[-1][idx] == min_distance]
        for occurrence in occurrences:
            x = occurrence[0]
            y = len(matrix) - 1
            while y > 0:
                diag_distance = matrix[y - 1, x - 1]
                ver_distance = matrix[y - 1, x]
                hor_distance = matrix[y, x-1]
                minimum = min((diag_distance, ver_distance, hor_distance))
                if minimum == diag_distance:
                    y -= 1
                    x -= 1
                elif minimum == ver_distance:
                    y -= 1
                else:
                    x -= 1
            occurrence.insert(0, x + 1)
                    
        return occurrences, min_distance
    
    def compliment(self, reverse = True):
        compliment_dict = None
        if self.type == 'DNA':
            compliment_dict = GLOBALS.DNA_COMPLIMENTS
        elif self.type == 'RNA':
            compliment_dict = GLOBALS.RNA_COMPLIMENTS
        else:
            raise TypeError(f'Reverse compliments can only be generated from RNA or DNA sequences, not {self.type} sequences')
        
        if reverse:
            sequence = self.sequence[::-1]
        else:
            sequence = self.sequence
            
        new_seq = ''
        for nucleotide in sequence:
            new_seq += compliment_dict[nucleotide]
        
        return Sequence(sequence = new_seq, name = self.name + ' reverse_compliment', seq_type = self.type)
    
    def write_fasta(self, filename: str = None, return_str: bool = False):
        fasta = '>' + self.name + ' size=' + str(self.length)
        
        for x in range(int(self.length / GLOBALS.FASTA_OUTPUT_LINELENGTH) + 1):
            fasta += '\n'
            fasta += self.sequence[x * GLOBALS.FASTA_OUTPUT_LINELENGTH: (x * GLOBALS.FASTA_OUTPUT_LINELENGTH) + GLOBALS.FASTA_OUTPUT_LINELENGTH]
        
        if return_str:
            return fasta
        open(filename, 'w').write(fasta)
        
    def to_rna(self):
        if not self.type == 'DNA':
            raise TypeError(f"Tried to convert a sequence of type {self.type} to RNA. Only DNA can be converted to RNA")
        return Sequence(sequence = self.sequence.replace('T', 'U'), name = self.name + ' RNA', seq_type = 'RNA')
    
    def to_dna(self):
        if not self.type == 'RNA':
            raise TypeError(f"Tried to convert a sequence of type {self.type} to DNA. Only RNA can be converted to DNA")
        return Sequence(sequence = self.sequence.replace('U', 'T'), name = self.name, seq_type = 'DNA')
        
    def to_protein(self, offset: int = 0):
        """
        Converts the nucleotide sequence to a protein sequence

        Parameters
        ----------
        offset : int, optional
            Which position in the sequence to start from. Warning: this will affect the reading frame.
            The default is 0.

        Returns
        -------
        A Sequence.

        """
        sequence = self
        if self.type == 'RNA':
            sequence = self.to_dna()
        protein_seq = ''
        for codon_start in range(offset, len(sequence.sequence), 3):
            try:
                codon = sequence.sequence[codon_start:codon_start + 3]
                protein_seq += GLOBALS.PROTEIN_CODES[codon]
            except KeyError:
                break
        
        return Sequence(sequence = protein_seq, name = sequence.name + ' protein', seq_type = 'protein')
    
    def gene_to_protein(self, gff_file):
        gff_lines = open(gff_file, 'r').readlines()
        gene_line = [x for x in gff_lines if '\tgene\t' in x][0]
        gene_start = int(gene_line.split('\t')[3])
        exons = [x for x in gff_lines if '\texon\t' in x]
        exon_indices = []
        sequence = ''
        for exon in exons:
            fields = exon.split('\t')
            sense = fields[6]
            start = int(fields[3]) - gene_start
            end = int(fields[4]) - gene_start
            offset = fields[7]
            if offset == '.':
                offset = 0
            else:
                offset = int(offset)
            exon_indices.append((start, end, sense))
            if sense == '+':
                sequence += self[start, end + 1].to_protein(offset = offset).sequence
            else:
                print(self[start : end + 1].reverse().to_protein(offset = offset).sequence)
                sequence += self[start : end + 1].reverse().to_protein(offset = offset).sequence
        return Sequence(sequence = sequence, name = self.name + ' protein', seq_type='protein')
