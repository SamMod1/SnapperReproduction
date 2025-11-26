from .globals import GLOBALS
from .Sequence import Sequence
import random
from multiprocessing import Pool, cpu_count


class MultiSequence:
    def __init__(self, fasta: str = None, sequences: dict = None):
        self.pattern = None
        self.length = 0
        if fasta:
            sequences = fasta.split('>')
            self.sequences = {}
            for sequence in sequences:
                if sequence:
                    sequence = Sequence(sequence)
                    self.length += len(sequence)
                    self.sequences[sequence.name] = sequence
                    
        elif sequences:
            self.sequences = sequences
            for sequence in self.sequences:
                self.length += len(self.sequences[sequence])
                
        else:
            self.sequences = {}
            
    def __len__(self):
        return len(self.sequences)
    
    def __iter__(self):
        return (self[x] for x in self.sequences.keys())
    
    def __str__(self):
        length = len(self)
        string = str(length) + " total sequences:\n"
        sequences = list(self.sequences.keys())
        if length > GLOBALS.SEQ_MAX_PRINT:
            sequences = sequences[:GLOBALS.SEQ_MAX_PRINT]
            sequences.append(f'... omitting {length - GLOBALS.SEQ_MAX_PRINT} sequence names')
        for sequence in sequences:
            string += str(sequence) + ' '
        
        return string
    
    def __getitem__(self, sequence_name: str):
        return self.sequences[sequence_name]
    
    def __setitem__(self, name, sequence: Sequence):
        sequence.name = name
        self.sequences[name] = sequence
        
    def full_len(self):
        length = 0
        for sequence in self.sequences:
            length += len(self.sequences[sequence])
        return length
        
    def remove(self, key):
        del self.sequences[key]
        
    def n_nucleotides(self):
        return self.length
    
    def to_rna(self):
        for sequence in self.sequences:
            self.sequences[sequence] = self.sequences[sequence].to_rna()
    
    def to_dna(self):
        for sequence in self.sequences:
            self.sequences[sequence] = self.sequences[sequence].to_dna()
    
    def to_protein(self):
        for sequence in self.sequences:
            self.sequences[sequence] = self.sequences[sequence].to_protein()
        
    def write_fasta(self, filename, return_str = False):
        fasta = []
        
        for sequence in self.sequences:
            fasta.append(self.sequences[sequence].write_fasta(return_str = True))
        
        fasta = '\n'.join(fasta)
        
        if not return_str:
            open(filename, 'w').write(fasta)
        else:
            return fasta
        
    def get_multiple(self, sequences):
        new_sequences = {}
        for sequence in sequences:
            new_sequences[sequence] = self.sequences[sequence]
        return MultiSequence(sequences = new_sequences)
    
    def exact_match(self, pattern):
        """Returns the indices of all exact matches in all sequences"""
        if type(pattern) == Sequence:
            pattern = pattern.sequence
        matches = [(self.sequences[x].name, self.sequences[x].exact_match(pattern)) for x in self.sequences]
        matches = []
        for sequence in self.sequences:
            these_matches = self.sequences[sequence].exact_match(pattern)
            if these_matches:
                matches.append((self.sequences[sequence].name, these_matches))
        
        return matches
        
    
    def approximate_match(self, pattern):
        if type(pattern) == Sequence:
            pattern = pattern.sequence
        matches = [[x] + self.sequences[x].approximate_match(pattern) for x in self.sequences]
        closest_match = min(matches, key = lambda x: x[1])
        closest_matches = [x for x in matches if x[1] == closest_match]
        matched_sequences = []
        for sequence_match in closest_matches:
            for match in sequence_match[1]:
                matched_sequences.append(self.sequences[sequence_match[0]][match[0]:match[1]])
                
        return matched_sequences
    
    def _assign_cpus(self):
        """
        This breaks up the list of sequences efficiently between each cpu. This is intended to
        be used in conjuncture with multiprocessing.map
        
        This function tries to distrubute sequences evenly in terms of both cpu load and memory load

        Returns
        -------
        new_keys : iterable
            DESCRIPTION.

        """
        keys = list(self.sequences.keys())
        keys = sorted(keys, key = lambda x: len(self.sequences[x]))
        new_keys = []
        cpus = cpu_count()
        jobs = [[] for _ in range(cpus)]
        
        for idx, key in enumerate(keys):
            cpu = (idx + 1) % (cpus * 2)
            if cpu > cpus:
                cpu = (cpus + 1) - (cpu - cpus)
            if not cpu:
                cpu = 1
            jobs[cpu - 1].append(key)
        
        for cpu in range(cpus):
            if not cpu % 2 == 0:
                jobs[cpu].reverse()
            new_keys.extend(jobs[cpu])
        
        return new_keys
    
    def _parallel_match(self, sequence):
        return (sequence,) + self.sequences[sequence].approximate_match(self.pattern)
    
    def parallel_approximate_match(self, pattern):
        self.pattern = pattern
        if type(pattern) == Sequence:
            pattern = pattern.sequence
        with Pool(processes = cpu_count()) as pool:
            matches = pool.map(self._parallel_match, self._assign_cpus())
        closest_match = min(matches, key = lambda x: x[2])[2]
        closest_matches = [x for x in matches if x[2] == closest_match]
        matched_sequences = MultiSequence()
        count = 1
        for sequence_match in closest_matches:
            for match in sequence_match[1]:
                matched_sequences[
                    f'{sequence_match[0]}-edit_distance:{closest_match}_number:{count}'
                    ] = self.sequences[sequence_match[0]][match[0]:match[1]]
                count += 1
        
        self.pattern = None
        
        return matched_sequences, matches
            