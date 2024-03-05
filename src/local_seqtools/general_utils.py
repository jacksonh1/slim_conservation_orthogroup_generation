import re
from pathlib import Path

from Bio import Align, AlignIO, Seq, SeqIO
from Bio.SeqRecord import SeqRecord


def import_fasta(fasta_path, output_format='list') -> list[SeqRecord]|dict[str, SeqRecord]:
    """import fasta file into a list or dictionary of SeqRecord objects

    Parameters
    ----------
    fasta_path : str
        file path to fasta file
    output_format : str, optional
        output_format of output. Either 'list' or 'dict'. by default 'list'

    Returns
    -------
    list or dictionary
        list or dictionary of SeqRecord objects for each sequence in the fasta file
    """
    allowed_formats = ['list', 'dict']
    with open(fasta_path) as handle:
        if output_format == 'list':
            seqs = list(SeqIO.parse(handle, 'fasta'))
        elif output_format == 'dict':
            seqs = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        else:
            raise ValueError(f"Invalid output format - {output_format}. Expected one of: {allowed_formats}")
    return seqs



class FastaImporter:
    """import fasta file and return seqrecord objects in various formats

    Parameters
    ----------
    fasta_path : str
        file path to fasta file
    """
    def __init__(self, fasta_path: str|Path):
        self.fasta_path = fasta_path

    def import_as_list(self) -> list[SeqRecord]:
        """return list of SeqRecord objects for each sequence in the fasta file

        Returns
        -------
        List[SeqRecord]
            list of SeqRecord objects
        """        
        with open(self.fasta_path) as handle:
            return list(SeqIO.parse(handle, 'fasta'))

    def import_as_dict(self) -> dict[str, SeqRecord]:
        """return dictionary of SeqRecord objects for each sequence in the fasta file

        Returns
        -------
        dict[str, SeqRecord]
            dictionary of SeqRecord objects, keys are the sequence ids and values are the SeqRecord objects
        """        
        with open(self.fasta_path) as handle:
            return SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        
    def import_as_alignment(self) -> Align.MultipleSeqAlignment:
        """return multiple sequence alignment object

        Returns
        -------
        Align.MultipleSeqAlignment
            multiple sequence alignment object
        """        
        with open(self.fasta_path) as handle:
            return AlignIO.read(handle, 'fasta')


def split_uniprot(prot_id):
    j = re.compile(r"^[st][pr]\|(.+)\|(.+)")
    prot = j.findall(prot_id)[0]
    accession= prot[0]
    name=prot[1]
    return name, accession


def get_regex_matches(regex_pattern: str, seq_str: str):
    """searches for all matches of a regex pattern in a sequence string
    returns a generator object that yields the match sequence, start index, and end index

    Parameters
    ----------
    regex_pattern : str
        regular expression pattern
    seq_str : str
        string to search for matches

    Yields
    ------
    tuple
        (match sequence, start index, end index)
    """    
    p = re.compile(regex_pattern)
    for m in p.finditer(seq_str):
        if m.start() == m.end():
            # even if there are groups in the lookahead, the first group should be the full match b/c that group surrounds the entire regex
            # so this will work whether or not there are groups in the lookahead
            match_seq = m.groups()[0]
        else:
            match_seq = seq_str[m.start() : m.end()]
        yield match_seq, m.start(), m.start() + len(match_seq)-1
