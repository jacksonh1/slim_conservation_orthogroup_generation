from Bio import SeqIO


def import_fasta(fasta_path, output_format='list') -> list[SeqIO.SeqRecord]|dict[str, SeqIO.SeqRecord]:
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

