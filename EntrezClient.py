# ### Documentation:
# # Entrez: https://biopython.org/docs/1.76/api/Bio.Entrez.html
# # Entrez tutorial: https://www.youtube.com/watch?v=tl4xqdfIBh0&t=55s
# # SeqIO docs: https://biopython.org/wiki/SeqIO
# # SeqIO : https://biopython.org/docs/1.85/Tutorial/chapter_seqio.html

from BaseEntrezClient import BaseEntrezClient
from Bio import Entrez, SeqIO
import config

class EntrezClient(BaseEntrezClient):

    def fetch_data(self, db=None, id_list=None, rettype=None, retmode=None, verbose=False):
        """Fetch the full record from the given IDs."""
        if not id_list:
            raise ValueError("fetch_data() requires a list of IDs, but nothing was provided.")

        db = db or config.database
        rettype = rettype or config.return_type
        retmode = retmode or config.return_mode

        # Dictionary to map rettypes to accepted format by SeqIO. There's many rettypes, but only fasta and gb are accepted by SeqIO.
        parser_dict = {"gb":"gb", "genbank":"gb", "fasta":"fasta", "fasta_cds_na":"fasta", "fasta_cds_aa":"fasta"}

        parser_format = parser_dict.get(rettype.lower())
        if parser_format is None:
            allowed_format = ", ".join(parser_dict.keys())
            raise ValueError(f"The rettype: {rettype} is not parseable by SeqIO. Use one of these: {allowed_format}")

        with Entrez.efetch(db=db, id=id_list, rettype=rettype, retmode=retmode) as handle:
            records = list(SeqIO.parse(handle, parser_format))

        if verbose:
            print(f"Total records fetched:{len(records)}")
            # Make sure there are at least three records before printing.
            if len(records) > 3:
                print("Structure of records (first three shown):")
                for r in records[:3]:
                    print(f"    Name: {r.name}, ID: {r.id}, head: {r.seq[:10]}..., length: {len(r.seq)}, description: {r.description}, dbx_identifiers: {r.dbxrefs}")
                print()

        return records

    def explore_record(self, record):
        """Get the attributes from one SeqRecord using python vars().
        vars(record) is endless because some atributes have nested list with thousands of gene names, annotations , etc,
        therefore use vars(record).keys() or .items() to get only the top level attributes."""
        for key, value in vars(record).items():
            try:
                length = len(value)
            except TypeError:
                length = None
            print(f"Attribute: {key}, type: {type(value)}, length: {length}")

# TODO:
# def seq_to_fasta()
#SeqIO write to fasta outputs
    