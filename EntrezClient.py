# ### Documentation:
# # Entrez: https://biopython.org/docs/1.76/api/Bio.Entrez.html
# # Entrez tutorial: https://www.youtube.com/watch?v=tl4xqdfIBh0&t=55s
# # SeqIO docs: https://biopython.org/wiki/SeqIO
# # SeqIO : https://biopython.org/docs/1.85/Tutorial/chapter_seqio.html


# from Bio import Entrez, SeqIO
# import config
# import csv

# class EntrezClient:
#     def __init__(self):
#         Entrez.email = config.email
#         Entrez.api_key = config.api_key

#     def show_databases(self):
#         """Get info of all databases.
#         Each connection is called a 'handle' and the information is saved as a 'record'."""
#         with Entrez.einfo() as handle:
#             return Entrez.read(handle)

#     def search_database(self, db = None, term = None, retmax = None, verbose = False):
#         """Get ID list of any database search (config file by default).
#         Retmax gives the desired maximum number of records, the default is 20."""
#         db = db or config.database
#         term = term or config.search_term
#         retmax = retmax or config.maximum_returns

#         with Entrez.esearch(db = db, term = term, retmax = retmax) as handle:
#             record = Entrez.read(handle)
            
#         # Return summary of accessible options with verbose.
#         if verbose:
#             print(f"All accessible items of the search {config.search_term}:")
#             for item in record.keys():
#                 print(item)
#             print(f"Total entries: {record['Count']}")
#             print(f"List of IDs: {record['IdList']}\n")

#         return record["IdList"]

#     def fetch_data(self, db = None, id_list = None, rettype = None, retmode = None, verbose = False):
#         """Fetch the full record from the given IDs."""
#         if not id_list:
#             raise ValueError("fetch_data() requires a list of IDs, but nothing was provided. Run search_database() to get IDs or supply them manually.")

#         db = db or config.database
#         rettype = rettype or config.return_type
#         retmode = retmode or config.return_mode

#         # Dictionary to map rettypes to accepted format by SeqIO. There's many rettypes, but only fasta and gb are accepted by SeqIO.
#         parser_dict = {"gb": "gb", "genbank": "gb", "fasta": "fasta", "fasta_cds_na": "fasta", "fasta_cds_aa": "fasta"}
#         parser_format = parser_dict.get(rettype.lower())
#         if parser_format is None:
#             allowed_format = ", ".join(parser_dict.keys())
#             raise ValueError(f"The rettype: {rettype} is not parseable by SeqIO. Use one of these: {allowed_format}")

#         with Entrez.efetch(db=db, id=id_list, rettype=rettype, retmode=retmode) as handle:
#             records = list(SeqIO.parse(handle, parser_format))

#         if verbose:
#             print(f"Total records fetched:{len(records)}")
#             # Make sure there are at least three records before printing.
#             if len(records) > 3:
#                 print("Structure of records (first three shown):")
#                 for record in records[:3]:
#                     print(f"    Name: {record.name}, ID: {record.id}, head: {record.seq[:10]}..., length: {len(record.seq)}, description: {record.description}, dbx_identifiers: {record.dbxrefs}")
#                 print("\n")
                
#         return records
    
#     def explore_record(self, record):
#         """Get the attributes from one SeqRecord using python vars().
#         vars(record) is endless because some atributes have nested list with thousands of gene names, annotations , etc,
#         therefore use vars(record).keys() or .items() to get only the top level attributes."""
        
#         for key, value in vars(record).items():
#             # Add length to attributes without length before printing.
#             try:
#                 length = len(value)
#             except TypeError:
#                 length = None
            
#             print(f"Attribute: {key}, type: {type(value)}, length: {length}")


from Bio import Entrez, SeqIO
from BaseEntrezClient import BaseEntrezClient
import config

class EntrezClient(BaseEntrezClient):

    def fetch_data(self, db=None, id_list=None, rettype=None, retmode=None, verbose=False):
        """Fetch the full record from the given IDs."""
        if not id_list:
            raise ValueError("fetch_data() requires a list of IDs, but nothing was provided.")

        db = db or config.database
        rettype = rettype or config.return_type
        retmode = retmode or config.return_mode

        # Map rettypes that SeqIO can parse
        parser_dict = {
            "gb": "gb",
            "genbank": "gb",
            "fasta": "fasta",
            "fasta_cds_na": "fasta",
            "fasta_cds_aa": "fasta"
        }

        parser_format = parser_dict.get(rettype.lower())
        if parser_format is None:
            allowed_format = ", ".join(parser_dict.keys())
            raise ValueError(
                f"The rettype: {rettype} is not parseable by SeqIO. "
                f"Use one of these: {allowed_format}"
            )

        with Entrez.efetch(db=db, id=id_list, rettype=rettype, retmode=retmode) as handle:
            records = list(SeqIO.parse(handle, parser_format))

        if verbose:
            print(f"Total records fetched:{len(records)}")
            if len(records) > 3:
                print("Structure of records (first three shown):")
                for r in records[:3]:
                    print(f"    Name: {r.name}, ID: {r.id}, head: {r.seq[:10]}..., length: {len(r.seq)}, description: {r.description}, dbx_identifiers: {r.dbxrefs}")
                print()

        return records

    def explore_record(self, record):
        """Explore top-level attributes of a SeqRecord."""
        for key, value in vars(record).items():
            try:
                length = len(value)
            except TypeError:
                length = None
            print(f"Attribute: {key}, type: {type(value)}, length: {length}")

# NEXT
# def seq_to_fasta()
#SeqIO write to fasta outputs
    