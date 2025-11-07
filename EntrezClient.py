### Documentation:

# Entrez: https://biopython.org/docs/1.76/api/Bio.Entrez.html
# Entrez youtube tutorial: https://www.youtube.com/watch?v=tl4xqdfIBh0&t=55s
# SeqIO documentation: https://biopython.org/wiki/SeqIO

# Import modules.
from Bio import Entrez, SeqIO # NCBI access.
import config # import configuration file

class EntrezClient:

    def __init__(self):
        Entrez.email = config.email
        Entrez.api = config.api_key

    def show_databases(self):
        """ Return all available databases with Entrez. """
        with Entrez.einfo() as handle:
            return Entrez.read(handle)
        
    def search_database(self, db = None, term = None, retmax = None):
        """ Search the database of choice witht the desired search term.
            If nothing is specified in main.py, it will use the defaults from the config file """
        db = db or config.database
        term = term or config.search_term
        retmax = retmax or config.maximum_returns
        
        with Entrez.esearch(db = db, term = term, retmax = retmax) as handle:
            record = Entrez.read(handle)


        # Explore the record variable, it is a dictionary.
        print("All accessible items of the search:")
        for i in record.keys():
            print(i)
        print("\n")
        print(f"Get how many IDs are found with the query: {record['Count']}") # same as len(record["IdList"]). This gets ALL appearances of the search. 
        print(f"Print the IDs of the query: {record['IdList']}") # this gives a list of IDs, default is 20 and can be changed with retmax.
        print("\n")
    
        return record["Idlist"]
    
    def fetch_data(self, db = None, id = None, rettype = None, retmode = None):
        """ Fetch the full record for the given IDs. """

        if not id:
            raise ValueError("No IDs have been provided. Run a search with search_database() or provide them manually")
        
        db = db or config.database
        rettype = rettype or config
        retmode = retmode or config

        with Entrez.efetch(db = db, id = id, rettype = rettype, retmode = retmode) as handle:
            records = list(SeqIO.parse(handle, ))
    






    # def fetch_records(self, db, id_list, rettype="gb", retmode="text"):
    #     """Fetch full records for given IDs."""
    #     with Entrez.efetch(db=db, id=id_list, rettype=rettype, retmode=retmode) as handle:
    #         records = list(SeqIO.parse(handle, "gb"))
    #     return records