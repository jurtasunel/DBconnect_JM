from Bio import Entrez
import config

class BaseEntrezClient:
    def __init__(self):
        Entrez.email = config.email
        Entrez.api_key = config.api_key

    def show_databases(self):
        """Get info of all databases."""
        with Entrez.einfo() as handle:
            return Entrez.read(handle)

    def search_database(self, db=None, term=None, retmax=None, verbose=False):
        """Generic Entrez search."""
        db = db or config.database
        term = term or config.search_term
        retmax = retmax or config.maximum_returns

        with Entrez.esearch(db=db, term=term, retmax=retmax) as handle:
            record = Entrez.read(handle)

        if verbose:
            print(f"All accessible items of the search '{term}' in the database '{db}':")
            for k in record.keys():
                print(k)
            print(f"Total entries: {record['Count']}")
            print(f"List of IDs: {record['IdList']}\n")

        return record["IdList"]
