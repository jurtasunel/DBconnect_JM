from Bio import Entrez
import config

class BaseEntrezClient:
    def __init__(self):
        Entrez.email = config.email
        Entrez.api_key = config.api_key

    def show_databases(self):
        """Get info of all databases.
        Each connection is called a 'handle' and the information is saved as a 'record'."""
        with Entrez.einfo() as handle:
            return Entrez.read(handle)

    def search_database(self, db=None, term=None, retmax=None, verbose=False):
        """Get ID list of any database search (config file by default).
        Retmax gives the desired maximum number of records, the default is 20."""
        db = db or config.database
        term = term or config.search_term
        retmax = retmax or config.maximum_returns

        with Entrez.esearch(db=db, term=term, retmax=retmax) as handle:
            record = Entrez.read(handle)

        # Return summary of accessible options with verbose.
        if verbose:
            print(f"All accessible items of the search '{term}' in the database '{db}':")
            for key in record.keys():
                print(key)
            print(f"Total entries: {record['Count']}")
            print(f"List of IDs: {record['IdList']}\n")

        return record["IdList"]
