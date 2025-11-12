### Documentation:
# Entrez: https://biopython.org/docs/1.76/api/Bio.Entrez.html
# Entrez tutorial: https://www.youtube.com/watch?v=tl4xqdfIBh0&t=55s
# SeqIO docs: https://biopython.org/wiki/SeqIO

from Bio import Entrez, SeqIO
import config

class EntrezClient:
    def __init__(self):
        Entrez.email = config.email
        Entrez.api_key = config.api_key

    def show_databases(self):
        """Return all available databases."""
        with Entrez.einfo() as handle:
            return Entrez.read(handle)

    def search_database(self, db=None, term=None, retmax=None, verbose=False):
        """Search database and return ID list."""
        db = db or config.database
        term = term or config.search_term
        retmax = retmax or config.maximum_returns

        with Entrez.esearch(db=db, term=term, retmax=retmax) as handle:
            record = Entrez.read(handle)

        if verbose:
            print("All accessible items of the search:")
            for k in record.keys():
                print(k)
            print()
            print(f"Found {record['Count']} entries.")
            print(f"IDs: {record['IdList']}\n")

        return record["IdList"]

    def fetch_data(self, db=None, id_list=None, rettype=None, retmode=None, verbose=False):
        """Fetch the full record for given IDs."""
        if not id_list:
            raise ValueError("No IDs provided. Run search_database() or supply manually.")

        db = db or config.database
        rettype = rettype or config.return_type
        retmode = retmode or config.return_mode

        with Entrez.efetch(db=db, id=id_list, rettype=rettype, retmode=retmode) as handle:
            records = list(SeqIO.parse(handle, rettype))

        if verbose:
            print(f"Fetched {len(records)} records.\n")
            for record in records:
                print(record.name, len(record.seq), record.description)
            print()

        return records
