from Bio import Entrez
from BaseEntrezClient import BaseEntrezClient
import csv
from io import StringIO

class SRAClient(BaseEntrezClient):

    def fetch_runinfo(self, id_list=None, verbose=False):
        """Fetch SRA runinfo (CSV metadata)."""
        if not id_list:
            raise ValueError("fetch_runinfo() requires an ID or list of IDs.")

        # Convert single ID into a list
        if isinstance(id_list, str):
            id_list = [id_list]

        # EFetch runinfo
        with Entrez.efetch(db="sra", id=",".join(id_list), rettype="runinfo", retmode="text") as handle:
            raw = handle.read()

        # Handle bytes or str safely
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")

        reader = csv.DictReader(StringIO(raw))
        records = list(reader)


        if verbose:
            print(f"RunInfo records fetched: {len(records)}")
            if len(records) > 0:
                print("Example metadata:")
                print(records[0])
                print()

        return records

# NEXT:
# def download_fastqs()
# SEND desired IDS and metadata to another file that DOWNLOADS RAW FASTQS
# USE sra-toolkit and faster-dump with python subprocess 