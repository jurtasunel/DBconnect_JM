

#     # Do custom search modifiying some config attributes.
#     new_id_list = client.search_database(db = "sra", term = "SRR390728", verbose = True)
#     new_records = client.fetch_data(id_list = new_id_list, db = "sra", verbose = True, rettype = "runinfo", retmode = "text")
    

from EntrezClient import EntrezClient
from SRAClient import SRAClient

def main():
    # Initialize connection clients.
    ncbi_sequences = EntrezClient() # for complete assembled sequences.
    ncbi_sra = SRAClient() # for raw sequencing data.

    # Show available databases.
    print("List of entrez databases:\n", ncbi_sequences.show_databases(), "\n")

    # Get the IDs from the default database search from the config file.
    id_list = ncbi_sequences.search_database(verbose = True)
    # Get full records of each ID.
    records = ncbi_sequences.fetch_data(id_list = id_list, verbose = True)

    # Build dictionary to access individual record by their ID.
    records_by_id = {id: record for id, record in zip(id_list, records)}

    # Explore a single record.
    id_of_choice = id_list[0]
    print(f"Attributes of NCBI record {id_of_choice}:")
    ncbi_sequences.explore_record(records_by_id[id_of_choice])

    # Explore SRA records.
    print("\nSearching SRA...")
    sra_ids = ncbi_sra.search_database(db = "sra", term = "SRR390728", verbose = True)
    runinfo = ncbi_sra.fetch_runinfo(sra_ids, verbose=True)
    
    
    print("all ran ok")
    
if __name__ == "__main__":
    main()
