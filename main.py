# from EntrezClient import EntrezClient

# def main():
#     # Initialize the client.
#     client = EntrezClient()

#     # Show available databases.
#     databases = client.show_databases()
#     print(f"List of entrez databases: {databases}\n")

#     # Get the IDs from the default database search from the config file.
#     id_list = client.search_database(verbose = True)

#     # Fetch data from the IDs.
#     if id_list:
#         records = client.fetch_data(id_list = id_list, verbose = True)

#     # Build dictionary to access individual record by ID.
#     records_by_id = {unique_id : record for unique_id, record in zip(id_list, records)}

#     # Example: inspect a record by ID.
#     id_of_choice = id_list[0]
#     print(f"Attributes of record {id_of_choice}:")
#     client.explore_record(records_by_id[id_of_choice])

#     # Do custom search modifiying some config attributes.
#     new_id_list = client.search_database(db = "sra", term = "SRR390728", verbose = True)
#     new_records = client.fetch_data(id_list = new_id_list, db = "sra", verbose = True, rettype = "runinfo", retmode = "text")
    
# if __name__ == "__main__":
#     main()

from EntrezClient import EntrezClient
from SRAClient import SRAClient

def main():
    # Initialize clients
    ncbi_sequences = EntrezClient()
    ncbi_sra = SRAClient()

    # Show DB list
    print("Databases:\n", ncbi_sequences.show_databases(), "\n")

    # Default search (nucleotide)
    id_list = ncbi_sequences.search_database(verbose=True)
    records = ncbi_sequences.fetch_data(id_list=id_list, verbose=True)

    # Build UID → record dictionary
    records_by_id = {uid: rec for uid, rec in zip(id_list, records)}

    # Explore a single record
    first_id = id_list[0]
    print(f"Attributes of record {first_id}:")
    ncbi_sequences.explore_record(records_by_id[first_id])

    # SRA example
    print("\nSearching SRA...")
    sra_ids = ncbi_sra.search_database(db="sra", term="SRR390728", verbose=True)
    runinfo = ncbi_sra.fetch_runinfo(sra_ids, verbose=True)
    
    
    print("all ran ok")
    
if __name__ == "__main__":
    main()
