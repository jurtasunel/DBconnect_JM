from EntrezClient import EntrezClient

def main():
    # Initialize the client
    client = EntrezClient()

    # Show available databases
    databases = client.show_databases()
    print("Available databases:")
    print(databases.keys())
    print()

    # Search the default database and term from config
    id_list = client.search_database()
    print(f"Found {len(id_list)} sequence IDs")

    # Fetch data for these IDs
    if id_list:
        records = client.fetch_data(id_list=id_list)
        print(f"Fetched {len(records)} records\n")

        # Access individual genomes: print first 10 nucleotides of each
        for record in records:
            print(record.id, record.seq[:10], len(record.seq), record.description)

if __name__ == "__main__":
    main()
