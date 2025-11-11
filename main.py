<<<<<<< HEAD

import config
from EntrezClient import EntrezClient

def main():
    client = EntrezClient()
    
    ids = client.search_database(verbose=True)
    records = client.fetch_data(id_list = ids, verbose=True)


    print(records)

if __name__ == "__main__":
    main()
=======
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
        records = client.fetch_data(id=id_list)
        print(f"Fetched {len(records)} records")

        # Loop thorhg the records
        for record in records:
            print(record.name, len(record.seq), record.id, record.seq[:10], record.description)

if __name__ == "__main__":
    main()
>>>>>>> 62ebc2d (added main)
