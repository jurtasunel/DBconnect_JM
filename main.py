
import config
from EntrezClient import EntrezClient

def main():
    client = EntrezClient()
    
    ids = client.search_database(verbose=True)
    records = client.fetch_data(id_list = ids, verbose=True)


    print(records)

if __name__ == "__main__":
    main()