### This contains setup variables and credentials.

# Documentation:
# Rettype and retmode for atabase: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?


# Email credentials.
email = "jmurtasun.94@gmail.com"
api_key = "1c5b74af4f785579858076ca3d607d927008"

# Search information.
database = "nucleotide"
#search_term = 'esxi[Gene name] AND "Mycobacterium bovis"[Organism]'
search_term = '"coesin[Gene name] Fragaria"[Organism]'
maximum_returns = 21

# Data fetching information.
return_type = "gb"
return_mode = "text"

# TODO:
# migrate config file to .yaml format