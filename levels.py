# Script to save the level of each GO term to a file

# Since we have already goatools installed, we just check the goatools version
from importlib.metadata import version
version('goatools')
from goatools.obo_parser import GODag

# Load the GO DAG from the OBO file
obo_file = "go-basic.obo"
godag = GODag(obo_file)

# Open a text file for writing
output_file = "go_terms_with_levels.txt"
with open(output_file, 'w') as f:
    # Write the header
    f.write("GO Term\tLevel\n")
    
    # Iterate through all GO terms in the GODag
    for go_id, go_term in godag.items():
        level = go_term.level if hasattr(go_term, 'level') else 'N/A'
        f.write(f"{go_id}\t{level}\n")

print(f"GO terms and levels written to {output_file}")
