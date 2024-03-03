import sys 
from project_model import * 

# Initialize database 
init()

if len(sys.argv) < 3:
    print("Please enter two proteins to be matched!", file=sys.stderr)
    sys.exit(1)

# Enter accesssion number 1
protein_1 = sys.argv[1]
# Enter acession number 2 
protein_2 = sys.argv[2]

# Lookup proteins by name 
for protein in Proteins.select(Proteins.q.accession == protein_1):
    protein_1_key = protein

for protein in Proteins.select(Proteins.q.accession == protein_2):
    protein_2_key = protein

# Boolean to check if alignment_exists 
align_bool = 0 
# Fetch forward alignment
try:
    align_1 =  Alignments.select(Alignments.q.query_protein == protein_1_key
                                & Alignments.q.ref_protein == protein_2_key)[0]
except NameError:
    print("The protein accessions can not be found in the database!")
    sys.exit(1)
except IndexError:
    align_bool = 1

# Fetch reverse alignment 
try:
    align_2 = Alignments.select(Alignments.q.query_protein == protein_2_key
                                & Alignments.q.ref_protein == protein_1_key)[0]
except NameError:
    print("The protein accessions can not be found in the database!")
    sys.exit(1)
except IndexError:
    align_bool = 1

# Print Alignments
if(align_bool == 1):
    print("No alignment found between proteins!")
else:
    print("The alignment between protein 1", protein_1, "and protein 2", protein_2, ":")
    print("Alignment Length:", align_1.alignment_length)
    print("Alignment Score:", align_1.alignment_score)
    print("Eval:", align_1.eval)
    print("Identity:", align_1.identity)

    print("")

    print("The alignment between protein 2", protein_2, "and protein 1", protein_1, ":")
    print("Alignment Length:", align_2.alignment_length)
    print("Alignment Score:", align_2.alignment_score)
    print("Eval:", align_2.eval)
    print("Identity:", align_2.identity)


