# Special modules for running blast
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import sys 
import os
from project_model import *

init(new = True)

# Messages
blast_error = "ERROR: Failure to blast! Check your parameters again!"

# Data Structs
# Store the query proteins
query_proteins = dict()
# Store the reference proteins 
ref_proteins = dict()

# Functions
# Function to get the parameters for blast search from input file 
def get_args():
    input_file_error = "ERROR: The input file is invalid! Try again."
    paramList = []
    # open input file
    input_handle = open("input.txt", 'r')
    # seperate values from headers
    for line in input_handle:
        linesplit = line.replace(" ", "").strip("\n").split("=")
        try: 
            paramList.append(linesplit[1])
        except IndexError:
            print(input_file_error)
            sys.exit(1)
    # terminate if missing or extra parameters
    if(len(paramList) != 3):
        print(input_file_error)
        sys.exit(1)
    # get arguments from list 
    blast_command = paramList[0]
    blast_query = paramList[1]
    blast_db_path = paramList[2]
    # close file 
    input_handle.close()
    return blast_command, blast_query, blast_db_path


# Get the blast arguments 
blast_command, blast_query, blast_db_path = get_args()

# Build the command-line for forward blast
cmdline = NcbiblastpCommandline(cmd=blast_command,
                                query=blast_query,
                                db=blast_db_path,
                                outfmt=5,
                                out="blast_results.xml")
# ...and execute.
try:
    stdout, stderr = cmdline()
except:
    print(blast_error)
    sys.exit(1)

# Build the cmdline for reverse blast
cmdlinereverse = NcbiblastpCommandline(cmd = blast_command, 
                                       query = blast_db_path,
                                       db = blast_query,
                                       outfmt = 5,
                                       out="blast_results_reverse.xml"
                                       )
# and continue
try:
    stdout, stderr = cmdlinereverse()
except:
    print("Failure to blast! Check your parameters again!")
    sys.exit(1)


# Open forward blast file 
result_handle = open('blast_results.xml')
for blast_result in NCBIXML.parse(result_handle):
    # Get the query protein handle
    query_term = blast_result.query.split('|')
    # Extract the query protein info
    try:
        protein_id = query_term[4].strip(']').split('[')
    except IndexError: 
        continue
    # Update proteins table 
    p1 = Proteins(accession = query_term[3], protein_name = protein_id[0], organism = protein_id[1])

    # for each alignment 
    for align in blast_result.alignments:
        # Get the reference protein handle
        ref_term = align.title.split('|')
        # Get the reference protein idenitfier 
        try:
            ref_id = ref_term[6].strip(']').split('[')
        except IndexError:
            continue
        # if protein exists, the protein is pulled from the table. if not it is added.
        try:
            p2 = Proteins(accession = ref_term[5], protein_name = ref_id[0], organism = ref_id[1])
        except dberrors.DuplicateEntryError:
            p2 = Proteins.byAccession(ref_term[5])
        # Get the best high scoring pair
        best_hsp = align.hsps[0]
        # Add alignment to Alignments table 
        a = Alignments(query_protein = p1, ref_protein = p2,
                       alignment_length = align.length, alignment_score = best_hsp.score,
                       eval = best_hsp.expect, identity = best_hsp.identities)
# Close file handle 
result_handle.close()
        
# Open reverse blast file
result_handle = open("blast_results_reverse.xml")
for blast_result in NCBIXML.parse(result_handle):
    # Get the query protein handle
    query_term = blast_result.query.split('|')
    # Get the query protein identifier 
    try:
        protein_id = query_term[4].strip(']').split('[')
    except IndexError:
        continue
    # if protein exists, the protein is pulled from the table. if not it is added.
    try:
        p1 = Proteins(accession = query_term[3], protein_name = protein_id[0], organism = protein_id[1])
    except dberrors.DuplicateEntryError:
        p1 = Proteins.byAccession(query_term[3])
    for align in blast_result.alignments:
        # Get the reference protein handle
        ref_term = align.title.split('|')
        # Get the reference protein identifier 
        try:
            ref_id = ref_term[6].strip(']').split('[')
        except IndexError:
            continue
        # if protein exists, the protein is pulled from the table. if not it is added.
        try:
            p2 = Proteins(accession = ref_term[5], protein_name = ref_id[0], organism = ref_id[1])
        except dberrors.DuplicateEntryError:
            p2 = Proteins.byAccession(ref_term[5])
        # Get best high scoring pair
        best_hsp = align.hsps[0]
        # Add alignment a to alignments table
        a = Alignments(query_protein = p1, ref_protein = p2,
                       alignment_length = align.length, alignment_score = best_hsp.score,
                       eval = best_hsp.expect, identity = best_hsp.identities)
result_handle.close()










