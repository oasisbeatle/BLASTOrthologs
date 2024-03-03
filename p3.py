from project_model import * 

# Initialize project model 
init() 
# Set that has all unique proteins 
protein_set = set()
# Set of proteins explored by ortholog finder
explored_proteins = set()
# List of alignments in for a blast query 
alignments_list = []
# List of orthologs 
ortholog_list = []
        
# Get all proteins from the protein table
for p in Proteins.select():
    protein_set.add(p)

# for each protein in the protein set 
for p in protein_set:
    # if the protein is already explored move to the next protein 
    if p in explored_proteins:
        continue
    
    # Add current protein to explored list 
    explored_proteins.add(p)
    
    # Clear alignments for previous protein 
    alignments_list.clear()
    # For all alignments with protein p as a query sequence
    for align in Alignments.select(Alignments.q.query_protein == p):
        # Add all alignments to the list 
        alignments_list.append(align)
        # Sort by e-value 
        sorted_list = sorted(alignments_list, key = lambda x: (x.eval, -(x.alignment_length)))
        # Return the best if exists 
        try:
            best_hit = sorted_list[0]
        except IndexError:
            continue   

   # Get the best hit protein 
    p2 = best_hit.ref_protein

    # Clear alignments list 
    alignments_list.clear() 
    # Search for all the best hit alignments
    for align in Alignments.select(Alignments.q.query_protein == p2):
         alignments_list.append(align)
    # Sort the alignments of the best hit
    sorted_list.clear()
    sorted_list = sorted(alignments_list, key = lambda x: (x.eval, -(x.alignment_length)))

    # Get the reverse best hit protein 
    try:
        best_reverse_hit = sorted_list[0]
    except IndexError:
        continue
    reverse_protein = best_reverse_hit.ref_protein

    # If there is a mutual best hit then we add to ortholog list
    if(p.accession == reverse_protein.accession):
        ortholog_list.append((p, p2))

# Print orthologs
for ortholog in ortholog_list:
    print("The protein", ortholog[0].protein_name, "with acession", ortholog[0].accession, "in organism", ortholog[0].organism,
          "is orthologous to protein", ortholog[1].protein_name, "with acession", ortholog[1].accession, "in organism", 
          ortholog[1].organism)