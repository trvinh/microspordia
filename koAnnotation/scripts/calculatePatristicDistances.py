import os, sys
import dendropy
from dendropy import treecalc

# Calculate the patristic distances for the gene tree and store it in a temporary file
# 	INPUT	1) Tree
#		2) Output file
#

gene_tree = os.path.abspath(sys.argv[1])
output = os.path.abspath(sys.argv[2])

fnew = open(output, 'w')
tree = dendropy.Tree.get(path=gene_tree, schema="newick")
pdm = tree.phylogenetic_distance_matrix()
for idx1, taxon1 in enumerate(tree.taxon_namespace):
	for taxon2 in tree.taxon_namespace[idx1+1:]:
		# Weighted == Sum of branch lengths
		weighted_patristic_distance = pdm.patristic_distance(taxon1, taxon2)
		# Unweighted == Edge counts between 2 species
        	unweighted_patristic_distance = pdm.path_edge_count(taxon1, taxon2)
   		fnew.write(taxon1.label.replace(' ', '_') + '\t' + taxon2.label.replace(' ', '_') + '\t' + str(weighted_patristic_distance) + '\t' + str(unweighted_patristic_distance) + '\n')

fnew.close()
