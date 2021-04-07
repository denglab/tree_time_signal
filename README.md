# Usage
scan a phylogenetic tree and find subtrees with strong temporal siganls

# input
a rooted phylogenetic tree

# parameters
size indicates the minimum number of leaves within a internal node;

threshold indicates the minimum squared coefficient (R2) of either the Spearman's or the Pearson's correlation;

sources indicates if wants to calculate simpson index of sources within a internal node;

simpson_threhold indicates the minimum value of simpson index.

# example
input: SE_SNP_tree_msa_phyml.tree

output: SE_SNP_tree_msa_phyml_time_signals.pdf

