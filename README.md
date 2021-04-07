# Usage
Scan a phylogenetic tree and find subtrees with strong temporal siganls

# Input
A rooted phylogenetic tree

# Output
Subtrees with strong temporal signals

# Parameters
"size" indicates the minimum number of leaves within a internal node;

"threshold" indicates the minimum squared coefficient (R2) of either the Spearman's or the Pearson's correlation;

"sources" indicates if wants to calculate simpson index of sources within a internal node, default is "none";

"simpson_threhold" indicates the minimum value of simpson index.

# Example
input: SE_SNP_tree_msa_phyml.tree

output: SE_SNP_tree_msa_phyml_time_signals.pdf

