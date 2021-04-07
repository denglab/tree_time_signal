#-------------------------------------------------------------------------------
# python 2.7
# Name: screen_subtrees_with_time_signal
# Author: Shaoting Li, Shaokang Zhang
# Usage: scan a phylogenetic tree and find subtrees with strong temporal siganls
# input: a rooted phylogenetic tree
# parameters: size indicates the minimum number of leaves within a internal node;
#             threshold indicates the minimum squared coefficient (R2) of either the Spearman's or the Pearson's correlation;
#             sources indicates if wants to calculate simpson index of sources within a internal node;
#             simpson_threhold indicates the minimum value of simpson index.
#-------------------------------------------------------------------------------

from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace,TextFace
import os

filepath = os.getcwd()

def simpson(counts):
  import numpy as np
  #calculate simpson index based on counts, e.g. [3,2,1]
  counts=np.array(counts)
  d=float(sum(counts*(counts-1)))/((counts.sum())*(counts.sum()-1))
  return (1-d)

def Uniq(L,sort_on_fre="none"):
  #return the uniq list and the count number
  Old=L
  L.sort()
  L = [L[i] for i in range(len(L)) if L[i] not in L[:i]]
  count=[]
  for j in range(len(L)):
    y=0
    for x in Old:
      if L[j]==x:
        y+=1
    count.append(y)
  if sort_on_fre!="none":
    d=zip(*sorted(zip(count, L)))
    L=d[1]
    count=d[0]
  return (L,count)

def layout(node):
##    if node.is_leaf() and 'showname' in node.features:
##        # Add node name to laef nodes
##        N = AttrFace("name", fsize=14, fgcolor="black")
##        faces.add_face_to_node(N, node, 0)
    if "Rsize" in node.features:
        # Creates a sphere face whose size is proportional to node's
        # feature "Rsize"
        C = CircleFace(radius=node.Rsize*3, color="RoyalBlue", style="sphere")
        # make the sphere transparent
        C.opacity = 0.9
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")
    if "Rpearsize" in node.features:
        # Creates a sphere face whose size is proportional to node's
        # feature "Rpearsize"
        C = CircleFace(radius=node.Rpearsize*3, color="lightgreen", style="sphere")
        # make the sphere transparent
        C.opacity = 0.9
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")


def scan_internals_pearR(tree,size,threshold,sources="none",simpson_threhold=0.4):
    global t
    #sources is defaulted to be "none"
    import math,seaborn
    import numpy as np
    from scipy.stats import pearsonr,spearmanr
    R_list=[]
    R2_list=[]
    S_index_list=[]
    tree_path=os.path.join(filepath, tree)
    t=Tree(tree_path,format=0)
    internals_dict={}
    internal_nodes=[]
    avoid_sources=["Unknown"] # sources to be omitted
    i=0
    path_trees='time_signal_trees'
    if not os.path.exists(path_trees):
        os.mkdir(path_trees)
    for node in t.traverse():
        if len(node) >= size:
            internal_nodes.append(node)
            dist_list=[]
            year_list=[]
            if sources!="none":
                source_list=[]
            node.add_features(nodetype='internal')
            conf=node.support
            for leaf in node: #may change with different label format
              ###time
                internal_dist=node.get_distance(leaf)
                year_list.append(leaf.name.split('_')[1])
                dist_list.append(internal_dist)
              ###end of time
              ###sources
                if sources!="none":
                    z=leaf.name
                    s=z.split("_")[4]
                    if s not in avoid_sources:
                        source_list.append(s)
              ####end of sources
            len_leaves=len(year_list)
            x_years=np.asarray(year_list).astype(np.int)
            y_dists=np.asarray(dist_list)
            R,P=spearmanr(x_years,y_dists)
            Rpear,Ppear=pearsonr(x_years,y_dists)
            ###for sources
            if sources!="none":
                source_names,source_fre=Uniq(source_list)
                s_index=simpson(source_fre)
                S_index_list.append(s_index)
            ###end of sources
            if math.isnan(R)!=True:
                if R*R >= threshold:
                    i+=1
                    nodetree=str(i)+'_R2_'+str(round(R*R,2))+'.tree'
                    node.write(outfile=filepath+'/'+path_trees+'/'+nodetree,format=0)
                    nt=Tree(filepath+'/'+path_trees+'/'+nodetree)
                    leaves=[leaf.name.replace("'","") for leaf in nt]
                    leaves_num=len(leaves)
                    leave_first=leaves[0].split('_')[0]
                    leavesfile=open(filepath+'/'+path_trees+'/'+nodetree+'.'+str(leaves_num)+'.'+leave_first+'.leaves.txt','w')
                    leavesfile.write("\n".join(leaves))
                    internals_dict[node]=R,P
                    node.add_features(Rsize=int(R*R*50))
                    R2_text=TextFace('R2='+str(round(R*R,2)))
                    #node.add_face(R2_text,column=0,position='branch-top')
                    leaves_text=TextFace('Leaves='+str(len_leaves))
                    #node.add_face(leaves_text,column=0,position='branch-bottom')
                    R2_list.append(R*R)
                    for leaf in node:
                        leaf.add_features(showname=True)
                elif Rpear*Rpear >= threshold:
                    i+=1
                    nodetree=str(i)+'_R2_'+str(round(Rpear*Rpear,2))+'.tree'
                    node.write(outfile=filepath+'/'+path_trees+'/'+nodetree,format=0)
                    nt=Tree(filepath+'/'+path_trees+'/'+nodetree)
                    leaves=[leaf.name.replace("'","") for leaf in nt]
                    leaves_num=len(leaves)
                    leave_first=leaves[0].split('_')[0]
                    leavesfile=open(filepath+'/'+path_trees+'/'+nodetree+'.'+str(leaves_num)+'.'+leave_first+'.leaves.txt','w')
                    leavesfile.write("\n".join(leaves))
                    internals_dict[node]=Rpear,Ppear
                    node.add_features(Rpearsize=int(Rpear*Rpear*50))
                    R2_text=TextFace('R2='+str(round(Rpear*Rpear,2)))
                    #node.add_face(R2_text,column=0,position='branch-top')
                    leaves_text=TextFace('Leaves='+str(len_leaves))
                    #node.add_face(leaves_text,column=0,position='branch-bottom')
                    R2_list.append(Rpear*Rpear)
                    for leaf in node:
                        leaf.add_features(showname=True)
                else:
                    internals_dict[node]=R,P
                    R2_list.append(R*R)
            ###for sources
            if sources!="none":
                if s_index<=simpson_threhold:#more clonal, low diversity
                    nstyle["hz_line_color"] = "blue"
                    node.set_style(nstyle)
                    source_text=TextFace('S='+str(round(s_index,2)),fgcolor="blue",fsize=15)
                    node.add_face(source_text,column=1,position='branch-bottom')
                else:
                    nstyle["hz_line_color"] = "green"
                    node.set_style(nstyle)
                    source_text=TextFace('S='+str(round(s_index,2)),fgcolor="green",fsize=15)
                    node.add_face(source_text,column=1,position='branch-bottom')
            ###end of sources
    ###for time
##    seaborn.set(style="white", palette="muted", color_codes=True)
##    sns_plot=seaborn.distplot(np.array(R2_list),rug=True)
##    fig = sns_plot.get_figure()
##    fig.savefig(os.path.join(filepath, tree.rsplit('.')[0]+"_R2_distribution.png"))
##    sns_plot.clear()
    ###end of time
    ###for source
    if sources!="none":
        seaborn.set(style="white", palette="muted", color_codes=True)
        sns_plot2=seaborn.distplot(np.array(S_index_list),rug=True)
        fig2 = sns_plot2.get_figure()
        fig2.savefig(filepath+tree.rsplit('.')[0]+"_simpson_index_distribution.png")
        sns_plot2.clear()
    ###end of source
    ns = NodeStyle()
    ns["vt_line_width"] = 2
    ns["hz_line_width"] = 2
    ns["size"] = 0
    for node in t.traverse():
        node.set_style(ns)
    ts=TreeStyle()
    ts.layout_fn = layout
    ts.mode = "c"
    ts.scale =180
    ts.show_leaf_name = False
    ts.force_topology=True
    ts.allow_face_overlap=True
    #ts.branch_vertical_margin=2
    #t.render(filepath+tree.rsplit('.')[0]+"_time_signals.png",dpi=300,tree_style=ts)
    outpath=os.path.join(filepath, tree.rsplit('.')[0]+"_time_signals.pdf")
    t.render(outpath,tree_style=ts)

scan_internals_pearR(tree='SE_SNP_tree_msa_phyml.tree',size=20,threshold=0.4,sources="none",simpson_threhold=0.4)

