import os
import numpy as np
import re
from Bio import Phylo
from io import StringIO
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from screeninfo import get_monitors
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from ete3 import Tree, TreeStyle, NodeStyle
import matplotlib.colors as mcolors
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import mplcursors


def read_data(folder):
    #read all trees and node_weights from input folder
    #fasta_files_list, newick_files_list
    def read_seqname_and_weights_from_fastafile(file_path):
        #file_path: A fasta file containing sequences.
        #output (fasta_dict): A dict containing extracted sequence_names with their weights(abundances)
        with open(file_path, 'r') as file:
            fasta_dict = {}
            for line in file:
                if line.startswith('>'):
                    stripline = line.strip()[1:].split('@')
                    current_key = stripline[0]
                    current_val = 1 if len(stripline) == 1 else float(stripline[1])
                    fasta_dict[current_key] = current_val
        return fasta_dict
    
    fasta_files = [f for f in os.listdir(folder) if f.endswith(".fasta") or f.endswith(".fa")]
    fasta_files = sorted(fasta_files, key=lambda x: int(''.join(filter(str.isdigit, x))))
    newick_files = [f for f in os.listdir(folder) if f.endswith(".nk")]
    newick_files = sorted(newick_files, key=lambda x: int(''.join(filter(str.isdigit, x))))
    fasta_files_list = []
    newick_files_list = []
    for filename in fasta_files:
        fastafilename = folder+filename
        file_dict = read_seqname_and_weights_from_fastafile(fastafilename)
        fasta_files_list.append(file_dict)
    for filename in newick_files:
        newickfilename = folder+filename
        with open(newickfilename, 'r') as file:
            count = 0
            line = file.read().strip()
            while '):' in line:
                line = line.replace( '):', f')uv{count}:', 1)
                count += 1
        newick_files_list.append(line)
    return fasta_files_list, newick_files_list


def simplify_tree_by_node_weights(newick_string, node_weights, weight_threshold):
    tree = Phylo.read(StringIO(newick_string), "newick")

    def process_node(clade, parent=None):
        # Check if the node is below the threshold
        if (clade.name in node_weights and clade.name != 'naive' and node_weights[clade.name] < weight_threshold) or (clade.name is not None and 'uv' in clade.name):
            # If the node is below the threshold, remove it and reconnect children to parent
            if parent:
                for child in clade.clades:
                    child.branch_length = (child.branch_length or 0) + (clade.branch_length or 0)
                    parent.clades.append(child)
            return None  # Mark node for removal
        else:
            # Process children recursively
            new_clades = []
            for child in clade.clades:
                processed_child = process_node(child, clade)
                if processed_child:
                    new_clades.append(processed_child)
            clade.clades = new_clades
            return clade

    # Start processing from the root
    tree.root = process_node(tree.root)

    # Convert the simplified tree back to Newick format
    output = StringIO()
    Phylo.write(tree, output, "newick")
    return output.getvalue().strip()


def simplify_trees(newick_trees, all_node_weights, thresh_method = 'mean'):
    new_newick_trees = []
    n = len(newick_trees)
    for i in range(n):
        tree = newick_trees[i]
        node_weights = all_node_weights[i]
        if thresh_method == 'most':
            k = 5
            if len(node_weights) > k:
                thresh = sorted(node_weights.values())[::-1][k-1]
            else:
                thresh = 1
        elif thresh_method == 'mean':
            thresh = np.mean(list(node_weights.values()))
        elif thresh_method == 'IQR' or thresh_method == 'iqr':
            q1, q3 = np.percentile(list(node_weights.values()), [25, 75])
            iqr = q3 - q1
            thresh = q3 + 1.5 * iqr
        elif thresh_method == 'zscore':
            z = 2
            thresh = z * np.std(list(node_weights.values())) + np.mean(list(node_weights.values()))
        else:
            thresh = np.percentile(list(node_weights.values()),95)
        temp_new_newick = simplify_tree_by_node_weights(tree, node_weights, thresh)
        new_newick = re.sub(r'\)[^)]*?;', ');', temp_new_newick)
        new_newick_trees.append(new_newick)
    return new_newick_trees

def Prepare_results_for_clustering(newicktrees, nodes_and_weights):
    filtered_nodes_and_weights=[]
    phylo_trees = []
    for i in range(len(newicktrees)):
        tree = Phylo.read(StringIO(newicktrees[i]), "newick")
        phylo_trees.append(tree)
        node_names = set(match[0] for match in re.findall(r"([a-zA-Z0-9]+):([\d.]+)", newicktrees[i]))
        filtered_nodeweights = {k: v for k, v in nodes_and_weights[i].items() if k in node_names}
        filtered_nodes_and_weights.append(filtered_nodeweights)
    return phylo_trees, filtered_nodes_and_weights


def create_ancestor_dict(tree):
    ancestor_dict = {}
    for clade in tree.find_clades(order="preorder"):
        # Get path to the current node
        path_to_node = tree.get_path(clade)

        # If the node is not the root, assign its parent as ancestor
        if len(path_to_node) > 1:
            parent = path_to_node[-2]  # Second to last in the path is the parent
            ancestor_dict[clade.name] = parent.name
        else:
            # Root node has no ancestor
            ancestor_dict[clade.name] = None
    if ((next(iter(ancestor_dict))=='uv0' or next(iter(ancestor_dict))==None) and ancestor_dict[next(iter(ancestor_dict))]==None):
        ancestor_dict.pop(next(iter(ancestor_dict)))
    return ancestor_dict


def tree_levels_to_dict(tree):
    #tree: A rooted tree with labeled nodes.
    #output: A dictionary mapping levels to lists of node labels.
    def traverse(node, current_level):
        levels[current_level].append(node.name)
        for child in node.clades:
            traverse(child, current_level + 1)
    levels = defaultdict(list)
    # Start traversal from the root of the tree
    traverse(tree.root, 0)
    alllevels = dict(levels)
    alllevels.pop(0)
    return alllevels


def construct_sortedshapeMat_for_tree(tree, ancDict):
    levels_info = tree_levels_to_dict(tree)
    max_nrow = len(levels_info)
    max_ncol = max([len(val) for k, val in levels_info.items()])
    shmat = np.zeros((max_nrow,max_ncol))
    for k, vals in levels_info.items():
        tempL = []
        valsL = len(vals)
        if valsL==1:
            shmat[k-1, 0] = ancDict.loc[vals[0]].w
        elif (valsL >=2 and len(levels_info[k-1])==1):
            tempL = [ancDict.loc[i].w for i in vals]
            sort_inds = np.argsort(tempL)[::-1]
            tempL = list(np.array(tempL)[sort_inds])
            levels_info[k] = list(np.array(levels_info[k])[sort_inds])
            shmat[k-1,:valsL] = tempL
        else:
            y = []
            z = []
            for node in levels_info[k-1]:
                x = ancDict.index[ancDict.anc == node].tolist()
                tempL = list(ancDict.w[x])
                sort_inds = np.argsort(tempL)[::-1]
                tempL = list(np.array(tempL)[sort_inds])
                z.extend(tempL)
                x = list(np.array(x)[sort_inds])
                y.extend(x)
            levels_info[k] = y
            shmat[k-1,:len(y)] = z
    return shmat, levels_info


def info_from_trees(phylotree_list, fasta_list):
    shapemat_list = []
    ancDict_list = []
    n = len(fasta_list)
    max_score = 0
    for i in range(n):
        wlist = []
        ancDict = create_ancestor_dict(phylotree_list[i])
        for k in ancDict.keys():
            w = fasta_list[i][k]
            wlist.append(w)
        df = pd.DataFrame(list(zip(ancDict.values(),wlist)),columns=['anc','w'],index=list(ancDict.keys()))
        ancDict_list.append(df)
        if len(ancDict) > max_score:
            max_score = len(ancDict)
        shapemat, _ = construct_sortedshapeMat_for_tree(phylotree_list[i], df)
        shapemat_list.append(shapemat)
    return shapemat_list, ancDict_list, max_score


def convert_shmat_to_scoremat(shapemat_list, max_score):
    final_mat_list = []
    for mat in shapemat_list:
        non_zero_elements = mat[np.nonzero(mat)]
        sorted_elements = np.sort(non_zero_elements)[::-1]
        replacement_map = {value: i for i, value in enumerate(sorted_elements, start=1)}
        for value, new_value in replacement_map.items():
            mat[mat == value] = max_score - new_value + 1
        final_mat_list.append(mat)
    return final_mat_list


def disMat_of_mats(mat_list):
    n = len(mat_list)
    dis_mat = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            shape_i = np.shape(mat_list[i])
            shape_j = np.shape(mat_list[j])
            A =np.zeros((max(shape_i[0],shape_j[0]),max(shape_i[1],shape_j[1])))
            B =np.zeros((max(shape_i[0],shape_j[0]),max(shape_i[1],shape_j[1])))
            A[:shape_i[0], :shape_i[1]] = mat_list[i]
            B[:shape_j[0], :shape_j[1]] = mat_list[j]
            dis_mat[i, j] = dis_mat[j, i] = np.sqrt(np.sum((A-B)**2))
    return dis_mat


def rectangular_layout(tree):
    positions = {}
    x_offset = 0
    y_offset = 0
    level_spacing = 50  # Vertical spacing between levels
    sibling_spacing = 100  # Horizontal spacing between siblings

    def assign_positions(node, x, y):
        nonlocal x_offset
        if node.is_leaf():
            positions[node.name] = (x_offset, y)
            x_offset += sibling_spacing
        else:
            child_positions = []
            for child in node.children:
                assign_positions(child, x, y - level_spacing)
                child_positions.append(positions[child.name][0])
            positions[node.name] = (sum(child_positions) / len(child_positions), y)

    assign_positions(tree, x_offset, y_offset)
    return positions


def interactive_plot_bcellTree(newick, node_weights, title_name ,title_num = 0):
    node_info = {}
    for key, val in node_weights.items():
        node_info[key] = f'Node_name: {key}\n\nAbundancy: {int(val)}'
    t = Tree(newick, format=1)
    tree_nodes = {node.name for node in t.traverse() if node.name}
    node_info.update({key: 'Node_name: ''\n\nAbundancy: 1' for key in tree_nodes if key not in node_info})
    node_weights.update({key: 1 for key in tree_nodes if key not in node_weights})
    max_weight = max(node_weights.values())
    node_sizes = {node: (weight / max_weight) * 1000 for node, weight in node_weights.items()}
    weights = np.array(list(node_weights.values()))
    normalized_weights = (weights - weights.min()) / (weights.max() - weights.min())  # Scale between 0 and 1
    colors = cm.viridis(normalized_weights)
    
    # Create a mapping of nodes to colors
    node_colors = {node: colors[idx] for idx, node in enumerate(node_weights.keys())}
    pos = rectangular_layout(t)
    
    # Draw the tree with rectangular edges
    fig = plt.figure(figsize=(8, 5))
    fig.canvas.manager.set_window_title(f'{title_name} T{title_num+1}')

    x_coords = [x for x, y in pos.values()]
    y_coords = [y for x, y in pos.values()]
    x_range = max(x_coords) - min(x_coords)
    y_range = max(y_coords) - min(y_coords)
    x_margin = x_range * 0.3 if x_range > 0 else 1.0
    y_margin = y_range * 0.3 if y_range > 0 else 1.0
    plt.gca().set_xlim(min(x_coords) - x_margin, max(x_coords) + x_margin)
    plt.gca().set_ylim(min(y_coords) - y_margin, max(y_coords) + y_margin)
    
    # Draw nodes
    scatter = plt.scatter([], [], s=[], alpha=0.9, color=[], edgecolor="black", zorder=2)
    scatter.set_offsets([pos[node] for node in tree_nodes])
    scatter.set_sizes([node_sizes[node] for node in tree_nodes])
    scatter.set_color([node_colors[node] for node in tree_nodes])

    for i, node in enumerate(tree_nodes):
        if node == "naive":
            plt.scatter(
                pos[node][0], pos[node][1],
                s=20,
                alpha=0.9,
                color="black",
                edgecolor="black",
                zorder=2, 
                marker="^" 
            )
    
    # Draw edges
    for node in t.traverse("postorder"):
        if not node.is_root():
            parent = node.up
            x_start, y_start = pos[parent.name]
            x_end, y_end = pos[node.name]
            # Draw horizontal line from parent to child
            plt.plot([x_start, x_end], [y_start, y_start], color="black", lw=1, zorder=1)
            # Draw vertical line down to child
            plt.plot([x_end, x_end], [y_start, y_end], color="black", lw=1, zorder=1)

    # Draw labels
    # for node, (x, y) in pos.items():
    #     plt.text(x, y + 10, node, fontsize=8, ha="center", va="center")

    # Display the plot
    plt.axis("off")
    
    # Add mplcursors for interactive annotations
    cursor = mplcursors.cursor(scatter, hover=True)

    @cursor.connect("add")
    def on_add(sel):
        node_name = list(tree_nodes)[sel.index]
        info_text = node_info.get(node_name, "No information available")
        # width = 20  # Adjust the width as per your preference
        # centered_text = "\n".join(line.center(width) for line in info_text.splitlines())
        # sel.annotation.set_text(centered_text)
        sel.annotation.set_text(info_text)
        sel.annotation.set_multialignment('left') 
        bbox_properties = dict(
        boxstyle="round,pad=0.7",
        edgecolor="black",
        facecolor="yellow",
        linewidth=1,
        alpha=0.7
        )
        sel.annotation.set_bbox(bbox_properties)
        sel.annotation.arrowprops = None
        sel.annotation.get_bbox_patch()
        
        @sel.annotation.axes.figure.canvas.mpl_disconnect
        def remove_annotation(event):
            sel.annotation.set_visible(False)
            sel.annotation.axes.figure.canvas.draw_idle()

    plt.show()


def plotTrees(newick_list, node_weights_list, input_trees, all_node_weights, figsize=(12, 12)):
    """
    Plot multiple Newick trees in a single figure and trigger a function when a subplot is clicked.

    Args:
        newick_list (list): List of Newick strings representing trees.
        node_weights_list (list): List of dictionaries with node weights for each tree.
        figsize (tuple): Figure size for the entire plot.
        on_subplot_click (function): Function to execute when a subplot is clicked.
                                      Receives the subplot index as an argument.
    """

    n_trees = len(newick_list)
    cols = int(np.ceil(np.sqrt(n_trees)))
    rows = int(np.ceil(n_trees / cols))

    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    fig.canvas.manager.set_window_title("Minified trees")

    axes = axes.flatten()

    def on_click(event):
        # Check if the click happened within any subplot
        for i, ax in enumerate(axes):
            if ax == event.inaxes:  # Check if the clicked axis matches the current axis
                interactive_plot_bcellTree(newick_list[i], node_weights_list[i], title_name = 'Minified', title_num = i)
                interactive_plot_bcellTree(input_trees[i], all_node_weights[i], title_name = 'Original', title_num = i)
                break

    # Connect the click event to the handler
    fig.canvas.mpl_connect("button_press_event", on_click)
    for i, (newick, node_weights) in enumerate(zip(newick_list, node_weights_list)):
        ax = axes[i]
        t = Tree(newick, format=1)
        max_weight = max(node_weights.values())
        node_sizes = {node: (weight / max_weight) * 500 for node, weight in node_weights.items()}
        weights = np.array(list(node_weights.values()))
        normalized_weights = (weights - weights.min()) / (weights.max() - weights.min())
        colors = cm.viridis(normalized_weights)
        node_colors = {node: colors[idx] for idx, node in enumerate(node_weights.keys())}
        pos = rectangular_layout(t)

        x_coords = [x for x, y in pos.values()]
        y_coords = [y for x, y in pos.values()]
        x_range = max(x_coords) - min(x_coords)
        y_range = max(y_coords) - min(y_coords)
        x_margin = x_range * 0.3 if x_range > 0 else 1.0
        y_margin = y_range * 0.3 if y_range > 0 else 1.0
        ax.set_xlim(min(x_coords) - x_margin, max(x_coords) + x_margin)
        ax.set_ylim(min(y_coords) - y_margin, max(y_coords) + y_margin)

        #ax.set_title(f"Tree {i+1}", fontsize=10)
        for node, (x, y) in pos.items():
            if node != "":
                size = node_sizes[t.search_nodes(name=node)[0].name]
                color = node_colors[t.search_nodes(name=node)[0].name]
                if node == "naive":
                    ax.scatter(x, y, s=12, alpha=0.9, color="black", edgecolor="black", zorder=2, marker="^")
                else:
                    ax.scatter(x, y, s=size, alpha=0.9, color=color, edgecolor="black", zorder=2)

        for node in t.traverse("postorder"):
            if not node.is_root():
                parent = node.up
                x_start, y_start = pos[parent.name]
                x_end, y_end = pos[node.name]
                ax.plot([x_start, x_end], [y_start, y_start], color="black", lw=1, zorder=1)
                ax.plot([x_end, x_end], [y_start, y_end], color="black", lw=1, zorder=1)
        ax.text(0.5, -0.2, f"Tree {i+1}", transform=ax.transAxes, ha='center', fontsize=7)
        ax.axis("off")

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    return fig


def hierclust0(tree, num_clusters):
    all_nodes = [node for node in tree.traverse() if node != tree.get_tree_root()]
    node_labels = [node.name if node.name else f"Node_{i}" for i, node in enumerate(all_nodes)]
    num_nodes = len(all_nodes)
    distance_matrix = np.zeros((num_nodes, num_nodes))
    for i, node1 in enumerate(all_nodes):
        for j, node2 in enumerate(all_nodes):
            distance_matrix[i, j] = tree.get_distance(node1, node2)
    distance_matrix = (distance_matrix + distance_matrix.T) / 2
    condensed_matrix = squareform(distance_matrix)
    linkage_matrix = linkage(condensed_matrix, method="average")
    cluster_assignments = fcluster(linkage_matrix, num_clusters, criterion='maxclust')
    clusters = {}
    for node_label, cluster_id in zip(node_labels, cluster_assignments):
        clusters.setdefault(cluster_id, []).append(node_label)
    leaf_clusters = {}
    for cluster_id, nodes in clusters.items():
        leaves_only = [node for node in nodes if tree.search_nodes(name=node) and tree.search_nodes(name=node)[0].is_leaf()]
        if leaves_only:
            leaf_clusters[cluster_id] = leaves_only
    leaf_clusters = {i + 1: value for i, value in enumerate(leaf_clusters.values())}
    return leaf_clusters


def hierclust(tree, num_clusters):
    internal_nodes = 0
    for node in tree.traverse():
        if not node.is_leaf():
            if not all(child.dist == 0 for child in node.children):
                internal_nodes += 1
    if num_clusters > internal_nodes:
        num_clusters = internal_nodes
    k = num_clusters
    while True:
        result = hierclust0(tree, k)
        if len(result) == num_clusters:
            return result
        k += 1


def plot_with_ete3(tree, num_clusters, output_cluster_file, show_clusters):

    def contract_zero_length(tree):
        for node in tree.traverse("postorder"):  # Postorder ensures children are processed before parents
            if not node.is_leaf() and node.dist == 0:
                parent = node.up
                if parent is not None:  # Ensure the parent exists
                    for child in node.get_children():
                        child.detach()
                        parent.add_child(child)
                    node.detach()

    # Remove clade names from internal nodes
    for clade in tree.get_nonterminals():
        clade.name = None

    # Convert to Newick and load into ete3 Tree
    newick_str = tree.format("newick")
    newick_str = newick_str[:newick_str.find(',Outgroup')]+newick_str[newick_str.rfind(')'):]
    ete_tree = Tree(newick_str)

    # Contract zero-length edges
    contract_zero_length(ete_tree)
    clusters = hierclust(ete_tree, num_clusters)
    k = len(clusters)
    leaf_num = len(set().union(*clusters.values()))
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_scale = False
    # if leaf_num > 20:
    #     ts.rotation = 90
    # Assign unique styles to each cluster
    # colors = ["red", "blue", "green", "purple", "orange"]
    cmap = plt.get_cmap("tab10")  # Use the colormap directly
    colors = [mcolors.rgb2hex(cmap(i / k)) for i in range(k)]
    for key, vals in clusters.items():
        for val in vals:
            node = ete_tree.search_nodes(name = val)[0]
            ns = NodeStyle()
            ns["fgcolor"] = colors[key % len(colors)]  # Assign a color
            ns["size"] = 4
            node.set_style(ns)

    for node in ete_tree.traverse():
        if not node.is_leaf() and not hasattr(node, 'style_applied'):  # Check if already styled
            ns = NodeStyle()
            ns["fgcolor"] = 'black'
            ns["size"] = 1.5
            node.set_style(ns)
    with open(output_cluster_file, 'w') as file:
        for key, values in clusters.items():
            file.write(f"{key}:\t{', '.join(values)}\n")
    #ete_tree.render('fig_nj_clustering.png', tree_style=ts, w=500, h=800, units="px")
    #ete_tree.render("fig_nj_clustering.svg", tree_style=ts)
    ete_tree.render("fig_nj_clustering.pdf", tree_style=ts)
    if show_clusters:
        ete_tree.show(tree_style=ts)


def nj_clustering_with_outgroup(distance_matrix, num_clusters, output_cluster_file, show_clusters, plot_method, out_dist):
    # Prepare the distance matrix
    list_distance_matrix = distance_matrix.tolist()
    names = ['T' + str(i) for i in range(1, np.shape(distance_matrix)[0] + 1)]
    outgroup_label = "Outgroup"
    names.append(outgroup_label)

    # Set distances to outgroup
    outgroup_distances = [out_dist] * len(list_distance_matrix)
    for i in range(len(list_distance_matrix)):
        list_distance_matrix[i].append(outgroup_distances[i])
    outgroup_distances.append(0)
    list_distance_matrix.append(outgroup_distances)

    lower_tri_dist_matrix = [[list_distance_matrix[i][j] for j in range(i + 1)] for i in range(len(list_distance_matrix))]
    matrix = DistanceMatrix(names, lower_tri_dist_matrix)

    # Perform Neighbor Joining
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(matrix)
    tree.root_with_outgroup(outgroup_label)

    if plot_method == 'ete':
        plot_with_ete3(tree, num_clusters, output_cluster_file, show_clusters)
        return None
    else:
        # Default matplotlib plot
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(1, 1, 1)
        def label_func(clade):
            return clade.name if clade.is_terminal() else ""
        tree.prune(outgroup_label)
        Phylo.draw(tree, label_func=label_func, do_show=False, axes=ax)
        return fig


def get_monitor_sizes():
    monitor = get_monitors()[0]
    screen_width, screen_height = monitor.width, monitor.height
    adjusted_width = screen_width / 100
    adjusted_height = screen_height / 100
    return adjusted_width, adjusted_height


def miniclust(folder, output_dismat_file, output_cluster_file, num_clusters = 4, show_clusters = False, show_simp_trees = False, b = 0):
    all_node_weights, input_trees =read_data(folder)
    new_trees = simplify_trees(input_trees, all_node_weights, thresh_method = 'mean')
    phylo_new_trees, new_trees_nodes_weights = Prepare_results_for_clustering(new_trees, all_node_weights)
    shapemat_list, _, max_score = info_from_trees(phylo_new_trees, new_trees_nodes_weights)
    scoremat_list = convert_shmat_to_scoremat(shapemat_list, max_score + b)
    distance_mat = disMat_of_mats(scoremat_list)
    normalized_distance_mat = distance_mat/np.max(distance_mat)
    dismat_df = pd.DataFrame(np.round(normalized_distance_mat,5), index = range(1,np.shape(normalized_distance_mat)[0]+1), columns = range(1,np.shape(normalized_distance_mat)[0]+1))
    dismat_df.to_csv(output_dismat_file, index=True, sep='\t')
    if show_simp_trees == True:
        width, height = get_monitor_sizes()
        fig_all = plotTrees(new_trees, new_trees_nodes_weights, input_trees, all_node_weights, figsize = (width, height))
        fig_all.savefig("fig_all.png")
        plt.show()
    fig_nj_clustering = nj_clustering_with_outgroup(normalized_distance_mat, num_clusters, output_cluster_file, show_clusters, plot_method='ete', out_dist=1.01)
    if fig_nj_clustering is not None: fig_nj_clustering.savefig("fig_nj_clustering.png")
    return new_trees, new_trees_nodes_weights