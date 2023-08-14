import pandas as pd
from tqdm import tqdm
import ete3


def load_tree(newick):
    """
    Loads the tree and returns the list of nodes and the list of names of the nodes.

    Parameters
    ----------
    t : ete3.TreeNode
        Tree to be analyzed.

    Returns
    -------
    t : ete3.TreeNode
        Tree to be analyzed.
    node_list : list
        List of nodes in the tree.
    names_list : list
        List of names of the nodes in the tree.

    """
    t = ete3.Tree(newick)
    
    node_list = []
    names_list=[]


    n=0
    for node in t.traverse("postorder"):
        n=n+1
        node_list.append(node)
        if node.name=='':
            node.name=n
            names_list.append(n)
        else:
            names_list.append(node.name)

    return t, node_list, names_list


def get_distance_matrix(t, node_list, names_list):
    """
    Returns the distance matrix between all pairs of nodes in the tree.
    
    Parameters
    ----------
    t : ete3.TreeNode
        Tree to be analyzed.
    node_list : list
        List of nodes in the tree.
    names_list : list
        List of names of the nodes in the tree.

        
    Returns
    -------
    ist_matrix : pandas.DataFrame
        Distance matrix between all pairs of nodes in the tree.
        
        """
    dist_matrix = pd.DataFrame(index=names_list, columns=names_list)


    for i in tqdm(range(len(node_list))):
        for j in range(i+1, len(node_list)):
            nodei=node_list[i]
            nodej=node_list[j]
            distance = t.get_distance(nodei, nodej)
            dist_matrix.loc[names_list[i], names_list[j]]=distance

    dist_matrix.to_csv("dist_matrix_all_incl_internal.csv")

    return dist_matrix


def find_nodes_with_descendants_from_matrix(leaves, tree, dist_matrix):
    """
    Returns the list of nodes with all descendants in the input leaf set.
    
    Parameters
    ----------
    leaves : list
        List of leaves in the tree which include a feature.
    tree : ete3.TreeNode
        Tree to be analyzed.
    dist_matrix : pandas.DataFrame
        Distance matrix between all pairs of nodes in the tree.
        
    Returns
    -------
    sum_dist : float
        Sum of distances between the ancestor leaves.
    """
    # Find the set of leaf names
    leaf_set = set(leaves)
    # Initialize a list to store the nodes with all descendants in the input leaf set
    nodes_with_descendants = []
    # Traverse the tree and check if each of node's descendants is in the input leaf set
    nodes=[]
    
    for node in tree.traverse("levelorder"):
        if set(node.get_leaf_names()).issubset(leaf_set):
            #print(set(node.get_leaf_names()))
            nodes_with_descendants.append(node.name)
            leaf_set=leaf_set-set(node.get_leaf_names())
            
    #return nodes_with_descendants
    nodes_with_descendants = [str(i) for i in nodes_with_descendants]
    # Calculate the distances between the ancestors
    distances = dist_matrix.loc[nodes_with_descendants, nodes_with_descendants]
    sum_dist = distances.sum().sum()

    
    return sum_dist

