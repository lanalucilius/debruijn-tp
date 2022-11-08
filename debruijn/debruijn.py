#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Lucilius"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Lucilius"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Lucilius"
__email__ = "lucilius.lana@hotmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file,"r") as monfich:
        for i in monfich:
            yield next(monfich).strip("\n")
            next(monfich)
            next(monfich)
    

def cut_kmer(read, kmer_size):
    for i in range(0,len(read)-kmer_size+1):
        yield read[i:i+kmer_size]
        


def build_kmer_dict(fastq_file, kmer_size):
    dict={}
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read,kmer_size):
            if kmer in dict:
                dict[kmer]+=1
            else:
                dict[kmer]=1
    return dict
                
                


def build_graph(kmer_dict):
    G= nx.DiGraph()
    for kmer in kmer_dict:
        G.add_edge(kmer[:-1],kmer[1:],weight=kmer_dict[kmer])
    return G
    


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node==True and delete_sink_node==True:
            graph.remove_nodes_from(path)
        elif delete_entry_node==False and delete_sink_node==False:
            graph.remove_nodes_from(path[1:len(path)-1])
        elif delete_entry_node==False and delete_sink_node==True:
            graph.remove_nodes_from(path[1:len(path)])
        elif delete_entry_node==True and delete_sink_node==False:
            graph.remove_nodes_from(path[:len(path)-1])
    return graph


    


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    # Calcule de l'ecart type des poids moyens des chemins
    ecart_weight=statistics.stdev(weight_avg_list)
    #si il est superieur a 0: on selectionne le chemin avec le poids le plus elevé 
    if(ecart_weight>0):
        path_list.pop(weight_avg_list.index(max(weight_avg_list)))
        remove_paths(graph,path_list,delete_entry_node,delete_sink_node)
    #Calcule de l'ecart type pour la longueur
    ecart_length=statistics.stdev(path_length)
    #si l'écart des poids est nul on vérifie ecart type pour la longueur
    if(ecart_weight==0):
        if(ecart_length==0):
            path_list.pop(randint(0,len(path_length)))
            remove_paths(graph,path_list,delete_entry_node,delete_sink_node)
        if(ecart_length>0):
            path_list.pop(path_length.index(max(path_length)))
            remove_paths(graph,path_list,delete_entry_node,delete_sink_node)
    
    return graph
    

def path_average_weight(graph, path):
    """Compute the weight of a path"""
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list=[]
    weight_avg_list=[]
    path_length=[]
    for path in nx.all_simple_paths(graph,ancestor_node,descendant_node):
        path_list.append(path)
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph,path))
    clean_graph=select_best_path(graph,path_list,path_length,weight_avg_list)
    return clean_graph 


    pass

def simplify_bubbles(graph):
    bubble= False
    for noeud in graph:
        liste_predecesseurs= list(graph.predecessors(noeud))
        if len(liste_predecesseurs) > 1:
            for noeud2 in graph:
                noeud_ancetre= nx.lowest_common_ancestor(graph,noeud,noeud2)
                if noeud_ancetre !=None:
                    bubble=True
                    break
        if bubble == True:
            break
    if bubble:
        graph=simplify_bubbles(solve_bubble(graph,noeud_ancetre,noeud2))
    return graph
    

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    liste_starting_node=[]
    for node in graph:
        if len(list(graph.predecessors(node)))==0:
            liste_starting_node.append(node)
    return liste_starting_node
    

def get_sink_nodes(graph):
    liste_sink_node=[]
    for node in graph:
        if len(list(graph.successors(node)))==0:
            liste_sink_node.append(node)
    return liste_sink_node
    

def get_contigs(graph, starting_nodes, ending_nodes):
    tuple=[]
    for node1 in starting_nodes:
        for node2 in ending_nodes:
            if nx.has_path(graph,node1,node2):
                for path in nx.all_simple_paths(graph,node1,node2):
                    contig=""
                    contig +=path[0][0]
                    for node in path:
                        contig += node[-1]
                    tuple.append([contig, len(contig)])
    
    return tuple

def save_contigs(contigs_list, output_file):
    n=0
    
    with open(output_file, "w") as f:
        
        for contigs in contigs_list:
            
            f.write(f">contig_{n} len={contigs[1]}\n")
            
            f.write(f"{textwrap.fill(contigs[0],width=80)}\n")
            
            n+=1
    pass


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    read_fastq()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == '__main__':
    main()
