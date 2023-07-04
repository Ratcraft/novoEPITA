from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import argparse

function = "debruijn_build"
description = "Builds a de Bruijn graph from a set of reads."
arguments = ["reads", "k"]
outputs = ["graph"]

def debruijn_build(reads, k):
    """
    Builds a de Bruijn graph from a set of reads.
    :param reads: the reads to build the graph from
    :param k: the k-mer size
    :return: the de Bruijn graph
    """
    # Initialize the graph
    graph = {}
    # Iterate over the reads
    reads = SeqIO.parse(reads, "fasta")

    for read in reads:
        print(1)
        # Iterate over the k-mers of the read
        for i in range(len(read.seq) - k + 1):
            # Get the k-mer
            kmer = str(read.seq[i:i + k])
            # Check if the k-mer is already in the graph
            if kmer not in graph:
                # Add the k-mer to the graph
                graph[kmer] = [0, [], []]
            # Increment the multiplicity of the k-mer
            graph[kmer][0] += 1
            # Check if the k-mer is not the first k-mer of the read
            if i > 0:
                # Get the predecessor
                pred = str(read.seq[i - 1:i + k - 1])
                # Add the predecessor to the k-mer
                graph[kmer][1].append(pred)
            # Check if the k-mer is not the last k-mer of the read
            if i < len(read.seq) - k:
                # Get the successor
                succ = str(read.seq[i + 1:i + k + 1])
                # Add the successor to the k-mer
                graph[kmer][2].append(succ)
    # Normalize the graph
    #graph = debruijn_normalize(graph)
    # Return the graph
    return graph

def main():
    """
    Builds a de Bruijn graph from a set of reads.
    """
    # Parse the arguments
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("reads", help="the reads to build the graph from")
    parser.add_argument("k", type=int, help="the k-mer size")
    parser.add_argument("graph", help="the output graph")
    args = parser.parse_args()
    # Build the graph
    print("Building the de Bruijn graph...")
    graph = debruijn_build(args.reads, args.k)
    # Write the graph to the output file
    with open(args.graph, "w") as file:
        for kmer in graph:
            file.write(kmer + " -> " + ",".join(graph[kmer]) + "\n")

if __name__ == "__main__":
    main()