==============================================================================
    INFLUMAX: Influence Maximization in Continuous Time Diffusion Networks
==============================================================================

The problem of finding the optimal set of source nodes in a diffusion
network that maximizes the spread of information, influence, and diseases in
a limited amount of time depends dramatically on the underlying temporal
dynamics of the network. However, this still remains largely unexplored to
date. To this end, we have developed INFLUMAX, an efficient approximation
algorithm with provable near-optimal performance to find the most
influential set of source nodes in the continuous time influence
maximization problem. The algorithm computes analytically the influence
using CTMCs and exploits submodularity for the optimization.

For more information about the procedure see:
  Influence Maximization in Continuous Time Diffusion Networks
  Manuel Gomez-Rodriguez, Bernhard Schölkopf
  http://www.stanford.edu/~manuelgr/influmax/
  
In order to compile on MacOS: 'make' OR 'make opt'.
In order to compile in Linux: 'make linux' OR 'make opt_linux'. 
You will need to install gfortran to make it work.

For networks operations, our code uses the library SNAP (http://snap.stanford.edu), developed by Jure Leskovec.
For computing dominator trees, we hacked the code at www.cs.princeton.edu/~rwerneck/dominators/,
developed by L. Georgiadis, R. E. Tarjan, and R. F. Werneck. 
For computing matrix exponentials, we hacked the Expokit (http://www.maths.uq.edu.au/expokit), developed by R. B. Sidje</a>.

/////////////////////////////////////////////////////////////////////////////
Parameters:

   -n:Input ground-truth network (one file) (default:'example-network.txt')
   -o:Output file name(s) prefix (default:'network')
   -s:Number of sources (default:1)
   -t:Time window (default:10)
   -e:Matrix exponential. 0:Parlett, 1:Krylov, 2:Hessenberg, 3:Krylov(big)+Hessenberg(small) (default:3)
   -mr:Maximum length for minimum reachable path from source(s) (default:-1, no maximum length)
   -mp:Maximum length for diffusion path from source(s) (default:-1, no maximum length)
   -b:Baselines
    0:INFLUMAX, 1:random sources, 2:outdegree, 3:complete search, 4:given source list, 5:given source file (default:0)
   -ir:Repetitions random sources selection if -b:1 (default:100)
   -ls:List of sources (id1;id2;id3) if -b:4 or source filename if -b:5 (default:'')


/////////////////////////////////////////////////////////////////////////////
Usage:

./influmax -i:example_network.txt

It outputs two text files: influence-average-*.txt and influence-info-*.txt, where * is the argument given with -o. The format of 
the output files is described in the files themselves.

/////////////////////////////////////////////////////////////////////////////

Format gound truth:

The ground truth input file should have two blocks separated by a blank line
- A first block with a line per node. The format of every line is <id>,<name>
- A second block with a line per directed edge. The format of every line is <id1>,<id2>,<alpha value>

/////////////////////////////////////////////////////////////////////////////
Additional Tool:

In addition, we provide:
- generate_nets: It allows to build Kronecker and Forest-Fire networks. Please, run 
without any argument to find out how to use it.
- find_sets: It computes self-dominant sets given a network. Run without arguments to find out how to use it.
- dominator_trees: It computes dominator trees given source set and sink node. Run without arguments to find out how to use it.