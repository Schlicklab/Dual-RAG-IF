# RNA Inverse Folding Using Dual Graph Representations (Schlick Lab)

Dual-RAG-IF finds RNA sequences with the minimal mutations that fold onto a (closely related) target 2D fold described as a dual graph. In dual graphs, double-stranded stems are represented as vertices, and single-stranded loops are edges. For information on RAG (RNA-As-Graphs) using a tree and dual graphs, see our website [RNA-As-Graph](http://www.biomath.nyu.edu/?q=rag/home) and also papers whose pdf files are downloadable from our publication page [Publications](http://www.biomath.nyu.edu/?q=content/schlick-group-publication-list): [Gan et al 2003](https://academic.oup.com/nar/article/31/11/2926/1220259) and other more recent papers on dual graphs. This algorithm in RAG-IF is based on the 2020 initial RAG-IF program [Jain et al., J Struct. Biol., 2020](https://www.sciencedirect.com/science/article/pii/S1047847719302679?via%3Dihub), which uses a genetic algorithm to mutate an initial sequence, twist its 2D folding at each step, and thereby achieve collectively a large number of sequences, which fold onto the desired fold. All such sequences are then sorted by increasing number of mutations. See details in [Jain et al 2020](https://www.sciencedirect.com/science/article/pii/S1047847719302679?via%3Dihub). We have successfully applied Dual-RAG-IF to find structure-alternating mutations for the SARS-CoV-2 frameshifting element [Schlick et al., Biophys. J., 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7575535/); [Schlick et al., JACS, 2021](https://pubs.acs.org/doi/10.1021/jacs.1c03003). See also commentary in [Biophis. J. 2021](https://www.sciencedirect.com/science/article/pii/S0006349521000898).

Dual-RAG-IF has three steps:
(1) Identify mutation regions, which can be done manually or automatically.
For manual design, the input file needs to contain a target 2D structure (in dot-bracket format) compatible with the target dual graph, and a sequence specifying the mutation regions by writing the residues as 'N'.
For automatic design, only the target dual graph is needed.
(2) Produce candidate sequences with mutations using a genetic algorithm (GA). In each step, we test the output fold by applying two 2D folding programs that can handle pseudoknots.
(3) Optimize the candidate sequences by sorting and retaining only minimal or essential mutations.

## Software requirements

RNA 2D folding prediction (command-line version): 

* [NUPACK (v3.2.2)](http://nupack.org/downloads)
* [IPknot](https://github.com/satoken/ipknot)
* [PKNOTS (v1.2)](http://eddylab.org/software.html)

RNA 2D structure format conversion:

[RNAstructure](https://rna.urmc.rochester.edu/Text/index.html)


## Usage

Dual-RAG-IF can be run with:
```
python dualRAGIF.py originalRNA.ct designMethod target [templateSeq=?.out k=?]
```

Input:

* originalRNA.ct: the original sequence's 2D structure in ct format
* designMethod: 1 for automatic design, or 2 for manual design
* target: a dual graph ID (see our tree and dual graph libraries ordered by graph IDs in [Dual Graphs Database](http://www.biomath.nyu.edu/?q=dual_vertices.php) for automatic design, or a two-line design file containing a target 2D structure (in dot-bracket format) and a sequence with mutation regions in 'N'
* templateSeq [optional]: a template sequence
* k [optional]: folding prediction program (1 for PKNOTS, 2 for NUPACK, 3 for IPknot)

Output:
[graphID]min_mut: the final optimized mutation file

The expected run time is 2-12 hours, for sequences of size about 80 nt.

## Examples

Two sample runs are provided for the 77nt SARS-CoV-2 frameshifting element in ./Example/, one using automatic design for 3_3 target dual graph with command:
```
python dualRAGIF.py 77nt_FSE.ct 1 3_3 tmpf=tmpf.out k=3
```
and the other using manual design for 3_5 target dual graph with command:
```
python dualRAGIF.py 77nt_FSE.ct 2 3_5inpf tmpf=tmpf.out k=3
```

Since dual graph 3_3 corresponds to two different RNA motifs, two result files are generated: 3_3_1min_mut, 3_3_2min_mut. To read the results, each line writes a possible mutation combination, for example: [70G-U] means mutating residue 70G to U. For target 3_5 dual graph, only one result file 3_5min_mut. The 3_5 mutant results are published in [Figure 2 in Schlick et al. Biophys J. 2021](https://ars.els-cdn.com/content/image/1-s2.0-S0006349520308146-gr2_lrg.jpg) ![Figure 2 in Schlick et al. Biophys J. 2021](https://ars.els-cdn.com/content/image/1-s2.0-S0006349520308146-gr2_lrg.jpg).


## References

Gan H.H., Pasquali S., Schlick T. Exploring the repertoire of RNA secondary motifs using graph theory; implications for RNA design. *Nucleic Acids Res.* 2003;31:2926–2943.

Jain S., Tao Y., Schlick T. Inverse folding with RNA-As-Graphs produces a large pool of candidate sequences with target topologies. *J. Struct. Biol.* 2020;209:107438. 

Schlick, T., Zhu, Q., Jain, S. Yan, S. Structure-altering mutations of the SARS-CoV-2 frameshifting RNA element. *Biophys. J.* 2021; 120, 1040–1053.

Schlick, T. et al. To knot or not to knot: Multiple conformations of the SARS-CoV-2 frameshifting RNA element. *J. Amer. Chem. Soc.* 2021; 143, 11404–11422.

Chen, S. Graph, pseudoknot, and SARS-CoV-2 genomic RNA: A biophysical synthesis. *Biophys. J.*, 2021; 120(6), 980–982. 


