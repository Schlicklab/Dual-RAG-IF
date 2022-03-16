# Dual-RAG-IF
RNA inverse folding using dual graph representations (Schlick Lab)

Dual-RAG-IF finds the minimal mutations that make an RNA sequence fold onto a (closely related) target dual graph, where double-stranded stems are represented as vertices and single-stranded loops as edges. This algorithm is based on the 2020 initial RAG-IF program (Jain et al., J Struct. Biol., 2020). We have successfully applied Dual-RAG-IF to find structure-alternating mutations for the SARS-CoV-2 frameshifting element (Schlick et al., Biophys. J., 2021; Schlick et al., JACS, 2021).

Dual-RAG-IF has three steps: \
(1) Identify mutation regions, which can be done manually or automatically. \
For manual design, the input file needs to contain a target 2D structure (in dot-bracket format) compatible with the target dual graph, and a sequence specifying the mutation regions by writing the residues as 'N'. \
For automatic design, only the target dual graph is needed. \
(2) Produce candidate sequences with mutations using a genetic algorithm (GA). \
(3) Optimize the candidate sequences by sorting and retaining only minimal or essential mutations.

The 2D folding prediction programs used are PKNOTS, NUPACK, and IPknot. The user needs to install command-line version NUPACK and IPknot (installing PKNOTS is optional): \
NUPACK: http://nupack.org/downloads (version 3.2.2 used) \
IPknot: https://github.com/satoken/ipknot \
PKNOTS: http://eddylab.org/software.html (version v1.2 used)

An RNA 2D structure format conversion program is also needed: \
RNAstructure command-line version: https://rna.urmc.rochester.edu/Text/index.html

The command needed to run Dual-RAG-IF is: \
&nbsp;&nbsp;&nbsp;&nbsp;% python dualRAGIF.py originalRNA.ct designMethod target [templateSeq=?.out k=?] \
&nbsp;&nbsp;&nbsp;&nbsp;@ originalRNA.ct contains the original sequence's 2D structure in ct format \
&nbsp;&nbsp;&nbsp;&nbsp;@ designMethod is 1 for automatic design, or 2 for manual design \
&nbsp;&nbsp;&nbsp;&nbsp;@ target is a dual graph ID if automatic design used, or a two-line design file containing a target 2D structure (in dot-bracket format) and a sequence with mutation regions in 'N' \
&nbsp;&nbsp;&nbsp;&nbsp;@ [optional] templateSeq=a file containing a template sequence \
&nbsp;&nbsp;&nbsp;&nbsp;@ [optional] k=folding prediction program (1 for PKNOTS, 2 for NUPACK, 3 for IPknot)

The final optimized mutation results are written in a file named after the target dual graph: [graphID]min_mut. \
The expected run time is 2-12 hours.
      
Two example runs are provided for the 77nt SARS-CoV-2 frameshifting element in the "Example" folder, one using automatic design for 3_3 target dual graph with command: \
&nbsp;&nbsp;&nbsp;&nbsp;% python dualRAGIF.py 77nt_FSE.ct 1 3_3 tmpf=tmpf.out k=3 \
and the other using manual design for 3_5 target dual graph with command: \
&nbsp;&nbsp;&nbsp;&nbsp;% python dualRAGIF.py 77nt_FSE.ct 2 3_5inpf tmpf=tmpf.out k=3
      
Since dual graph 3_3 corresponds to two different RNA motifs, two result files are generated: 3_3_1min_mut, 3_3_2min_mut. To read the results, each line writes a possible mutation combination, for example: [70G-U] means mutating residue 70G to U. For target 3_5 dual graph, only one result file 3_5min_mut.


