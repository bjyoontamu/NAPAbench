NAPAbench Ver 2.0
Contact: bjyoon@ece.tamu.edu, larcwind@tamu.edu

Input
	General Parameters
		1. The number of networks: The number of networks in each dataset. It should be positive decimal, and greater or equal to 2.
		2. The number of dataset: The number of sample dataset. It should be positive decimal, and greater or equal to 1. 
		3. Phylogenetic tree file path: The underlying phylogenetic tree of the networks in Newick format. *
		4. Output directory where the generated dataset will be placed.
	Network Growth Model Parameters
		1. Network growth mode: DMR, DMR, CG, STICK
		2. Network growth model parameters
		3. The ancestral network size: the size of seed network. 

Output **
	A.net, ..        Network files. These files define the structure of each of generated networks.
	A.fo, ..       Functional annotation Network files. These files include the functional ortholgy group of each network node.
	A-B.sim, ..   The similarity scores of nodes across the networks.

Useful tips
	Tuning parameters for the network growth model
		
** Output file format
	Network files: [e.g A.net]
              a1	a2
              a3	a1
              a4	a2
              a2	a3
	Functional annotation files: [e.g. A.fo]
              a1	FO:1
              a2	FO:2
              a1	FO:2
              a4	FO:3
	Similarity score files: [e.g. A-B.sim]
              a1	b1	153
              a1	b3	55
              a1	b7	49
              a2	b3	444
              a3	b3	211
              a3	b4	122
              a4	b5	251
              a4	b8	71