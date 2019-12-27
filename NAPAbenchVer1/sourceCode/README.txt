------------------------------------------------------------------------
Network Family Geneating algorithm
------------------------------------------------------------------------
This algorithm generates a set of evolutionary related networks based on
the algorithm presented in the following paper:

- S.M.E. Sahraeian and B.J. Yoon, “A Network Synthesis Model for
Generating Protein Interaction Network Families”, PLoS ONE, 2012, to
appear.

The algorithm is developed by S.M.E. Sahraeian and Byung-Jun Yoon 
in GSP lab, Department of Electrical and Computer Engineering, 
Texas A&M University.

------------------------------------------------------------------------

The source code has been made freely available as PUBLIC DOMAIN
software and hence is not subject to copyright in the  United
States.  This system and/or any portion of  the  source  code
may be used, modified, or redistributed without restrictions.  
The software is distributed WITHOUT WARRANTY, express or implied.
The authors accept NO LEGAL LIABILITY OR  RESPONSIBILITY  for
loss due to reliance on the program. 

------------------------------------------------------------------------
The source code is available in Matlab and Python. (Hint: Matlab 
implementation is more efficient.)

Matlab source code: GenerateNetworkFamily.m
Python surce code: GenerateNetworkFamily.py

--------------------------Matlab implementaion--------------------------

	Usage:
	GenerateNetworkFamily(M,tree_file,Na,growth_model,out_path, grow_parameters)

	inputs  M               Number of networks
			tree_file       The underlying phylogentic tree of the networks
							in Newick format
			Na              Number of nodes in the ancestral network
			growth_model    The duplication algorithm to be used to grow
							the networks. The options are 'DMC', 'DMR', and
							'CG' schemes.
			out_path        Output directory where the generated networks
							files will be placed.
			grow_parameters The parameters of the growth model: 
							for DMC: grow_parameters=[q_con, q_mod] 
							for DMR: grow_parameters=[q_new, q_del] 
							for CG: grow_parameters=delta
	output files:
			A.net,..        Network files. These files define the
							structure of each of generated networks.
			A.fo,..       	Functional annotation Network files. These
							files include the functional ortholgy group of
							each network node.
			A-B.sim,..   	The similarity scores of nodes across the
							networks.

	Example:
		GenerateNetworkFamily(5,'test\my_tree.txt',200,'DMC','test\out')

--------------------------Python implementaion--------------------------

	inputs  <M>             	Number of networks
			<tree.txt>      	The underlying phylogentic tree of the networks
			<Na>            	Number of nodes in the ancestral network
			<growth_model>  	The duplication algorithm to be used to grow
								the networks. The options are 'DMC', 'DMR', and
								'CG' schemes
			<output_path>   	Output directory where the generated networks
								files will be placed
			[options]	--q_con	q_con parametre in DMC growth model (sould be between 0 and 1). [Default: 0.1]
						--q_mod	q_mod parametre in DMC growth model (sould be between 0 and 1). [Default: 0.6]
						--q_new	q_new parametre in DMR growth model (sould be between 0 and 1). [Default: 0.12]
						--q_del	q_del parametre in DMR growth model (sould be between 0 and 1). [Default: 0.635]
						--delta	delta parametre in CG  growth model (should be an integer). [Default: 4]
						-h		Help. Print Usage.
	output files
			A.net,...     	Network files. These files define the
							structure of each of generated networks.
			A.fo,...      	Functional annotation Network files. These
							files include the functional ortholgy group of
							each network node.
			A-B.sim, ...    The similarity scores of nodes across the
							networks

	Example: 
		GenerateNetworkFamily.py 5 test\my_tree.txt 200 DMC test\out

------------------------------------------------------------------------  
	Input tree file format (Newick format):

	For genrating the following tree:

               /\
           100/  \50
             A   /\
             150/  \110
               /    \
              B     /\
                 60/  \70
                  C   /\
                   30/  \80
                    D    E

		The tree file should be written as:
		(A:100,(B:150,(C:60,(D:30,E:80):70):110):50)
		Then, if the Number of nodes in the ancestral network (Na) is 200,
		the networks will be of sizes:
		|A|=300, |B|=400, |C|=420, |D|=460, |E|=510

	Onput file format:

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

------------------------------------------------------------------------  
For more information on the algorithms, please see:

S.M.E. Sahraeian and B.J. Yoon, “A Network Synthesis Model for Generating
Protein Interaction Network Families”, PLoS ONE, 2012, to appear.

By Sayed Mohammad Ebrahim Sahraeian and Byung-Jun Yoon
July 2012
Contact: msahraeian@tamu.edu, bjyoon@ece.tamu.edu
-------------------------------------------------------------------------
