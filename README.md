------------------------------------------------------------------------
NAPAbench Ver.2 : Network Alignment Performance Assessment benchmark
------------------------------------------------------------------------

Network synthesis models in NAPAbench provide effective means to generate synthetic network families that can be used to rigorously assess the performance of network alignment algorithms. In recent years, the protein-protein-interaction (PPI) databases have been significantly updated, hence the network synthesis models in NAPAbench need to be updated to be able to create synthetic network families whose characteristics are close to those of real PPI networks. In NAPAbench Ver.2, we have updated all parameters so that the our benchmark dataset NAPAbench Ver. 2 can mimic the latest real PPI networks. 

NAPAbench has been made freely available as PUBLIC DOMAIN
software and hence is not subject to copyright in the  United
States.  This system and/or any portion of  the  source  code
may be used, modified, or redistributed without restrictions.  
The software is distributed WITHOUT WARRANTY, express or implied.
The authors accept NO LEGAL LIABILITY OR  RESPONSIBILITY  for
loss due to reliance on the program. 



NAPAbench consists of three suites of datasets:

2way(pairwise) alignment dataset
5way alignment dataset
8way alignment dataset


Each each suite has three subcategories: DMR, DMC, CG, STICKY which are named 
based on the network growth model used to construct that set.
In each category, we have 10 independently generated network family sets.

2way(pairwise) alignment dataset: 
The netwok families in the "2way" set consists of two networks generated 
from an ancerstral network  of size 2000 along the following tree:
(A:1000,B:2000)
so the number of nodeos in each network is as follows:
|A|=3000, |B|=4000

5way alignment dataset: 
The netwok families in the "5way" set consists of five networks generated 
from an ancerstral network  of size 1000 along the following tree:
(A:250,(B:250,(C:250,(D:250,E:250):250):250):250)
so the number of nodeos in each network is as follows:
|A|=1250, |B|=1500, |C|=1750, |D|=2000, |E|=2000

8way alignment dataset: 
The netwok families in the "8way" set consists of eight networks generated 
from an ancerstral network  of size 700 along the following tree:
(((A:100,B:100):100,(C:100,D:100):100):100,((E:100,F:100):100,(G:100,H:100):100):100)
so the number of nodeos in each network is as follows:
|A|=1000, |B|=1000, |C|=1000, |D|=1000, |E|=1000, |F|=1000, |G|=1000, |H|=1000



Each Family consists of the following data files:

Network files: A.net, B.net, ...
	These files define the structure of each of generated networks.
	
Functional annotation files: A.fo, B.fo, ...
	These files include the functional ortholgy group of each network node.
	
Similarity score files: A-B.sim, A-C.sim, B-C.sim, ...
	These files incude the similarity scores of nodes across the networks.

Data file format:

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

Contact: bjyoon@ece.tamu.edu