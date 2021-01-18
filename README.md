# LAYLA: a method using beam search to design the lay-ups of composite laminates with many plies

--------------------------------------------------------------------------

The files correspond to the article in the Journal Composites Part C: Open Access: 
https://doi.org/10.1016/j.jcomc.2020.100072 (2021)

LAYLA is a deterministic method for optimising composite laminate lay-ups satisfying 
lay-up design guidelines.

--------------------------------------------------------------------------
Requirements:

1. A python IDE running python 3.7+

2. The following libraries should accessible:

	- matplotlib
	- pandas
	- numpy

---------------------------------------------------------------------------
In order to use it:

1. clone or download the repository in your desired folder.

2. Set up your environment with the appropriate libraries.

3. Change the settings and run one of the files used to test LAYLA: 
run_LAYLA_V02.py, run_LAYLA_V02_SSpop.py, run_LAYLA_V02_loop.py, run_LAYLA_vs_BBK.py 

Refer to the documentation for more details.
--------------------------------------------------------------------------
Folder Structure

- src and subfolders contain the files needed to run the code

- populations contains the files storing the stacking sequence populations used for 
testing LAYLA

- results-paper-LAYLA-V02 contains the results and analyses generated for the paper.

- run_LAYLA_V02.py is used for to run LAYLA for optimising one composite laminate lay-up.

- run_LAYLA_V02_SSpop.py is used for testing LAYLA capacity at retrieving populations 
of laminate lay-ups based on their lamination parameters and ply counts.

- run_LAYLA_V02_loop.py is used for testing LAYLA capacity at retrieving laminate 
lay-ups based on their lamination parameters and increasing ply counts.

- run_LAYLA_vs_BBK.py is used for comparing the optimiser LAYLA and the 
'branch-and-bound'-based optimiser of Liu Xiaoyang by running both optimiser for a 
population of lay-up lamination parameters.

--------------------------------------------------------------------------
Version 1.0.0

--------------------------------------------------------------------------
License

This project is licensed under the MIT License. See the LICENSE for details.

--------------------------------------------------------------------------
Acknowledgments

- Terence Macquart, Paul Weaver and Alberto Pirrera, my PhD supervisors.

--------------------------------------------------------------------------
Author:

Noémie Fedon

For any questions or feedback don't hesitate to contact me at: noemiefedon@gmail.com
or through github at https://github.com/noemiefedon/RELAY
