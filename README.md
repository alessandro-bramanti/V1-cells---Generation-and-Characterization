**********************************************************
*							                                           *
*		             Dendritic Computation 		               *
* 	  and the Fine Structure of Receptive Fields:        *
*		          A Model of V1 Neurons                      *
*	                                                       *
*		                A.P.Bramanti                         *
*							                                           *
**********************************************************

The MATLAB Scripts in this repository reproduce the main result
of the paper mentioned above - arXiv preprint available

Usage is documented in comments inside each file.
All scripts have been tested with MATLAB R2022b (9.13.0.2049777)


MAIN FILES (they may call ancillary scripts as described in comments):

generate_V1_cells.m 
	generates a cell set and saves it in two .mat files 
	cells and angular distributions. 
full_characterization.m
	performs the complete characterization of the cell set
	including spot-mapping and response to gratings.
	other files are recalled for each characterizations step.
	unneeded steps can be skipped commenting out the corresponding 
	section inside the file

V1_cell_set.mat
	set of 100 cells analyzed in the paper
	(format of the variables described in generate_V1_cell.m)
V1_angular_dists.mat
	angular distributions of excitatory and
	inhibitory inputs in the set of 100 cells
	analyzed in the paper
	(format described in generate_V1_cell.m)
V1_pinwheels_300_300.mat
	orientation "pinwheel" map used for 
	the paper's cell set
	 
