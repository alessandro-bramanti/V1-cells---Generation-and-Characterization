"Dendritic Computation and the Fine Structure of Receptive Fields: a Model of V1 Neurons"
 
A.P.Bramanti                         



The MATLAB Scripts in this repository reproduce the main result
of the paper mentioned above - arXiv preprint available

Usage is documented in comments inside each file.
All scripts have been tested with MATLAB R2022b (9.13.0.2049777)

DISCLAIMER: The code is provided "as is. Developer makes no warranties, 
express or implied, and hereby disclaims all implied warranties, 
including any warranty of merchantability and warranty of fitness for a particular purpose.

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


V1_pinwheels_300_300.mat
	orientation "pinwheel" map used for 
	the paper's cell set
	 
