%{
****************************************************
*							                       *
*		   Dendritic Computation 		           *
* 	and the Fine Structure of Receptive Fields:    *
*		A Model of V1 Neurons			           *
*                       					       *
*		      A.P.Bramanti			               *
*							                       *
****************************************************

The MATLAB Scripts in this repository reproduce the main result
of the paper mentioned above - arXiv preprint available.
Please cite the paper whenever using the results of these scripts or
modified versions.

DISCLAIMER: The code is provided "as is. Developer makes no warranties, 
express or implied, and hereby disclaims all implied warranties, 
including any warranty of merchantability and warranty of fitness for a particular purpose.

Generates a new V1 cell set

%}

%close all
%clear variables

fprintf('Setting and saving basic constants...\n');

pixel_fact=1;

X_size=100;
Y_size=100;

angular_steps=10;

dummy_pixel_x=1;
dummy_pixel_y=1;
dummy_pixel_z=1;

%% ******************** LGN cells ************************

fprintf('Setting and saving LGN cells...\n');

global LGN_collinear_Gabor

inner_diameters_ON=5; %7*2; % optimum diameter (max response) from LGN spatial feedback and properties, Andolina et al. 2013: the diameter is 0.7 degs in the paper, so there is factor of 20
outer_diameters_ON=10; %20*2; % see comment above: 2 is the maximum diameter estimated from the figure;

C50=0.15;

LGN_contrast_gain=5; % coefficient in logistic function
LGN_contrast_threshold=0.18; 

LGN_margin=ceil(max(outer_diameters_ON/pixel_fact));

Gabor_transverse_sfs=0.5;
Gabor_collinear_radii=6; % 2
Gabor_collinear_sfs=8;
Gabor_transverse_radii=4; % 2

load LGN_cells.mat

LGN_collinear_Gabor=generate_LGN_to_V1_Gabor_filters(Gabor_collinear_radii,Gabor_collinear_sfs,Gabor_transverse_radii,angular_steps);
                                                                                                    
%% ************* V1 cell basic parameters *****************

global g_L g_a g_NMDA g_SOMA
global dendrite_number
global dendrite_n_comp
global ex_syn_threshold

ex_syn_threshold=0; % threshold of curve in the last comparment; if the sum of active ex connections on a dendrite is smaller, the dendrite will not be computed and its potential will be set to 0 to speed up computation

dendrite_number=5; % dendrites per cell
dendrite_n_comp=4; % compartments per dendrite

g_a=1e-7*[0.8150    0.6525    0.3710    0.2109]; %transverse conductances from soma to periphery
g_L=1e-7*[0.1017    0.0433    0.0246    0.0281]; % leakage conductances from soma to periphery
g_NMDA=1e-6*[0.1730    0.0868    0.0436    0.0219]; % conductance of excitatory NMDA connection
g_SOMA=1.6568e-07;

% size of pinwheel map
X_size_gen=300;
Y_size_gen=300;

% load pinwheel map
load(strcat('V1_pinwheels_',num2str(X_size_gen),'_',num2str(Y_size_gen)));

V1_min_space_radius=11;
V1_max_space_radius=round(V1_min_space_radius*sqrt(2)); 

V1_lower_X_margin=0;
V1_lower_Y_margin=0;

V1_X_margin=0;
V1_Y_margin=0;

V1_overload_factor=3;

% *** the following are from "Synaptic basis..." ***

V1_ex_angle_min_std=(40/180)*angular_steps; 
V1_ex_angle_max_std=(50/180)*angular_steps;
V1_in_angle_max_std=(70/180)*angular_steps;
V1_min_centr_ex_to_in_ratio=0.3; 
V1_max_centr_ex_to_in_ratio=0.7;

% ***********************************************

V1_central_radius=1.5; % to choose the cell's preferred orientation according to the position within the orientation map

V1_tot_conns=300; % input connections per cell

V1_conns_per_dendr=V1_tot_conns/dendrite_number;

V1_min_pw_radius=dendrite_n_comp*angular_steps/V1_grad;

V1_pw_to_space_coeff=V1_min_pw_radius/V1_min_space_radius;

V1_branch_sigma=1/5*(V1_min_space_radius+V1_max_space_radius)/2;

dendr_x_coord=[60 120 180 240]; % [microns]

V1_ex_space_fracs=exp(-dendr_x_coord/90);
V1_ex_space_fracs=V1_ex_space_fracs/sum(V1_ex_space_fracs);

V1_in_space_fracs=exp(-dendr_x_coord/90);
V1_in_space_fracs=V1_in_space_fracs/sum(V1_in_space_fracs);

%% ******************** GENERATION *************************

fprintf('Setting and saving V1 cells...\n');

% ******************************************************************
% choose a seed and uncomment to make random generation reproducible
% seed=...
% rng(seed); 
% ******************************************************************


[V1_ex_conn_inds,V1_in_conn_inds,V1_x0,V1_y0,V1_space_radii,V1_ex_angle_std_devs,V1_in_angle_std_devs]=generate_cells (V1_pinwheels,V1_grad,dendrite_number,dendrite_n_comp,V1_tot_conns,V1_conns_per_dendr,V1_min_space_radius,V1_max_space_radius,V1_pw_to_space_coeff,V1_ex_space_fracs,V1_in_space_fracs,V1_branch_sigma,V1_central_radius,V1_min_centr_ex_to_in_ratio,V1_max_centr_ex_to_in_ratio,V1_ex_angle_min_std,V1_ex_angle_max_std,V1_in_angle_max_std,X_size,Y_size,V1_lower_X_margin,V1_lower_Y_margin,V1_X_margin,V1_Y_margin,V1_overload_factor,1,angular_steps);

% V1_ex_conn_inds and V1_in_conn_inds contain the indexes of the excitatory
% and inhibitory connections respectively for all cells, mapped in the
% Gabor filter space, which is X_size x Y_size x (2*angular_steps). The
% 2 factor is because for each direction there are two filters of opposite
% polarities. Index 1 is for dummy, because the pixel in position (1,1,1)
% in the Gabor filter is set to zero by convention. The two variables are
% 3D matrices of size (X_size*Y_size*dendrite_number) x dendrite_n_comp x
% m, with m of arbitrary length. The first index scans the dendritic
% branches. The dendrites belonging to the n-th cell have indexes
% dendrite_number*(n-1)+1:1:dendrite_number*n

V1_cell_number=length(V1_x0);

clear V1_pinwheels

save V1_cell_set 

% calculate effective angular distributions of ex and in conns

[V1_cell_angle_dists_ex,V1_cell_angle_dists_in]=generate_V1_cells_dists (V1_ex_conn_inds,V1_in_conn_inds,dendrite_number,X_size,Y_size,angular_steps);

save V_angle_dists V1_cell_angle_dists_ex V1_cell_angle_dists_in -v7.3



