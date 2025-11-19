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

This scripts characterizes the generated cell set throughç
    - gratings and calculation of correction coefficients
    - spot-mapping
    - end-stopping
%}

%***** how many cells do you want to pick up? ****

cell_number=100;

%***** load cell set *****

cell_set_filename='V1_cell_set.mat';
load(cell_set_filename);

%***** set file for saving simulation *****

save_filename='V1_characterization.mat';

%**** choose whether to pick up the same cells as in the paper ***********
% (works if the generated set is the same and cell_number=100)

paper_set=1; % if 0 new cells will be picked up even from the same set

%***** parameters *******

% avoid picking up peripheral cells
X_margin=30;
Y_margin=30;

% gratings
phase_step=2*pi/20;
min_period=6;
max_period=50;
period_step=2;

periods=min_period:period_step:max_period;
phases=0:phase_step:2*pi-phase_step;

[y_unrot,x_unrot]=ndgrid(1:1:Y_size,1:1:X_size);

%***** pick up cells to be characterized *****

if paper_set
    test_cell_inds=[62870 32573 36177 50794 38495 32305 ...
        55704 37356 31309 49681 47812 53889 49664 27716 43967 ...
        29853 53520 54693 29309 57485 44208 29578 57141 38879 31108 ...
        55703 48128 44896 39469 30731 61077 46920 58401 59592 30448 ...
        34387 43318 46050 40979 28933 58397 27763 40073 34940 46362 ...
        62814 36126 35525 29217 48771 62843 47251 46903 62205 56233 ...
        32223 28689 62240 27193 37968 60107 38553 60117 48431 31086 ...
        36460 62850 28105 55594 50897 55050 28616 61398 58397 40042 ...
        51157 35211 43605 46908 38191 46303 43105 46053 30107 38238 ...
        56269 40340 51804 45442 42149 42438 59285 38260 49317 60714 ...
        32557 29302 56846 31102 58041];
else
    cells_left=cell_number; %#ok<UNRCH> 
    test_cell_inds=[];
    while cells_left>0
        new_test_cell_inds=ceil(rand(1,cells_left)*length(V1_x0));
        new_test_cell_inds=new_test_cell_inds(V1_x0(new_test_cell_inds)>=min_x0 & V1_x0(new_test_cell_inds)<=max_x0 & ...
                                              V1_y0(new_test_cell_inds)>=min_y0 & V1_y0(new_test_cell_inds)<=max_y0);
        test_cell_inds=[test_cell_inds new_test_cell_inds]; %#ok<AGROW>
        cells_left=cells_left-length(new_test_cell_inds);
    end
end

test_dendr_inds=zeros(dendrite_number*cell_number,1);
for i=1:cell_number
    test_dendr_inds((i-1)*dendrite_number+1:1:(i-1)*dendrite_number+dendrite_number,1)=[(test_cell_inds(i)-1)*dendrite_number+1:1:(test_cell_inds(i)-1)*dendrite_number+dendrite_number]';
end

%************ CHARACTERIZATION ********************


%% *********** PSTH ********************
% calculate PSTH based on gratings
% calculate correction factors

fprintf('Initializing grating input stimuli...\n');
init_PSTH_inputs;
fprintf('done.\n');

fprintf('Computing angular responses (contrast: 100%%)...\n');
PSTH_full_angle_soma_response;

% calculate correction factors
calculate_correction_factors_PSTH_grid_mapping    
for i=1:cell_number
    V1_soma_resp_map(i,:,:,:)=V1_soma_resp_map(i,:,:,:)*correction_factors(i); %#ok<SAGROW> 
end

% V1_soma_resp_map is a 4D matrix of dimensions:
% cell_number x length(periods) x length(phases) x angular_steps
% squeeze(V1_soma_resp_map(n,:,:,ang) is the response of cell n
% to all the gratings oriented along angle ang 
% ang=1:1:angular_steps; in the paper, angular_steps=10
% meaning angles (in degrees) 0, 18, 36... 162.
% somatic potentials are in V 
% the results can be represented as in Figure SF-8 of the paper except that
% no thresholds are applied automatically

save(save_filename,"V1_soma_resp_map","correction_factors");

%% *********** spot-mapping *******************
% spot mapping
% mOSInom and mOSIeff

fprintf('Spotwise mapping of receptive fields...\n');
init_spotwise_input;
PSTH_spot_map;  

% calculates effective RF's radii
eff_r=zeros(1,cell_number);
for i=1:cell_number
    p=pos_field_map(:,:,i);
    n=neg_field_map(:,:,i);
    s=sum(p(:)>0.01*max(p(:)) | n(:)>0.01*max(n(:)));
    eff_r(i)=sqrt(s/pi)/2;
end

% retrieves nominal maximum orientation
cell_angular_dists_filename='V1_angular_dists.mat';
load(cell_angular_dists_filename);
V1_max_angle=zeros(cell_number,1);
for i=1:cell_number
    [~,V1_max_angle(i)]=max(V1_cell_angle_dists_ex(test_cell_inds(i),:));
end
clear V1_cell_angle_dists_ex V1_cell_angle_dists_in

% each RF spot map is composed of two maps obtained with the positive and
% negative spot respectively, pos_field_map and neg_field_map
% both are Y_size x X_size x cell_number in size
% they are already corrected (correction factors) and are in V
% no thresholds are applied automatically
% show_spot_mapping_results shows all maps annotated, as explained inside
% use show_field_map for individual RF maps

show_spot_mapping_results

% calculate mOSInom and mOSIeff
calculate_mOSI_nom_eff

save(save_filename,"pos_field_map","neg_field_map","mOSInom","mOSIeff","-append");

%% *********** spot-mapping *******************
% WARNING: both previous sections must have been run in advance 
% to calculate the effective radii and the mOSI coefficients

fprintf('End-stopping characterization...\n');
end_stopping_complete_set_test

% display the results like in Figure SF-7 of the paper
show_end_stopping_test_results

save(save_filename,"V1_end_stopping_resps","-append");


