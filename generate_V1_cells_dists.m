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

ancillary script

%}

function [V1_cell_angle_dists_ex,V1_cell_angle_dists_in]=generate_V1_cells_dists (V1_ex_conn_inds,V1_in_conn_inds,dendrite_number,X_size,Y_size,angular_steps)


V1_cell_number=size(V1_ex_conn_inds,1)/dendrite_number;

V1_cell_angle_dists_ex=zeros(V1_cell_number,angular_steps);
V1_cell_angle_dists_in=zeros(V1_cell_number,angular_steps);

for i=1:V1_cell_number
    dendr_inds=dendrite_number*(i-1)+1:1:dendrite_number*i;

    conn_inds=V1_ex_conn_inds(dendr_inds,:,:);
    conn_inds=conn_inds(conn_inds(:)>1);

    [~,~,conn_inds]=ind2sub([Y_size X_size 2*angular_steps],conn_inds);
    conn_inds(conn_inds>angular_steps)=conn_inds(conn_inds>angular_steps)-angular_steps;

    dist=histcounts(conn_inds,0.5:1:angular_steps+0.5);
    dist=dist/sqrt(sum(dist.^2));

    V1_cell_angle_dists_ex(i,:)=dist;

    conn_inds=V1_in_conn_inds(dendr_inds,:,:);
    conn_inds=conn_inds(conn_inds(:)>1);

    [~,~,conn_inds]=ind2sub([Y_size X_size 2*angular_steps],conn_inds);
    conn_inds(conn_inds>angular_steps)=conn_inds(conn_inds>angular_steps)-angular_steps;

    dist=histcounts(conn_inds,0.5:1:angular_steps+0.5);
    dist=dist/sqrt(sum(dist.^2));

    V1_cell_angle_dists_in(i,:)=dist;

end

