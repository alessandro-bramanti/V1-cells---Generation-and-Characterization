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

spot-maps the receptive field of the cells

%}

pos_field_map=zeros(Y_size,X_size,cell_number);
neg_field_map=zeros(Y_size,X_size,cell_number);

spot_radius=2;

tot_pixels=(X_size-2*spot_radius)*(Y_size-2*spot_radius);

pix_ind=0;
part_pix_ind=0;

t0=tic;

for y_ind=spot_radius+1:Y_size-spot_radius-1
    shift_ind_y=-Y_size/2+y_ind;
    for x_ind=spot_radius+1:X_size-spot_radius-1
        
        shift_ind_x=-X_size/2+x_ind;
        
        I_conv=circshift(pos_spot_resp_0,[shift_ind_y shift_ind_x 0]);   
        
        if X_size/2+shift_ind_x<=2*spot_radius
            I_conv(:,X_size-spot_radius-1:1:X_size,:)=0;
        elseif X_size/2+shift_ind_x>=X_size-2*spot_radius
            I_conv(:,1:1:spot_radius+1,:)=0;
        end
        
        if Y_size/2+shift_ind_y<=2*spot_radius
            I_conv(Y_size-spot_radius-1:1:Y_size,:,:)=0;
        elseif Y_size/2+shift_ind_y>=Y_size-2*spot_radius
            I_conv(1:1:spot_radius+1,:,:)=0;
        end
        
        V1_soma_resp_map_pos=V1_soma_response(cell_number,V1_ex_conn_inds(test_dendr_inds,:,:),V1_in_conn_inds(test_dendr_inds,:,:),I_conv);  

        for i=1:cell_number
            pos_field_map(y_ind,x_ind,i)=pos_field_map(y_ind,x_ind,i)+V1_soma_resp_map_pos(i);        
        end
     
        I_conv=circshift(neg_spot_resp_0,[shift_ind_y shift_ind_x 0]);                    
        
        if X_size/2+shift_ind_x<=2*spot_radius
            I_conv(:,X_size-spot_radius-1:1:X_size,:)=0;
        elseif X_size/2+shift_ind_x>=X_size-2*spot_radius
            I_conv(:,1:1:spot_radius+1,:)=0;
        end
        
        if Y_size/2+shift_ind_y<=2*spot_radius
            I_conv(Y_size-spot_radius-1:1:Y_size,:,:)=0;
        elseif Y_size/2+shift_ind_y>=Y_size-2*spot_radius
            I_conv(1:1:spot_radius+1,:,:)=0;
        end
        
        V1_soma_resp_map_neg=V1_soma_response(cell_number,V1_ex_conn_inds(test_dendr_inds,:,:),V1_in_conn_inds(test_dendr_inds,:,:),I_conv);    
        
        for i=1:cell_number
            neg_field_map(y_ind,x_ind,i)=neg_field_map(y_ind,x_ind,i)+V1_soma_resp_map_neg(i);
        end
        
        pix_ind=pix_ind+1;
        part_pix_ind=part_pix_ind+1;
        
        if part_pix_ind>=0.1*tot_pixels
            part_pix_ind=0;
            fprintf('Completed %d%% of mapping, time elapsed: %d minutes.\n',round(pix_ind/tot_pixels*100),round(toc(t0)/60));
        end
    end
end

RF_areas=zeros(1,length(test_cell_inds));
ov_indexes=zeros(1,length(test_cell_inds));
mean_pos_vals=zeros(1,length(test_cell_inds));
max_pos_vals=zeros(1,length(test_cell_inds));

RF_x_center=zeros(1,length(test_cell_inds));
RF_y_center=zeros(1,length(test_cell_inds));
RF_radius=zeros(1,length(test_cell_inds));

pos_field_map=pos_field_map*correction_factors(i);
neg_field_map=neg_field_map*correction_factors(i);

set(0,'DefaultFigureVisible','on')





