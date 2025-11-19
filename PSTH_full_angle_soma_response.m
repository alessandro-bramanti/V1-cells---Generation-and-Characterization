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

map response of cells to gratings
%}

V1_soma_resp_map=zeros(cell_number,length(periods),length(phases),angular_steps);

inp_ind=1;

t0=tic;

for angle_val=1:angular_steps
    for period_ind=1:1:length(periods)
        for phase_ind=1:1:length(phases)

             V1_soma_resp_map(:,period_ind,phase_ind,angle_val)=V1_soma_response(cell_number,V1_ex_conn_inds(test_dendr_inds,:,:),V1_in_conn_inds(test_dendr_inds,:,:),LGN_responses{inp_ind});    

            inp_ind=inp_ind+1;
            
        end
    end
    
    fprintf('Examined %d angular orientations out of %d; %d minutes elapsed.\n',angle_val,angular_steps,round(toc(t0)/60));
        
end

