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

initialize LGN responses to gratings
%}

I_raws=cell(angular_steps*length(periods)*length(phases),1);
LGN_responses=cell(angular_steps*length(periods)*length(phases),1);

inp_ind=1;

for angle_val=1:angular_steps

    theta=(angle_val-1)/angular_steps*pi;
    x_rot=cos(theta)*x_unrot+sin(theta)*y_unrot;
    y_rot=sin(theta)*x_unrot-cos(theta)*y_unrot;

    for period_ind=1:1:length(periods)
        for phase_ind=1:1:length(phases)

            I_raw=(1/2*cos(2*pi*(y_rot-Y_size/2)/periods(period_ind)+phases(phase_ind))+1/2);

            I_raws{inp_ind}=I_raw;
            
            preliminary_LGN_response;
            
            initial_step=1;

            LGN_response;

            LGN_responses{inp_ind}=I_conv;
                        
            inp_ind=inp_ind+1;
            
        end    
    end
end