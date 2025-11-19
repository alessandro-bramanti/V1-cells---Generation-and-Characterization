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

analyzed dendrite by dendrite the end-stopping curve 
relative to cell# cell_ind for a given angle=1, 2,... angular_steps
%}

cell_ind=2;
angle=3;

figure
subplot(dendrite_number+1,1,1);
h=plot(V1_end_stopping_resps{cell_ind}(angle,:)*1000);
set(h,'LineWidth',2,'Color',line_color(angle,:));
M=0;
for k=1:dendrite_number
    M=max([M; squeeze(V1_comp_1_ex_inputs{cell_ind}(k,angle,:)); squeeze(V1_comp_1_in_inputs{cell_ind}(k,angle,1))]);
end
M=M*1.3;
for k=1:dendrite_number
    subplot(dendrite_number+1,1,k+1);
    hold on
    h=plot(squeeze(V1_comp_1_ex_inputs{cell_ind}(k,angle,:)),'r');
    set(h,'LineWidth',2);
    h=plot(squeeze(V1_comp_1_in_inputs{cell_ind}(k,angle,:)),'b');
    set(h,'LineWidth',2);    
    v=axis;
    axis([v(1) v(2) 0 M])
end