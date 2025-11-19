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

shows the RF spot-maps for all the cells
annotated with nominal orientation and sigma_rel
as in Figure SF-7 of the paper

iteratively calls show_field_map

%}

close all

x0=0.772;
y0=0.222;

rx=0.03;
ry=0.045;

for i=1:length(test_cell_inds)
    show_field_map(pos_field_map(:,:,i),neg_field_map(:,:,i),0,0.0001);
    axis square
    set(gca,'Color','k')
    title(sprintf('%d %f',(V1_max_angle(i)),(V1_in_angle_std_devs(test_cell_inds(i))/V1_ex_angle_std_devs(test_cell_inds(i)))));
    hold on

    cell_str=strcat('C',num2str(i));
    sigma_rel_str=num2str(V1_in_angle_std_devs(test_cell_inds(i))/V1_ex_angle_std_devs(test_cell_inds(i)),'%.2f');

    annotation('textbox',[0.7, 0.79, 0.1, 0.1],'String',sigma_rel_str,'FontSize',14,'BackgroundColor','black','EdgeColor','none','Color','white','FontWeight','bold');
    angle=(V1_max_angle(i)-1)*pi/angular_steps;
    annotation('line',[x0-rx*cos(angle), x0+rx*cos(angle)],[y0-ry*sin(angle), y0+ry*sin(angle)],'Color','white','LineWidth',5)

end