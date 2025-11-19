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

represents the results of the end-stopping test as in Figure SF-7 of the paper

%}

cm=colormap;
close all

line_color=zeros(2*angular_steps,3);
line_color(1,:)=cm(1,:);
for k=2:2*angular_steps
    cm_pos=round(size(cm,1)/(2*angular_steps-1)*(k-1));
    line_color(k,:)=cm(cm_pos,:);
end

char_size=25;
x_axis_size=size(V1_end_stopping_resps{1},2);
x_margin=round(0.05*x_axis_size);

for i=1:length(test_cell_inds)
    figure;
    hold on;
    title(strcat('mOSI_{nom} ',num2str(mOSInom(i)),' \sigma_{rel} ',num2str(V1_in_angle_std_devs(test_cell_inds(i))/V1_ex_angle_std_devs(test_cell_inds(i)))))
    for angle=1:2*angular_steps
        if max(V1_end_stopping_resps{i}(angle,:))>0.015
            h=plot(V1_end_stopping_resps{i}(angle,:)*1000);
            set(h,'LineWidth',5,'Color',line_color(angle,:));
        end
    end
    h=plot([1:size(V1_end_stopping_resps{i},2)],15*ones(1,size(V1_end_stopping_resps{i},2)),'r--');
    set(h,'LineWidth',5);

    mOSI_str=num2str(mOSInom(i),'%.2f');

    annotation('textbox',[0.15, 0.125, 0.1, 0.1],'String',mOSI_str,'FontSize',20,'BackgroundColor','white','EdgeColor','none','Color','red','FontWeight','bold')
    annotation('rectangle',[0.685, 0.125, 0.2 0.25],'EdgeColor','black','FaceColor','white')

    M=max(V1_end_stopping_resps{i}(:));

    if max(V1_end_stopping_resps{i}(1,:))>0.015
        f=max(V1_end_stopping_resps{i}(1,:))/M;
        annotation('arrow',[0.785 0.785+0.08*f],[0.25 0.25],'LineWidth',3,'Color',line_color(1,:));
    end

    for k=2:2*angular_steps
        if max(V1_end_stopping_resps{i}(k,:))>0.015
            ang=(k-1);
            cm_pos=round(size(cm,1)/(2*angular_steps-1)*(k-1));
            line_color(k,:)=cm(cm_pos,:);
            f=max(V1_end_stopping_resps{i}(k,:))/M;
            annotation('arrow',[0.785 0.785+0.08*f*cos((k-1)/angular_steps*pi)],[0.25 0.25+0.1*f*sin((k-1)/angular_steps*pi)],'LineWidth',3,'Color',line_color(k,:));
        end
    end


end

figure
hold on
h=quiver([0 1],[0 0]);
set(h,'LineWidth',5,'Color',line_color(1,:));

for k=2:2*angular_steps
    ang=(k-1);
    cm_pos=round(size(cm,1)/(2*angular_steps-1)*(k-1));
    h=quiver([0 cos((k-1)/angular_steps*pi)],[0 sin((k-1)/angular_steps*pi)]);
    set(h,'LineWidth',5,'Color',line_color(k,:));
end
