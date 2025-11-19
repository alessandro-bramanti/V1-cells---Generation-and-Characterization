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

tests all the cells for end-stopping effect
%}

close all

% draw arrow legenda

V1_end_stopping_resps=cell(1,length(test_cell_inds));
V1_end_stopping_dendr_resps=cell(1,length(test_cell_inds));

cell_count=1;                

part_cell_count=1;  

t0=tic;

for i=1:length(test_cell_inds)

    dendr_inds=dendrite_number*(test_cell_inds(i)-1)+1:1:dendrite_number*test_cell_inds(i);

    % initialize
    
    cell_ind=test_cell_inds(i);
    
    ex_conns=V1_ex_conn_inds(dendr_inds,:,:);
    [ex_y,ex_x,ex_angle]=ind2sub([Y_size X_size 2*angular_steps],ex_conns);
    ex_angle(ex_angle>angular_steps)=ex_angle(ex_angle>angular_steps)-angular_steps;
    
    in_conns=V1_in_conn_inds(dendr_inds,:,:);
    [in_y,in_x,in_angle]=ind2sub([Y_size X_size 2*angular_steps],in_conns);
    in_angle(in_angle>angular_steps)=in_angle(in_angle>angular_steps)-angular_steps;
    
    % check end-stopping properties
    
    M=zeros(1,angular_steps);
    M_ind=zeros(1,angular_steps);
    
    for k=1:angular_steps
        [M(k)]=sum(sum(squeeze(V1_soma_resp_map(i,:,:,k))));
    end
    
    single_test_dendr_inds=dendrite_number*(test_cell_inds(i)-1)+1:1:dendrite_number*test_cell_inds(i);
    
    contrast=1;
    r_steps=20;


    V1_end_stopping_resps{i}=zeros(2*angular_steps,2*r_steps-1);
    V1_end_stopping_dendr_resps{i}=zeros(dendrite_number,2*angular_steps,2*r_steps-1);

    bar_width=eff_r(i);
    
    for angle=1:2*angular_steps
    
        if angle>angular_steps
            angle_val=angle-angular_steps;
            up_to_down_or_right_to_left=1;
        else
            angle_val=angle;
            up_to_down_or_right_to_left=0;
        end
    
        theta=(angle_val-1)/angular_steps*pi;
        
        if theta~=pi/2 
            m=tan(theta);
            q1=V1_y0(test_cell_inds(i))+bar_width/(2*cos(theta))-m*V1_x0(test_cell_inds(i));
            q2=V1_y0(test_cell_inds(i))-bar_width/(2*cos(theta))-m*V1_x0(test_cell_inds(i));
            full_input=zeros(size(x_unrot))+contrast*((x_unrot-V1_x0(test_cell_inds(i))).^2+((y_unrot-V1_y0(test_cell_inds(i))).^2)<=eff_r(i)^2).*(y_unrot<=m*x_unrot+q1 & y_unrot>=m*x_unrot+q2);
        else
            q1=V1_x0(test_cell_inds(i))+bar_width/(2*cos(theta));
            q2=V1_x0(test_cell_inds(i))-bar_width/(2*cos(theta));
            full_input=zeros(size(x_unrot))+contrast*((x_unrot-V1_x0(test_cell_inds(i))).^2+((y_unrot-V1_y0(test_cell_inds(i))).^2)<=eff_r(i)^2).*(x_unrot<=q1 & y_unrot>=q2);
        end
    
       
        if up_to_down_or_right_to_left
            rs=eff_r(i):-eff_r(i)/(r_steps-1):-eff_r(i); %#ok<*IJCL> 
        else
            rs=-eff_r(i):eff_r(i)/(r_steps-1):eff_r(i); %#ok<*IJCL>
        end
    
        saf_coeff=1.2;
          
        if theta~=0 && theta~=pi/2
            m_norm=-1/m;
            if up_to_down_or_right_to_left
                q_lim=(V1_y0(test_cell_inds(i))+saf_coeff*eff_r*sin(theta))-m_norm*(V1_x0(test_cell_inds(i))+saf_coeff*eff_r*cos(theta));
            else
                q_lim=(V1_y0(test_cell_inds(i))-saf_coeff*eff_r*sin(theta))-m_norm*(V1_x0(test_cell_inds(i))-saf_coeff*eff_r*cos(theta));
            end
            ind=1;
            for r_val=rs
                I_raw=full_input;
                q_lim_2=(V1_y0(test_cell_inds(i))+r_val*sin(theta))-m_norm*(V1_x0(test_cell_inds(i))+r_val*cos(theta));
                if up_to_down_or_right_to_left
                    I_raw=I_raw.*(y_unrot>=m_norm*x_unrot+q_lim_2 & y_unrot<=m_norm*x_unrot+q_lim);
                else
                    I_raw=I_raw.*(y_unrot<=m_norm*x_unrot+q_lim_2 & y_unrot>=m_norm*x_unrot+q_lim);
                end
                preliminary_LGN_response;
                attention_focused=1;
                initial_step=1;
                LGN_response;
                
                [V1_end_stopping_resps{i}(angle,ind),V1_end_stopping_dendr_resps{i}(:,angle,ind)]=V1_soma_response(1,V1_ex_conn_inds(single_test_dendr_inds,:,:),V1_in_conn_inds(single_test_dendr_inds,:,:),I_conv);    
                V1_end_stopping_resps{i}(angle,ind)=V1_end_stopping_resps{i}(angle,ind)*correction_factors(i);
                V1_end_stopping_dendr_resps{i}(:,angle,ind)=V1_end_stopping_dendr_resps{i}(:,angle,ind)*correction_factors(i);
                ind=ind+1;
            end
        elseif theta==pi/2
            if up_to_down_or_right_to_left
                q_lim=(V1_x0(test_cell_inds(i))+saf_coeff*eff_r(i)); %#ok<*UNRCH> 
            else
                q_lim=(V1_x0(test_cell_inds(i))-saf_coeff*eff_r(i));
            end
            ind=1;
            for r_val=rs
                I_raw=full_input;
                q_lim_2=V1_x0(test_cell_inds(i))+r_val;
                if up_to_down_or_right_to_left
                    I_raw=I_raw.*(x_unrot>=q_lim_2 & x_unrot<=q_lim);
                else
                    I_raw=I_raw.*(x_unrot<=q_lim_2 & x_unrot>=q_lim);
                end
                preliminary_LGN_response;
                attention_focused=1;
                initial_step=1;
                LGN_response;
                
                [V1_end_stopping_resps{i}(angle,ind),V1_end_stopping_dendr_resps{i}(:,angle,ind)]=V1_soma_response(1,V1_ex_conn_inds(single_test_dendr_inds,:,:),V1_in_conn_inds(single_test_dendr_inds,:,:),I_conv);    
                V1_end_stopping_resps{i}(angle,ind)=V1_end_stopping_resps{i}(angle,ind)*correction_factors(i);
                V1_end_stopping_dendr_resps{i}(:,angle,ind)=V1_end_stopping_dendr_resps{i}(:,angle,ind)*correction_factors(i);
                ind=ind+1;
            end
        else % theta=0
            if up_to_down_or_right_to_left
                q_lim=(V1_x0(test_cell_inds(i))+saf_coeff*eff_r(i)); %#ok<*UNRCH> 
            else
                q_lim=(V1_x0(test_cell_inds(i))-saf_coeff*eff_r(i));
            end
            ind=1;
            for r_val=rs
                I_raw=full_input;
                q_lim_2=V1_x0(test_cell_inds(i))+r_val;
                if up_to_down_or_right_to_left
                    I_raw=I_raw.*(x_unrot>=q_lim_2 & x_unrot<=q_lim);
                else
                    I_raw=I_raw.*(x_unrot<=q_lim_2 & x_unrot>=q_lim);
                end
                preliminary_LGN_response;
                attention_focused=1;
                initial_step=1;
                LGN_response;
                
                [V1_end_stopping_resps{i}(angle,ind),V1_end_stopping_dendr_resps{i}(:,angle,ind)]=V1_soma_response(1,V1_ex_conn_inds(single_test_dendr_inds,:,:),V1_in_conn_inds(single_test_dendr_inds,:,:),I_conv);    
                V1_end_stopping_resps{i}(angle,ind)=V1_end_stopping_resps{i}(angle,ind)*correction_factors(i);
                V1_end_stopping_dendr_resps{i}(:,angle,ind)=V1_end_stopping_dendr_resps{i}(:,angle,ind)*correction_factors(i);
                ind=ind+1; 
            end
    
        end
    end

    cell_count=cell_count+1;
    part_cell_count=part_cell_count+1;
    
    if part_cell_count>=0.1*cell_number
        part_cell_count=0;
        fprintf('Completed %d%% of characterization, time elapsed: %d minutes.\n',round(cell_count/cell_number*100),round(toc(t0)/60));
    end


end
