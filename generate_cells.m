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

function [ex_conn_inds,in_conn_inds,x0,y0,cell_space_radii,ex_angle_std_devs,in_angle_std_devs]=generate_cells (pinwheels,pw_grad,dendrite_number,dendrite_n_comp,tot_conns,conns_per_dendr,min_space_radius,max_space_radius,pw_to_space_coeff,ex_space_fracs,in_space_fracs,branch_sigma,central_radius,min_centr_ex_to_in_ratio,max_centr_ex_to_in_ratio,ex_angle_min_std,ex_angle_max_std,in_angle_max_std,X_size,Y_size,lower_X_margin,lower_Y_margin,X_margin,Y_margin,overload_factor,lower_overload_factor,angular_steps)


% coordinates in the pinwheel map space

[pw_y,pw_x]=ndgrid(1:size(pinwheels,1),1:size(pinwheels,2));

[y0,x0]=ndgrid(Y_margin+1/overload_factor:1/overload_factor:Y_size-Y_margin,X_margin+1/overload_factor:1/overload_factor:X_size-X_margin);

% central coordinates of the cells (within the previous, possibly with
% fractional values in case overload_factor>1

x0=reshape(x0,[],1);
y0=reshape(y0,[],1);

cell_number=length(x0);

cell_space_radii=zeros(cell_number,1);

ex_angle_std_devs=zeros(cell_number,1);
in_angle_std_devs=zeros(cell_number,1);

% ex and in connection inds per cell (each cell is an element of the cell
% array, because each cell has different numbers of connections

cell_ex_inds=cell(cell_number,1);
cell_in_inds=cell(cell_number,1);

tot_ex_comp_length=0;
tot_in_comp_length=0;

% the minimum space radius is calculated for subsequent use

min_angle_radius=min_space_radius*pw_to_space_coeff;

pw_X_size=size(pinwheels,2);
pw_Y_size=size(pinwheels,1);

pw_numel=pw_X_size*pw_Y_size;

central_x=round(X_size/2);
central_y=round(Y_size/2);

pw_central_x=round(pw_X_size/2);
pw_central_y=round(pw_Y_size/2);

pw_central_inds_0=find((pw_x-pw_central_x).^2+(pw_y-pw_central_y).^2<=central_radius^2);

for n_cell=1:cell_number % main loop: cell-by-cell synthesis
        
    % choose the space radius (in the visual coordinates) and ex and in
    % connection numbers subsequently
    
    space_radius=min_space_radius+(max_space_radius-min_space_radius)*rand(1);
    
    cell_space_radii(n_cell)=space_radius;
    
    % central coordinate of the current cell in the pinwheel system

    
    pw_x0=pw_central_x+round((x0(n_cell)-central_x)/(pw_grad)*pw_to_space_coeff);
    pw_x0=min(pw_x0,pw_X_size);
    pw_y0=pw_central_y+round((y0(n_cell)-central_y)/(pw_grad)*pw_to_space_coeff);
    pw_y0=min(pw_y0,pw_Y_size);

    % indexes (in the pinwheel map space) of the central area, shaping the
    % cell's response
    
    central_ex_to_in_ratio=min_centr_ex_to_in_ratio+(max_centr_ex_to_in_ratio-min_centr_ex_to_in_ratio)*rand(1);
    ex_std=ex_angle_min_std+rand(1)*(ex_angle_max_std-ex_angle_min_std); % standard deviation (on average 45 degress for V1 cells, 0.25 for upper-than-V1

    a=1;                           % integral of n_ex connections at a minimum radius up to a multiplicative constant (so fixed at 1 by convention)
    b=1/central_ex_to_in_ratio*a*in_angle_max_std/ex_std; % (a/sigma_ex)/(b/sigma_in)=central_ex_to_in_ratio, ex-to-in distribution rate at the peak
                                       % small cells tend to be simple,
                                       % and so, by assumption, have
                                       % maximum inhib. std. dev.
    
    ex_conn_num_min=tot_conns*a/(a+b);    % # of ex connections at the smallest radius
    
    in_conn_num_max=tot_conns-ex_conn_num_min;
    
    angle_radius=space_radius*pw_to_space_coeff;
    
    ex_conn_num=round(ex_conn_num_min*(angle_radius/min_angle_radius)^2); % # of ex connections at the true radius
    a_large=a*(space_radius/min_space_radius)^2; % integral of ex conn distribution at true radius up to a multiplicative constant
    in_conn_num=tot_conns-ex_conn_num;    % # of inhibitory connections
    b_large=b*in_conn_num/in_conn_num_max; % b at real radius, considering that b/(b+a)*n_tot=n_in and that tot_conns-in_conn_num=ex_conn_num
    in_std=central_ex_to_in_ratio*ex_std*b_large/a_large; % standard deviation of in connection at large radius, recalling that (a/sigma_ex)/(b/sigma_in)=0.5

    % V1 cells inherit the orientations in the layer below
        
    pw_central_inds=pw_central_inds_0+(pw_y0-pw_central_y)+pw_Y_size*(pw_x0-pw_central_x);
    
    pw_central_inds=pw_central_inds(pw_central_inds>0 & pw_central_inds<=pw_numel);
    pw_central_inds=pw_central_inds(pw_central_inds>0 & pw_central_inds<=pw_numel);
    
    % cut out the relevant part of the pinwheel map

    loc_pw_x=round(pw_x0-angle_radius);
    loc_pw_y=round(pw_y0-angle_radius);

    pw_x2=round(pw_x0+angle_radius);
    pw_y2=round(pw_y0+angle_radius);

    loc_pw_x=max(loc_pw_x,1);
    loc_pw_x=min(loc_pw_x,size(pinwheels,2));        

    pw_x2=max(pw_x2,1);
    pw_x2=min(pw_x2,size(pinwheels,2));        

    loc_pw_y=max(loc_pw_y,1);
    loc_pw_y=min(loc_pw_y,size(pinwheels,1));        

    pw_y2=max(pw_y2,1);
    pw_y2=min(pw_y2,size(pinwheels,1));    

    local_pw=pinwheels(loc_pw_y:pw_y2,loc_pw_x:pw_x2);

    % compute map coordinates

    loc_pw_x=pw_x(loc_pw_y:pw_y2,loc_pw_x:pw_x2);
    loc_pw_y=pw_y(loc_pw_y:pw_y2,loc_pw_x:pw_x2);       
    
    central_pw_set=pinwheels(pw_central_inds); % central connections of V1 cell, shaping its response
    
    central_angle=mode(central_pw_set(:));


    ex_angle_std_devs(n_cell)=ex_std;    
    in_angle_std_devs(n_cell)=in_std;      

    
    % refer the map to a fictitious central "zero" angle

    local_pw=local_pw-central_angle;
    local_pw(local_pw>angular_steps/2)=local_pw(local_pw>angular_steps/2)-angular_steps;    
    local_pw(local_pw<=-angular_steps/2)=local_pw(local_pw<=-angular_steps/2)+angular_steps;        

    % compute the distribution of ex and in connections centred on
    % fictitious 0, partitioned between spatial areas
    % consider one sub-distribution per angle within the central connection
    % radius

    ex_conn_dist=round(randn(1,ex_conn_num)*ex_std);
    in_conn_dist=round(randn(1,in_conn_num)*in_std);
            
    ex_conn_dist(ex_conn_dist<=-angular_steps/2)=ex_conn_dist(ex_conn_dist<=-angular_steps/2)+angular_steps;
    ex_conn_dist(ex_conn_dist>angular_steps/2)=ex_conn_dist(ex_conn_dist>angular_steps/2)-angular_steps;

    in_conn_dist(in_conn_dist<=-angular_steps/2)=in_conn_dist(in_conn_dist<=-angular_steps/2)+angular_steps;
    in_conn_dist(in_conn_dist>angular_steps/2)=in_conn_dist(in_conn_dist>angular_steps/2)-angular_steps;

    conn_vals=-angular_steps/2+1:1:angular_steps/2;
    
    ex_conn_freqs=histcounts(ex_conn_dist,[conn_vals-0.5 max(conn_vals)+0.5]);
    in_conn_freqs=histcounts(in_conn_dist,[conn_vals-0.5 max(conn_vals)+0.5]);
 
        
    % find angle values included in the pinwheel submap cut out (in
    % principle all values will be included)
    
    % shifted_angle_values=unique(local_pw(~isnan(local_pw)));

    % *_angle_ex and *_angle_in contain the x and y coordinates of the ex
    % and in connections respectively, in the pinwheel space
    
    x_angle_ex=zeros(1,ex_conn_num);
    y_angle_ex=zeros(1,ex_conn_num);

    angle_ex_ind=1;
    
    x_angle_in=zeros(1,in_conn_num);
    y_angle_in=zeros(1,in_conn_num);

    angle_in_ind=1;
    
    layer_angle_ex=zeros(1,ex_conn_num);
    layer_angle_in=zeros(1,in_conn_num);    
     
    conn_angle_ex=zeros(length(conn_vals),length(ex_space_fracs));
    conn_angle_in=zeros(length(conn_vals),length(in_space_fracs));
    
  
    for j=1:length(conn_vals)

        conn_angle_ex(j,:)=round(ex_conn_freqs(j)*ex_space_fracs);
        if sum(conn_angle_ex(j,:))<ex_conn_freqs(j)
            conn_angle_ex(j,1)=conn_angle_ex(j,1)+(ex_conn_freqs(j)-sum(conn_angle_ex(j,:)));
        else
            j2=length(ex_space_fracs);
            while (sum(conn_angle_ex(j,:))>ex_conn_freqs(j))
                conn_angle_ex(j,j2)=max([conn_angle_ex(j,j2)-1 0]);
                j2=j2-1;
                if j2<1
                    j2=length(ex_space_fracs);
                end
            end
        end

        conn_angle_in(j,:)=round(in_conn_freqs(j)*in_space_fracs);
        if sum(conn_angle_in(j,:))<in_conn_freqs(j)
            conn_angle_in(j,1)=conn_angle_in(j,1)+(in_conn_freqs(j)-sum(conn_angle_in(j,:)));
        else
            j2=length(in_space_fracs);
            while (sum(conn_angle_in(j,:))>in_conn_freqs(j))
                conn_angle_in(j,j2)=max([conn_angle_in(j,j2)-1 0]);
                j2=j2-1;
                if j2<1
                    j2=length(in_space_fracs);
                end
            end
        end

    end    
    
    
    % for each of the angle values in the submap choose a correct number of
    % ex and in connections corresponding to the value itself and the space
    % region 

    for j=1:length(ex_space_fracs)

        for k=1:1:length(conn_vals)

            ang_inds=find(local_pw==conn_vals(k) & ((loc_pw_x-pw_x0).^2+(loc_pw_y-pw_y0).^2>=((j-1)/length(ex_space_fracs)*angle_radius)^2 & (loc_pw_x-pw_x0).^2+(loc_pw_y-pw_y0).^2<=(j/length(ex_space_fracs)*angle_radius)^2));

            if ~isempty(ang_inds)

                ang_space_inds=ceil(rand(1,conn_angle_ex(k,j))*length(ang_inds));
                x_angle_ex(1,angle_ex_ind:1:angle_ex_ind+length(ang_space_inds)-1)=loc_pw_x(ang_inds(ang_space_inds));
                y_angle_ex(1,angle_ex_ind:1:angle_ex_ind+length(ang_space_inds)-1)=loc_pw_y(ang_inds(ang_space_inds));
                angle_ex_ind=angle_ex_ind+length(ang_space_inds);

                ang_space_inds=ceil(rand(1,conn_angle_in(k,j))*length(ang_inds));
                x_angle_in(1,angle_in_ind:1:angle_in_ind+length(ang_space_inds)-1)=loc_pw_x(ang_inds(ang_space_inds));
                y_angle_in(1,angle_in_ind:1:angle_in_ind+length(ang_space_inds)-1)=loc_pw_y(ang_inds(ang_space_inds));
                angle_in_ind=angle_in_ind+length(ang_space_inds);

            end

        end
        
    end

    % renormalize ex and in connections with central angle value

    ex_conn_dist=ex_conn_dist+central_angle;
    ex_conn_dist(ex_conn_dist>angular_steps)=ex_conn_dist(ex_conn_dist>angular_steps)-angular_steps;
    ex_conn_dist(ex_conn_dist<1)=ex_conn_dist(ex_conn_dist<1)+angular_steps;

    in_conn_dist=in_conn_dist+central_angle;
    in_conn_dist(in_conn_dist>angular_steps)=in_conn_dist(in_conn_dist>angular_steps)-angular_steps;
    in_conn_dist(in_conn_dist<1)=in_conn_dist(in_conn_dist<1)+angular_steps;
    

    % find the angle formed by the coordinates in the pinwheel space for
    % each ex or in connection - NOT correlated to the input angle value

    ex_space_angles=atan2(y_angle_ex-pw_y0,x_angle_ex-pw_x0)+pi;
    in_space_angles=atan2(y_angle_in-pw_y0,x_angle_in-pw_x0)+pi;

    % coordinate angles in pw 

    space_angles=[ex_space_angles in_space_angles];

    % corresponding input angles

    space_angle_conn_angle=[ex_conn_dist in_conn_dist];

    % connection type: 1 = ex, 0 = in

    space_angle_type=[ones(size(ex_space_angles)), zeros(size(in_space_angles))];

    % corresponding x and y coordinates

    x_angle=[x_angle_ex x_angle_in];
    y_angle=[y_angle_ex y_angle_in];
    
    layer_angle=[layer_angle_ex layer_angle_in];
    
    % distribute the connections among dendrites: choose a starting angle
    % in pw coordinates, at random

    start_angle=2*pi*rand(1);
    space_angles(space_angles<start_angle)=space_angles(space_angles<start_angle)+2*pi;
    space_angles=space_angles-start_angle;

    % sort the space angles, and the other vectors relatedly: 

    [space_angles,inds]=sort(space_angles);
    space_angle_type=space_angle_type(inds);
    space_angle_conn_angle=space_angle_conn_angle(inds);
    x_angle=x_angle(inds);
    y_angle=y_angle(inds);
    layer_angle=layer_angle(inds);
    
    start_angle_ind=1;

    ex_comp_length=0;
    in_comp_length=0;



    for dendr_num=1:1:dendrite_number % dendrite cycle
        
        % connections with space angles in pw coordinates between
        % start_angle_ind and end_angle_ind affer to the dendr_num-th
        % dendrite

        end_angle_ind=min(start_angle_ind+conns_per_dendr-1,length(x_angle));

        angle_inds=start_angle_ind:1:end_angle_ind;

        % isolate the relevant portions of the vectors above
        
        x_dendr_angle=x_angle(angle_inds);
        y_dendr_angle=y_angle(angle_inds);    
        dendr_conn_angle=space_angle_conn_angle(angle_inds);
        dendr_conn_type=space_angle_type(angle_inds);
        dendr_space_angles=space_angles(angle_inds);
        pw_dendr_radii=sqrt((x_dendr_angle-pw_x0).^2+(y_dendr_angle-pw_y0).^2);
        
        % if the cell is in V1 it is fed by inputs from a packet of
        % layers below (2*angular_step Gaussian filters,
        % corresponding to the angles with both orientations of pos
        % and neg lobes)
        
        dendr_conn_angle=dendr_conn_angle+angular_steps*(rand(1,length(dendr_conn_angle))>0.5);   
        
  %      [pw_radii,sort_inds]=sort(pw_radii);
  %      dendr_conn_angle=dendr_conn_angle(sort_inds);
  %      dendr_type=dendr_type(sort_inds);
  %      dendr_space_angles=dendr_space_angles(sort_inds);
        
        for comp_num=1:dendrite_n_comp % dendrite compartment cycle
        

            % ** calculate inds of ex connections for each dendrite
        
            % convert ex connection coordinates from pinwheel space to visual
            % coordinates, compartment by compartment; compartments are
            % linearly divided according to the radii
            
            comp_inds=find(pw_dendr_radii>(comp_num-1)/dendrite_n_comp*angle_radius & pw_dendr_radii<=comp_num/dendrite_n_comp*angle_radius & dendr_conn_type==1);
            x_dendr_space=x0(n_cell)+pw_dendr_radii(comp_inds)/angle_radius*space_radius.*cos(dendr_space_angles(comp_inds));
            y_dendr_space=y0(n_cell)+pw_dendr_radii(comp_inds)/angle_radius*space_radius.*sin(dendr_space_angles(comp_inds));

            % add some small Gaussian radial displacement (tuft of cells
            % below spread around)

            dendr_space_shift_radii=branch_sigma*randn(1,length(x_dendr_space));
            dendr_space_shift_angles=2*pi*rand(1,length(x_dendr_space));

            x_dendr_space=round(x_dendr_space+dendr_space_shift_radii.*cos(dendr_space_shift_angles));
            y_dendr_space=round(y_dendr_space+dendr_space_shift_radii.*sin(dendr_space_shift_angles));

            % find values out of range and set them to 1                         
            % (dummy value = 0)             

            inds=find(x_dendr_space>lower_X_margin & x_dendr_space<=X_size-lower_X_margin & y_dendr_space>lower_Y_margin & y_dendr_space<=Y_size-lower_Y_margin);   
            cell_ex_inds{n_cell}(dendr_num,comp_num,1:length(inds))=sub2ind([Y_size*lower_overload_factor X_size*lower_overload_factor 2*angular_steps],y_dendr_space(inds),x_dendr_space(inds),dendr_conn_angle(comp_inds(inds)));

            ex_comp_length=max(ex_comp_length,length(x_dendr_space));




            % find max length of ex compartments in this cell (cycle
            % through all dendrites)
            
            

            % ** do the same for in connections
        
            comp_inds=find(pw_dendr_radii>(comp_num-1)/dendrite_n_comp*angle_radius & pw_dendr_radii<=comp_num/dendrite_n_comp*angle_radius & dendr_conn_type==0);

            x_dendr_space=x0(n_cell)+pw_dendr_radii(comp_inds)/angle_radius*space_radius.*cos(dendr_space_angles(comp_inds));
            y_dendr_space=y0(n_cell)+pw_dendr_radii(comp_inds)/angle_radius*space_radius.*sin(dendr_space_angles(comp_inds));

            dendr_space_shift_radii=branch_sigma*randn(1,length(x_dendr_space));
            dendr_space_shift_angles=2*pi*rand(1,length(x_dendr_space));

            x_dendr_space=round(x_dendr_space+dendr_space_shift_radii.*cos(dendr_space_shift_angles));
            y_dendr_space=round(y_dendr_space+dendr_space_shift_radii.*sin(dendr_space_shift_angles));

            inds=find(x_dendr_space>lower_X_margin & x_dendr_space<=X_size-lower_X_margin & y_dendr_space>lower_Y_margin & y_dendr_space<=Y_size-lower_Y_margin);   

            cell_in_inds{n_cell}(dendr_num,comp_num,1:length(inds))=sub2ind([Y_size*lower_overload_factor X_size*lower_overload_factor 2*angular_steps],y_dendr_space(inds),x_dendr_space(inds),dendr_conn_angle(comp_inds(inds)));

            in_comp_length=max(in_comp_length,length(x_dendr_space));

                
        end % end of dendrite compartment cycles
            
        start_angle_ind=end_angle_ind+1;
        
    end % end of dendrite cycles
    
    % output the cell angular distribution 
    % (for correlation to upper levels)
    
    tot_ex_comp_length=max(tot_ex_comp_length,ex_comp_length);
    tot_in_comp_length=max(tot_in_comp_length,in_comp_length);

end % end cell loop


% set all connections in big matrixes for efficiency (no cell vars, no long for
% loops)

ex_conn_inds=zeros(cell_number*dendrite_number,dendrite_n_comp,tot_ex_comp_length);
in_conn_inds=zeros(cell_number*dendrite_number,dendrite_n_comp,tot_in_comp_length);

dendr_ind=1;

for n_cell=1:cell_number
  if size(cell_ex_inds{n_cell},1)<dendr_num
      cell_ex_inds{n_cell}=[cell_ex_inds{n_cell}; zeros(dendr_num-size(cell_ex_inds{n_cell},1),size(cell_ex_inds{n_cell},2),size(cell_ex_inds{n_cell},3))];
  end          
  if size(cell_in_inds{n_cell},1)<dendr_num
      cell_in_inds{n_cell}=[cell_in_inds{n_cell}; zeros(dendr_num-size(cell_in_inds{n_cell},1),size(cell_in_inds{n_cell},2),size(cell_in_inds{n_cell},3))];
  end      
  for dendr_num=1:1:dendrite_number % dendrite cycle
    for comp_num=1:dendrite_n_comp % dendrite compartment cycle
        ex_conn_inds(dendr_ind,comp_num,1:size(cell_ex_inds{n_cell},3))=cell_ex_inds{n_cell}(dendr_num,comp_num,:);
        in_conn_inds(dendr_ind,comp_num,1:size(cell_in_inds{n_cell},3))=cell_in_inds{n_cell}(dendr_num,comp_num,:);        
    end
    dendr_ind=dendr_ind+1;
  end
end

% set all 0's to 1

ex_conn_inds(ex_conn_inds==0)=1;
in_conn_inds(in_conn_inds==0)=1;

ex_conn_inds=smart_cast(ex_conn_inds,dendrite_number);
in_conn_inds=smart_cast(in_conn_inds,dendrite_number);       





