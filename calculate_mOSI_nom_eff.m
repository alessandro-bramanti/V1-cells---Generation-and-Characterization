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


calculate mOSInom and mOSIeff
%}

s_rel=V1_in_angle_std_devs(test_cell_inds(:))./V1_ex_angle_std_devs(test_cell_inds(:));

mOSInom=zeros(1,cell_number);
mOSIeff=zeros(1,cell_number);

s=zeros(1,angular_steps);
for i=1:cell_number
    a=reshape(V1_soma_resp_map(i,:,:,:),1,[]);
    M=max(a(:));
    for k=1:angular_steps
        grid_mean_resp=reshape(V1_soma_resp_map(i,:,:,k),1,[]);
        s(k)=mean(grid_mean_resp); %#ok<AGROW> 
    end
    orth_angle=V1_max_angle(i)+angular_steps/2*(V1_max_angle(i)<=angular_steps/2)-angular_steps/2*(V1_max_angle(i)>angular_steps/2);
    mOSInom(i)=(s(V1_max_angle(i))-s(orth_angle))/(s(V1_max_angle(i))+s(orth_angle));
    [~,eff_max_angle]=max(s);
    orth_angle=eff_max_angle+angular_steps/2*(eff_max_angle<=angular_steps/2)-angular_steps/2*(eff_max_angle>angular_steps/2);
    mOSIeff(i)=(s(eff_max_angle)-s(orth_angle))/(s(eff_max_angle)+s(orth_angle));
end

figure;
scatter(s_rel,mOSInom,'r^')
hold
scatter(s_rel,mOSIeff,'ks')
legend('mOSI_{nom}','mOSI_{eff}')
