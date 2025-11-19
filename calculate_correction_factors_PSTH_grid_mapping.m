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


calculate correction factors
%}

grid_80_perc=zeros(1,cell_number);

s_rel=V1_in_angle_std_devs(test_cell_inds(:))./V1_ex_angle_std_devs(test_cell_inds(:));

for i=1:cell_number
    a=reshape(V1_soma_resp_map(i,:,:,:),1,[]);
    grid_80_perc(i)=prctile(a,80)*1000;
end

figure;
scatter(s_rel,grid_80_perc,'k')

correction_factors=15./grid_80_perc;

figure;
scatter(s_rel,correction_factors,'k')

s_rel=s_rel';

fun=@(x)sum((correction_factors-x(1)-x(2).*(s_rel-x(3))).^2);
options=optimset('MaxFunEvals',10000);

x0=[1 0.4 0];

[c,fval,exitflag,output]=fminsearch(fun,x0, options);

figure;
scatter(s_rel,correction_factors,'k')
hold on
scatter(s_rel,c(1)+c(2).*(s_rel-c(3)))

