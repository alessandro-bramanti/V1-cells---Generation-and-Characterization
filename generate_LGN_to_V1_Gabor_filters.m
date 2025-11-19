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

function [LGN_collinear_Gabor]=generate_LGN_to_V1_Gabor_filters(Gabor_collinear_radii,Gabor_collinear_sfs,Gabor_transverse_radii,angular_steps)

LGN_collinear_Gabor=cell(length(Gabor_collinear_radii)*length(Gabor_collinear_sfs),angular_steps);

for r=1:length(Gabor_transverse_radii)
    
    [y,x]=ndgrid(-2*Gabor_collinear_radii(r):1:2*Gabor_collinear_radii(r),-2*Gabor_collinear_radii(r):1:2*Gabor_collinear_radii(r));

    for i=1:angular_steps
        
        rot_y=-x*sin(pi/angular_steps*(i-1))+y*cos(pi/angular_steps*(i-1));

        ind=1;

        for m=1:length(Gabor_collinear_radii)
            for k=1:length(Gabor_collinear_sfs)
                gab=(abs(rot_y)<=Gabor_transverse_radii).*(1.*(rot_y>0)-1.*(rot_y<0));
                gab=gab.*(x.^2+y.^2<=Gabor_collinear_radii(r)^2);

                LGN_collinear_Gabor{ind,i}=gab/sum(gab(gab(:)>0));
                ind=ind+1;
            end
        end
    end
end