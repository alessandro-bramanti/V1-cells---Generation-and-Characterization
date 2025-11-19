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

I_raw(I_raw<=0)=0;

I=I_raw;

[y_margin_grid,x_margin_grid]=ndgrid(1:1:Y_size,1:1:X_size);

w_plus_ON=conv2(I,LGN_cells_ON_con_center{1},'same');
w_minus_ON=conv2(I,LGN_cells_ON_con_surround{1},'same');

av_lum_ON=conv2(I,LGN_cells_ON_average{1},'same');

C=(w_plus_ON+w_minus_ON)./av_lum_ON; % Weber contrast

C(C<0)=0;
C(av_lum_ON==0)=0;

C(x_margin_grid<=LGN_margin)=0;
C(x_margin_grid>=X_size-LGN_margin)=0;
C(y_margin_grid<=LGN_margin)=0;
C(y_margin_grid>=Y_size-LGN_margin)=0;

I_LGN_ON_con=(C.^2./(C.^2+C50.^2));

C_ON=C;

