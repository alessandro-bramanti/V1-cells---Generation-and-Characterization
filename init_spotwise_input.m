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

initializes the spot stimuli (bright and dark)
spot_radius in pixels, 0<spot_intensity<=1
%}

spot_radius=2;
spot_intensity=1;

I_raw=zeros(size(x_unrot));
I_raw((x_unrot-X_size/2).^2+(y_unrot-Y_size/2).^2<=spot_radius^2)=spot_intensity; %#ok<NASGU>

preliminary_LGN_response;

LGN_response;

pos_spot_resp_0=I_conv;

pos_spot_resp_0(pos_spot_resp_0<0.05)=0;

I_raw=spot_intensity*ones(size(x_unrot));
I_raw((x_unrot-X_size/2).^2+(y_unrot-Y_size/2).^2<=spot_radius^2)=0;

preliminary_LGN_response;

LGN_response;

neg_spot_resp_0=I_conv;

neg_spot_resp_0(:,1:X_margin,:)=0;
neg_spot_resp_0(:,X_size-X_margin+1:X_size,:)=0;

neg_spot_resp_0(1:Y_margin,:,:)=0;
neg_spot_resp_0(Y_size-Y_margin+1:Y_size,:,:)=0;

neg_spot_resp_0(neg_spot_resp_0<0.05)=0;