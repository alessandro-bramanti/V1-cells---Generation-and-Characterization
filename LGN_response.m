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
of the paper mentioned above - arXiv preprint available

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

I_LGN_con=(C.^2./(C.^2+C50.^2));

I_conv=zeros(Y_size,X_size,angular_steps*2);    
    
for m=1:angular_steps
    for n=1:1:size(LGN_collinear_Gabor,1)
%        conv_surf=conv2(I_LGN_con,LGN_collinear_Gabor{n,m}.*(LGN_collinear_Gabor{n,m}>0),'same');
%        conv_surf_2=conv2(I_LGN_coff,-LGN_collinear_Gabor{n,m}.*(LGN_collinear_Gabor{n,m}<0),'same');
        Gabor_conv=conv2(I_LGN_con,LGN_collinear_Gabor{n,m},'same');

        %Gabor_conv=conv_surf+conv_surf_2;

        Gabor_conv(1:2*Gabor_transverse_radii+2,:)=0;
        Gabor_conv(size(Gabor_conv,1)-2*Gabor_transverse_radii-2:1:size(Gabor_conv,1),:)=0;            
        Gabor_conv(:,1:2*Gabor_transverse_radii+2)=0; 
        Gabor_conv(:,size(Gabor_conv,2)-2*Gabor_transverse_radii-2:1:size(Gabor_conv,2))=0;            

%        Gabor_conv=sign(Gabor_conv).*(abs(Gabor_conv)>0.1);      
        
        I_conv(:,:,m)=Gabor_conv.*(Gabor_conv>0);
        I_conv(:,:,m)=I_conv(:,:,m).*(I_conv(:,:,m)>0);
        I_conv(:,:,m+angular_steps)=-Gabor_conv.*(Gabor_conv<0);      
        I_conv(:,:,m+angular_steps)=I_conv(:,:,m+angular_steps).*(I_conv(:,:,m+angular_steps)>0);

    end
end

I_conv(dummy_pixel_y,dummy_pixel_x,dummy_pixel_z)=0;

