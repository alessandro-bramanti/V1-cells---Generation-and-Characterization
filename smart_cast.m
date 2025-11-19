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

function cast_var = smart_cast (input_var,varargin)

if length(varargin)==1
    f=varargin{1}(1);
else
    f=1;
end;

M=double(max(input_var(:)))*f;

if M<=2^8-1
    cast_var=uint8(input_var);
elseif M<=2^16-1
    cast_var=uint16(input_var);
elseif M<=2^32-1
    cast_var=uint32(input_var);
else
    cast_var=input_var;
end;