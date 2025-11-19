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

short for surf with some settings
%}

function show_surf (varargin)


if length(varargin)>1 && ~isempty(varargin{2}(1))
    surf_vals=zeros(2*varargin{3}(1)+size(varargin{1},1),2*varargin{2}(1)+size(varargin{1},2));
    surf_vals(varargin{3}(1)+1:1:varargin{3}(1)+size(varargin{1},1),varargin{2}(1)+1:1:varargin{2}(1)+size(varargin{1},2))=varargin{1};
else
    surf_vals=varargin{1};
end

if length(varargin)>3
    surf_vals(1,1)=varargin{4}(1);
    surf_vals(1,2)=varargin{5}(1);
end

figure;
h=surf(surf_vals);
set(h,'EdgeAlpha',0);
view([0 -90]);
colorbar;
axis equal;