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

represents spot_field_map (ON and OFF regions in red and blue,
respectively)
receives up to 4 arguments:
    pos_field_map (mandatory) calculated for the n-th cell 
    neg_field_map (mandatory) calculated for the n-th cell 
    crisp (optional) =0 for a smooth color scale, =1 for a crisp color scale (mandatory)
    thresh, minimum somatic voltage represented. cuts off asymptotic tails,
    0.0001 can be a good choice, remember that the voltages are in volts
%}

function show_field_map(varargin)

pos_field_map=varargin{1};
neg_field_map=varargin{2};
crisp=varargin{3}(1);
if nargin>3
    thresh=varargin{4}(1);
else
    thresh=0;
end

field_map=pos_field_map-neg_field_map;

margin=5;

M=max(abs(field_map));

x0=find(max(abs(field_map))>thresh*M,1,'first');
x0=max(x0-margin,2);

x1=find(max(abs(field_map))>thresh*M,1,'last');
x1=min(x1+margin,size(field_map,2));

y0=find(max(abs(field_map'))>thresh*M,1,'first');
y0=max(y0-margin,1);

y1=find(max(abs(field_map'))>thresh*M,1,'last');
y1=min(y1+margin,size(field_map,1));

red_field_map=field_map(y0:y1,x0:x1);
red_pos_field_map=pos_field_map(y0:y1,x0:x1);
red_neg_field_map=neg_field_map(y0:y1,x0:x1);

if length(varargin)>4
    red_pos_field_map=red_pos_field_map/varargin{5}(1);    
    red_neg_field_map=red_neg_field_map/varargin{6}(1);    
else
    red_pos_field_map=(red_pos_field_map-min(abs(red_pos_field_map(:))))/(max(abs(red_pos_field_map(:)))-min(abs(red_pos_field_map(:))));    
    red_neg_field_map=(red_neg_field_map-min(abs(red_neg_field_map(:))))/(max(abs(red_neg_field_map(:)))-min(abs(red_neg_field_map(:))));    
end

if crisp
    cmap=[reshape(1.0*(red_pos_field_map>red_neg_field_map),[],1),zeros(numel(red_pos_field_map),1),reshape(1.0*(red_pos_field_map<red_neg_field_map),[],1)];
else
    cmap=[reshape(red_pos_field_map,[],1),zeros(numel(red_pos_field_map),1),reshape(red_neg_field_map,[],1)];
end
c_mat=reshape(1:size(cmap,1),size(red_pos_field_map,1),[]);

figure;
h=surf(red_field_map,c_mat);
colormap(cmap);
view([0 90]);
set(h,'EdgeAlpha',0);

    



