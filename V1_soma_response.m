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

Solves the dendritic circuit
%}

function [V_soma,V_dendr_out] = V1_soma_response (n_neurons,V_ex_conn_inds,V_in_conn_inds,I_conv,varargin)

global g_a g_L g_NMDA g_SOMA

I_conv(1)=0; % for safety, in case it's forgotten outside

E_m=-70e-3; % rest membrane potential

SOMA_fact=g_a(1)/(g_SOMA+g_a(1));

dendrite_number=5;

dendrite_n_comp=4;

g_GABA=g_NMDA/2;

V_vals=E_m:0.1e-3:0;

ex_syn_threshold=0;

NMDA_exp_shift=23.7e-3;
NMDA_exp_const=12.5e-3;

V_soma=zeros(n_neurons,1);
V_dendr_out=zeros(dendrite_number,1);

V_step=abs(V_vals(2)-V_vals(1));

dendr_inds=[1:n_neurons*dendrite_number]'; %#ok<NBRAK>

N_ex_syn=5*sum(I_conv(V_ex_conn_inds),3);
dendr_inds=dendr_inds(max(N_ex_syn,[],2)>ex_syn_threshold);

if isempty(dendr_inds)
    return;
end

V=zeros(length(dendr_inds),dendrite_n_comp);
I_leak=zeros(length(dendr_inds),dendrite_n_comp); 
I_NMDA=zeros(length(dendr_inds),dendrite_n_comp); 
I_ax=zeros(length(dendr_inds),dendrite_n_comp); 
V_dendr=zeros(length(dendr_inds),1);


continue_search=ones(size(dendr_inds));

N_ex_syn=sum(I_conv(V_ex_conn_inds(dendr_inds,:,:)),3);
N_in_syn=sum(I_conv(V_in_conn_inds(dendr_inds,:,:)),3);

ind=length(V_vals)*ones(size(dendr_inds));
exit_loop=0;
while ~exit_loop

    V(:,1)=V_vals(ind); 
    I_leak(:,1)=(g_SOMA*g_a(1)/(g_SOMA+g_a(1))+g_GABA(1)*N_in_syn(:,1)+g_L(1)).*(V(:,1)-E_m);
    I_NMDA(:,1)=-N_ex_syn(:,1).*V(:,1)*g_NMDA(1)./(1+exp(-(V(:,1)+NMDA_exp_shift)/NMDA_exp_const));
    I_ax(:,2)=I_NMDA(:,1)-I_leak(:,1); % outflowing

    for k=2:dendrite_n_comp-1
        V(:,k)=V(:,k-1)-I_ax(:,k)/g_a(k);
        I_leak(:,k)=(g_L(k)+g_GABA(k)*N_in_syn(:,k)).*(V(:,k)-E_m);
        I_NMDA(:,k)=-N_ex_syn(:,k).*V(:,k)*g_NMDA(k)./(1+exp(-(V(:,k)+NMDA_exp_shift)/NMDA_exp_const));
        I_ax(:,k+1)=I_NMDA(:,k)-I_leak(:,k)+I_ax(:,k); 
    end

    V(:,dendrite_n_comp)=V(:,dendrite_n_comp-1)-I_ax(:,dendrite_n_comp)/g_a(dendrite_n_comp);
    I_leak(:,dendrite_n_comp)=(g_L(dendrite_n_comp)+g_GABA(dendrite_n_comp)*N_in_syn(:,dendrite_n_comp)).*(V(:,dendrite_n_comp)-E_m);
    I_NMDA(:,dendrite_n_comp)=-N_ex_syn(:,dendrite_n_comp).*V(:,dendrite_n_comp)*g_NMDA(dendrite_n_comp)./(1+exp(-(V(:,dendrite_n_comp)+NMDA_exp_shift)/NMDA_exp_const));

    diff_bal=(sum(I_NMDA(:,:),2)-sum(I_leak(:,:),2));

    if ind<length(V_vals)-1
        growth_sign=sign(diff_bal);
        if sum(growth_sign~=prev_growth_sign)>0
            V_dendr(growth_sign~=prev_growth_sign)=V(growth_sign~=prev_growth_sign,1)-E_m+V_step*abs(abs(diff_bal(growth_sign~=prev_growth_sign))./(abs(prev_diff_bal(growth_sign~=prev_growth_sign))+abs(diff_bal(growth_sign~=prev_growth_sign))));
            continue_search=continue_search.*(growth_sign==prev_growth_sign & ind>1);
            prev_growth_sign=growth_sign;
            if max(continue_search)==0
                exit_loop=1;
            end
        end
        ind=ind-continue_search;
        continue_search(ind==0)=0;
        ind(ind==0)=1;
        if max(continue_search)==0
            V_dendr(continue_search==1)=0;
            exit_loop=1;
        end
    elseif ind==length(V_vals)-1
        prev_growth_sign=sign(diff_bal);
        ind=ind-continue_search;
    else
        ind=ind-continue_search;
    end
    
    prev_diff_bal=diff_bal;

end


neuron_inds=ceil(dendr_inds/(dendrite_number));

for i=1:length(neuron_inds)
    V_soma(neuron_inds(i))=V_soma(neuron_inds(i))+V_dendr(i);
end

V_dendr_out(dendr_inds)=V_dendr;

V_soma=V_soma*SOMA_fact;

