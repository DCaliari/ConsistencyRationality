function [a, b_u] = E_fun_nonpar_simplified(pi_par,Ns,F)

%
F_u_y=F.*repmat(pi_par',Ns,1);
%
a = F_u_y./sum(F_u_y,2);
b_u=mean(a,1)';
%
end
