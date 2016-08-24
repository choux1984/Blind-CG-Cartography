function out_Rep = myRepThm(alpha,evl_pnt,input,myKfunc)
% Evaluation of function with kernels

% INPUT:
%   alpha               2-by-Nc coefficients of kerner inputs
%   evl_pnt             2-by-(t*Ng) collection of phi's to evaluate kernels
%   input               2-by-1 kernel input 
%   myKfunc             Kernel function
% OUTPUT:
%   out_Rep             function value evaluated with the input


% num_para := total number of parameter 
num_para = size(alpha,1);
vec_K = zeros(1,num_para);

num_mea = size(evl_pnt,3);
Nc = size(evl_pnt,2);
% cnt = 1;
% for i = 1 : num_mea
%     for j = 1 : num_grid
%       
%         % vec_K := represent vector
% %         vec_K(cnt) = exp(-(2 * K_var)^(-1) * dist^2);
%         vec_K(cnt) = myKfunc(input,evl_pnt(:,j,i));
%         cnt = cnt + 1;
%     end
% end
% out_Rep = vec_K * alpha;

% vec_K = myKfunc(input*ones(1,Nc),evl_pnt(:,:));

for i = 1 : Nc
    vec_K(1,i) = myKfunc(input,evl_pnt(:,i));
end

out_Rep = vec_K * alpha;

end