function MSE_w = estimate_error_w(w_fun,w_est,phi_col)
% Approximate MSE of functions evalulated only at phi_col
% INPUT:
%   w_fun               Original function with two inputs phi_1, phi_2
%   w_est               Estimated function with vector input [phi_1;phi_2]
%   phi_col             2-by-(n_measurements*Ng): collection of all feature vectors
%                            
%
% OUTPUT:
%    MSE_w              Approximated MSE of estimated weight function

N_K  = size(phi_col,2);

if N_K > 30000
    N_K = 30000;
end

eval_w = w_fun(phi_col(1,1:N_K), phi_col(2,1:N_K));
eval_w_est = zeros(1,N_K);
for k = 1 : N_K
    eval_w_est(k) = w_est(phi_col(:,k));
end
sum_err = sum((eval_w-eval_w_est).^2);

MSE_w = sum_err/(norm(eval_w,2)^2);



end