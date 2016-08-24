function alpha = est_alpha(K,idx_phi,f,s_check,mu_w)
    % Estimator for alpha
    
    % INPUT:
    %   K              Kernel matrix
    %   idx_phi        (n_measurements * Ng) x 1: Cluster indices where phi_i's belong to 
    %                  idx_phi(1) means that the location of 1 on the 1st row of R matrix 
    %   f              (Ng-by-1) initialization of the spatial loss field (SLF)
    %   s_check        Noisy received signals
    %   lambda         Parameter for Kernel smoothness regularizer
    %
    % OUTPUT:
    %   alpha          Estimate of alpha
    
    t = size(s_check,1);
    Nc = size(K,2);

    temp_fR = myfRMat(t,Nc,f,idx_phi);    
    K1 = temp_fR * K;
    
%     A = (K1' * K1 + lambda * t * K);
%     alpha = A\(K1' * s_check);
    A = (temp_fR' * K1 + mu_w * t * eye(Nc));
    alpha = A\(temp_fR' * s_check);    
end


function fR_Mat = myfRMat2(t,Nc,f,idx_phi)
% fR_mat = ( I \kron f^T ) R

fR_Mat = zeros(t,Nc);
Ng = size(f,1);



% for i = 1 : Ng
%     for tau = 1 : t
%         ind_R = (tau - 1) * Ng + i;
%         ind_op = idx_phi(ind_R);
%         fR_Mat(tau,ind_op) = fR_Mat(tau,ind_op) + f(i,1); 
%     end
% end
for tau = 1 : t
    ind_R = (tau - 1) * Ng;
    for i = 1 : Ng
        ind_op = idx_phi(ind_R+i);
        fR_Mat(tau,ind_op) = fR_Mat(tau,ind_op) + f(i,1); 
    end
end

end

function fR_Mat = myfRMat(t,Nc,f,idx_phi)
% fR_mat = ( I \kron f^T ) R
% i) R_tilde:=[R_1' ; ... ; R_T'] * f
% ii) Rf_vec := sum(R_tilde,2)
% iii) fR_Mat = unvec(Rf_vec)'

Ng = size(f,1);

idx_phi_mat = reshape(idx_phi,Ng,t);
row_inds = idx_phi_mat + Nc*ones(Ng,1)*(0:t-1);
col_inds = repmat(1:Ng,1,t)';
entry_vals = repmat(f,t,1);

Rf_vec = full(sum(sparse(row_inds(:),col_inds,entry_vals,Nc*t,Ng,length(idx_phi)),2));

fR_Mat = reshape(Rf_vec,Nc,t)';


end

