function RK = compute_nRK(R,K)
% Compute R*K. RK(i,:)= K(idx_phi(i),:) where R = Ng*t-by-Nc repetition
% matrix (cluster indicator)
%
% INPUT:
%   K             Nc-by-Nc Kernel matrix
%   idx_phi       Ng*t-by-1 Cluster indices that features belong to
%
% OUTPUT:
%   RK            R * K

RK = R*K;

end