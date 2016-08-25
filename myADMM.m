function est_f = myADMM(A,D,P,T,rho_1,rho_2,mu_f,s_check,num_itr)

%%% INPUT
%     D               N_y(N_x-1)-by-Ng discrete gradient matrix
%     P               Permutation matrix to make row-major vector of 'f'

[N_Dx,N_g] = size(D);
N_Px = size(P,1);

til_A = (2/T) * A' * A;
til_D = rho_1 * D'*D;
til_P = rho_2 * P'*P;
til_H = (til_A + til_D + til_P);
til_s = (2/T) * A' * s_check;
til_V = til_D + rho_2 * eye(N_g);

% Initialize
f = randn(N_g,1);
v = P * f;
gamma_x = zeros(N_Dx,1);
gamma_y = zeros(N_Dx,1);
gamma_z = zeros(N_Px,1);
d_x = zeros(N_Dx,1);
d_y = zeros(N_Dx,1);

for k = 1 : num_itr
    k
    % [S1] Dual ascent for aux variables
    gamma_x = gamma_x + rho_1 * (D*f - d_x);
    gamma_y = gamma_y + rho_1 * (D*v - d_y);
    gamma_z = gamma_z + rho_2 * (P*f - v);
    % [S2] Updates of aux variables d_x and d_y
    d_x = wthresh(D*f + gamma_x./rho_1,'s', mu_f/rho_1);
    d_y = wthresh(D*v + gamma_y./rho_1,'s', mu_f/rho_1);

    % [S3] Updates of SLF
    f = ( til_H ) \(D' *(rho_1 * d_x - gamma_x) + P' *(rho_2 *v -gamma_z) + til_s);
    
    % [S4] Updates of aux variable v
    v = til_V \ ( D' *(rho_1 * d_y - gamma_y) + rho_2 * P * f + gamma_z);
    
    cost(k) = (1/T) * norm(A*f - s_check,2)^2;
end
est_f = f;
figure
plot(cost);
end