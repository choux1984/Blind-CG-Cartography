% Simulations for blind channel-gain cartography
% This is the main file

filename = 'Map_15_15.csv';

% estimation parameters

mu_w_vec = 0.26;
mu_f_vec = 1e-4;
K_std_vec = 0.15; %0.19
Nc_vec = 10;
rho_vec = 1e-3;
t_slot_vec = 10;

% mu_w_vec = 0.2;
% mu_f_vec = 1e-4;
% K_std_vec = 0.15; %0.17
% Nc_vec = 10;
% t_slot_vec = 10;
reg_f_type = 'TV'; % choose regularizer type for the SLF: 'T' = Tikhonov; 'S' = l_1; 'TV' = total variation 

size_mu_w = length(mu_w_vec);
size_mu_f = length(mu_f_vec);
size_K = length(K_std_vec);
size_Nc = length(Nc_vec);
size_t = length(t_slot_vec);
size_rho = length(rho_vec);
tic
for i = 1:size_mu_w
    for l = 1 : size_mu_f
        for j = 1 : size_K
            for k = 1 : size_Nc
                for t = 1 : size_t
                    for r = 1 : size_rho
                        mu_w = mu_w_vec(i);
                        mu_f = mu_f_vec(l);
                        K_std = K_std_vec(j);
                        Nc = Nc_vec(k);
                        t_slots = t_slot_vec(t);
                        rho = rho_vec(r);
                        simulation_01(mu_w,mu_f,K_std,Nc,t_slots,reg_f_type,filename,rho);
                    end
                end
            end
        end
    end
end

toc