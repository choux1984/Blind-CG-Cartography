function compare_W(w_fun,w_est,N_x,N_y,lambda_W,eps_W,K_std,mu_w,Nc,idx_sign,rho)
%   Display the original and estimated weighted functions
% 
% INPUT:
%   w             Original function
%   w_est         Estimated function
%   N_x/N_y       Number of rows and columns of a spatial loss field
%


    % Evaluation of w and w_est
    min_phi1 = 0.1;
    max_phi1 = sqrt((1-N_x)^2 + (1-N_y)^2);
    max_phi2 = (N_x-1) + (N_y-1);
    max_phi2 = max_phi2; 
    
%     rng_phi1 = min_phi1:0.5:max_phi1;
%     rng_phi2 = min_phi1:0.1:max_phi2;  
    
    rng_phi1 = 16:0.15:18;
    rng_phi2 = 16:0.01:18;  
    
%     rng_phi1 = 10:0.3:15;
%     rng_phi2 = 12:0.01:16;  
%     
    len_phi1 = size(rng_phi1,2);
    len_phi2 = size(rng_phi2,2);
    for i = 1 : len_phi1
        for j = 1 : len_phi2
            if rng_phi1(i)+eps_W > rng_phi2(j) && rng_phi1(i) <= rng_phi2(j)            
%                 hat_w(i,j) = NaN;
%                 w_o(i,j) = NaN;
                hat_w(i,j) = w_est([rng_phi1(1,i); rng_phi2(1,j)]);
                w_o(i,j) = myinvElip(w_fun,rng_phi1(1,i),rng_phi2(1,j),lambda_W,eps_W);
            elseif rng_phi1(i) > rng_phi2(j)
                hat_w(i,j) = NaN;
                w_o(i,j) = NaN;                
            else
                hat_w(i,j) = w_est([rng_phi1(1,i); rng_phi2(1,j)]);
                w_o(i,j) = myinvElip(w_fun,rng_phi1(1,i),rng_phi2(1,j),lambda_W,eps_W);

            end
        end
    end
    
    h = figure
    plot(rng_phi2(1:1:end)',w_o(:,1:1:end)');
    hold on
    plot(rng_phi2(1:1:end)',idx_sign .* hat_w(:,1:1:end)','--');
    title(sprintf('Reconstructed w with K-std=%g, mu_w=%g, rho=%g, Nc = %g', K_std, mu_w,rho,Nc))
    file_name=sprintf('est_w_K_std_%g_mu_w_%g_rho_%g_Nc_%g.fig',K_std, mu_w,rho,Nc);
    saveas(h,file_name)
% 
%     % 2-D plot while fixing phi_1
%     figure
%     num_plot = 4;
%     itv_plot = floor(len_phi1/num_plot); %interval of phi_1 for the subplots
%     for n = 1 : num_plot
%         subplot(2,2,n);
%         vec_inds = (n-1)*itv_plot + 1:n*itv_plot + 1;
%         plot(rng_phi2',w_o(vec_inds,:)')
%         for k = 1:length(vec_inds)
%             leg{k} = sprintf('\\phi_1 = %g', rng_phi1(vec_inds(k)));
%         end
%         legend(leg,'Location','southwest');
% 
%         hold on
%         plot(rng_phi2',hat_w(vec_inds,:)','--')
%         hold off
%         xlabel('\phi_2')
%         ylabel('w')      
%     end
% 
    % 3D plot of Original and Estimated functions
%     figure
%     surf(rng_phi2,rng_phi1,w_o)
%     title('Ground truth w')
%     ylabel('\phi_1')
%     xlabel('\phi_2')  
%     figure
%     
%     surf(rng_phi2,rng_phi1,idx_sign .* hat_w)
%     title('Estimated w')
%     ylabel('\phi_1')
%     xlabel('\phi_2')    
    

end

function w = myinvElip(w_fun,rng_phi1,rng_phi2,lambda_W,eps_W)

    rej_ind = rng_phi2 < rng_phi1 + eps_W;
    maxMu = rej_ind .* (rng_phi1 + eps_W);
    new_phi2 = maxMu + rng_phi2.* mod(rej_ind+1,2);    
    
    til_mu = 0.5 * sqrt(new_phi2.^2 - rng_phi1.^2);
    
    temp_W = w_fun(rng_phi1, til_mu );
    rej_max_ind = new_phi2 < rng_phi1 + 0.5 * lambda_W ;
    w = rej_max_ind .* temp_W;

end
