function [s_check,Tx_pos,Rx_pos,test_W] = myRxSigCompare(t_slots,F,sigma_2,w_fun,lambda_W,eps_W)  
    %   Received signal generator
    % 
    % INPUT:
    %   num_CR_set    Number of pairs of sensors per time slot
    %   t_slots       Number of time slots
    %   F             Spatial loss field
    %   sigma_2       Variance of noise
    %   w_fun         Weight function
    %   lambda_W      Parameter to determine the semi-minor axis of the first Fresnel zone
    %   eps_W         Parameter to determine the threshold for the small
    %
    % OUTPUT:
    %   s_check       (num_CR_set * t_slots)-by-1 vector of corrupted noise signals 
    %   Tx_pos        2-by-(num_CR_set * t_slots) collections of Tx coordinates   
    %   Rx_pos        2-by-(num_CR_set * t_slots) collections of Rx coordinates

    % pre allocate for spped
    num_CR_set = 1;
    s = NaN(num_CR_set,1); % noiseless measurements
    [N_x,N_y] = size(F);
     
    % collect the locations of all the measurements
    Tx_pos = zeros(2,t_slots);
    Rx_pos = zeros(2,t_slots);
% 
%     itv = 0.1;
%     
%     pool1 = 
    
    test_W = zeros(N_x*N_y,t_slots);
    
    for tau = 1 : t_slots
%         [W,coll_posCR1,coll_posCR2] = inv_ellipW(N_x,N_y,num_CR_set,w_fun,lambda_W,eps_W); 
        [W,norm_W,coll_posCR1,coll_posCR2] = test_inv_ellipW(N_x,N_y,num_CR_set,w_fun,lambda_W,eps_W); 
        test_W(:,tau) = norm_W(:);
%         test_W(:,tau) = W(:);
        
        
        % embed the locations of Tx and Rx
        ini_point = 1 + (tau - 1) * num_CR_set;
        Tx_pos(:,ini_point:ini_point + (num_CR_set -1)) = coll_posCR1; %Tx location
        Rx_pos(:,ini_point:ini_point + (num_CR_set -1)) = coll_posCR2; %Rx location

        for j = 1 : num_CR_set  
            s(j) = sum(sum(W(:,:,j).*F)); % signal w/o noise
        end

        k = (tau-1)*num_CR_set;
        s_check(k + 1 : k + num_CR_set,1) = s + sqrt(sigma_2)*randn(num_CR_set,1); % signal vector with noise
    end
end

function [W,norm_W,posCR1,posCR2] = test_inv_ellipW(sizeX,sizeY,numCR,w_fun,lambda_W,eps_W)
    % Return W with its weight d1 / (d2 + d3) where d1: main axix, d2 % d3
    % distances from tx and rx to the arbitrary point
    assert(numCR==1);

    rows = repmat((1:sizeX)', 1, sizeY);  % first coordinate of grid point
    cols = repmat(1:sizeY, sizeX, 1);     % second coordinate of grid point
        
    posCR1 = round(diag([sizeX-1,sizeY-1])*rand(2,1)*100)/100+1;
    posCR2 = round(diag([sizeX-1,sizeY-1])*rand(2,1)*100)/100+1;
    
    phi1 = norm(posCR1-posCR2);
    phi2 = sqrt( (rows-posCR1(1)).^2 + (cols-posCR1(2)).^2 ) + sqrt( (rows-posCR2(1)).^2 + (cols-posCR2(2)).^2 );
        
    rej_mu_minMat = phi2 < phi1 + eps_W;
    minMu = rej_mu_minMat .* (phi1 + eps_W);
    new_phi2 = minMu + phi2.* mod(rej_mu_minMat+1,2);    
    
    til_mu = 0.5 * sqrt(new_phi2.^2 - phi1.^2);
    
    temp_W = w_fun(phi1, til_mu );
    rej_mu_maxMat = new_phi2 < phi1 + 0.5 * lambda_W ;


    W = rej_mu_maxMat .* temp_W;
    
%     % Test weight
    temp_W_norm = ones(size(W,1),size(W,2)) ./ sqrt(phi1);
    norm_W = rej_mu_maxMat .* temp_W_norm ;
    
end