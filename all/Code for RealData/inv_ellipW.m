function [W,posCR1,posCR2] = inv_ellipW(sizeX,sizeY,posCR1,posCR2,w_fun,lambda_W,eps_W,resolution)
    % Return W with its weight d1 / (d2 + d3) where d1: main axix, d2 % d3
    % distances from tx and rx to the arbitrary point

    grid_X = sizeX:-1/resolution:0;
    grid_Y = 0:1/resolution:sizeY;
    
    y_axis = repmat(grid_X', 1, size(grid_Y,2));  % first coordinate of grid point
    x_axis = repmat(grid_Y, size(grid_X,2), 1);     % second coordinate of grid point
        
    phi1 = norm(posCR1-posCR2).* 0.3048;
    phi2 = sqrt( (x_axis-posCR1(1)).^2 + (y_axis-posCR1(2)).^2 ).* 0.3048 + sqrt( (x_axis-posCR2(1)).^2 + (y_axis-posCR2(2)).^2 ).* 0.3048;
        
    rej_mu_minMat = phi2 < phi1 + eps_W;
    minMu = rej_mu_minMat .* (phi1 + eps_W);
    new_phi2 = minMu + phi2.* mod(rej_mu_minMat+1,2);    
    
    til_mu = 0.5 * sqrt(new_phi2.^2 - phi1.^2);
    
    temp_W = w_fun(phi1, til_mu );
    rej_mu_maxMat = new_phi2 < phi1 + 0.5 * lambda_W ;


    W = rej_mu_maxMat .* temp_W;
    
    
end