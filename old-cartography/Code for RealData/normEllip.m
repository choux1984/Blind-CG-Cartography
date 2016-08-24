function W = normEllip(sizeX,sizeY,posCR1,posCR2,w_fun,lambda_W,resolution)
    % Return W with its weight 1 /sqrt(d1) where d1: main axix, d2 % d3
    % distances from tx and rx to the arbitrary point
    
    grid_X = sizeX:-1/resolution:0;
    grid_Y = 0:1/resolution:sizeY;
    
    y_axis = repmat(grid_X', 1, size(grid_Y,2));  % first coordinate of grid point
    x_axis = repmat(grid_Y, size(grid_X,2), 1);     % second coordinate of grid point
        
    phi1 = norm(posCR1-posCR2).* 0.3048;
    phi2 = sqrt( (x_axis-posCR1(1)).^2 + (y_axis-posCR1(2)).^2 ).* 0.3048 + sqrt( (x_axis-posCR2(1)).^2 + (y_axis-posCR2(2)).^2 ).* 0.3048;
       
    temp_W = w_fun(phi1);
    rej_mu_maxMat = phi2 <= phi1 + 0.5 * lambda_W ;

    W = rej_mu_maxMat .* temp_W;    
end