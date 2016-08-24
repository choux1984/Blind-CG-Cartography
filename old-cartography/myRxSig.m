function [s_check,Tx_pos,Rx_pos] = myRxSig(t_slots,F,sigma_2,w_fun,lambda_W,eps_W)  
    
end


function [W,posCR1,posCR2] = inv_ellipW(sizeX,sizeY,numCR,w_fun,lambda_W,eps_W)
    % Return W (btw rx and tx at rbitrary positions), following the inverse area elliptical model
    % INPUT:
    %   sizeX         Number of rows of W
    %   sizeY         Number of columns of W
    %
    % OUTPUT:
    %   W            (SizeX-by-SizeY) weight matrix
    %   posCR1       2-by-1 vector of the tx location
    %   posCR2       2-by-1 vector of the rx location
    
    assert(numCR==1);

    
end

function [W,posCR1,posCR2] = inv_ellipW_bound(sizeX,sizeY,numCR,w_fun,lambda_W,eps_W)
    % Return W (btw rx and tx at rbitrary positions on a boundary of the area), following the inverse area elliptical model   
    assert(numCR==1);

    rows = repmat((1:sizeX)', 1, sizeY);  % first coordinate of grid point
    cols = repmat(1:sizeY, sizeX, 1);     % second coordinate of grid point
    
    rnd_idx = randperm(4);
   
    switch rnd_idx(1)
        case 1
            posCR1 = [sizeX*rand(1)+1,1];       
        case 2
            posCR1 = [1,sizeY*rand(1)+1];           
        case 3
            posCR1 = [sizeX,sizeY*rand(1)+1];        
        case 4
            posCR1 = [sizeX*rand(1)+1,sizeY];                
    end
    
    switch rnd_idx(2)
        case 1
            posCR2 = [sizeX*rand(1)+1,1];       
        case 2
            posCR2 = [1,sizeY*rand(1)+1];           
        case 3
            posCR2 = [sizeX,sizeY*rand(1)+1];        
        case 4
            posCR2 = [sizeX*rand(1)+1,sizeY];                
    end
    
    
    phi1 = norm(posCR1-posCR2);
    phi2 = sqrt( (rows-posCR1(1)).^2 + (cols-posCR1(2)).^2 ) + sqrt( (rows-posCR2(1)).^2 + (cols-posCR2(2)).^2 );
        
    rej_mu_minMat = phi2 < phi1 + eps_W;
    minMu = rej_mu_minMat .* (phi1 + eps_W);
    new_phi2 = minMu + phi2.* mod(rej_mu_minMat+1,2);    
    
    til_mu = 0.5 * sqrt(new_phi2.^2 - phi1.^2);
    
    temp_W = w_fun(phi1, til_mu );
    rej_mu_maxMat = new_phi2 < phi1 + 0.5 * lambda_W ;

    W = rej_mu_maxMat .* temp_W;
    
    
end