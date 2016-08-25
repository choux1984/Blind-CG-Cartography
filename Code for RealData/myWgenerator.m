function [W_col,W] = myWgenerator(Tx_Pos_total,Rx_Pos_total,numCR,t,sizeX,sizeY,w_fun,sizeElip,resolution,beta,eps_W)
% Generate (Nx*Ny-by-numCR*t) matrix which contains collection of vec(W)s
% Require normEllip function

for nS = 1: t
    % Ini_Node
    Ini_Node = (nS-1)*numCR;
    for nCR = 1 : numCR
        posCR1 = Tx_Pos_total(Ini_Node + nCR,:);
        posCR2 = Rx_Pos_total(Ini_Node + nCR,:);
        W(:,:,nCR,nS) = normEllip(sizeX,sizeY,posCR1,posCR2,w_fun,sizeElip,resolution);
%         W(:,:,nCR,nS) = inv_ellipW(sizeX,sizeY,posCR1,posCR2,w_fun,sizeElip,eps_W,resolution);
    end
end

a = W(:,:,1,5);
[sizeX_F,sizeY_F] = size(a);

% for i = 59 : 60
%     for j = 1 : 10
%         W(:,:,j,i) = zeros(sizeX_F,sizeY_F);
%     end
% end

for i = 6
    for j = 81 : 100
        W(:,:,j,i) = zeros(sizeX_F,sizeY_F);
    end
end

for tau = 1 : t

    for j = 1 : numCR  
        w_col = reshape(W(:,:,j,tau),sizeX_F*sizeY_F,1);
        W_col(:, (tau-1)*numCR +j) = sqrt(beta^(t-tau))*w_col;
    end
        
end

end