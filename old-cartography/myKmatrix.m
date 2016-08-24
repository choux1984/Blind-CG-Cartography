function K = myKmatrix(phi_i,myKfunc)
    % Kernel matrix generator

    % INPUT:
    %   phi_i               2-by-Nc matrix of kerner inputs
    %   myKfunc             Kernel function
    %
    % OUTPUT:
    %       K               Nc x Nc kernel matrix 
    
    Nc = size(phi_i,2);
    K = zeros(Nc);
    
    for m0 = 1 : Nc
        for m1 = 1 : Nc
             K(m0,m1) = myKfunc(phi_i(:,m0),  phi_i(:,m1) );
        end
    end
    min_eig = min(eig(K));
    if min_eig <= 0
        K = K - min_eig* eye(size(K,1));
    end
    K = K + 10^8*eps* eye(Nc);
    
    if rank(K) < Nc
        Nc
        rank_K = rank(K)
    end
end