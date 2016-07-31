function [alpha,evl_pnt,mK,R] = estimate_alpha(N_x,N_y,total_Tx_pos,total_Rx_pos,K_var,f,lambda,s_check,clustering,Nc,myKfunc)

    t = size(total_Tx_pos,2);
    Ng = N_x * N_y;

    [phi_i, phi_col] = myPhi(N_x,N_y,t,total_Tx_pos,total_Rx_pos,Ng);

    if clustering
        % Rotate phi to feet with matlab function 'kmeans'
        phi = phi_col';
        % idx_phi := cluster index where the entries of phi belong to / C := a set of centroids
        [idx_phi,C] = kmeans(phi,Nc,'Replicates',2,'MaxIter',500,'start','uniform','emptyaction','singleton'); 
        [tempC, ia, ic] = unique(C,'rows'); % Rmove redundant centroids
        idx_phii = ic(idx_phi);
        C = tempC;
        Nc = size(ia,1);
        idx_phi = idx_phii;

        R = zeros(t*Ng,Nc); 
        for i = 1 : Nc
            R(:,i)= (idx_phi==i);
        end

        C = C';
        K = myKmatrix(clustering,Nc,K_var,C,myKfunc);

        evl_pnt = C; % evl_pnt := evluation point of kerenl.
        
        figure;
        plot(C(:,1),C(:,2),'rx')
        hold all
        for i = 1 : Nc
            plot(phi(idx_phi==i,1),phi(idx_phi==i,2),'.')
        end

    else  % no clustering

        K = myKmatrix(t,Ng,K_var,phi_i,myKfunc);
        R = eye(t * Ng);
        evl_pnt = phi_i; % evl_pnt := evluation point of kerenl.

    end
    
    dim_eye = size(K,1);
    mK = K - eye(dim_eye)*(-10^5*eps+min(eig(K)));

    K1 = kron(eye(t),f') * R * mK;
    alpha = (R'*kron(eye(t),f) * K1 + lambda * t * eye(dim_eye))\(R'*(kron(eye(t),f) * s_check));
end


function [phi_i, phi_col] = myPhi(N_x,N_y,t,total_Tx_pos,total_Rx_pos,Ng)
    % idx_mat := matrix containing the coordinates of the area
    idx_rows = repmat((1:N_x)', 1, N_y);
    idx_cols = repmat(1:N_y, N_x, 1);
    v_idx_rows = reshape(idx_rows, N_x*N_y,1);
    v_idx_cols = reshape(idx_cols, N_x*N_y,1);
    v_idx = [v_idx_rows, v_idx_cols]';
    
    %t = size(total_Tx,2);
    phi_i = zeros(2,N_x*N_y,t);

    for m0 = 1 : t
        Tx_m0 = total_Tx_pos(:,m0);
        Rx_m0 = total_Rx_pos(:,m0);

        for i = 1 : N_x * N_y
            phi_i(:,i,m0) = [norm(Tx_m0 - Rx_m0); norm(Tx_m0 - v_idx(:,i)) +  norm(v_idx(:,i) - Rx_m0)];
        end
        % phi_col := collection of possible phi
        phi_col(:,1+(m0-1)*Ng : m0*Ng) = phi_i(:,:,m0);
    end    
end


function K = myKmatrix(t,dim_subMat,K_var,phi_i,myKfunc)

K = zeros(t * dim_subMat, t * dim_subMat); % Initialize K matrix
K_mm = zeros(dim_subMat,dim_subMat); % Initialize the block gram matirx K_mm
size_feat = size(phi_i,1); % size_feat := size of the input feature vector for kernel
phi_m1 = phi_i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construction of K %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m0 = 1 : t
    for m1 = 1 : t
        % dist_idx_m0(1) := dummy matrices to calculate the distances among
        % all possible coordinates
        dist_idx_m0 = repmat(reshape(phi_i(:,:,m0),size_feat * dim_subMat,1), 1, dim_subMat);
        dist_idx_m1 = repmat(phi_m1(:,:,m1), dim_subMat, 1);
        dist2_idx = (dist_idx_m0 - dist_idx_m1).*(dist_idx_m0 - dist_idx_m1);
        % calculate K_{m,m'}
        for i = 1 : dim_subMat
            % calculate the (squared) distance row-wisely
            dist2_row = sum(dist2_idx(1 + size_feat * (i -1):size_feat*i,:));
            dist_row = sqrt(dist2_row);
            K_mm(i,:) = myKfunc(K_var,dist_row);                    
        end
        K(1 + (m0-1) * dim_subMat : m0 * dim_subMat, 1 + (m1-1) * dim_subMat : m1 * dim_subMat) = K_mm;
    end
end

end
