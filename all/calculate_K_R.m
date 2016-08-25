function [evl_pnts,mK,R,idx_phi,phi_col] = calculate_K_R(N_x,N_y,Tx_pos,Rx_pos,Nc,myKfunc)
    % Generator for a Kernal matrix K and a mapping matrix R

    % INPUT:
    %   Tx_pos         2-by-(num_CR_set * t_slots) collections of Tx coordinates   
    %   Rx_pos         2-by-(num_CR_set * t_slots) collections of Rx coordinates
    %   myKfunc        Kernel function
    %   Nc             Number of clusters
    %                  Nc = [] means no clustering (Nc = n_measurements * Ng)
    % OUTPUT:
    %   evl_pnts       2-by-Nc mat: collections of centroids
    %   mK             Kernel matrix
    %                  [Clustering] Nc-by-Nc
    %                  [Non-clustering] (n_measurements * Ng)-by- (n_measurements * Ng)
    %   R              (n_measurements * Ng) x Nc: Mapping matrix from phi_i to corresponding bar_phi
    %   idx_phi        (n_measurements * Ng) x 1: Cluster indicies where phi_i's belong to
    %   phi_col        (n_measurements*N_g-by-2): collection of all feature
    %                   vectors
    
    % Varialbles:
    n_measurements = size(Tx_pos,2); % Number of measurements = t_slots * num_CR_set(=1)
    
    % Dependent variables
    Ng = N_x * N_y; % Number of grid points
    
    

    [phi_i, phi_col] = myPhi(N_x,N_y,Tx_pos,Rx_pos);  % each column of phi_col contains a vector phi
    
    clustering = ~isempty(Nc);
    if clustering
        % idx_phi := cluster index where the entries of phi belong to C := a set of centroids
        [idx_phi,C] = kmeans(phi_col',Nc,'Replicates',2,'MaxIter',500,'start','sample','emptyaction','singleton'); 
        
        if  (sum(sum(isnan(C)))>0) || ( sum(sum(C==Inf))>0 )
           error('kmeans returned NaN'); 
        end
        
        [C, ia, ic] = unique(C,'rows'); % Rmove redundant centroids
        idx_phi = ic(idx_phi);
        Nc = size(ia,1);
       
        Rt= zeros(Nc,n_measurements*Ng); 
        Rt( idx_phi+Nc*(0:length(idx_phi)-1)') = 1; 
        R = Rt'; 
        
        C = C';
        K = myKmatrix(clustering,Nc,C,myKfunc);

        evl_pnts = C; % evl_pnt := evluation points of kerenl.

    else  % no clustering
       

        K = myKmatrix(n_measurements,Ng,phi_i,myKfunc);
        R = eye(n_measurements * Ng);
        evl_pnts = phi_i; % evl_pnt := evluation point of kerenl.

    end
    
    dim_eye = size(K,1);

    mK = K + eye(dim_eye)*((1e5)*eps+min(eig(K)));
    
end


function [phi_i, phi_col] = myPhi(N_x,N_y,Tx_pos,Rx_pos)
    % Computation of all feature vectors associated with the sensor
    % locations

    % INPUT:
    %   Tx_pos     2-by-(num_CR_set * t_slots) collections of Tx coordinates   
    %   Rx_pos     2-by-(num_CR_set * t_slots) collections of Rx coordinates        
    % OUTPUT:
    %   phi_i      2-by-Ng-by-(num_CR_set * t_slots) fibers along 1st
    %              dimension contain feature vectors
    %   phi_col    2-by-(Ng * (num_CR_set * t_slots)) matricization of phi_i

    Ng = N_x * N_y;
    % idx_mat := matrix containing the coordinates of the area
    idx_rows = repmat((1:N_x)', 1, N_y);
    idx_cols = repmat(1:N_y, N_x, 1);
    grid_points_coord = [idx_rows(:), idx_cols(:)]';
    
    n_measurements = size(Tx_pos,2);
    phi_i = zeros(2,Ng,n_measurements);

    for m0 = 1 : n_measurements
        Tx_m0 = Tx_pos(:,m0);
        Rx_m0 = Rx_pos(:,m0);

        for i = 1 : Ng
            phi_i(:,i,m0) = [norm(Tx_m0 - Rx_m0); norm(Tx_m0 - grid_points_coord(:,i)) +  norm(grid_points_coord(:,i) - Rx_m0)];
        end
        
        phi_col(:,1+(m0-1)*Ng : m0*Ng) = phi_i(:,:,m0);
    end    
end


function K = myKmatrix(n_measurements,dim_subMat,phi_i,myKfunc)
    % Kernel matrix generator

    % INPUT:
    %   n_measurements      Number of measurements  
    %   dim_subMat          Ng-by-Ng Dim. of a submatrix K_i
    %                       [Clustering] Ng = Nc 
    %   phi_i               2-by-Ng-by-n_measurements tensor of kerner inputs
    %   myKfunc             Kernel function
    %
    % OUTPUT:
    %       K               (n_measurements * dim_subMat)^2 kernel matrix 
    
    K = zeros(n_measurements * dim_subMat, n_measurements * dim_subMat); % Initialize K matrix
    K_mm = zeros(dim_subMat,dim_subMat); % Initialize the block gram matirx K_mm
    size_feat = size(phi_i,1); % size_feat := size of the input feature vector for kernel
    phi_m1 = phi_i;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Construction of K %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for m0 = 1 : n_measurements
        for m1 = 1 : n_measurements
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
                K_mm(i,:) = myKfunc(dist_row);                    
            end
            K(1 + (m0-1) * dim_subMat : m0 * dim_subMat, 1 + (m1-1) * dim_subMat : m1 * dim_subMat) = K_mm;
        end
    end

end