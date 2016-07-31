function [evl_pnts,idx_phi,phi_col] = generate_feature_vecs(N_x,N_y,Tx_pos,Rx_pos,Nc,clustering_type,lambda_W)
    % Generator for a Kernal matrix K and a mapping matrix R

    % INPUT:
    %   Tx_pos         2-by-(num_CR_set * t_slots) collections of Tx coordinates   
    %   Rx_pos         2-by-(num_CR_set * t_slots) collections of Rx coordinates
    %   Nc             Number of clusters
    %                  Nc = [] means no clustering (Nc = n_measurements * Ng)
    %   clustering_type
    %                  'kmeans'
    %                  'random'
    % OUTPUT:
    %   evl_pnts       2-by-Nc mat: collections of centroids
    %   idx_phi        (n_measurements * Ng) x 1: Cluster indices where phi_i's belong to
    %   phi_col        2-by-(n_measurements*Ng): collection of all feature
    %                   vectors
    
    debug = 0;
    
    % Varialbles:
    n_measurements = size(Tx_pos,2); % Number of measurements = t_slots * num_CR_set(=1)
    
    % Dependent variables
    Ng = N_x * N_y; % Number of grid points
    
%     [phi_col] = myPhi(N_x,N_y,Tx_pos,Rx_pos);  % each column of phi_col contains a vector phi
    [phi_col] = myPhiEllipseOnly(N_x,N_y,Tx_pos,Rx_pos,lambda_W); % only for the features within the ellipse 
    
    clustering = ~isempty(Nc);
    if clustering
        switch clustering_type
            case 'kmeans'
                % idx_phi := cluster index where the entries of phi belong to C := a set of centroids
                [idx_phi,C] = kmeans(phi_col',Nc,'Replicates',2,'MaxIter',500,'start','sample','emptyaction','singleton');
                
                if  (sum(sum(isnan(C)))>0) || ( sum(sum(C==Inf))>0 )
                    error('kmeans returned NaN');
                end                
                [C, ia, ic] = unique(C,'rows'); % Rmove redundant centroids
                idx_phi = ic(idx_phi);
                Nc = size(ia,1);
              
                evl_pnts = C'; % evl_pnt := evluation points of kerenl.
            case 'random'
                % choose Nc feature vectors (cols of phi_col) --> centroids   
                temp_phi_col = phi_col';
                uni_phi_col = unique(temp_phi_col,'rows')';
                
                idx_centroids = randperm(size(uni_phi_col,2),Nc);
                rand_centroids = uni_phi_col(:,idx_centroids);                
                             
                % assign each col of phi_col to a centroid --> idx_phi
                matD = pdist2(phi_col',rand_centroids');% matD := distance matrix between coordinates and centroids
                [~,idx_phi] = min(matD');
                idx_phi = idx_phi';
                
                evl_pnts = rand_centroids;
        end
        if debug
             plot_clusters(evl_pnts,phi_col,idx_phi);
        end
        plot_clusters(evl_pnts,phi_col,idx_phi);
    else  % no clustering
        R = eye(n_measurements * Ng);
        evl_pnts = phi_col; % evl_pnt := evluation point of kerenl.
    end
    
    
   
end

function [phi_col] = myPhi(N_x,N_y,Tx_pos,Rx_pos)
    % Computation of all feature vectors associated with the sensor
    % locations

    % INPUT:
    %   Tx_pos     2-by-(num_CR_set * t_slots) collections of Tx coordinates   
    %   Rx_pos     2-by-(num_CR_set * t_slots) collections of Rx coordinates        
    % OUTPUT:
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

function [phi_col] = myPhiEllipseOnly(N_x,N_y,Tx_pos,Rx_pos,lambda_W)
    % Computation of all feature vectors associated only within the ellipse
    % defined by the sensors 

    % INPUT:
    %   Tx_pos     2-by-(num_CR_set * t_slots) collections of Tx coordinates   
    %   Rx_pos     2-by-(num_CR_set * t_slots) collections of Rx coordinates        
    % OUTPUT:
    %   phi_col    2-by-(# of grids within an ellipse * (num_CR_set * t_slots)) matricization of phi_i

    Ng = N_x * N_y;
    % idx_mat := matrix containing the coordinates of the area
    idx_rows = repmat((1:N_x)', 1, N_y);
    idx_cols = repmat(1:N_y, N_x, 1);
    grid_points_coord = [idx_rows(:), idx_cols(:)]';
    
    n_measurements = size(Tx_pos,2);

    cnt = 1;
    for m0 = 1 : n_measurements
        Tx_m0 = Tx_pos(:,m0);
        Rx_m0 = Rx_pos(:,m0);
        
        for i = 1 : Ng
            phi1 = norm(Tx_m0 - Rx_m0);
            phi2 = norm(Tx_m0 - grid_points_coord(:,i)) +  norm(grid_points_coord(:,i) - Rx_m0);
            if phi2 <= phi1 + 0.5 * lambda_W
                phi_col(:,cnt) = [phi1; phi2];
            else
                phi_col(:,cnt) = [phi1; phi1 + 0.5 * lambda_W];
            end
            cnt = cnt + 1;

        end
        
%         phi_col(:,1+(m0-1)*Ng : m0*Ng) = phi_i(:,:,m0);
    end    
end


function K = old_myKmatrix(n_measurements,dim_subMat,phi_i,myKfunc)
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

function plot_clusters(evl_pnt,phi_col,idx_phi)


    % Plot phi's and clusters
    figure;
    plot(evl_pnt(1,:),evl_pnt(2,:),'rx')
    hold all
    for i = 1 : size(evl_pnt,2)
        %phi = phi';
        plot(phi_col(1,idx_phi==i),phi_col(2,idx_phi==i),'.')
    end

    xlabel('\phi_1')
    ylabel('\phi_2')
    

end

