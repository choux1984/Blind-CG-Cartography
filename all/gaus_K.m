function K = gaus_K(N_x,N_y,gamma,t,num_CR_set,total_col_pos1,total_col_pos1)

% Find coordinate matrix and its vectorization
sizeX = N_x;
sizeY = N_y;
vec_gx = [1:sizeX]';
vec_gy = [1:sizeY];
for i = 1 : sizeX
    g(:,i,1) = vec_gx;
    g(i,1:sizeY,2) = vec_gy;
end

g_vec_x = reshape(g(:,:,1),sizeX*sizeY,1);
g_vec_y = reshape(g(:,:,2),sizeX*sizeY,1);

for m_1 = 1 : t * num_CR_set
    for m_2 = 1 : t * num_CR_set
        phi_11 = norm(total_col_pos1(:,m_1) - total_col_pos2(:,m_1),2);
        phi_21 = norm(total_col_pos1(:,m_2) - total_col_pos2(:,m_2),2);
    end    
end

end
