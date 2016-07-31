function R = myiCovMat(sizeX,sizeY,var_shadow)

% Find the inverse of spatial covariance matrix R
vec_gx = [1:sizeX]';
vec_gy = [1:sizeY];
for i = 1 : sizeX
    g(:,i,1) = vec_gx;
    g(i,1:sizeY,2) = vec_gy;
end

g_vec_x = reshape(g(:,:,1),sizeX*sizeY,1);
g_vec_y = reshape(g(:,:,2),sizeX*sizeY,1);
%
u = 1;
for i = 1: sizeX * sizeY
    for j = 1 : sizeX * sizeY
        d_R(u,1) = sqrt((g_vec_x(j,1) - g_vec_x(i,1))^2 + (g_vec_y(j,1) - g_vec_y(i,1))^2);
        u = u + 1;
    end
end
temp_R = reshape(d_R,sizeX*sizeY,sizeX*sizeY);
R = var_shadow * exp(-1*temp_R);
R = inv(R);

end