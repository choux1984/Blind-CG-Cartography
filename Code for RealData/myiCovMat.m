function inv_R = myiCovMat(sizeX,sizeY,resolution)

% Find the inverse of spatial covariance matrix R
grid_X = sizeX:-1/resolution:0;
grid_Y = 0:1/resolution:sizeY;

y_axis = repmat(grid_X', 1, size(grid_Y,2));  % first coordinate of grid point
x_axis = repmat(grid_Y, size(grid_X,2), 1);     % second coordinate of grid point

grid_y_vec = y_axis(:);
grid_x_vec = x_axis(:);

N_x = size(grid_X,2);
N_y = size(grid_Y,2);

u = 1;
for i = 1: N_x * N_y
    for j = 1 : N_x * N_y
        vec_R(u,1) = sqrt((grid_x_vec(j,1) - grid_x_vec(i,1))^2 + (grid_y_vec(j,1) - grid_y_vec(i,1))^2)* 0.3048;
        u = u + 1;
    end
end
temp_R = reshape(vec_R,N_x*N_y,N_x*N_y);
R = exp(-1*temp_R);
inv_R = inv(R);

end