function CM = myCovMat(sizeX,sizeY)

% sizeX(Y) := size of the coordinate on X(Y) axis

% vec_gx(gy) := vectorize the every coordinate along the axis x (y)
vec_gx = [1:sizeX]';
vec_gy = [1:sizeY];

for i = 1 : sizeX
    g(:,i,1) = vec_gx;
    g(i,1:sizeY,2) = vec_gy;
end

g_vec_x = reshape(g(:,:,1),sizeX*sizeY,1);
g_vec_y = reshape(g(:,:,2),sizeX*sizeY,1);

cnt = 1;
for i = 1: sizeX * sizeY
    for j = 1 : sizeX * sizeY
        % comlumn-wisely calculate covariance matrix 
        vec_CM(cnt,1) = sqrt((g_vec_x(j,1) - g_vec_x(i,1))^2 + (g_vec_y(j,1) - g_vec_y(i,1))^2) ;
        cnt = cnt + 1;
    end
end
temp_CM = reshape(vec_CM,sizeX*sizeY,sizeX*sizeY);
CM = exp(-1*temp_CM);

end