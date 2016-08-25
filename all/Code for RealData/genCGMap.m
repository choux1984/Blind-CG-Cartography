function estCGMap = genCGMap(F,x_rx,sizeX,sizeY,resolution,gamma,w_fun,sizeElip)
% function to generate CG map with the rx located at x_rx
%%% INPUT
%   F             N_x-by-N_y spatial loss field
%%% OUTPUT
%

[N_x,N_y] = size(F);
estCGMap = zeros(N_x,N_y);

grid_Y = sizeX:-1/resolution:0; % Origin is at bottom left
grid_X = 0:1/resolution:sizeY;
y_axis = repmat(grid_Y', 1, size(grid_X,2));  % first coordinate of grid point
x_axis = repmat(grid_X, size(grid_Y,2), 1);     % second coordinate of grid point

Gain_Mat = 64 * ones(N_x,N_y);
PathLoss_Mat = myPathloss(x_axis,y_axis,x_rx,gamma);
Shadow_Mat = myShadow(grid_X,grid_Y,sizeX,sizeY,x_rx,w_fun,sizeElip,resolution,F);

figure
imagesc(Shadow_Mat);

a = imagesc(Shadow_Mat);
imcontrast(a)

estCGMap = Gain_Mat - PathLoss_Mat + Shadow_Mat;

figure
imagesc(estCGMap);

b = imagesc(estCGMap);
imcontrast(b)

% keyboard



end

function PathLoss_Mat = myPathloss(x_axis,y_axis,x_rx,gamma)
% function to generate path-loss map
%%% INPUT
%   gamma           Pathloss exponent
%%% OUTPUT
%   PathLoss_Mat

    phi1 = sqrt( (x_axis-x_rx(1)).^2 + (y_axis-x_rx(2)).^2 ).* 0.3048;
    PathLoss_Mat = gamma .* 10 * log(phi1)./log(10);

end

function Shadow_Mat = myShadow(grid_X,grid_Y,sizeX,sizeY,x_rx,w_fun,sizeElip,resolution,F)

size_i = size(grid_X,2);
size_j = size(grid_Y,2);

Shadow_Mat = zeros(size_j,size_i);

for i = 1 : size_i
    for j = 1 : size_j
        x_test = [grid_X(i);grid_Y(j)];
        W = normEllip(sizeX,sizeY,x_rx,x_test,w_fun,sizeElip,resolution);
        Shadow_Mat(i,j) = sum(sum(W.*F));
    end
end

end