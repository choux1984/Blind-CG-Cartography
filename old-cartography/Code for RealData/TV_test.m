function W = TV_test(X)
% close all
% X = magic(5)
[N_x,N_y] = size(X);
N_g = N_x * N_y;

% TV-1
sub_col_1 = [1,-1,zeros(1,N_x-2)];
sub_Mat_1 = zeros(N_x-1,N_x);
for i = 1 : N_x-1
    sub_Mat_1(i,:) = circshift(sub_col_1,[1,i-1]);
end
TV_1 = kron(eye(N_y),sub_Mat_1);

% TV-2
sub_col_2 = [1,zeros(1,N_x-1),-1,zeros(1,N_g - (N_x +1))];
TV_2 = zeros(N_g-N_x,N_g);
for i = 1 : N_g-N_x
    TV_2(i,:) = circshift(sub_col_2,[1,i-1]);
end

W = [TV_1;TV_2];

% keyboard
end