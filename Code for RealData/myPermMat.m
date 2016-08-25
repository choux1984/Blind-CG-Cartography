function P = myPermMat(N_x,N_y)
% generate a permuation matrix to change column-major form vector to
% row-major form

N_g = N_x * N_y;
seedVec = [1,zeros(1,N_g-1)];
for i = 1 : N_y
    seedMat(i,:) = circshift(seedVec,[0,(i-1)*N_x]);
end

P = seedMat;
for i = 2 : N_x
    P = [P;circshift(seedMat,[0,(i-1)])];
end

end