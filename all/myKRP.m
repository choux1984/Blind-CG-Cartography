function C = myKRP(A,B)
    % Column dim. of A and B should be the same
    [n,l] = size(A);
    [m] = size(B,1);

    C = zeros(n*m,l); %initialization
    for i = 1 : l
       ab = B(:,i) * A(:,i)';
       C(:,i) = ab(:); % vectorization
    end
end