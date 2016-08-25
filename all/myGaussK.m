function outK = myGaussK(var,input,parameter)
    % dist := Distance between input and parameter
    dist = norm(input - parameter,2);
    outK = exp(-(2 * var)^(-1) * dist^2);
end