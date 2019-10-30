function [ erc ] = ERC( Phi,lambda )
% [ erc ] = ERC( Phi,lambda )
%   compute the exact recovery coefficient
%   ERC = 
%    Inputs
%       Phi : dictionary matrix [D x N], D is the number of dimension, and
%       N is the number of atoms
%       lambda : subset indices of columns of Phi [1 x N] array of binary
%       elements or [1 x n] (n<N) integer indices
%    Outputs
%       erc: scalar, exact recovery coefficient

[D,N] = size(Phi);

if islogical(lambda)
    lambda = find(lambda);
    lambda = reshape(lambda,1,length(lambda));
end

if D<length(lambda)
    error(['Subset size is too large']);
end

nlambda = setdiff(1:size(Phi,2),lambda);

erc = 1 - operatorNorm(pinv(Phi(:,lambda))*Phi(:,nlambda),1,1);

end

