function [ perc ] = PERC( A,Gamma )
% [ perc ] = PERC( Phi,lambda )
%   compute the positive exact recovery coefficient
%   PERC = 
%    Inputs
%       Phi : dictionary matrix [D x N], D is the number of dimension, and
%       N is the number of atoms
%    Outputs
%       lambda : subset indices of columns of Phi [1 x N] array of binary
%       elements or [1 x n] (n<N) integer indices

[D,N] = size(A);

if islogical(Gamma)
    Gamma = find(Gamma);
    Gamma = reshape(Gamma,1,length(Gamma));
end

if D<length(Gamma)
    error(['Subset size is too large']);
end

cGamma = setdiff(1:N,Gamma);
Ag = A(:,Gamma);
Agcmp = A(:,cGamma);
AgtAg = Ag'*Ag;
% AgtAginv = inv(AgtAg);
Agpinv = pinv(Ag);

perc = 1 - max(sum(pinv(A(:,Gamma))*A(:,cGamma),1),[],2);

% perc2 = ( 1 - max(sum(Agcmp' * Agpinv', 2)) );


end

