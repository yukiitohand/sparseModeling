function [ perc ] = PERC( A,Gamma )
% [ perc ] = PERC( A,Gamma )
%   compute the positive exact recovery coefficient (PERC)
%     PERC(Gamma) =      min      1-sum(pinv(AG)*Aj)
%                    j \in Gcmp  |__________________|
%                                         ||
%                                         \/
%                                     PSC(Gamma;j)
%   where PSC is (positive subset coherence).
%
%    Inputs
%       A : Matrix [L x N] of atoms, L is the number of dimension, and
%           N is the number of atoms
%       Gamma: the set of indices (boolean with size [1,N] or integers 
%              in the range[1,N])
%    Outputs
%       perc: scalar, positive exact recovery coefficient
%
%   REFERENCE
%     Y. Itoh, M. F. Duarte and M. Parente, "Perfect Recovery Conditions 
%     for Non-negative Sparse Modeling," in IEEE Transactions on Signal 
%     Processing, vol. 65, no. 1, pp. 69-80, 1 Jan.1, 2017.
%     doi: 10.1109/TSP.2016.2613067
%   
%   Copyright © 2019 Yuki Itoh


[L,N] = size(A);

if islogical(Gamma)
    Gamma = find(Gamma);
    Gamma = reshape(Gamma,1,length(Gamma));
end

if L<length(Gamma)
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

