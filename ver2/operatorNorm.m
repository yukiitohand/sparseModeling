function [ val ] = operatorNorm( A,p,q )
% [ val ] = operatorNorm( A,p,q )
%   compute the operator norm ||A||_{p,q}
%    Inputs
%       A : matrix
%       p : domain of the operator norm [1,2,'inf']
%       q : co-domain of the operator norm ['1','2','inf']
%    Outputs
%       val : double precision scalar

if isnumeric(p)
    p = num2str(p);
end

if isnumeric(q)
    q = num2str(q);
end

if (strcmp(p,'1') && strcmp(q,'1'))
    sumAcol = sum(abs(A),1);
    val = max(sumAcol);
elseif (strcmp(p,'1') && strcmp(q,'2'))
    sumAcol2nrm = sqrt(sum(A.^2,1));
    val = max(sumAcol2nrm);
elseif (strcmp(p,'1') && strcmp(q,'inf'))
    val = max(max(abs(A)));
elseif (strcmp(p,'2') && strcmp(q,'2'))
    val = svds(A,1);
elseif (strcmp(p,'2') && strcmp(q,'inf'))
    sumArow2nrm = sqrt(sum(A.^2,2));
    val = max(sumArow2nrm);
elseif (strcmp(p,'inf') && strcmp(q,'inf'))
    sumArow = sum(abs(A),2);
    val = max(sumArow);
else
    error(['Unrecognized domain: p=' p ' and q=' q]);
end

end
