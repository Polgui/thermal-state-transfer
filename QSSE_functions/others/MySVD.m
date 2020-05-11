function [U,L1,V, projection_norm]=MySVD(C21, maxSchmidtrank, Eps)
% [U,L1,V, projection_norm]=MySVD(C21, maxSchmidtrank, Eps)
% Performs the SVD of the matrix C21, but retains only the maxSchmidtrank
% higher singular values larger than Eps
%size(C21)
[U,S,V]=svd(C21,'econ');
L1=diag(S);                     % Singular values vector

% Truncation of the singular values
Bond=length(L1);
for j=length(L1)-1:-1:1
    if L1(j:end)'*L1(j:end)/(L1'*L1) < Eps
        Bond=j;
    end
end

Bondlength1=min(Bond,maxSchmidtrank);
if Bondlength1==maxSchmidtrank
    warning('mytest:maxrank','Bond dimension saturated for MAXSCHMIDTRANK=%i. This might yield unprecise results.',maxSchmidtrank)
    warning('off','mytest:maxrank')
end
L1=L1(1:Bondlength1);
projection_norm=norm(L1);
L1=L1./norm(L1);
U=U(:,1:Bondlength1);
V=V(:,1:Bondlength1);


end

