function UMPS = UnitaryTensor(gamma,Omega,Nbins, dimS, dimB,dt )

% sm=[0 1;
%     0 0];

a=diag(sqrt(1:dimB-1),1);   % Ladder
sm=diag(sqrt(1:dimS-1),1);   % Ladder

idB=eye(dimB);

Op=0;

for j=1:Nbins
Op=Op+sqrt(gamma*dt)*...
    (KroneckerProduct(sm,KroneckerProduct(idB,j-1),a',KroneckerProduct(idB,Nbins-j))...
    -KroneckerProduct(sm',KroneckerProduct(idB,j-1),a,KroneckerProduct(idB,Nbins-j)));
end

Op=Op+(Omega/2)*KroneckerProduct(sm-sm',KroneckerProduct(idB,Nbins))*dt + 0*1i*KroneckerProduct(sm'*sm,KroneckerProduct(idB,Nbins))*dt;

UMPS=expm(Op);
UMPS=reshape(UMPS,[ones(1,Nbins)*dimB,dimS,ones(1,Nbins)*dimB,dimS]); 
UMPS=permute(UMPS,[(Nbins:-1:1),Nbins+1,Nbins+1+(Nbins:-1:1),2*Nbins+2]); 
    

end

