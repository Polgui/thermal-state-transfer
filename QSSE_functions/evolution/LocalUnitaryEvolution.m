function [Gamma,Lambda,ErrorTracker] = LocalUnitaryEvolution(Gamma,Lambda,siteV,Unitary_V,maxSchmidtrank, Eps)

% Physical dimensions
dimB=size(Gamma{siteV+1},3);
dimS=size(Gamma{siteV},3);

Nbins=length(size(Unitary_V))/2 -1; %Number of time bins interacting

%% Merge everything

Gamma_tot=Gamma{siteV};
for l=1:Nbins
Gamma_tot=Lambda_multiplication(Gamma_tot,Lambda{siteV+(l-1)},l+1);
Gamma_tot=tensor_contraction(Gamma_tot,Gamma{siteV+l},l+1,1);
end

Gamma_tot=permute(Gamma_tot,[1 Nbins+2 3:Nbins+1 Nbins+3 2]);

BondlengthL=size(Gamma_tot,1); % left MPS dimension of Gamma_feedback 
BondlengthR=size(Gamma_tot,2); % right MPS dimension of Gamma_feedback 

if siteV > 1
    Gamma_tot=Lambda_multiplication(Gamma_tot,Lambda{siteV-1},1);
end

Gamma_tot=Lambda_multiplication(Gamma_tot,Lambda{siteV+Nbins},2);


%% Apply the unitary

Gamma_tot=tensor_contraction(Gamma_tot,Unitary_V,(3:length(size(Gamma_tot))),((length(size(Unitary_V))/2 +1):length(size(Unitary_V))));
Gamma_tot=permute(Gamma_tot,[1 (3:length(size(Gamma_tot))) 2]);

%% Bring everything back in the original form
ErrorTracker=ones(1,2);
for l=1:Nbins
    Gamma_tot=reshape(Gamma_tot,[BondlengthL*dimB,dimB^(Nbins-l)*dimS*BondlengthR]); % reshape into a matrix

    % Apply the SVD
    [U,L1,V, projection_norm]=MySVD(Gamma_tot,maxSchmidtrank,Eps);
    ErrorTracker(1)=ErrorTracker(1)*projection_norm;
    ErrorTracker(2)=min(ErrorTracker(2),projection_norm);
    
    Bondlength=length(L1);
    Lambda{siteV+(l-1)}=L1;
    % New state
    AL_feedback=permute(reshape(U,[BondlengthL,dimB,Bondlength]),[1,3,2]);
    
    if siteV+(l-2)>0
        AL_feedback=Lambda_multiplication(AL_feedback,1./Lambda{siteV+(l-2)},1);
    end
    Gamma{siteV+(l-1)}=AL_feedback;
    % All the rest
    V=V(:,1:Bondlength);
    
    Gamma_tot=diag(L1)*V';
    BondlengthL=Bondlength;

end

% New state for the system
Gamma{siteV+Nbins}=Lambda_multiplication(permute(reshape(V',[BondlengthL,dimS,BondlengthR]),[1,3,2]),1./Lambda{siteV+Nbins},2);

