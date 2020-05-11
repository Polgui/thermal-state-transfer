function [Gamma,Lambda, ExcitationsS1, ExcitationsS2,rhoS1,rhoS2, rho_final]...
        = FullEvolution(Gamma0,Lambda0,Gammain, Lambdain, gammaS1, gammaS2, N_wav, dimB, dimS1, dimS2, dt, Nt, MS1S2, maxSchmidtrank, Eps, showprogress)
    
warning('on','mytest:maxrank')


if showprogress
    prog_bar = waitbar(0);
    titleHandle = get(findobj(prog_bar,'Type','axes'),'Title');
    set(titleHandle,'FontSize',20);
    waitbar(0,prog_bar,sprintf('0.00%',0));
    tic
end

siteS1=length(Gamma0)-1;
siteS2=siteS1-2*MS1S2-1;

Gamma=Gamma0;
Lambda=Lambda0;

ExcitationsS1=zeros(1,Nt);
ExcitationsS2=zeros(1,Nt);

aS1=diag(sqrt(1:dimS1-1),1);
aS2=diag(sqrt(1:dimS2-1),1);


%% Time evolution  
for k=1:Nt 
    
     for j=1:N_wav
        Gamma{siteS1+j+1}=Gammain{j};
        Lambda{siteS1+j+1}=Lambdain{j};
     end

     Unitary_S1 = UnitaryTensor(gammaS1(k),0,1, dimS1, dimB,dt);
     Unitary_S2 = UnitaryTensor(gammaS2(k),0,1, dimS2, dimB,dt);
     
     [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,siteS1+1,siteS1+3,maxSchmidtrank, Eps);
     
     [Gamma,Lambda] = LocalUnitaryEvolution(Gamma,Lambda,siteS1,Unitary_S1,maxSchmidtrank, Eps);
     siteS1=siteS1+1;
     
     [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,siteS1,siteS1+1,maxSchmidtrank, Eps);
     siteS1=siteS1+1;
     
     [Gamma,Lambda] = LocalUnitaryEvolution(Gamma,Lambda,siteS2,Unitary_S2,maxSchmidtrank, Eps);
     %Lambda{siteS2}
     siteS2=siteS2+1;
     
     [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,siteS2,siteS2+1,maxSchmidtrank, Eps);
     siteS2=siteS2+1;
     
     tens=Lambda_multiplication(Gamma{siteS1}, Lambda{siteS1-1},1);
     rhoS1=tensor_contraction(tens, conj(tens),[1 2],[1 2]);
     ExcitationsS1(k)=trace(rhoS1*(aS1)'*aS1);

     tens=Lambda_multiplication(Lambda_multiplication(Gamma{siteS2},Lambda{siteS2},2), Lambda{siteS2-1},1);
     rhoS2=tensor_contraction(tens, conj(tens),[1 2],[1 2]);
     ExcitationsS2(k)=trace(rhoS2*(aS2)'*aS2);

     
     % Display progression
     if showprogress
         waitbar(k/Nt,prog_bar,sprintf('%3.2f%%',100*k/Nt));
     end
end

     [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,siteS1+1,siteS2+1,maxSchmidtrank, Eps);
     siteS1=siteS1+1;

     tens=Lambda_multiplication(Lambda_multiplication(Gamma{siteS2}, Lambda{siteS2-1},1),Lambda{siteS2},2);
     tens=tensor_contraction(tens,Gamma{siteS2+1},2,1);
     tens=Lambda_multiplication(tens,Lambda{siteS2+1},3);
     tens=permute(tens,[1 3 2 4]);
     tens=reshape(tens,[size(tens,1) size(tens,2) size(tens,3)*size(tens,4)]);
     rho_final=tensor_contraction(tens, conj(tens),[1 2],[1 2]);

if showprogress
    toc
    close(prog_bar)
end