%% SETTINGS
clearvars
addpath('QSSE_functions/analysis/')
addpath('QSSE_functions/evolution/')
addpath('QSSE_functions/others/')
addpath('QSSE_functions/unitaries/')

disp('CURRENT PID:') 
feature getpid
 
c2 = parcluster('local');
c2.NumWorkers = 12;
mypool=parpool(c2.NumWorkers); 
  
Nth=0.00001;
  
gamma_m=1;
t0=5;

dt=0.05/gamma_m;          % adaptative Timestep
DelayS1S2=0.1;

Asynchronicity=0;

Totaltime=20 ;        % Total simulation time
t = (dt:dt:Totaltime);

maxSchmidtrank=100; % Maximum Schmidt rank for the SVD 
Eps=1e-12;

dimB=10;
dimS1=round(dimB);
dimS2=dimS1;

showprogress=1;
showtransfer=1;
saveresults=0;

N_wav=2;                                % Number of waveguides

Nt=round(Totaltime/dt);    % Total number of timesteps
MS1S2=round(DelayS1S2/dt);

gammaS1=(t<t0).*exp(gamma_m*(t-t0))./(2-exp(gamma_m*(t-t0))) + (t>=t0).*(t<=2*t0)*gamma_m;

    GammainB = reshape(eye(dimB),1,dimB,dimB);
    GammainA = reshape(eye(dimB),dimB,1,dimB);
    LambdaBA=((Nth+1)/Nth).^(-(0:dimB-1)/2);
    LambdaBA=LambdaBA'./norm(LambdaBA); 
        
     Gammain=cell(1,N_wav);
    Gammain{1}=GammainB;
    Gammain{2}=GammainA;
    
    Lambdain=cell(1,N_wav);
    Lambdain{1}=LambdaBA;
    Lambdain{2}=1;     
    Gamma0=cell(1,1);
    Lambda0=cell(1,1);
    
    Gamma0{1}=reshape([1 zeros(1,dimS2-1)],[1,1,dimS2]);
    Lambda0{1}=Lambdain{2};
    
    for j=1:MS1S2
        Gamma0{(j-1)*2+2}=Gammain{1};
        Lambda0{(j-1)*2+2}=Lambdain{1};
        Gamma0{(j-1)*2+3}=Gammain{2};
        Lambda0{(j-1)*2+3}=Lambdain{2};
    end
    
    RealS1(:,:,1)=[1 0];
    RealS1(:,:,2)=[0 1];
    for iter=3:dimS1
        RealS1(:,:,iter)=zeros(1,2);
    end
    Gamma0{2*MS1S2+2}= RealS1;
    Lambda0{2*MS1S2+2}=[1/sqrt(2); 1/sqrt(2)];


    AuxS1(:,:,1)=[1; 0];
    AuxS1(:,:,2)=[0; 1];
    Gamma0{2*MS1S2+3}= AuxS1;
    Lambda0{2*MS1S2+3}=1; 
    
    
Psi_init=1/sqrt(2)*(KroneckerProduct( [1; 0],[1; zeros(dimS2-1,1)]) - KroneckerProduct([0; 1],[0;1; zeros(dimS2-2,1)])); 
    
%%

Fidelity=zeros(length(Asynchronicity),1);

for iter=1:length(Asynchronicity)
iter/length(Asynchronicity)

    Async=Asynchronicity(iter);
    gammaS2=(t-Async>=t0+DelayS1S2).*(t-Async<= 2*t0+DelayS1S2).*exp(-gamma_m*(t-Async-t0-DelayS1S2))./(2-exp(-gamma_m*(t-Async-t0-DelayS1S2))) + (t-Async<t0+DelayS1S2).*(t-Async>= DelayS1S2)*gamma_m;

         [Gamma,Lambda, Exc1, Exc2,rhoS1,rhoS2, rho_final]...
            = FullEvolution(Gamma0,Lambda0,Gammain, Lambdain, gammaS1, gammaS2, N_wav, dimB, dimS1, dimS2, dt, Nt, MS1S2, maxSchmidtrank, Eps, showprogress);

        if showtransfer==1

            figure
            plot((dt:dt:Totaltime), Exc1, '-b','linewidth', 2)
            hold on
            plot((dt:dt:Totaltime), Exc2, '-r','linewidth', 2)
            plot(t, gammaS1/max(gammaS1), '--b','markersize', 10)
            plot(t, gammaS2/max(gammaS1), '--r', 'markersize', 10)
            xlabel('\Gamma_{max}t')
            axis([0 Totaltime 0 1.05])
            title('State transfer')
            set(gca, 'fontsize', 20)
            legend('|e_1|^2', '|e_2|^2', '\gamma_1(t)/\Gamma_{max}', '\gamma_2(t)/\Gamma_{max}')

        end  
        
Fidelity(iter)= Psi_init'*rho_final*Psi_init;     
    
end    
%%
    
delete(mypool)

if saveresults==1
    save([datestr(datetime('now')) '.mat'], 'Asynchronicity','Nth', 'Totaltime', 'gamma_m', 't0', 'Fidelity')
end

