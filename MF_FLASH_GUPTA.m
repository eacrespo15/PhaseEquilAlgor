function [Beta,Thita,Xfinal,K,phasetype,rho_phase,val,time]=MF_FLASH_GUPTA(P,T,zfeed,guess_x_initial,guess_beta,phasetype,saftparam)
%Multiphase flash algorithm by Gupta et al.
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Objective: Given the feed composition zfeed, temperature T and pressure P, and a set of possible
%phases defined by guess_x, guess_beta and phase type, performs a multiphase 
%calculation to determine which phases are stable, phase fractions and compositions
%
%References: Martin et al. 2011 - 10.1016/j.ece.2011.08.003
%Gupta et al., Fluid Phase Equilibria, 63 (1991), 65-89
%
%Input Parameters:
%P: Pressure (MPa)
%T: Temperature (K)
%zfeed: Feed composition (1xNC)
%guess_x: Initial guess for phase compositions (NF x NC)
%guess_beta: Initial guess for phase fraction (NF x 1)
%phase: Type of phase expected. (NFx1) e.g. [1 1 2] for liquid liquid vapor
%
%Output Parameters:
%Beta:  Calculated phase fractions (NFx1)
%Thita: Calculated phase stability parameters (NFx1)
%x:     Calculated phase compositions (NF x NC)
%K:     Distribution coefficients (NF x NC) K(Phase i, component j)= x(Phase i, component j) / x(phase 1,component j)
%phasetype: phase order may change from the initial input
%rho_phase: Density of the different phases (NF x 1) / mol.L^-1
%val:   final value of the objective function
%time:  time required for the calculations




t=cputime;                      %Start timer
NF=length(guess_beta);          %Number of Phases
NC=length(zfeed);               %Number of Components

guess_x=guess_x_initial;
% Initial guess for XBeta
Xini=guess_beta(2:NF);
for i=2:NF
    if Xini (i-1)==0
        Xini(i-1)=-1;
    end
end

%Options for fsolve of objective function
OPTIONS=optimset('Display','off','MaxIter',200,'MaxFunEvals',200,'TolX',1e-9,'TolFun',1e-6);

%Tolerance and maximum number of iterations
tol=1e-4;
val=1;
iter=0; 
maxiter=100;
epsilon=1e-6;


% While cycle that continues iterations until differences between calculated compositions in two 
%successive iterations is lower than tol,or until the maximum number of iterations is reached

while val>tol && iter<maxiter
    
    iter=iter+1;
    guess_x_old=guess_x;
    
    %Determination of fugacities
    fi=zeros(NF,NC);
    rho_phase=zeros(NF,1);
    for j=1:NF
        zphase=guess_x(j,:);
        [f,rho_phase(j,1)]=simplex_fug_saft(zphase,T,P,phasetype(j),saftparam);
        fug=exp(f);
        for i=1:NC
            fi(j,i)=fug(i);
        end
    end
    
    %Distribution coefficients
    K=zeros(NF,NC);
    for j=1:NF
        for i=1:NC
            K(j,i)=fi(1,i)/fi(j,i);
        end
    end
    
    %Minimization of Objective Function
    XBeta=fsolve(@(XBeta) obj_MF_Flash_GUPTA(XBeta,zfeed,K,NF,NC),Xini,OPTIONS);
    
    %Store phase fractions and stability results
    Beta=zeros(NF,1);
    Thita=zeros(NF,1);
    
    for j=2:NF
        if XBeta(j-1)<0
            Beta(j)=0;
            Thita(j)=-XBeta(j-1);
        else
            Beta(j)=XBeta(j-1);
            Thita(j)=0;
        end
    end
    
    %First phase is the reference phase so Beta should be >0 and thita=0
    Thita(1)=0;
    Beta(1)=1-sum(Beta);
    Beta=abs(Beta);
    
    %Rounds Beta and thita values smaller than epsilon
    for i=1:NF
        if Beta(i)<=epsilon
            Beta(i)=0;
        end
        if Thita(i)<=epsilon
            Thita(i)=0;
        end
    end
    
    %Normalization of Beta values
    Beta=Beta/sum(Beta);
    
    %Prepare initial guess for following iteration
    for j=2:NF
        if Thita(j)>0
            Xini(j-1)=-Thita(j);
        else
            Xini(j-1)=Beta(j);
        end
    end
    
    %Updating compositions
    for i=1:NC
        sumat=0;
        for j=2:NF
            sumat=sumat+(K(j,i)*exp(Thita(j))-1)*Beta(j); 
        end
        guess_x(1,i)=zfeed(i)/(1+sumat); %Composition of reference phase
    end
    
    for i=1:NC
        for j=2:NF
            guess_x(j,i)=guess_x(1,i)*K(j,i)*exp(Thita(j)); %Composition of remaining phases
        end
    end
    
    for j=1:NF
        sumat=sum(abs(guess_x(j,:)));
        for i=1:NC
            guess_x(j,i)=abs(guess_x(j,i))/sumat; %Normalizing the compositions to 1
        end
    end
    
    %Evaluate difference between successive iterations
    val=0;
    for j=1:NF
        for i=1:NC
            val=val+abs(guess_x(j,i)-guess_x_old(j,i));
        end
    end
    
    Xfinal=guess_x; %Return the calculated compositions
end

if (iter == maxiter)
    s = sprintf('%f',val);
    warning('MATLAB:EoS', ['Convergence error in MultiFlash. Maximum number of iterations reached. Final value of objective function: ' s '.']);
end


%Shift the phases order to try for convergence
if val>tol
    order=input('Do you want to try shifting the order of the phases? Yes (Y) No (N)','s');
    if order=='Y'
        shift=1;
        guess_beta_shift=zeros(1,NF);
        guess_x_shift=zeros(NF,NC);
        while val>tol && shift<NF
            for j=1:NF
                j_shift=j+shift;
                if j_shift>NF
                    j_shift=j_shift-NF;
                end
                guess_beta_shift(j_shift)=guess_beta(j);
                phasetype(j_shift)=phasetype(j);
                for i=1:NC
                    guess_x_shift(j_shift,i)=guess_x_initial(j,i);
                end
            end
            [Beta_shift,Thita_shift,Xfinal_shift,K_shift,phasetype_shift,rho_phase_shift,val_shift,~]=MF_FLASH_GUPTA(P,T,zfeed,guess_x_shift,guess_beta_shift,phasetype,saftparam);
            shift=shift+1;
        end
        if val_shift <tol
            Beta=Beta_shift;
            Thita=Thita_shift;
            Xfinal=Xfinal_shift;
            K=K_shift;
            phasetype=phasetype_shift;
            rho_phase=rho_phase_shift;
            val=val_shift;
        end
    end
end


%Merges multiple occurencies of identical phases
for j = 1:NF-1
    if Thita(j) == 0
        for k = j+1:NF
            if (abs(K(k,1) - K(j,1)) < 1e-3)
                Beta(j) = Beta(j) + Beta(k);
                Beta(k) = 0;
            end
        end
    end
end


%Transforms K into format used by other programs
for j = 1:NF
    for i = 1:NC
        K(j,i) = 1/K(j,i);
    end
end

Beta=Beta';
Thita=Thita';
rho_phase=rho_phase';

time=cputime-t;

end