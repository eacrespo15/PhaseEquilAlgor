function[Betas,XPhase,NPhases,rho_phase,error,niter]=MF_FLASH_HEIDEMANN(P,T,zfeed,saftparam)
%Multiphase flash algorithm based on the works of R. Heidemann - Univ. Calgary
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Objective:Multiphase flash algorithm
%
%Reference: The works of Robert R. Heidemann  - Univ. Calgary
%
%Input parameters:
%T-         Temperature (K)
%P-         Pressure (MPa)
%zfeed-     Feed composition (NCx1)
%
%Output parameters:
%Betas      - Phase fractions (NFx1)
%XPhase     - Phase compositions (NFxNC)
%NPhases    - Number of stable phases in equilibrium
%rhophase   - Density of the phases in equilibrium (NFx1)

NC=length(zfeed); %Number of components
NF=NC+1;          %Start the algorithm with NC pure liquid phases and a vapour phase
error=1;          %Initialize convergence criterium
tol=1e-10;        %Tolerance for convergence
niter=0;          %Initialize iteration counter
maxiter=500;      %Maximum number of iterations

%Initial estimates of phase fractions
guess_beta=zeros(NF,1);
guess_beta(1:NF-1)=zfeed;
guess_beta=guess_beta*(NF-1)/NF;
guess_beta(NC+1)=1;
for i=1:NC
    guess_beta(NC+1)=guess_beta(NC+1)-guess_beta(i);
end

%Initial estimates for mole fractions. Vapor phase composition from mass balance
guess_x=zeros(NF,NC);
for i=1:NC
    for j=1:NF-1
        if i==j
            guess_x(j,i)=1;
            guess_x(NF,i)=guess_x(NF,i)-guess_x(j,i)*guess_beta(j);
        end
    end
    guess_x(NF,i)=(zfeed(i)+guess_x(NF,i))/guess_beta(NF);
end

%Mass balance on each phase
sumx=zeros(NF,1);
for j=1:NF
    for i=1:NC
        sumx(j)=sumx(j)+guess_x(j,i);
    end
end
if any(sumx>=1.001) ||  any(sumx<=0.999) 
    error('Sum of initial mole fractions in at least one of the phases different than 1');
end

while error>tol && niter<maxiter
    niter=niter+1;
    guess_x_old=guess_x;
    
    %Reference phase is the one with maximum beta
    ref=1;
    for j=1:NF
        if guess_beta(j)/guess_beta(ref)>1
            ref=j;
        end
    end
    
    %Determination of fugacities
    fi=zeros(NF,NC);
    rho_phase=zeros(NF,1);
    for j=1:NF
        zphase=guess_x(j,:);
        [f,rho_phase(j,1)]=simplex_fug_saft(zphase,T,P,0,saftparam);
        fug=exp(f);
        for i=1:NC
            fi(j,i)=fug(i);
        end
    end
    
    %Distribution coefficients
    K=zeros(NF,NC);
    for j=1:NF
        for i=1:NC
            K(j,i)=fi(ref,i)/fi(j,i);
        end
    end
    
    V=1./rho_phase;
    V=V/1000;
    
    %Call the stability analysis
    [NF,guess_beta,K]=chkphs(NF,NC,guess_beta,K,V);
    
    %Call the routine for material balances of the flash: determination of
    %phase amounts and phase compositions
    [guess_beta,guess_x]=fracts(NF,NC,zfeed,guess_beta,K);
    
    %Evaluation of the convergence criterium: Difference in compositions
    error=0;
    for j=1:NF
        for i=1:NC
            error=error+abs(guess_x(j,i)-guess_x_old(j,i));
        end
    end
    Betas=guess_beta;
    XPhase=guess_x;
    NPhases=NF;
    rho_phase=rho_phase';
end

