function [XPhase,rho_phase,delta,niter]=LLE_Binary(P,T,guess_ini,phasetype,saftparam)
%Binary LLE calculations using the isofugacity criterion
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%NF=2;NC=2
%Inputs:
%P - Pressure (MPa)
%T - Temperature (K)
%guess_ini - Initial estimate for phase compositions (NFXNC)
%phasetype - Type of phases to look for (1-Liquid; 2- Vapour;0- Auto)(NFx1)
%saftparam - EoS parameters
%
%Outputs:
%XPhase     - Composition of the Two-Phases in equilibrium (NF x NC)
%rho_phase  - Density of the phases in equilibrium (NFx1)
%delta      - Final value of the convergence criterium
%niter      - Final number of iterations
%
%Objective:Algorithm for binary LLE calculations. Can be easily adapted to 
%VLE or activity coefficient models byreplacing fi by the correspondent F1 and F2 functions.
%
%Reference:R. Privat et al. / Computers and Chemical Engineering 50 (2013) 139– 151

%1)Set the iteration counter to 0 and tol to 1e-6
niter=0;
tol=1e-6;
delta=1;
guess_x=guess_ini;
rho_phase=zeros(2,1);
fi=zeros(2,2);

while delta>tol && niter<1000
%1) Update iteration counter
niter=niter+1;    

%2) Provide initial estimates for the unknowns
zphase=[guess_x(1) 1-guess_x(1); guess_x(2) 1-guess_x(2)];

%3) Determine the fugacity coefficients which are the Fi functions for
%Phi-Phi approach
for j=1:2
    [f,rho_phase(j,1)]=simplex_fug_saft(zphase(j,:),T,P,phasetype(j),saftparam);
    fug=exp(f);
    fi(j,:)=fug;
end

%4) Define a convergence criterion delta:
delta=abs(zphase(2,1)*fi(2,1)-zphase(1,1)*fi(1,1))+abs(zphase(2,2)*fi(2,2)-zphase(1,2)*fi(1,2));

%5) Determine A and B
A=[fi(1,1) -fi(2,1);-fi(1,2) fi(2,2)];
b=[0;fi(2,2)-fi(1,2)];

%6) Update the unknowns
guess_x=A\b;
end

XPhase=zphase;
if niter==1000
    disp('Warning: The method did not converge');
    XPhase=NaN;
end
end