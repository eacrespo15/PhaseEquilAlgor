function[Betas,XPhase,rho_phase,niter,time]=MF_FLASH_MICHELSEN(P,T,zfeed,guess_x,guess_beta,phasetype,saftparam)
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Objective:Multiphase flash algorithm
%
%Reference: Notes taken in DTU PhD. Summer course on Thermodynamic models;
%Michelsen & Mollerup - Thermodynamic models - Fundamentals and
%Computational Aspects
%
%Input Variables:
%P          - Pressure/MPa
%T          - Temperature/K
%zfeed      - Feed molar compositions (NCx1)
%guess_x    - Guess of phase compositions (NFxNC)
%guess_beta - Guess of phase fractions (NFx1)
%phasetype  - Type of phases expected (NFx1) (0) - Most stable (1) liquid
%(2) vapour
%saftparam  - EoS parameters

t=cputime;                  %Start timer
NC=length(zfeed);           %Number of components
NF=length(guess_beta);      %Number of Phases (assumed)
fi=zeros(NF,NC);            %Initialize the matrix of fugacity coefs.
rho_phase=zeros(1,NF);      %Initialize the matrix for phase densities

delta_betas=1;
niter=0;
while delta_betas>1e-6 && niter<100
    niter=niter+1;
    for j=1:NF
        [f,rho_phase(1,j)]=simplex_fug_saft(guess_x(j,:),T,P,phasetype(j),saftparam);
        fug=exp(f);
        fi(j,:)=fug;
    end
    [Betas,guess_x]=MPh_Rachford_Rice(NF,zfeed,fi,guess_beta);
    delta_betas=norm(Betas-guess_beta);
    guess_beta=Betas;
end
XPhase=guess_x;

time=t+cputime;

end
