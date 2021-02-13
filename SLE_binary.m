function[XSLE,fval,exitflag,ACTCOEF,rhomix]=SLE_binary(P,T,xguess,index,melting,saftparam)
%Algorithm for binary SLE calculations
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Objective:Algorithm for the solid-liquid equilibria of binary mixtures
%For a given temperature determines the solute solubility in the liquid phase
%
%Outputs:
%
%XSLE       - solute solubility (mol/mol)
%fval       - Objective function value
%exitflag   - Convergence status
%ACTCOEF    - Activity coefficients (1 X 2)
%rhomix     - Mixture's density (mol/L)
%
%Inputs:
%P          - System's Pressure (MPa)
%T          - Temperature (K)
%xguess     - Mole fraction of the solute (guess)
%index      - Freeze out component index
%MELTING    - Melting properties previously read from input file [Tm;Hm;dCpm,Ntrs,Ttrs,Htrs]

%1) Define options settings for the fzero subroutine
options=optimset('Display','off','MaxIter',100);

%2) Define the melting properties
Tm=melting.Tm(index);
Hm=melting.Hm(index);
Cpm=melting.dCp(index);
NT=melting.Ntrs(index);
if NT>0
    Tt=melting.Ttrs(1:NT,index);
    Ht=melting.Htrs(1:NT,index);
else
    Tt=0;
    Ht=0;
end

%3) Define a function handle
fsolubility=@(x)find_solubility(x,T,P,Tm,Hm,Cpm,NT,Tt,Ht,index,saftparam);

%4) Solve the SLE equation for the zero
[XSLE,~,fval,exitflag]=lsqnonlin(fsolubility,xguess,0,1,options);
% A possible improvement would be to define the maximum x to the previous x
% on the following points

%5) Calculate the fugacity coefficients in the liquid mixture
if index==2
    zmix=[1-XSLE XSLE];
else
    zmix=[XSLE 1-XSLE];
end
[LFUG_Mix,rhomix]=simplex_fug_saft(zmix,T,P,1,saftparam);

%6) Calculate the fugacity coefficients in the pure compounds
[LFUG_Pur1,~]=simplex_fug_saft([1 0],T,P,1,saftparam);
[LFUG_Pur2,~]=simplex_fug_saft([0 1],T,P,1,saftparam);
LFUG_Pure=[LFUG_Pur1(1) LFUG_Pur2(2)];
%Calculate the activity coefficients in the liquid phase
ACTCOEF=exp(LFUG_Mix-LFUG_Pure);
end

function [fobj]=find_solubility(x,T,P,Tm,Hm,Cpm,NT,Tt,Ht,id,saftparam)
%Define R
R=8.3144;

%Calculate psi(k)
psi=Hm*(T-Tm)/(R*Tm*T)+Cpm/R*(log(T/Tm)-(T-Tm)/T);

%Accounting for solid-solid transitions
if NT>0
    for m=1:NT
        if T<Tt(m)
            psi=psi+Ht(m)*(T-Tt(m))/(R*Tt(m)*T);
        end
    end
end

%Calculate the fugacity coefficient in the liquid mixture
if id==2
    zmix=[1-x x];
else
    zmix=[x 1-x];
end

%Calculate the fugacity coefficients in the mixture
[LFUG_Mix,~]=simplex_fug_saft(zmix,T,P,1,saftparam);

%Calculate the fugacity coefficients in the pure compounds
[LFUG_Pur1,~]=simplex_fug_saft([1 0],T,P,1,saftparam);
[LFUG_Pur2,~]=simplex_fug_saft([0 1],T,P,1,saftparam);
LFUG_Pure=[LFUG_Pur1(1) LFUG_Pur2(2)];

%Calculate the activity coefficients in the liquid phase
ACTCOEF=exp(LFUG_Mix-LFUG_Pure);

%Calculate the objective function
fobj=log(x*ACTCOEF(id))-psi;
end
