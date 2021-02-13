function[output_variable,y,niter,rhol,rhov]=bubble_point(given_variable,z,spec,guess,critical,saftparam)
%Algorithm for Bubble Point calculations without using explicit derivatives of LnFUG
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%References: Notes taken during the DTU PhD. Summer course on Thermodynamic models
%Michelsen & Mollerup - Thermodynamic models: Fundamentals and Computational Aspects
%
%Input variables:
%given_variable - Fixed T(K)for Bubble P calculations or fixed P(MPa) for Bubble T calculations
%z              - Feed composition in molar fractions [x1,x2,...,xN]
%spec           - (1) - Bubble T calculation (2)- Bubble P calculation
%guess          - Guess for the unknown variable. If spec=1 a T_guess must be given; if spec=2 a P_guess must be given
%critical       - Critical Properties (structure read previously from the input file)
%
%Output Variables:
%output_variable    - Bubble Point Pressure/MPa (spec=2) or Temperature/K(spec=1)
%y                  - Vapour Phase composition (1xNC)
%niter              - Number of iterations until convergence
%rhol               - Density of the liquid phase (mol/L)
%rhov               - Density of the vapour phase (mol/L)


%Load the critical properties and initialize K-values
K=zeros(1,length(z));
Tc=critical.Tc;
Pc=critical.Pc;
omega=critical.w;

if spec==2
    
    %Spec==2 Means Bubble Pressure calculation (T is given)
    T=given_variable;
    
    %Calculate initial estimate for Pressure if guess is not given
    P=0;
    if ~exist('guess','var') || isempty(guess)
        for i=1:length(z)
            P=P+z(i)*Pc(i)*exp(5.373*(1+omega(i))*(1-Tc(i)/T));
        end
    else
        P=guess;
    end
    
    %Calculate startup K-Factors (through Wilson)
    for i=1:length(z)
        K(i)=Pc(i)/P*exp(5.373*(1+omega(i))*(1-Tc(i)/T));
    end
    
    %Initialize variable to converge (sum of vapour compositions=1) and iteration counter
    sumy=2;
    niter=0;
    yold=K.*z;
    y=yold;
    
    %Cycle of Convergence
    while abs(sumy-1)>1e-6 && niter<1000
        niter=niter+1;
        niter_inner=0;
        dy=ones(length(y),1);
        while norm(dy)>1e-6 && niter_inner<100
            niter_inner=niter_inner+1;
            [LCFUGL,rhol]=simplex_fug_saft(z,T,P,1,saftparam);
            [LCFUGV,rhov]=simplex_fug_saft(y,T,P,2,saftparam);
            K=exp(LCFUGL-LCFUGV);
            y=z.*K;
            dy=y-yold;
            yold=y;
            y=y/sum(y);
        end
        sumy=sum(yold);
        tempy(niter)=sumy;
        tempp(niter)=P;
        if niter<=2
            if sumy>1
                P=P*1.05;
            else
                P=P*0.95;
            end
        else
            Pnew=tempp(niter-2)+(1-tempy(niter-2))*(tempp(niter-1)-tempp(niter-2))/(tempy(niter-1)-tempy(niter-2));
            if Pnew<0
                P=0.75*P;
            else
                P=Pnew;
            end
        end
    end
    output_variable=P;
elseif spec==1
    
    %Spec==1 Means Bubble Temperature calculation (P is given)
    P=given_variable;
    
    %Calculate initial estimate for Temperature and K-Factors
    if ~exist('guess','var') || isempty(guess)
        T=350;
        dT=1;
        niter=0;
        while abs(dT)>1e-6 && niter<100
            for i=1:length(z)
                K(i)=Pc(i)/P*exp(5.373*(1+omega(i))*(1-Tc(i)/T));
            end
            f=z*K'-1;
            dfdt=1/T^2*sum(5.373*(1+omega).*Tc.*K.*z);
            dT=-f/dfdt;
            T=T+dT;
            niter=niter+1;
        end
    else
        T=guess;
        for i=1:length(z)
            K(i)=Pc(i)/P*exp(5.373*(1+omega(i))*(1-Tc(i)/T));
        end
    end
    
    %Initialize variable to converge (sum of vapour compositions=1) and iteration counter 
    sumy=2;
    niter=0;
    yold=K.*z;
    y=yold;
    
    %Cycle of Convergence
    while abs(sumy-1)>1e-6 && niter<1000
        niter=niter+1;
        niter_inner=0;
        dy=ones(length(y),1);
        while norm(dy)>1e-6 && niter_inner<100
            niter_inner=niter_inner+1;
            [LCFUGL,rhol]=simplex_fug_saft(z,T,P,1,saftparam);
            [LCFUGV,rhov]=simplex_fug_saft(y,T,P,2,saftparam);
            K=exp(LCFUGL-LCFUGV);
            y=z.*K;
            dy=y-yold;
            yold=y;
            y=y/sum(y);
        end
        sumy=sum(yold);
        tempy(niter)=sumy;
        tempt(niter)=T;
        if niter<=2
            if sumy<1
                T=T*1.05;
            else
                T=T*0.95;
            end
        else
            T=tempt(niter-2)+(1-tempy(niter-2))*(tempt(niter-1)-tempt(niter-2))/(tempy(niter-1)-tempy(niter-2));
        end
    end
    output_variable=T;
else
    error('Invalid Specification for Bubble Point Calculation');
end
end


