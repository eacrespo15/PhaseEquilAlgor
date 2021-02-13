function[viscosity]=SSViscosity(T,P,z,rho,critical,viscoparam,dilute_term,mrule)
%Viscosity calculations using the Free Volume Theory
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Objective:Viscosity calculations by coupling the Free Volume theory to the soft-SAFT EoS
%
%References: 
%
% Input variables:
% T-Temperature
% P-Pressure
% z-Mixture Composition (1xNC)
% rho- System's density =f(T,P,z) previously calculated with soft-SAFT (mol/L)
% critical - critical properties necessary for the diluted term and Mw
% (structure read from the input file)
% viscoparam - parameters for FVT theory (structure read from the input file)
% diluted_term - Enable or disable the diluted term
% mrule - set of mixing rules for the dense term (1) - Linear (2)- Polishuk and
% Yitzhak (IECR 53,2014,959-971) (3)- Baylaucq et al.
%
% Output variables:
% viscosity - Viscosity in mPa.s /cP

%Read pure-component properties
Mw=critical.Mw;
TC=critical.Tc;
VC=critical.Vc;
OMEGA=critical.w;

%Read viscosity parameters
alfa=viscoparam.alfa;
B=viscoparam.B;
Lv=viscoparam.Lv;
dipolerorig=viscoparam.dipole;
kappa=viscoparam.kappa;

%Load R value and initiliaze variables
physics=module_physics;
Bmix=0;
Lvmix=0;
alfamix=0;
Mwmix=0;
NC=length(z);
viscosity=0;

%Reduce the dipole moment
dipoler=zeros(1,NC);
for i=1:NC
    dipoler(i)=(131.3*dipolerorig(i))/sqrt(1000*VC(i)*TC(i));
end

%Dilute Gas Viscosity (Ind. Eng. Chem. Res., Vol 27, No.4, 1988)
if dilute_term==1
   Fc=zeros(1,NC);
   collision=zeros(1,NC);
   N0=zeros(1,NC);
   N0mix=0;
   for i=1:NC
       Fc(i)=1-0.2756*OMEGA(i)+0.059035*power(dipoler(i),4)+kappa(i);
       tvstar=1.2593*T/TC(i);
       collision(i)=(1.16145/power(tvstar,0.14874))+0.52487/exp(0.77320*tvstar)+2.16178/exp(2.43787*tvstar)-0.0006435*power(tvstar,0.14874)*sin(18.0323*power(tvstar,-0.76830)-7.27371);
       N0(i)=sqrt(Mw(i)*T)*Fc(i)*0.0040785/(collision(i)*power(VC(i)*1000,2/3));
   end
   %Sum the contribution of each component to the diluted term
   %(https://doi.org/10.1063/1.1747673)
   phi=zeros(NC,NC);
   for i=1:NC
       sum1=0;
       for j=1:NC
           phi(i,j)=power(1+power(N0(i)/N0(j),0.5)*power(Mw(j)/Mw(i),1/4),2)/(4/sqrt(2)*power(1+Mw(i)/Mw(j),0.5));
           if i~=j
               sum1=sum1+z(j)*phi(i,j);
           end
       end
       N0mix=N0mix+N0(i)/(1+sum1/z(i));
   end
   viscosity=viscosity+N0mix;
end

% Dense term % (J. Phys. Chem. B 2013, 117,8159-8171)
if mrule==1 %Linear mixing Rule
    for i=1:NC
        Mwmix=Mwmix+Mw(i)*z(i);
        Bmix=Bmix+B(i)*z(i);
        Lvmix=Lvmix+Lv(i)*z(i);
        alfamix=alfamix+alfa(i)*z(i);
    end
    dN=Lvmix*(0.1*P+0.0001*alfamix*power(rho,2)*Mwmix)*sqrt(0.001*Mwmix/(3*physics.R*T))*exp(Bmix*power((1000*P+alfamix*power(rho,2)*Mwmix)/(rho*physics.R*T),3/2));
elseif mrule==2 %Mixing rule of Polishuk and Yitzhak (IECR 53,2014,959-971)
    dummyB=0;
    for i=1:NC
        Mwmix=Mwmix+Mw(i)*z(i);
        alfamix=alfamix+alfa(i)*z(i);
        Lvmix=Lvmix+Lv(i)*z(i);
        dummyB=dummyB+z(i)/B(i);
    end
    Bmix=1/dummyB;
    dN=Lvmix*(0.1*P+0.0001*alfamix*power(rho,2)*Mwmix)*sqrt(0.001*Mwmix/(3*physics.R*T))*exp(Bmix*power((1000*P+alfamix*power(rho,2)*Mwmix)/(rho*physics.R*T),3/2));
elseif mrule==3 %Mixing rule of Baylaucq
    fvi=zeros(1,NC);
    dummyB=0;
    dummyfv=0;
    for i=1:NC
        fvi(i)=power((rho*physics.R*T)/(alfa(i)*power(rho,2)*Mw(i)+1000*P),3/2);
        Mwmix=Mwmix+Mw(i)*z(i);
        Lvmix=Lvmix+Lv(i)*z(i);
        dummyB=dummyB+z(i)/B(i);
        dummyfv=dummyfv+z(i)/fvi(i);
    end
    Bmix=1/dummyB;
    fvmix=1/dummyfv;
    dN=0.0001*rho*Lvmix*sqrt(0.001*Mwmix*physics.R*T/3)*power(fvmix,-2/3)*exp(Bmix/fvmix);
end

% Summation of the two contributions
viscosity=viscosity+dN;
end


