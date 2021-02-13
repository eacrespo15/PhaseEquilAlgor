function [tml,tmv,LCFUGWV,LCFUGWL,LCFUGZ]=stability_analysis(T,P,z,K,saftparam)
%Stability analysis for the Two-Phases PT-Flash
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Objective:Performs the stability analysis for the two phases PT-Flash
%
%Reference: Notes taken in DTU PhD. Summer course on Thermodynamic models;
%Michelsen & Mollerup - Thermodynamic models - Fundamentals and
%Computational Aspects

[LCFUGZ]=simplex_fug_saft(z,T,P,0,saftparam);
Di=LCFUGZ+log(z);

WV=K.*z; %Vapour-Like estimate
WL=z./K; %Liquid-Like estimate

WL=WL/sum(WL);
WV=WV/sum(WV);

criterio=1;
iteracao_sl=0;
while max(abs(criterio))>1e-9 && iteracao_sl<200
    iteracao_sl=iteracao_sl+1;
    [LCFUGWL]=simplex_fug_saft(WL,T,P,0,saftparam);
    LWL_new=Di-LCFUGWL;
    criterio=exp(LWL_new)-WL;
    WL=exp(LWL_new);
    WL=WL/sum(WL);
end
tml=1+sum(WL.*(log(WL)+LCFUGWL-Di-1));

criterio=1;
iteracao_sv=0;
while max(abs(criterio))>1e-9 && iteracao_sv<200
    iteracao_sv=iteracao_sv+1;
    [LCFUGWV]=simplex_fug_saft(WV,T,P,0,saftparam);
    LWV_new=Di-LCFUGWV;
    criterio=exp(LWV_new)-WV;
    WV=exp(LWV_new);
    WV=WV/sum(WV);
end
tmv=1+sum(WV.*(log(WV)+LCFUGWV-Di-1));