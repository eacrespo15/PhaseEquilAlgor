function [Beta,x,y,iter_total,iter_acel,iter_extrap]=VLE_PTFLASH(T,P,Kguess,z,approach,check_stab,rrsolver_type,saftparam)
%Two-Phase PT Flash
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Objective:Two-Phase PT FLASH With/Without Acceleration through DEVM
%DEVM - Dominant Eigen Value Method
%
%Reference: Notes taken in DTU PhD. Summer course on Thermodynamic models;
%Michelsen & Mollerup - Thermodynamic models - Fundamentals and
%Computational Aspects
%
% Output variables:
% Beta          - Ratio between V and F
% x             - Liquid composition (NCx1)
% y             - Vapour composition (NCx1)
% iter_total    - Total number of iterations
% iter_acel     - Number of accelerated steps
% iter_extrap   - Number of extrapolated steps
%
% Input variables:
% T                 - Temperature (K)
% P                 - Pressure (MPa)
% Kguess            - K-factors guess (NCx1)
% z                 - Feed composition (NCx1)
% approach          - Acceleration method: (0)- Unaccelerated Flash (1)- After every 5 steps of Succesive Substition (2)- When the eigen value is changing only within a tolerance
% check_stab        - Stability analysis:  (0)- Disabled (1)- Enabled (2) - Full Negative Flash (Autom. enabled if stability analysis suggests one phase)
% rrsolver_type     - Solver of Rachford-Rice Equations: (0)-Full RRSolver (1)-Negative RRSolver (2)-Sloppy RRSolver


%% Initialization of variables
iteracao=0;     %Counter of iterations to perform acceleration
iter_total=0;   %Counter of the total number of iterations
iter_acel=0;    %Counter of acceleration trials
iter_extrap=0;  %Counter of successful extrapolations
Beta=0.5;       %Initial Beta=0.5
delta_LK=1;     %Initialization of the convergence criterium
K=Kguess;       %K_guess are initial guesses for K (e.g. Wilson K-Factors)

%% Check for the optional arguments
%a) Acceleration Approach:
%0- Unaccelerated Flash
%1- After every 5 steps of Succesive Substition
%2- When the eigen value is changing only within a tolerance
if ~exist('approach','var') || isempty(approach)
    approach=0;
end
%b) Stability Analysis:
%0-Disabled
%1-Enabled
%2-Full Negative Flash (Enabled automatically if stability analysis suggests one-phases)
if ~exist('check_stab','var') || isempty(check_stab)
    check_stab=1;
end
%c) RRSolver:
%0-Full RRSolver
%1-Negative RRSolver
%2-Sloppy RRSolver
if ~exist('rrsolver_type','var') || isempty(rrsolver_type)
    rrsolver_type=0;
end

%% STABILITY ANALYSIS
if check_stab==1
    [tml,tmv,CFUGWV,CFUGWL,CFUGZ]=stability_analysis(T,P,z,K,saftparam);
    if ~(abs(tmv)<1e-5 || abs(tml)<1e-5) && abs(tml-tmv)>1e-5
        % if any tm<0 then we are in the two-phase region
        if tml<0 && tmv<0
            K=exp(CFUGWL-CFUGWV);
        elseif tml<0
            K=exp(CFUGWL-CFUGZ);
        elseif tmv<0
            K=exp(CFUGZ-CFUGWV);
        else
            %Perform a Full-Negative Flash
            check_stab=2;
        end
    else
        %Perform a Full-Negative Flash
        check_stab=2;
    end
end

%% FLASH CALCULATIONS
if check_stab==1 || check_stab==0
    if approach==0
        while delta_LK>1e-7 && iter_total<500
            iter_total=iter_total+1;                        %Iteration Counter
            [Beta,x,y]=RR_solver(K,z,rrsolver_type);        %Call to RR_solver with actual K's
            [LCFUGL]=simplex_fug_saft(x,T,P,1,saftparam);     %Calculate Liquid Phase Fugacities
            [LCFUGV]=simplex_fug_saft(y,T,P,2,saftparam);     %Calculate Vapor Phase Fugacities
            LK(iter_total,:)=LCFUGL-LCFUGV;                 %Calculate Ln(K)
            K=exp(LK(end,:));                               %Update K-Factors
            if iter_total>1
                delta_LK=max(abs(LK(iter_total,:)-LK(iter_total-1,:))); %Update Convergence criterium
            end
        end
    elseif approach==1
        while delta_LK>1e-7 && iter_total<500
            
            %Successive Substitution Steps (cycles of 5 iterations)
            while delta_LK>1e-7 && iteracao<5
                iteracao=iteracao+1;
                iter_total=iter_total+1;
                [Beta,x,y]=RR_solver(K,z,rrsolver_type);    %Calls RR_solver for actual K-Factors
                [LCFUGL]=simplex_fug_saft(x,T,P,1,saftparam);     %Calculate Liquid Phase Fugacities
                [LCFUGV]=simplex_fug_saft(y,T,P,2,saftparam);     %Calculate Vapor Phase Fugacities
                LK(iter_total,:)=LCFUGL-LCFUGV;             %Calculate Ln(K)
                K=exp(LK(iter_total,:));                    %Update K-Factors
                if iter_total>1
                    delta_LK=max(abs(LK(iter_total,:)-LK(iter_total-1,:))); %Update Convergence Criterium
                end
            end
            
            if delta_LK>1e-7
                G_SS=(1-Beta)*sum(x.*(LCFUGL+log(x)))+Beta*sum(y.*(LCFUGV+log(y))); %Calculate Gibbs energy at each 5th step
                iter_acel=iter_acel+1;                                              %Count one additional trial for extrapolation
                Dk=LK(end-1,:)-LK(end-2,:);
                Dk1=LK(end,:)-LK(end-1,:);
                lambda=sum(Dk1.^2)/sum(Dk.*Dk1);                                    % Eq.33 Ch.10 - I don't understand why it is not Dk1.*Dk/Dk.^2
                LNK=LK(end,:)+Dk1*lambda/(1-lambda);                                %Extrapolation
                K_ext=exp(LNK);                                                     %Calculate Extrapolated K-Factors
                [Beta,x,y]=RR_solver(K_ext,z);                                      %Call RR_Solver to obtain new compositions
                [LCFUGL]=simplex_fug_saft(x,T,P,1,saftparam);     %Calculate Liquid Phase Fugacities
                [LCFUGV]=simplex_fug_saft(y,T,P,2,saftparam);     %Calculate Vapor Phase Fugacities
                G_A=(1-Beta)*sum(x.*(LCFUGL+log(x)))+Beta*sum(y.*(LCFUGV+log(y)));  %Calculate Gibbs energy after accelerated step
                
                %If Gibbs energy is decreased accept extrapolated K_values
                if G_A-G_SS<0
                    K=K_ext;
                    iter_extrap=iter_extrap+1;                                  %Count a successful extrapolation
                    LK(iter_total,:)=LNK;
                    delta_LK=max(abs(LK(iter_total,:)-LK(iter_total-1,:)));     %Update Convergence Criterium
                end
            end
            iteracao=0;
        end
    elseif approach==2
        dfactor=2;
        while delta_LK>1e-7 && iter_total<500
            %Successive Substitution Steps (Until dfactor>1)
            while delta_LK>1e-7 && dfactor>1 && iter_total<500
                iteracao=iteracao+1;
                iter_total=iter_total+1;
                [Beta,x,y]=RR_solver(K,z,rrsolver_type);    %Calls RR_solver for actual K-Factors
                [LCFUGL]=simplex_fug_saft(x,T,P,1,saftparam);     %Calculate Liquid Phase Fugacities
                [LCFUGV]=simplex_fug_saft(y,T,P,2,saftparam);     %Calculate Vapor Phase Fugacities
                LK(iter_total,:)=LCFUGL-LCFUGV;             %Calculate Ln(K)
                K=exp(LK(iter_total,:));                    %Update K-Factors
                if iter_total>2
                    delta_LK=max(abs(LK(iter_total,:)-LK(iter_total-1,:)));                  %Update Convergence Criterium
                    G_SS=(1-Beta)*sum(x.*(LCFUGL+log(x)))+Beta*sum(y.*(LCFUGV+log(y)));      %Calculate Gibbs energy at each SS step
                    Dk=LK(end-1,:)-LK(end-2,:);
                    Dk1=LK(end,:)-LK(end-1,:);
                    lambda=sum(Dk1.^2)/sum(Dk.*Dk1);                                         % Eq.33 Ch.10
                    factor(iter_total-2)=lambda/(1-lambda);
                end
                if iter_total>3
                    dfactor=abs(factor(iter_total-2)-factor(iter_total-3))*100/abs(factor(iter_total-3));
                end
            end
            
            iter_acel=iter_acel+1;                      %Count one additional trial for extrapolation
            LNK=LK(end,:)+Dk1*lambda/(1-lambda);        %Extrapolation
            K_ext=exp(LNK);                             %Calculate Extrapolated K-Factors
            [Beta,x,y]=RR_solver(K_ext,z);              %Call RR_Solver to obtain new compositions
            [LCFUGL]=simplex_fug_saft(x,T,P,1,saftparam);     %Calculate Liquid Phase Fugacities
            [LCFUGV]=simplex_fug_saft(y,T,P,2,saftparam);     %Calculate Vapor Phase Fugacities
            G_A=(1-Beta)*sum(x.*(LCFUGL+log(x)))+Beta*sum(y.*(LCFUGV+log(y)));%Calculate Gibbs energy after accelerated step
            %If Gibbs energy is decreased accept extrapolated K_values
            if G_A-G_SS<0
                K=K_ext;
                iter_extrap=iter_extrap+1;              %Count a successful extrapolation
                LK(iter_total,:)=LNK;
                delta_LK=max(abs(LK(iter_total,:)-LK(iter_total-1,:))); %Update Convergence Criterium
            end
            dfactor=2;
        end
    else
        Error ('Invalid Approach selected for Acceleration in the Two-Phase PT-FLASH');
    end
elseif check_stab==2
    if approach==0
        while (delta_LK)>1e-7 && iter_total<500
            iter_total=iter_total+1;                    %Iteration Counter
            [Beta,x,y]=RR_solver(K,z,1);                %Call to RR_solver with actual K's
            [LCFUGL]=simplex_fug_saft(x,T,P,1,saftparam);     %Calculate Liquid Phase Fugacities
            [LCFUGV]=simplex_fug_saft(y,T,P,2,saftparam);     %Calculate Vapor Phase Fugacities
            LK(iter_total,:)=LCFUGL-LCFUGV;             %Calculate Ln(K)
            K=exp(LK(end,:));                           %Update K-Factors
            if iter_total>1
                delta_LK=max(abs(LK(iter_total,:)-LK(iter_total-1,:))); %Update Convergence criterium
            end
        end
    elseif approach==1
        while delta_LK>1e-7 && Beta>=0 && Beta<=1 && iter_total<500
            %Successive Substitution Steps (cycles of 5 iterations)
            while delta_LK>1e-7 && iteracao<5
                iteracao=iteracao+1;
                iter_total=iter_total+1;
                [Beta,x,y]=RR_solver(K,z,1);                %Calls RR_solver for actual K-Factors
                [LCFUGL]=simplex_fug_saft(x,T,P,1,saftparam);     %Calculate Liquid Phase Fugacities
                [LCFUGV]=simplex_fug_saft(y,T,P,2,saftparam);     %Calculate Vapor Phase Fugacities
                LK(iter_total,:)=LCFUGL-LCFUGV;             %Calculate Ln(K)
                K=exp(LK(iter_total,:));                    %Update K-Factors
                if iter_total>1
                    delta_LK=max(abs(LK(iter_total,:)-LK(iter_total-1,:))); %Update Convergence Criterium
                end
            end
            % Extrapolation step
            if delta_LK>1e-7 && Beta>=0 && Beta<=1
                G_SS=(1-Beta)*sum(x.*(LCFUGL+log(x)))+Beta*sum(y.*(LCFUGV+log(y))); %Calculate Gibbs energy at each 5th step
                iter_acel=iter_acel+1;                                              %Count one additional trial for extrapolation
                Dk=LK(end-1,:)-LK(end-2,:);
                Dk1=LK(end,:)-LK(end-1,:);
                lambda=sum(Dk1.^2)/sum(Dk.*Dk1);                                    % Eq.33 Ch.10
                LNK=LK(end,:)+Dk1*lambda/(1-lambda);                                %Extrapolation
                K_ext=exp(LNK);                                                     %Calculate Extrapolated K-Factors
                [Beta,x,y]=RR_solver(K_ext,z,1);                                    %Call RR_Solver to obtain new compositions
                [LCFUGL]=simplex_fug_saft(x,T,P,1,saftparam);     %Calculate Liquid Phase Fugacities
                [LCFUGV]=simplex_fug_saft(y,T,P,2,saftparam);     %Calculate Vapor Phase Fugacities
                G_A=(1-Beta)*sum(x.*(LCFUGL+log(x)))+Beta*sum(y.*(LCFUGV+log(y)));  %Calculate Gibbs energy after accelerated step
                %If Gibbs energy is decreased accept extrapolated K_values
                if G_A-G_SS<0
                    K=K_ext;
                    iter_extrap=iter_extrap+1;                                      %Count a successful extrapolation
                    LK(iter_total,:)=LNK;
                    delta_LK=max(abs(LK(iter_total,:)-LK(iter_total-1,:)));         %Update Convergence Criterium
                end
            end
            iteracao=0;
        end
    elseif approach==2
        dfactor=2;
        while delta_LK>1e-7 && Beta>=0 && Beta<=1 && iter_total<500
            
            %Successive Substitution Steps (Until dfactor>1)
            while delta_LK>1e-7 && dfactor>1 && Beta>=0 && Beta<=1 && iter_total<500
                iter_total=iter_total+1;
                [Beta,x,y]=RR_solver(K,z,1);                                            %Calls RR_solver for actual K-Factors
                [LCFUGL]=simplex_fug_saft(x,T,P,1,saftparam);     %Calculate Liquid Phase Fugacities
                [LCFUGV]=simplex_fug_saft(y,T,P,2,saftparam);     %Calculate Vapor Phase Fugacities
                LK(iter_total,:)=LCFUGL-LCFUGV;                                         %Calculate Ln(K)
                K=exp(LK(iter_total,:));                                                %Update K-Factors
                if iter_total>2
                    delta_LK=max(abs(LK(iter_total,:)-LK(iter_total-1,:)));             %Update Convergence Criterium
                    G_SS=(1-Beta)*sum(x.*(LCFUGL+log(x)))+Beta*sum(y.*(LCFUGV+log(y))); %Calculate Gibbs energy at each 5th step
                    Dk=LK(end-1,:)-LK(end-2,:);
                    Dk1=LK(end,:)-LK(end-1,:);
                    lambda=sum(Dk1.^2)/sum(Dk.*Dk1);                                    % Eq.33 Ch.10
                    factor(iter_total-2)=lambda/(1-lambda);
                end
                if iter_total>3
                    dfactor=abs(factor(iter_total-2)-factor(iter_total-3))*100/abs(factor(iter_total-3));
                end
            end
            
            %Acceleration trial
            if delta_LK>1e-7 && Beta>=0 && Beta<=1
                iter_acel=iter_acel+1;                                              %Count one additional trial for extrapolation
                LNK=LK(end,:)+Dk1*lambda/(1-lambda);                                %Extrapolation
                K_ext=exp(LNK);                                                     %Calculate Extrapolated K-Factors
                [Beta,x,y]=RR_solver(K_ext,z,1);                                    %Call RR_Solver to obtain new compositions
                [LCFUGL]=simplex_fug_saft(x,T,P,1,saftparam);     %Calculate Liquid Phase Fugacities
                [LCFUGV]=simplex_fug_saft(y,T,P,2,saftparam);     %Calculate Vapor Phase Fugacities
                G_A=(1-Beta)*sum(x.*(LCFUGL+log(x)))+Beta*sum(y.*(LCFUGV+log(y)));  %Calculate Gibbs energy after accelerated step
                %If Gibbs energy is decreased accept extrapolated K_values
                if G_A-G_SS<0
                    K=K_ext;
                    iter_extrap=iter_extrap+1;                                      %Count a successful extrapolation
                    LK(iter_total,:)=LNK;
                    delta_LK=max(abs(LK(iter_total,:)-LK(iter_total-1,:)));         %Update Convergence Criterium
                end
            end
            dfactor=2;
        end
    else
        Error ('Invalid Approach selected for Acceleration in the Two-Phase PT-FLASH');
    end
else
    Error ('Invalid Method choosen for Stability Analysis')
end