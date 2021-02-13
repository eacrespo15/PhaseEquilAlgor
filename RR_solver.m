function [Beta,x,y]=RR_solver(K,z,solvertype)
%Rachford-Rice Equation Solver
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code revised in: November 2019
%
%Objective:Solver for the Rachford-Rice Equations to be used in Flash Calculations
%
%References: Notes taken in DTU PhD. Summer course on Thermodynamic models;
%Michelsen & Mollerup - Thermodynamic models - Fundamentals and
%Computational Aspects - p.253-254; Chapter 10
%
%Solver types:
%Full RRSolver(0); Negative RR_Solver (1); Sloppy RRSolver(2)

if ~exist('solvertype','var') || isempty(solvertype)
    solvertype=0;
end


if solvertype==0
    % Calculate g(0) and g(1);
    g_0=Rachford_Rice(0,z,K); %g for Beta=0
    g_1=Rachford_Rice(1,z,K); %g for Beta=1
    
    % Check if it will flash - If it does not than we are in the presence
    % of a subcooled liquid (Beta<0) or a superheated vapour (Beta>1)
    if g_0<0 ||  g_1>0
        error ('Not flashing - Single Phase Region') % Eqs. 4 and 5 of Ch.10
    end
    
    % Define upper and lower bounds for Beta
    Beta_min=0;
    Beta_max=1;
    
    % Add a better lower bound (if possible)
    % Great when there is a very volatile compound in small amount in the feed
    if any(K>1)==1
        u_LB=find(K>1);
        Rig_Beta_min=zeros(1,length(u_LB));
        for i=1:length(u_LB)
            Rig_Beta_min(i)=(K(u_LB(i))*z(u_LB(i))-1)/(K(u_LB(i))-1); %Eq.7 Ch10
        end
        if any(Rig_Beta_min>0)==1
            Beta_min=max(Rig_Beta_min);
        end
    end
    
    % Add a better upper bound (if possible)
    % Great when there is a non-volatile compound in small amount in the feed
    if any(K<1)==1
        u_UB=find(K<1);
        Rig_Beta_max=zeros(1,length(u_UB));
        for i=1:length(u_UB)
            Rig_Beta_max(i)=(1-z(u_UB(i)))/(1-K(u_UB(i))); %Eq.8 Ch10
        end
        if any(Rig_Beta_max>0)==1
            Beta_max=min(Rig_Beta_max);
        end
    end
    
    %Calculate first value of g for Beta=0.5*(Beta_min+Beta_max)
    Beta_old=0.5*(Beta_min+Beta_max);
    Beta=Beta_old;
    
    g_ini=Rachford_Rice(Beta_old,z,K);
    g=g_ini;
    
    %Initialization of dBeta and counter of iterations
    deltaBeta=1;
    counter=0;
    while abs(g)>1e-6 || abs(deltaBeta)>1e-6
        
        %Verify if for the current beta, g is >0 or <0 and redefine Beta bounds
        if g>0
            Beta_min=Beta;
        else
            %Change to liquid fraction as the iterative variable
            Beta_max=Beta;
        end
        
        % Take a Newton step
        [g_b,dg_b]=Rachford_Rice(Beta,z,K);
        deltaBeta=-g_b/dg_b;
        Beta=Beta_old+deltaBeta;
        
        % Check if Beta is between the interval prior to acceptal
        if Beta>=Beta_min && Beta<=Beta_max
            Beta=Beta;
        else
            deltaBeta= 0.5*(Beta_min+Beta_max)-Beta_old;
            Beta=Beta_old+deltaBeta;
        end
        
        % Calculate g and dg for the new Beta and update Beta_old
        [g,~]=Rachford_Rice(Beta,z,K);
        Beta_old=Beta;
        counter=counter+1;
    end
    % Calculate compositions after convergence (Eq.6 - Ch.10)
    x=z./(1-Beta+Beta*K);
    y=K.*x;
elseif solvertype==1
    %Define Lower and Upper bound for Beta
    Beta_min=-1/(max(K)-1);
    Beta_max=1/(1-min(K));
    
    % Add a better lower bound (if possible)
    % Great when there is a very volatile compound in small amount in the feed
    if any(K>1)==1
        u_LB=find(K>1);
        Rig_Beta_min=zeros(1,length(u_LB));
        for i=1:length(u_LB)
            Rig_Beta_min(i)=(K(u_LB(i))*z(u_LB(i))-1)/(K(u_LB(i))-1);
        end
        if any(Rig_Beta_min>0)==1
            Beta_min=max(Rig_Beta_min);
        end
    end
    
    % Add a better upper bound (if possible)
    % Great when there is a non-volatile compound in small amount in the feed
    if any(K<1)==1
        u_UB=find(K<1);
        Rig_Beta_max=zeros(1,length(u_UB));
        for i=1:length(u_UB)
            Rig_Beta_max(i)=(1-z(u_UB(i)))/(1-K(u_UB(i)));
        end
        if any(Rig_Beta_max>0)==1
            Beta_max=min(Rig_Beta_max);
        end
    end
    
    %Calculate first value of g for Beta=0.5*(Beta_min+Beta_max)
    Beta_old=0.5*(Beta_min+Beta_max);
    Beta=Beta_old;
    g_ini=Rachford_Rice(Beta_old,z,K);
    g=g_ini;
    
    %Initialization
    deltaBeta=1;
    counter=0;
    
    while abs(g)>1e-6 || abs(deltaBeta)>1e-6
        
        %Verify if for the current beta, g is >0 or <0 and redefine Beta bounds
        if g>0
            Beta_min=Beta;
        else
            Beta_max=Beta;
        end
        
        % Take a Newton step
        [g_b,dg_b]=Rachford_Rice(Beta,z,K);
        deltaBeta=-g_b/dg_b;
        Beta=Beta_old+deltaBeta;
        
        % Check if Beta is between the interval prior to acceptal
        if Beta>=Beta_min && Beta<=Beta_max
            Beta=Beta;
        else
            deltaBeta= 0.5*(Beta_min+Beta_max)-Beta_old;
            Beta=Beta_old+deltaBeta;
        end
        
        % Calculate g and dg for the new Beta and update Beta_old
        [g,~]=Rachford_Rice(Beta,z,K);
        Beta_old=Beta;
        counter=counter+1;
    end
    % Calculate compositions after convergence (Eq.6 - Ch.10)
    %if Beta>=0 && Beta<=1
        x=z./(1-Beta+Beta*K);
        y=K.*x;
    %else
        %x=z;
        %y=z;
    %end
elseif solvertype==2
    
    % Calculate g(0) and g(1);
    g_0=Rachford_Rice(0,z,K); %g for Beta=0
    g_1=Rachford_Rice(1,z,K); %g for Beta=1
    
    % Check if it will flash - If it does not than we are in the presence
    % of a subcooled liquid (Beta<0) or a superheated vapour (Beta>1)
    if g_0<0 ||  g_1>0
        error ('Not flashing') % Eqs. 4 and 5 of Ch.10
    end
    
    % Define upper and lower bounds for Beta
    Beta_min=0;
    Beta_max=1;
    
    % Add a better lower bound (if possible)
    % Great when there is a very volatile compound in small amount in the feed
    if any(K>1)==1
        u_LB=find(K>1);
        Rig_Beta_min=zeros(1,length(u_LB));
        for i=1:length(u_LB)
            Rig_Beta_min(i)=(K(u_LB(i))*z(u_LB(i))-1)/(K(u_LB(i))-1);%Eq.7 Ch10
        end
        if any(Rig_Beta_min>0)==1
            Beta_min=max(Rig_Beta_min);
        end
    end
    
    % Add a better upper bound (if possible)
    % Great when there is a non-volatile compound in small amount in the feed
    
    if any(K<1)==1
        u_UB=find(K<1);
        Rig_Beta_max=zeros(1,length(u_UB));
        for i=1:length(u_UB)
            Rig_Beta_max(i)=(1-z(u_UB(i)))/(1-K(u_UB(i)));%Eq.8 Ch10
        end
        if any(Rig_Beta_max>0)==1
            Beta_max=min(Rig_Beta_max);
        end
    end
    
    % Calculate first value of g for Beta=0.5*(Beta_min+Beta_max)
    Beta_old=0.5*(Beta_min+Beta_max);
    Beta=Beta_old;
    
    % Take a Newton step
    [g_b,dg_b]=Rachford_Rice(Beta,z,K);
    deltaBeta=-g_b/dg_b;
    Beta=Beta+deltaBeta;
    
    % Calculate compositions(Eq.6 - Ch.10)
    x=z./(1-Beta+Beta*K);
    y=K.*x;
else
    error('Invalid RR_Solver selected')
end
end


function [gsub,dgsub]=Rachford_Rice(Betasub,zsub,Ksub)
%Rachford-rice equation
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code revised in: November 2019
%
%Objective:Subroutine that given a vapour fraction, composition and K-factors 
%return the values of g and dg. The term sub indicates that are the
%variables from the subroutine.
%
%Reference: Notes taken in DTU PhD. Summer course on Thermodynamic models;
%Michelsen & Mollerup - Thermodynamic models - Fundamentals and
%Computational Aspects

gsub=0;
dgsub=0; %May be an optional argument
for i=1:length(zsub)
    gsub=gsub+zsub(i)*(Ksub(i)-1)/(1-Betasub+Betasub*Ksub(i)); % Eq 2 of Ch.10
    dgsub=dgsub-zsub(i)*(Ksub(i)-1)^2/(1-Betasub+Betasub*Ksub(i))^2; %Eq 3 of Ch.10
end
end

