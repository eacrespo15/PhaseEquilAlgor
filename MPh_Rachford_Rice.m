function [Betas,y]=MPh_Rachford_Rice(NF,z,FUGs,Betas_guess)
%Rachford-Rice multiphase equations solver
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Objective:Solver for the Rachford-Rice Equations to be used in 
%Multiphase Flash Calculations
%
%Reference: Notes taken in DTU PhD. Summer course on Thermodynamic models;
%Michelsen & Mollerup - Thermodynamic models - Fundamentals and
%Computational Aspects

flag=0;
count_total=1;
count1=1;
count_active=0;


%Step 1 - Define all NF Phases with Beta >0 as active
Betas=Betas_guess;
Logical_Betas=zeros(NF,1);
Logical_Betas(Betas>0)=1;

%Step 2 - Evaluate Q, calculate g and H
E=Betas'*(1./FUGs);
Q=sum(Betas)-z*log(E)';

%Gradient = dQ/dBetas
g=zeros(1,NF);
for k=1:NF
    g(k)=1-sum(z./(E.*FUGs(k,:)));
end

%Hessian Matrix= dg/dBetas
H=zeros(NF,NF);
for k=1:NF
    for l=1:NF
        H(k,l)=sum(z./(E.^2.*FUGs(l,:).*FUGs(k,:)));
        if k==l
            H(k,l)=H(k,l)+1e-10;
        end
    end
end

% Save a copy of Q function , gradient g and betas
g_copy(count1,:)=g; 
q_copy(count1,:)=Q;
betas_copy(count1,:)=Betas;


while flag==0 && count_active<5
    ERROR=1;
    count1=1;
    g_copy(2:end,:)=[];
    q_copy(2:end,:)=[];
    
    while ERROR>1e-5 && count1<100
        count_total=count_total+1;
        count1=count1+1;
        
        %If a phase j is inactive set g=0 and Hjk=Hkj=0 for any k phase and Hjj=1;
        g(Logical_Betas==0)=0;
        H(Logical_Betas==0,:)=0;
        H(:,Logical_Betas==0)=0;
        H(Logical_Betas==0,Logical_Betas==0)=1;
        
        %Step 3 - Calculate the Newton Step
        w=-H\g';
        ERROR=norm(w);
        
        % Step 4 - Calculate the new Betas!
        alpha=1;
        epsilon=1e-10;
        Betas_new=Betas+alpha*w;
        
        % Step 5- Check for negative Betas - Proceed to deactivate!
        if any(Betas_new<0)==1
            %Determine alpha so that all Betas become non-negative
            index_min_beta=find(Betas_new==min(Betas_new));
            alpha=-Betas_new(index_min_beta)/w(index_min_beta);
            Betas_new=Betas_new+alpha*w;
            %Evaluate the Q function for the new Betas
            E=Betas_new'*(1./FUGs);
            Qnew=sum(Betas_new)-z*log(E)';
            if Qnew<Q+epsilon
                Logical_Betas(index_min_beta)=0;
            else
                alpha=alpha/2;
                Betas_new=Betas+alpha*w;
            end
        end
        
        %Step 5/6
        E=Betas_new'*(1./FUGs);
        Qnew=sum(Betas_new)-z*log(E)';
        
        while Qnew>Q+epsilon
            alpha=alpha/2;
            Betas_new=Betas+alpha*w;
            Betas=Betas_new;
            E=Betas_new'*(1./FUGs);
            Qnew=sum(Betas_new)-z*log(E)';
        end
        
        Q=Qnew;
        Betas=Betas_new;
        
        % Atualizar E, Q, g e H
        
        E=Betas'*(1./FUGs);
        Q=sum(Betas)-z*log(E)';
        
        g=zeros(1,NF);
        for k=1:NF
            g(k)=1-sum(z./(E.*FUGs(k,:)));
        end
        
        H=zeros(NF,NF);
        for k=1:NF
            for l=1:NF
                H(k,l)=sum(z./(E.^2.*FUGs(l,:).*FUGs(k,:)));
                if k==l
                    H(k,l)=H(k,l)+1e-10;
                end
            end
        end
        
        g_copy(count1,:)=g;
        q_copy(count1,:)=Q;
        betas_copy(count1,:)=Betas;
        
    end
    
    %Step 7 - Check if it is necessary to activate/reactivate elements
    flag=1;
    if any(Logical_Betas==0)==1
       uzeros=find(Logical_Betas==0);
       if any(g_copy(end,uzeros)<=0)==1
           flag=0;
           count_active=count_active+1;
           Logical_Betas(uzeros(g_copy(end,uzeros)==min(g_copy(end,uzeros))))=1;
       end
    end 
    
end

%Step 8 - Calculate mole fractions in all phases
reped_z=repmat(z,NF,1);
reped_E=repmat(E,NF,1);
y=reped_z./reped_E./FUGs;
