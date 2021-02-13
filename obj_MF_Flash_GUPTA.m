function [fval]=obj_MF_Flash_GUPTA(XBeta,z,K,NF,NC)
%Objective function for Gupta multiphase flash algorithm
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Objective: objective function for MultiFlash algorithm, modification of Rachford-Rice equation 
%for multiple phases accounting for phase stability
%
%References: Martin et al. 2011 - 10.1016/j.ece.2011.08.003
%Gupta et al., Fluid Phase Equilibria, 63 (1991), 65-89


%Rounding value of beta and thita
epsilon = 1e-6; 

%Initialize Betas and Thitas
beta = zeros(1,NF);
thita = zeros(1,NF);

%Construct Beta and thita vectors from X
%If a phase is not present then its fraction is 0 and it must be unstable
%so it should get a positive thita. On the contrary if a certain phase is
%present is Beta=Betaguess="X" and thita must be 0 = stable
for i = 2:NF
    if XBeta(i-1) < 0
        beta(i) = 0;
        thita(i) = -XBeta(i-1);
    else
        beta(i) = XBeta(i-1);
        thita(i) = 0;
    end
end

%The first phase is taken as reference and always must have Beta>0 and Thita=0
thita(1)=0;
beta(1)=1-sum(beta);
beta=abs(beta);

%Round beta and thita values smaller than epsilon to 0
for i=1:NF
    if beta(i)<=epsilon
        beta(i)=0;
    end
    if thita(i)<=epsilon
        thita(i)=0;
    end
end

%Normalization of Beta values
beta=beta/sum(beta);

%Calculation of objective function
fval=zeros(1,NF-1);

for j=2:NF
    sumat=0;
    for i=1:NC
        sum_int=0;
        for k=2:NF
            sum_int=sum_int+beta(k)*(K(k,i)*exp(thita(k))-1); % Di-1
        end
        sumat=sumat+z(i)*(K(j,i)*exp(thita(j))-1)/(1+sum_int); %Equation for the mole fraction summation
    end
    fval(j-1)=sumat*1e8;
end


