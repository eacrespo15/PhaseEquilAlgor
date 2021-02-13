function [NF,guess_beta,K]=chkphs(NF,NC,guess_beta,K,V)
%Stability analysis for Heidemann's Flash Routine
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Reference: Robert R. Heidemann works - Univ. Calgary
%
%Input variables:
%NF- Test number of phases
%NC- Number of components
%guess_beta- A guess of the phase fractions
%K- Equilibrium K-factors
%V- Molar volumes=(1/rho) in m3/mol
%
%Output variables:
%NF - Final number of phases
%guess_beta - Phase fractions
%K- Equilibrium K-factors


if NF~=1
    kph=1;
    while kph<NF
        jph=kph;
        while jph <NF
            jph=jph+1;
            count1=0;
            count2=0;
            rmax=0;
            for i=1:NC
                tmp=abs(1-K(jph,i)/K(kph,i));
                if tmp<5e-2
                    count1=count1+1;
                end
                if tmp>rmax
                    rmax=tmp;
                end
            end
            
            tmp=abs(1-V(jph)/V(kph));
            if tmp<5e-2
                count2=count2+1;
            end
            if tmp>rmax
                rmax=tmp;
            end
            
            if rmax<5e-2 || (count1>=NC-1 && count2==1)
                guess_beta(kph)=guess_beta(kph)+guess_beta(jph);
                for j=jph:NF-1
                    guess_beta(j)=guess_beta(j+1);
                    for i=1:NC
                        K(j,i)=K(j+1,i);
                    end
                end
                NF=NF-1;
            end
        end
           kph=kph+1;
    end
end
    
end

