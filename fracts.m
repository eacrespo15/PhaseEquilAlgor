function [guess_beta,guess_x]=fracts(NF,NC,zfeed,guess_beta,K)
%Material balances for Heidemann's flash routine
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Objective:Material balances for Heidemann's flash rotine - determination
%of phase fractions and compositions.
%
%Reference: Robert R. Heidemann works - Univ. Calgary
%
%Input variables:
%NF-    Number of phases (guess)
%NC-    Number of components
%z-     feed composition
%guess_beta- Phase fractions (guess) (1XNF)
%K- Equilibrium K-factors (NFxNC)
%
%Output variables:
%guess_beta - Phase fractions (1XNF)
%guess_x-     Phase compositions (NFxNC)


delmax=1;
maxiter=50;
tmax=1e-10;
flag=0;
alp=1e-7;
func=0;
dbetas=zeros(NF,1);
A=zeros(NF,NF);
ipiv=0;

niter=0;
while flag==0
    niter=niter+1;
    
    den=zeros(1,NC);
    for i=1:NC
        for j=1:NF
            den(i)=den(i)+guess_beta(j)*K(j,i);
        end
    end
    
    s=zeros(NF,1);
    guess_x=zeros(NF,NC);
    for j=1:NF
        for i=1:NC
            guess_x(j,i)=zfeed(i)*K(j,i)/den(i);
            s(j)=s(j)+guess_x(j,i);
        end
    end
    
    func0=func;
    func=0;
    for j=1:NF
        func=func+guess_beta(j);
    end
    
    for i=1:NC
        func=func-zfeed(i)*log(den(i));
    end
    
    s0=zeros(NF,1);
    beta0=zeros(NF,1);
    x0=zeros(NF,NC);
    if niter==1
        func0=func+1e-3;
        for j=1:NF
            s0(j)=s(j);
            beta0=guess_beta(j);
            for i=1:NC
                x0(j,i)=guess_x(j,i);
            end
        end
    end
    
    if func>func0+1e-8
        for j=1:NF
            s(j)=s0(j);
            guess_beta(j)=beta0(j);
            for i=1:NC
                guess_x(j,i)=x0(j,i);
            end
        end
        delmax=delmax/3;
    end
    
    for j=1:NF
        beta0(j)=guess_beta(j);
        s0(j)=s(j);
        for i=1:NC
            x0(j,i)=guess_x(j,i);
        end
    end
    
    test=0;
    dummy=0;
    for j=1:NF
        dummy=dummy+guess_beta(j);
    end
    
    for j=1:NF
        guess_beta(j)=guess_beta(j)/dummy;
    end
    
    for j=1:NF
        if ~(guess_beta(j)==0 && s(j)<1)
            dummy=abs(1-s(j));
            if dummy>test
                test=dummy;
            end
        end
    end
    
    if test<=tmax || niter>maxiter
        for j=1:NF
            for i=1:NC
                guess_x(j,i)=guess_x(j,i)/s(j);
            end
        end
        flag=1;
        return
    end
    
    for j=1:NF
        dbetas(j)=-(1-s(j));
        for k=j:NF
            dummy=0;
            for i=1:NC
                dummy=dummy+guess_x(j,i)*guess_x(k,i)/zfeed(i);
            end
            if j==k
                dummy=dummy+alp;
            end
            A(j,k)=dummy;
            A(k,j)=dummy;
        end
    end
    
    for j=1:NF
        if guess_beta(j)<=0
            for k=1:NF
                if j~=k
                    A(k,j)=0;
                end
            end
        end
    end
    
    [ipiv,A]=mat_lowerupper(NF,ipiv,A);
    [dbetas]=back(NF,ipiv,A,dbetas);
    
    dummy=0;
    fa=1;
    fb=0;
    
    for j=1:NF
        if beta0(j)>0
            tem=-dbetas(j)/beta0(j);
            if tem>fb
                fb=tem;
            end
        end
    end
    
    fb=0.9999*fb;
    if fb>1
        fa=delmax/fb;
    end
    
    for j=1:NF
        guess_beta(j)=beta0(j)+fa*dbetas(j);
        if guess_beta(j)<0
            guess_beta(j)=0;
        end
    end
    
    for j=1:NF
        dummy=dummy+guess_beta(j);
    end
    
    for j=1:NF
        guess_beta(j)=guess_beta(j)/dummy;
    end
end
    for j=1:NF
        for i=1:NC
            guess_x(j,i)=guess_x(j,i)/s(j);
        end
    end
end

function [IPIV,A]=mat_lowerupper(N,IPIV,A)
%% Lower-Upper factorization of matrices

IPIV(N)=N;
N1=N-1;
for i=1:N1
    x=abs(A(i,i));
    IPIV(i)=i;
    i1=i+1;
    for j=i1:N
        y=A(j,i);
        if(y<0)
            y=-y;
        end
        if y>x
            x=y;
            IPIV(i)=j;
        end
    end
    
    if IPIV(i)~=i
        k=IPIV(i);
        for j=i:N
            x=A(i,j);
            A(i,j)=A(k,j);
            A(k,j)=x;
        end
    end
    
    for j=i1:N
        x=-A(j,i)/A(i,i);
        A(j,i)=x;
        for k=i1:N
            A(j,k)=A(j,k)+x*A(i,k);
        end
    end
end
end

function [v]=back(n,ipiv,a,v)
n1=n-1;
for i=1:n1
    i1=i+1;
    k=ipiv(i);
    if k~=i
        x=v(i);
        v(i)=v(k);
        v(k)=x;
    end
    for j=i1:n
        v(j)=v(j)+a(j,i)*v(i);
    end
end

v(n)=v(n)/(a(n,n)+1e-20);
for ii=2:n
    i=n+1-ii;
    i1=i+1;
    for j=i1:n
        v(i)=v(i)-a(i,j)*v(j);
    end
    v(i)=v(i)/(a(i,i)+1e-20);
end
end



