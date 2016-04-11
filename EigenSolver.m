function [ flux,k,iter ] = EigenSolver( a,h )
%Brandon Kuhnert
%Inputs: 

%Parameters

D=1; %cm    Diffusion Coefficient
Sab=0.7; %cm^-1     Macroscopic xsec for Absorption
S=8; %n/(cm^3*s)    Source
vSf=0.6; %cm^-1      Macroscopic xsec for fission * Nu
L=sqrt(D/Sab); %diffusion length

x=(-a+h):h:(a-h); %cell-centered
n=2*a/h-1;

%Build matrix

mdiag=((2*D/(h^2))+Sab).*ones(n,1);
sdiag=((-D)/(h^2)).*ones(n-1,1);
A=diag(mdiag)+diag(sdiag,-1)+diag(sdiag,1);

%Define tolerance

tolk=10^(-4);
tolflux=10^(-4);

%guess

flux0=ones(n,1);
k=1;

b=flux0.*vSf;
kres=10;

while kres>=tolk
    kold=k;
    [phi,iter]=GaussSeidel(A,b,tolflux);
    k=(sum(phi*vSf)/sum(flux0*vSf));
    kres=abs(k-kold);
    flux=(1/k).*phi;
end

plot(x,flux,'ro')
xlabel('x cm')
ylabel('flux')
title('Eigenvalue Criticality')

end

