function [ x,flux,n,A ] = FSDEdifference( a,h,D,Sab,S )
%Brandon Kuhnert
%Inputs: a boundary conditions, h mesh size, D diffusion constant,
%Sab macroscopic absorption xsec, S source function
%Output: x corresponding x's to flux valuses, flux vector

L=sqrt(D/Sab); %diffusion length
n=2*a/h;
x= (-a+h):h:(a-h);

%Build the matrix and vector to solve

mdiag=((2*D/(h^2))+Sab).*ones(n,1);
sdiag=((-D)/(h^2)).*ones(n-1,1);
A=diag(mdiag)+diag(sdiag,-1)+diag(sdiag,1);

s=zeros(n-1,1);
for k=1:(n-1)
    s(k)=S(x(k));
end

%Solve using Thomas Algorithm
%define constants

N=length(s);
t=zeros(N,1);
y=t;
d=diag(A);
q=diag(A,1);
a=diag(A,-1);
w=d(1);
y(1)=s(1)/w;

for i=2:N
    t(i-1)=q(i-1)/w;
    w=d(i)-a(i-1)*t(i-1);
    y(i)=(s(i)-a(i-1)*y(i-1))/w;
end
for j=(N-1):(-1):1
    y(j)=y(j)-(t(j)*y(j+1));
end

flux=y;

%plot

plot(x,flux,'r');

end

