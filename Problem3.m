function  Problem3
%Inputs
a=4;
D=1;
Sab=0.2;
S=8;
L=sqrt(D/Sab);

h=[1 0.5 0.25 0.1 0.05 0.01];
N=length(h);

for i=1:N
    %Numerical Expression
    [x,flux,n(i)]=FSDEdifferenceProblem3(4,h(i),1,0.2,(@(x) 8));
    %Analytical Expression
    fluxa=((-S/Sab)/exp(a/L)+exp(-a/L)).*(exp(x./L)+exp(-x./L))+(S/Sab);
    %relative error
    err(i)=max((fluxa-flux')./flux');
end

%plot
loglog(n,err)
xlabel('Mesh Size')
ylabel('Maximum Relative Error')
title('Error vs. Mesh Size')
end

