For 0 source

[x,flux,n]=FSDEdifference(4,0.1,1,0.2,(@(x) 0));
hold on
fluxa=zeros(1,n-1);
plot(x,fluxa,'bo');
legend('Numerical','Analytical')
xlabel('x cm')
ylabel('Flux')
title('No Source')
hold off

For Constant Source 

[x,flux,n]=FSDEdifference(4,0.1,1,0.2,(@(x) 8));
hold on
L=sqrt(1/0.2);
fluxa=((-8/0.2)/(exp(-4/L)+exp(4/L))).*(exp(-x./L)+exp(x./L))+(8/0.2);
plot(x,fluxa,'bo')
legend('Numerical','Analytical')
xlabel('x cm')
ylabel('Flux')
title('Constant Source')
hold off

For Cosine 

[x,flux,n]=FSDEdifference(4,0.1,1,0.2,(@(x) cos(x)));
hold on
L=sqrt(1/0.2);
D=1;
fluxa=(-1/(exp(-4/L)+exp(4/L)))*(cos(4)/(D*(1+1/(L^2)))).*(exp(-x./L)+exp(x./L))+cos(x)./(D*(1+1/L^2));
plot(x,fluxa,'bo')
legend('Numerical','Analytical')
xlabel('x cm')
ylabel('Flux')
title('Cosine')