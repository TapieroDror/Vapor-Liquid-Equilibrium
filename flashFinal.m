%solving Rashford-Rice equations
function f=flashFinal(x,z,K)
comp=length(z);
% the system of five algebraic equations
for i=1:comp
    f(i)=x(i+comp)-K(i)*x(i);
end
for i=1:comp
    f(i+comp)=x(i)-z(i)/(1+x(comp*2+1)*(K(i)-1));
end
f(comp*2+1)=0;
for i=1:comp
    f(comp*2+1)=f(comp*2+1)+z(i)*(K(i)-1)/(1+x(comp*2+1)*(K(i)-1));
end