function f=flash_T_CONST(x,prop1,prop2,T,index)
global z
phi=0;
% critical temperature and pressure and acentric factor
propNum=[prop1 prop2];
p1=criticalProperties(prop1);
p2=criticalProperties(prop2);
Pc=[p1(4) p2(4)]*10;%[bar]
Tc=[p1(3) p2(3)];%[K]
w=[p1(6) p2(6)];
% reduced temperature and pressure
Tr=T./Tc;
Pr=x(5)./Pc;
% alpha, Ai, Bi, Aij, A, B for the PR EOS
alpha=alphaFunctions(index,w,Tr,propNum);
Ai=0.45724.*alpha.*Pr./Tr.^2;
Bi=0.07780.*Pr./Tr;
for i=1:2
    for j=1:2
        Aij(i,j)=(Ai(i)*Ai(j))^0.5;
    end
end
Av=0;
for i=1:2
    for j=1:2
        Av=Av+x(i+2)*x(j+2)*Aij(i,j);
    end
end
Bv=0;
for i=1:2
    Bv=Bv+x(i+2)*Bi(i);
end
Bl=0;
for i=1:2
    Bl=Bl+x(i)*Bi(i);
end
Al=0;
for i=1:2
    for j=1:2
        Al=Al+x(i)*x(j)*Aij(i,j);
    end
end
Alsum=[0 0];
for i=1:2
    for j=1:2
        Alsum(i)=Alsum(i)+x(j)*Aij(i,j);
    end
end
Avsum=[0 0];
for i=1:2
    for j=1:2
        Avsum(i)=Avsum(i)+x(j+2)*Aij(i,j);
    end
end
% liquid and gas phase compressibility factors
Zv=max(roots([1 -1+Bv Av-3*Bv^2-2*Bv -Av*Bv+Bv^2+Bv^3]));
Zl=min(roots([1 -1+Bl Al-3*Bl^2-2*Bl -Al*Bl+Bl^2+Bl^3]));
% vapor and liquid phase fugacity coefficients 
phiv=exp((Zv-1).*Bi/Bv-log(Zv-Bv)-Av/(2*sqrt(2)*Bv)*log((Zv+(1+sqrt(2))*Bv)/(Zv+(1-sqrt(2))*Bv)).*(2.*Avsum./Av-Bi./Bv));
phil=exp((Zl-1).*Bi/Bl-log(Zl-Bl)-Al/(2*sqrt(2)*Bl)*log((Zl+(1+sqrt(2))*Bl)/(Zl+(1-sqrt(2))*Bl)).*(2.*Alsum./Al-Bi./Bl));
% equilibrium constant
K=phil./phiv;
% the system of five algebraic equations
for i=1:2
    f(i)=x(i+2)-K(i)*x(i);
end
for i=1:2
    f(i+2)=x(i)-z(i)/(1+phi*(K(i)-1));
end
f(5)=0;
for i=1:2
    f(5)=f(5)+z(i)*(K(i)-1)/(1+phi*(K(i)-1));
end