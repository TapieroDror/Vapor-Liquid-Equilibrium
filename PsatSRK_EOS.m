function P_sat = PsatSRK_EOS(prop_num,T)
% Parameters: T,P,Tc,Pc,w
% T: Temperature [K]                                          
% P: Pressure [Pa]                                             
% Tc: critical temperature [K]                               
% Pc: critical pressure [Pa]                                   
% w: accentic factor

prop=criticalProperties(prop_num);
Tc=prop(3); 
Pc=prop(4)*10^6;
w=prop(6);
P=Pc*0.1;

% Reduced variables
Tr = T/Tc ;

% Parameters of the EOS for a pure component
f = 0.480 + 1.574*w - 0.176*w^2;
alfa = (1 + f*(1 - sqrt(Tr)))^2;
A=0.42747*alfa*(P/Pc)/(T/Tc)^2;
B = 0.08664*(P/Pc)/(T/Tc);

% Compressibility factor
Z = roots([1 -1 (A-B-B^2) -(A*B)]);
ZR = [];
while (isreal(Z(1))==0||isreal(Z(2))==0||isreal(Z(3))==0)
    P=P+10;
    A=0.42747*alfa*(P/Pc)/(T/Tc)^2;
    B = 0.08664*(P/Pc)/(T/Tc);
    Z = roots([1 -1 (A-B-B^2) -(A*B)]);
    ZR = [];
end
for i = 1:3
   if isreal(Z(i))
      ZR(i) = Z(i);     
   end
end
Z1 = min(ZR);   
Z2 = max(ZR);
    
% Fugacity coefficient
fhi1 = exp(Z1 - 1 - log(Z1-B) - (A/B)*log((Z1+B)/Z1));
fhi2 = exp(Z2 - 1 - log(Z2-B) - (A/B)*log((Z2+B)/Z2));

% Finally calculating the saturated pressure
while abs((fhi1-fhi2))>0.00001
    if(fhi1<fhi2)
        P=P-1;
    else
        P=P+1;
    end
    A=0.42747*alfa*(P/Pc)/(T/Tc)^2;
    B = 0.08664*(P/Pc)/(T/Tc);
    Z = roots([1 -1 (A-B-B^2) -(A*B)]);
    ZR = [];
    for i = 1:3
       if isreal(Z(i))
          ZR(i) = Z(i);    
       end
    end
    Z1 = min(ZR);   
    Z2 = max(ZR);
    fhi1 = exp(Z1 - 1 - log(Z1-B) - (A/B)*log((Z1+B)/Z1));
    fhi2 = exp(Z2 - 1 - log(Z2-B) - (A/B)*log((Z2+B)/Z2));
end
P_sat=P;
