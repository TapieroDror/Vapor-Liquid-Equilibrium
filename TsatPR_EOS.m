function T_sat = TsatPR_EOS(prop_num,P)
% Parameters: T,P,Tc,Pc,w
% T: Temperature [K]                                          
% P: Pressure [Pa]                                             
% Tc: critical temperature [K]                               
% Pc: critical pressure [Pa]                                   
% w: accentic factor
R = 8.3144621; % gas constant J/(mol K)

prop=criticalProperties(prop_num);
Tc=prop(3); 
Pc=prop(4)*10^6;
w=prop(6);
T=Tc*0.0001;

% Reduced variables
Tr = T/Tc ;

% Parameters of the EOS for a pure component
f = 0.37464 + 1.54226*w - 0.26992*w^2;% simplification term

% w is the acentric factor
alfa = (1 + f*(1 - sqrt(Tr)))^2;
a = 0.45724*(R*Tc)^2/Pc*alfa;% attraction parameter
b = 0.0778*R*Tc/Pc;% repulsion parameter
A = a*P/(R*T)^2;
B = b*P/(R*T);

% Compressibility factor
Z = roots([1 -(1-B) (A-3*B^2-2*B) -(A*B-B^2-B^3)]);
ZR = [];
while (isreal(Z(1))==0||isreal(Z(2))==0||isreal(Z(3))==0)
    T=T+0.001;
    Tr = T/Tc ;
    alfa = (1 + f*(1 - sqrt(Tr)))^2;
    a = 0.45724*(R*Tc)^2/Pc*alfa;
    b = 0.0778*R*Tc/Pc;
    A = a*P/(R*T)^2;
    B = b*P/(R*T);
    Z = roots([1 -(1-B) (A-3*B^2-2*B) -(A*B-B^2-B^3)]);
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
fhi1 = exp(Z1 - 1 - log(Z1-B) - A/(2*B*sqrt(2))*log((Z1+(1+sqrt(2))*B)/(Z1+(1-sqrt(2))*B)));
fhi2 = exp(Z2 - 1 - log(Z2-B) - A/(2*B*sqrt(2))*log((Z2+(1+sqrt(2))*B)/(Z2+(1-sqrt(2))*B)));

% Finally calculating the saturated pressure
while abs((fhi1-fhi2))>0.0001
    if(fhi1>fhi2)
        T=T-0.001;
    else
        T=T+0.001;
    end
    Tr = T/Tc ;
    alfa = (1 + f*(1 - sqrt(Tr)))^2;
    a = 0.45724*(R*Tc)^2/Pc*alfa;
    b = 0.0778*R*Tc/Pc;
    A = a*P/(R*T)^2;
    B = b*P/(R*T);
    Z = roots([1 -(1-B) (A-3*B^2-2*B) -(A*B-B^2-B^3)]);
    ZR = [];
    for i = 1:3
       if isreal(Z(i))
        ZR(i) = Z(i);  
       end
    end
    Z1 = min(ZR);   
    Z2 = max(ZR);
    fhi1 = exp(Z1 - 1 - log(Z1-B) - A/(2*B*sqrt(2))*log((Z1+(1+sqrt(2))*B)/(Z1+(1-sqrt(2))*B)));
    fhi2 = exp(Z2 - 1 - log(Z2-B) - A/(2*B*sqrt(2))*log((Z2+(1+sqrt(2))*B)/(Z2+(1-sqrt(2))*B)));
end
T_sat=T;


