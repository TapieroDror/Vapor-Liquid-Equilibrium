% Liquid and Vapor Mole Fraction calculator using Peng-Robindon Equation of
% State and van der Waals Mixing Rules
% Do not Copy - All Rights Reserved to Tapiero Dror
%--------------------------------------------------------------------------
%critical constants and acentric factor
%Mw - molar mass [kg/kmol]
%R - gas constant [kJ/kg*K]
%Pc - critical pressure [MPa]
%Tc - critical temperature [K]
%Vc - critical volume [m^3/kmol]
%w - acentric factor [-] 
%Hf - Enthalpy of Formation [J/mol] (at 25[C], 100[KPa]) page 239
%--------------------------------------------------------------------------
%Substances Table:
%[number    substance                          formula  ]
%[1         Air                                -        ]
%[2         Ammonia                            NH3      ]
%[3         Argon                              Ar       ]
%[4         Benzene                            C6H6     ]
%[5         Bromine                            Br2      ]
%[6         n-Butane                           C4H10    ]
%[7         Carbon dioxide                     CO2      ]
%[8         Carbon monoxide                    CO       ]
%[9         Carbon tetrachloride               CCl4     ]
%[10        Chlorine                           Cl2      ]
%[11        Chloroform                         CHCl3    ]
%[12        Dichlorodifluoromethane(R-12)      CCl2F2   ]
%[13        Dichlorofluoromethane(R-21)        CHCl2F   ]
%[14        Ethane                             C2H6     ]
%[15        Ethyl alcohol (Ethanol)*           C2H5OH   ]
%[16        Ethylene (Ethene)                  C2H4     ]
%[17        Helium                             He       ]
%[18        n-Hexane                           C6H14    ]
%[19        Hydrogen(normal)                   H2       ]
%[20        Krypton                            Kr       ]
%[21        Methane                            CH4      ]
%[22        Methyl alcohol (Methanol)*         CH3OH    ]
%[23        Methyl chloride                    CH3Cl    ]
%[24        Neon                               Ne       ]
%[25        Nitrogen                           N2       ]
%[26        Nitrous oxide                      N2O      ]
%[27        Oxygen                             O2       ]
%[28        n-Pentane                          C5H12    ]
%[29        Propane                            C3H8     ]
%[30        Propylene (Propene)                C3H6     ]
%[31        Sulfur dioxide                     SO2      ]
%[32        Tetrafluoroethane(R-134a)          CF3CH2F  ]
%[33        Toluene                            C7H8     ]
%[34        Trichlorofluoromethane(R-11)       CCl3F    ]
%[35        Water*                             H2O      ]
%[36        Xenon                              Xe       ]
%*[Hf for gas                      ]
% [for liq Ethanol  (Hf(gas)-42380)]
% [for liq Methanol (Hf(gas)-37920)]
% [for liq Water    (Hf(gas)-44004)]
%--------------------------------------------------------------------------
%indexes for alpha function:
%(1) Soave (ORIGINAL)
%(2) Milind Joshipura
%(3) Khaled A. M. Gasem
%(4) Chorng H. Twu
%(5) Hamid Saffari and Alireza Zahedi (very good for Natural Gas)
%(6) Georges A. Melhem
%(7) Peyman Mahmoodi and Maryam Sedigh (proposed model 1)
%(8) Peyman Mahmoodi and Maryam Sedigh (proposed model 2)
%--------------------------------------------------------------------------
clc;clear;close all;
%---------------------------------inputs:----------------------------------
T=260;%[K]
P=1.5;%[MPa]
propNum=[6 7 14];
index=1;
nmoles=[20 50 30];
%--------------------------------------------------------------------------
%Display the names of the Substances
disp(["index = " num2str(index)])
options = optimset('Display','off');
fprintf('    ')
for i=1:length(nmoles)
    fprintf(PropertieNames(propNum(i))+'   ')
end
fprintf('\n')

%Calculate mole fractions
z=zeros(1,length(nmoles));
nsum=sum(nmoles);
for i=1:length(nmoles)
    z(i)=nmoles(i)/nsum;
end
disp('z: ')
disp(z)

%Guess for the molar fractions of the liquid for each component (will be used to solve the Rachford-Rice equations)
x=z;%guess for x

%Peng-Robinson Equation of State
comp=length(x);
for i=1:comp
    prop=criticalProperties(propNum(i));
    p(:,i)=prop';
end
Pc=p(4,:);%[MPa]
Tc=p(3,:);%[K]
w=p(6,:);
Tr=T./Tc;
Pr=P./Pc;
alpha=alphaFunctions(index,w,Tr,propNum);
% k=0.37464 + 1.54226.*w-0.26992.*w.^2;
% alpha=(1+k.*(1-Tr.^0.5)).^2;
Ai=0.45724.*alpha.*Pr./Tr.^2;
Aj=Ai;
Bi=0.07780.*Pr./Tr;

%Van der Waals Mixing Rules
kij=Kij(propNum);
%kij=zeros(length(propNum));
Aij=zeros(comp);
for i=1:comp
    for j=1:comp
    Aij(i,j)=(Ai(i)*Aj(j))^0.5;
    end
end
Aij=(1-kij).*Aij;
a=zeros(comp);
for i=1:comp
    for j=1:comp
        a(i,j)=x(i)*x(j)*Aij(i,j);
    end
end
Al=sum(sum(a));
Av=Al;
b=zeros(1,comp);
for i=1:comp
    b(i)=x(i)*Bi(i);
end
Alsum=rand(1,comp)*0;
for i=1:comp
    for j=1:comp
        Alsum(i)=Alsum(i)+x(j)*Aij(i,j);
    end
end
Avsum=rand(1,comp)*0;
for i=1:comp
    for j=1:comp
        Avsum(i)=Avsum(i)+x(j)*Aij(i,j);
    end
end
Bl=sum(b);
Bv=Bl;

% liquid and gas phase compressibility factors
Zv=max(roots([1 -1+Bv Av-3*Bv^2-2*Bv -Av*Bv+Bv^2+Bv^3]));
Zl=min(roots([1 -1+Bl Al-3*Bl^2-2*Bl -Al*Bl+Bl^2+Bl^3]));

% vapor and liquid phase fugacity coefficients 
phiv=exp((Zv-1).*Bi/Bv-log(Zv-Bv)-Av/(2*sqrt(2)*Bv)*log((Zv+(1+sqrt(2))*Bv)/(Zv+(1-sqrt(2))*Bv)).*(2.*Avsum./Av-Bi./Bv));
phil=exp((Zl-1).*Bi/Bl-log(Zl-Bl)-Al/(2*sqrt(2)*Bl)*log((Zl+(1+sqrt(2))*Bl)/(Zl+(1-sqrt(2))*Bl)).*(2.*Alsum./Al-Bi./Bl));

% equilibrium constant
K=phil./phiv;

% Solve Rachford-Rice Equations
x_y=[x 1-x 0.5];%the last numper is the guess of psi
[X]=fsolve(@flashFinal,x_y,options,x,K);
x=zeros(1,comp);
y=zeros(1,comp);
for i=1:(length(X)-1)/2
    x(i)=X(i);
end
for i=((length(X)-1)/2+1):(length(X)-1)
    y(i-((length(X)-1)/2+1)+1)=X(i);
end 

%Find the moles of each components in vapor and liquid
N=ones(1,comp*2);
options=optimoptions('fsolve','Display','off','MaxFunEvals',100,'Algorithm','levenberg-marquardt');
[n]=fsolve(@moleFinder,N,options,x,y,z,nmoles);
nMat=zeros(comp,2);
for i=1:comp
    nMat(i,1)=n(i);
    nMat(i,2)=n(i+comp);
end
check=sum(find(nMat<0));
if sum(check)~=0||isreal(x)==0||isreal(y)==0
    error('saturated mixture can not exist in this ratio')
end

%Display the moles fraction for each components in vapor and liquid
disp('Saturated Mixture: ')
disp('x: ')
disp(x)
disp('y: ')
disp(y)

%Display the moles for each components in vapor and liquid
names=string(zeros(1,length(propNum)));
for i=1:length(propNum)
    names(i)=PropertieNames(propNum(i));
end
names=names';
disp('nMat: ')
phase=string(zeros(1,3));
phase(1)="       ";
phase(2)="n(liq)";
phase(3)="n(vap)";
disp([phase;names num2str(nMat(:,1)) num2str(nMat(:,2))])

