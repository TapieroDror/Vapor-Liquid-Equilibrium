% Liquid and Vapor Mole Fraction calculator using Peng-Robindon Equation of
% State and Van der Waals Mixing Rules
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
%(2) Khaled A. M. Gasem
%(3) Chorng H. Twu
%(4) Hamid Saffari and Alireza Zahedi (very good for Natural Gas)
%(5) Georges A. Melhem
%(6) Peyman Mahmoodi and Maryam Sedigh (proposed model 1)
%(7) Peyman Mahmoodi and Maryam Sedigh (proposed model 2)
%--------------------------------------------------------------------------
clc;clear;close all;
%------------------------------Article Data:-------------------------------
%Conditions: T[K]
%260[K]
% A=[3.15	0.0000	0.0000
% 3.62	0.1700	0.0600
% 5.02	0.4520	0.1680
% 6.34	0.6120	0.2830
% 9.37	0.7840	0.5010
% 11.56	0.8990	0.6720
% 14.10	0.9680	0.8550
% 15.37	0.9800	0.9290
% 16.99	1.0000	1.0000
% ];
%A(:,1)=A(:,1)*1.01325; atm->bar
%--------------------------------------------------------------------------
%---------------------------------inputs:----------------------------------
P=2.5; %[MPa]
Pbar=P*10; %[bar]
propertie1=21; %number from table (propertie 1)
propertie2=14; %number from table (propertie 2)
index=3;
%--------------------------------------------------------------------------
% flash calculation using fsolve and a zeroorder collocation method
global z % Define global variable.
clear sol % sol is a variable of type double
% flash calculation using fsolve and a zeroorder collocation method

Trange=(TsatPR_EOS(propertie1,Pbar*100*10^3)+TsatPR_EOS(propertie2,Pbar*100*10^3))/2;
z=[0.0001 0.9999];
options = optimset('Display','off');
args=[0.01 0.9 0.01 0.9 Trange];
[X]=fsolve(@flash_P_CONST,args,options,propertie1,propertie2,Pbar,index);
x0=X;
sol(1,1)=X(1);
sol(2,1)=X(3);
sol(3,1)=X(5);
for i=1:100
    z=[0.01*i 1-0.01*i];
    [X]=fsolve(@flash_P_CONST,x0,options,propertie1,propertie2,Pbar,index);
    x0=X;
    sol(1,i+1)=X(1);
    sol(2,i+1)=X(3);
    sol(3,i+1)=X(5);
end
%-----------------------Plot the Binary VLE Diagram-----------------------
%--------------------------plotting bubble curve---------------------------
switch index
    case 1
        h(1)=plot(sol(1,:),sol(3,:),'bo',LineWidth=1.5);
    case 2
        h(2)=plot(sol(1,:),sol(3,:),'bx',LineWidth=1.5);
    case 3
        h(3)=plot(sol(1,:),sol(3,:),'b+',LineWidth=1.5);
    case 4
        h(4)=plot(sol(1,:),sol(3,:),'b*',LineWidth=1.5);
    case 5
        h(5)=plot(sol(1,:),sol(3,:),'bs',LineWidth=1.5);
    case 6
        h(6)=plot(sol(1,:),sol(3,:),'bd',LineWidth=1.5);
    case 7
        h(7)=plot(sol(1,:),sol(3,:),'bv',LineWidth=1.5);
end

    hold on
%---------------------------plotting due curve-----------------------------
switch index
    case 1
        h(8)=plot(sol(2,:),sol(3,:),'ro',LineWidth=1.5);
    case 2
        h(9)=plot(sol(2,:),sol(3,:),'rx',LineWidth=1.5);
    case 3
        h(10)=plot(sol(2,:),sol(3,:),'r+',LineWidth=1.5);
    case 4
        h(11)=plot(sol(2,:),sol(3,:),'r*',LineWidth=1.5);
    case 5
        h(12)=plot(sol(2,:),sol(3,:),'rs',LineWidth=1.5);
    case 6
        h(13)=plot(sol(2,:),sol(3,:),'rd',LineWidth=1.5);
    case 7
        h(14)=plot(sol(2,:),sol(3,:),'rv',LineWidth=1.5);
end

% Plot the Article Results
% h(15)=plot(A(:,2),A(:,1),'Ok',LineWidth=2);
% h(16)=plot(A(:,3),A(:,1),'Ok',LineWidth=2);


switch index
    case 1
        legend([h(1),h(8)],{'Original (LL)','Original (VL)'},FontName='times')
    case 2
        legend([h(2),h(9)],{'Khaled A. M. Gasem (LL)','Khaled A. M. Gasem (VL)'},FontName='times')
    case 3
        legend([h(3),h(10)],{'Chorng H. Twu (LL)','Chorng H. Twu (VL)'},FontName='times')
    case 4
        legend([h(4),h(11)],{'Saffari and Zahedi (LL)','Saffari and Zahedi (VL)'},FontName='times')
    case 5
        legend([h(5),h(12)],{'Georges A. Melhem (LL)','Georges A. Melhem (VL)'},FontName='times')
    case 6
        legend([h(6),h(13)],{'Mahmoodi and Sedigh PM2 (LL)','Mahmoodi and Sedigh PM2 (VL)'},FontName='times')
    case 7
        legend([h(7),h(14)],{'Mahmoodi and Sedigh PM3 (LL)','Mahmoodi and Sedigh PM3 (VL)'},FontName='times')
end
axis tight
prop1=PropertieNames(propertie1);
prop2=PropertieNames(propertie2);
xlabel([prop1+' vapor or liquid mole fraction'])
ylabel('temperature in K')
Pbar=int2str(Pbar);
title(["isobaric VLE diagram for "+prop1+"/"+prop2+" mixture at "+Pbar+" bar"])
grid on
shg