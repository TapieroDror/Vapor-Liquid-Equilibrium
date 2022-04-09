% Liquid and Vapor Ternary Diagram calculator using Peng-Robinson
% Equation of State and Van der Waals Mixing Rules
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% <><><><><>Do not Copy - All Rights Reserved to Tapiero Dror ©<><><><><><>
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
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

%*[___________Hf for gas___________]
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
%-------------------------------plot 1 guide-------------------------------
%
%                                   1  0
%                                 0.9  0.1
%                               0.8      0.2
%                             0.7          0.3
%                           0.6              0.4
%                         0.5                  0.5
%                       0.4                      0.6
%                     0.3                          0.7
%                   0.2                              0.8
%                 0.1                                  0.9
%                0.0                                     1.0
%                1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.0
%
% 1: ->
%    |
% 2: v
% 3: <-

%-------------------------------plot 2 guide-------------------------------
%
%                                   0  1
%                                 0.1  0.9
%                               0.2      0.8
%                             0.3          0.7
%                           0.4              0.6  
%                         0.5                  0.5
%                       0.6                      0.4
%                     0.7                          0.3
%                   0.8                              0.2
%                 0.9                                  0.1
%                1.0                                     0.0
%                0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
%
% 1: <-
%    |
% 2: v
% 3: ->
%--------------------------------------------------------------------------
clc;clear;close all;
%---------------------------------inputs:----------------------------------
propNum=[7 29 25];
P=2;%[MPa]
T=240;%[K]
index=1;
plotNum=1;
k=7;%number of iterations -> k^3
%------------------------------Article Data:-------------------------------
%Conditions: T[K] and P[MPa]
Ax=[0.0294 0.0048 0.9658
    0.0284 0.3539 0.6177
    0.0265 0.4229 0.5506
    0.0238 0.5375 0.4387
    0.0211 0.6317 0.3472
    0.0175 0.7413 0.2412
    0.0266 0.9238 0.0496
    0.0130 0.9870 0.0000
    ];
Ay=[0.8953 0.0094 0.0951
    0.5205 0.4037 0.0758
    0.4849 0.4424 0.0727
    0.4325 0.4968 0.0707
    0.4030 0.5374 0.0596
    0.3507 0.6092 0.0401
    0.8216 0.0859 0.0925
    0.3186 0.6814 0.0000
    ];
%--------------------------------------------------------------------------
count=1; %Each time we find a result, count will increase by one
resultx=ones(1,3); %result matrix for x
resulty=ones(1,3); %result matrix for y
resultz=ones(1,3); %result matrix for z
%-------------------Display the names of the Substances--------------------
options = optimset('Display','off');
for i=1:length(propNum)
    fprintf(PropertieNames(propNum(i))+'   ')
end
fprintf('\n')
%--------------------------------------------------------------------------
%All possible mole fraction combinations
for n1=1:k
    for n2=1:k
        for n3=1:k
            nmoles=[n1 n2 n3];

            %Calculate mole fractions
             z=zeros(1,length(nmoles));
             nsum=sum(nmoles);
             for i=1:length(nmoles)
                 z(i)=nmoles(i)/nsum;
             end
                
                %Guess for the molar fractions of the liquid for each 
                % component (will be used to solve the Rachford-Rice equations)
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
                if ~isreal(x(1))||~isreal(x(2))||~isreal(x(3))||~isreal(y(1))||~isreal(y(2))||~isreal(y(3))
                    continue
                else
                    Ix=sum(resultx(:,1)==x(1)&resultx(:,2)==x(2)&resultx(:,3)==x(3));
                    Iy=sum(resulty(:,1)==y(1)&resulty(:,2)==y(2)&resulty(:,3)==y(3));
                    if(Ix>0||Iy>0)
                        continue
                    end
                    resultz(count,1)=z(1);
                    resultz(count,2)=z(2);
                    resultz(count,3)=z(3);
                    resultx(count,1)=x(1);
                    resultx(count,2)=x(2);
                    resultx(count,3)=x(3);
                    resulty(count,1)=y(1);
                    resulty(count,2)=y(2);
                    resulty(count,3)=y(3);
                    count=count+1;
                end
        end
    end
end

%-----------------------Plot the Ternary VLE Diagram-----------------------
switch plotNum
    case 1
        terplot1
        terlabel1(PropertieNames(propNum(1)),PropertieNames(propNum(2)),PropertieNames(propNum(3)));
        ternarycx(resultx(2:end,1),resultx(2:end,2),resultx(2:end,3));
        ternarycy(resulty(2:end,1),resulty(2:end,2),resulty(2:end,3));
%         ternarycArticle(Ax(:,2),Ax(:,3),Ax(:,1));
%         ternarycArticle(Ay(:,2),Ay(:,3),Ay(:,1));
    case 2
        terplot2
        terlabel2(PropertieNames(propNum(1)),PropertieNames(propNum(2)),PropertieNames(propNum(3)));
        ternarycx(resultx(2:end,1),resultx(2:end,2),resultx(2:end,3));
        ternarycy(resulty(2:end,1),resulty(2:end,2),resulty(2:end,3));
%         ternarycArticle(Ax(:,3),Ax(:,1),Ax(:,2));
%         ternarycArticle(Ay(:,3),Ay(:,1),Ay(:,2));
end
%===================================END====================================