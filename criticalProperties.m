%critical constants and acentric factor
%Mw - molar mass [kg/kmol]
%R - gas constant [kJ/kg*K]
%Pc - critical pressure [MPa]
%Tc - critical temperature [K]
%Vc - critical volume [m^3/kmol]
%w - acentric factor [-] 
%Hf - Enthalpy of Formation [J/mol] (at 25[C], 100[KPa]) page 239
%chart:
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
%[22        Methyl alcohol(Methanol)*          CH3OH    ]
%[23        Methyl chloride (Chloromethane)    CH3Cl    ]
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
function prop=criticalProperties(numberInChart)
propertiesMatrix=...
[28.970 0.2870 132.5 3.77 0.0883 0.0335 0;% Air
17.030 0.4882 405.5 11.28 0.0724 0.2530 -45898;% Ammonia
39.948 0.2081 151.0 4.860 0.0749 0.0010 0;% Argon
78.115 0.1064 562.0 4.920 0.2603 0.2090 82880;% Benzene
159.81 0.0520 584.0 10.34 0.1355 0.1280 30910;% Bromine
58.124 0.1430 425.2 3.800 0.2547 0.1970 -125790;% n-Butane
44.010 0.1889 304.2 7.390 0.0943 0.2240 -393510;% Carbon dioxide
28.011 0.2968 133.0 3.500 0.0930 0.0480 -110530;% Carbon monoxide 
153.82 0.0541 556.4 4.560 0.2759 0.1910 0;% Carbon tetrachloride  
70.906 0.1173 417.0 7.710 0.1242 0.0730 0;% Chlorine
119.38 0.0696 536.6 5.470 0.2403 0.2280 0;% Chloroform
120.91 0.0688 384.7 4.010 0.2179 0.2040 0;% Dichlorodifluoromethane(R-12)
102.92 0.0808 451.7 5.170 0.1973 0.2060 0;% Dichlorofluoromethane(R-21)
30.070 0.2765 305.5 4.872 0.1480 0.0980 -84740;% Ethane 
46.070 0.1805 516.0 6.380 0.1673 0.6640 -235000;% Ethyl alcohol (Ethanol)
28.054 0.2964 282.4 5.120 0.1242 0.0860 52467;% Ethylene (Ethene)
4.0030 2.0769 5.300 0.230 0.0578 -0.388 0;% Helium 
86.179 0.0965 507.9 3.032 0.3677 0.3000 -167300;% n-Hexane
2.0160 4.1240 33.30 1.300 0.0649 -0.215 0;% Hydrogen(normal)
83.800 0.0992 209.4 5.500 0.0924 0.0050 0;% Krypton
16.043 0.5182 191.1 4.640 0.0993 0.0110 -74873;% Methane
32.042 0.2595 513.2 7.950 0.1180 0.5560 -201300;% Methyl alcohol (Methanol)
50.488 0.1647 416.3 6.680 0.1430 0.1530 0;% Methyl chloride 
20.183 0.4119 44.50 2.730 0.0417 -0.038 0;% Neon 
28.013 0.2968 126.2 3.390 0.0899 0.0390 0;% Nitrogen 
44.013 0.1889 309.7 7.270 0.0961 0.1430 82050;% Nitrous oxide 
31.999 0.2598 154.8 5.080 0.0780 0.0200 0;% Oxygen 
72.150 0.1152 469.6 3.375 0.3150 0.2540 -146500;% n-Pentane 
44.097 0.1885 370.0 4.260 0.1998 0.1530 -103900;% Propane 
42.081 0.1976 365.0 4.620 0.1810 0.1370 20430;% Propylene (Propene)
64.063 0.1298 430.7 7.880 0.1217 0.2440 -296842;% Sulfur dioxide
102.03 0.0815 374.2 4.059 0.1993 0.3270 0;% Tetrafluoroethane(R-134a)
92.141 0.0902 591.8 4.100 0.3140 0.2620 0;% Toluene
137.37 0.0605 471.2 4.380 0.2478 0.1890 0;% Trichlorofluoromethane(R-11)
18.015 0.4615 647.1 22.06 0.0560 0.0080 -241826;% Water  
131.30 0.0633 289.8 5.880 0.1186 0.3430 0];%Xenon
Mw=propertiesMatrix(numberInChart,1);
R=propertiesMatrix(numberInChart,2);
Tc=propertiesMatrix(numberInChart,3);
Pc=propertiesMatrix(numberInChart,4);
Vc=propertiesMatrix(numberInChart,5);
w=propertiesMatrix(numberInChart,6);
Hf=propertiesMatrix(numberInChart,7);
prop=[Mw R Tc Pc Vc w Hf];