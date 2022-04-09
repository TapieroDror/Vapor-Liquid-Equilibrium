function [L,M,N] = Twu (propNum)
param=[
    0.000000 0.000000 0.00000;%Air
    0.000000 0.000000 0.00000;%Ammonia
    0.036512 0.935460 3.97643;%Argon
    0.141223 0.844404 2.47802;%Benzene
    0.000000 0.000000 0.00000;%Bromine
    0.351718 0.843261 1.47039;%n-Butane
    0.000000 0.000000 0.00000;%Carbon dioxide
    0.000000 0.000000 0.00000;%Carbon monoxide
    0.000000 0.000000 0.00000;%Carbon tetrachloride
    0.000000 0.000000 0.00000;%Chlorine 
    0.000000 0.000000 0.00000;%Chloroform 
    0.000000 0.000000 0.00000;%(R-12)
    0.000000 0.000000 0.00000;%(R-21)
    0.311041 0.866279 1.29869;%Ethane
    0.000000 0.000000 0.00000;%Ethanol
    0.000000 0.000000 0.00000;%Ethene
    0.000000 0.000000 0.00000;%Helium
    0.088790 0.901528 5.08963;%n-Hexane
    0.000000 0.000000 0.00000;%Hydrogen
    0.000000 0.000000 0.00000;%Krypton
    0.081043 0.915696 2.61622;%Methane
    0.000000 0.000000 0.00000;%Methanol
    0.000000 0.000000 0.00000;%Methyl chloride
    0.000000 0.000000 0.00000;%Neon
    0.000000 0.000000 0.00000;%Nitrogen
    0.000000 0.000000 0.00000;%Nitrous oxide
    0.000000 0.000000 0.00000;%Oxygen
    0.343038 0.810237 1.55034;%n-Pentane
    0.283800 0.858285 1.59420;%Propane
    0.000000 0.000000 0.00000;%Propene
    0.000000 0.000000 0.00000;%Sulfur dioxide
    0.000000 0.000000 0.00000;%(R-134a)
    0.000000 0.000000 0.00000;%Toluene
    0.000000 0.000000 0.00000;%(R-11)
    0.000000 0.000000 0.00000;%Water
    0.000000 0.000000 0.00000;%Xenon
    0.169318 0.844860 3.03151;%n-Heptane 
    0.428179 0.805321 1.71898;%n-Octane 
    1.971720 1.628290 0.37765;%n-Nonane 
    0.315303 0.811589 2.42062;%n-Decane 
    0.246874 0.821219 3.03978;%n-Undecane 
    0.411574 0.802000 2.25086;%n-Dodecane 
    0.509689 0.822738 2.10994;%n-Tridecane 
    0.236127 0.821870 3.61671;%n-Tetradecane 
    0.229708 0.823364 3.85809;%n-Pentadecane 
    0.340390 0.810309 3.11841;%n-Hextadecane 
    0.134135 0.874867 6.76784;%n-Heptadecane 
    0.480721 0.783810 2.47646;%n-Octadecane 
    0.543018 0.772678 2.27430;%n-Nonadecane 
    2.731330 2.209970 0.30780;%Eicosane 
    0.076117 0.882790 3.80751;%Cyclohexane
    ];
L=param(propNum,1);
L=L';
M=param(propNum,2);
M=M';
N=param(propNum,3);
N=N';