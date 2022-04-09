function [m,n] = Melhem (propNum)
param=[
 0.0000  0.0000%AIR
 0.7393  0.2837%AMMONIA
 0.3804  0.2109%ARGON
 0.6670  0.4720%BENZENE
 0.5082  0.5011%BROMINE
 0.6543  0.4308%NORMAL BUTANE
 0.6877  0.3813%CARBON DIOXIDE
 0.4489  0.3207%CARBON MONOXIDE
 0.6433  0.4615%CARBON TETRACBLORIDE
 0.5014  0.3624%CHLORINE
 0.6705  0.4636%CHLOROFORM
 0.0000  0.0000%(R-12)
 0.0000  0.0000%(R-21)
 0.5336  0.2431%ETHANE
 1.2417  0.1381%ETHANOL
 0.4873  0.4570%ETHYLENE
 0.0000  0.0000%Helium
 0.7939  0.4116%NORMAL HEXANE
 0.0829 -0.4780%HYDROGEN
 0.0000  0.0000%Krypton
 0.4045  0.1799%METHANE
 1.2107 -0.5292%METHANOL
 0.6068  0.3011%METHYL CHLORIDE
 0.0000  0.0000%Neon
 0.4274  0.3359%NITROGEN
 0.5981  0.4349%NITROUS OXIDE
 0.4137  0.2784%OXYGEN
 0.7362  0.3513%NORMAL PENTANE
 0.5909  0.3997%PROPANE
 0.5743  0.4425%PROPYLENE
 0.7238  0.4537%SULFUR DIOXIDE
 0.0000  0.0000%(R-134a)
 0.7563  0.3125%TOLUENE
 0.0000  0.0000%(R-11)
 0.8893  0.0151%WATER
 0.0000  0.0000%Xenon (END)
 0.8313 -0.5701%1-1-DICHLOROETHANE
 0.7311  0.3594%1-2-DICHLOROETHANE
 0.6644  0.2885%1-3-BUTADIENE
 0.6396  0.4639%1-BUTENE
 0.7303  0.3486%2-3-DIMETHYL BUTANE
 1.0223  1.7081%2-BUTANOL
 0.7825  0.1751%2-METHYL PENTANE
 0.7366  0.0216%3-METHYL 1-BUTENE
 0.8317 -0.0149%ACETALDEHYDE
 1.0791 -0.5313%ACETIC ACID
 2.1183 -6.6355%ACETIC ANHYDRIDE
 0.8283  0.1495%ACETONE
 0.9813  0.9529%ACETONITRILE
 0.6640  0.2062%ACETYLENE
 0.9016 -0.2409%ACRYLONITRILE
 0.9227  0.1891%ANILINE
 0.7979  0.4790%BENZALDEHYDE
 0.7357  1.8414%BENZYL ALCOHOL
 0.7497  0.2442%BROMOBENZENE
 0.5924  0.1033%CARBON DISULFIDE
 0.7511  0.1902%CHLOROBENZENE
 0.6651  0.4814%CYCLOHEXANE
 0.6366  0.5328%CYCLOPENTANE
 0.6955  0.2458%DICHLORO METHANE
 2.4837 -6.5968%DIETHYLENE GLYCOL
 0.7712  0.6816%DIETHYL AMINE
 0.6572  0.3478%DIMETIYL ETHER
 2.2459 -4.9603%ETHYLENE GLYCOL
 0.6499  0.4672%ETHYLENE OXIDE
 0.8758  0.4057%ETHYL ACETATE
 0.8084  0.3129%ETHYL BENZENE
 0.7740  0.3761%ETHYL ETHER
 0.6495  0.3559%ETHYL MERCAPTAN
 0.4516  0.3038%FLUORINE
-0.1603 -0.9866%HELIUM-4
 0.6088  0.6640%HYDROGEN BROMIDE
 0.6018 -0.0507%HYDROGEN CHLORIDE
 1.1635 -2.0723%HYDROGEN CYANIDE
 1.1982 -2.8266%HYDROGEN FLUORIDE
 0.5193 -0.1840%HYDROGEN IODIDE
 0.5050  0.4719%HYDROGEN SULFID
 0.5024  0.5590%IDOINE
 0.6384  0.3959%ISOBUTANE
 0.6576  0.3844%ISOBUTYLENE
 0.7052  0.3361%ISOPENTANE
 1.1979  0.8456%ISOPROPANOL
 0.8718  0.1339%ISOPROPYL ETHER
 1.0553  1.5978%ISO BUTANOL
 0.8555  0.1447%METHYL ACETATE
 0.7699  0.5226%METHYL AMINE
 0.8455  0.2145%METHYL ETHYL KETONE
 0.6730  0.2581%METHYL FLUORIDE
 0.7202  0.5485%METHYL FORMATE
 0.5815  0.3936%METHYL IODIDE
 9.8477  0.2164%M-XYLENE
 1.3541 -1.5199%NITRIC OXIDE
 1.8386 -3.6754%NITROGEN DIOXIDE
 0.9585  0.6643%NITROMETHANE
 1.0705  1.3402%NORMAL BUTANOL
 1.0763  0.0164%NORMAL DECANE
 1.1778  0.0531%NORMAL DODECANE
 0.8620  0.3973%NORMAL HEPTANE
 1.0203  0.0654%NORMAL NONANE
 0.9273  0.3748%NORMAL OCTANE
 1.2243  0.1019%NORMAL TRIDECANE
 1.1394 -0.0157%NORMAL UNDECANE
 0.6838  0.2261%OZONE
 0.8283  0.2267%O-XYLENE
 1.0205 -0.0250%PHENOL
 1.1505  0.8075%PROPANOL
 0.7630  0.1883%PROPIONALDEHYDE
 0.7124  0.4424%PYRIDINE
 0.8478  0.1643%P-XYLENE
 0.4072  6.1502%SULFUR TRIOXIDE
 1.0636  1.7938];%TERT-BUTANOL
 
m=param(propNum,1);
m=m';
n=param(propNum,2);
n=n';