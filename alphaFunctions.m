function alpha = alphaFunctions (index,w,Tr,propNum)
%References:
%[1]    Soave, Giorgio. "Equilibrium constants from a modified Redlich-Kwong equation of state." Chemical engineering science 27.6 (1972): 1197-1203.
%[2] 
%[6]&[7]    Mahmoodi, Peyman, and Maryam Sedigh. "Second derivative of alpha functions in cubic equations of state." The Journal of Supercritical Fluids 120 (2017): 191-206.‏
switch index
    case 1 % Soave (ORIGINAL)‏
        k=0.37464 + 1.54226.*w-0.26992.*w.^2;
        alpha=(1+k.*(1-Tr.^0.5)).^2;
%     case 2 % Milind Joshipura
%         J=1.252*w+0.4754;
%         alpha=exp(J.*(1-Tr));
    case 2 % Khaled A. M. Gasem
        G3=0.134+0.508*w-0.0467*w.^2;
        alpha=exp((2+0.836.*Tr).*(1-Tr.^G3));
    case 3 % Chorng H. Twu
        [L,M,N]=Twu(propNum);
        sum1=sum(L~=0);
        sum2=sum(M~=0);
        sum3=sum(N~=0);
        if((sum1+sum2+sum3)~=3*length(propNum))
            error('For one or more of the materials, there are no parameters for this alpha function.')
        else
            alpha=Tr.^(N.*(M-1)).*exp(L.*(1-Tr.^(N.*M)));
        end
    case 4 % Hamid Saffari and Alireza Zahedi (very good for Natural Gas)
        k1=0.013145*w+0.003091;
        k2=0.482173*w-0.006487;
        k3=3.586161*w+0.721306;
        alpha=exp(k1.*Tr+k2.*log(Tr)+k3.*(1-sqrt(Tr)));
    case 5 % Georges A. Melhem
        [m,n]= Melhem(propNum);
        sum1=sum(m~=0);
        sum2=sum(n~=0);
        if((sum1+sum2)~=2*length(propNum))
            error('For one or more of the materials, there are no parameters for this alpha function.')
        else
            alpha=exp(m.*(1-Tr)+n.*(1-sqrt(Tr)).^2);
        end
    case 6 % Peyman Mahmoodi and Maryam Sedigh (proposed model 2)
        [~,~,~,~,~,C1,C2,~,~,~] = MahmoodiAndSedigh(propNum);
        sum1=sum(C1~=0);
        sum2=sum(C2~=0);
        if((sum1+sum2)~=2*length(propNum))
            error('For one or more of the materials, there are no parameters for this alpha function.')
        else
            alpha=exp(2*C1.*(1-sqrt(Tr))-(C2.*(1-sqrt(Tr))).^2+2/3*(C1.*(1-sqrt(Tr))).^3);
        end
    case 7 % Peyman Mahmoodi and Maryam Sedigh (proposed model 3)
        [~,~,~,~,~,~,~,C1,C2,C3] = MahmoodiAndSedigh(propNum);
        sum1=sum(C1~=0);
        sum2=sum(C2~=0);
        sum3=sum(C3~=0);
        if((sum1+sum2+sum3)~=3*length(propNum))
            error('For one or more of the materials, there are no parameters for this alpha function.')
        else
            alpha=exp(2*C1.*(1-sqrt(Tr))-(C2.*(1-sqrt(Tr))).^2+2/3*(C3.*(1-sqrt(Tr))).^3);
        end
end
