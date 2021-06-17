% Different gas composition and population ratio
% ...Wang group: May 2021...

clearvars;
close all;
clc
%% Experimental data (Gas composition exp.)(gDCW/L)

tt{1} = [0 11.375 23.29166667 35.41666667 48.5 71.91666667 96]; % time
tt{2} = [0 11.375 23.29166667 35.41666667 48.5 71.91666667 96];
tt{3} = [0 11.375 23.29166667 35.41666667 48.5 71.91666667 96];
tt{4} = [0 22.75 46.75 70.75 94]; %time

mm{1} = [0.017679 0.030483036 0.059890803 0.109108617 0.173973065 0.275477808 0.352663564]; % A % Methanotroph
mm{2} = [0.017679 0.03507124 0.075205667 0.127206728 0.193787397 0.308194375 0.380449242]; % B
mm{3} = [0.017679 0.032734772 0.089722134 0.19267495 0.297810934 0.436616768 0.5474978]; % C
mm{4} = [0.017679 0.110207105 0.209767561 0.315560355 0.385826233]; %D

mmStd{1} = [0 0.003547599 0.001233902 0.00275207 0.001255231 0.000298207 0.008551141];
mmStd{2} = [0 0.000752102 0.000342922 0.003890072 0.002626602 0.009344893 0.005478602];
mmStd{3} = [0 0.002635469 0.004257134 0.009899394 0.002460139 0.005217304 0.010789763];
mmStd{4} = [0 0.001951332 0.000894164 0.001885122 0.003564229];

aa{1} = [0.217225 0.315288749 0.418570692 0.576994989 0.79137755 1.120871373 1.383908997]; % A % Cyanobacteria
aa{2} = [0.217225 0.335948401 0.455636509 0.611512083 0.813753923 1.129469478 1.331724222]; % B
aa{3} = [0.217225 0.331142056 0.424598929 0.617572045 0.815117028 1.104361786 1.302238687]; % C
aa{4} = [0.217225 0.465019243 0.802410871 1.114492792 1.344836658]; %D

aaStd{1} = [0 0.00413915 0.018987169 0.011051697 0.015039582 0.01697042 0.011273322];
aaStd{2} = [0 0.016615659 0.011722354 0.013407165 0.017944794 0.019807368 0.021525667];
aaStd{3} = [0 0.005134249 0.021626812 0.046575595 0.012142364 0.017289911 0.062296375];
aaStd{4} = [0 0.01587751 0.015792079 0.035117572 0.040838102];

% Amount of gas consumption/production (every 24 hours)(mmol/L):
CO2pr24{1}=[0.096895855 0.169619144 0.231261414 0.285883926 0.407424635 0.298013742];
CO2pr24{2}=[0.105387013 0.188984094 0.224136357 0.320749619 0.453724417 0.330804162];
CO2pr24{3}=[0.184951283 0.297450856 0.436546881 0.46983158 0.545997593 0.443819152];
CO2pr24{4}=[0.323971462 0.525665813 0.488920498 0.301382765];

CO2co24{1}=[0.348962161 0.349111335 0.496276634 0.603729162 0.998255286 0.749910798];
CO2co24{2}=[0.376892309 0.448108991 0.522984833 0.748415778 1.022023639 0.771876378];
CO2co24{3}=[0.328613 0.331659517 0.56729103 0.57123617 0.78484893 0.530407332];
CO2co24{4}=[0.741957491 1.206873153 1.012279923 0.506164635];

O2pr24{1}=[0.44781785 0.456389169 0.645345792 0.784661452 1.288495824 0.965311359];
O2pr24{2}=[0.451115386 0.531203607 0.633731042 0.905656304 1.257389019 0.891242946];
O2pr24{3}=[0.41076625 0.422567516 0.669718577 0.687638246 0.926557764 0.626175323];
O2pr24{4}=[0.88544738 1.425935099 1.3069639 0.817014026];

O2co24{1}=[0.270822164 0.456389169 0.645345792 0.784661452 1.126094934 0.760917484];
O2co24{2}=[0.321531988 0.531203607 0.633731042 0.905656304 1.257389019 0.891242946];
O2co24{3}=[0.495274576 0.80309686 1.20409507 1.253675443 1.496289263 1.19501109];
O2co24{4}=[0.87531485 1.425935099 1.3069639 0.817014026];

CH4co24{1}=[0.20058726 0.34634061 0.483687126 0.597103488 0.86443306 0.588266706];
CH4co24{2}=[0.24395142 0.393716862 0.466950744 0.668228373 0.966092535 0.689175338];%mmol
CH4co24{3}=[0.38186893 0.60841941 0.907365888 0.965474721 1.161647342 0.9337899];
CH4co24{4}=[0.671328	1.0942595	1.0056785	0.627658];

% Gas measurement for every 24 hours(mmol/L):
measT = [0 11.375 11.375 23.29166667 23.29166667 35.41666667 35.41666667 48.5 48.5 71.91666667 71.91666667 96];
measCH4{1}= [7.64157 6.22137 7.654125 5.40516 7.68396 4.238895 7.675185 3.476145 7.66182 1.653645 7.625235 3.05589];
measCH4{2}= [24.36267 22.437165 24.40803 21.521325 24.50874 21.038835 24.70017 20.000955 24.62295 17.70538 24.613095 19.95917];
measCH4{3}= [24.425715 22.04466 24.243465 20.15479 24.60378 18.20775 24.69639 18.399045 24.615795 16.70493 24.52926 18.3495];
measCH4{4}= [32.3604  28.2228  32.3554	24.6505	32.3719	26.0772	32.3631	27.6395];

measCO2{1} = [4.71936445073079,3.01833383848589,4.73411935071058,2.67623149457803,4.73818966794639,2.43718188859702,4.71554852832222,2.12789647740284,4.71440375159965,2.12840526705732,4.70588152488718,2.13667309894255];% after adjusting the pH
measCO2{2} = [6.43831029837678,4.10989314339597,6.42113864753822,4.15100074089042,6.39837031050044,4.17644022361420,6.43373119148650,4.19984454772008,6.38679534586112,4.41612581666330,6.39391840102378,4.54947437192699];
measCO2{3} = [6.43182323028221,4.18417643294942,6.54146740082171,4.17423220179161,6.56220057924160,4.14767077524079,6.40027827170472,4.18926432949417,6.35993385869199,4.56583717249276,6.33780150872230,4.41523543476797];
measCO2{4} = [5.51979524483061,4.12851215733819,5.41365258974877,3.22846366269280,5.40767158348488,4.01567319997306,5.49963359601266,4.33173705125615];

measO2{1} = [0 0.846637904 0 0 0 0 0 0 0 0.2826726 0 0.4292925];
measO2{2} = [0 0.86388932 0 0 0 0 0 0 0 0 0 0];
measO2{3} = [4.13364 4.636917824 4.065264 0.34965 4.06371 0.292152 4.049724 0.254856 4.074588 0 4.068372 0];
measO2{4} = [0 0.134154 0 0 0 0 0 0];

%% parameters
% .....
% Insert Gas composition: % Rest of it will be filled with inert gas (Nitrogen)
N = 3; %Which experimental data set (1=20%CH4,10%Co2,; 2=60%CH4,30%CO2; 3=60%CH4,30%CO2,10%O2; 4=80%CH4,20%CO2)


if N==1
    CH4percentage = 0.2; %CH4 in gas
    CO2percentage = 0.1; %CO2 in gas
    O2percentage = 0.0; %O2 in gas
elseif N==2
    CH4percentage = 0.6; %CH4 in gas
    CO2percentage = 0.3; %CO2 in gas
    O2percentage = 0.0; %O2 in gas
elseif N==3
    CH4percentage = 0.6; %CH4 in gas
    CO2percentage = 0.3; %CO2 in gas
    O2percentage = 0.1; %O2 in gas
elseif N==4
    CH4percentage = 0.8; %CH4 in gas
    CO2percentage = 0.2; %CO2 in gas
    O2percentage = 0.0; %O2 in gas.
end
% Insert light intensity (umol/m2/s):
I0 = 180;
% Insert volume of liquid and gas phase (L):
Vl= 0.10;
Vg = 0.15; 

%microalgae Yields
YXpO2 = 24.54/1000; %
YO2CO2 = 1.28; %1.28-1.4
YXpCO2 = YXpO2*YO2CO2;

%methanotroph yields
YXmCO2 = 24.38/1000; 
YCO2CH4 = 0.48;%0.4-0.5
if N==3
YCO2CH4 = 0.40;
end
YXmCH4 = YXmCO2*YCO2CH4;
YO2CH4 = 1.4; %1.3-1.4
YXmO2 = YXmCH4/YO2CH4;

%Monod parameters
uPmax = 0.034; 
uMmax = 0.145; 
KsPco2 = 0.24;
KsMo2 = 0.005; 
KsMch4 = 0.028; 

% effective Henry's constants
eHCH4 = 0.03412; %effective Henry for Methane
eHO2 = 0.03168; %effective Henry for Oxygen
eHCO2 = 1.6120; %effective Henry for Carbon dioxide

% kLa
klaCH4=100;%80-200 (1/hr)
klaO2=1.17* klaCH4; % ref: (Yu, Ramsay and Ramsay 2006)
klaCO2=0.90*klaO2; % ref: (Babcock et al., 2002)

%A:Determining the effective light intensity
n = 1;
Ksi = 4.33; % when we use a constant Ksi

p = [CH4percentage CO2percentage O2percentage I0 Vl Vg]; % p(1)-p(6)
p = [p YXpO2 YO2CO2 YXpCO2]; % add microalgae yield parameters: p(7)-p(9)
p = [p YXmCO2 YCO2CH4 YXmCH4 YO2CH4 YXmO2]; % add methanotroph yield parameters: p(10)-p(14)
p = [p uPmax uMmax KsPco2 KsMo2 KsMch4]; % add Monod parameters: p(15)-p(19)
p = [p eHCH4 eHO2 eHCO2]; % effective Henry's constants: p(20)-p(22)
p = [p klaCH4 klaO2 klaCO2]; % kLa: p(23)-p(25)
p = [p n Ksi]; % light intensity: p(26)-p(27)
% time(hours)
ttt = [0 11.375 23.29166667 35.41666667 48.5 71.91666667 96];
%% initialization
% Insert initial individual biomass (inoculation):
X0P = 0.218;
X0M = 0.018;
% Insert total dissoved inorganic carbon in the liquid (mmol/L)(if any):
TIC = 200.7;%(mmol/L)
%Saturation Constant:
Hch4 = 0.0014; %Henry's constant (M/atm)
CH4l = (Hch4*CH4percentage*1.00*1000); % mmol/L
% CH4g = (CH4percentage*1.00*1000)/(0.08206*297); %mmol/L
CH4g = measCH4{N}(1);
Hco2 = 0.035; %Henry's constant (M/atm)
% % dissociation of CO2 due to pH change (%)
eta1=1;
CO2g = eta1* measCO2{N}(1);
CO2l = eta1*(Hco2*CO2percentage*1.00*1000); % mmol/L
% CO2g calculation using IC content and %CO2 in the gas phase:
% [CO2]g and TIC relationship for pH 8.7-9: TIC=alpha+beta[CO2]g
alpha = 170.81;
beta = 5.9235;
CO2g1=CO2g;
TIC1=TIC;
CO2g = (TIC1*Vl+CO2g1*Vg-alpha*Vl)/(beta*Vl+Vg);
Ho2 = 0.0013; %Henry's constant (M/atm)
O2l = (Ho2*O2percentage*1.00*1000);
O2g = (O2percentage*1.00*1000)/(0.08206*297);

% Initial conditions
% Y1 = Biomass Algae (gDW/L)
% Y2 = Biomass Methanotroph (gDW/L)
% Y3 = CH4g (mmol/L)
% Y4 = O2g (mmol/L) 
% Y5 = CO2g (mmol/L) 
% Y6 = CH4l (mmol/L) 
% Y7 = O2l (mmol/L) 
% Y8 = CO2l (mmol/L) 
Y0 = [X0P X0M CH4g O2g CO2g CH4l O2l CO2l]; %initial conditions 

nPara=21;
%% 0-11.375 hr
tspan=[0:0.001:ttt(2)];
[T1,Y1] = ode45(@(t,y)functionV7(t,y,p),tspan,Y0);%No2
Xp=Y1(:,1); % Y1 = Biomass Algae (gDW/L)
Xm=Y1(:,2);% Y2 = Biomass Methanotroph (gDW/L)
% Y3 = CH4g (mmol/L)
% Y4 = O2g (mmol/L) 
% Y5 = CO2g (mmol/L) 
CH4l=Y1(:,6);% Y6 = CH4l (mmol/L) 
O2l=Y1(:,7); % Y7 = O2l (mmol/L) 
CO2l=Y1(:,8);% Y8 = CO2l (mmol/L)
% methanotroph
uM_1 =uMmax.*(O2l./(KsMo2+O2l)).*(CH4l./(KsMch4+CH4l));
vMo2_1 = uM_1./YXmO2;
vMch4_1 = uM_1./YXmCH4;
vMco2_1 = uM_1./YXmCO2;
% photoautotroph
m = (-0.0175*I0+6.40);
Ia = I0.*exp(-m.*(Xp+Xm));
uP_1 = uPmax.*(CO2l./(KsPco2+CO2l)).*(Ia.^n./(Ksi.^n+Ia.^n));
vPco2_1 = uP_1./YXpCO2;
vPo2_1 = uP_1./YXpO2;

conPco2_1 = cumtrapz(T1,vPco2_1.*Xp)*Vl; % mmol
proPo2_1 = cumtrapz(T1,vPo2_1.*Xp)*Vl;
proMco2_1 = cumtrapz(T1,vMco2_1.*Xm)*Vl;
conMo2_1 = cumtrapz(T1,vMo2_1.*Xm)*Vl;
conMch4_1= cumtrapz(T1,vMch4_1.*Xm)*Vl;
CH4g_1 = (Y1(1,3)*Vg+Y1(1,6)*Vl)-conMch4_1;

conPco2{N}(1) = trapz(T1,vPco2_1.*Xp)*Vl; %mmol
proPo2{N}(1) = trapz(T1,vPo2_1.*Xp)*Vl;
proMco2{N}(1) = trapz(T1,vMco2_1.*Xm)*Vl;
conMo2{N}(1) = trapz(T1,vMo2_1.*Xm)*Vl;
conMch4{N}(1) = trapz(T1,vMch4_1.*Xm)*Vl;
conMch4_1p(1) = (Y1(1,3)-Y1(end,3))*Vg; 
proMco2_1p = (Y1(1,3)-Y1(end,3))*Vg*YCO2CH4;

%% 11.375-23.29
tspan=[0:0.001:ttt(3)-ttt(2)];
Y0(1)=Y1(end,1);
Y0(2)=Y1(end,2);
[T2,Y2] = ode45(@(t,y)functionV7(t,y,p),tspan,Y0);%No2
T2=T2+T1(end);
Xp=Y2(:,1); % Y1 = Biomass Algae (gDW/L)
Xm=Y2(:,2);% Y2 = Biomass Methanotroph (gDW/L)
% Y3 = CH4g (mmol/L)
% Y4 = O2g (mmol/L) 
% Y5 = CO2g (mmol/L) 
CH4l=Y2(:,6);% Y6 = CH4l (mmol/L) 
O2l=Y2(:,7); % Y7 = O2l (mmol/L) 
CO2l=Y2(:,8);% Y8 = CO2l (mmol/L)
% methanotroph
uM_2 =uMmax.*(O2l./(KsMo2+O2l)).*(CH4l./(KsMch4+CH4l));
vMo2_2 = uM_2./YXmO2;
vMch4_2 = uM_2./YXmCH4;
vMco2_2 = uM_2./YXmCO2;
% photoautotroph
m = (-0.0175*I0+6.40); 
Ia = I0.*exp(-m.*(Xp+Xm));
uP_2 = uPmax.*(CO2l./(KsPco2+CO2l)).*(Ia.^n./(Ksi.^n+Ia.^n));
vPco2_2 = uP_2./YXpCO2;
vPo2_2 = uP_2./YXpO2;

conPco2{N}(2) = trapz(T2,vPco2_2.*Xp)*Vl; %mmol
proPo2{N}(2) = trapz(T2,vPo2_2.*Xp)*Vl;
proMco2{N}(2) = trapz(T2,vMco2_2.*Xm)*Vl;
conMo2{N}(2) = trapz(T2,vMo2_2.*Xm)*Vl;
conMch4{N}(2) = trapz(T2,vMch4_2.*Xm)*Vl;
conMch4_1p(2) = (Y2(1,3)-Y2(end,3))*Vg; 

%% 24-36
tspan=[0:0.001:ttt(4)-ttt(3)];
Y0(1)=Y2(end,1);
Y0(2)=Y2(end,2);
[T3,Y3] = ode45(@(t,y)functionV7(t,y,p),tspan,Y0);%No2
T3=T3+T2(end);
Xp=Y3(:,1); % Y1 = Biomass Algae (gDW/L)
Xm=Y3(:,2);% Y2 = Biomass Methanotroph (gDW/L)
% Y3 = CH4g (mmol/L)
% Y4 = O2g (mmol/L) 
% Y5 = CO2g (mmol/L) 
CH4l=Y3(:,6);% Y6 = CH4l (mmol/L) 
O2l=Y3(:,7); % Y7 = O2l (mmol/L) 
CO2l=Y3(:,8);% Y8 = CO2l (mmol/L)
% methanotroph
uM_3 =uMmax.*(O2l./(KsMo2+O2l)).*(CH4l./(KsMch4+CH4l));
vMo2_3 = uM_3./YXmO2;
vMch4_3 = uM_3./YXmCH4;
vMco2_3 = uM_3./YXmCO2;
% photoautotroph
m = (-0.0175*I0+6.40); 
Ia = I0.*exp(-m.*(Xp+Xm));
uP_3 = uPmax.*(CO2l./(KsPco2+CO2l)).*(Ia.^n./(Ksi.^n+Ia.^n));
vPco2_3 = uP_3./YXpCO2;
vPo2_3 = uP_3./YXpO2;

conPco2{N}(3) = trapz(T3,vPco2_3.*Xp)*Vl; %mmol
proPo2{N}(3) = trapz(T3,vPo2_3.*Xp)*Vl;
proMco2{N}(3) = trapz(T3,vMco2_3.*Xm)*Vl;
conMo2{N}(3) = trapz(T3,vMo2_3.*Xm)*Vl;
conMch4{N}(3) = trapz(T3,vMch4_3.*Xm)*Vl;
conMch4_1p(3) = (Y3(1,3)-Y3(end,3))*Vg; 

%% 36-48
paramInBase_i=zeros(1,nPara);
paramInBase=zeros(1,nPara);
tspan=[0:0.001:ttt(5)-ttt(4)];
Y0(1)=Y3(end,1);
Y0(2)=Y3(end,2);
[T4,Y4] = ode45(@(t,y)functionV7(t,y,p),tspan,Y0);%No2
T4=T4+T3(end);
Xp=Y4(:,1); % Y1 = Biomass Algae (gDW/L)
Xm=Y4(:,2);% Y2 = Biomass Methanotroph (gDW/L)
% Y3 = CH4g (mmol/L)
% Y4 = O2g (mmol/L) 
% Y5 = CO2g (mmol/L) 
CH4l=Y4(:,6);% Y6 = CH4l (mmol/L) 
O2l=Y4(:,7); % Y7 = O2l (mmol/L) 
CO2l=Y4(:,8);% Y8 = CO2l (mmol/L)
% methanotroph
uM_4 =uMmax.*(O2l./(KsMo2+O2l)).*(CH4l./(KsMch4+CH4l));
vMo2_4 = uM_4./YXmO2;
vMch4_4 = uM_4./YXmCH4;
vMco2_4 = uM_4./YXmCO2;
% photoautotroph
m = (-0.0175*I0+6.40); 
Ia = I0.*exp(-m.*(Xp+Xm));
uP_4 = uPmax.*(CO2l./(KsPco2+CO2l)).*(Ia.^n./(Ksi.^n+Ia.^n));
vPco2_4 = uP_4./YXpCO2;
vPo2_4 = uP_4./YXpO2;

conPco2{N}(4) = trapz(T4,vPco2_4.*Xp)*Vl; %mmol
proPo2{N}(4) = trapz(T4,vPo2_4.*Xp)*Vl;
proMco2{N}(4) = trapz(T4,vMco2_4.*Xm)*Vl;
conMo2{N}(4) = trapz(T4,vMo2_4.*Xm)*Vl;
conMch4{N}(4) = trapz(T4,vMch4_4.*Xm)*Vl;
conMch4_1p(4) = (Y4(1,3)-Y4(end,3))*Vg; 

%% 48-72
paramInBase_i=zeros(1,nPara);
paramInBase=zeros(1,nPara);
tspan=[0:0.001:ttt(6)-ttt(5)];
Y0(1)=Y4(end,1);
Y0(2)=Y4(end,2);
[T5,Y5] = ode45(@(t,y)functionV7(t,y,p),tspan,Y0);%No2
T5=T5+T4(end);
Xp=Y5(:,1); % Y1 = Biomass Algae (gDW/L)
Xm=Y5(:,2);% Y2 = Biomass Methanotroph (gDW/L)
% Y3 = CH4g (mmol/L)
% Y4 = O2g (mmol/L) 
% Y5 = CO2g (mmol/L) 
CH4l=Y5(:,6);% Y6 = CH4l (mmol/L) 
O2l=Y5(:,7); % Y7 = O2l (mmol/L) 
CO2l=Y5(:,8);% Y8 = CO2l (mmol/L)
% methanotroph
uM_5 =uMmax.*(Y5(:,7)./(KsMo2+Y5(:,7))).*(Y5(:,6)./(KsMch4+Y5(:,6)));
vMo2_5 = uM_5./YXmO2;
vMch4_5 = uM_5./YXmCH4;
vMco2_5 = uM_5./YXmCO2;
% photoautotroph
m = (-0.0175*I0+6.40); 
Ia = I0.*exp(-m.*(Xp+Xm));
uP_5 = uPmax.*(CO2l./(KsPco2+CO2l)).*(Ia.^n./(Ksi.^n+Ia.^n));
vPco2_5 = uP_5./YXpCO2;
vPo2_5 = uP_5./YXpO2;

conPco2{N}(5) = trapz(T5,vPco2_5.*Xp)*Vl; %mmol
proPo2{N}(5) = trapz(T5,vPo2_5.*Xp)*Vl;
proMco2{N}(5) = trapz(T5,vMco2_5.*Xm)*Vl;
conMo2{N}(5) = trapz(T5,vMo2_5.*Xm)*Vl;
conMch4{N}(5) = trapz(T5,vMch4_5.*Xm)*Vl;
conMch4_1p(5) = (Y5(1,3)-Y5(end,3))*Vg; 

%% 72-96
paramInBase_i=zeros(1,nPara);
paramInBase=zeros(1,nPara);
tspan=[0:0.001:ttt(7)-ttt(6)];
Y0(1)=Y5(end,1);
Y0(2)=Y5(end,2);
[T6,Y6] = ode45(@(t,y)functionV7(t,y,p),tspan,Y0);%No2
T6=T6+T5(end);
Xp=Y6(:,1); % Y1 = Biomass Algae (gDW/L)
Xm=Y6(:,2);% Y2 = Biomass Methanotroph (gDW/L)
% Y3 = CH4g (mmol/L)
% Y4 = O2g (mmol/L) 
% Y5 = CO2g (mmol/L) 
CH4l=Y6(:,6);% Y6 = CH4l (mmol/L) 
O2l=Y6(:,7); % Y7 = O2l (mmol/L) 
CO2l=Y6(:,8);% Y8 = CO2l (mmol/L)
% methanotroph
uM_6 =uMmax.*(Y6(:,7)./(KsMo2+Y6(:,7))).*(Y6(:,6)./(KsMch4+Y6(:,6)));
vMo2_6 = uM_6./YXmO2;
vMch4_6 = uM_6./YXmCH4;
vMco2_6 = uM_6./YXmCO2;
% photoautotroph
m = (-0.0175*I0+6.40); 
Ia = I0.*exp(-m.*(Xp+Xm));
uP_6 = uPmax.*(CO2l./(KsPco2+CO2l)).*(Ia.^n./(Ksi.^n+Ia.^n));
vPco2_6 = uP_6./YXpCO2;
vPo2_6 = uP_6./YXpO2;

conPco2{N}(6) = trapz(T6,vPco2_6.*Xp)*Vl; %mmol
proPo2{N}(6) = trapz(T6,vPo2_6.*Xp)*Vl;
proMco2{N}(6) = trapz(T6,vMco2_6.*Xm)*Vl;
conMo2{N}(6) = trapz(T6,vMo2_6.*Xm)*Vl;
conMch4{N}(6) = trapz(T6,vMch4_6.*Xm)*Vl;
conMch4_1p(6) = (Y6(1,3)-Y6(end,3))*Vg; 

T=[T1;T2;T3;T4;T5;T6];
Y=[Y1;Y2;Y3;Y4;Y5;Y6];

%% Population Ratio Calculation
ratio = Y(:,1)./Y(:,2);
ttr{1} = [0 11.375 23.29166667 35.41666667 48.5 71.91666667 96];
ttr{2} = [0 11.375 23.29166667 35.41666667 48.5 71.91666667 96];
ttr{3} = [0 11.375 23.29166667 35.41666667 48.5 71.91666667 96];
rr{1}=[12.28717688 10.34308898 6.988897664 5.288262331 4.548850986 4.068826374 3.924162112];
rr{2}=[12.28717688 9.579028284 6.058539564 4.707230657 4.099209739 3.564795887 3.400399199];
rr{3}=[12.28717688 10.11591151 4.532376611 3.092064578 2.612857175 2.489802986 2.326416667];

%% Plots
% biomass prediction
figure;
plot(T,Y(:,1),'b',T,Y(:,2),'r')
hold on
plot(tt{N},aa{N},'bs',tt{N},mm{N},'ro')
errorbar(tt{N},aa{N},aaStd{N}, 'LineStyle','none')
errorbar(tt{N},mm{N},mmStd{N}, 'LineStyle','none')
% plot(tt2,aaa3,'bs',tt2,mmm2,'ro')
legend('Estimated Photoautotroph', 'Estimated Methanotroph','Measured Photoautotroph', 'Measured Methanotroph')
ylim([0,1.5])
ylabel('Biomass(g/L)')
xlabel('Time(hr)')

% population ratio
figure;
plot(T,ratio(:,1))
hold on
plot(ttr{N},rr{N},'v')
ylabel('Population Ratio (P/M)')
xlabel('Time(hr)')


% Change in Gas Phase
figure
plot(measT,measCH4{N},'-.o',measT,measCO2{N},'-.o',measT,measO2{N},'-.o');
hold on
plot(T,Y(:,3),T,Y(:,5),T,Y(:,4))%
ylabel('Concentration (mmol/L)')
xlabel('Time(hr)')
% legend('CH4g' ,'O2g', 'CO2g')% 
 
legend('measured CH4g','measured CO2g','measured O2g','predicted CH4g' ,'predicted CO2g', 'predicted O2g')% 


% plot for comparison of production/consumption every 12/24 hours
figure
% Bar(s)
yy = [CH4co24{N};conMch4{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-12','12-24','24-36','36-48','48-72','72-96'})
legend('Measured','Predicted')
ylabel('CH4 mmol')


figure
yy = [O2co24{N};conMo2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-12','12-24','24-36','36-48','48-72','72-96'})
legend('Measured','Predicted')
ylabel('O2 Consumed mmol')


figure
yy = [O2pr24{N};proPo2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-12','12-24','24-36','36-48','48-72','72-96'})
legend('Measured','Predicted')
ylabel('O2 Produced mmol')


figure
yy = [CO2pr24{N};proMco2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-12','12-24','24-36','36-48','48-72','72-96'})
legend('Measured','Predicted')
ylabel('CO2 Produced mmol')


figure
yy = [CO2co24{N};conPco2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-12','12-24','24-36','36-48','48-72','72-96'})
legend('Measured','Predicted')
ylabel('CO2 Consumed mmol')
% 
