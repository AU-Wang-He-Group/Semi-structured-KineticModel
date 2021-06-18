% Coculture vs single culture experiment/Model
% ...Wang group: May 2021...

clearvars;
close all;
clc
%% Experimental data (coculture vs single culture)(gDCW/L)({1}coculture, {2}single photoautotroph, {3}single methanotroph)

tt{1} = [0 23.5 47.125 71.625 94.16666667];
tt{2} = [0 23.5 47.125 71.625 94.16666667];
tt{3} = [0 23.5 47.125 71.625 94.16666667];
mm{1} = [0.017679 0.110207105 0.220661725 0.337445477 0.431468526]; % methanotroph in coculture
aa{1} = [0.217225 0.465019243 0.802410871 1.114492792 1.334836658]; % photoautotroph in coculture
% Biomass for single culture:
aa{2}=[0.217225 0.3985425 0.6346875 0.93385 1.19047 ]; %photoautotroph single culture 
mm{2}=[0 0 0 0 0];
mm{3}=[0.017679 0.07331575 0.148876083 0.224379 0.292360333];% methanotroph single culture
aa{3}=[0 0 0 0 0];

aaStd{1} = [0 0.02987751 0.037920792 0.071175724 0.098381021];
mmStd{1} = [0 0.004951332 0.00894164 0.009885122 0.015642293];
aaStd{2} = [0 0.011489313 0.018526938 0.016904796 0.018245443];
mmStd{2} = [0 0 0 0 0];
aaStd{3} = [0 0 0 0 0];
mmStd{3} = [0 0.003042542 0.0056857 0.008250038 0.007957544];
% Amount of gas consumption/production (every 24 hours)(mmol):
CO2pr24{1} = [0.366061955 0.514858516 0.511918214 0.33825995];
CO2pr24{2} = [];
CO2pr24{3} = [];

CO2co24{1} = [0.8043522256 1.215077062 1.196634677 0.653420507];
CO2co24{2} = [];
CO2co24{3} = [];

O2pr24{1} = [1.088040282 1.394903514 1.388308027 0.890541176];
O2pr24{2} = [];
O2pr24{3} = [];

O2co24{1} = [0.9940440282 1.394903514 1.388308027 0.890541176];
O2co24{2} = [];
O2co24{3} = [];

CH4co24{1} = [0.768462406 1.072621908 1.0649628 0.68970823];
CH4co24{2} = [];%mmol
CH4co24{3} = [];

% Gas measurement for every 24 hours(mmol/L):
measT = [0 23.5 23.5 47.125 47.125 71.625 71.625 94.16666667];
measCH4{1} = [29.36982 24.11619 29.2143 22.00431 29.11048 22.46308 29.35205 24.65254];
measCH4{2} = [29.12007 28.6839102 28.74099 28.64991 28.32996 28.57323 28.41465 28.03125];
measCH4{3} = [23.67658785 19.5429 23.45412963 18.88293 23.1760299 18.97176 23.1874905 18.76386];

measCO2{1} = [6.04288940526706,4.88119869340608,6.05221721492557,4.60504984171887,6.05866188455580,4.63091331575402,6.03220482252307,4.75353162248266];% after adjusting the pH
measCO2{2} = [6.059193 4.220185667 6.00535633 3.843325 6.037251 3.971323333 5.906315 4.578161333];
measCO2{3} = [6.002546 4.325434511 6.08452123 4.958744 6.119547 4.955841227 6.091489 5.987564238];

measO2{1} = [0 0.223776 0 0 0 0 0 0];
measO2{2} = [0 4.80409776 0 6.285412 0 7.84252 0 6.67184];
measO2{3} = [5.01015816 0.0518 5.646073608 0 6.4299858 0.08806 5.8531928 0];

%% parameters
N = 1; %Which experimental data set (1=coculture; 2=photoautotroph; 3=methanotroph)

% Insert Gas composition: % Rest of it will be filled with inert gas (Nitrogen)
CH4percentage = 0.7; %CH4 in gas
CO2percentage = 0.3; %CO2 in gas
O2percentage = 0.0; %O2 in gas
% Insert light intensity (umol/m2/s):
I0 = 180;
% Insert volume of liquid and gas phase (L):
Vl= 0.10;
Vg = 0.15; 

%photoautotroph Yields
YXpO2 = 24.54/1000; %
YO2CO2 = 1.28; %1.1-1.4
YXpCO2 = YXpO2*YO2CO2;

%Methanotroph yields
YXmCO2 = 24.38/1000; 
YCO2CH4 = 0.48;
YXmCH4 = YXmCO2*YCO2CH4;
YO2CH4 = 1.35; %1.3-1.4
YXmO2 = YXmCH4/YO2CH4;

%Monod parameters
uPmax = 0.034; % light intensity of photoautotroph in coculture exp=0.034, single culture=0.024
if N==2
uPmax = 0.024;
end
uMmax = 0.145; % light intensity of methanotroph in coculture exp=0.145, single culture=0.117
if N==3
uPmax = 0.024;
uMmax = 0.117;
end
KsPco2 = 0.24;
KsMo2 = 0.005; 
KsMch4 = 0.028; 

% effective Henry's constants
eHCH4 = 0.03412; %effective Henry for Methane
eHO2 = 0.03168; %effective Henry for Oxygen
eHCO2 = 1.6120; %0.8530; %effective Henry for Carbon dioxide

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
ttt = [0 23.5 47.125 71.625 94.16666667];
%% initialization

% Insert initial individual biomass (inoculation):
X0P = 0.218;
X0M = 0.018;
if N==2
X0M = 0;
end
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
alpha = 170.26;
beta = 5.3595;
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
%% 0-23.5 hr
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
% metahnotroph
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

%% 23.5-47.125
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

%% 47.125-71.625
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

%% 71.625-94.16666667
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


T=[T1;T2;T3;T4];
Y=[Y1;Y2;Y3;Y4];

%% Population Ratio Calculation
ratio(:,1) = Y(:,1)./Y(:,2);
ttr{1} = [0 23.5 47.125 71.625 94.16666667];
ttr{2} = [0 23.5 47.125 71.625 94.16666667];
ttr{3} = [0 23.5 47.125 71.625 94.16666667];
rr{1}=[12.28717688 4.219503298 3.636384478 3.302734418 3.093705747];
rr{2}=[0 0 0 0 0];
rr{3}=[0 0 0 0 0];


%% Computing the predicted gas productiion/consumption (Overall)
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
plot(T,Y(:,3),T,Y(:,4),T,Y(:,5))%
ylabel('Concentration (mmol/L)')
xlabel('Time(hr)')
% legend('CH4g' ,'O2g', 'CO2g')% 
hold on 
plot(measT,measCH4{N},'-.o',measT,measCO2{N},'-.o',measT,measO2{N},'-.o');
legend('predicted CH4g' ,'predicted O2g', 'predicted CO2g', 'measured CH4g','measured CO2g','measured O2g')% 


% % plot for comparison of production/consumption every 24 hours
figure
% Bar(s)
yy = [CH4co24{N};conMch4{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-24','24-48','48-72','72-95'})
legend('Measured','Predicted')
ylabel('CH4 Consumed mmol')


figure
yy = [O2co24{N};conMo2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-24','24-48','48-72','72-95'})
legend('Measured','Predicted')
ylabel('O2 Consumed mmol')


figure
yy = [O2pr24{N};proPo2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-24','24-48','48-72','72-95'})
legend('Measured','Predicted')
ylabel('O2 Produced mmol')


figure
yy = [CO2pr24{N};proMco2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-24','24-48','48-72','72-95'})
legend('Measured','Predicted')
ylabel('CO2 Produced mmol')


figure
yy = [CO2co24{N};conPco2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-24','24-48','48-72','72-95'})
legend('Measured','Predicted')
ylabel('CO2 Consumed mmol')
% 
