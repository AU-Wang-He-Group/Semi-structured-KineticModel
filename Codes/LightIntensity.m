%% Light intensity experiment/Model
% ...Wang group: May 2021...

clearvars;
% close all;
% clc
%% Experimental data (Light Intensity exp.)(gDCW/L)
% tt” shows the time, “mm” shows methanotroph biomass, and “aa” shows photoautotroph biomass 
tt{1} = [0 22.75 46.75 70.75 94];
tt{2} = [0 22.75 46.75 70.75 94];
tt{3} = [0 22.75 46.75 70.75 94];
tt{4} = [0 22.75 46.75 70.75 94]; %time
mm{1} = [0.0180235 0.051664965 0.096299958 0.151554201 0.20806785]; % methanotroph biomass (gDCW/L) at=60 umol/m2/hr
mm{2} = [0.0180235 0.068770359 0.147220883 0.219770961 0.284556418]; % methanotroph biomass (gDCW/L) at=100 umol/m2/hr
mm{3} = [0.0180235 0.09662902 0.206249947 0.308093869 0.349184709];      % methanotroph biomass (gDCW/L) at=140 umol/m2/hr
mm{4} = [0.017679 0.110207105 0.220661725 0.337445477 0.431468526]; % methanotroph biomass (gDCW/L) at=180 umol/m2/hr

mmStd{1} = [0 0.000457587 0.005560311 0.006315208 0.004431002]; %standard deviation
mmStd{2} = [0 0.002172932 0.002637890 0.001087146 0.011093641]; %standard deviation
mmStd{3} = [0 0.002811037 0.005277562 0.003303156 0.020272729]; %standard deviation
mmStd{4} = [0 0.001951332 0.000894164 0.001885122 0.035642293]; %standard deviation

aa{1} = [0.2220025 0.401825842 0.534477604 0.680298288 0.781827768]; % Photoautotroph biomass (gDCW/L) at=60 umol/m2/hr. 
aa{2} = [0.2220025 0.440423680 0.638798864 0.842569244 0.984226541]; % Photoautotroph biomass (gDCW/L) at=100 umol/m2/hr.
aa{3} = [0.2220025 0.469621049 0.747693536 1.006626101 1.16958353]; % Photoautotroph biomass (gDCW/L) at=140 umol/m2/hr.
aa{4} = [0.217225 0.465019243 0.802410871 1.114492792 1.334836658]; % Photoautotroph biomass (gDCW/L) at=180 umol/m2/hr.

aaStd{1} = [0 0.007558861 0.021086318 0.025406247 0.036286293]; %standard deviation
aaStd{2} = [0 0.007480141 0.000296108 0.018024975 0.055135276]; %standard deviation
aaStd{3} = [0 0.020103624 0.015858815 0.033594477 0.079984467]; %standard deviation
aaStd{4} = [0 0.029877510 0.037920792 0.071175724 0.108381021 ]; %standard deviation


% Amount of gas consumption/production (every 24 hours)(mmol):
CO2pr24{1}=[0.245468465	0.248454860	0.18548912	0.14123111];
CO2pr24{2}=[0.262900229 0.357714406 0.24304152 0.205368819];
CO2pr24{3}=[0.276481506 0.46135173 0.377612248 0.289897651];
CO2pr24{4}=[0.323971462 0.525665813 0.488920498 0.301382765];

CO2co24{1}=[0.520282178	0.518209234	0.358969085	0.223184078];
CO2co24{2}=[0.624455642 0.749615655 0.529131615 0.385122351];
CO2co24{3}=[0.701221346 0.969366723 0.749641883 0.5216605495];
CO2co24{4}=[0.741957491 1.206873153 1.012279923 0.506164635];

O2pr24{1}=[0.660901336 0.66033521 0.47869901 0.3894201];
O2pr24{2}=[0.719857489 0.977500352 0.6538711 0.550659056];
O2pr24{3}=[0.758654622 1.25317674 1.01534448 0.781587144];
O2pr24{4}=[0.88544738 1.425935099 1.3069639 0.817014026];

O2co24{1}=[0.652132155 0.66033521 0.47869901 0.3894201];
O2co24{2}=[0.711792334 0.977500352 0.6538711 0.550659056];
O2co24{3}=[0.74758775  1.25317674 1.01534448 0.781587144];
O2co24{4}=[0.87531485 1.425935099 1.3069639 0.817014026];

CH4co24{1}=[0.5032675	0.511325	0.366907	0.2934675];
CH4co24{2}=[0.547088	0.7473385	0.505007	0.4194205];
CH4co24{3}=[0.57501	    0.9644644	0.779703	0.600468];
CH4co24{4}=[0.671328	1.0942595	1.0056785	0.627658];
%Net change of gas consumption/production
%[CO2photo.,CO2meth.,O2photo.,O2meth.,CH4meth.]
%[CO2photo.,O2photo.,CH4meth., O2meth.,CO2meth.]
netgaschange{1} = [1.620644575	0.820643555 2.189355656	2.180586475 1.67496];
netgaschange{2} = [2.288325263	1.069024974 2.901887997	2.893822842 2.218854];
netgaschange{3} = [2.941890502	1.405343135 3.808762986	3.797696114 2.9196454];
netgaschange{4} = [3.467275202  1.639940538 4.435360405	4.425227875 3.398924];

% % Gas measurement (Headspace) for every 24 hours(mmol/L):
measT = [0 22.75 22.75 46.75 46.75 70.75 70.75 94 ]; %time of each sampling points
measCH4{1}= [32.3639  28.6654  32.3633	28.6483	32.3763	29.4769	32.3671	30.1286];
measCH4{2}= [32.3677  28.8381  32.3693	27.9026	32.3524	28.933	32.362	29.8109];
measCH4{3}= [32.3407  28.5987  32.3569  26.0721 32.3555 26.3929 32.3608 28.0352];
measCH4{4}= [32.3604  28.2228  32.3554	24.6505	32.3719	26.0772	32.3631	27.6395];
measCO2{1} = [5.59376978514178,4.12743315147841,5.64016299589143,4.16784535596417,5.56774432545295,4.34858894052671,5.54285040748973,4.42553377786758];
measCO2{2} = [5.5979524483061,4.12743315147841,5.57792820098337,4.00263549655454,5.55600458004984,4.34179969017310,5.54921532969624,4.45721694618441];
measCO2{3} = [5.53111066208661,4.12743315147841,5.55869199164815,3.96784535596417,5.50098336364249,4.13048427291709,5.56095507509935,4.40177140162996];
measCO2{4} = [5.51979524483061,4.12851215733819,5.41365258974877,3.22846366269280,5.40767158348488,4.01567319997306,5.49963359601266,4.33173705125615];
measO2{1} = [0 0.062288 0 0 0 0 0 0];
measO2{2} = [0 0.084565 0 0 0 0 0 0];
measO2{3} = [0 0.105643 0 0 0 0 0 0];
measO2{4} = [0 0.134154 0 0 0 0 0 0];


%% ODE

%% parameters
% Insert Gas composition: % Rest of it will be filled with inert gas (Nitrogen)
CH4percentage = 0.8; %CH4 in gas
CO2percentage = 0.2; %CO2 in gas
O2percentage = 0.0; %O2 in gas

%...
% Insert light intensity (umol/m2/s)(60/100/140/180):
N = 3; %Which experimental data set (light intensities: 1=60; 2=100; 3=140, 4=180)
if N==1
    I0=60;
elseif N==2
    I0=100;
elseif N==3
    I0=140;
elseif N==4
    I0=180;
end
% Insert volume of liquid and gas phase (L):
Vl= 0.10;
Vg = 0.15; 

%microalgae Yields
YXpO2 = 24.54/1000; %
YO2CO2 = 1.28; %1.1-1.4
YXpCO2 = YXpO2*YO2CO2;

%methanotroph yields
YXmCO2 = 24.38/1000; 
YCO2CH4 = 0.48;
YXmCH4 = YXmCO2*YCO2CH4;
YO2CH4 = 1.35; %1.3-1.4
YXmO2 = YXmCH4/YO2CH4;

%Monod parameters
uPmax = 0.034; % light intensity exp=0.034
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
Ksi = 3.33; % when we use a constant Ksi

p = [CH4percentage CO2percentage O2percentage I0 Vl Vg]; % p(1)-p(6)
p = [p YXpO2 YO2CO2 YXpCO2]; % add microalgae yield parameters: p(7)-p(9)
p = [p YXmCO2 YCO2CH4 YXmCH4 YO2CH4 YXmO2]; % add methanotroph yield parameters: p(10)-p(14)
p = [p uPmax uMmax KsPco2 KsMo2 KsMch4]; % add Monod parameters: p(15)-p(19)
p = [p eHCH4 eHO2 eHCO2]; % effective Henry's constants: p(20)-p(22)
p = [p klaCH4 klaO2 klaCO2]; % kLa: p(23)-p(25)
p = [p n Ksi]; % light intensity: p(26)-p(27)
% time(hours)
ttt = [0 22.75 46.75 70.75 94];
%% initialization
% Insert initial individual biomass (inoculation):
X0P = 0.222;
X0M = 0.018;
% Insert total dissoved inorganic carbon in the liquid (mmol/L)(if any):
TIC = 200.7;%(mmol/L)
%Saturation Constant:
Hch4 = 0.0014; %Henry's constant (M/atm)
CH4l = (Hch4*CH4percentage*1.00*1000); % mmol/L
% CH4g = (CH4percentage*1.00*1000)/(0.08206*297); %mmol/L
CH4g = measCH4{N}(1);
Hco2 = 0.035; %Henry's constant (M/atm)
% CO2g = (CO2percentage*1.00*1000)/(0.08206*297);
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
%% 0-22.75 hr
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

%% 22.75-46.75
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

%% 46.75-70.75
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

%% 70.75-94
paramInBase_i=zeros(1,nPara);
paramInBase=zeros(1,nPara);
tspan=[0:0.001:ttt(5)-ttt(4)];
Y0(1)=Y3(end,1);
Y0(2)=Y3(end,2);
[T4,Y4] = ode45(@(t,y)functionV7(t,y,p),tspan,Y0);
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

NetGas=[sum(conPco2{1,N}),sum(proMco2{1,N}),sum(proPo2{1,N}),sum(conMo2{1,N}),sum(conMch4{1,N})];
T=[T1;T2;T3;T4];
Y=[Y1;Y2;Y3;Y4];

%% Population Ratio Calculation
ratio(:,1) = Y(:,1)./Y(:,2);
% Population Ratio, Experimental data
ttr{1} = [0 22.75 46.75 70.75 94];
ttr{2} = [0 22.75 46.75 70.75 94];
ttr{3} = [0 22.75 46.75 70.75 94];
ttr{4} = [0 22.75 46.75 70.75 94];

rr{1}=[12.31739118	7.777530518	5.550133303	4.488811815	3.757561622];
rr{2}=[12.31739118	6.404266117	4.339050622	3.833851571	3.458809848];
rr{3}=[12.31739118	4.860041518	3.625181712	3.267270797	3.449702977];
rr{4}=[12.28717688	4.219503298	3.636384475	3.302734421	3.093705746];

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
legend('Estimated photoautotroph', 'Estimated Methanotroph','Measured photoautotroph', 'Measured Methanotroph')
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
 

legend('measured CH4g','measured CO2g','measured O2g','predicted CH4g' , 'predicted CO2g','predicted O2g')% 


% plot for comparison of production/consumption every 24 hours
figure
% Bar(s)
yy = [CH4co24{N};conMch4{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-23','23-47','47-71','71-94'})
legend('Measured','Predicted')
ylabel('CH4 mmol')


figure
yy = [O2co24{N};conMo2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-23','23-47','47-71','71-94'})
legend('Measured','Predicted')
ylabel('O2 Consumed mmol')


figure
yy = [O2pr24{N};proPo2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-23','23-47','47-71','71-94'})
legend('Measured','Predicted')
ylabel('O2 Produced mmol')


figure
yy = [CO2pr24{N};proMco2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-23','23-47','47-71','71-94'})
legend('Measured','Predicted')
ylabel('CO2 Produced mmol')


figure
yy = [CO2co24{N};conPco2{N}]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'0-23','23-47','47-71','71-94'})
legend('Measured','Predicted')
ylabel('CO2 Consumed mmol')
% 
%Gas Net Change
figure
yy = [netgaschange{N};NetGas]; 
b = bar(yy','FaceColor','flat');%
set(gca, 'XTickLabel',{'CO_2pro.','CO_2con.','O_2pro.','O_2con.','CH_4con.'})
legend('coculture (Mea.)','coculture (Pre.)')
ylabel('Consumption/production (mmol)')