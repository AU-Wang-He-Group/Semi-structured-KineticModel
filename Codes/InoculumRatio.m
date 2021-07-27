%% Inoculum Ratio experiment/Model
% ...Wang group: May 2021...

clearvars;
close all;
clc
%% Experimental data (Initial inoculum Ratio exp.)(gDCW/L)
% tt” shows the time, “mm” shows methanotroph biomass, and “aa” shows photoautotroph biomass 
tt{1} = [0 22.75 46.75 70.75 94];
tt{2} = [0 22.75 46.75 70.75 94];
tt{3} = [0 22.75 46.75 70.75 94];
tt{4} = [0 22.75 46.75 70.75 94];%time
mm{1} = [0.01819575 0.027123679 0.041503693 0.097494607 0.196882834]; % methanotroph biomass (gDCW/L) at ratio 1.5:1
mm{2} = [0.01819575 0.052293095 0.100581669 0.211556291 0.314165181];      % methanotroph biomass (gDCW/L) at ratio 4:1
mm{2} = [0.01819575 0.091072733 0.178128773 0.306920756 0.405720018]; % methanotroph biomass (gDCW/L) at ratio 8.5:1
mm{2} = [0.0180235  0.111752842 0.238617699 0.382078957 0.505804087]; % methanotroph biomass (gDCW/L) at ratio 12.5:1

mmStd{1} = [0 9.21854E-05 0.004898803 0.013932441 0.018542928]; %standard deviation
mmStd{2} = [0 0.002708581 0.004922449 0.002794844 0.002980078]; %standard deviation
mmStd{3} = [0 0.0184383   0.040245925 0.023223035 0.024588569]; %standard deviation
mmStd{4} = [0 0.01951332  0.015941648 0.01985122  0.033642293]; %standard deviation

aa{1} = [0.023395  0.0676678   0.123095789 0.269983025 0.555524497]; % Photoautotroph biomass (gDCW/L) at ratio 1.5:1 
aa{2} = [0.07117   0.181577513 0.335740593 0.627240323 0.925206505]; % Photoautotroph biomass (gDCW/L) at ratio 4:1
aa{3} = [0.1510225 0.356493527 0.601992805 0.909339909 1.197337925]; % Photoautotroph biomass (gDCW/L) at ratio 8.5:1
aa{4} = [0.2220025 0.486175566 0.835248658 1.170625584 1.40166692]; % Photoautotroph biomass (gDCW/L) at ratio 12.5:1

aaStd{1} = [0 0.003124351 0.032510945 0.06064455  0.070965952]; %standard deviation
aaStd{2} = [0 0.008954408 0.026135291 0.027521763 0.0144062]; %standard deviation
aaStd{3} = [0 0.00458466  0.057897276 0.042496708 0.066276703]; %standard deviation
aaStd{4} = [0 0.029877510 0.037920792 0.071175724 0.098381021 ]; %standard deviation


% Amount of gas consumption/production (every 24 hours)(mmol/L):
CO2pr24{1}=[0.03548	    0.080521003	0.19540433	0.340703418];
CO2pr24{2}=[0.11422265	0.233961015	0.436684868	0.450880338];
CO2pr24{3}=[0.3070674	0.415481265	0.591893618	0.431410725];
CO2pr24{4}=[0.323971462 0.525665813 0.488920498 0.301382765];

CO2co24{1}=[0.07565452	0.170313649	0.428535222	0.8275192];
CO2co24{2}=[0.3018682	0.532935761	0.952011506 1.008954703];
CO2co24{3}=[0.693878038	0.936986112	1.334828177	0.972909953];
CO2co24{4}=[0.741957491 1.206873153 1.012279923 0.506164635];

O2pr24{1}=[0.09691344	0.218171785	0.53672362	0.932032095];
O2pr24{2}=[0.3155321	0.64922071	1.174867395	1.24032097];
O2pr24{3}=[0.888857767	1.20027921	1.709914895	1.24629765];
O2pr24{4}=[0.88544738   1.425935099 1.3069639   0.817014026];

O2co24{1}=[0.0967214    0.218171785	0.53672362	0.932032095];
O2co24{2}=[0.3155321	0.64922071	1.174867395	1.24032097];
O2co24{3}=[0.8870836	1.20027921	1.709914895	1.24629765];
O2co24{4}=[0.87531485   1.425935099 1.3069639   0.817014026];

CH4co24{1}=[0.0744	    0.16782445	0.4097874	0.72156315];
CH4co24{2}=[0.242717	0.4932467	0.90374415	0.95640075];
CH4co24{3}=[0.682372	0.9232917	1.31531915	0.9586905];
CH4co24{4}=[0.671328	1.0942595	1.0056785	0.627658];
% % Gas measurement (Headspace) for every 24 hours(mmol/L):
measT = [0 22.75 22.75 46.75 46.75 70.75 70.75 94 ]; %time of each sampling points
measCH4{1}= [32.3083  31.8283	32.3313 31.0483	32.30883  29.50375	32.3636  26.4729];
measCH4{2}= [32.3395  30.4381	32.3109	29.2163	32.34783  26.6269	32.3613  25.9007];
measCH4{3}= [32.3359  28.2335	32.3643	26.2383	32.31403  23.8281	32.3662	 26.1811];
measCH4{4}= [32.3604  28.2228   32.3554	24.6505	32.3719	  26.0772	32.3631	 27.6395];
measCO2{1} = [5.353282352  3.880612649	5.302881249	3.918608204	5.234793215	4.088543322	5.211387953	4.160886858];
measCO2{2} = [5.294370244  3.956423188	5.226282211	3.918608204	5.172024558	3.771521313	5.228409962	3.651287997];
measCO2{3} = [5.283731489  3.91923513	5.277956165	3.875601536	5.272332823	3.581615943	5.283579507	4.060739176];
measCO2{4} = [5.519795244  4.128512157  5.413652589 3.228463662 5.407671583 4.015673199 5.499633596 4.3317370512];
measO2{1} = [0 0 0 0 0 0 0 0];
measO2{2} = [0 0 0 0 0 0 0 0];
measO2{3} = [0 0 0 0 0 0 0 0];
measO2{4} = [0 0.134154 0 0 0 0 0 0];

%% ODE

%% parameters
% Insert Gas composition: % Rest of it will be filled with inert gas (Nitrogen)
CH4percentage = 0.8; %CH4 in gas
CO2percentage = 0.2; %CO2 in gas
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

%methanotroph yields
YXmCO2 = 24.38/1000; 
YCO2CH4 = 0.48;
YXmCH4 = YXmCO2*YCO2CH4;
YO2CH4 = 1.35; %1.3-1.4
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
%...
N =1; %Which experimental data set (population gDCW Ratios ( photoautotroph:methanotroph )): 1=1.5:1; 2=3:1; 3=4:1, 4=8.5:1, 5=12.5:1)

% Insert initial individual biomass (inoculation):
if N==1 %1.5:1
    X0P = 0.0233;
    X0M = 0.018;
elseif N==2 %4:1
    X0P = 0.071;
    X0M = 0.018;
elseif N==3 %8.5:1
    X0P = 0.151;
    X0M = 0.018;
elseif N==4 %12.5:1
    X0P = 0.222;
    X0M = 0.018;
end

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
CO2g =  measCO2{N}(1);
CO2l = (Hco2*CO2percentage*1.00*1000); % mmol/L
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
% Population Ratio, Experimental data
ttr{1} = [0 22.75 46.75 70.75 94];
ttr{2} = [0 22.75 46.75 70.75 94];
ttr{3} = [0 22.75 46.75 70.75 94];
ttr{4} = [0 23.5 47.125 71.625 94.167];
rr{1}=[1.285739802 2.494786927 2.965899664 2.769209847 2.821599461];
rr{2}=[3.911352926 3.472303803 3.337989868 2.964886182 2.944968317];
rr{3}=[8.299877719 3.914382673 3.379537146 2.962784017 2.95114333];
rr{4}=[12.31739118 4.350453709 3.500363393 3.063831603 3.071165668];

%% Computing the predicted gas productiion/consumption (Overall)
%% Plots
% biomass prediction
figure;
plot(T,Y(:,1),'b',T,Y(:,2),'r')
hold on
plot(tt{N},aa{N},'bs',tt{N},mm{N},'ro')
errorbar(tt{N},aa{N},aaStd{N}, 'LineStyle','none')
errorbar(tt{N},mm{N},mmStd{N}, 'LineStyle','none')
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
ylim([0,14])

% Change in Gas Phase
figure
plot(T,Y(:,3),T,Y(:,4),T,Y(:,5))%
ylabel('Concentration (mmol/L)')
xlabel('Time(hr)')
% legend('CH4g' ,'O2g', 'CO2g')% 
hold on 
plot(measT,measCH4{N},'-.o',measT,measCO2{N},'-.o',measT,measO2{N},'-.o');
legend('predicted CH4g' ,'predicted O2g', 'predicted CO2g', 'measured CH4g','measured CO2g','measured O2g')% 


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
