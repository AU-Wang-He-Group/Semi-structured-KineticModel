% Semi-structured Kinetic model/Function to take care of in situ production/consumption of the strains in the
% coculture system
% ...Wang group May 2021...

function dydt=f(t,y,p)

Xp = y(1);
Xm = y(2);
%Gas:
CH4g = y(3);
O2g = y(4);
CO2g = y(5);
%Liquid:
CH4l = y(6);
O2l = y(7);
CO2l = y(8);
%parameters:
CH4percentage = p(1);
CO2percentage = p(2);
O2percentage = p(3);

%light intensity:
I0 = p(4); %(umol/m2/s)

%Volume of gas and liquid phase:
Vl=p(5); %L
Vg=p(6); %L

%microalgae Yields
YXpO2 = p(7);
YO2CO2 = p(8);
YXpCO2 = p(9);

%methanotroph yields
YXmCO2 = p(10);
YCO2CH4 = p(11);
YXmCH4 = p(12);
YO2CH4 = p(13);
YXmO2 = p(14);

%Monod parameters
uPmax = p(15);
uMmax = p(16);
KsPco2 = p(17);
KsMo2 = p(18);
KsMch4 = p(19);

% %Michaelis menten parameters
% vPmaxco2= 3.97;% 
% vMmaxo2= 1.9;%
% vMmaxch4= 0.039;%
% KmPco2 = 58.63;%
% KmMo2 = 0.037;%
% KmMch4 = 0.00072;

%Saturation Constant:
Vratio = Vl/Vg;
eHch4 = p(20); %effective Henry for Methane
eHo2 = p(21); %effective Henry for Oxygen
eHco2 = p(22); %effective Henry for Carbon dioxide

%kLa's
klaCH4 = p(23);%80-200 (1/hr)
klaO2 = p(24); % ref: (Yu, Ramsay and Ramsay 2006)
klaCO2 = p(25); % ref: (Babcock et al., 2002)

%A:Determining the effective light intensity
n = p(26);
Ksi = p(27); % when we use a constant Ksi
m = (-0.0175*I0+6.40); %
Ia = I0*exp(-m*(Xp+Xm));%attenuated light intensity

%Determining the photoautotroph growth
uP = ((uPmax*CO2l)/(KsPco2+CO2l))*(Ia^n/(Ksi^n+Ia^n));
vPco2 = uP/YXpCO2;
vPo2 = uP/YXpO2;
%comparison of carbon dioxide consumption between Michaelis-Menten eq. and through cell growth
% vPco22 = (vPmaxco2*CO2l)/(KmPco2+CO2l);
% if vPco22>vPco2
%     vPco2=vPco22;
% end
%Detirmine Methanotroph growth
uM =((uMmax*O2l)/(KsMo2+O2l))*(CH4l/(KsMch4+CH4l));
vMo2 = uM/YXmO2;
vMch4 = uM/YXmCH4;
vMco2 = uM/YXmCO2;
%comparison of oxygen consumption between Michaelis-Menten eq. and through cell growth
% vMo22 = (vMmaxo2*O2l)/(KmMo2+O2l);
% if vMo22>vMo2
%     vMo2=vMo22;
% end
% vMch4 = (vMmaxch4*CH4l)/(KmMch4+CH4l);

% Update the variables:
dydt = zeros(length(y),1);    % a column vector
dydt(1) = uP*Xp;
dydt(2) = uM*Xm;
dydt(3) = -klaCH4*((eHch4*CH4g)-CH4l)*Vratio; %CH4g
dydt(4) = -klaO2*((eHo2*O2g)-O2l)*Vratio; %O2g
dydt(5) = -klaCO2*((eHco2*CO2g)-CO2l)*Vratio; %CO2g
    
dydt(6) = klaCH4*((eHch4*CH4g)-CH4l)- vMch4*Xm; %CH4l
dydt(7) = klaO2*((eHo2*O2g)-O2l) - vMo2*Xm + vPo2*Xp; %O2l
dydt(8) = klaCO2*((eHco2*CO2g)-CO2l) - vPco2*Xp + vMco2*Xm; %CO2l

end



