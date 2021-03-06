function [Vis] = GlycerolViscosity(Pw,Temp)
%function [Vis,Pv] = GlycerolViscosity(Pw,Temp)

% function [Vis,Pw] = GlycerolViscosity(Pw,Pv,Temp);
%  Pv - percentage glycerol by volume
%  Pw - percentage glycerol by weight
%  Temp - temperature
%  Vis - viscosity
% from: Handbook of Chemistry and Physics (The Chemical Rubber Publishing Co., Cleveland OH, 1960, 42nd Edition) p2212

%if nargin<3; Temp=20; end;

% percentage by weight is provided
rho=1.26108;
Pv=1/(rho/Pw+1-Pw);

% if isempty(Pw);  % percentage by volume was provided
%   rho=1.26108; Pv=Pv/100;
%   Pw=rho*Pv./(rho*Pv+(1-Pv)); Pw=Pw*100;
% end;

T=[0 10 20 30 40 50 60 70 80 90 100];
bf=[...
0 1.792 1.308 1.005 0.8007 0.6560 0.5494 0.4688 0.4061 0.3565 0.3165 0.2838
10 2.44 1.74 1.31 1.03 0.826 0.680 0.575 0.500 0 0 0
20 3.44 2.41 1.76 1.35 1.07 0.879 0.731 0.635  0 0 0
30 5.14 3.49 2.50 1.87 1.46 1.16 0.956 0.816 0.690 0 0
40 8.25 5.37 3.72 2.72 2.07 1.62 1.30 1.09 0.918 0.763 0.668
50 14.6 9.01 6.00 4.21 3.10 2.37 1.86 1.53 1.25 1.05 0.910
60 29.9 17.4 10.8 7.19 5.08 3.76 2.85 2.29 1.84 1.52 1.28
65 45.7 25.3 15.2 9.85 6.80 4.89 3.66 2.91 2.28 1.86 1.55
67 55.5 29.9 17.7 11.3 7.73 5.50 4.09 3.23 2.50 2.03 1.68
70 76 38.8 22.5 14.1 9.40 6.61 4.86 3.78 2.90 2.34 1.93
75 132 65.2 35.5 21.2 13.6 9.25 6.61 5.01 3.80 3.00 2.43
80 255 116 60.1 33.9 20.8 13.6 9.42 6.94 5.13 4.03 3.18
85 540 223 109 58 33.5 21.2 14.2 10.0 7.28 5.52 4.24
90 1310 498 219 109 60.0 35.5 22.5 15.5 11.0 7.93 6.00
91 1590 592 259 127 68.1 39.8 25.1 17.1 11.9 8.62 6.40
92 1950 729 310 147 78.3 44.8 28.0 19.0 13.1 9.46 6.82
93 2400 860 367 172 89 51.5 31.6 21.2 14.4 10.3 7.54
94 2930 1040 437 202 105 58.4 35.4 23.6 15.8 11.2 8.19
95 3690 1270 523 237 121 67.0 39.9 26.4 17.5 12.4 9.08
96 4600 1580 624 281 142 77.8 45.4 29.7 19.6 13.6 10.1
97 5770 1950 765 340 166 88.9 51.9 33.6 21.9 15.1 10.9
98 7370 2460 939 409 196 104 59.8 38.5 24.8 17.0 12.2
99 9420 3090 1150 500 235 122 69.1 43.6 27.8 19.0 13.2
100 12070 3900 1410 612 284 142 81.3 50.6 31.9 21.3 14.8];

Pt=bf(:,1)'; Vt=bf(:,2:end);

Vis=interp2(Pt,T,Vt',Pw,Temp);