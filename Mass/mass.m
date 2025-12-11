%[text] # Constants

G=6.67430e-11; % Gravitational constant 
Q=1.602176634e-19; % Elemetrary charge
r=0.877e-15; % Proton Radius
eo = 8.8541878188e-12; %Vacuum permittivity 
c=physconst('LightSpeed');
alpha = 0.0072973525643 % Fine-structure constant %[output:5519051d]


%%
%[text] ## Speed variance and beta correction
% Estimation of the electron velocity  based on the fine-struture constant
SeRel = alpha % Standard deviation of Reletative Electron velocity  %[output:2120f4a9]
abscoef = 2*G*Q/(4*pi()*eo*(c*r)^3)/(SeRel^2) %[output:2ac141a4]
b= abscoef*Q*r/(2*c) %[output:5a1b6c18]
%%
%[text] ## Force between atoms

syms c v m Q r real positive;
assume(v>0)
assume(m>0)
assume(Q>0)
assume(r>0)
assume(c>0)

cx= (c^2-v^2)^(1/2) %[output:1324b646]
beta = 1/2*m*Q*r/cx %[output:3599f470]
Fnet = simplify (beta*((beta-1)+(cx/c)^3),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:1658c6ac]
Ftay = simplify (taylor(Fnet,v,ExpansionPoint=0,Order=4)) %[output:1cef324b]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":33.6}
%---
%[output:5519051d]
%   data: {"dataType":"textualVariable","outputData":{"name":"alpha","value":"0.0073"}}
%---
%[output:2120f4a9]
%   data: {"dataType":"textualVariable","outputData":{"name":"SeRel","value":"0.0073"}}
%---
%[output:2ac141a4]
%   data: {"dataType":"textualVariable","outputData":{"name":"abscoef","value":"1.9861e+05"}}
%---
%[output:5a1b6c18]
%   data: {"dataType":"textualVariable","outputData":{"name":"b","value":"4.6543e-38"}}
%---
%[output:1324b646]
%   data: {"dataType":"symbolic","outputData":{"name":"cx","value":"\\sqrt{c^2 -v^2 }"}}
%---
%[output:3599f470]
%   data: {"dataType":"symbolic","outputData":{"name":"beta","value":"\\frac{Q\\,m\\,r}{2\\,\\sqrt{c^2 -v^2 }}"}}
%---
%[output:1658c6ac]
%   data: {"dataType":"symbolic","outputData":{"name":"Fnet","value":"\\frac{Q^2 \\,m^2 \\,r^2 }{4\\,{\\left(c^2 -v^2 \\right)}}-\\frac{Q\\,m\\,r}{2\\,\\sqrt{c^2 -v^2 }}+\\frac{Q\\,m\\,r\\,{\\left(c^2 -v^2 \\right)}}{2\\,c^3 }"}}
%---
%[output:1cef324b]
%   data: {"dataType":"symbolic","outputData":{"name":"Ftay","value":"\\frac{Q\\,m\\,r\\,{\\left(Q\\,c^2 \\,m\\,r+Q\\,m\\,r\\,v^2 -3\\,c\\,v^2 \\right)}}{4\\,c^4 }"}}
%---
