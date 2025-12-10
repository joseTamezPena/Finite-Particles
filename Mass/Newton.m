% Constants
electron_mass = 9.10938356e-31; % kg
electron_charge = 1.602176634e-19; % C
speed_of_light = 299792458; % m/s
electric_permeability = 8.854187817e-12; % F/m
electron_radius = 1.0e-20 %[output:40354b7d]

% Display the constants
fprintf('Electron Mass: %.5e kg\n', electron_mass); %[output:8b80d836]
fprintf('Electron Charge: %.5e C\n', electron_charge); %[output:042ae10c]
fprintf('Speed of Light: %.5e m/s\n', speed_of_light); %[output:91863228]
fprintf('Electric Permeability: %.5e F/m\n', electric_permeability); %[output:2373cc35]
%%
% Classical electron radius
classical_electron_radius = (electron_charge^2) / (4 * pi * electric_permeability * electron_mass * speed_of_light^2); % m

% Display the classical electron radius
fprintf('Classical Electron Radius: %.5e m\n', classical_electron_radius); %[output:2a1248ff]
%%
%[text] ## The Newton force based on finite particle model
% External field accting on a finite particle
% Abs_int = Abs_ext

% q_e emitted photons per second
% m_o absorbed photons per second
% c speed of light
% r_o particle size(radius)
% L_o separation between two particle elements
% a particle aceleration
% Dt crossing time of light 

syms q_e m_o c r_o L_o a
v = a*L_o/c %[output:7b314032]
Abs_int_l = q_e*m_o*(c+v)*r_o^3/(16*pi*(L_o*c)^2) %[output:5a7bc871]
Abs_int_r = q_e*m_o*(c-v)*r_o^3/(16*pi*(L_o*c)^2) %[output:1f4fcd8e]
Abs_int = simplify(Abs_int_l-Abs_int_r) %[output:78206d3d]

% Newton force is the speed of light times the interal absorved
Force = c*Abs_int %[output:9338cc03]
% The particle mass is:
mass = Force/a %[output:00977042]
% The particle energy is proportional to q^2/L_o. Hence E=mc^2
%%
%[text] ## The number of absorbed elements per second
%particle mass

mass %[output:283dd844]
syms tot_ets; % Total events We joint the emited:q_e and absorbed m_o
mass = subs(mass,q_e*m_o,tot_ets) %[output:7202a2c5]
mass = subs(mass,r_o,electron_radius) % Use the estimated (observed) electron radius %[output:1e9d2580]
mass = subs(mass,L_o,0.5*classical_electron_radius) % L_o is assumed half the classical electron radius %[output:1bf0037d]
mass = subs(mass,c,speed_of_light) % Subs the speed of light %[output:68076807]
% Solve for the number of interactions per second
tot_p_second = solve(mass==electron_mass) %[output:13fac332]
vpa(tot_p_second) %[output:0545720b]
fprintf('Estimated number of absorbed photons per second: %.5e\n', vpa(tot_p_second)); %[output:7258661e]

%%
%[text] ## Interactions based on Coulomb's law
% Elements per second based on electron charge

ColumbsF = (electron_charge^2)/(4*pi*electric_permeability) %[output:38f85a59]
abscharge = tot_ets*r_o^3/(4*pi) % Based total absorbed between two electrons %[output:9e8c9127]
abscharge = subs(abscharge,r_o,electron_radius) % The electron size %[output:596c1116]

Qtot_p_second = solve(ColumbsF==abscharge) %[output:1e21467c]
fprintf('Estimated number of absorbed photons per second: %.5e \n', vpa(Qtot_p_second)) %[output:1833f1e8]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:40354b7d]
%   data: {"dataType":"textualVariable","outputData":{"name":"electron_radius","value":"1.0000e-20"}}
%---
%[output:8b80d836]
%   data: {"dataType":"text","outputData":{"text":"Electron Mass: 9.10938e-31 kg\n","truncated":false}}
%---
%[output:042ae10c]
%   data: {"dataType":"text","outputData":{"text":"Electron Charge: 1.60218e-19 C\n","truncated":false}}
%---
%[output:91863228]
%   data: {"dataType":"text","outputData":{"text":"Speed of Light: 2.99792e+08 m\/s\n","truncated":false}}
%---
%[output:2373cc35]
%   data: {"dataType":"text","outputData":{"text":"Electric Permeability: 8.85419e-12 F\/m\n","truncated":false}}
%---
%[output:2a1248ff]
%   data: {"dataType":"text","outputData":{"text":"Classical Electron Radius: 2.81794e-15 m\n","truncated":false}}
%---
%[output:7b314032]
%   data: {"dataType":"symbolic","outputData":{"name":"v","value":"\\frac{L_o \\,a}{c}"}}
%---
%[output:5a7bc871]
%   data: {"dataType":"symbolic","outputData":{"name":"Abs_int_l","value":"\\frac{m_o \\,q_e \\,{r_o }^3 \\,{\\left(c+\\frac{L_o \\,a}{c}\\right)}}{16\\,{L_o }^2 \\,c^2 \\,\\pi }"}}
%---
%[output:1f4fcd8e]
%   data: {"dataType":"symbolic","outputData":{"name":"Abs_int_r","value":"\\frac{m_o \\,q_e \\,{r_o }^3 \\,{\\left(c-\\frac{L_o \\,a}{c}\\right)}}{16\\,{L_o }^2 \\,c^2 \\,\\pi }"}}
%---
%[output:78206d3d]
%   data: {"dataType":"symbolic","outputData":{"name":"Abs_int","value":"\\frac{a\\,m_o \\,q_e \\,{r_o }^3 }{8\\,L_o \\,c^3 \\,\\pi }"}}
%---
%[output:9338cc03]
%   data: {"dataType":"symbolic","outputData":{"name":"Force","value":"\\frac{a\\,m_o \\,q_e \\,{r_o }^3 }{8\\,L_o \\,c^2 \\,\\pi }"}}
%---
%[output:00977042]
%   data: {"dataType":"symbolic","outputData":{"name":"mass","value":"\\frac{m_o \\,q_e \\,{r_o }^3 }{8\\,L_o \\,c^2 \\,\\pi }"}}
%---
%[output:283dd844]
%   data: {"dataType":"symbolic","outputData":{"name":"mass","value":"\\frac{m_o \\,q_e \\,{r_o }^3 }{8\\,L_o \\,c^2 \\,\\pi }"}}
%---
%[output:7202a2c5]
%   data: {"dataType":"symbolic","outputData":{"name":"mass","value":"\\frac{{r_o }^3 \\,{\\textrm{tot}}_{\\textrm{ets}} }{8\\,L_o \\,c^2 \\,\\pi }"}}
%---
%[output:1e9d2580]
%   data: {"dataType":"symbolic","outputData":{"name":"mass","value":"\\frac{{\\textrm{tot}}_{\\textrm{ets}} }{8000000000000000000000000000000000000000000000000000000000000\\,L_o \\,c^2 \\,\\pi }"}}
%---
%[output:1bf0037d]
%   data: {"dataType":"symbolic","outputData":{"name":"mass","value":"\\frac{549755813888\\,{\\textrm{tot}}_{\\textrm{ets}} }{6196716403399281543340482159010207396931946277618408203125\\,c^2 \\,\\pi }"}}
%---
%[output:68076807]
%   data: {"dataType":"symbolic","outputData":{"name":"mass","value":"\\frac{137438953472\\,{\\textrm{tot}}_{\\textrm{ets}} }{139233273967962276116475909047168375842762344518632744438946247100830078125\\,\\pi }"}}
%---
%[output:13fac332]
%   data: {"dataType":"symbolic","outputData":{"name":"tot_p_second","value":"\\frac{724088025030325759890700151428021481751867582552124880113098015499417670071125030517578125\\,\\pi }{784637716923335095479473677900958302012794430558004314112}"}}
%---
%[output:0545720b]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\text{2.8991591545042659355792865333043e+33}"}}
%---
%[output:7258661e]
%   data: {"dataType":"text","outputData":{"text":"Estimated number of absorbed photons per second: 2.89916e+33\n","truncated":false}}
%---
%[output:38f85a59]
%   data: {"dataType":"textualVariable","outputData":{"name":"ColumbsF","value":"2.3071e-28"}}
%---
%[output:9e8c9127]
%   data: {"dataType":"symbolic","outputData":{"name":"abscharge","value":"\\frac{{r_o }^3 \\,{\\textrm{tot}}_{\\textrm{ets}} }{4\\,\\pi }"}}
%---
%[output:596c1116]
%   data: {"dataType":"symbolic","outputData":{"name":"abscharge","value":"\\frac{{\\textrm{tot}}_{\\textrm{ets}} }{4000000000000000000000000000000000000000000000000000000000000\\,\\pi }"}}
%---
%[output:1e21467c]
%   data: {"dataType":"symbolic","outputData":{"name":"Qtot_p_second","value":"\\frac{1115634247822495060675773714820024906657636165618896484375\\,\\pi }{1208925819614629174706176}"}}
%---
%[output:1833f1e8]
%   data: {"dataType":"text","outputData":{"text":"Estimated number of absorbed photons per second: 2.89916e+33 \n","truncated":false}}
%---
