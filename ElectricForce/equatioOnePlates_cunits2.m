%[text] # One Charged Plate and a Point Particle
clear;
% The variables
% a : angle of point in plate
% d : distance from plate to particle
% v : relative particle velocity [0...1]
% v_a: change in velocity due to aceleration v_a = a*Dt = a*R/C
% r : radial distance between axis to point in plate
syms a d v v_a r real positive;
assume(d>0)
assume(r>=0)
assume(v>=0)
assumeAlso (v < 1)
assume(a>=0)
taylorOrder = 5;

% Particle velocity
vp = [v,0,0];
%%
%[text] ## The vectors
%[text] d= distance from the charge to the plate
%[text] r = radius of a circle in the plate
%[text] The speed of light is 1 (c=1)
%[text] 
R = (d^2+r^2)^(1/2) % total distance from a point in the plate to the charge %[output:0c987dd1]

% Get two points per circle. 
% a is the angle
vec_l1 = [d,r*cos(a),r*sin(a)]/R  % velocity vector of light from a point 1 in left positive plate at angle a %[output:7595b9fe]

vec_l2 = [d,-r*cos(a),-r*sin(a)]/R  % velocity vector of light from a point 2 in left positive plate at 180+a %[output:6ebffa30]
r_vec_l1 = (vec_l1 - vp) % left plate relative velocity of 1 particle (v in light speed fraction =v/c) %[output:269a14cd]
r_vec_l2 = (vec_l2 - vp) % left plate relative velocity of 2 particle (v in light speed fraction =v/c) %[output:273c8f2d]

%%
%[text] ## The Anomality 
% 1 particle disparity between photon angle and trayectory


%psin = simplify(sqrt(1 - dot(vp,vp))*sqrt(dot(r_vec_l1,r_vec_l1)));
psin = simplify(sqrt(1 - dot(vp,vp)));


%%
%[text] ## The Relative Net Force of the Two Points
% Force of the two point in plate and moving particle
dforce = (psin*vec_l1 + psin*vec_l2 )/(2*R^2); % net Force in the x direction first points.

dforce = simplify(expand(dforce),'Steps',100) % Two charges %[output:16ef0998]
dforce_x = dforce(1) %[output:1862f7b5]
%dforce_x = simplify(taylor(expand(subs(dforce_x/2,v,v+v_a) + dforce_x/2),v_a,order=2))
%%
%[text] ## Check at v = 0, r =0
df_v_zero = simplify(subs(dforce_x,[v,r],[0,0]),'Steps',100) % Due to the two point charges %[output:3d0726da]
df_v_one = simplify(subs(dforce_x,[v,r],[v,0]),'Steps',100) % Due to the two point charges %[output:69e50e72]
%%
%[text] ## Integrate all the angles
%dforce_xT = expand(taylor(dforce_x,v,Order=taylorOrder));
dforce_xT = dforce_x;

dcircle = int(dforce_xT,a,0,pi); % From zero to pi i.e., the full circle
%%
dcircle = simplify(dcircle,steps=100); % Net force of the two full circles with charge density
df_v_zero = simplify(subs(dcircle,[v,r],[v,0]),'Steps',100) % Net force of the full zero radius circle (point) %[output:4bf04a0c]
%[text] ## 
%%
%[text] ## Integrate for all values of r
syms R_o real
assumeAlso (R_o>0)
dtotf = r*dcircle/pi  %[output:7b7123d6]
dtotf = simplify(expand(taylor(dtotf,v,Order=taylorOrder))); % Taylor
%dtotf = simplify(subs(dtotf,v,0)) %static
totInt = int(dtotf,r,0,R_o); % Add all the circles
%%
totf = simplify(totInt,'Criterion','preferReal','Steps',100) % The force on a charge due a charge surface density %[output:3b43bbb2]

%%
%[text] ## Set the Plates Radius to Infinity
sympref('PolynomialDisplayStyle', 'ascend');
totfInf = simplify(limit(totf,R_o,Inf)) % Net force by plate charge density %[output:34b522cd]
%totfInf = subs(totfInf,v,-v)

%%
figure(1)
plotrange = [0,0.99];
fplot(totfInf,plotrange) %[output:81e51b41]
hold on %[output:81e51b41]
relEq = sqrt(1-v^2) %[output:4739c8cb]
relEqApx = taylor(relEq,v,Order=taylorOrder) %[output:38ff51e6]
fplot(relEq,plotrange) %[output:81e51b41]
fplot(relEqApx,plotrange) %[output:81e51b41]
xlabel('v/c')  %[output:81e51b41]
ylabel('ratio')  %[output:81e51b41]
legend({'Force Carriers',"Lorentz Mass","Lorentz Apx"},'Location','northeast') %[output:81e51b41]
hold off %[output:81e51b41]



%[text] ## Set the Plates Radius to 0.1\*d
totf_d = simplify(limit(totf,R_o,0.1*d)) % Net force by plate charge density %[output:42cc0d02]
totf_dV0 = subs(totf_d,v=0) %[output:931e7217]

totf_dRel = simplify(totf_d/totf_dV0) %[output:1f52ffa5]
totf_dRel0 = subs(totf_dRel,v_a,0) %[output:75135088]

%%
%[text] ## Plots
figure(1)
plotrange = [0,0.99];
fplot(totfInf,plotrange) %[output:91191ae5]
hold on %[output:91191ae5]
relEq = sqrt(1-v^2) %[output:7166c1d5]
relEqApx = taylor(relEq,v,Order=taylorOrder) %[output:5d71f32e]
fplot(relEq,plotrange) %[output:91191ae5]
fplot(relEqApx,plotrange) %[output:91191ae5]
fplot(totf_dRel0,plotrange) %[output:91191ae5]
xlabel('v/c')  %[output:91191ae5]
ylabel('ratio')  %[output:91191ae5]
legend({'Force Carriers',"Lorentz Mass","Lorentz Apx","Fd"},'Location','northeast') %[output:91191ae5]
hold off %[output:91191ae5]
%[text] ## 
%%
%[text] ## 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":22}
%---
%[output:0c987dd1]
%   data: {"dataType":"symbolic","outputData":{"name":"R","value":"\\sqrt{r^2 +d^2 }"}}
%---
%[output:7595b9fe]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_l1","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:6ebffa30]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_l2","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:269a14cd]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l1","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }}-v & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:273c8f2d]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l2","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }}-v & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:16ef0998]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce","value":"\\left(\\begin{array}{ccc}\n\\frac{d\\,\\sqrt{1-v^2 }}{{{\\left(r^2 +d^2 \\right)}}^{3\/2} } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:1862f7b5]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce_x","value":"\\frac{d\\,\\sqrt{1-v^2 }}{{{\\left(r^2 +d^2 \\right)}}^{3\/2} }"}}
%---
%[output:3d0726da]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\frac{1}{d^2 }"}}
%---
%[output:69e50e72]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_one","value":"\\frac{\\sqrt{1-v^2 }}{d^2 }"}}
%---
%[output:4bf04a0c]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\frac{\\pi \\,\\sqrt{1-v^2 }}{d^2 }"}}
%---
%[output:7b7123d6]
%   data: {"dataType":"symbolic","outputData":{"name":"dtotf","value":"\\frac{d\\,r\\,\\sqrt{1-v^2 }}{{{\\left(r^2 +d^2 \\right)}}^{3\/2} }"}}
%---
%[output:3b43bbb2]
%   data: {"dataType":"symbolic","outputData":{"name":"totf","value":"\\frac{{\\left(\\frac{d}{\\sqrt{{R_o }^2 +d^2 }}-1\\right)}\\,{\\left(-8+4\\,v^2 +v^4 \\right)}}{8}"}}
%---
%[output:34b522cd]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInf","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}"}}
%---
%[output:4739c8cb]
%   data: {"dataType":"symbolic","outputData":{"name":"relEq","value":"\\sqrt{1-v^2 }"}}
%---
%[output:38ff51e6]
%   data: {"dataType":"symbolic","outputData":{"name":"relEqApx","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}"}}
%---
%[output:81e51b41]
%   data: {"dataType":"image","outputData":{"dataUri":"data:,","height":0,"width":0}}
%---
%[output:42cc0d02]
%   data: {"dataType":"symbolic","outputData":{"name":"totf_d","value":"{\\left(\\frac{10\\,\\sqrt{101}}{101}-1\\right)}\\,{\\left(-1+\\frac{v^2 }{2}+\\frac{v^4 }{8}\\right)}"}}
%---
%[output:931e7217]
%   data: {"dataType":"symbolic","outputData":{"name":"totf_dV0","value":"1-\\frac{10\\,\\sqrt{101}}{101}"}}
%---
%[output:1f52ffa5]
%   data: {"dataType":"symbolic","outputData":{"name":"totf_dRel","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}"}}
%---
%[output:75135088]
%   data: {"dataType":"symbolic","outputData":{"name":"totf_dRel0","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}"}}
%---
%[output:7166c1d5]
%   data: {"dataType":"symbolic","outputData":{"name":"relEq","value":"\\sqrt{1-v^2 }"}}
%---
%[output:5d71f32e]
%   data: {"dataType":"symbolic","outputData":{"name":"relEqApx","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}"}}
%---
%[output:91191ae5]
%   data: {"dataType":"image","outputData":{"dataUri":"data:,","height":0,"width":0}}
%---
