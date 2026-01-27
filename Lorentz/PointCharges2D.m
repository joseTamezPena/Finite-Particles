%[text] # Electric and magnetic field due to a charge moving at constant velocity
clear;
syms  q  real 
syms t epsilon0 mu0 real positive
syms x y c v_x v_y c_x c_y real positive
assumeAlso (v_x^2 + v_y^2 < c^2)
assumeAlso (c^2 - v_x^2 > 0)
assumeAlso (c^2 - v_y^2 > 0)

syms vp_x vp_y real positive

assumeAlso (vp_x^2+vp_y^2 < c^2)
assumeAlso (c^2 - vp_x^2 > 0)
assumeAlso (c^2 - vp_y^2 > 0)
assumeAlso (c^2*x - v_y^2*x - v_x*sqrt(c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2) + v_x*v_y*y >0 )
assumeAlso (c^2*y - v_x^2*y - v_y*sqrt(c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2) + v_x*v_y*x >0 )


% The emission charge is at zero at t=0
r_o = [0,0,0];
% Define the position of the space where the field will be measured
r_p = [x, y, 0];

% Define the velocity vector of the moving charge
v = [v_x, v_y, 0];

% Light vector
cv = [c_x, c_y, 0];

% The velocity of the reciever charged particle at r_p

vp = [vp_x,vp_y,0];


%%
%[text] ## Compute the time to travel

% the two reach the same point
eq_1=(v+cv)*t==r_p %[output:12620cf2]
eq_2 = dot(cv,cv)==c^2 %[output:223fbfc8]


[ccx,ccy,tt] = solve([eq_1,eq_2],[c_x,c_y,t],"Real",true,"IgnoreAnalyticConstraints",true); %[output:3ad81ded]
ccx = simplify(ccx(2));
ccy = simplify(ccy(2));
cvec = [-ccx,ccy,0];
assumeAlso (c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2 > 0)
simplify(dot(cvec,cvec)) %[output:7e43424c]

tt = simplify(tt(1));
simplify(subs(ccx,y,0)), simplify(subs(ccy,y,0)), simplify(subs(tt,y,0)) %[output:5489c3d5] %[output:19fe85d4] %[output:56a589a1]
%%

r_c = simplify(expand(v*tt)) %[output:560d68f8]

simplify(subs(r_c,y,0)) %[output:8c55f7d6]
simplify(subs(r_c,[y,v_x],[0,0])) %[output:62a4a62b]

%%
%[text] ## Electric and Magnetic Field
assumeAlso (v_x*x - sqrt(c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2) + v_y*y>0)


% Electric field E due to a moving charge
pvec = r_p - r_c;
simplify(subs(pvec,[y,v_x],[0,0])) %[output:2b2ea048]

E = (1/(4*pi*epsilon0)) * (q * pvec) / norm(pvec)^3;
E = simplify(subs(E,y=0));
% Magnetic field B due to a moving charge
B = simplify(subs((1.0/(c^2)) * cross(v,E),y,0));

% Display the equations
simplify(subs(E,[y,v_x],[0,0])), simplify(subs(B,[y,v_x],[0,0])) %[output:95efb827] %[output:721998ca]
%%
%[text] ## The Lorenz force



% Define the Lorentz force F acting on a second charge with vp velocity
F = simplify(q * (E + cross(vp, B)),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

Fx = simplify(subs(F,[y,v_x],[0,0]));


Ft = simplify(taylor(Fx,v_y,0,order=2),'Criterion','preferReal','Steps',100) %[output:26331b5d]

F_static = subs(Ft,v_y,0) %[output:798825f9]
Ft_neutral = simplify((Ft - F_static),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:1675db2f]
ft_magnetic = simplify(q * cross(vp, B),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = taylor(ft_magnetic,[v_x,v_y],[0,0],order=2);
ft_magnetic = simplify(taylor(ft_magnetic,[vp_x,vp_y],[0,0],order=2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = simplify(subs(ft_magnetic,[v_x,y],[0,0])) %[output:346cae46]
%[text] ## The corpuscle force
%[text] 
%[text] The corpuscle force is:
%[text] $&dollar&;&dollar&; \n\\mathbf{F} = \\frac{k q\_1 q\_2}{4 \\pi r^2} \\left( \\frac{\\|\\mathbf{c} + \\mathbf{v\_1}- \\mathbf{v\_2}\\|^2}{\\mathbf{c} \\cdot ( \\mathbf{c} + \\mathbf{v\_1} - \\mathbf{v\_2} )} \\right)  \\hat{o\_1}, \n&dollar&;&dollar&;$
%[text] 
%%
k= q^2/(4*pi*epsilon0);
assumeAlso(c^2 - v_y^2>0)
assumeAlso(c^2 - v_x^2>0)
assumeAlso(c - v_x>0)
assumeAlso(c - v_y>0)
assumeAlso(2*c - v_x>0)
assumeAlso(2*c - v_y>0)
assumeAlso(c - vp_x>0)
assumeAlso(c - vp_y>0)
assumeAlso(c^2*y^2 - v_x^2*y^2 + c^2*x^2 >0)
assumeAlso(c^2 - vp_x^2 - vp_y^2 >0)
assumeAlso(c^2 - vp_y^2 + 2*c*vp_x - vp_x^2 > 0)

R2 = dot(r_p,r_p);

ov = simplify(cvec/sqrt(dot(cvec,cvec)));

vs = cvec+v;
rvel = vs-vp;
vprel = v-vp;

% Dot products
rvelm2 = dot(rvel,rvel);
vrmag2 = dot(vprel,vprel);
rvels = dot(vs,vs);


% Assume at y=0 and v_x=0

s_sta2 = simplify(subs(rvels,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %Speed^2 
c_sta2 = simplify(subs(rvelm2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %Speed^2 
cmag = sqrt(c_sta2); %Speed 
smag = sqrt(s_sta2); %Speed 
dpmp_n = simplify(subs(dot(ov,rvel/cmag),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %raw product
dpmp_sr = simplify(subs(dot(cvec,vs),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %raw product
dpmp_sn = simplify(subs(dot(ov,vs/smag),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %raw product
dpmp_nn = simplify(subs(dot(rvel/cmag,vs/smag),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %raw product

%% Distance Factors
p1 = simplify(subs(k/R2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % Density at propagation surface
p2 = simplify(s_sta2/c^2,'Criterion','preferReal','Steps',100); % velocity correction for distance at emission

%% Action Direction in velocity vector

ov0 = simplify(subs(ov,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % The action orientation

%% Action Omegas 

% === Works for magnetic force===
%p3 = simplify(subs((cmag/c)/dpmp_n,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % density and velocity correction
p3 = simplify(subs((cmag*c)/(s_sta2*dpmp_n),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % density and velocity correction
st_force = simplify(subs(p1*p2*p3*ov0,v_y,0),'Criterion','preferReal','Steps',100);
%===== 

% === Works Imitate Lorenz force===
%vr_sta2 = simplify(subs(dot(vprel,vprel)/c^2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %Speed
%fact = simplify(subs(sqrt(dot(v/c,ov)^2),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100);
%p3 = simplify(subs(sqrt(1-vr_sta2) + fact,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % density and velocity correction
%p4 = simplify(subs((dot(vp,ov0)*v - 2*dot(vp,v)*ov0)/c^2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100);
%st_force = simplify(subs(p1*p2*(p3*ov0+p4),v_y,0),'Criterion','preferReal','Steps',100)
%===== 



st_forceT = simplify(taylor(st_force,vp_y,0,order=2));
st_forceT = simplify(taylor(st_forceT,vp_x,0,order=2)) %[output:9aff099c]

F2_c = simplify(p1*p2*p3*ov0,'Criterion','preferReal','Steps',100);
F2T_c = simplifyFraction(simplify(taylor(F2_c,v_y,0,order=2),'Criterion','preferReal','Steps',100),'Expand',true);

% The net Force due to corpuscles
F2T_net = simplify(taylor(F2T_c,vp_y,0,order=2));
F2T_net = simplify(taylor(F2T_net,vp_x,0,order=2),'Criterion','preferReal','Steps',100),Ft %[output:0bf16652] %[output:8c554a06]
%F2T_net = simplify(taylor(F2T_c,[vp_x,vp_y],[0,0],order=2))
% The net magnetic Force due to corpuscles
F2T_m=simplify(expand(F2T_net-st_forceT),'Criterion','preferReal','Steps',100) %[output:60146e59]
F2_m=simplify(expand(F2_c-st_force),'Criterion','preferReal','Steps',100);

F2T_mnet = simplify(taylor(F2_m,v_y,0,order=2),'Criterion','preferReal','Steps',100);
F2T_mnet = simplify(taylor(F2T_mnet,vp_x,0,order=2),'Criterion','preferReal','Steps',100);
F2T_mnet = simplify(taylor(F2T_mnet,vp_y,0,order=2),'Criterion','preferReal','Steps',100) %[output:33432fd0]

%%
% The Lorentz total force
ft_magnetic %[output:912e680e]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:12620cf2]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_1","value":"\\left(\\begin{array}{ccc}\nt\\,{\\left(c_x +v_x \\right)}=x & t\\,{\\left(c_y +v_y \\right)}=y & 0=0\n\\end{array}\\right)"}}
%---
%[output:223fbfc8]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_2","value":"{c_y }^2 +{c_x }^2 =c^2"}}
%---
%[output:3ad81ded]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:7e43424c]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"c^2"}}
%---
%[output:5489c3d5]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-\\sqrt{c^2 -{v_y }^2 }"}}
%---
%[output:19fe85d4]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-v_y"}}
%---
%[output:56a589a1]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\frac{x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 }"}}
%---
%[output:560d68f8]
%   data: {"dataType":"symbolic","outputData":{"name":"r_c","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{v_x \\,\\sigma_1 }{-c^2 +{v_y }^2 +{v_x }^2 } & \\frac{v_y \\,\\sigma_1 }{-c^2 +{v_y }^2 +{v_x }^2 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =v_x \\,x-\\sqrt{c^2 \\,y^2 -{v_x }^2 \\,y^2 +2\\,v_x \\,v_y \\,x\\,y+c^2 \\,x^2 -{v_y }^2 \\,x^2 }+v_y \\,y\n\\end{array}"}}
%---
%[output:8c55f7d6]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{v_x \\,x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 } & \\frac{v_y \\,x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 } & 0\n\\end{array}\\right)"}}
%---
%[output:62a4a62b]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & \\frac{v_y \\,x}{\\sqrt{c^2 -{v_y }^2 }} & 0\n\\end{array}\\right)"}}
%---
%[output:2b2ea048]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\nx & -\\frac{v_y \\,x}{\\sqrt{c^2 -{v_y }^2 }} & 0\n\\end{array}\\right)"}}
%---
%[output:95efb827]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{q\\,{{\\left(c^2 -{v_y }^2 \\right)}}^{3\/2} }{4\\,c^3 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q\\,v_y \\,{\\left(c^2 -{v_y }^2 \\right)}}{4\\,c^3 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:721998ca]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & 0 & -\\frac{q\\,v_y \\,{{\\left(c^2 -{v_y }^2 \\right)}}^{3\/2} }{4\\,c^5 \\,\\varepsilon_0 \\,x^2 \\,\\pi }\n\\end{array}\\right)"}}
%---
%[output:26331b5d]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:798825f9]
%   data: {"dataType":"symbolic","outputData":{"name":"F_static","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:1675db2f]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_neutral","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:346cae46]
%   data: {"dataType":"symbolic","outputData":{"name":"ft_magnetic","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:9aff099c]
%   data: {"dataType":"symbolic","outputData":{"name":"st_forceT","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:0bf16652]
%   data: {"dataType":"symbolic","outputData":{"name":"F2T_net","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y +c\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:8c554a06]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:60146e59]
%   data: {"dataType":"symbolic","outputData":{"name":"F2T_m","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:33432fd0]
%   data: {"dataType":"symbolic","outputData":{"name":"F2T_mnet","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:912e680e]
%   data: {"dataType":"symbolic","outputData":{"name":"ft_magnetic","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
