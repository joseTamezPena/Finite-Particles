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

Ft = taylor(Fx,v_y,0,order=2) %[output:3d8c6337]

%Ft = Fx;
Ft = taylor(Ft,[vp_x,vp_y],[0,0],order=2);

Ft = simplify(Ft,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:6718b8ac]
F_static = subs(Ft,v_y,0) %[output:4e0b658b]
Ft_neutral = simplify((Ft - F_static),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:1058ae39]
ft_magnetic = simplify(q * cross(vp, B),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = taylor(ft_magnetic,[v_x,v_y],[0,0],order=2);
ft_magnetic = simplify(taylor(ft_magnetic,[vp_x,vp_y],[0,0],order=2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = simplify(subs(ft_magnetic,[v_x,y],[0,0])) %[output:922d6e6a]
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
pvec2 = dot(pvec,pvec);
rvels = dot(vs,vs);



% Assume at y=0 and v_x=0
cprod = simplify(subs(cross(cvec,rvel),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100);
mcprod = simplify(dot(cprod,cprod)/c^2,'Criterion','preferReal','Steps',100);

vr_sta2 = simplify(subs(dot(vprel,vprel)/c^2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %Speed
s_sta2 = simplify(subs(rvels,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %Speed^2 
c_sta2 = simplify(subs(rvelm2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %Speed^2 
cmag = simplify(subs(sqrt(rvelm2),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %Speed 
vsmag = simplify(subs(sqrt(rvels),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %Speed 
dpmp_r = simplify(subs(dot(cvec,rvel),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %raw product
dpmp_n = simplify(subs(dot(ov,rvel/cmag),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %Cosine
dotrel_c = simplify(subs(dot(vprel/c,ov),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %Cosine;

%% Action Factors
p1 = simplify(subs(k/R2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % distance factor
p2 = simplify(s_sta2/c^2,'Criterion','preferReal','Steps',100); % density and velocity correction
%% Action Direction in velocity vector
ov0 = simplify(subs(ov,[v_x,y],[0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100); % The action orientation
%% Action direction in relative velocity
p3 = simplify(sqrt(1-vr_sta2),'Criterion','preferReal','Steps',100); % density and velocity correction
%p4 = subs(dotrel_c*vprel/c,[v_x,y],[0,0]);
p4 = simplify(subs((dot(vp,ov0)*v - 2*dot(vp,v)*ov0)/c^2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100);

st_force = simplify(subs(p1*p2*p3*ov0,[v_x,v_y,y],[0,0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
st_force = simplify(taylor(st_force,vp_y,0,order=2));
st_force = simplify(taylor(st_force,vp_x,0,order=2)) %[output:9fb2f05d]

F2_c = simplify(p1*p2*p3*(ov0+p4),'Criterion','preferReal','Steps',100);
F2T_c = simplify(taylor(F2_c,v_y,0,order=2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

% The net Force due to corpuscles
F2T_net = simplify(taylor(F2T_c,vp_x,0,order=2));
F2T_net = simplify(taylor(F2T_net,vp_y,0,order=2),'Criterion','preferReal','Steps',100),Ft %[output:10b606d1] %[output:23b0c75e]
%F2T_net = simplify(taylor(F2T_c,[vp_x,vp_y],[0,0],order=2))
% The net magnetic Force due to corpuscles
F2T_m=simplify(expand(F2T_net-st_force),'Criterion','preferReal','Steps',100) %[output:5a895bcf]

%%
dp1 = subs(dot(vp,ov0)*v - 2*dot(vp,v)*ov0,v_x,0);
dp1 = simplify(taylor(dp1,v_y,0,order=2));
dp1 = simplify(taylor(dp1,vp_x,0,order=2));
simplify(taylor(dp1,vp_y,0,order=2)) %[output:0c3b538c]
dp1 = subs(dot(vprel,ov0)*c*ov0,[v_x,y],[0,0]);
dp1 = simplify(taylor(dp1,v_y,0,order=2));
dp1 = simplify(taylor(dp1,vp_x,0,order=2));
simplify(taylor(dp1,vp_y,0,order=2)) %[output:6af7a300]
dp1 = subs(dot(vprel,ov0)*vprel,[v_x,y],[0,0]);
dp1 = simplify(taylor(dp1,v_y,0,order=2));
dp1 = simplify(taylor(dp1,vp_x,0,order=2));
simplify(taylor(dp1,vp_y,0,order=2)) %[output:987668a1]
dp1 = subs(dot(vprel,ov0)*(vprel-c*ov0),[v_x,y],[0,0]);
dp1 = simplify(taylor(dp1,v_y,0,order=2));
dp1 = simplify(taylor(dp1,vp_x,0,order=2));
simplify(taylor(dp1,vp_y,0,order=2)) %[output:7fd3ad77]
%%
% The Lorentz total force
ft_magnetic %[output:1ed65c91]

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
%[output:3d8c6337]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi }-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(-c^3 +c^2 \\,{\\textrm{vp}}_x \\right)}}{4\\,c^4 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:6718b8ac]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:4e0b658b]
%   data: {"dataType":"symbolic","outputData":{"name":"F_static","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:1058ae39]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_neutral","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:922d6e6a]
%   data: {"dataType":"symbolic","outputData":{"name":"ft_magnetic","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:9fb2f05d]
%   data: {"dataType":"symbolic","outputData":{"name":"st_force","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:10b606d1]
%   data: {"dataType":"symbolic","outputData":{"name":"F2T_net","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:23b0c75e]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:5a895bcf]
%   data: {"dataType":"symbolic","outputData":{"name":"F2T_m","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:0c3b538c]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-2\\,v_y \\,{\\textrm{vp}}_y  & v_y \\,{\\textrm{vp}}_x  & 0\n\\end{array}\\right)"}}
%---
%[output:6af7a300]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\nv_y \\,{\\textrm{vp}}_y -c\\,{\\textrm{vp}}_x  & v_y \\,{\\textrm{vp}}_x  & 0\n\\end{array}\\right)"}}
%---
%[output:987668a1]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{v_y \\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y }{c} & -{\\textrm{vp}}_x \\,{\\left(v_y -{\\textrm{vp}}_y \\right)} & 0\n\\end{array}\\right)"}}
%---
%[output:7fd3ad77]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\nc\\,{\\textrm{vp}}_x -{\\textrm{vp}}_y \\,{\\left(v_y +\\frac{v_y \\,{\\textrm{vp}}_x }{c}\\right)} & -{\\textrm{vp}}_x \\,{\\left(2\\,v_y -{\\textrm{vp}}_y \\right)} & 0\n\\end{array}\\right)"}}
%---
%[output:1ed65c91]
%   data: {"dataType":"symbolic","outputData":{"name":"ft_magnetic","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
