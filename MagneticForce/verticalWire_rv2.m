%[text] # A single vertical wire with current
clear;

syms  q  real 
syms t c epsilon0 mu0 real positive
syms x y v_y c_x real 
syms c_y real 
syms vp_x vp_y real 


assumeAlso (t > 0)
assumeAlso (x > 0)
assumeAlso (y > 0)
assumeAlso (c > 0)
assumeAlso (v_y > 0)
assumeAlso (c_x > 0)
assumeAlso (c_y < 0)
assumeAlso (vp_x > 0)
assumeAlso (vp_y > 0)
assumeAlso (c - vp_x > 0)
assumeAlso (c - vp_y > 0)
assumeAlso (c^2 - v_y^2 > 0)
assumeAlso (c^2 - vp_y^2 - vp_y^2 > 0)
assumeAlso (c^2*y^2 + c^2*x^2 - v_y^2*x^2 > 0)


% Define the position of the space where the force will be measured
r_p = [x, 0, 0];

% The velocity of the charged particle at r_p

vp = [vp_x,vp_y,0];



%%
%[text] ## Solving for the light vector reaching the charge at x
% Light vector leaving the wire
cv = [c_x, c_y, 0];

%arbritary point in the wire

r_w = [0,y,0];

% the velocity of the wire charges
v_w = [0,v_y,0];

% the light wire reaching the measurement point
eq_1 = r_w + (v_w+cv)*t - r_p == [0,0,0] %[output:4770e46b]
eq_2 = dot(cv,cv)-c^2==0 %[output:0dd89255]

[ccx,ccy,tt,par,cond] = solve([eq_1(1),eq_1(2),eq_2],[c_x,c_y,t],"IgnoreAnalyticConstraints",true,ReturnConditions=true);
%%

c_x = simplify(ccx(2)), c_y = simplify(ccy(2)),tt = simplify(tt(2)) %[output:53386cb5] %[output:2c6ee76c] %[output:5970a22c]
cvec=[c_x,c_y,0];
%%
subs(cvec,[y,x],[0,1]) %[output:771f98c2]

%%
%[text] ## The Force
%[text] 
k_e = q^2/(4*pi*epsilon0);
v_rel1 = cvec + v_w;
v_rel2 = cvec + v_w - vp;
ov = cvec/c;
v_rel1m = dot(v_rel1,v_rel1);
v_rel2m = dot(v_rel2,v_rel2);
R = r_w-r_p;
R2 = dot(R,R);


mag1 = simplify(sqrt(v_rel1m),'Criterion','preferReal','Steps',100); %Speed 
mag2 = simplify(sqrt(v_rel2m),'Criterion','preferReal','Steps',100); %Speed 

p1 = simplify(1/R2,'Criterion','preferReal','Steps',100); % distance factor
p2 = simplify((v_rel1m/c^2)*(1.0/dot(v_rel1,ov)),'Criterion','preferReal','Steps',100); % Density correction
p3 = simplify((v_rel2m/(dot(v_rel2,ov))*(c/sqrt(v_rel1m))),'Criterion','preferReal','Steps',100); % Absortion correction factor

%p4 = simplify((dot(vp,ov)*v_w - dot(vp,v_w)*ov)/c^2,'Criterion','preferReal','Steps',100);
p4=0;
F_a = simplify(p1*p2*p3*(ov+p4),'Criterion','preferReal','Steps',100);



%%
%[text] ## Integrate the static charges
F_n = simplify(subs(F_a,v_y,0),'Criterion','preferReal','Steps',100);
F_nt = simplify(taylor(F_n,vp_x,0,order=2),'Criterion','preferReal','Steps',100);
F_nt = simplify(taylor(F_nt,vp_y,0,order=2),'Criterion','preferReal','Steps',100);
simplify(subs(F_a,[v_y,y,vp_y],[0,0,0]),'Criterion','preferReal','Steps',100) %[output:88e9cb12]
simplify(subs(F_nt,[y,vp_y],[0,0])) %[output:1d65d78b]
syms L real positive

Ft_n = int(F_nt,y,-L,L);
Ft_nInf = simplify(limit(Ft_n,L,Inf),'Criterion','preferReal','Steps',100) %[output:609707b7]


%%
%[text] ## Integrating the current
% Define the velocity vector of the moving charge
F_ct = simplify(taylor(F_a,v_y,0,order=2),'Criterion','preferReal','Steps',100);
F_ct = simplify(taylor(F_ct,vp_y,0,order=2),'Criterion','preferReal','Steps',100);
F_ct = simplify(taylor(F_ct,vp_x,0,order=2),'Criterion','preferReal','Steps',100);

simplify(subs(F_ct,y,0)),[-(q^2*(- c^2 + v_y*vp_y + c*vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c - vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:47d80b9f] %[output:3f3c5c0a]

Ft_c = int(F_ct,y,-L,L);
subs(Ft_c,[vp_x,vp_y],[0,0]) %[output:57bbc25e]
limit(Ft_c,L,0) %[output:5c284957]
%Set L to Infinity
Ft_cInf = simplify(limit(Ft_c,L,Inf),'Criterion','preferReal','Steps',100) %[output:45fd90e0]
%[text] ## The net force
F_atot = k_e*simplify(Ft_cInf - Ft_nInf,'Criterion','preferReal','Steps',100) %[output:23424762]


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:4770e46b]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_1","value":"\\left(\\begin{array}{ccc}\nc_x \\,t-x=0 & y+t\\,{\\left(c_y +v_y \\right)}=0 & 0=0\n\\end{array}\\right)"}}
%---
%[output:0dd89255]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_2","value":"{c_y }^2 +{c_x }^2 -c^2 =0"}}
%---
%[output:53386cb5]
%   data: {"dataType":"symbolic","outputData":{"name":"c_x","value":"-\\frac{x\\,{\\left(v_y \\,y-\\sqrt{c^2 \\,y^2 +c^2 \\,x^2 -{v_y }^2 \\,x^2 }\\right)}}{y^2 +x^2 }"}}
%---
%[output:2c6ee76c]
%   data: {"dataType":"symbolic","outputData":{"name":"c_y","value":"-\\frac{y\\,\\sqrt{c^2 \\,y^2 +c^2 \\,x^2 -{v_y }^2 \\,x^2 }+v_y \\,x^2 }{y^2 +x^2 }"}}
%---
%[output:5970a22c]
%   data: {"dataType":"symbolic","outputData":{"name":"tt","value":"\\frac{v_y \\,y+\\sqrt{c^2 \\,y^2 +c^2 \\,x^2 -{v_y }^2 \\,x^2 }}{c^2 -{v_y }^2 }"}}
%---
%[output:771f98c2]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\sqrt{c^2 -{v_y }^2 } & -v_y  & 0\n\\end{array}\\right)"}}
%---
%[output:88e9cb12]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{\\frac{{\\textrm{vp}}_x }{c}-1}{x^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:1d65d78b]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{c-{\\textrm{vp}}_x }{c\\,x^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:609707b7]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n\\frac{4\\,c-\\pi \\,{\\textrm{vp}}_x }{2\\,c\\,x} & -\\frac{{\\textrm{vp}}_y \\,{\\left(3\\,\\pi \\,c+8\\,{\\textrm{vp}}_x \\right)}}{6\\,c^2 \\,x} & 0\n\\end{array}\\right)"}}
%---
%[output:47d80b9f]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{-c^2 +v_y \\,{\\textrm{vp}}_y +c\\,{\\textrm{vp}}_x }{c^2 \\,x^2 } & -\\frac{v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{c^2 \\,x^2 } & 0\n\\end{array}\\right)"}}
%---
%[output:3f3c5c0a]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y +c\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:57bbc25e]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{2\\,L^{11} \\,c^2 +10\\,L^9 \\,c^2 \\,x^2 +20\\,L^7 \\,c^2 \\,x^4 +20\\,L^5 \\,c^2 \\,x^6 +10\\,L^3 \\,c^2 \\,x^8 +2\\,L\\,c^2 \\,x^{10} }{{\\left(L^2 \\,c^2 \\,x+c^2 \\,x^3 \\right)}\\,{{\\left(L^2 +x^2 \\right)}}^{9\/2} } & -\\frac{L\\,v_y }{L^2 \\,c+c\\,x^2 }-\\frac{L\\,v_y \\,x^2 }{L^2 \\,c\\,x^2 +c\\,x^4 } & 0\n\\end{array}\\right)"}}
%---
%[output:5c284957]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:45fd90e0]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_cInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{-4\\,c^2 +4\\,v_y \\,{\\textrm{vp}}_y +\\pi \\,c\\,{\\textrm{vp}}_x }{2\\,c^2 \\,x} & -\\frac{3\\,\\pi \\,c\\,{\\textrm{vp}}_y -12\\,v_y \\,{\\textrm{vp}}_x +8\\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y }{6\\,c^2 \\,x} & 0\n\\end{array}\\right)"}}
%---
%[output:23424762]
%   data: {"dataType":"symbolic","outputData":{"name":"F_atot","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,c^2 \\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
