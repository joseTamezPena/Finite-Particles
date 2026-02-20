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
ov = cvec/c;
v_rel1m = dot(v_rel1,v_rel1);
R = r_w-r_p;
R2 = dot(R,R);


cmag = simplify(sqrt(v_rel1m),'Criterion','preferReal','Steps',100); %Speed 

p1 = simplify(1/R2,'Criterion','preferReal','Steps',100); % distance factor
p2 = simplify(cmag/c,'Criterion','preferReal','Steps',100); % from source factor

p4 = simplify((dot(vp,ov)*v_w - dot(vp,v_w)*ov)/c^2,'Criterion','preferReal','Steps',100);
F_a = simplify(p1*p2*(ov+p4),'Criterion','preferReal','Steps',100);



%%
%[text] ## Integrate the static charges
F_n = simplify(subs(F_a,v_y,0),'Criterion','preferReal','Steps',100);
F_nt = simplify(taylor(F_n,vp_x,0,order=2),'Criterion','preferReal','Steps',100);
F_nt = simplify(taylor(F_nt,vp_y,0,order=2),'Criterion','preferReal','Steps',100);
simplify(subs(F_a,[v_y,y,vp_y],[0,0,0]),'Criterion','preferReal','Steps',100) %[output:69794517]
simplify(subs(F_nt,[y,vp_y],[0,0])) %[output:08035b3f]
syms L real positive

Ft_n = int(F_nt,y,-L,L);
Ft_nInf = simplify(limit(Ft_n,L,Inf),'Criterion','preferReal','Steps',100) %[output:8e92f737]


%%
%[text] ## Integrating the current
% Define the velocity vector of the moving charge
F_ct = simplify(taylor(F_a,v_y,0,order=2),'Criterion','preferReal','Steps',100);
F_ct = simplify(taylor(F_ct,vp_y,0,order=2),'Criterion','preferReal','Steps',100);
F_ct = simplify(taylor(F_ct,vp_x,0,order=2),'Criterion','preferReal','Steps',100);

simplify(subs(F_ct,y,0)),[-(q^2*(- c^2 + v_y*vp_y + c*vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c - vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:2a5ce6ac] %[output:9d84fd44]

Ft_c = int(F_ct,y,-L,L);
subs(Ft_c,[vp_x,vp_y],[0,0]) %[output:0e3db015]
limit(Ft_c,L,0) %[output:3c283673]
%Set L to Infinity
Ft_cInf = simplify(limit(Ft_c,L,Inf),'Criterion','preferReal','Steps',100) %[output:2e7d4978]
%[text] ## The net force
F_atot = k_e*simplify(Ft_cInf - Ft_nInf,'Criterion','preferReal','Steps',100) %[output:14e192bb]


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
%[output:69794517]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{1}{x^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:08035b3f]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{1}{x^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:8e92f737]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n\\frac{2}{x} & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:2a5ce6ac]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{-c^2 +v_y \\,{\\textrm{vp}}_y }{c^2 \\,x^2 } & -\\frac{v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{c^2 \\,x^2 } & 0\n\\end{array}\\right)"}}
%---
%[output:9d84fd44]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y +c\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:0e3db015]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{2\\,L}{x\\,\\sqrt{L^2 +x^2 }} & -\\frac{2\\,L\\,v_y }{c\\,{\\left(L^2 +x^2 \\right)}} & 0\n\\end{array}\\right)"}}
%---
%[output:3c283673]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:2e7d4978]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_cInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{2\\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y \\right)}}{c^2 \\,x} & \\frac{2\\,v_y \\,{\\textrm{vp}}_x }{c^2 \\,x} & 0\n\\end{array}\\right)"}}
%---
%[output:14e192bb]
%   data: {"dataType":"symbolic","outputData":{"name":"F_atot","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,c^2 \\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
