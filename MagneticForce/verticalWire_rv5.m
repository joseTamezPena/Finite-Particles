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
eq_1 = r_w + (v_w+cv)*t - r_p == [0,0,0] %[output:6dac617e]
eq_2 = dot(cv,cv)-c^2==0 %[output:587d2e34]

[ccx,ccy,tt,par,cond] = solve([eq_1(1),eq_1(2),eq_2],[c_x,c_y,t],"IgnoreAnalyticConstraints",true,ReturnConditions=true);
%%

c_x = simplify(ccx(2)), c_y = simplify(ccy(2)),tt = simplify(tt(2)) %[output:12118750] %[output:1e40c523] %[output:97cb076c]
cvec=[c_x,c_y,0];
%%
    subs(cvec,[y,x],[0,1]) %[output:3f487f1c]

%%
%[text] ## The Force
%[text] 
k_e = q^2/(4*pi*epsilon0);
vrel = (v_w - vp)/c;
v_rel1 = cvec + v_w;
v_rel2 = cvec + v_w - vp;
ov = cvec/c;
vrelm = sqrt(dot(vrel,vrel));
v_rel1m = dot(v_rel1,v_rel1);
v_rel2m = dot(v_rel2,v_rel2);
R = r_w-r_p;
R2 = dot(R,R);


mag1 = simplify(sqrt(v_rel1m),'Criterion','preferReal','Steps',100); %Speed 
mag2 = simplify(sqrt(v_rel2m),'Criterion','preferReal','Steps',100); %Speed 
urvel = v_rel2/mag2;
dots = simplify(dot(v_rel1,ov),'Criterion','preferReal','Steps',100);
coss = simplify(dot(urvel,ov),'Criterion','preferReal','Steps',100);

p1 = simplify(1/R2,'Criterion','preferReal','Steps',100); % distance factor
p2 = simplify((v_rel1m/c^2),'Criterion','preferReal','Steps',100); % Density correction
p3 = simplify(mag2/dots,'Criterion','preferReal','Steps',100);
%p4 = simplify(dot((vp-v_w)/mag2,ov)*ov + dot(vp/c,v_w/c)*vp/c + p3*dot(vp/c,ov)*v_w/c,'Criterion','preferReal','Steps',100);

%p4 = simplify( -dot(vrel,urvel)*urvel,'Criterion','preferReal','Steps',100);
%[(q^2*v_y*vp_y*(16*c - 9*sym(pi)*vp_x))/(48*c^3*epsilon0*x*sym(pi)), (q^2*v_y*vp_x)/(3*c^2*epsilon0*x*sym(pi)), sym(0)]

%p4 = simplify( -dot(vrel,urvel)*ov,'Criterion','preferReal','Steps',100);
%[(q^2*v_y*vp_y*(64*c - 3*sym(pi)*vp_x))/(96*c^3*epsilon0*x*sym(pi)), (q^2*v_y*vp_x)/(6*c^2*epsilon0*x*sym(pi)), sym(0)]

%p4 = simplify( -dot(vrel,ov)*urvel,'Criterion','preferReal','Steps',100);
%[-(q^2*v_y*vp_y*(32*c + 3*sym(pi)*vp_x))/(96*c^3*epsilon0*x*sym(pi)), (2*q^2*v_y*vp_x)/(3*c^2*epsilon0*x*sym(pi)), sym(0)]

%p4 = simplify( -dot(vrel,ov)*ov,'Criterion','preferReal','Steps',100);
[-(q^2*v_y*vp_x*vp_y)/(16*c^3*epsilon0*x), (q^2*v_y*vp_x)/(2*c^2*epsilon0*x*sym(pi)), sym(0)] %[output:07d4b63b]

%p4 = 0;
[-(q^2*v_y*vp_y*(32*c + 3*sym(pi)*vp_x))/(96*c^3*epsilon0*x*sym(pi)), (q^2*v_y*(3*sym(pi)*c + 4*vp_x))/(24*c^2*epsilon0*x*sym(pi)), sym(0)] %[output:05284df3]
%p4 = simplify( -dot(vrel,ov)*ov + dot(vrel,urvel)*ov/2,'Criterion','preferReal','Steps',100);
%[-(q^2*v_y*vp_y*(8*c + sym(pi)*vp_x))/(16*c^3*epsilon0*x*sym(pi)), (q^2*v_y*(sym(pi)*c + 8*vp_x))/(16*c^2*epsilon0*x*sym(pi)), sym(0)]

%p4 = simplify( dot(vrel,v_rel2/c)*v_rel2/c,'Criterion','preferReal','Steps',100);
%[-(q^2*v_y*vp_y*(5*c - 3*sym(pi)*vp_x))/(3*c^3*epsilon0*x*sym(pi)), -(q^2*v_y*(- 3*sym(pi)*c + 8*vp_x))/(12*c^2*epsilon0*x*sym(pi)), sym(0)]


%p4 = simplify( dot(vp/c,urvel)*ov,'Criterion','preferReal','Steps',100);
%[(q^2*v_y*vp_y*(16*c - 3*sym(pi)*vp_x))/(96*c^3*epsilon0*x*sym(pi)), (q^2*v_y*(3*sym(pi)*c + 4*vp_x))/(24*c^2*epsilon0*x*sym(pi)), sym(0)]

%p4 = simplify( dot(vp/c,ov)*ov,'Criterion','preferReal','Steps',100);
%[-(q^2*v_y*vp_y*(16*c + 9*sym(pi)*vp_x))/(96*c^3*epsilon0*x*sym(pi)), (q^2*v_y*(3*sym(pi)*c + 8*vp_x))/(24*c^2*epsilon0*x*sym(pi)), sym(0)]

%p4 = simplify( dot(vp/c,ov)*urvel,'Criterion','preferReal','Steps',100);
%[-(q^2*v_y*vp_y*(32*c + 3*sym(pi)*vp_x))/(96*c^3*epsilon0*x*sym(pi)), (q^2*v_y*(3*sym(pi)*c + 16*vp_x))/(24*c^2*epsilon0*x*sym(pi)), sym(0)]

p4 = simplify( vrelm*ov,'Criterion','preferReal','Steps',100);

F_a = simplify(p1*p2*p3*(ov+p4),'Criterion','preferReal','Steps',100);

%%

%%
%[text] ## Integrate the static charges
F_n = simplify(subs(F_a,v_y,0),'Criterion','preferReal','Steps',100);
F_nt = simplify(taylor(F_n,vp_x,0,order=2),'Criterion','preferReal','Steps',100);
F_nt = simplify(taylor(F_nt,vp_y,0,order=2),'Criterion','preferReal','Steps',100);
simplify(subs(F_a,[v_y,y,vp_y],[0,0,0]),'Criterion','preferReal','Steps',100) %[output:1a603e1c]
simplify(subs(F_nt,[y,vp_y],[0,0])) %[output:42fad958]
syms L real positive

Ft_n = int(F_nt,y,-L,L);
Ft_nInf = simplify(limit(Ft_n,L,Inf),'Criterion','preferReal','Steps',100) %[output:5ed0cbab]


%%
%[text] ## Integrating the current
% Define the velocity vector of the moving charge
F_ct = simplify(taylor(F_a,v_y,0,order=2),'Criterion','preferReal','Steps',100);
F_ct = simplify(taylor(F_ct,vp_y,0,order=2),'Criterion','preferReal','Steps',100);
F_ct = simplify(taylor(F_ct,vp_x,0,order=2),'Criterion','preferReal','Steps',100); %[output:1e755d69]

simplify(subs(k_e*F_ct,y,0)),q^2*[-((- c^2 + v_y*vp_y + c*vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), -(v_y*(c - vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)]

Ft_c = int(F_ct,y,-L,L);
subs(Ft_c,[vp_x,vp_y],[0,0]);
limit(Ft_c,L,0)
%Set L to Infinity
Ft_cInf = simplify(limit(Ft_c,L,Inf),'Criterion','preferReal','Steps',100);
%[text] ## The net force
F_atot = k_e*simplify(Ft_cInf - Ft_nInf,'Criterion','preferReal','Steps',100)


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:6dac617e]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_1","value":"\\left(\\begin{array}{ccc}\nc_x \\,t-x=0 & y+t\\,{\\left(c_y +v_y \\right)}=0 & 0=0\n\\end{array}\\right)"}}
%---
%[output:587d2e34]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_2","value":"{c_y }^2 +{c_x }^2 -c^2 =0"}}
%---
%[output:12118750]
%   data: {"dataType":"symbolic","outputData":{"name":"c_x","value":"-\\frac{x\\,{\\left(v_y \\,y-\\sqrt{c^2 \\,y^2 +c^2 \\,x^2 -{v_y }^2 \\,x^2 }\\right)}}{y^2 +x^2 }"}}
%---
%[output:1e40c523]
%   data: {"dataType":"symbolic","outputData":{"name":"c_y","value":"-\\frac{y\\,\\sqrt{c^2 \\,y^2 +c^2 \\,x^2 -{v_y }^2 \\,x^2 }+v_y \\,x^2 }{y^2 +x^2 }"}}
%---
%[output:97cb076c]
%   data: {"dataType":"symbolic","outputData":{"name":"tt","value":"\\frac{v_y \\,y+\\sqrt{c^2 \\,y^2 +c^2 \\,x^2 -{v_y }^2 \\,x^2 }}{c^2 -{v_y }^2 }"}}
%---
%[output:3f487f1c]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\sqrt{c^2 -{v_y }^2 } & -v_y  & 0\n\\end{array}\\right)"}}
%---
%[output:07d4b63b]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y }{16\\,c^3 \\,\\varepsilon_0 \\,x} & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,c^2 \\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:05284df3]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y \\,{\\left(32\\,c+3\\,\\pi \\,{\\textrm{vp}}_x \\right)}}{96\\,c^3 \\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(3\\,\\pi \\,c+4\\,{\\textrm{vp}}_x \\right)}}{24\\,c^2 \\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:1a603e1c]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{c^2 -{{\\textrm{vp}}_x }^2 }{c^2 \\,x^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:42fad958]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{c-{\\textrm{vp}}_x }{c\\,x^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:5ed0cbab]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n\\frac{{\\left(c+{\\textrm{vp}}_y \\right)}\\,{\\left(4\\,c-\\pi \\,{\\textrm{vp}}_x \\right)}}{2\\,c^2 \\,x} & -\\frac{{\\textrm{vp}}_y \\,{\\left(3\\,\\pi \\,c+4\\,{\\textrm{vp}}_x \\right)}}{6\\,c^2 \\,x} & 0\n\\end{array}\\right)"}}
%---
%[output:1e755d69]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('sym\/taylor', 'C:\\Program Files\\MATLAB\\R2025b\\toolbox\\symbolic\\symbolic\\@sym\\taylor.m', 127)\" style=\"font-weight:bold\">sym\/taylor<\/a> (<a href=\"matlab: opentoline('C:\\Program Files\\MATLAB\\R2025b\\toolbox\\symbolic\\symbolic\\@sym\\taylor.m',127,0)\">line 127<\/a>)\nUnable to compute Taylor expansion. Try 'series' for more general expansion."}}
%---
