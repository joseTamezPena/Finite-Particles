%[text] # A single vertical wire with current
clear;

syms  q  real 
syms t epsilon0 mu0 real positive
syms x y v_y c_x real 
syms c_y real 
syms vp_x vp_y real 


assumeAlso (t > 0)
assumeAlso (x > 0)
assumeAlso (y > 0)

assumeAlso (v_y > 0)
assumeAlso (c_x > 0)
assumeAlso (c_y < 0)
assumeAlso (vp_x > 0)
assumeAlso (vp_y > 0)
assumeAlso (1 - vp_x > 0)
assumeAlso (1 - vp_y > 0)
assumeAlso (1 - v_y^2 > 0)
assumeAlso (1 - vp_y^2 - vp_y^2 > 0)
assumeAlso (y^2 + x^2 - v_y^2*x^2 > 0)


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
eq_1 = r_w + (v_w+cv)*t - r_p == [0,0,0] %[output:0e0af0c6]
eq_2 = dot(cv,cv)-1^2==0 %[output:1fcaf330]

[ccx,ccy,tt,par,cond] = solve([eq_1(1),eq_1(2),eq_2],[c_x,c_y,t],"IgnoreAnalyticConstraints",true,ReturnConditions=true);
%%

c_x = simplify(ccx(2)),c_y = simplify(ccy(2)),tt = simplify(tt(2)) %[output:4ca87760] %[output:637e3edd] %[output:63b52d71]


cvec=[c_x,c_y,0];
%%
simplify(subs(cvec,[x,y,v_y],[5,5,0])),simplify(subs(cvec,[x,y,v_y],[5,0,0.1])) %[output:4fd0f84b] %[output:1610fcc1]

%%
%[text] ## The Force
k_e = q^2/(4*pi*epsilon0);
R = r_p-r_w;
R2 = dot(R,R);

v_rel1 = cvec + v_w;
v_rel2 = v_rel1 - vp;
ov = cvec;
v_rel1m = dot(v_rel1,v_rel1);
v_rel2m = dot(v_rel2,v_rel2);


p1 = simplify(v_rel2m/R2,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
p2 = simplify(ov/dot(cvec,v_rel2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
F_a = k_e*simplify(p1*p2,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);


%%
%[text] ## Integrate the static charges
F_n = simplify(subs(F_a,v_y,0));
F_nt = taylor(F_n,[vp_x,vp_y],[0,0],order=2);
%%
syms L real positive

Ft_n = int(F_nt,y,-L,L);
Ft_nInf = simplify(limit(Ft_n,L,Inf)) %[output:5c4da99c]
%%
%[text] ## Net forces
netf = simplify(F_a-F_n);

Fn_ct = simplify(taylor(netf,[vp_x,vp_y],[0,0],order=2));
Fn_ct = simplify(taylor(Fn_ct,v_y,0,order=2));
subs(Fn_ct,y,0) %[output:3717c36f]
Fnt_c = int(Fn_ct,y,-L,L);
Fnt_cInf = simplify(limit(Fnt_c,L,Inf)) %[output:2187b570]
%%
%[text] ## Integrating the current
% Define the velocity vector of the moving charge
F_ct = simplify(taylor(F_a,v_y,0,order=2));
F_ct = simplify(taylor(F_ct,[vp_x,vp_y],[0,0],order=2));
Ft_c = int(F_ct,y,-L,L);
%%
%Set L to Infinity
Ft_cInf = simplify(limit(Ft_c,L,Inf)) %[output:33008a39]
%%
%[text] ## The net force
F_atot = simplify(Ft_cInf - Ft_nInf) %[output:72f17789]
sol=[-(q^2*v_y*vp_y)/(2*epsilon0*x*sym(pi)), +(q^2*v_y*vp_x)/(2*epsilon0*x*sym(pi)), sym(0)] %[output:767a0fa9]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:0e0af0c6]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_1","value":"\\left(\\begin{array}{ccc}\nc_x \\,t-x=0 & y+t\\,{\\left(c_y +v_y \\right)}=0 & 0=0\n\\end{array}\\right)"}}
%---
%[output:1fcaf330]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_2","value":"-1+{c_y }^2 +{c_x }^2 =0"}}
%---
%[output:4ca87760]
%   data: {"dataType":"symbolic","outputData":{"name":"c_x","value":"-\\frac{x\\,{\\left(v_y \\,y-\\sqrt{y^2 +x^2 -{v_y }^2 \\,x^2 }\\right)}}{y^2 +x^2 }"}}
%---
%[output:637e3edd]
%   data: {"dataType":"symbolic","outputData":{"name":"c_y","value":"-\\frac{y\\,\\sqrt{y^2 +x^2 -{v_y }^2 \\,x^2 }+v_y \\,x^2 }{y^2 +x^2 }"}}
%---
%[output:63b52d71]
%   data: {"dataType":"symbolic","outputData":{"name":"tt","value":"-\\frac{v_y \\,y+\\sqrt{y^2 +x^2 -{v_y }^2 \\,x^2 }}{-1+{v_y }^2 }"}}
%---
%[output:4fd0f84b]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{\\sqrt{2}}{2} & -\\frac{\\sqrt{2}}{2} & 0\n\\end{array}\\right)"}}
%---
%[output:1610fcc1]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{3\\,\\sqrt{11}}{10} & -\\frac{1}{10} & 0\n\\end{array}\\right)"}}
%---
%[output:5c4da99c]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-4+\\pi \\,{\\textrm{vp}}_x \\right)}}{8\\,\\varepsilon_0 \\,x\\,\\pi } & -\\frac{q^2 \\,{\\textrm{vp}}_y }{8\\,\\varepsilon_0 \\,x} & 0\n\\end{array}\\right)"}}
%---
%[output:3717c36f]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(-x^3 +{\\textrm{vp}}_x \\,x^3 \\right)}}{4\\,\\varepsilon_0 \\,x^5 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:2187b570]
%   data: {"dataType":"symbolic","outputData":{"name":"Fnt_cInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:33008a39]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_cInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-4+4\\,v_y \\,{\\textrm{vp}}_y +\\pi \\,{\\textrm{vp}}_x \\right)}}{8\\,\\varepsilon_0 \\,x\\,\\pi } & -\\frac{q^2 \\,{\\left(\\pi \\,{\\textrm{vp}}_y -4\\,v_y \\,{\\textrm{vp}}_x \\right)}}{8\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:72f17789]
%   data: {"dataType":"symbolic","outputData":{"name":"F_atot","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:767a0fa9]
%   data: {"dataType":"symbolic","outputData":{"name":"sol","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
