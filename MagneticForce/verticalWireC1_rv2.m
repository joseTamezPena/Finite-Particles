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
assumeAlso (y^2 + x^2 - v_y^2*y^2 > 0)
assumeAlso (sqrt(y^2 + x^2 - v_y^2*x^2) > v_y*y)


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
%[text] ## The Force
k_e = q^2/(4*pi*epsilon0);
R = r_w-r_p;
R2 = dot(R,R);

v_rel1 = cvec + v_w;
v_rel1m = simplify(dot(v_rel1,v_rel1),'Criterion','preferReal','Steps',100);

v_rel = cvec + v_w - vp;
v_relm = simplify(dot(v_rel,v_rel),'Criterion','preferReal','Steps',100);

p1 = simplify(1/R2,'Criterion','preferReal','Steps',100); % distance factor0
p2 = simplify(v_rel1m/dot(v_rel1,cvec),'Criterion','preferReal','Steps',100); % density and velocity correction
%p2 = 1;
p3 = simplify((v_relm/dot(v_rel,cvec)*(1.0/sqrt(v_rel1m))),'Criterion','preferReal','Steps',100); % Abs correction
%p4 = simplify((dot(vp,cvec)*v_w - dot(vp,v_w)*cvec),'Criterion','preferReal','Steps',100);
p4 = 0;
F_a = k_e*simplify(p1*p2*p3*(cvec+p4),'Criterion','preferReal','Steps',100);


%%
%[text] ## Integrate the static charges
F_n = simplify(subs(F_a,v_y,0));
F_nt = simplify(taylor(F_n,vp_y,0,order=3));
F_nt = simplify(taylor(F_nt,vp_x,0,order=3));
%%
simplify(subs(F_nt,y,0)) %[output:86cdccd6]
syms L real positive

Ft_n = int(F_nt,y,-L,L);
Ft_nInf = simplify(limit(Ft_n,L,Inf)) %[output:095c183f]
Ft_nInf = simplify(taylor(Ft_nInf,vp_y,0,order=2));
Ft_nInf = simplify(taylor(Ft_nInf,vp_x,0,order=2)) %[output:422aa895]
%%
%[text] ## Net forces
netf = simplify(F_a-F_n,'Criterion','preferReal','Steps',100);

Fn_ct = simplify(taylor(netf,v_y,0,order=2));
Fn_ct = simplify(taylor(Fn_ct,vp_y,0,order=2));
Fn_ct = simplify(taylor(Fn_ct,vp_x,0,order=2));
subs(Fn_ct,y,0);
Fnt_c = int(Fn_ct,y,-L,L);
Fnt_cInf = simplify(limit(Fnt_c,L,Inf),'Criterion','preferReal','Steps',100) %[output:782be1b2]

%%
%[text] ## Integrating the current

F_ct = simplify(taylor(F_a,v_y,0,order=2));
F_ct = simplify(taylor(F_ct,vp_x,0,order=2));
F_ct = simplify(taylor(F_ct,vp_y,0,order=2));

Ft_c = int(F_ct,y,-L,L);
%%
%Set L to Infinity
Ft_cInf = simplify(limit(Ft_c,L,Inf)) %[output:59a0a403]
%%
%[text] ## The net force
% The corpuscular force
F_atot = simplify(Ft_cInf - Ft_nInf) %[output:2de8960d]
%The Lorentz force
sol=[-(q^2*v_y*vp_y)/(2*epsilon0*x*sym(pi)), +(q^2*v_y*vp_x)/(2*epsilon0*x*sym(pi)), sym(0)] %[output:0ad9e7ae]

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
%[output:86cdccd6]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(1+{{\\textrm{vp}}_y }^2 -{\\textrm{vp}}_x +{\\textrm{vp}}_x \\,{{\\textrm{vp}}_y }^2 +{{\\textrm{vp}}_x }^2 \\,{{\\textrm{vp}}_y }^2 \\right)}}{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:095c183f]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(240+160\\,{{\\textrm{vp}}_y }^2 -60\\,\\pi \\,{\\textrm{vp}}_x +15\\,\\pi \\,{\\textrm{vp}}_x \\,{{\\textrm{vp}}_y }^2 +80\\,{{\\textrm{vp}}_x }^2 +48\\,{{\\textrm{vp}}_x }^2 \\,{{\\textrm{vp}}_y }^2 \\right)}}{480\\,\\varepsilon_0 \\,x\\,\\pi } & -\\frac{q^2 \\,{\\textrm{vp}}_y \\,{\\left(12\\,\\pi +32\\,{\\textrm{vp}}_x -3\\,\\pi \\,{{\\textrm{vp}}_x }^2 \\right)}}{96\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:422aa895]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-4+\\pi \\,{\\textrm{vp}}_x \\right)}}{8\\,\\varepsilon_0 \\,x\\,\\pi } & -\\frac{q^2 \\,{\\textrm{vp}}_y \\,{\\left(3\\,\\pi +8\\,{\\textrm{vp}}_x \\right)}}{24\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:782be1b2]
%   data: {"dataType":"symbolic","outputData":{"name":"Fnt_cInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:59a0a403]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_cInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-4+4\\,v_y \\,{\\textrm{vp}}_y +\\pi \\,{\\textrm{vp}}_x \\right)}}{8\\,\\varepsilon_0 \\,x\\,\\pi } & -\\frac{q^2 \\,{\\left(3\\,\\pi \\,{\\textrm{vp}}_y -12\\,v_y \\,{\\textrm{vp}}_x +8\\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y \\right)}}{24\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:2de8960d]
%   data: {"dataType":"symbolic","outputData":{"name":"F_atot","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:0ad9e7ae]
%   data: {"dataType":"symbolic","outputData":{"name":"sol","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
