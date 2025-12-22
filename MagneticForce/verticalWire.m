%[text] # A single vertical wire with current

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
simplify(subs(cvec,[x,y,v_y],[1,-1,0])) %[output:5763abb0]
simplify(subs(cvec,[x,y,v_y],[1,0,-0.01*c])) %[output:039885b1]
%%
%[text] ## The Force
%[text] 
k_e = q^2/(4*pi*epsilon0);
v_rel1 = cvec + v_w;
v_rel2 = v_rel1 - vp;
ov = cvec/c;
v_rel1m = sqrt(dot(v_rel1,v_rel1));
v_rel2m = sqrt(dot(v_rel2,v_rel2));
p1 = v_rel1m/c;
p2 = (1.0/(dot(ov,v_rel2/v_rel2m)*dot(ov,v_rel1/v_rel1m)));
R = r_p-r_w;
R2 = dot(R,R);
F_a = simplify(p1*p2/R2*ov);
simplify(subs(F_a,[v_y,y,x],[0,0,x])) %[output:32a8f242]
simplify(subs(F_a,[v_y,y,x],[0,0,-x])) %[output:0112436b]

%%
%[text] ## Integrate the static charges
F_n = simplify(subs(F_a,v_y,0));
F_nt = taylor(F_n,[vp_x,vp_y],[0,0],order=2);
simplify(subs(F_a,[v_y,y,vp_y],[0,0,0])) %[output:2cdd2279]
simplify(subs(F_nt,[y,vp_y],[0,0])) %[output:839ce010]
%%
syms L real positive

Ft_n = int(F_nt,y,-L,L);
Ft_nInf = simplify(limit(Ft_n,L,Inf))
%%
%[text] ## Integrating the current
% Define the velocity vector of the moving charge
F_ct = simplify(taylor(F_a,v_y,0,order=2));
F_ct = simplify(taylor(F_ct,[vp_x,vp_y],[0,0],order=2));
simplify(subs(F_a,[y,vp_y],[0,0]));
simplify(subs(F_ct,[y,vp_y],[0,0]))
Ft_c = int(F_ct,y,-L,L);
%%
%Set L to Infinity
Ft_cInf = simplify(limit(Ft_c,L,Inf))
%%
%[text] ## The net force
F_atot = k_e*simplify(Ft_cInf - Ft_nInf)


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
%[output:5763abb0]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{\\sqrt{2}\\,c}{2} & \\frac{\\sqrt{2}\\,c}{2} & 0\n\\end{array}\\right)"}}
%---
%[output:039885b1]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{3\\,\\sqrt{1111}\\,c}{100} & \\frac{c}{100} & 0\n\\end{array}\\right)"}}
%---
%[output:32a8f242]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{\\sqrt{c^2 +{{\\textrm{vp}}_y }^2 -2\\,c\\,{\\textrm{vp}}_x +{{\\textrm{vp}}_x }^2 }}{x^2 \\,{\\left(c-{\\textrm{vp}}_x \\right)}} & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:0112436b]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{\\sqrt{c^2 +{{\\textrm{vp}}_y }^2 +2\\,c\\,{\\textrm{vp}}_x +{{\\textrm{vp}}_x }^2 }}{x^2 \\,{\\left(c+{\\textrm{vp}}_x \\right)}} & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:2cdd2279]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{1}{x^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:839ce010]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{1}{x^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
