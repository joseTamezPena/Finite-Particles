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
vpa(subs(cvec,[y,x,v_y,c],[0,10,-0.01,1])) %[output:30db109a]
%%
%[text] ## The Force
%[text] 
k_e = q^2/(4*pi*epsilon0);
v_rel1 = cvec + v_w;
v_rel2 = v_rel1 - vp;
ov = cvec/c;
v_rel1m = dot(v_rel1,v_rel1);
v_rel2m = dot(v_rel2,v_rel2);
R = r_p-r_w;
R2 = dot(R,R);

%--- Working
%p1 = v_rel2m/R2;
%p2 = (ov/dot(cvec,v_rel2));


%p1= simplify(ov/R2);
%p2 = (sqrt(v_rel2m)/(c^2))*sqrt(v_rel2m)/dot(ov,v_rel2);

dpmp_r = dot(ov,v_rel2);

%p1 = (ov/R2)*(v_rel1m/v_rel2m)/(c^2);
p1 = (ov/R2)*(1/v_rel2m);
p2 = dpmp_r*(sqrt(v_rel2m)-dpmp_r)+v_rel2m;

F_a = simplify(p1*p2,'Criterion','preferReal','Steps',100);



simplify(subs(F_a,[v_y,y,x],[0,0,x])) %[output:22845809]
simplify(subs(F_a,[v_y,y,x],[0,0,-x])) %[output:434155ba]
simplify(subs(F_a,[y,vp_x],[0,0])) %[output:797813ec]

%%
%[text] ## Integrate the static charges
F_n = simplify(subs(F_a,v_y,0));
F_nt = simplify(taylor(F_n,[vp_x,vp_y],[0,0],order=2));
simplify(subs(F_a,[v_y,y,vp_y],[0,0,0])) %[output:87993c3e]
simplify(subs(F_nt,[y,vp_y],[0,0])) %[output:25cc57d5]
%%
syms L real positive

Ft_n = int(F_nt,y,-L,L);
Ft_nInf = simplify(limit(Ft_n,L,Inf)) %[output:62333062]
%%
%[text] ## Integrating the current
% Define the velocity vector of the moving charge
F_ct = simplify(taylor(F_a,v_y,0,order=2));
F_ct = simplify(taylor(F_ct,[vp_x,vp_y],[0,0],order=2));
simplify(subs(F_a,[y,vp_x],[0,0])) %[output:25a6b11e]
simplify(subs(F_ct,y,0)) %[output:4d035445]
sol1=[-(- c^2 + v_y*vp_y + c*vp_x)/(c^2*x^2), -(v_y*(c - vp_x))/(c^2*x^2), sym(0)] %[output:24c40381]

Ft_c = int(F_ct,y,-L,L);
%%
%Set L to Infinity
Ft_cInf = simplify(limit(Ft_c,L,Inf)) %[output:5323fa24]
%%
%[text] ## The net force
F_atot = k_e*simplify(Ft_cInf - Ft_nInf) %[output:79bf4b54]


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
%[output:30db109a]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0.99994999874993749609347654199058 & 0.01 & 0\n\\end{array}\\right)"}}
%---
%[output:22845809]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{c\\,\\sqrt{c^2 +{{\\textrm{vp}}_y }^2 -2\\,c\\,{\\textrm{vp}}_x +{{\\textrm{vp}}_x }^2 }-{\\textrm{vp}}_x \\,\\sqrt{c^2 +{{\\textrm{vp}}_y }^2 -2\\,c\\,{\\textrm{vp}}_x +{{\\textrm{vp}}_x }^2 }+{{\\textrm{vp}}_y }^2 }{x^2 \\,{\\left(c^2 +{{\\textrm{vp}}_y }^2 -2\\,c\\,{\\textrm{vp}}_x +{{\\textrm{vp}}_x }^2 \\right)}} & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:434155ba]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{c\\,\\sqrt{c^2 +{{\\textrm{vp}}_y }^2 +2\\,c\\,{\\textrm{vp}}_x +{{\\textrm{vp}}_x }^2 }+{\\textrm{vp}}_x \\,\\sqrt{c^2 +{{\\textrm{vp}}_y }^2 +2\\,c\\,{\\textrm{vp}}_x +{{\\textrm{vp}}_x }^2 }+{{\\textrm{vp}}_y }^2 }{x^2 \\,{\\left(c^2 +{{\\textrm{vp}}_y }^2 +2\\,c\\,{\\textrm{vp}}_x +{{\\textrm{vp}}_x }^2 \\right)}} & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:797813ec]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{\\sqrt{c^2 -{v_y }^2 }\\,{\\left(c^2 -{v_y }^2 +{{\\textrm{vp}}_y }^2 -\\frac{{\\left(c^2 -{v_y }^2 +v_y \\,{\\textrm{vp}}_y \\right)}\\,{\\left(v_y \\,{\\textrm{vp}}_y -c\\,\\sqrt{\\sigma_1 }+c^2 -{v_y }^2 \\right)}}{c^2 }\\right)}}{c\\,x^2 \\,\\sigma_1 } & -\\frac{v_y \\,{\\left(2\\,{v_y }^3 \\,{\\textrm{vp}}_y -{v_y }^2 \\,{{\\textrm{vp}}_y }^2 +c^3 \\,\\sqrt{\\sigma_1 }-{v_y }^4 +c^2 \\,{v_y }^2 +c^2 \\,{{\\textrm{vp}}_y }^2 -2\\,c^2 \\,v_y \\,{\\textrm{vp}}_y -c\\,{v_y }^2 \\,\\sqrt{\\sigma_1 }+c\\,v_y \\,{\\textrm{vp}}_y \\,\\sqrt{\\sigma_1 }\\right)}}{c^3 \\,x^2 \\,\\sigma_1 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =c^2 -{v_y }^2 +{{\\textrm{vp}}_y }^2 \n\\end{array}"}}
%---
%[output:87993c3e]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{1}{x^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:25cc57d5]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{1}{x^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:62333062]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n\\frac{2}{x} & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:25a6b11e]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{\\sqrt{c^2 -{v_y }^2 }\\,{\\left(c^2 -{v_y }^2 +{{\\textrm{vp}}_y }^2 -\\frac{{\\left(c^2 -{v_y }^2 +v_y \\,{\\textrm{vp}}_y \\right)}\\,{\\left(v_y \\,{\\textrm{vp}}_y -c\\,\\sqrt{\\sigma_1 }+c^2 -{v_y }^2 \\right)}}{c^2 }\\right)}}{c\\,x^2 \\,\\sigma_1 } & -\\frac{v_y \\,{\\left(2\\,{v_y }^3 \\,{\\textrm{vp}}_y -{v_y }^2 \\,{{\\textrm{vp}}_y }^2 +c^3 \\,\\sqrt{\\sigma_1 }-{v_y }^4 +c^2 \\,{v_y }^2 +c^2 \\,{{\\textrm{vp}}_y }^2 -2\\,c^2 \\,v_y \\,{\\textrm{vp}}_y -c\\,{v_y }^2 \\,\\sqrt{\\sigma_1 }+c\\,v_y \\,{\\textrm{vp}}_y \\,\\sqrt{\\sigma_1 }\\right)}}{c^3 \\,x^2 \\,\\sigma_1 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =c^2 -{v_y }^2 +{{\\textrm{vp}}_y }^2 \n\\end{array}"}}
%---
%[output:4d035445]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{-c^2 +v_y \\,{\\textrm{vp}}_y }{c^2 \\,x^2 } & -\\frac{v_y }{c\\,x^2 } & 0\n\\end{array}\\right)"}}
%---
%[output:24c40381]
%   data: {"dataType":"symbolic","outputData":{"name":"sol1","value":"\\left(\\begin{array}{ccc}\n-\\frac{-c^2 +v_y \\,{\\textrm{vp}}_y +c\\,{\\textrm{vp}}_x }{c^2 \\,x^2 } & -\\frac{v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{c^2 \\,x^2 } & 0\n\\end{array}\\right)"}}
%---
%[output:5323fa24]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_cInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{-6\\,c^2 +4\\,v_y \\,{\\textrm{vp}}_y }{3\\,c^2 \\,x} & \\frac{v_y \\,{\\left(-3\\,\\pi \\,c+4\\,{\\textrm{vp}}_x \\right)}}{6\\,c^2 \\,x} & 0\n\\end{array}\\right)"}}
%---
%[output:79bf4b54]
%   data: {"dataType":"symbolic","outputData":{"name":"F_atot","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{3\\,c^2 \\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(-3\\,\\pi \\,c+4\\,{\\textrm{vp}}_x \\right)}}{24\\,c^2 \\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
