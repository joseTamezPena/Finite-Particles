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
eq_1 = r_w + (v_w+cv)*t - r_p == [0,0,0] %[output:1ff24391]
eq_2 = dot(cv,cv)-c^2==0 %[output:6bbb3b47]

[ccx,ccy,tt,par,cond] = solve([eq_1(1),eq_1(2),eq_2],[c_x,c_y,t],"IgnoreAnalyticConstraints",true,ReturnConditions=true);
%%

c_x = simplify(ccx(2)), c_y = simplify(ccy(2)),tt = simplify(tt(2)) %[output:09d3b159] %[output:00c1f76f] %[output:10a1a446]
cvec=[c_x,c_y,0];
%%
simplify(subs(cvec,[x,y,v_y],[5,5,0])) %[output:5432fe5e]

%%
%[text] ## The Force
%[text] 
k_e = q^2/(4*pi*epsilon0);
v_rel1 = cvec + v_w;
v_rel2 = v_rel1 - vp;
ov = cvec/c;
v_rel1m = dot(v_rel1,v_rel1);
v_rel2m = dot(v_rel2,v_rel2);
p1 = (v_rel1m/v_rel2m)/c;
p2 = dot(ov,v_rel2);
R = r_p-r_w;
R2 = dot(R,R);
F_a = simplify(p1*p2/R2*ov) %[output:123e3fac]
%%
%[text] ## Integrate the static charges
F_n = simplify(subs(F_a,v_y,0)) %[output:11596177]
F_nt = taylor(F_n,[vp_x,vp_y],[0,0],order=2) %[output:706e7afd]
%%
syms L real positive

Ft_n = int(F_nt,y,-L/2,L/2) %[output:42fde108]
Ft_nInf = simplify(limit(Ft_n,L,Inf)) %[output:2e4fc24a]
%%
%[text] ## Integrating the current
% Define the velocity vector of the moving charge
F_ct = simplify(taylor(F_a,v_y,0,order=2)) %[output:88a8756f]
F_ct = simplify(taylor(F_ct,[vp_x,vp_y],[0,0],order=2)) %[output:9207b679]
Ft_c = int(F_ct,y,-L,L) %[output:30586e83]
%%
%Set L to Infinity
Ft_cInf = simplify(limit(Ft_c,L,Inf)) %[output:361bdefe]
%%
%[text] ## The net force
F_atot = -k_e*simplify(Ft_cInf - Ft_nInf) %[output:9c2f28e2]


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:1ff24391]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_1","value":"\\left(\\begin{array}{ccc}\nc_x \\,t-x=0 & y+t\\,{\\left(c_y +v_y \\right)}=0 & 0=0\n\\end{array}\\right)"}}
%---
%[output:6bbb3b47]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_2","value":"{c_y }^2 +{c_x }^2 -c^2 =0"}}
%---
%[output:09d3b159]
%   data: {"dataType":"symbolic","outputData":{"name":"c_x","value":"-\\frac{x\\,{\\left(v_y \\,y-\\sqrt{c^2 \\,y^2 +c^2 \\,x^2 -{v_y }^2 \\,x^2 }\\right)}}{y^2 +x^2 }"}}
%---
%[output:00c1f76f]
%   data: {"dataType":"symbolic","outputData":{"name":"c_y","value":"-\\frac{y\\,\\sqrt{c^2 \\,y^2 +c^2 \\,x^2 -{v_y }^2 \\,x^2 }+v_y \\,x^2 }{y^2 +x^2 }"}}
%---
%[output:10a1a446]
%   data: {"dataType":"symbolic","outputData":{"name":"tt","value":"\\frac{v_y \\,y+\\sqrt{c^2 \\,y^2 +c^2 \\,x^2 -{v_y }^2 \\,x^2 }}{c^2 -{v_y }^2 }"}}
%---
%[output:5432fe5e]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{\\sqrt{2}\\,c}{2} & -\\frac{\\sqrt{2}\\,c}{2} & 0\n\\end{array}\\right)"}}
%---
%[output:123e3fac]
%   data: {"dataType":"symbolic","outputData":{"name":"F_a","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n-\\frac{x\\,{\\left(v_y \\,y-\\sigma_4 \\right)}\\,\\sigma_3 \\,\\sigma_1 }{\\sigma_2 } & -\\frac{{\\left(y\\,\\sigma_4 +v_y \\,x^2 \\right)}\\,\\sigma_3 \\,\\sigma_1 }{\\sigma_2 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =c^2 \\,x^2 -{v_y }^2 \\,x^2 +c^2 \\,y^2 -{\\textrm{vp}}_x \\,x\\,\\sigma_4 -v_y \\,y\\,\\sigma_4 +{\\textrm{vp}}_y \\,y\\,\\sigma_4 +v_y \\,{\\textrm{vp}}_y \\,x^2 +v_y \\,{\\textrm{vp}}_x \\,x\\,y\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =c^3 \\,{{\\left(y^2 +x^2 \\right)}}^4 \\,{\\left(\\frac{{{\\left(y\\,\\sigma_4 +{\\textrm{vp}}_y \\,x^2 -v_y \\,y^2 +{\\textrm{vp}}_y \\,y^2 \\right)}}^2 }{{{\\left(y^2 +x^2 \\right)}}^2 }+{{\\left({\\textrm{vp}}_x +\\frac{x\\,{\\left(v_y \\,y-\\sigma_4 \\right)}}{y^2 +x^2 }\\right)}}^2 \\right)}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 ={v_y }^2 \\,y^2 -{v_y }^2 \\,x^2 +c^2 \\,x^2 +c^2 \\,y^2 -2\\,v_y \\,y\\,\\sigma_4 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 =\\sqrt{c^2 \\,y^2 +c^2 \\,x^2 -{v_y }^2 \\,x^2 }\n\\end{array}"}}
%---
%[output:11596177]
%   data: {"dataType":"symbolic","outputData":{"name":"F_n","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{c\\,x\\,{\\left(c\\,x^2 +c\\,y^2 -{\\textrm{vp}}_x \\,x\\,\\sigma_2 +{\\textrm{vp}}_y \\,y\\,\\sigma_2 \\right)}}{\\sigma_1 } & -\\frac{c\\,y\\,{\\left(c\\,x^2 +c\\,y^2 -{\\textrm{vp}}_x \\,x\\,\\sigma_2 +{\\textrm{vp}}_y \\,y\\,\\sigma_2 \\right)}}{\\sigma_1 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 ={{\\left(y^2 +x^2 \\right)}}^{5\/2} \\,{\\left(\\frac{{{\\left({\\textrm{vp}}_x \\,x^2 +{\\textrm{vp}}_x \\,y^2 -c\\,x\\,\\sigma_2 \\right)}}^2 }{{{\\left(y^2 +x^2 \\right)}}^2 }+\\frac{{{\\left({\\textrm{vp}}_y \\,x^2 +{\\textrm{vp}}_y \\,y^2 +c\\,y\\,\\sigma_2 \\right)}}^2 }{{{\\left(y^2 +x^2 \\right)}}^2 }\\right)}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =\\sqrt{y^2 +x^2 }\n\\end{array}"}}
%---
%[output:706e7afd]
%   data: {"dataType":"symbolic","outputData":{"name":"F_nt","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{c\\,{\\textrm{vp}}_y \\,x\\,\\sigma_2 }{\\sigma_1 }-\\frac{c\\,{\\textrm{vp}}_x \\,x\\,\\sigma_3 }{\\sigma_1 }+\\frac{c\\,x\\,{\\left(c\\,y^2 +c\\,x^2 \\right)}}{\\sigma_1 \\,\\sigma_4 } & \\frac{c\\,{\\textrm{vp}}_x \\,y\\,\\sigma_3 }{\\sigma_1 }-\\frac{c\\,{\\textrm{vp}}_y \\,y\\,\\sigma_2 }{\\sigma_1 }-\\frac{c\\,y\\,{\\left(c\\,y^2 +c\\,x^2 \\right)}}{\\sigma_1 \\,\\sigma_4 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 ={{\\left(y^2 +x^2 \\right)}}^{5\/2} \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =\\frac{y\\,\\sqrt{y^2 +x^2 }}{\\sigma_4 }-\\frac{2\\,c\\,y\\,{\\left(c\\,y^2 +c\\,x^2 \\right)}}{\\sqrt{y^2 +x^2 }\\,{\\sigma_4 }^2 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =\\frac{x\\,\\sqrt{y^2 +x^2 }}{\\sigma_4 }-\\frac{2\\,c\\,x\\,{\\left(c\\,y^2 +c\\,x^2 \\right)}}{\\sqrt{y^2 +x^2 }\\,{\\sigma_4 }^2 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 =\\frac{c^2 \\,x^2 }{y^2 +x^2 }+\\frac{c^2 \\,y^2 }{y^2 +x^2 }\n\\end{array}"}}
%---
%[output:42fde108]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_n","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{{\\textrm{vp}}_x \\,\\sigma_2 \\,{\\sigma_1 }^{5\/2} +2\\,L\\,c\\,{\\sigma_1 }^2 +2\\,L\\,{\\textrm{vp}}_x \\,x\\,{\\sigma_1 }^{3\/2} }{x\\,{\\left(4\\,x^2 \\,{\\left(2\\,L^2 \\,c\\,\\sqrt{\\sigma_1 }+4\\,c\\,x^2 \\,\\sqrt{\\sigma_1 }\\right)}+L^4 \\,c\\,\\sqrt{\\sigma_1 }\\right)}} & \\frac{{\\textrm{vp}}_y \\,\\sigma_2 }{c\\,x}-\\frac{2\\,L\\,{\\textrm{vp}}_y }{c\\,\\sigma_1 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =L^2 +4\\,x^2 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =\\mathrm{atan}\\left(\\frac{L}{2\\,x}\\right)\n\\end{array}"}}
%---
%[output:2e4fc24a]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n\\frac{2\\,c+\\frac{\\pi \\,{\\textrm{vp}}_x }{2}}{c\\,x} & \\frac{\\pi \\,{\\textrm{vp}}_y }{2\\,c\\,x} & 0\n\\end{array}\\right)"}}
%---
%[output:88a8756f]
%   data: {"dataType":"symbolic","outputData":{"name":"F_ct","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{c\\,x\\,\\sigma_1 }{\\sigma_5 }-\\frac{v_y \\,x\\,{\\left(c\\,{\\left(c^2 \\,\\sigma_8 \\,{\\left(\\frac{y}{\\sigma_7 }-\\frac{2\\,c\\,y\\,\\sigma_3 }{\\sigma_6 }\\right)}+\\frac{2\\,c^2 \\,y\\,\\sigma_8 }{\\sigma_7 }\\right)}\\,\\sigma_1 -\\frac{c^3 \\,{\\sigma_8 }^{3\/2} \\,\\sigma_2 }{\\sigma_7 }\\right)}}{\\sigma_4 } & -\\frac{v_y \\,{\\left(c\\,{\\left(c^2 \\,\\sigma_8 \\,{\\left(\\frac{x^2 }{\\sigma_7 }+\\frac{2\\,c\\,y^2 \\,\\sigma_3 }{\\sigma_6 }\\right)}-\\frac{2\\,c^2 \\,y^2 \\,\\sigma_8 }{\\sigma_7 }\\right)}\\,\\sigma_1 +\\frac{c^3 \\,y\\,{\\sigma_8 }^{3\/2} \\,\\sigma_2 }{\\sigma_7 }\\right)}}{\\sigma_4 }-\\frac{c\\,y\\,\\sigma_1 }{\\sigma_5 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =c\\,x^2 +c\\,y^2 -{\\textrm{vp}}_x \\,x\\,\\sqrt{\\sigma_8 }+{\\textrm{vp}}_y \\,y\\,\\sqrt{\\sigma_8 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 ={\\textrm{vp}}_y \\,x^2 -c\\,y\\,\\sqrt{\\sigma_8 }+{\\textrm{vp}}_x \\,x\\,y\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =c\\,\\sqrt{\\sigma_8 }-{\\textrm{vp}}_x \\,x+{\\textrm{vp}}_y \\,y\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 =c^3 \\,{\\sigma_8 }^4 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_5 ={\\sigma_8 }^{5\/2} \\,\\sigma_7 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_6 =\\sqrt{\\sigma_8 }\\,{\\sigma_7 }^2 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_7 =\\frac{{{\\left({\\textrm{vp}}_x \\,x^2 +{\\textrm{vp}}_x \\,y^2 -c\\,x\\,\\sqrt{\\sigma_8 }\\right)}}^2 }{{\\sigma_8 }^2 }+\\frac{{{\\left({\\textrm{vp}}_y \\,x^2 +{\\textrm{vp}}_y \\,y^2 +c\\,y\\,\\sqrt{\\sigma_8 }\\right)}}^2 }{{\\sigma_8 }^2 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_8 =y^2 +x^2 \n\\end{array}"}}
%---
%[output:9207b679]
%   data: {"dataType":"symbolic","outputData":{"name":"F_ct","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{x\\,{\\left(c^2 \\,\\sigma_1 +c\\,{\\textrm{vp}}_x \\,x-2\\,c\\,v_y \\,y-c\\,{\\textrm{vp}}_y \\,y+v_y \\,{\\textrm{vp}}_y \\,\\sigma_1 \\right)}}{\\sigma_2 } & -\\frac{c^2 \\,y\\,\\sigma_1 +c\\,v_y \\,x^2 -c\\,v_y \\,y^2 -c\\,{\\textrm{vp}}_y \\,y^2 +c\\,{\\textrm{vp}}_x \\,x\\,y+v_y \\,{\\textrm{vp}}_x \\,x\\,\\sigma_1 }{\\sigma_2 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\sqrt{y^2 +x^2 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =c^2 \\,{{\\left(y^2 +x^2 \\right)}}^2 \n\\end{array}"}}
%---
%[output:30586e83]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_c","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{L\\,{\\left(2\\,c^2 \\,\\sigma_2 +2\\,v_y \\,{\\textrm{vp}}_y \\,\\sigma_2 +c\\,{\\textrm{vp}}_x \\,x\\right)}}{\\sigma_1 }+\\frac{{\\textrm{vp}}_x \\,\\mathrm{atan}\\left(\\frac{L}{x}\\right)}{c\\,x} & -\\frac{2\\,L\\,v_y \\,{\\textrm{vp}}_x \\,\\sigma_2 +2\\,L\\,c\\,v_y \\,x+L\\,c\\,{\\textrm{vp}}_y \\,x-c\\,{\\textrm{vp}}_y \\,\\mathrm{atan}\\left(\\frac{L}{x}\\right)\\,{\\left(L^2 +x^2 \\right)}}{\\sigma_1 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =x\\,{\\left(L^2 \\,c^2 +c^2 \\,x^2 \\right)}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =\\sqrt{L^2 +x^2 }\n\\end{array}"}}
%---
%[output:361bdefe]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_cInf","value":"\\left(\\begin{array}{ccc}\n\\frac{4\\,c^2 +4\\,v_y \\,{\\textrm{vp}}_y +\\pi \\,c\\,{\\textrm{vp}}_x }{2\\,c^2 \\,x} & -\\frac{-\\frac{\\pi \\,c\\,{\\textrm{vp}}_y }{2}+2\\,v_y \\,{\\textrm{vp}}_x }{c^2 \\,x} & 0\n\\end{array}\\right)"}}
%---
%[output:9c2f28e2]
%   data: {"dataType":"symbolic","outputData":{"name":"F_atot","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,c^2 \\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
