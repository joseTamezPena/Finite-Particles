%[text] # A single vertical wire with current

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
eq_1 = r_w + (v_w+cv)*t - r_p == [0,0,0] %[output:833e62db]
eq_2 = dot(cv,cv)-1^2==0 %[output:5326bc5f]

[ccx,ccy,tt,par,cond] = solve([eq_1(1),eq_1(2),eq_2],[c_x,c_y,t],"IgnoreAnalyticConstraints",true,ReturnConditions=true);
%%

c_x = simplify(ccx(2)),c_y = simplify(ccy(2)),tt = simplify(tt(2)) %[output:95186bb8] %[output:0b4e5fe9] %[output:0266cab0]



%%
simplify(subs(cvec,[x,y,v_y],[5,5,0])) %[output:44192c96]

%%
%[text] ## The Force
k_e = q^2/(4*pi*epsilon0);
v_rel1 = cvec + v_w;
v_rel2 = v_rel1 - vp;
ov = cvec;
v_rel1m = sqrt(dot(v_rel1,v_rel1));
v_rel2m = sqrt(dot(v_rel2,v_rel2));
p1 = v_rel1m;
%p2 = dot(ov,v_rel2/v_rel2m)*dot(ov,v_rel2/v_rel2m);
%p2 = dot(ov,v_rel2/v_rel2m)*dot(v_rel1/v_rel1m,v_rel2/v_rel2m);
p2 = dot(ov,v_rel2/v_rel2m);
%p3 = dot(ov,v_rel2/v_rel2m)-dot(ov,v_rel1/v_rel1m)/2;
%p3 = dot(ov,v_rel1/v_rel1m);
%p3 = v_rel1m;
%p3 = (dot(ov,v_rel2/v_rel2m)+v_rel1m)/2;
%p3 = (dot(ov,v_rel2/v_rel2m)+dot(ov,v_rel1/v_rel1m))/2;
%p3 = dot(ov,v_rel2/v_rel2m)*v_rel2m;
%p3 = v_rel2m;
%p3 = v_rel2m/v_rel1m;
%p3 = dot(ov,v_rel2/v_rel2m)/v_rel1m;
p3 = v_rel1m/v_rel2m;
R = r_p-r_w;
R2 = dot(R,R);
F_a = simplify(p1*p2*p3/R2*ov) %[output:61a42253]
%%
%[text] ## Integrate the static charges
F_n = simplify(subs(F_a,v_y,0)) %[output:5f906a2d]
F_nt = taylor(F_n,[vp_x,vp_y],[0,0],order=2) %[output:9e56b267]
%%
syms L real positive

Ft_n = int(F_nt,y,-L,L) %[output:4d7da4b9]
Ft_nInf = simplify(limit(Ft_n,L,Inf)) %[output:2dc71a1d]
%%
%[text] ## Integrating the current
% Define the velocity vector of the moving charge
F_ct = simplify(taylor(F_a,v_y,0,order=2)) %[output:9ab85781]
F_ct = simplify(taylor(F_ct,[vp_x,vp_y],[0,0],order=2)) %[output:55cf698f]
Ft_c = int(F_ct,y,-L,L) %[output:808fcaaf]
%%
%Set L to Infinity
Ft_cInf = simplify(limit(Ft_c,L,Inf)) %[output:468b879b]
%%
%[text] ## The net force
F_atot = -k_e*simplify(Ft_cInf - Ft_nInf) %[output:0d9d4bc9]


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:833e62db]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_1","value":"\\left(\\begin{array}{ccc}\nc_x \\,t-x=0 & y+t\\,{\\left(c_y +v_y \\right)}=0 & 0=0\n\\end{array}\\right)"}}
%---
%[output:5326bc5f]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_2","value":"-1+{c_y }^2 +{c_x }^2 =0"}}
%---
%[output:95186bb8]
%   data: {"dataType":"symbolic","outputData":{"name":"c_x","value":"-\\frac{x\\,{\\left(v_y \\,y-\\sqrt{y^2 +x^2 -{v_y }^2 \\,x^2 }\\right)}}{y^2 +x^2 }"}}
%---
%[output:0b4e5fe9]
%   data: {"dataType":"symbolic","outputData":{"name":"c_y","value":"-\\frac{y\\,\\sqrt{y^2 +x^2 -{v_y }^2 \\,x^2 }+v_y \\,x^2 }{y^2 +x^2 }"}}
%---
%[output:0266cab0]
%   data: {"dataType":"symbolic","outputData":{"name":"tt","value":"-\\frac{v_y \\,y+\\sqrt{y^2 +x^2 -{v_y }^2 \\,x^2 }}{-1+{v_y }^2 }"}}
%---
%[output:44192c96]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{\\sqrt{2}}{2} & -\\frac{\\sqrt{2}}{2} & 0\n\\end{array}\\right)"}}
%---
%[output:61a42253]
%   data: {"dataType":"symbolic","outputData":{"name":"F_a","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n-\\frac{x\\,{\\left(v_y \\,y-\\sigma_4 \\right)}\\,\\sigma_3 \\,\\sigma_2 }{\\sigma_1 } & -\\frac{{\\left(y\\,\\sigma_4 +v_y \\,x^2 \\right)}\\,\\sigma_3 \\,\\sigma_2 }{\\sigma_1 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 ={{\\left(y^2 +x^2 \\right)}}^2 \\,{\\left({{\\left(y\\,\\sigma_4 +{\\textrm{vp}}_y \\,x^2 -v_y \\,y^2 +{\\textrm{vp}}_y \\,y^2 \\right)}}^2 +{{\\left({\\textrm{vp}}_x \\,x^2 -x\\,\\sigma_4 +{\\textrm{vp}}_x \\,y^2 +v_y \\,x\\,y\\right)}}^2 \\right)}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =x^2 -{v_y }^2 \\,x^2 +y^2 +v_y \\,{\\textrm{vp}}_y \\,x^2 -{\\textrm{vp}}_x \\,x\\,\\sigma_4 -v_y \\,y\\,\\sigma_4 +{\\textrm{vp}}_y \\,y\\,\\sigma_4 +v_y \\,{\\textrm{vp}}_x \\,x\\,y\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 ={v_y }^2 \\,y^2 -{v_y }^2 \\,x^2 +x^2 +y^2 -2\\,v_y \\,y\\,\\sigma_4 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 =\\sqrt{y^2 +x^2 -{v_y }^2 \\,x^2 }\n\\end{array}"}}
%---
%[output:5f906a2d]
%   data: {"dataType":"symbolic","outputData":{"name":"F_n","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{x\\,{\\left(x^2 +y^2 -{\\textrm{vp}}_x \\,x\\,\\sqrt{y^2 +x^2 }+{\\textrm{vp}}_y \\,y\\,\\sqrt{y^2 +x^2 }\\right)}}{\\sigma_1 } & -\\frac{y\\,{\\left(x^2 +y^2 -{\\textrm{vp}}_x \\,x\\,\\sqrt{y^2 +x^2 }+{\\textrm{vp}}_y \\,y\\,\\sqrt{y^2 +x^2 }\\right)}}{\\sigma_1 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\sqrt{y^2 +x^2 }\\,{\\left({{\\left({\\textrm{vp}}_x \\,x^2 -x\\,\\sqrt{y^2 +x^2 }+{\\textrm{vp}}_x \\,y^2 \\right)}}^2 +{{\\left(y\\,\\sqrt{y^2 +x^2 }+{\\textrm{vp}}_y \\,x^2 +{\\textrm{vp}}_y \\,y^2 \\right)}}^2 \\right)}\n\\end{array}"}}
%---
%[output:9e56b267]
%   data: {"dataType":"symbolic","outputData":{"name":"F_nt","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{x\\,\\sigma_4 }{\\sigma_3 }-\\frac{{\\textrm{vp}}_x \\,x\\,\\sigma_2 }{\\sigma_4 }+\\frac{{\\textrm{vp}}_y \\,x\\,\\sigma_1 }{\\sigma_4 } & \\frac{{\\textrm{vp}}_x \\,y\\,\\sigma_2 }{\\sigma_4 }-\\frac{y\\,\\sigma_4 }{\\sigma_3 }-\\frac{{\\textrm{vp}}_y \\,y\\,\\sigma_1 }{\\sigma_4 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{y\\,\\sigma_4 }{\\sigma_3 }-\\frac{2\\,y\\,{{\\left(y^2 +x^2 \\right)}}^{5\/2} }{{\\sigma_3 }^2 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =\\frac{x\\,\\sigma_4 }{\\sigma_3 }-\\frac{2\\,x\\,{{\\left(y^2 +x^2 \\right)}}^{5\/2} }{{\\sigma_3 }^2 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =x^2 \\,{\\left(y^2 +x^2 \\right)}+y^2 \\,{\\left(y^2 +x^2 \\right)}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 =\\sqrt{y^2 +x^2 }\n\\end{array}"}}
%---
%[output:4d7da4b9]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_n","value":"\\left(\\begin{array}{ccc}\n\\frac{L\\,{\\textrm{vp}}_x }{L^2 +x^2 }+\\frac{{\\textrm{vp}}_x \\,\\mathrm{atan}\\left(\\frac{L}{x}\\right)}{x}+\\frac{2\\,L}{x\\,\\sqrt{L^2 +x^2 }} & \\frac{{\\textrm{vp}}_y \\,\\mathrm{atan}\\left(\\frac{L}{x}\\right)}{x}-\\frac{L\\,{\\textrm{vp}}_y }{L^2 +x^2 } & 0\n\\end{array}\\right)"}}
%---
%[output:2dc71a1d]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n\\frac{4+\\pi \\,{\\textrm{vp}}_x }{2\\,x} & \\frac{\\pi \\,{\\textrm{vp}}_y }{2\\,x} & 0\n\\end{array}\\right)"}}
%---
%[output:9ab85781]
%   data: {"dataType":"symbolic","outputData":{"name":"F_ct","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{x\\,\\sigma_2 }{\\sqrt{\\sigma_5 }\\,\\sigma_1 }-\\frac{v_y \\,x\\,{\\left({\\left(\\sigma_5 \\,{\\left(\\frac{y}{\\sigma_1 }-\\frac{2\\,y\\,{\\sigma_5 }^{3\/2} \\,\\sigma_4 }{{\\sigma_1 }^2 }\\right)}+\\frac{2\\,y\\,\\sigma_5 }{\\sigma_1 }\\right)}\\,\\sigma_2 -\\frac{{\\sigma_5 }^{3\/2} \\,\\sigma_3 }{\\sigma_1 }\\right)}}{{\\sigma_5 }^2 } & -\\frac{v_y \\,{\\left({\\left(\\sigma_5 \\,{\\left(\\frac{x^2 }{\\sigma_1 }+\\frac{2\\,y^2 \\,{\\sigma_5 }^{3\/2} \\,\\sigma_4 }{{\\sigma_1 }^2 }\\right)}-\\frac{2\\,y^2 \\,\\sigma_5 }{\\sigma_1 }\\right)}\\,\\sigma_2 +\\frac{y\\,{\\sigma_5 }^{3\/2} \\,\\sigma_3 }{\\sigma_1 }\\right)}}{{\\sigma_5 }^2 }-\\frac{y\\,\\sigma_2 }{\\sqrt{\\sigma_5 }\\,\\sigma_1 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 ={{\\left({\\textrm{vp}}_x \\,x^2 -x\\,\\sqrt{\\sigma_5 }+{\\textrm{vp}}_x \\,y^2 \\right)}}^2 +{{\\left(y\\,\\sqrt{\\sigma_5 }+{\\textrm{vp}}_y \\,x^2 +{\\textrm{vp}}_y \\,y^2 \\right)}}^2 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =x^2 +y^2 -{\\textrm{vp}}_x \\,x\\,\\sqrt{\\sigma_5 }+{\\textrm{vp}}_y \\,y\\,\\sqrt{\\sigma_5 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 ={\\textrm{vp}}_y \\,x^2 -y\\,\\sqrt{\\sigma_5 }+{\\textrm{vp}}_x \\,x\\,y\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 ={\\textrm{vp}}_y \\,y-{\\textrm{vp}}_x \\,x+\\sqrt{\\sigma_5 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_5 =y^2 +x^2 \n\\end{array}"}}
%---
%[output:55cf698f]
%   data: {"dataType":"symbolic","outputData":{"name":"F_ct","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{x\\,{\\left(x^2 +y^2 +v_y \\,{\\textrm{vp}}_y \\,x^2 +v_y \\,{\\textrm{vp}}_y \\,y^2 +{\\textrm{vp}}_x \\,x\\,\\sqrt{\\sigma_1 }-2\\,v_y \\,y\\,\\sqrt{\\sigma_1 }-{\\textrm{vp}}_y \\,y\\,\\sqrt{\\sigma_1 }\\right)}}{{\\sigma_1 }^{5\/2} } & -\\frac{x^2 \\,y+y^3 +v_y \\,x^2 \\,\\sqrt{\\sigma_1 }-v_y \\,y^2 \\,\\sqrt{\\sigma_1 }-{\\textrm{vp}}_y \\,y^2 \\,\\sqrt{\\sigma_1 }+v_y \\,{\\textrm{vp}}_x \\,x^3 +{\\textrm{vp}}_x \\,x\\,y\\,\\sqrt{\\sigma_1 }+v_y \\,{\\textrm{vp}}_x \\,x\\,y^2 }{{\\sigma_1 }^{5\/2} } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =y^2 +x^2 \n\\end{array}"}}
%---
%[output:808fcaaf]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_c","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{2\\,L^3 }{\\sigma_2 }+\\frac{2\\,L\\,x}{{\\sigma_3 }^{3\/2} }+\\frac{{\\textrm{vp}}_x \\,\\mathrm{atan}\\left(\\frac{L}{x}\\right)}{x}+\\frac{L\\,{\\textrm{vp}}_x \\,x^2 }{\\sigma_1 }+\\frac{2\\,L^3 \\,v_y \\,{\\textrm{vp}}_y }{\\sigma_2 }+\\frac{2\\,L\\,v_y \\,{\\textrm{vp}}_y \\,x}{{\\sigma_3 }^{3\/2} } & \\frac{{\\textrm{vp}}_y \\,\\mathrm{atan}\\left(\\frac{L}{x}\\right)}{x}-\\frac{L\\,{\\textrm{vp}}_y }{\\sigma_3 }-\\frac{L\\,v_y }{\\sigma_3 }-\\frac{L\\,v_y \\,x^2 }{\\sigma_1 }-\\frac{2\\,L^3 \\,v_y \\,{\\textrm{vp}}_x }{\\sigma_2 }-\\frac{2\\,L\\,v_y \\,{\\textrm{vp}}_x \\,x}{{\\sigma_3 }^{3\/2} } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =L^2 \\,x^2 +x^4 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =x\\,{\\sigma_3 }^{3\/2} \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =L^2 +x^2 \n\\end{array}"}}
%---
%[output:468b879b]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_cInf","value":"\\left(\\begin{array}{ccc}\n\\frac{4+4\\,v_y \\,{\\textrm{vp}}_y +\\pi \\,{\\textrm{vp}}_x }{2\\,x} & \\frac{\\pi \\,{\\textrm{vp}}_y -4\\,v_y \\,{\\textrm{vp}}_x }{2\\,x} & 0\n\\end{array}\\right)"}}
%---
%[output:0d9d4bc9]
%   data: {"dataType":"symbolic","outputData":{"name":"F_atot","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
