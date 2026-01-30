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
simplify(subs(cvec,[x,y,v_y],[5,5,0])),simplify(subs(cvec,[x,y,v_y],[5,0,0.1])),subs(tt,[x,y,v_y],[5,5,0]) %[output:4fd0f84b] %[output:1610fcc1] %[output:3a11bbd6]

%%
%[text] ## The Force
k_e = q^2/(4*pi*epsilon0);
R = r_w-r_p;
R2 = dot(R,R);

v_rel1 = cvec + v_w;
v_rel2 = cvec + v_w - vp;
v_rel1m = simplify(dot(v_rel1,v_rel1),'Criterion','preferReal','Steps',100);
v_rel2m = simplify(dot(v_rel2,v_rel2),'Criterion','preferReal','Steps',100);

%p1 = simplify(1/R2,'Criterion','preferReal','Steps',100); % Distance to source
%p2 = simplify(v_rel1m,'Criterion','preferReal','Steps',100); % Ajust for number of corpuscles per unit of volume
%p3 = simplify(sqrt(v_rel2m)/dot(v_rel2,cvec),'Criterion','preferReal','Steps',100); % Ajust for number of corpuscles per unit of volume
%p4 = simplify(dot(vp,v_w)*vp + dot(vp,cvec)*v_w,'Criterion','preferReal','Steps',100);
%F_a = k_e*simplify(p1*p2*p3*(cvec+p4),'Criterion','preferReal','Steps',100);


ov = cvec;
vprel = v_w-vp;
%v_rel1m = dot(v_rel1,v_rel1);
%smag = sqrt(v_rel2m);
%cos_t = dot(v_rel2/smag,ov);
%deltac = sqrt(vr_sta2);

%p1= simplify(v_rel1m/R2,'Criterion','preferReal','Steps',100); % Distance to source
%p2 = simplify(smag/cos_t,'Criterion','preferReal','Steps',100); % Ajust for number of corpuscles per unit of volume
%p3 = simplify(1/v_rel1m,'Criterion','preferReal','Steps',100); %Adjust for capture
%F_a = k_e*simplify(p1*p2*p3*ov,'Criterion','preferReal','Steps',100);
%===============

vr_sta2 = simplify(dot(vprel,vprel),'Criterion','preferReal','Steps',100); %Speed
p1 = simplify(1/R2,'Criterion','preferReal','Steps',100); % distance factor
p2 = simplify(v_rel1m,'Criterion','preferReal','Steps',100); % density and velocity correction
p3 = simplify(sqrt(1 - vr_sta2) + sqrt(dot(v_w,ov)^2),'Criterion','preferReal','Steps',100); % density and velocity correction
p4 = simplify((dot(vp,ov)*v_w - 2*dot(vp,v_w)*ov),'Criterion','preferReal','Steps',100);
F_a = k_e*simplify(p1*p2*p3*(ov+p4),'Criterion','preferReal','Steps',100);




%%
%[text] ## Integrate the static charges
%F_a
F_n = simplify(subs(F_a,v_y,0));
F_nt = simplify(taylor(F_n,vp_y,0,order=3));
F_nt = simplify(taylor(F_nt,vp_x,0,order=3));
%%
simplify(subs(F_nt,y,0)) %[output:6bc9e4ac]
syms L real positive

Ft_n = int(F_nt,y,-L,L);
Ft_nInf = simplify(limit(Ft_n,L,Inf)) %[output:0cb5770f]
Ft_nInf = simplify(taylor(Ft_nInf,vp_y,0,order=2));
Ft_nInf = simplify(taylor(Ft_nInf,vp_x,0,order=2)) %[output:309d3a05]
%%
%[text] ## Net forces
netf = simplify(F_a-F_n);

Fn_ct = simplify(taylor(netf,v_y,0,order=3));
Fn_ct = simplify(taylor(Fn_ct,vp_y,0,order=3));
Fn_ct = simplify(taylor(Fn_ct,vp_x,0,order=3));
subs(Fn_ct,y,0);
Fnt_c = int(Fn_ct,y,-L,L);
Fnt_cInf = simplify(limit(Fnt_c,L,Inf),'Criterion','preferReal','Steps',100) %[output:2c392e9c]
Fnt_cInf = simplify(taylor(Fnt_cInf,v_y,0,order=2));
Fnt_cInf = simplify(taylor(Fnt_cInf,vp_y,0,order=2));
Fnt_nInf = simplify(taylor(Fnt_cInf,vp_x,0,order=2)) %[output:204b1d78]

%%
%[text] ## Integrating the current

F_ct = simplify(taylor(F_a,v_y,0,order=2));
F_ct = simplify(taylor(F_ct,vp_x,0,order=2));
F_ct = simplify(taylor(F_ct,vp_y,0,order=2));
simplify(subs(F_ct,y,0)),[-(q^2*(- 1 + v_y*vp_y + vp_x))/(4*epsilon0*x^2*sym(pi)), -(q^2*v_y*(1 - vp_x))/(4*epsilon0*x^2*sym(pi)), sym(0)] %[output:6ed6db58] %[output:24e4a871]
Ft_c = int(F_ct,y,-L,L);
%%
%Set L to Infinity
Ft_cInf = simplify(limit(Ft_c,L,Inf)) %[output:23ffe4f2]
%%
%[text] ## The net force
% The corpuscular force
F_atot = simplify(Ft_cInf - Ft_nInf) %[output:42ef8dc6]
%The Lorentz force
sol=[-(q^2*v_y*vp_y)/(2*epsilon0*x*sym(pi)), +(q^2*v_y*vp_x)/(2*epsilon0*x*sym(pi)), sym(0)] %[output:9c1e7aa3]

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
%[output:3a11bbd6]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\sqrt{50}"}}
%---
%[output:6bc9e4ac]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-4+2\\,{{\\textrm{vp}}_y }^2 +2\\,{{\\textrm{vp}}_x }^2 +{{\\textrm{vp}}_x }^2 \\,{{\\textrm{vp}}_y }^2 \\right)}}{16\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:0cb5770f]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-4+2\\,{{\\textrm{vp}}_y }^2 +2\\,{{\\textrm{vp}}_x }^2 +{{\\textrm{vp}}_x }^2 \\,{{\\textrm{vp}}_y }^2 \\right)}}{8\\,\\varepsilon_0 \\,x\\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:309d3a05]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:2c392e9c]
%   data: {"dataType":"symbolic","outputData":{"name":"Fnt_cInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\left(20\\,v_y +24\\,{\\textrm{vp}}_y +66\\,v_y \\,{{\\textrm{vp}}_y }^2 +6\\,v_y \\,{{\\textrm{vp}}_x }^2 -36\\,{{\\textrm{vp}}_x }^2 \\,{\\textrm{vp}}_y +51\\,v_y \\,{{\\textrm{vp}}_x }^2 \\,{{\\textrm{vp}}_y }^2 \\right)}}{48\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x \\,{\\left(2+2\\,v_y \\,{\\textrm{vp}}_y -{{\\textrm{vp}}_y }^2 \\right)}}{4\\,\\varepsilon_0 \\,x\\,\\pi }-\\frac{q^2 \\,v_y \\,{\\left(-4\\,v_y \\,{\\textrm{vp}}_y +2\\,{{\\textrm{vp}}_y }^2 +2\\,{{\\textrm{vp}}_x }^2 -4\\,v_y \\,{{\\textrm{vp}}_x }^2 \\,{\\textrm{vp}}_y +{{\\textrm{vp}}_x }^2 \\,{{\\textrm{vp}}_y }^2 \\right)}}{32\\,\\varepsilon_0 \\,x} & 0\n\\end{array}\\right)"}}
%---
%[output:204b1d78]
%   data: {"dataType":"symbolic","outputData":{"name":"Fnt_nInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:6ed6db58]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-1+v_y \\,{\\textrm{vp}}_y \\right)}}{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(-1+{\\textrm{vp}}_x \\right)}}{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:24e4a871]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-1+v_y \\,{\\textrm{vp}}_y +{\\textrm{vp}}_x \\right)}}{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(-1+{\\textrm{vp}}_x \\right)}}{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:23ffe4f2]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_cInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-1+v_y \\,{\\textrm{vp}}_y \\right)}}{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:42ef8dc6]
%   data: {"dataType":"symbolic","outputData":{"name":"F_atot","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:9c1e7aa3]
%   data: {"dataType":"symbolic","outputData":{"name":"sol","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
