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
v_rel3 = cvec - v_w + vp;
ov = cvec;
v_rel1m = dot(v_rel1,v_rel1);
v_rel2m = simplify(dot(v_rel2,v_rel2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
%v_rel2m = simplify(dot(v_rel2,v_rel3),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)
%v_rel2m = simplify(sqrt(dot(v_rel2,v_rel2))*sqrt(dot(v_rel3,v_rel3)),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);


%p1 = simplify(v_rel2m/R2,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
%p2 = simplify(1.0/dot(cvec,v_rel2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
%F_a = k_e*simplify(p1*p2*ov,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);


p1= simplify(1/R2);
%p2 = (sqrt(dot(v_rel1,v_rel1))*v_rel2m)/dot(ov,v_rel2);
%p2 = (sqrt(v_rel1m*v_rel2m))/(dot(ov,v_rel1/sqrt(v_rel1m)));
%p2 = (sqrt(v_rel1m)/(1*1))*(sqrt(v_rel2m));
%p2 = (sqrt(v_rel1m)/(1*1))*(1/dot(ov,v_rel2/sqrt(v_rel2m)));
%p2 = (sqrt(v_rel1m)/(1*1))*(1/dot(ov,v_rel2/sqrt(v_rel1m)));
%p2 = (sqrt(v_rel1m)/(1*1))*(1/dot(ov,v_rel2));
%p2 = (sqrt(v_rel2m)/(1*1))*(1/dot(ov,v_rel2/sqrt(v_rel2m)));
%p2 = (sqrt(v_rel1m)/(1*1))*(sqrt(v_rel2m)/(sqrt(v_rel1m)*dot(ov,v_rel2/sqrt(v_rel2m))));

%p2 = (sqrt(v_rel1m)/(1*1))*(1.0/dot(ov,v_rel2/sqrt(v_rel2m))^(3/2));


%%%Basics
p2a = simplify(sqrt(v_rel1m)/(1^2),'Criterion','preferReal','Steps',100); % For source at release
cmag = sqrt(v_rel2m); % courpuscle speed
rveln = v_rel2/cmag; % courpuscle unitary vector
p2 = p2a*1/dot(rveln,ov); 


F_a = k_e*simplify(p1*p2*ov,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);


%%
%[text] ## Integrate the static charges
F_n = simplify(subs(F_a,v_y,0));
F_nt = taylor(F_n,[vp_x,vp_y],[0,0],order=2);
%%
syms L real positive

Ft_n = int(F_nt,y,-L,L);
Ft_nInf = simplify(limit(Ft_n,L,Inf)) %[output:1124f7d4]
%%
%[text] ## Net forces
netf = simplify(F_a-F_n);

Fn_ct = simplify(taylor(netf,[vp_x,vp_y],[0,0],order=2));
Fn_ct = simplify(taylor(Fn_ct,v_y,0,order=2));
subs(Fn_ct,y,0) %[output:2644ec52]
Fnt_c = int(Fn_ct,y,-L,L);
Fnt_cInf = simplify(limit(Fnt_c,L,Inf)) %[output:2252b924]
%%
%[text] ## Integrating the current

F_ct = simplify(taylor(F_a,v_y,0,order=2));
F_ct = simplify(taylor(F_ct,[vp_x,vp_y],[0,0],order=2));
Ft_c = int(F_ct,y,-L,L);
%%
%Set L to Infinity
Ft_cInf = simplify(limit(Ft_c,L,Inf)) %[output:63e6a5f1]
%%
%[text] ## The net force
% The corpuscular force
F_atot = simplify(Ft_cInf - Ft_nInf) %[output:8ddb0029]
%The Lorentz force
sol=[-(q^2*v_y*vp_y)/(2*epsilon0*x*sym(pi)), +(q^2*v_y*vp_x)/(2*epsilon0*x*sym(pi)), sym(0)] %[output:7c06eb8d]

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
%[output:1124f7d4]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_nInf","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:2644ec52]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:2252b924]
%   data: {"dataType":"symbolic","outputData":{"name":"Fnt_cInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{3\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{6\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:63e6a5f1]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_cInf","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-3+2\\,v_y \\,{\\textrm{vp}}_y \\right)}}{6\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{6\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:8ddb0029]
%   data: {"dataType":"symbolic","outputData":{"name":"F_atot","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{3\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{6\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:7c06eb8d]
%   data: {"dataType":"symbolic","outputData":{"name":"sol","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,\\varepsilon_0 \\,x\\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,\\varepsilon_0 \\,x\\,\\pi } & 0\n\\end{array}\\right)"}}
%---
