%[text] # Electric and magnetic field due to a charge moving at constant velocity
clear;
syms  q  real 
syms t epsilon0 mu0 real positive
syms x y c v_x v_y c_x c_y real positive
assumeAlso (v_x^2 + v_y^2 < c^2)
assumeAlso (c^2 - v_x^2 > 0)
assumeAlso (c^2 - v_y^2 > 0)

syms vp_x vp_y real positive

assumeAlso (vp_x^2+vp_y^2 < c^2)
assumeAlso (c^2 - vp_x^2 > 0)
assumeAlso (c^2 - vp_y^2 > 0)
assumeAlso (c^2*x - v_y^2*x - v_x*sqrt(c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2) + v_x*v_y*y >0 )
assumeAlso (c^2*y - v_x^2*y - v_y*sqrt(c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2) + v_x*v_y*x >0 )


% The emission charge is at zero at t=0
r_o = [0,0,0];
% Define the position of the space where the field will be measured
r_p = [x, y, 0];

% Define the velocity vector of the moving charge
v = [v_x, v_y, 0];

% Light vector
cv = [c_x, c_y, 0];

% The velocity of the reciever charged particle at r_p

vp = [vp_x,vp_y,0];


%%
%[text] ## Compute the time to travel

% the two reach the same point
eq_1=(v+cv)*t==r_p %[output:17824102]
eq_2 = dot(cv,cv)==c^2 %[output:6ceb1ed0]


[ccx,ccy,tt] = solve([eq_1,eq_2],[c_x,c_y,t],"Real",true,"IgnoreAnalyticConstraints",true); %[output:3889e994]
ccx = simplify(ccx(2));
ccy = simplify(ccy(2));
cvec = [-ccx,ccy,0];
assumeAlso (c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2 > 0)
simplify(dot(cvec,cvec)) %[output:86ec241a]

tt = simplify(tt(1));
simplify(subs(ccx,y,0)), simplify(subs(ccy,y,0)), simplify(subs(tt,y,0)) %[output:4e05b108] %[output:673b68f0] %[output:51d43245]
%%

r_c = simplify(expand(v*tt)) %[output:06cf47f5]

simplify(subs(r_c,y,0)) %[output:5c22f2ee]
simplify(subs(r_c,[y,v_x],[0,0])) %[output:1b99f8a5]

%%
%[text] ## Electric and Magnetic Field
assumeAlso (v_x*x - sqrt(c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2) + v_y*y>0)


% Electric field E due to a moving charge
pvec = r_p - r_c;
simplify(subs(pvec,[y,v_x],[0,0])) %[output:5fb7952a]

E = (1/(4*pi*epsilon0)) * (q * pvec) / norm(pvec)^3;
E = simplify(subs(E,y=0));
% Magnetic field B due to a moving charge
B = simplify(subs((1.0/(c^2)) * cross(v,E),y,0));

% Display the equations
simplify(subs(E,[y,v_x],[0,0])), simplify(subs(B,[y,v_x],[0,0])) %[output:3f2b2aca] %[output:92cfaf83]
%%
%[text] ## The Lorenz force



% Define the Lorentz force F acting on a second charge with vp velocity
F = simplify(q * (E + cross(vp, B)),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

Fx = simplify(subs(F,[y,v_x],[0,0]));

Ft = taylor(Fx,v_y,0,order=2) %[output:2b29e235]

%Ft = Fx;
Ft = taylor(Ft,[vp_x,vp_y],[0,0],order=2);

Ft = simplify(Ft,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:50771bda]
F_static = subs(Ft,v_y,0) %[output:8786377a]
Ft_neutral = simplify((Ft - F_static),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:8eec7d0b]
ft_magnetic = simplify(q * cross(vp, B),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = taylor(ft_magnetic,[v_x,v_y],[0,0],order=2);
ft_magnetic = simplify(taylor(ft_magnetic,[vp_x,vp_y],[0,0],order=2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = simplify(subs(ft_magnetic,[v_x,y],[0,0])) %[output:62f2fcba]
%[text] ## The corpuscle force
%[text] 
%[text] The corpuscle force is:
%[text] $&dollar&;&dollar&; \n\\mathbf{F} = \\frac{k q\_1 q\_2}{4 \\pi r^2} \\left( \\frac{\\|\\mathbf{c} + \\mathbf{v\_1}- \\mathbf{v\_2}\\|^2}{\\mathbf{c} \\cdot ( \\mathbf{c} + \\mathbf{v\_1} - \\mathbf{v\_2} )} \\right)  \\hat{o\_1}, \n&dollar&;&dollar&;$
%[text] 
%%
k= q^2/(4*pi*epsilon0);
assumeAlso(c^2 - v_y^2>0)
assumeAlso(c^2 - v_x^2>0)
assumeAlso(c - v_x>0)
assumeAlso(c - v_y>0)
assumeAlso(c^2*y^2 - v_x^2*y^2 + c^2*x^2 >0)

R2 = dot(r_p,r_p);

ov = simplify(cvec/sqrt(dot(cvec,cvec)));

vs = cvec+v;
rvel = vs-vp;

% Dot products
rvelm2 = dot(rvel,rvel);
pvec2 = dot(pvec,pvec);
vsmag2 = dot(vs,vs);


% Assume at y=0 and v_x=0

% corpuscle Vectors

vsmag = simplify(subs(sqrt(vsmag2),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % courpuscle source speed
cmag = simplify(subs(sqrt(rvelm2),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % courpuscle raw speed
rveln = simplify(subs(rvel/cmag,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % courpuscle unitary vector
rvelc = simplify(subs(rvel/c,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % courpuscle relative vector
rvelvs = simplify(subs(rvel/vsmag,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % courpuscle relative vector

dpmp_r = simplify(subs(dot(ov,rvel),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %raw product
dpmp_c = simplify(subs(dot(ov,rvelc),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %c-relative product
dpmp_n = simplify(subs(dot(ov,rveln),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %cosine
dpmp_vs = simplify(subs(dot(ov,rvelvs),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %cosine
dpmp_s = simplify(sqrt(1.0-dpmp_n^2),'Criterion','preferReal','Steps',100); %sine
dpmp_sn = simplify(sqrt(1.0-dpmp_c^2),'Criterion','preferReal','Steps',100); %c-rel sine

% Action Factors
p1 = simplify(subs(k/R2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % distance factor
p2 = simplify(subs((vsmag/(c*cmag))^2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % source velocity correction
ov0 = simplify(subs(ov,[v_x,y],[0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100); % The action orientation

%p3 = simplify(subs(rvelm/dot(cvec,rvel),[v_x,y],[0,0])/p2,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
%E=[(q^2*(c - vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[-(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c - vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)]
%%
%[text] ## 
%[text] ## 
%[text] ## Testing different options for absortion correction
%p3 = cmag^2;
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[sym(0), -(q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:86784a1f] %[output:0d364586]

%p3 = c^2;
E=[(q^2*(c + 2*vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[sym(0), -(q^2*v_y*(c + 2*vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:3bdf1515] %[output:2f9819a8]

%p3 = c*cmag*dpmp_n; 
E=[(q^2*(c + vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c + vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:603b98fe] %[output:20d55730]

%p3 = c*cmag/dpmp_n; 
E=[(q^2*(c + vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[-(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c + vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:62b6d338] %[output:59ffd085]
%p3 = (c/dpmp_n)^2;
E=[(q^2*(c + 2*vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[-(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c + 2*vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:6c161213] %[output:45fda936]

%p3 = c*dpmp_r;
E=[(q^2*(c + vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c + vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:73106f40] %[output:2562a565]

%p3 = cmag*dpmp_r; 
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:96d6e8d9] %[output:89b04bb0]

%p3 = c*cmag*dpmp_s;
E=[(q^2*vp_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[-(q^2*v_y*(c + vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:2aa9fce0] %[output:06875754]

%p3 = cmag^2*dpmp_s^2; 
E=[sym(0), sym(0), sym(0)],M=[-(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), sym(0), sym(0)] %[output:64acbeb4] %[output:754de36e]

%p3 = cmag^2*dpmp_sn^2; 
E=[(q^2*vp_x)/(2*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[-(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*vp_x)/(2*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:0adc392a] %[output:75cfe504]

%p3 = (cmag*dpmp_n)^2;
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:09c60a23] %[output:7084384e]

%p3 = (cmag*dpmp_c)^2;
E=[(q^2*(c - 2*vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c - 2*vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:0cb3f0a4] %[output:5fbeadbe]


%p3=(cmag/dpmp_n)^2;
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[-(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:57903419] %[output:97ef7d5b]
%p3=c/dpmp_c; 

%p3=c/(dpmp_n^(3/2)); 

%%
%[text] ## Combinations
%%p3=cmag*dpmp_r+cmag^2*dpmp_s^2;
%%p3=cmag*(dpmp_r+cmag*dpmp_s^2);
%%p3=cmag*dpmp_r+cmag^2-dpmp_r^2;
%%p3=cmag*(dpmp_r+cmag)-dpmp_r^2;
%p3=dpmp_r*(cmag-dpmp_r)+cmag^2;
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[-(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:6ccba2e0] %[output:174f29a2]

p3=2*cmag^2/dpmp_n-c^2;



%%
eforce = simplify(subs((k/pvec2)*p2*p3*ov0,[v_x,v_y,y],[0,0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
eforce = simplify(taylor(eforce,[vp_x,vp_y],[0,0],order=2)) %[output:682515f6]


F2_c = simplify(p1*p2*p3*ov0,'Criterion','preferReal','Steps',100);
F2T_c = simplify(taylor(F2_c,v_y,0,order=2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

% The net Force due to corpuscles
F2T_net = simplify(taylor(F2T_c,[vp_x,vp_y],[0,0],order=2)) %[output:352684c3]
% The net magnetic Force due to corpuscles
F2T_m=simplify(F2T_net-eforce) %[output:65aa9d2a]


%%
% The Lorentz magnetic forces
Ft_neutral,ft_magnetic %[output:862974f8] %[output:5de0f273]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:17824102]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_1","value":"\\left(\\begin{array}{ccc}\nt\\,{\\left(c_x +v_x \\right)}=x & t\\,{\\left(c_y +v_y \\right)}=y & 0=0\n\\end{array}\\right)"}}
%---
%[output:6ceb1ed0]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_2","value":"{c_y }^2 +{c_x }^2 =c^2"}}
%---
%[output:3889e994]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:86ec241a]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"c^2"}}
%---
%[output:4e05b108]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-\\sqrt{c^2 -{v_y }^2 }"}}
%---
%[output:673b68f0]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-v_y"}}
%---
%[output:51d43245]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\frac{x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 }"}}
%---
%[output:06cf47f5]
%   data: {"dataType":"symbolic","outputData":{"name":"r_c","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{v_x \\,\\sigma_1 }{-c^2 +{v_y }^2 +{v_x }^2 } & \\frac{v_y \\,\\sigma_1 }{-c^2 +{v_y }^2 +{v_x }^2 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =v_x \\,x-\\sqrt{c^2 \\,y^2 -{v_x }^2 \\,y^2 +2\\,v_x \\,v_y \\,x\\,y+c^2 \\,x^2 -{v_y }^2 \\,x^2 }+v_y \\,y\n\\end{array}"}}
%---
%[output:5c22f2ee]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{v_x \\,x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 } & \\frac{v_y \\,x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 } & 0\n\\end{array}\\right)"}}
%---
%[output:1b99f8a5]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & \\frac{v_y \\,x}{\\sqrt{c^2 -{v_y }^2 }} & 0\n\\end{array}\\right)"}}
%---
%[output:5fb7952a]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\nx & -\\frac{v_y \\,x}{\\sqrt{c^2 -{v_y }^2 }} & 0\n\\end{array}\\right)"}}
%---
%[output:3f2b2aca]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{q\\,{{\\left(c^2 -{v_y }^2 \\right)}}^{3\/2} }{4\\,c^3 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q\\,v_y \\,{\\left(c^2 -{v_y }^2 \\right)}}{4\\,c^3 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:92cfaf83]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & 0 & -\\frac{q\\,v_y \\,{{\\left(c^2 -{v_y }^2 \\right)}}^{3\/2} }{4\\,c^5 \\,\\varepsilon_0 \\,x^2 \\,\\pi }\n\\end{array}\\right)"}}
%---
%[output:2b29e235]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi }-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(-c^3 +c^2 \\,{\\textrm{vp}}_x \\right)}}{4\\,c^4 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:50771bda]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:8786377a]
%   data: {"dataType":"symbolic","outputData":{"name":"F_static","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:8eec7d0b]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_neutral","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:62f2fcba]
%   data: {"dataType":"symbolic","outputData":{"name":"ft_magnetic","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:86784a1f]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:0d364586]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n0 & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:3bdf1515]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c+2\\,{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:2f9819a8]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n0 & -\\frac{q^2 \\,v_y \\,{\\left(c+2\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:603b98fe]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c+{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:20d55730]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c+{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:62b6d338]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c+{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:59ffd085]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c+{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:6c161213]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c+2\\,{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:45fda936]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c+2\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:73106f40]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c+{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:2562a565]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c+{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:96d6e8d9]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:89b04bb0]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:2aa9fce0]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\textrm{vp}}_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:06875754]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\left(c+{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:64acbeb4]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n0 & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:754de36e]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:0adc392a]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\textrm{vp}}_x }{2\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:75cfe504]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:09c60a23]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:7084384e]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:0cb3f0a4]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-2\\,{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:5fbeadbe]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-2\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:57903419]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:97ef7d5b]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:6ccba2e0]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:174f29a2]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:682515f6]
%   data: {"dataType":"symbolic","outputData":{"name":"eforce","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-2\\,{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:352684c3]
%   data: {"dataType":"symbolic","outputData":{"name":"F2T_net","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +2\\,v_y \\,{\\textrm{vp}}_y +2\\,c\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-2\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:65aa9d2a]
%   data: {"dataType":"symbolic","outputData":{"name":"F2T_m","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-2\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:862974f8]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_neutral","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:5de0f273]
%   data: {"dataType":"symbolic","outputData":{"name":"ft_magnetic","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
