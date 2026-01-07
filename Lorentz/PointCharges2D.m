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
eq_1=(v+cv)*t==r_p %[output:49f885b4]
eq_2 = dot(cv,cv)==c^2 %[output:2aeb5593]


[ccx,ccy,tt] = solve([eq_1,eq_2],[c_x,c_y,t],"Real",true,"IgnoreAnalyticConstraints",true); %[output:160968dd]
ccx = simplify(ccx(2));
ccy = simplify(ccy(2));
cvec = [-ccx,ccy,0];
assumeAlso (c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2 > 0)
simplify(dot(cvec,cvec)) %[output:88484ec5]

tt = simplify(tt(1));
simplify(subs(ccx,y,0)), simplify(subs(ccy,y,0)), simplify(subs(tt,y,0)) %[output:62dfa894] %[output:7ff2eb3c] %[output:95bd8d4b]
%%

r_c = simplify(expand(v*tt)) %[output:14169921]

simplify(subs(r_c,y,0)) %[output:0ff0c22f]
simplify(subs(r_c,[y,v_x],[0,0])) %[output:9c15096d]

%%
%[text] ## Electric and Magnetic Field
assumeAlso (v_x*x - sqrt(c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2) + v_y*y>0)


% Electric field E due to a moving charge
pvec = r_p - r_c;
simplify(subs(pvec,[y,v_x],[0,0])) %[output:5a8ed240]

E = (1/(4*pi*epsilon0)) * (q * pvec) / norm(pvec)^3;
E = simplify(subs(E,y=0));
% Magnetic field B due to a moving charge
B = simplify(subs((1.0/(c^2)) * cross(v,E),y,0));

% Display the equations
simplify(subs(E,[y,v_x],[0,0])), simplify(subs(B,[y,v_x],[0,0])) %[output:73749778] %[output:9963bf93]
%%
%[text] ## The Lorenz force



% Define the Lorentz force F acting on a second charge with vp velocity
F = simplify(q * (E + cross(vp, B)),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

Fx = simplify(subs(F,[y,v_x],[0,0]));

Ft = taylor(Fx,v_y,0,order=2) %[output:715ae53c]

%Ft = Fx;
Ft = taylor(Ft,[vp_x,vp_y],[0,0],order=2);

Ft = simplify(Ft,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:6a6a9ea7]
F_static = subs(Ft,v_y,0) %[output:363b55de]
Ft_neutral = simplify((Ft - F_static),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:51db289c]
ft_magnetic = simplify(q * cross(vp, B),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = taylor(ft_magnetic,[v_x,v_y],[0,0],order=2);
ft_magnetic = simplify(taylor(ft_magnetic,[vp_x,vp_y],[0,0],order=2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = simplify(subs(ft_magnetic,[v_x,y],[0,0])) %[output:40305650]
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
rvelm = dot(rvel,rvel);
pvec2 = dot(pvec,pvec);
vsmag = dot(vs,vs);


% Assume at y=0 and v_x=0

% corpuscle Vectors

cmag = simplify(subs(sqrt(rvelm),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % courpuscle raw speed
rveln = simplify(subs(rvel/cmag,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % courpuscle unitary vector
rvelc = simplify(subs(rvel/c,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % courpuscle relative vector

dpmp_r = simplify(subs(dot(ov,rvel),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %raw product
dpmp_c = simplify(subs(dot(ov,rvelc),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %c-relative product
dpmp_n = simplify(subs(dot(ov,rveln),[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); %cosine
dpmp_s = simplify(sqrt(1.0-dpmp_n^2),'Criterion','preferReal','Steps',100); %sine
dpmp_sn = simplify(sqrt(1.0-dpmp_c^2),'Criterion','preferReal','Steps',100); %c-rel sine

% Action Factors
p1 = simplify(subs(k/R2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % distance factor
p2 = simplify(subs(sqrt(vsmag)/c^2,[v_x,y],[0,0]),'Criterion','preferReal','Steps',100); % source velocity correction
ov0 = simplify(subs(ov,[v_x,y],[0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100); % The action orientation

p3 = simplify(subs(rvelm/dot(cvec,rvel),[v_x,y],[0,0])/p2,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
E=[(q^2*(c - vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[-(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c - vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:213e5620] %[output:876f187a]
%%
%[text] ## 
%[text] ## 
%[text] ## Testing different options for absortion correction
%p3 = cmag;
E=[(q^2*(c - vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[sym(0), -(q^2*v_y*(c - vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:0758925b] %[output:169b257a]
%p3 = c;
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[sym(0), -(q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:102f85cd] %[output:4edeb801]
%p3 = c*dpmp_n; 
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:2eee45e2] %[output:94668975]
%p3 = c/dpmp_n; %or cmag/dpmp_c;
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[-(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:0d82eeb8] %[output:40259810]


%p3 = dpmp_r; 
E=[(q^2*(c - vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c - vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:17072a30] %[output:4f41e5bd]
%p3 = cmag*dpmp_c; 
E=[(q^2*(c - 2*vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c - 2*vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:5087ec0b] %[output:7c96db09]
%p3 = cmag*dpmp_s;
E=[(q^2*vp_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[-(q^2*v_y*(c - vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:93ccdba8] %[output:7b1fe95b]
%p3 = cmag*dpmp_s^2; % or c*dpmp_s^2;
E=[sym(0), sym(0), sym(0)], M= [-(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), sym(0), sym(0)] %[output:71f0c583] %[output:5a873ac6]
%p3 = cmag*dpmp_sn^2; % or c*dpmp_sn^2
E=[(q^2*vp_x)/(2*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[-(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*vp_x)/(2*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:408a2ee3] %[output:39e2d897]

%p3 = cmag*dpmp_n^2;
E=[(q^2*(c - vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c - vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:87855c9a] %[output:27b05733]
%p3 = cmag*dpmp_c^2;
E=[(q^2*(c - 3*vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c - 3*vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:2d3f3cfc] %[output:8b087497]

p3=cmag/dpmp_n; 
E=[(q^2*(c - vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[-(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c - vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:12513be8] %[output:432de499]

%p3=c/dpmp_c; 
E=[(q^2*(c + vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[-(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c + vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:0c1642a5] %[output:2a7f71ff]
%p3=c/(dpmp_n^(3/2)); 
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[(q^2*(- 4*c^2 + 3*v_y*vp_y))/(8*c^2*epsilon0*x^2*sym(pi)), (q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:2dbfede1] %[output:44e46089]

%%
%[text] ## Combinations
%p3=cmag/dpmp_n+c/dpmp_c;
E=[q^2/(2*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[-(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y)/(2*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:30162d12] %[output:886607cb]
%p3=c*dpmp_n+c/dpmp_n;
E=[q^2/(2*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[sym(0), -(q^2*v_y)/(2*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:8e82a98e] %[output:117e9327]
%p3= c/dpmp_n - c*dpmp_n;
E=[sym(0), sym(0), sym(0)],M=[(q^2*v_y*vp_y)/(2*c^2*epsilon0*x^2*sym(pi)), sym(0), sym(0)] %[output:0640159d] %[output:4ee36d78]
%p3= dpmp_r - c*dpmp_n; %or cmag*dpmp_c - dpmp_r
E=[-(q^2*vp_x)/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[sym(0), (q^2*v_y*vp_x)/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:5d553da2] %[output:8a0065ff]
%p3=c/dpmp_c + dpmp_r - c*dpmp_n;
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[-(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:5d3a5a83] %[output:1648a7d3]
%p3=cmag*(dpmp_c^2-dpmp_n^2);
E=[-(q^2*vp_x)/(2*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[sym(0), (q^2*v_y*vp_x)/(2*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:9566a7c8] %[output:687d4529]
%p3=cmag*(dpmp_n^2+dpmp_sn^2/2);
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:27a1f563] %[output:394b0f03]
%p3=cmag-dpmp_r;
E=[sym(0), sym(0), sym(0)], M=[-(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), sym(0), sym(0)] %[output:08f088c5] %[output:73d49389]
%p3=cmag+cmag*dpmp_s^2/2;
%p3=cmag*(1+dpmp_s^2/2);
%p3=cmag*(1.5-dpmp_n^2/2);
E=[(q^2*(c - vp_x))/(4*c*epsilon0*x^2*sym(pi)), sym(0), sym(0)], M=[-(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y*(c - vp_x))/(4*c^2*epsilon0*x^2*sym(pi)), sym(0)] %[output:8d68895c] %[output:76ccc2f1]
%p3=2*dpmp_r-cmag*dpmp_c;
E=[q^2/(4*epsilon0*x^2*sym(pi)), sym(0), sym(0)],M=[(q^2*v_y*vp_y)/(4*c^2*epsilon0*x^2*sym(pi)), -(q^2*v_y)/(4*c*epsilon0*x^2*sym(pi)), sym(0)] %[output:9b323e14] %[output:2f218cbb]
%%
eforce = simplify(subs((k/pvec2)*p2*p3*ov0,[v_x,v_y,y],[0,0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
eforce = simplify(taylor(eforce,[vp_x,vp_y],[0,0],order=2)) %[output:4ee7b89f]


F2_c = simplify(p1*p2*p3*ov0,'Criterion','preferReal','Steps',100);
subs(F2_c,vp_x,0) %[output:09e56986]
F2T_c = simplify(taylor(F2_c,v_y,0,order=2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

% The net Force due to corpuscles
F2T_net = simplify(taylor(F2T_c,[vp_x,vp_y],[0,0],order=2)) %[output:820d5b48]
% The net magnetic Force due to corpuscles
F2T_m=simplify(F2T_net-eforce) %[output:7dbeb18f]


%%
% The Lorentz magnetic forces
Ft_neutral,ft_magnetic %[output:2a7d5919] %[output:8a477894]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:49f885b4]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_1","value":"\\left(\\begin{array}{ccc}\nt\\,{\\left(c_x +v_x \\right)}=x & t\\,{\\left(c_y +v_y \\right)}=y & 0=0\n\\end{array}\\right)"}}
%---
%[output:2aeb5593]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_2","value":"{c_y }^2 +{c_x }^2 =c^2"}}
%---
%[output:160968dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:88484ec5]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"c^2"}}
%---
%[output:62dfa894]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-\\sqrt{c^2 -{v_y }^2 }"}}
%---
%[output:7ff2eb3c]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-v_y"}}
%---
%[output:95bd8d4b]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\frac{x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 }"}}
%---
%[output:14169921]
%   data: {"dataType":"symbolic","outputData":{"name":"r_c","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{v_x \\,\\sigma_1 }{-c^2 +{v_y }^2 +{v_x }^2 } & \\frac{v_y \\,\\sigma_1 }{-c^2 +{v_y }^2 +{v_x }^2 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =v_x \\,x-\\sqrt{c^2 \\,y^2 -{v_x }^2 \\,y^2 +2\\,v_x \\,v_y \\,x\\,y+c^2 \\,x^2 -{v_y }^2 \\,x^2 }+v_y \\,y\n\\end{array}"}}
%---
%[output:0ff0c22f]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{v_x \\,x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 } & \\frac{v_y \\,x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 } & 0\n\\end{array}\\right)"}}
%---
%[output:9c15096d]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & \\frac{v_y \\,x}{\\sqrt{c^2 -{v_y }^2 }} & 0\n\\end{array}\\right)"}}
%---
%[output:5a8ed240]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\nx & -\\frac{v_y \\,x}{\\sqrt{c^2 -{v_y }^2 }} & 0\n\\end{array}\\right)"}}
%---
%[output:73749778]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{q\\,{{\\left(c^2 -{v_y }^2 \\right)}}^{3\/2} }{4\\,c^3 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q\\,v_y \\,{\\left(c^2 -{v_y }^2 \\right)}}{4\\,c^3 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:9963bf93]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & 0 & -\\frac{q\\,v_y \\,{{\\left(c^2 -{v_y }^2 \\right)}}^{3\/2} }{4\\,c^5 \\,\\varepsilon_0 \\,x^2 \\,\\pi }\n\\end{array}\\right)"}}
%---
%[output:715ae53c]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi }-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(-c^3 +c^2 \\,{\\textrm{vp}}_x \\right)}}{4\\,c^4 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:6a6a9ea7]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:363b55de]
%   data: {"dataType":"symbolic","outputData":{"name":"F_static","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:51db289c]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_neutral","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:40305650]
%   data: {"dataType":"symbolic","outputData":{"name":"ft_magnetic","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:213e5620]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:876f187a]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:0758925b]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:169b257a]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n0 & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:102f85cd]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:4edeb801]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n0 & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:2eee45e2]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:94668975]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:0d82eeb8]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:40259810]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:17072a30]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:4f41e5bd]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:5087ec0b]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-2\\,{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:7c96db09]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-2\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:93ccdba8]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\textrm{vp}}_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:7b1fe95b]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:71f0c583]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n0 & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:5a873ac6]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:408a2ee3]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\textrm{vp}}_x }{2\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:39e2d897]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:87855c9a]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:27b05733]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:2d3f3cfc]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-3\\,{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:8b087497]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-3\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:12513be8]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:432de499]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:0c1642a5]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c+{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:2a7f71ff]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c+{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:2dbfede1]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:44e46089]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(-4\\,c^2 +3\\,v_y \\,{\\textrm{vp}}_y \\right)}}{8\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:30162d12]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{2\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:886607cb]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{2\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:8e82a98e]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{2\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:117e9327]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n0 & -\\frac{q^2 \\,v_y }{2\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:0640159d]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n0 & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:4ee36d78]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:5d553da2]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\textrm{vp}}_x }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:8a0065ff]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n0 & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:5d3a5a83]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:1648a7d3]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:9566a7c8]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\textrm{vp}}_x }{2\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:687d4529]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n0 & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:27a1f563]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:394b0f03]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:08f088c5]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n0 & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:73d49389]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:8d68895c]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:76ccc2f1]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:9b323e14]
%   data: {"dataType":"symbolic","outputData":{"name":"E","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:2f218cbb]
%   data: {"dataType":"symbolic","outputData":{"name":"M","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:4ee7b89f]
%   data: {"dataType":"symbolic","outputData":{"name":"eforce","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:09e56986]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{q^2 \\,{\\left(c^2 -{v_y }^2 \\right)}\\,{\\left(c^2 -{v_y }^2 +{{\\textrm{vp}}_y }^2 \\right)}}{\\sigma_1 } & -\\frac{q^2 \\,v_y \\,\\sqrt{c^2 -{v_y }^2 }\\,{\\left(c^2 -{v_y }^2 +{{\\textrm{vp}}_y }^2 \\right)}}{\\sigma_1 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi \\,{\\left(c^2 -{v_y }^2 +v_y \\,{\\textrm{vp}}_y \\right)}\n\\end{array}"}}
%---
%[output:820d5b48]
%   data: {"dataType":"symbolic","outputData":{"name":"F2T_net","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y +c\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:7dbeb18f]
%   data: {"dataType":"symbolic","outputData":{"name":"F2T_m","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:2a7d5919]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_neutral","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:8a477894]
%   data: {"dataType":"symbolic","outputData":{"name":"ft_magnetic","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
