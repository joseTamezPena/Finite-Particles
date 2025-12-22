%[text] # Electric and magnetic field due to a charge moving at constant velocity
clear;
syms  q  real 
syms t epsilon0 mu0 real positive
syms x c v_y c_x c_y real positive
syms  vp_y real positive

assumeAlso (c - v_y > 0)
assumeAlso (c - vp_y > 0)


% The emission charge is at zero at t=0
r_o = [0,0,0];
% Define the position of the space where the field will be measured
r_p = [x, 0, 0];

% Define the velocity vector of the moving charge
v = [0, v_y, 0];

% Light vector
cv = [c_x, c_y, 0];

% The velocity of the reciever charged particle at r_p

vp = [0,vp_y,0];


%%
%[text] ## Compute the time to travel

% the two reach the same point
eq_1=(v+cv)*t==r_p %[output:49f885b4]
eq_2 = dot(cv,cv)==c^2 %[output:2aeb5593]

[ccx,ccy,tt,par,con] = solve([eq_1(1),eq_1(2),eq_2],[c_x,c_y,t],"Real",true,"IgnoreAnalyticConstraints",true,ReturnConditions=true) %[output:23aa21c1]
%%
ccx = simplify(ccx(1)); %[output:36410db9]
ccy = simplify(ccy(1));
cvec = [-ccx,ccy,0];
assumeAlso (c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2 > 0)
simplify(dot(cvec,cvec))

tt = simplify(tt(1));
simplify(subs(ccx,y,0)), simplify(subs(ccy,y,0)), simplify(subs(tt,y,0))
%%

r_c = simplify(expand(v*tt))

simplify(subs(r_c,y,0))
simplify(subs(r_c,[y,v_x],[0,0]))

%%
%[text] ## Electric and Magnetic Field
assumeAlso (v_x*x - sqrt(c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2) + v_y*y>0)


% Electric field E due to a moving charge
pvec = r_p - r_c;
simplify(subs(pvec,[y,v_x],[0,0]))

E = (1/(4*pi*epsilon0)) * (q * pvec) / norm(pvec)^3;
E = simplify(E);
% Magnetic field B due to a moving charge
B = simplify((1.0/(c^2)) * cross(v,E));

% Display the equations
simplify(subs(E,[y,v_x],[0,0])), simplify(subs(B,[y,v_x],[0,0]))
%%
%[text] ## The Lorenz force

simplify(cross(vp,simplify(subs(B,[y,v_x],[0,0]))))

% Define the Lorentz force F acting on a second charge with vp velocity
F = simplify(q * (E + cross(vp, B)),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

Fx = simplify(subs(F,[y,v_x],[0,0]))

Ft = taylor(Fx,v_y,0,order=2)

%Ft = Fx;
Ft = taylor(Ft,[vp_x,vp_y],[0,0],order=2);

Ft = simplify(Ft,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)
F_static = subs(Ft,v_y,0)
Ft_neutral = simplify((Ft - F_static),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)
ft_magnetic = simplify(q * cross(vp, B),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = taylor(ft_magnetic,[v_x,v_y],[0,0],order=2);
ft_magnetic = simplify(taylor(ft_magnetic,[vp_x,vp_y],[0,0],order=2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = simplify(subs(ft_magnetic,[v_x,y],[0,0]))
%%
%[text] ## The corpuscle force
%[text] 
%[text] $&dollar&;&dollar&; \n\\mathbf{F} = \\frac{k q\_1 q\_2}{4 \\pi r^2} \\frac{\\|\\mathbf{c} + \\mathbf{v\_1}\\|}{\\|c\\|}  \\left( \\frac{\\mathbf{c} \\cdot ( \\mathbf{c} + \\mathbf{v\_1} - \\mathbf{v\_2} )}{\\|\\mathbf{c}\\| \\|\\mathbf{c} + \\mathbf{v\_1} - \\mathbf{v\_2}\\|} \\right)^2 \\hat{o\_1}, \n&dollar&;&dollar&;$
%[text] 
k= q^2/(4*pi*epsilon0);
assumeAlso(c^2 - v_y^2>0)
assumeAlso(c^2 - v_x^2>0)
assumeAlso(c - v_x>0)
assumeAlso(c - v_y>0)
R2 = dot(r_p,r_p);

ov = simplify(cvec/sqrt(dot(cvec,cvec)));

cv2 = cvec;

vs = cv2+v;
rvel= vs-vp;

p4=simplify(subs(ov,[v_x,y,vp_x],[0,0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)
simplify(subs(vs,[v_x,y,vp_x],[0,0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)
simplify(subs(rvel,[v_x,y,vp_x],[0,0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)

vsmag = dot(vs,vs);
rvelm = dot(rvel,rvel);
pvec2 = dot(pvec,pvec);
p1=simplify(subs((sqrt(vsmag)/c),[v_x,y,vp_x],[0,0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)
p2=simplify(subs((dot(ov,rvel/sqrt(rvelm))^2),[v_x,y,vp_x],[0,0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)
p3=simplify(subs((k/R2),[v_x,y,vp_x],[0,0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)

F_c = (k/R2)*(sqrt(vsmag)/c)*(dot(ov,rvel/sqrt(rvelm))^2)*ov;
simplify(subs(F_c,[v_x,y,vp_x],[0,0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)
F_c = simplify(p1*p2*p3*p4,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)

F_axis = simplify(subs(F_c,[v_x,y,vp_x],[0,0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
FT_c = taylor(F_axis,v_y,0,order=2);
FT_c = taylor(FT_c,vp_y,0,order=2)

%%


%F_c = simplify(q^2/(4*pi*epsilon0*R2*c)*(vsmag/rvelm)*dot(ov,rvel)*ov);
F_c = simplify(q^2/(4*pi*epsilon0*R2*c)*sqrt(vsmag)/dot(ov,rvel/sqrt(rvelm))^(3/2)*dot(ov,vs/sqrt(vsmag))^(3/2)*ov);

rv2 = cv2 - vp;
rv2m = dot(rv2,rv2);

Fx_c_static = simplify(q^2/(4*pi*epsilon0*pvec2)/(dot(ov,rv2/sqrt(rv2m))^(3/2))*ov);
%Fx_c_static = +simplify(q^2/(4*pi*epsilon0*R2)*r_p/sqrt(dot(r_p,r_p)));
F_c = F_c - Fx_c_static;

Fx_c = simplify(subs(F_c,[y,v_x],[0,0]));
Ft_c = taylor(Fx_c,v_y,0,order=2);
Ft_c = taylor(Ft_c,[vp_x,vp_y],[0,0],order=2);
Ft_c = simplify(Ft_c,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)

%%

Ft_neutral

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:49f885b4]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_1","value":"\\left(\\begin{array}{ccc}\nc_x \\,t=x & t\\,{\\left(c_y +v_y \\right)}=0 & 0=0\n\\end{array}\\right)"}}
%---
%[output:2aeb5593]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_2","value":"{c_y }^2 +{c_x }^2 =c^2"}}
%---
%[output:23aa21c1]
%   data: {"dataType":"text","outputData":{"text":" \nccx =\n \nEmpty sym: 0-by-1\n \n \nccy =\n \nEmpty sym: 0-by-1\n \n \ntt =\n \nEmpty sym: 0-by-1\n \n \npar =\n \nEmpty sym: 1-by-0\n \n \ncon =\n \nEmpty sym: 0-by-1\n \n","truncated":false}}
%---
%[output:36410db9]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Index exceeds array bounds.\n\nError in indexing (line 967)\n                R_tilde = builtin('subsref',L_tilde,Idx);\n                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"}}
%---
