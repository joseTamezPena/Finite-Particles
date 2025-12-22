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

Ft = taylor(Fx,v_y,0,order=2) %[output:162765a8]

%Ft = Fx;
Ft = taylor(Ft,[vp_x,vp_y],[0,0],order=2);

Ft = simplify(Ft,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:3d514c36]
F_static = subs(Ft,v_y,0) %[output:4f74780b]
Ft_neutral = simplify((Ft - F_static),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:62a1ecf1]
ft_magnetic = simplify(q * cross(vp, B),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = taylor(ft_magnetic,[v_x,v_y],[0,0],order=2);
ft_magnetic = simplify(taylor(ft_magnetic,[vp_x,vp_y],[0,0],order=2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
ft_magnetic = simplify(subs(ft_magnetic,[v_x,y],[0,0])) %[output:0882f2c3]
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
assumeAlso(c^2*y^2 - v_x^2*y^2 + c^2*x^2 >0)
R2 = dot(r_p,r_p);

ov = simplify(cvec/sqrt(dot(cvec,cvec)));


cv2 = cvec;

vs = cv2+v;
rvel= vs-vp;

p4=simplify(subs(ov,[v_x,y],[0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

vsmag = dot(vs,vs);
rvelm = dot(rvel,rvel);
pvec2 = dot(pvec,pvec);

eforce = simplify(subs(k/pvec2*p4*(dot(p4,rvel/sqrt(rvelm))^2),[v_x,v_y,y],[0,0,0]));
eforce = taylor(eforce,[vp_x,vp_y],[0,0],order=2) %[output:1d20e58e]

p1=simplify(subs((sqrt(vsmag)/c),[v_x,y],[0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
p2=simplify(subs((dot(p4,rvel/sqrt(rvelm))^2),[v_x,y],[0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
p3=simplify(subs((k/R2),[v_x,y],[0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

F_c = simplify(p1*p2*p3*p4,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

FT_c = taylor(F_c,v_y,0,order=2);
FT_c = taylor(FT_c,[vp_x,vp_y],[0,0],order=2) %[output:4bf374a5]

%%
%% Second Option

%p2=simplify(2-subs(1.0/dot(ov,rvel/sqrt(rvelm))^2,[v_x,y],[0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
dotp = simplify(dot(p4,rvel/sqrt(rvelm)));
%p2=simplify(subs(1.0/dotp,[v_x,y],[0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

p2=simplify(subs(sqrt(rvelm/vsmag)/dotp,[v_x,y],[0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
%p2a=simplify(subs(1.0/dot(vs/sqrt(vsmag),rvel/sqrt(rvelm)),[v_x,y],[0,0]),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);
%p2=simplify(p2*p2a,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)

simplify(taylor(subs(p4,[v_x,y],[0,0]),[vp_x,vp_y],[0,0],order=2)) %[output:0e15a269]
simplify(taylor(subs(rvel,[v_x,y],[0,0]),[vp_x,vp_y],[0,0],order=2)) %[output:0663d88f]
p2s = simplify(taylor(p2,v_y,0,order=2));
simplify(taylor(p2s,[vp_x,vp_y],[0,0],order=2)) %[output:92768fa2]

F2_c = simplify(p1*p2*p3*p4,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

F2T_c = taylor(F2_c,v_y,0,order=2);
F2T_c = simplify(taylor(F2T_c,[vp_x,vp_y],[0,0],order=2)-eforce) %[output:4a5087e7]


%%

Ft_neutral,ft_magnetic %[output:17012985] %[output:059e7f34]

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
%[output:162765a8]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi }-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(-c^3 +c^2 \\,{\\textrm{vp}}_x \\right)}}{4\\,c^4 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:3d514c36]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:4f74780b]
%   data: {"dataType":"symbolic","outputData":{"name":"F_static","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:62a1ecf1]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_neutral","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:0882f2c3]
%   data: {"dataType":"symbolic","outputData":{"name":"ft_magnetic","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:1d20e58e]
%   data: {"dataType":"symbolic","outputData":{"name":"eforce","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:4bf374a5]
%   data: {"dataType":"symbolic","outputData":{"name":"FT_c","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi }+\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y }{4\\,c\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:0e15a269]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{\\sqrt{c^2 -{v_y }^2 }}{c} & -\\frac{v_y }{c} & 0\n\\end{array}\\right)"}}
%---
%[output:0663d88f]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\sqrt{c^2 -{v_y }^2 }-{\\textrm{vp}}_x  & -{\\textrm{vp}}_y  & 0\n\\end{array}\\right)"}}
%---
%[output:92768fa2]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-\\frac{-c^2 +v_y \\,{\\textrm{vp}}_y +c\\,{\\textrm{vp}}_x }{c^2 }"}}
%---
%[output:4a5087e7]
%   data: {"dataType":"symbolic","outputData":{"name":"F2T_c","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(v_y \\,{\\textrm{vp}}_y +c\\,{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:17012985]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_neutral","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:059e7f34]
%   data: {"dataType":"symbolic","outputData":{"name":"ft_magnetic","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\textrm{vp}}_x }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
