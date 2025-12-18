%[text] # Electric and magnetic field due to a charge moving at constant velocity

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


% The charge is at zero at t=0
r_o = [0,0,0];
% Define the position of the space where the field will be measured
r_p = [x, y, 0];

% Define the velocity vector of the moving charge
v = [v_x, v_y, 0];

% Light vector
cv = [c_x, c_y, 0];

% The velocity of the charged particle at r_p

vp = [vp_x,vp_y,0];


%%
%[text] ## Compute the time to travel

% the two reach the same point
eq_1=(v+cv)*t==r_p %[output:49f885b4]
eq_2 = dot(cv,cv)==c^2 %[output:2aeb5593]


[ccx,ccy,tt] = solve([eq_1,eq_2],[c_x,c_y,t],"Real",true,"IgnoreAnalyticConstraints",true); %[output:6053e22f]
ccx = simplify(ccx(2));
ccy = simplify(ccy(2));
cvec = [-ccx,ccy,0];
assumeAlso (c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2 > 0)
simplify(dot(cvec,cvec));

tt = simplify(tt(1));
simplify(subs(ccx,y,0)) %[output:228a2de4]
simplify(subs(ccy,y,0)) %[output:2e1230b4]
simplify(subs(tt,y,0)) %[output:17e89be4]
%%

r_c = simplify(expand(v*tt)) %[output:199f6115]

simplify(subs(r_c,y,0)) %[output:62dd322f]
simplify(subs(r_c,[y,v_x],[0,0])) %[output:788d5e1d]

%%
%[text] ## Electric and Magnetic Field
assumeAlso (v_x*x - sqrt(c^2*y^2 - v_x^2*y^2 + 2*v_x*v_y*x*y + c^2*x^2 - v_y^2*x^2) + v_y*y>0)


% Electric field E due to a moving charge
pvec = r_p - r_c;
simplify(subs(pvec,[y,v_x],[0,0])) %[output:85c11572]

E = (1/(4*pi*epsilon0)) * (q * pvec) / norm(pvec)^3;
E = simplify(E);
% Magnetic field B due to a moving charge
B = simplify((1.0/(c^2)) * cross(v,E));

% Display the equations
simplify(subs(E,[y,v_x],[0,0])), simplify(subs(B,[y,v_x],[0,0])) %[output:51211e6e] %[output:996c96eb]
%%
%[text] ## The Lorenz force

simplify(cross(vp,simplify(subs(B,[y,v_x],[0,0])))) %[output:5187c7fb]

% Define the Lorentz force F acting on a second charge with vp velocity
F = simplify(q * (E - simplify(cross(vp, B))),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100);

Fx = simplify(subs(F,[y,v_x],[0,0])) %[output:21e0d8eb]

Ft = taylor(Fx,v_y,0,order=2) %[output:9f5127b9]
Ft = [Ft(1),-Ft(2),0] %[output:8a90037a]
%Ft = Fx;
Ft = taylor(Ft,[vp_x,vp_y],[0,0],order=2);

Ft = simplify(Ft,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:763f366c]
F_static = subs(Ft,v_y,0) %[output:3ec7afee]
Ft_neutral = simplify((Ft - F_static),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:3424ac65]
%%
%[text] ## The corpuscle force
%[text] 
%[text] $&dollar&;&dollar&; \n\\mathbf{F} = \\frac{k q\_1 q\_2}{4 \\pi r^2} \\frac{\\|\\mathbf{c} + \\mathbf{v\_1}\\|}{\\|c\\|}  \\left( \\frac{\\mathbf{c} \\cdot ( \\mathbf{c} + \\mathbf{v\_1} - \\mathbf{v\_2} )}{\\|\\mathbf{c}\\| \\|\\mathbf{c} + \\mathbf{v\_1} - \\mathbf{v\_2}\\|} \\right)^2 \\hat{o\_1}, \n&dollar&;&dollar&;$
%[text] 
assumeAlso(c^2 - v_y^2>0)
R2 = dot(r_p,r_p);

ov = simplify(pvec/sqrt(dot(pvec,pvec)));
cv2 = c*ov;
simplify(subs(cv2,[v_x,y],[0,0])) %[output:381f1329]
simplify(subs(cvec,[y,v_x],[0,0])) %[output:8da40636]


vs = cv2+v;
subs(vs,[v_x,y],[0,0]) %[output:30ca583a]
rvel= vs-vp;
vsmag = dot(vs,vs);
rvelm = dot(rvel,rvel);

pvec2 = dot(pvec,pvec);

F_c = simplify(q^2/(4*pi*epsilon0*R2*c)*(vsmag/rvelm)*dot(ov,rvel)*ov);

rv2 = cv2 - vp;
rv2m = dot(rv2,rv2);

Fx_c_static = +simplify(q^2/(4*pi*epsilon0*pvec2)*dot(ov,rv2)*(c/rv2m)*ov);
%Fx_c_static = +simplify(q^2/(4*pi*epsilon0*R2)*r_p/sqrt(dot(r_p,r_p)));
F_c = F_c - Fx_c_static;

Fx_c = simplify(subs(F_c,[y,v_x],[0,0]));
Ft_c = taylor(Fx_c,v_y,0,order=2);
Ft_c = taylor(Ft_c,[vp_x,vp_y],[0,0],order=2);
Ft_c = simplify(-Ft_c,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:5f6c51dd]

%%

Ft_neutral %[output:7b029db7]

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
%[output:6053e22f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:228a2de4]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-\\sqrt{c^2 -{v_y }^2 }"}}
%---
%[output:2e1230b4]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-v_y"}}
%---
%[output:17e89be4]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\frac{x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 }"}}
%---
%[output:199f6115]
%   data: {"dataType":"symbolic","outputData":{"name":"r_c","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{v_x \\,\\sigma_1 }{-c^2 +{v_y }^2 +{v_x }^2 } & \\frac{v_y \\,\\sigma_1 }{-c^2 +{v_y }^2 +{v_x }^2 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =v_x \\,x-\\sqrt{c^2 \\,y^2 -{v_x }^2 \\,y^2 +2\\,v_x \\,v_y \\,x\\,y+c^2 \\,x^2 -{v_y }^2 \\,x^2 }+v_y \\,y\n\\end{array}"}}
%---
%[output:62dd322f]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{v_x \\,x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 } & \\frac{v_y \\,x\\,{\\left(-\\sqrt{c^2 -{v_y }^2 }+v_x \\right)}}{-c^2 +{v_y }^2 +{v_x }^2 } & 0\n\\end{array}\\right)"}}
%---
%[output:788d5e1d]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & \\frac{v_y \\,x}{\\sqrt{c^2 -{v_y }^2 }} & 0\n\\end{array}\\right)"}}
%---
%[output:85c11572]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\nx & -\\frac{v_y \\,x}{\\sqrt{c^2 -{v_y }^2 }} & 0\n\\end{array}\\right)"}}
%---
%[output:51211e6e]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{q\\,{{\\left(c^2 -{v_y }^2 \\right)}}^{3\/2} }{4\\,c^3 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q\\,v_y \\,{\\left(c^2 -{v_y }^2 \\right)}}{4\\,c^3 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:996c96eb]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & 0 & -\\frac{q\\,v_y \\,{{\\left(c^2 -{v_y }^2 \\right)}}^{3\/2} }{4\\,c^5 \\,\\varepsilon_0 \\,x^2 \\,\\pi }\n\\end{array}\\right)"}}
%---
%[output:5187c7fb]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n-\\frac{q\\,v_y \\,{\\textrm{vp}}_y \\,{{\\left(c^2 -{v_y }^2 \\right)}}^{3\/2} }{4\\,c^5 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q\\,v_y \\,{\\textrm{vp}}_x \\,{{\\left(c^2 -{v_y }^2 \\right)}}^{3\/2} }{4\\,c^5 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:21e0d8eb]
%   data: {"dataType":"symbolic","outputData":{"name":"Fx","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,\\sigma_1 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y \\right)}}{\\sigma_2 } & \\frac{q^2 \\,v_y \\,{\\left({\\textrm{vp}}_x \\,\\sigma_1 -c^4 +c^2 \\,{v_y }^2 \\right)}}{\\sigma_2 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 ={{\\left(c^2 -{v_y }^2 \\right)}}^{3\/2} \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =4\\,c^5 \\,\\varepsilon_0 \\,x^2 \\,\\pi \n\\end{array}"}}
%---
%[output:9f5127b9]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi }-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(-c^4 +c^3 \\,{\\textrm{vp}}_x \\right)}}{4\\,c^5 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:8a90037a]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi }-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & -\\frac{q^2 \\,v_y \\,{\\left(-c^4 +c^3 \\,{\\textrm{vp}}_x \\right)}}{4\\,c^5 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:763f366c]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,{\\left(-c^2 +v_y \\,{\\textrm{vp}}_y \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:3ec7afee]
%   data: {"dataType":"symbolic","outputData":{"name":"F_static","value":"\\left(\\begin{array}{ccc}\n\\frac{q^2 }{4\\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:3424ac65]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_neutral","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
%[output:381f1329]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\sqrt{c^2 -{v_y }^2 } & -v_y  & 0\n\\end{array}\\right)"}}
%---
%[output:8da40636]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\sqrt{c^2 -{v_y }^2 } & -v_y  & 0\n\\end{array}\\right)"}}
%---
%[output:30ca583a]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{c\\,{\\left(c^2 \\,x-{v_y }^2 \\,x\\right)}}{\\sigma_1 } & v_y -\\frac{c\\,v_y \\,\\sqrt{c^2 \\,x^2 -{v_y }^2 \\,x^2 }}{\\sigma_1 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\sqrt{{{\\left(c^2 \\,x-{v_y }^2 \\,x\\right)}}^2 -{v_y }^2 \\,{\\left(-c^2 \\,x^2 +{v_y }^2 \\,x^2 \\right)}}\n\\end{array}"}}
%---
%[output:5f6c51dd]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_c","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{2\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:7b029db7]
%   data: {"dataType":"symbolic","outputData":{"name":"Ft_neutral","value":"\\left(\\begin{array}{ccc}\n-\\frac{q^2 \\,v_y \\,{\\textrm{vp}}_y }{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & \\frac{q^2 \\,v_y \\,{\\left(c-{\\textrm{vp}}_x \\right)}}{4\\,c^2 \\,\\varepsilon_0 \\,x^2 \\,\\pi } & 0\n\\end{array}\\right)"}}
%---
