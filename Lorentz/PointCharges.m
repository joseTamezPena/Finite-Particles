%[text] # Electric and magnetic field due to a charge moving at constat velocity

syms  q v_x v_y v_z c_x c_y c_z x y z real positive
syms t epsilon0 mu0 real positive


% The charge is at zero at t=0
r_o = [0,0,0];
% Define the position of the space where the field will be measured
r_p = [x, y, z];

% Define the velocity vector of the moving charge
v = [v_x, v_y, v_z];

% Light vector
cv = [c_x, c_y, c_z];
%%
%[text] ## Compute the time to travel
syms R C V real positive
syms vp_x vp_y cp_x cp_y real positive
assumeAlso (vp_x^2 + vp_y^2 < C)
assumeAlso (t > 0)
assumeAlso (R > 0)
assumeAlso (C > 0)

vp = [vp_x,vp_y];
cp = [cp_x,cp_y];

% the two reach the same point
eq_1=(vp+cp)*t==R %[output:0e96de90]
eq_2 = dot(cp,cp)==C %[output:4603ef49]


[ccx,ccy,tt] = solve([eq_1,eq_2],[cp_x,cp_y,t],"Real",true,"IgnoreAnalyticConstraints",true) %[output:732ed5e9] %[output:27ec514a] %[output:7b67f273] %[output:27dd32f0]
ccx = simplify(subs(ccx(2),vp_y^2,V-vp_x^2)) %[output:75881d02]
ccy = simplify(subs(ccy(2),vp_y^2,V-vp_x^2)) %[output:92cc4795]
tt = simplify(subs(tt(2),vp_y^2,V-vp_x^2)) %[output:531b8898]
%%

%%


%%
r_c = cv*t
% Electric field E due to a moving charge
E = (1/(4*pi*epsilon0)) * (q * (r_p - r_c)) / norm(r_p - r_c)^3;

% Magnetic field B due to a moving charge
B = simplify((1.0/(c^2)) * cross(v,E));

% Display the equations
E, B
%%
%[text] ## The Lorenz force
syms vp_x vp_y vp_z real
vp = [vp_x,vp_y,vp_z];
% Define the Lorentz force F acting on a second charge with v2 velocity
F = simplify(q * (E + cross(vp, B)),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)
%%
%[text] ## The corpuscle force
%[text] 
%[text] \\$\\$ 
%[text] \\mathbf{F} = \\frac{k q\_1 q\_2}{4 \\pi r^2} \\frac{\\|\\mathbf{c} + \\mathbf{v\_1}\\|}{\\|c\\|}  \\left( \\frac{\\mathbf{c} \\cdot ( \\mathbf{c} + \\mathbf{v\_1} - \\mathbf{v\_2} )}{\\|\\mathbf{c}\\| \\|\\mathbf{c} + \\mathbf{v\_1} - \\mathbf{v\_2}\\|} \\right)^2 \\hat{o\_1}, 
%[text] \\$\\$
%[text] 
R2 = dot(r_p,r_p);
cv = r_p/t-v;
ov = cv/c;
vs = cv+v;
rvel= vs-vp;
vsmag = norm(vs);
rvelm = norm(rvel);
F_c = simplify(q^2/(4*pi*epsilon0*R2)*(vsmag/c)*(dot(ov,rvel/rvelm)^2)*ov)
Ft_c = taylor(F_c,[v_x,v_y,v_z],[0,0,0],order=2);
Ft_c = simplify(Ft_c,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100)


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:0e96de90]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_1","value":"\\left(\\begin{array}{cc}\nt\\,{\\left({\\textrm{cp}}_x +{\\textrm{vp}}_x \\right)}=R & t\\,{\\left({\\textrm{cp}}_y +{\\textrm{vp}}_y \\right)}=R\n\\end{array}\\right)"}}
%---
%[output:4603ef49]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_2","value":"{{\\textrm{cp}}_y }^2 +{{\\textrm{cp}}_x }^2 =C"}}
%---
%[output:732ed5e9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:27ec514a]
%   data: {"dataType":"symbolic","outputData":{"name":"ccx","value":"\\left(\\begin{array}{c}\n{\\textrm{vp}}_y -\\frac{R\\,{\\textrm{vp}}_x +R\\,{\\textrm{vp}}_y +R\\,\\sqrt{2\\,C-{{\\textrm{vp}}_y }^2 +2\\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y -{{\\textrm{vp}}_x }^2 }}{2\\,R}\\\\\n{\\textrm{vp}}_y -\\frac{R\\,{\\textrm{vp}}_x +R\\,{\\textrm{vp}}_y -R\\,\\sqrt{2\\,C-{{\\textrm{vp}}_y }^2 +2\\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y -{{\\textrm{vp}}_x }^2 }}{2\\,R}\n\\end{array}\\right)"}}
%---
%[output:7b67f273]
%   data: {"dataType":"symbolic","outputData":{"name":"ccy","value":"\\left(\\begin{array}{c}\n{\\textrm{vp}}_x -\\frac{R\\,{\\textrm{vp}}_x +R\\,{\\textrm{vp}}_y +R\\,\\sqrt{2\\,C-{{\\textrm{vp}}_y }^2 +2\\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y -{{\\textrm{vp}}_x }^2 }}{2\\,R}\\\\\n{\\textrm{vp}}_x -\\frac{R\\,{\\textrm{vp}}_x +R\\,{\\textrm{vp}}_y -R\\,\\sqrt{2\\,C-{{\\textrm{vp}}_y }^2 +2\\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y -{{\\textrm{vp}}_x }^2 }}{2\\,R}\n\\end{array}\\right)"}}
%---
%[output:27dd32f0]
%   data: {"dataType":"symbolic","outputData":{"name":"tt","value":"\\left(\\begin{array}{c}\n\\frac{R\\,{\\textrm{vp}}_x +R\\,{\\textrm{vp}}_y +R\\,\\sqrt{2\\,C-{{\\textrm{vp}}_y }^2 +2\\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y -{{\\textrm{vp}}_x }^2 }}{-C+{{\\textrm{vp}}_y }^2 +{{\\textrm{vp}}_x }^2 }\\\\\n\\frac{R\\,{\\textrm{vp}}_x +R\\,{\\textrm{vp}}_y -R\\,\\sqrt{2\\,C-{{\\textrm{vp}}_y }^2 +2\\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y -{{\\textrm{vp}}_x }^2 }}{-C+{{\\textrm{vp}}_y }^2 +{{\\textrm{vp}}_x }^2 }\n\\end{array}\\right)"}}
%---
%[output:75881d02]
%   data: {"dataType":"symbolic","outputData":{"name":"ccx","value":"\\frac{{\\textrm{vp}}_y }{2}-\\frac{{\\textrm{vp}}_x }{2}+\\frac{\\sqrt{2\\,C-V+2\\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y }}{2}"}}
%---
%[output:92cc4795]
%   data: {"dataType":"symbolic","outputData":{"name":"ccy","value":"\\frac{{\\textrm{vp}}_x }{2}-\\frac{{\\textrm{vp}}_y }{2}+\\frac{\\sqrt{2\\,C-V+2\\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y }}{2}"}}
%---
%[output:531b8898]
%   data: {"dataType":"symbolic","outputData":{"name":"tt","value":"-\\frac{R\\,{\\left({\\textrm{vp}}_x +{\\textrm{vp}}_y -\\sqrt{2\\,C-V+2\\,{\\textrm{vp}}_x \\,{\\textrm{vp}}_y }\\right)}}{C-V}"}}
%---
