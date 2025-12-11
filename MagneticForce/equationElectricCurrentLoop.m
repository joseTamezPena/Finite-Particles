%[text] # Current on a loop
% The variables
syms r v vp c a real positive;
assume(r>0)
assume(vp>0)
assume(v>0)
assume(c>0)
assume(a>0)

%[text] ## The fixed protons on the wire
%%
C = -[cos(a), sin(a),0] % force particle arriving the center of the loop %[output:12ff9cd3]
vec_p = [vp,0,0]  % speed of particle at the center of the loop %[output:259580f3]
vec_r = C - vec_p  % relative speed of force particle %[output:6f40bd2f]
vec_rn = simplify(vec_r/sqrt(dot(vec_r,vec_r))) %[output:7c482a92]
dotCVp = simplify(dot(C,vec_rn)) %[output:66a38352]
cos2dot = simplify(dotCVp^2) %[output:32a410b5]
%%
%[text] ## Integrate the loop: X
cos2dotx_p= simplify(expand(-cos2dot*cos(a))) %[output:87a8fc5c]
allloopx = simplify(int(cos2dotx_p,a,IgnoreAnalyticConstraints=true,IgnoreSpecialCases=true),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:58259d11]
netatXul = simplify(subs(allloopx,a,2*pi)) %[output:90792206]
netatXll = subs(allloopx,a,0) %[output:2a4a8c54]
netXp = netatXul - netatXll %[output:6060a364]
%%
%[text] ## Integrate the loop: Y
cos2doty_p= simplify(expand(-cos2dot*sin(a)),'Steps',100) %[output:82fdec6a]
allloopy = simplify(int(cos2doty_p,a,IgnoreAnalyticConstraints=true,IgnoreSpecialCases=true),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:61027145]
netatYul = simplify(subs(allloopy,a,2*pi)) %[output:2fd429f9]
netatYll = simplify(subs(allloopy,a,0)) %[output:19052f23]
netYp = netatYul - netatYll %[output:9f72009d]
%%
%[text] ## Adding the current of the electrons of the loop 
%[text] The current goes counter clock wise.
%[text] The electrons are moving clock wise
vloop = v*[sin(a),-cos(a),0] %[output:01eb994f]
RotA = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,0]' %[output:58a202ae]
uv = [0,1,0]*RotA %[output:50be4e02]
C= [(1-v^2)^(1/2),v,0]*RotA %[output:555e49b4]
vec_r = C + vloop - vec_p  % relative speed of force particle %[output:10db4e1a]
vec_rn = simplify(vec_r/sqrt(dot(vec_r,vec_r)),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:5d9ec1de]
dotCVp = simplify(dot(C,vec_rn),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:50a991bb]

syms v real positive
assume(v>0)
assume(v^2<1)
cosdott = simplify(taylor(dotCVp,v,ExpansionPoint=0,Order=2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:07adacbe]

syms v real positive
assume(v > 0)
cos2dot = simplify(cosdott^2,'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:943cd90a]

assume(v>0)
cos2dott = simplify(taylor(cos2dot,vp,ExpansionPoint=0,Order=2),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:4ca31f5d]
%%
%[text] ## Integrating the current: X
syms r v vp c a real positive;
assume(r>0)
assume(vp>0)
assume(v>0)
assume(c>0)
assume(a>0)

cos2dotx= simplify(cos2dott*C(1),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:8aa3ec10]

allloopxI = simplify(int(cos2dotx,a,IgnoreAnalyticConstraints=true,IgnoreSpecialCases=true),'Criterion','preferReal','Steps',100) %[output:76575491]
%%
assume(vp>0)
netatxul = simplify(subs(allloopx,a,2*pi)) %[output:56d378c0]
netatxll = simplify(subs(allloopx,a,0)) %[output:4f14e003]
netX = simplify(netXp - (1-v^2)^(1/2)*(netatxul - netatxll)) %[output:1500c8ee]
%%
%[text] ## Integrating the current: Y
cos2doty= simplify(cos2dott*C(2),'Steps',100) %[output:7c821ac0]
allloopy = simplify(int(cos2doty,a,IgnoreAnalyticConstraints=true,IgnoreSpecialCases=true),'Criterion','preferReal','Steps',100) %[output:6e48af87]
netatYul = subs(allloopy,a,2*pi) %[output:034f9dc6]
netatYll = subs(allloopy,a,0) %[output:18875016]
netY = simplify(netYp - (1-v^2)^(1/2)*(netatYul-netatYll)) %[output:4328e22d]
%%
%[text] ## The total force of a moving particle in the center of the current loop
syms I Q_1 Q_2 R v e_o v_1 m_o real positive;
totF = Q_1*Q_2/(4*pi*e_o*R)*[netX,netY,0] % The total  %[output:34529a40]
totFT = taylor(totF,v,ExpansionPoint=0,Order=2) %[output:0df8079a]
v = I/(c*Q_1) %[output:067ed5e0]
e_o=1/(c^2*m_o) %[output:3cb9009c]
vp = v_1/c %[output:52bda58d]
totFT= simplify(subs(totFT),'IgnoreAnalyticConstraints',true,'Criterion','preferReal','Steps',100) %[output:64c59c17]


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":22}
%---
%[output:12ff9cd3]
%   data: {"dataType":"symbolic","outputData":{"name":"C","value":"\\left(\\begin{array}{ccc}\n-\\cos \\left(a\\right) & -\\sin \\left(a\\right) & 0\n\\end{array}\\right)"}}
%---
%[output:259580f3]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_p","value":"\\left(\\begin{array}{ccc}\n\\mathrm{vp} & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:6f40bd2f]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_r","value":"\\left(\\begin{array}{ccc}\n-\\mathrm{vp}-\\cos \\left(a\\right) & -\\sin \\left(a\\right) & 0\n\\end{array}\\right)"}}
%---
%[output:7c482a92]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_rn","value":"\\left(\\begin{array}{ccc}\n-\\frac{\\mathrm{vp}+\\cos \\left(a\\right)}{\\sqrt{2\\,\\mathrm{vp}\\,\\cos \\left(a\\right)+{\\mathrm{vp}}^2 +1}} & -\\frac{\\sin \\left(a\\right)}{\\sqrt{2\\,\\mathrm{vp}\\,\\cos \\left(a\\right)+{\\mathrm{vp}}^2 +1}} & 0\n\\end{array}\\right)"}}
%---
%[output:66a38352]
%   data: {"dataType":"symbolic","outputData":{"name":"dotCVp","value":"\\frac{\\mathrm{vp}\\,\\cos \\left(a\\right)+1}{\\sqrt{2\\,\\mathrm{vp}\\,\\cos \\left(a\\right)+{\\mathrm{vp}}^2 +1}}"}}
%---
%[output:32a410b5]
%   data: {"dataType":"symbolic","outputData":{"name":"cos2dot","value":"\\frac{{{\\left(\\mathrm{vp}\\,\\cos \\left(a\\right)+1\\right)}}^2 }{2\\,\\mathrm{vp}\\,\\cos \\left(a\\right)+{\\mathrm{vp}}^2 +1}"}}
%---
%[output:87a8fc5c]
%   data: {"dataType":"symbolic","outputData":{"name":"cos2dotx_p","value":"-\\frac{\\cos \\left(a\\right)\\,{{\\left(\\mathrm{vp}\\,\\cos \\left(a\\right)+1\\right)}}^2 }{2\\,\\mathrm{vp}\\,\\cos \\left(a\\right)+{\\mathrm{vp}}^2 +1}"}}
%---
%[output:58259d11]
%   data: {"dataType":"symbolic","outputData":{"name":"allloopx","value":"\\begin{array}{l}\n\\frac{{\\mathrm{vp}}^2 \\,\\sin \\left(a\\right)}{4}-\\frac{\\mathrm{vp}\\,\\sin \\left(2\\,a\\right)}{8}-\\frac{3\\,\\sin \\left(a\\right)}{4}-\\frac{\\frac{a}{8}+\\sigma_1 }{\\mathrm{vp}}-{\\mathrm{vp}}^3 \\,{\\left(\\frac{a}{8}-\\sigma_1 \\right)}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{\\mathrm{atan}\\left(\\frac{\\sin \\left(\\frac{a}{2}\\right)\\,{\\left(-1+\\mathrm{vp}\\right)}}{\\cos \\left(\\frac{a}{2}\\right)\\,{\\left(1+\\mathrm{vp}\\right)}}\\right)}{4}\n\\end{array}"}}
%---
%[output:90792206]
%   data: {"dataType":"symbolic","outputData":{"name":"netatXul","value":"-\\frac{\\pi \\,{\\left(1+{\\mathrm{vp}}^4 \\right)}}{4\\,\\mathrm{vp}}"}}
%---
%[output:2a4a8c54]
%   data: {"dataType":"symbolic","outputData":{"name":"netatXll","value":"0"}}
%---
%[output:6060a364]
%   data: {"dataType":"symbolic","outputData":{"name":"netXp","value":"-\\frac{\\pi \\,{\\left(1+{\\mathrm{vp}}^4 \\right)}}{4\\,\\mathrm{vp}}"}}
%---
%[output:82fdec6a]
%   data: {"dataType":"symbolic","outputData":{"name":"cos2doty_p","value":"-\\frac{\\sin \\left(a\\right)\\,{{\\left(\\mathrm{vp}\\,\\cos \\left(a\\right)+1\\right)}}^2 }{2\\,\\mathrm{vp}\\,\\cos \\left(a\\right)+{\\mathrm{vp}}^2 +1}"}}
%---
%[output:61027145]
%   data: {"dataType":"symbolic","outputData":{"name":"allloopy","value":"\\frac{\\mathrm{vp}\\,{\\cos \\left(a\\right)}^2 }{4}-\\cos \\left(a\\right)\\,{\\left(-\\frac{3}{4}+\\frac{{\\mathrm{vp}}^2 }{4}\\right)}+\\frac{\\log \\left(2\\,\\mathrm{vp}\\,\\cos \\left(a\\right)+{\\mathrm{vp}}^2 +1\\right)\\,{{\\left(-1+{\\mathrm{vp}}^2 \\right)}}^2 }{8\\,\\mathrm{vp}}"}}
%---
%[output:2fd429f9]
%   data: {"dataType":"symbolic","outputData":{"name":"netatYul","value":"\\frac{\\mathrm{vp}}{4}-\\frac{{\\mathrm{vp}}^2 }{4}+\\frac{\\log \\left(1+\\mathrm{vp}\\right)\\,{{\\left(-1+{\\mathrm{vp}}^2 \\right)}}^2 }{4\\,\\mathrm{vp}}+\\frac{3}{4}"}}
%---
%[output:19052f23]
%   data: {"dataType":"symbolic","outputData":{"name":"netatYll","value":"\\frac{\\mathrm{vp}}{4}-\\frac{{\\mathrm{vp}}^2 }{4}+\\frac{\\log \\left(1+\\mathrm{vp}\\right)\\,{{\\left(-1+{\\mathrm{vp}}^2 \\right)}}^2 }{4\\,\\mathrm{vp}}+\\frac{3}{4}"}}
%---
%[output:9f72009d]
%   data: {"dataType":"symbolic","outputData":{"name":"netYp","value":"0"}}
%---
%[output:01eb994f]
%   data: {"dataType":"symbolic","outputData":{"name":"vloop","value":"\\left(\\begin{array}{ccc}\nv\\,\\sin \\left(a\\right) & -v\\,\\cos \\left(a\\right) & 0\n\\end{array}\\right)"}}
%---
%[output:58a202ae]
%   data: {"dataType":"symbolic","outputData":{"name":"RotA","value":"\\left(\\begin{array}{ccc}\n\\cos \\left(a\\right) & \\sin \\left(a\\right) & 0\\\\\n-\\sin \\left(a\\right) & \\cos \\left(a\\right) & 0\\\\\n0 & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:50be4e02]
%   data: {"dataType":"symbolic","outputData":{"name":"uv","value":"\\left(\\begin{array}{ccc}\n-\\sin \\left(a\\right) & \\cos \\left(a\\right) & 0\n\\end{array}\\right)"}}
%---
%[output:555e49b4]
%   data: {"dataType":"symbolic","outputData":{"name":"C","value":"\\left(\\begin{array}{ccc}\n\\cos \\left(a\\right)\\,\\sqrt{1-v^2 }-v\\,\\sin \\left(a\\right) & \\sin \\left(a\\right)\\,\\sqrt{1-v^2 }+v\\,\\cos \\left(a\\right) & 0\n\\end{array}\\right)"}}
%---
%[output:10db4e1a]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_r","value":"\\left(\\begin{array}{ccc}\n\\cos \\left(a\\right)\\,\\sqrt{1-v^2 }-\\mathrm{vp} & \\sin \\left(a\\right)\\,\\sqrt{1-v^2 } & 0\n\\end{array}\\right)"}}
%---
%[output:5d9ec1de]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_rn","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n-\\frac{{{\\left(1-v^2 \\right)}}^{1\/4} \\,{\\left(\\mathrm{vp}-\\cos \\left(a\\right)\\,\\sqrt{1-v^2 }\\right)}}{\\sigma_1 } & \\frac{\\sin \\left(a\\right)\\,{{\\left(1-v^2 \\right)}}^{3\/4} }{\\sigma_1 } & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\sqrt{{\\mathrm{vp}}^2 \\,\\sqrt{1-v^2 }+\\left|1-v^2 \\right|\\,\\sqrt{1-v^2 }-\\mathrm{vp}\\,\\cos \\left(a\\right)+v^2 \\,\\mathrm{vp}\\,\\cos \\left(a\\right)-\\mathrm{vp}\\,\\cos \\left(a\\right)\\,\\left|1-v^2 \\right|}\n\\end{array}"}}
%---
%[output:50a991bb]
%   data: {"dataType":"symbolic","outputData":{"name":"dotCVp","value":"\\begin{array}{l}\n\\frac{\\sigma_2 -\\sigma_1 +v\\,\\mathrm{vp}\\,\\sin \\left(a\\right)\\,\\sigma_3 }{{{\\left(1-v^2 \\right)}}^{1\/4} \\,\\sqrt{{\\mathrm{vp}}^2 \\,\\sigma_3 +\\sigma_2 -\\mathrm{vp}\\,\\cos \\left(a\\right)+v^2 \\,\\mathrm{vp}\\,\\cos \\left(a\\right)-\\sigma_1 }}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\mathrm{vp}\\,\\cos \\left(a\\right)\\,\\left|1-v^2 \\right|\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =\\left|1-v^2 \\right|\\,\\sigma_3 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =\\sqrt{1-v^2 }\n\\end{array}"}}
%---
%[output:07adacbe]
%   data: {"dataType":"symbolic","outputData":{"name":"cosdott","value":"\\frac{v\\,\\mathrm{vp}\\,\\sin \\left(a\\right)-\\mathrm{vp}\\,\\cos \\left(a\\right)+1}{\\sqrt{{\\mathrm{vp}}^2 -2\\,\\mathrm{vp}\\,\\cos \\left(a\\right)+1}}"}}
%---
%[output:943cd90a]
%   data: {"dataType":"symbolic","outputData":{"name":"cos2dot","value":"\\frac{{{\\left(v\\,\\mathrm{vp}\\,\\sin \\left(a\\right)-\\mathrm{vp}\\,\\cos \\left(a\\right)+1\\right)}}^2 }{{\\mathrm{vp}}^2 -2\\,\\mathrm{vp}\\,\\cos \\left(a\\right)+1}"}}
%---
%[output:4ca31f5d]
%   data: {"dataType":"symbolic","outputData":{"name":"cos2dott","value":"2\\,v\\,\\mathrm{vp}\\,\\sin \\left(a\\right)+1"}}
%---
%[output:8aa3ec10]
%   data: {"dataType":"symbolic","outputData":{"name":"cos2dotx","value":"{\\left(\\cos \\left(a\\right)\\,\\sqrt{1-v^2 }-v\\,\\sin \\left(a\\right)\\right)}\\,{\\left(2\\,v\\,\\mathrm{vp}\\,\\sin \\left(a\\right)+1\\right)}"}}
%---
%[output:76575491]
%   data: {"dataType":"symbolic","outputData":{"name":"allloopxI","value":"\\sin \\left(a\\right)\\,\\sqrt{1-v^2 }+v\\,\\cos \\left(a\\right)-a\\,v^2 \\,\\mathrm{vp}+v^2 \\,\\mathrm{vp}\\,\\cos \\left(a\\right)\\,\\sin \\left(a\\right)-v\\,\\mathrm{vp}\\,{\\cos \\left(a\\right)}^2 \\,\\sqrt{1-v^2 }"}}
%---
%[output:56d378c0]
%   data: {"dataType":"symbolic","outputData":{"name":"netatxul","value":"-\\frac{\\pi \\,{\\left(1+{\\mathrm{vp}}^4 \\right)}}{4\\,\\mathrm{vp}}"}}
%---
%[output:4f14e003]
%   data: {"dataType":"symbolic","outputData":{"name":"netatxll","value":"0"}}
%---
%[output:1500c8ee]
%   data: {"dataType":"symbolic","outputData":{"name":"netX","value":"\\frac{\\pi \\,{\\left(1+{\\mathrm{vp}}^4 \\right)}\\,{\\left(\\sqrt{1-v^2 }-1\\right)}}{4\\,\\mathrm{vp}}"}}
%---
%[output:7c821ac0]
%   data: {"dataType":"symbolic","outputData":{"name":"cos2doty","value":"{\\left(\\sin \\left(a\\right)\\,\\sqrt{1-v^2 }+v\\,\\cos \\left(a\\right)\\right)}\\,{\\left(2\\,v\\,\\mathrm{vp}\\,\\sin \\left(a\\right)+1\\right)}"}}
%---
%[output:6e48af87]
%   data: {"dataType":"symbolic","outputData":{"name":"allloopy","value":"v\\,\\sin \\left(a\\right)-\\cos \\left(a\\right)\\,\\sqrt{1-v^2 }-v^2 \\,\\mathrm{vp}\\,{\\cos \\left(a\\right)}^2 +a\\,v\\,\\mathrm{vp}\\,\\sqrt{1-v^2 }-v\\,\\mathrm{vp}\\,\\cos \\left(a\\right)\\,\\sin \\left(a\\right)\\,\\sqrt{1-v^2 }"}}
%---
%[output:034f9dc6]
%   data: {"dataType":"symbolic","outputData":{"name":"netatYul","value":"2\\,\\pi \\,v\\,\\mathrm{vp}\\,\\sqrt{1-v^2 }-\\sqrt{1-v^2 }-v^2 \\,\\mathrm{vp}"}}
%---
%[output:18875016]
%   data: {"dataType":"symbolic","outputData":{"name":"netatYll","value":"-v^2 \\,\\mathrm{vp}-\\sqrt{1-v^2 }"}}
%---
%[output:4328e22d]
%   data: {"dataType":"symbolic","outputData":{"name":"netY","value":"2\\,\\pi \\,v\\,\\mathrm{vp}\\,{\\left(-1+v^2 \\right)}"}}
%---
%[output:34529a40]
%   data: {"dataType":"symbolic","outputData":{"name":"totF","value":"\\left(\\begin{array}{ccc}\n\\frac{Q_1 \\,Q_2 \\,{\\left(1+{\\mathrm{vp}}^4 \\right)}\\,{\\left(\\sqrt{1-v^2 }-1\\right)}}{16\\,R\\,e_o \\,\\mathrm{vp}} & \\frac{Q_1 \\,Q_2 \\,v\\,\\mathrm{vp}\\,{\\left(-1+v^2 \\right)}}{2\\,R\\,e_o } & 0\n\\end{array}\\right)"}}
%---
%[output:0df8079a]
%   data: {"dataType":"symbolic","outputData":{"name":"totFT","value":"\\left(\\begin{array}{ccc}\n0 & -\\frac{Q_1 \\,Q_2 \\,v\\,\\mathrm{vp}}{2\\,R\\,e_o } & 0\n\\end{array}\\right)"}}
%---
%[output:067ed5e0]
%   data: {"dataType":"symbolic","outputData":{"name":"v","value":"\\frac{\\textrm{I}}{Q_1 \\,c}"}}
%---
%[output:3cb9009c]
%   data: {"dataType":"symbolic","outputData":{"name":"e_o","value":"\\frac{1}{c^2 \\,m_o }"}}
%---
%[output:52bda58d]
%   data: {"dataType":"symbolic","outputData":{"name":"vp","value":"\\frac{v_1 }{c}"}}
%---
%[output:64c59c17]
%   data: {"dataType":"symbolic","outputData":{"name":"totFT","value":"\\left(\\begin{array}{ccc}\n0 & -\\frac{\\textrm{I}\\,Q_2 \\,m_o \\,v_1 }{2\\,R} & 0\n\\end{array}\\right)"}}
%---
