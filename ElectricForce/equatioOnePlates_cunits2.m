%[text] # One Charged Plate and a Point Particle
clear;
% The variables
% a : angle of point in plate
% d : distance from plate to particle
% v : relative particle velocity [0...1]
% v_a: change in velocity due to aceleration v_a = a*Dt = a*R/C
% r : radial distance between axis to point in plate
syms a d v v_a r real positive;
assume(d>0)
assume(r>=0)
assume(v>=0)
assumeAlso (v < 1)
assume(a>=0)
taylorOrder = 11;

% Particle velocity
vp = [v,0,0];
%%
%[text] ## The vectors
%[text] d= distance from the charge to the plate
%[text] r = radius of a circle in the plate
%[text] The speed of light is 1 (c=1)
%[text] 
R = (d^2+r^2)^(1/2) % total distance from a point in the plate to the charge %[output:0c987dd1]

% Get two points per circle. 
% a is the angle
vec_l1 = [d,r*cos(a),r*sin(a)]/R  % velocity vector of light from a point 1 in left positive plate at angle a %[output:7595b9fe]

vec_l2 = [d,-r*cos(a),-r*sin(a)]/R  % velocity vector of light from a point 2 in left positive plate at 180+a %[output:6ebffa30]
r_vec_l1 = (vec_l1 - vp) % left plate relative velocity of 1 particle (v in light speed fraction =v/c) %[output:269a14cd]
r_vec_l2 = (vec_l2 - vp) % left plate relative velocity of 2 particle (v in light speed fraction =v/c) %[output:273c8f2d]

%%
%[text] ## The Anomality 
% 1 particle disparity between photon angle and trayectory


psin = simplify( (dot(r_vec_l1,r_vec_l1)/dot(r_vec_l1,vec_l1)) ) %[output:3d0726da]


%%
%[text] ## The Relative Net Force of the Two Points
% Force of the two point in plate and moving particle
dforce = (psin*vec_l1 + psin*vec_l2 )/(2*R^2); % net Force in the x direction first points.

dforce = simplify(expand(dforce),'Steps',100) % Two charges %[output:69e50e72]
dforce_x = dforce(1) %[output:4bf04a0c]
%dforce_x = simplify(taylor(expand(subs(dforce_x/2,v,v+v_a) + dforce_x/2),v_a,order=2))
%%
%[text] ## Check at v = 0, r =0
df_v_zero = simplify(subs(dforce_x,[v,r],[0,0]),'Steps',100) % Due to the two point charges %[output:7b7123d6]
df_v_one = simplify(subs(dforce_x,[v,r],[v,0]),'Steps',100) % Due to the two point charges %[output:3b43bbb2]
%%
%[text] ## Integrate all the angles
%dforce_xT = expand(taylor(dforce_x,v,Order=taylorOrder));
dforce_xT = dforce_x;

dcircle = int(dforce_xT,a,0,pi); % From zero to pi i.e., the full circle
%%
dcircle = simplify(dcircle,steps=100); % Net force of the two full circles with charge density
df_v_zero = simplify(subs(dcircle,[v,r],[v,0]),'Steps',100) % Net force of the full zero radius circle (point) %[output:34b522cd]
%[text] ## 
%%
%[text] ## Integrate for all values of r
syms R_o real
assumeAlso (R_o>0)
dtotf = r*dcircle/pi  %[output:4739c8cb]
dtotf = simplify(expand(taylor(dtotf,v,Order=taylorOrder))); % Taylor
%dtotf = simplify(subs(dtotf,v,0)) %static
totInt = int(dtotf,r,0,R_o); % Add all the circles
%%
totf = simplify(totInt,'Criterion','preferReal','Steps',100); % The force on a charge due a charge surface density

%%
%[text] ## Set the Plates Radius to Infinity
sympref('PolynomialDisplayStyle', 'ascend');
totfInf = simplify(limit(totf,R_o,Inf)) % Net force by plate charge density %[output:079a3f68]
%totfInf = subs(totfInf,v,-v)

%%
figure(1) %[output:580f0635]
plotrange = [0,0.99];
fplot(totfInf,plotrange) %[output:580f0635]
hold on %[output:580f0635]
relEq = sqrt(1-v^2) %[output:50ce538e]
relEqApx = taylor(relEq,v,Order=taylorOrder) %[output:18f3f6fe]
fplot(relEq,plotrange) %[output:580f0635]
fplot(relEqApx,plotrange) %[output:580f0635]
xlabel('v/c')  %[output:580f0635]
ylabel('ratio')  %[output:580f0635]
legend({'Force Carriers',"Lorentz Mass","Lorentz Apx"},'Location','northeast') %[output:580f0635]
hold off %[output:580f0635]



hold off %[output:580f0635]
%[text] ## 
%%
%[text] ## 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":22}
%---
%[output:0c987dd1]
%   data: {"dataType":"symbolic","outputData":{"name":"R","value":"\\sqrt{r^2 +d^2 }"}}
%---
%[output:7595b9fe]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_l1","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:6ebffa30]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_l2","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:269a14cd]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l1","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }}-v & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:273c8f2d]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l2","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }}-v & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:3d0726da]
%   data: {"dataType":"symbolic","outputData":{"name":"psin","value":"\\frac{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 -2\\,d\\,v\\,\\sqrt{r^2 +d^2 }}{d^2 +r^2 -d\\,v\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:69e50e72]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n-\\frac{d^4 \\,{\\left(v-v^3 \\right)}-d^3 \\,{\\left(\\sigma_1 -\\sigma_1 \\,v^2 \\right)}+r^2 \\,{\\left(d^2 \\,{\\left(v-v^3 \\right)}-d\\,{\\left(\\sigma_1 +\\sigma_1 \\,v^2 \\right)}\\right)}}{{{\\left(r^2 +d^2 \\right)}}^2 \\,{\\left(r^2 +d^2 -d^2 \\,v^2 \\right)}} & 0 & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\sqrt{r^2 +d^2 }\n\\end{array}"}}
%---
%[output:4bf04a0c]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce_x","value":"\\begin{array}{l}\n-\\frac{d^4 \\,{\\left(v-v^3 \\right)}-d^3 \\,{\\left(\\sigma_1 -\\sigma_1 \\,v^2 \\right)}+r^2 \\,{\\left(d^2 \\,{\\left(v-v^3 \\right)}-d\\,{\\left(\\sigma_1 +\\sigma_1 \\,v^2 \\right)}\\right)}}{{{\\left(r^2 +d^2 \\right)}}^2 \\,{\\left(r^2 +d^2 -d^2 \\,v^2 \\right)}}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\sqrt{r^2 +d^2 }\n\\end{array}"}}
%---
%[output:7b7123d6]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\frac{1}{d^2 }"}}
%---
%[output:3b43bbb2]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_one","value":"-\\frac{-1+v}{d^2 }"}}
%---
%[output:34b522cd]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"-\\frac{\\pi \\,{\\left(-1+v\\right)}}{d^2 }"}}
%---
%[output:4739c8cb]
%   data: {"dataType":"symbolic","outputData":{"name":"dtotf","value":"\\begin{array}{l}\n\\frac{r\\,{\\left(\\pi \\,d\\,r^2 \\,{\\left(\\sigma_1 -d\\,v+v^2 \\,\\sigma_1 +d\\,v^3 \\right)}+\\pi \\,d^3 \\,{\\left(-1+v^2 \\right)}\\,{\\left(d\\,v-\\sigma_1 \\right)}\\right)}}{\\pi \\,{{\\left(r^2 +d^2 \\right)}}^2 \\,{\\left(r^2 +d^2 -d^2 \\,v^2 \\right)}}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\sqrt{r^2 +d^2 }\n\\end{array}"}}
%---
%[output:079a3f68]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInf","value":"1-\\frac{v}{2}+\\frac{2\\,v^2 }{3}+\\frac{v^3 }{4}+\\frac{2\\,v^4 }{15}+\\frac{v^5 }{12}+\\frac{2\\,v^6 }{35}+\\frac{v^7 }{24}+\\frac{2\\,v^8 }{63}+\\frac{v^9 }{40}+\\frac{2\\,v^{10} }{99}"}}
%---
%[output:50ce538e]
%   data: {"dataType":"symbolic","outputData":{"name":"relEq","value":"\\sqrt{1-v^2 }"}}
%---
%[output:18f3f6fe]
%   data: {"dataType":"symbolic","outputData":{"name":"relEqApx","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}-\\frac{v^6 }{16}-\\frac{5\\,v^8 }{128}-\\frac{7\\,v^{10} }{256}"}}
%---
%[output:580f0635]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQ1wXtV55584pCBCqC2gGywlIzbYNNPpGsik0gp2TLvp2ttZe3ZCiCzPdIjGk3oXQpUt0YehATsb25KFZ6IGMuOyXsfdjmw3gWlwuruE7gQKUdGmTuKZZVNwaBQjOdkAthNCZAp1dp6rHPno6v24933v1zn3d2c0tqRzz33O7zn3PX895znnvO0Xv\/jFL4QLAhCAAAQgAAEIOETgbQgYh7yFqRCAAAQgAAEIBAQQMHQECEAAAhCAAAScI4CAcc5lGAwBCEAAAhCAAAKGPgABCEAAAhCAgHMEEDDOuQyDIQABCEAAAhBAwNAHIAABCEAAAhBwjgACxjmXYTAEIAABCEAAAggY+gAEIAABCEAAAs4RQMA45zIMhgAEIAABCEAAAUMfgAAEIAABCEDAOQIIGOdchsEQgAAEIAABCCBg6AMQgAAEIAABCDhHAAHjnMswGAIQgAAEIAABBAx9AAIQgAAEIAAB5wggYJxzGQZDAAIQgAAEIICAoQ9AAAIQgAAEIOAcAQSMcy7DYAhAAAIQgAAEEDD0AQhAAAIQgAAEnCOAgHHOZRgMAQhAAAIQgAAChj4AAQhAAAIQgIBzBBAwzrkMgyEAAQhAAAIQQMDQByAAAQhAAAIQcI4AAsY5l2EwBCAAAQhAAAIImAh9YGZmRvSLCwIQgAAEIOAagbeu\/HW5+drlrpld114ETB1EKlwGBgZkamqqLkwKQAACEIAABIpE4PylV8pP\/82oPHbHDd6JGARMnZ727LPPyubNm2VsbEza2tqK1C+X2KIia3x83Alb1XjsTbc7wRe+hgB9oZx94aUz5+Q\/fuVlWfbzV+WxO66Xrq6udEFkXDsCJqKAmZiYKLzzjdhywVbFjr3pvu3wha8hQF8oZ1+489B35chT\/0cue2ZMDv+XPyn8GBbXSwgYBEzcPpNYeT5UE0NZsSL4whcBk24fKDLfQ9\/8kYw+\/n35\/Wt+Ig99+g5x5Q\/bOB5DwHgkYDRf55FHHpFbb71V2tvb4\/SDXMpib7rY4QtfQ4C+UK6+cPL0Obn+s38rQ+s6ZO2v\/r8gDQIBk24fKGTtrv0VW0iIGAUBCEAAApkQUPGiU0ea\/\/KdP\/6Xzk3Vx4FEBMajCEwcx1MWAhCAAAT8I6DTRjp9pKuO3tt6CQLGPxdHbxERmOisKAkBCEAAAvkRUOGi0ZeHet8vvR98d2CIz2MYERgiMPm9bTwZAhCAAAQSIWDyXlS4qIAxFwImEbxuVuKz8930CFZDAAIQgIBNIJz3Yv\/O5zGMCAwRGD4JIAABRwlwzImjjkvAbF1palabhvNeEDAJAPahCp\/Vqw\/+oQ0QKCsBjjkpq+fn293Z2Rnsuj597jLZ+IVvL8p7QcCUu2+UYv4QF0MAAu4ScOmYE3cpF9NyczTE6EP\/VXYeu1huet\/yRXkvCJhi+i1zq4jAZI6cB0IAAhEI8NkUAZKnRYzvf3bzoKxcvSbY76Xa5XM\/IQemTgf32fmevts0CwKlIMBnUyncXLGRxvetG++TR+\/bFOz3goApb38opfNxNwQg4C4B1wXM6dOnZcuWLXL8+PFFTtiwYYOMjIxIS0tLps4xPM1Dd+\/eLT09PU3ZMDc3J8PDw9Lb25voQYrGVp1Cuu3f3lLTRtf7Sa3GEYEhAtPUC8rNEIBAPgRcH5iMgBkaGkp0cG\/EG8ryU5\/6lBw4cEBWrVolRnh0d3c3LWIasafePXF8H6dsvecW7fcIGARM0fok9kAAAhEIuD4w1RMwRkQcPXo0oGEOIzxx4oQMDg4uENq\/f3\/wfxPNWbNmjejPWltbRcv29fXJqVOnxP65jbeaWNF7d+3aJXv37g3qOnLkiGzbti24deXKlYHY0WXMGmE5e\/asPP3008HPHn300SXfmwhMNXvUl6Ojo\/Lyyy\/LBz7wAbnvvvvkM5\/5jJi2h6NBcXwfp2yEbleoIggYBEyhOiTGQAAC0Qi4PjDVEzA6oKvw0OkknWYyERKlo6LkgQceCCI34WkaFRrT09Py8Y9\/PBA1JsJj12dPTxlBtGfPniD6UukyAsMII61Lrz\/8wz8MBIwKGn2OsSX8vQqY1atXV7XHbp\/aYNqgdSqnHTt2yP333x8IKb3i+D5O2Wg9rzilEDAImOL0RiyBAAQiEwgPTLobq34V9dJEUzvZtFIOjImSqMCwc0fsKMmNN94YRGCM4AhHSkz7w6KjWrlqP6\/F0QgMI2DMVFM4mmOLK61PhY8dHTIRnhdeeGHR72wBU01Qbd68eSEqVctWBExR34gM7PLZ+Rng4xEQgEBKBMKfTbob6+jj0yk9rflqh9Z1yNC6axYqqhWBqfQ7Hfw7OjpEBYw9tRMWKraA0UHevszUjx1piRKBCU9naZ1bt25diMCYKaJwNCgsYKrZ8+qrr8qhQ4cWJS9re\/ft2xeYb6bPwm0L\/7xZsdO8l7OtgQhMHd4ImGw7JE+DAASiEfAlAlMpibeaENBIR1jA1IrAhEVBJbLVcmDsqZsnnnhCJicnFwRGOAITVcBUs0d9We13asfdd98t99xzz8IUV5xxKU7ZaD2vOKUQMAiY4vRGLIEABCITcH1gaiYHxo7AhOtRcaFiQ6eZPvGJTyzkwOjPDx8+vDCFY4OutgrJ5LKYOjUfRwWP5tZo\/o2ZQooiYMI5MLY9OoVkCxgTbdJl3OTAVH8lEDAImMgfmBSEAASKQ8B3AVNrFZItYNQj1Vb32D+vNH0UFjH2FI9OEWl0SC87X0fzdDRB+PHHH19YLRRFwKjgqWZPOAITzg9iCqnye+elgLHVa7WPG3vTomrL6\/Re1z8kivNxiyUQgECSBPhsSpKmW3XF8X2csm5REPFOwJjEp1q7KBoVbJbh2eHB8O6PPjvftc6KvRCAwAUCfDaVtzfE8X2csq4R9UbAmJDb8uXLAx+sX7++6g6K9ZaoVQorRsn2ds352AsBCLhLwOeByV2vZGN5HN\/HKZuN9ck9xSsBc+bMmWBDId0\/oNoW0HG3iDbO7+\/vl87OzoC87r6oX1wQgAAE8iLg88CUF1NXnlvP9zMzM6Jfes3OzsrAwECkPWNcab+x0xsBYxpUT6CY369bt04efvjhYIfHKDkwtmNVzOgXFwQgAIG8CNQbxPKyi+emT6Ce78fHx0W\/7MvHWYTSCpiTJ08uLKertsW0Ot90lLGxMWlrayMCk\/67yRMgAIEIBOoNYhGqyLVIvWXUaRoXd\/dd84fvsWPHFg58NPaZM5KSOL06apvr+d6OwExNTQViBgETlW6O5aJGYOwpplo7MdbrKDk2lUdDAAIlJuD6Z5OLAkYPbfy93\/u9hfxKHW927twpzz33nGzatCmzk6vj+D5OWddep9JFYNRB4WXWtdS4z853rbNiLwQgcIGA659N9QSMvdXFhg0bgl1w9bJPf9aogqYA6M8qnVqt+8Vcfvnli35nNpQz6QNa92c\/+9lFXSscTTF\/GF9\/\/fXy4osvyr333iu6YlXHjoMHDwb3\/uZv\/mYgYMJ7uNj7ydjHA0T5ebX+Hsf3ccq69n6VUsCoQ+1DtczJombTItuJPjvftc6KvRCAQDkEjL3VhREoukAjfPqz+YNU\/9XPb\/28Dp9afddddwXCwk4V0CmW8GZ4Wkelbfv153bu5Je+9KWFbf11+kgvPf1az2nauHHjokUkdnTfPu9I79HIze233y7Vfl7tZGy9N864FKesa+9XKQRM+FwNuwPo\/426D+8BE7ejuOZ87IUABNwlEB6Y3np5Wt78cXEPc3zHr3XIRVd1LACvFYHRttlb65s\/Oh988MHgFGqTAhCuw\/6sv+KKKxadWm3XWUnA1LLHrvepp56StWvXBpEfI0IeffTRQMCoULIvWxBVOrDRjDFRzmxq9A9rBIy773jTlvvs\/KbhUAEEIJAbgfBn05m\/2CFnvrQ9N3vqPXjFbdtlxUfvjyRgwpuLmml+nepRAWO27g9P15jKdQqo0qnVRiiEBUzU3El9rl4qYj784Q8vRHF0RasRMPbUl5a1jzAwCb\/6c3uaqtrPqzGNMy7FKVvPh0X7vXcRmKQB++z8pFlRHwQgkB2BMkdgbAETPqnZeCCc21gtAqOR91p7h2l9dgRGc2h27NgRnAz9+uuvB1NXJq\/yd3\/3d4ODHs0J2\/WmpEw7jM2VZgsq9ag441Kcstn13mSehICpw9Fn5yfThagFAhDIg4Drn021pmzq5cDYA7+dw2jfp1NIdp5LNQGj0RO9KuVAVhMWJhnXLE2uJmA0svL5z38+WHr9rW99K8iV0eeY1UuaA1Pt5+TA1H+rEDAImPq9hBIQgEDhCPgiYHQ1kH0ZUVBrFZItYMKnVpupmVoRGL1HIyXnz5+Xq6++Wr72ta8tssFeIRSOwOip0ipMDh8+vGgvMTOFZE8HffrTn5bvfOc7wZRXeLWUsbOa\/bU6XBzfxylbuE5exyAEDALGtT6LvRCAQMyVKADzi0AcURKnrGuUEDAIGNf6LPZCAAIImFL3gTiiJE5Z16AiYBAwrvVZ7IUABBAwpe4DcURJnLKuQUXAIGBc67PYCwEIIGBK3QfiiJI4ZV2DioBBwLjWZ7EXAhBAwJS6D8QRJXHKugYVAYOAca3PYi8EIICAKXUfiCNK4pR1DSoCBgHjWp\/FXghAwAMBU+8wxzSdXOsA33rPVUFgzluqtVdLvXqa+X0cURKnbDM25XEvAgYBk0e\/45kQgECTBFwfmFwUMGYDOnWdOX26STc2dHsc38cp25AxOd6EgEHA5Nj9eDQEINAoAdcHpnoCptZGdmfPnpWnn35adNO78AZxZiM8E2W5\/PLL5ejRowFm\/Z0eBaCb2OkGenqvHuarZyzZl31Okf1zrfPgwYNy2223yec+9znZu3evtLa2Lhw1cP311web2506dWrhrCOzS68e9Bj1qIB6fSKO7+OUrffcov0eAYOAKVqfxB4IQCACAdcHpmaOEtADEs3W\/\/ZRAvb0jiLs6+uTu+66KzglWsupsBgZGZFKp1Fr+WpnFxl36C67epn69FRq3ZnXCJOTJ08GAkZPnh4cHAwOnmxvb184tVp\/rgdB1jq2IILrJY7v45SN8uwilUHAIGCK1B+xBQIQiEggPDCd\/\/mM\/NPPZyLenX2xt1\/aLssubV94cC0BY59bpIct6vcqQB588MFAFHR3dwciIlyHHeHQs5CMiNBcFbvOSgKmXkTIPr\/I1GfESKXTrO3Ii0ZuNMqzbNmyhahNMx6II0rilG3GpjzuRcAgYPLodzwTAhBokkB4YJp7\/nMy9\/x4k7Wmd3vLdf3Sct0nIwkYjXRMTk4G0RIVMGY6SEWAChj7NGozHWRbrlNAN954Y6TDHO0pICOMKlEwB0VqFMdcOgWlERdzonX4kElzPpKWtwVNs5TjiJI4ZZu1K+v7ETAImKz7HM+DAAQSIFDmCIwtYO6++2655557JLwiqNZhjnYExoiPWuLFCBD9157+MaJk48aNMjw8vBAZCkdk1Fdf\/epXA6\/rCdTNrl6KI0rilE2gW2ZaBQIGAZNph+NhEIBAMgRcH5iayYEJRzqMsDBRkgceeEB0CmnXrl0LUzbVppAefvjhwCG18lKq2WoiRffdd5985jOfCeoxOTZm+krt2LFjh9x\/\/\/1BbowmAd97771B1KbRK47v45Rt1J687kPAIGDy6ns8FwIQaIKA6wOTEQW6Gsi+zCqiWquQbAFjoh1mpZFZQVQrAqP36NTT+fPn5eqrr5avfe1ri2zYunXrIkFjcnB0ukinnMxlBNPOnTvl0UcfFU0u3rdvX\/Bre4WUba+Knunp6aYSeeP4Pk7ZJrpjLrciYBAwuXQ8HgoBCDRHwOeBqTky2d+d1PLoqJbH8X2cslGfX5RyCBgETFH6InZAAAIxCPg8MMXAUIiiCJh83ICAQcDk0\/N4KgQg0BQBBExT+Jy+OY7v45R1DQoCBgHjWp\/FXghAwDoLqb+\/Xzo7O2FSIgKzs7MyMDAQ5NnoRnq1LgSMYx0jznp7TcKyNzsKN9Vn5zvmVsyFAAQsAroUWAexqakpuJSQgIrWsbGxYKdfBIwnHUDFi2aBVzvLwm6mmbc8duyYHDhwoOLafASMJx2DZkDAQwIqYvSLK1kCL505J3ce+nu56X3LZWhdR7KVJ1SbCpd64kUf5fMY5s0UklmSt3z58qB7rF+\/Pthqup4yNedo6O6OlTYX8tn5Cb1HVAMBCEDAGwKHvvkjufPQdwPhMrTuGufb5fMY5pWAOXPmTLAO394RsVrvU8Gjmwvp+nwVMQgY599TGgABCECgKQKjj39fRh+f9ka8EIFpqjtkf3OlQ7UqWWFOFdXzMqLkwNiJclFDd9m3nidCAAIQgEAjBFS8aPSl94Pvdj7yYk8txkn4bYRbnvd4E4ExEKMIGE3cNds5q6OjCBjbSSpm9IsLAhCAAATcJnDy9Dk59M0feiNe1Bvj4+PBl31FWbHkmidLKWB0ymjt2rXB8rOoq5A027utrS3wLxEY17o59kIAAhBYSkDFi0ZevvHiWXlw0\/vl5mvncyhdv+wIjK5SUzGDgHHAq\/UiMNXO39CmVXKwzwlQDrgTEyEAAQikQkDFy8YvfDuo+7E7bpD3tl6SynPyrtTnMayUERi7Q0WNwPioXvN+sXg+BCAAgTwIlEW8KFsETB49rMFnVorA1DqnAgHTIGhugwAEIOAgAbPSSPd4OXrnDQ62IJ7JCJh4vLwq7bPzvXIUjYEABCBQg4BJ1vVtmXQ9p\/s8hnk3hVTPmXF\/77Pz47KgPAQgAAEXCdjJuj4sk47jA5\/HMARMnZ7gs\/PjvASUhQAEIOAigTLlu1Tyj89jGAIGAePiZxI2QwACEKhL4JnvnQ1WGukKI59XGtUCgYCp2038LeCz8\/31Gi2DAATKTsDHYwEa8anPYxgRGCIwjbwT3AMBCECgkAR83Fm3GdAImGboOX6vz8533DWYDwEIQGARgbLnu1TqDj6PYURgiMDwEQgBCEDAeQJmyqjM+S4IGOe7cbIN8Fm9JkuK2iAAAQhkT6Cs+7tEJe3zGEYEhghM1PeAchCAAAQKRcCeMirb\/i5RHYGAiUrKw3I+O99Dd9EkCECgJASYMormaJ\/HMCIwRGCivQWUggAEIFAAAkwZxXMCAiYeL69K++x8rxxFYyAAAe8JMGUU38U+j2FEYIjAxH8juAMCEIBAxgSYMmoMOAKmMW5e3OWz871wEI2AAAS8JlDmgxiTcKzPYxgRGCIwSbwj1AEBCEAgcQJ21OXBTe+Xm69dnvgzfK8QAeO7h2u0z2fnl9itNB0CECgwARJ1k3OOz2MYERgiMMm9KdQEAQhAoEkCeoL0Jw5\/N6iFvV2ahCkiCJjmGTpbg8\/Od9YpGA4BCHhHwI66qHB5qPf93rUxjwb5PIYRgSECk8c7xTMhAAEILBBgeXR6nQEBkx7bwtfss\/MLDx8DIQABrwnYUZeb3rc8iLroYYxcyRHweQwjAkMEJrk3hZogAAEIRCRw6Js\/El1lpBe5LhGhNVAMAdMANF9u8dn5vviIdkAAAu4QMPu6qIAxwoWoS3r+83kMIwJDBCa9N4eaIQABCPySgD1dpIKFfV2y6RoImGw4J\/aU0dFR6ejokJ6enop1nj59WrZs2SLHjx8Pfr9hwwYZGRmRlpaWJeV9dn5iwKkIAhCAQA0CJOnm1z18HsO8i8CoeNm3b5\/s3r27ooCZm5uT4eFh6e7uDn5vvl+5cqUMDQ0hYPJ7z3gyBCDgGYFw1OWxO24gSTdjHyNgMgbeyONMVGX58vmtptevX181AhOu\/8iRIzI5OVkxCuOz8xvhzD0QgAAEohAgSTcKpfTL+DyGeROBUQFz5swZ0UiKHWGJ0j0QMFEoUQYCEIBAfQIadbnz0HflGy+eJUm3Pq7USyBgUkec3APCU0T1ajaRm02bNlWM2Bjn9\/f3S2dnZ1Bde3t78MUFAQhAAALzBEjSLU5PmJmZEf3Sa3Z2VgYGBmRiYkK6urqKY2QClngTgTEs4ggYU1bvrZfEa7NWMaNfXBCAAATKTiAsXNjTJf8eMT4+LvplXwiY\/P1S14KoAiaKeNGHmQjM2NiYtLW1EYGp6wEKQAACZSFAnksxPW1HYKampgIxg4Appq8WWRVFwNRbeWRX6PP8oQPuxEQIQKCABOzN6HRPF1YXFdBJvzTJ5zGslFNIutT61KlTVaeNEDDFfRmxDAIQyI9AeLpoaN01QaIuV3EJIGCK65slllWKwJif9fb2yurVqxdtYmcqWLNmjezfv19aW1sX1emz8x1yK6ZCAAI5EiDPJUf4TT7a5zHMuwhMk75ecrvPzk+aFfVBAAJ+EUC4uO9Pn8cwBEyd\/umz891\/NWkBBCCQBgEjXDRJVy+dJur94NXsopsG7JTr9HkMQ8AgYFJ+fageAhBwhYAKF92AbvTx7wcm3\/S+5aJ5LpwW7YoHl9qJgHHXd01b7rPzm4ZDBRCAgBcEbOGi\/zd7uSBc3Hevz2MYERgiMO6\/obQAAhBoiEBYuJiIy83Xzp8px+U+AQSM+z5suAU+O79hKNwIAQg4TcDOcTERF81xQbg47daKxvs8hhGBIQLj3xtLiyAAgYoEwsm5GnFBuPjdWRAwfvu3Zut8dn6J3UrTIVAqAuHl0CTnlsf9Po9hRGCIwJTnTaalECgZAYRLyRxeobkImBL3AZ+dX2K30nQIeE0gfFYR+7h47e7SziIQgSECU943m5ZDwCMCZkXRxP\/+YbCXiy6BNsuhPWomTYlJwOc\/whEwCJiYrwPFIQCBIhEIryhS4cIhi0XyUL62IGDy5Z\/r0312fq5geTgEINAUgUr5Lawoagqplzf7PIYRgSEC4+VLS6Mg4COBatNEnFPko7eTaRMCJhmOTtbis\/OddAhGQ6CEBOxoizafaaISdoIGm+zzGEYEhghMg68Ft0EAAmkSqBZtuel9K9gxN03wntWNgPHMoXGa47Pz43CgLAQgkA0Boi3ZcC7LU3wew4jAEIEpy3tMOyFQWAJGtDzzvbOLlkATbSmsy5wxDAHjjKuSN9Rn5ydPixohAIGoBCpNEekW\/zdfuyLYv4ULAkkQ8HkMIwJDBCaJd4Q6IACBCASMaHnme2fk0Dd\/FNxBQm4EcBRpmAACpmF07t\/os\/Pd9w4tgEDxCahoCZ8CbXbJZflz8f3nuoU+j2FEYIjAuP5+Yj8ECkegWqRFp4jYbK5w7vLaIASM1+6t3TifnV9it9J0CCROANGSOFIqTICAz2NYoSIwR44ckW3bti1y2e7du6WnpycBNzZWhc\/Ob4wId0EAAoYAooW+UHQCPo9hhREwKl4OHz4s+\/fvl9bW1qBPnD59WrZs2SKbNm1KTcSMjo5KR0dH1fp9dn7RXzzsg0ARCYSXPKuNJqeFZc9F9Fi5bfJ5DCuEgDFCZWhoSLq6uhb1NoWvIsMWNkl1R6133759UivK47Pzk+JIPRDwmYCdhGsiLipY3rPiEtn8W1cH\/9587XKfEdA2hwn4PIaVUsAYwbR8+fyHzvr164nAOPyCYjoEkiYQ3qPFRFnYpyVp0tSXNgEETNqERSTLKSQVMGfOnJGVK1fK8PCwdHd31xUw\/f390tnZGZBob28PvrggAAE\/CFTKZWFqyA\/flrEVMzMzol96zc7OysDAgExMTCyZ4XCdTSEiMAZi1km8c3NzkQTM6h\/3NuXnZZfWFjvLWpb+\/u2he5Zd2rbEBrvecB3m\/nrPbqph3AwBRwlUmhYiyuKoMzF7CYHx8XHRL\/tCwHjWUaIKmD\/\/07vkw7+\/NdiM6pkXzy6hoPPh711xidx07fIgmc9c538+r4CrXed\/PlvxV\/9U4b7zcxfqqldvteeFxYwRPRfEzgWRZMouLUPkybPXoBTNsSMs5v9hwUIuSym6QikaaUdgpqamAjGDgPHM9VEFzObNmxc536xCUBx6+NpLZ+Z32jSXihidK5\/\/d0XqCX5hQWMLICN8wmWMeIpStpLbK0V\/VAjZkSItYwsgokGevUAFbg6CpcDOwbRMCZADkwJue+XR6tWrg+XSx48fr\/ikNWvWpLIKqVEBU8lIE5L+xotnAjFj\/5Vn\/tKz\/8LLQtg06zZb9Biho4LI\/nlYCIV\/X0v8qLgJR38QPc16rXz3m\/dN\/5CwzxgiwlK+vkCLlxJAwHjaK5IUMNUQGSFz8vRcEK35RmgKykw5mdNnXRA2UbuDETq1xE+139nPqDSdZSI9F10xv+x+PvrD9FZU37hazkQ69T1SsRKeDtJ2sVLIVe9idxoEEDBpULXqzGsfmCwETLVojf5cP4SNsKk0DWU+jIMcm9aWUuw3UUv0aLRHBU+tKI8tduwpLRPZQeik\/DInWH04slIpqmnECvkrCYKnKq8IIGBSdmdeAiZKs7J2vh2xqTQVZcLi4eko\/XkZN9OyBY8tbN58ZSpw71uvPlvRzbbQeceV88vjETlR3oh0ytjiJBxZMX0+y7yydFpJrRDInkDWY1iWLcx1GXWlZdOVGr9161bRXXrzuIrifDvHRjlUSh42H\/RljNzU6xu20DGipl5Ep5bIeceVi3eMrvd8fj9PIEpUxfTfm69dUYqoI30DAmkSKMoYlkYbcxUwpkG1IjBpNDpOnS44PyxuqkVujMCxozdlmZqK4vNK0ZwoIkeTkcNRnDILHDtPZV5sL85VCQtt\/R6xEqWHUgYC8Qm4MIbFb9X8HYUQMI0an8V9rjvfHkw032Z+QFm69NuwNGe8zOfdzO9po4nF84MM572oyAnn4eh0VaW8HBPBeccVXQvLyzXp2Ic8nHAkxURXqiWpm1wVLUe+ShafXDwDAvMEXB\/DavkRAVOnl\/vsfDPo6L8modgIHPOzSngqiRwTySm70AkLHBPBqZSLY3JuTPSmKOLGiF791yxNNn0lnGyuP7f7w7z\/56d+bBHMYAIBCORDwOcxrDAC5sSJE9LX1yenTp1a4uW09oGJ0p18dn4RE06NAAAet0lEQVSU9puBywxmJopj\/wVub+Jn12kiOOHBzBY7ZRvkjMAxgqZW9MaemlJx08y0lC1K1EfzQmRuYQPGC\/69sCGj8eWFSNx8BE4FCpGUqG8P5SCQLwGfx7BCCBh7OfPGjRuD84l6e3vFbHCnCbxdXfkkTfrs\/CRfq\/Bf7WGhc2HQXDpAhgdKFTzmL\/vF\/7YsDJz2PfbxDUm2Keu6KombcOQmPC2l38+eu1Kmzv76EkFSj3lYYBrWRmCWTVxm7W+eB4EsCPg8hhVCwISTeEdHR6WjoyM4IVrhHzp0SEZGRqSlZX4Ay\/Ly2flZcgw\/q5rg0XImumMGYPOzKPbaYsYIobBAskVRuE4dvOtd4Xrt8hrZqHUZYWfKhKNX4e+1PrPZX+fy56X9klek7ZJXgn\/1e\/uaOXelzJ67IhAz+n\/Ntbnois5fTuVcEH8Ik3oe5vcQ8IeAz2NYIQWMLq+enp4Olk4rfBU0+\/fvl9bW1sx7lc\/OzxxmAg+0B3jzfyMaaokD+76wyKg2BZaAuRWrqBQxskVR+PdmI0NTmR2hUiGjAseekqoUtTHTURfybvKJaKbFlHohAIHKBHwewwohYBS7ihS9wqLliSeekMnJSSIwvJ25EGhE3BRhSivqdJRZIdVsjk0uzuGhEIBAXQIImLqImi8Q3tZfBc2+fftk5cqVcuDAAVm1alXzD2mgBp+d3wAObnGcgC1sKq2QMjk2F7\/n1qClCBvHHY75pSfg8xhWiAiMES+auJtXsm61Xu6z80v\/ZgNggYAKmzdffTY4aVxXRlWahiJaQ4eBgHsEfB7DCiFg2InXvZcCi\/0nUG8aSqM1iBr\/+wEtdJsAAiYD\/2nibp65LkRgMnAyj\/CCgAqbN176ctCWcLQGUeOFi2mERwQQMCk700Rgjh8\/XvFJbGSXsgOoHgJNEkDUNAmQ2yGQEgEETEpgXajWZ+e7wB8b3SVQT9SQKOyub7HcHQI+j2GFyIEpclfw2flF5o5t\/hEwp33r9JOugDLTUNpS+1woVj7553talB8Bn8cwBEydfuWz8\/N7pXgyBOYJ2InClfJpTJTm4vd8JBA5XBCAQDwCPo9hCBgETLy3gdIQSJmAPfU09\/z4wtOI0qQMnuq9JICA8dKt0Rrls\/OjEaAUBPIlUC1KY2+6x7RTvj7i6cUl4PMYRgSGCExx3zwsg0AFAlGmnRA0dB0IzBNAwJS4J\/js\/BK7laZ7RsBMO1XLo0HQeOZwmhOZgM9jGBEYIjCRXwQKQsAVAggaVzyFnWkTQMCkTTih+nU3323btgW17d69W3p6eqrWbJfdsGFD1dOufXZ+QtipBgKFJ4CgKbyLMDAlAj6PYd5EYE6cOCGDg4OyZ8+eoBuY\/1c6xVodqqdd79+\/X1paWmR4eDg49XpoaGhJF\/LZ+Sm9L1QLgcITqCRo7FVOLNsuvAsxMCIBn8cwbwRM+CwlFSgdHR0VozDhsrXOYfLZ+RH7P8Ug4D2BN37wjJw\/93eLznZSQWP2oWm57pPeM6CBfhLweQzzRsCoYNHLRFHC39tds1IEpru7u6LY8dn5fr6utAoCzRHQ6Mybrz4bbLIX3odGBQ0Jwc3x5e5sCfg8hnklYOyIi0ZVpqenK04LaffRKae+vj45deqUTExMSFdXV8VeZZzf398vnZ2dQZn29vbgiwsCEPCfQLXpJqIz\/vve1RbOzMyIfuk1OzsrAwMDNcc5V9tZSgGj4ubw4cNBDkxra2uQD2NHb8LRms2bNy\/yr4oZ\/eKCAATKReCtl6flF2\/7kbz16rPy+rG9suyyiwIAZrqJ6Ey5+kNRWzs+Pi76ZV+1\/lAvajvq2eWVgLFFSDVRMjc3FyTt2lNGdgJwOOnXRGDGxsakra2NCEy9HsXvIVAiArWiM4iZEnWEgjXVjsBMTU0FYgYBUzAn2eaEp4yqJfE2KmB8dH6B3YlpEHCOgB2dsXNn5oVMp5AI7JxLvTCYHBgH3BhnGXWlKSTNhRkZGQmWVduXz853wK2YCAFnCZjozBsvPRIkBNtTTSzTdtatzhnu8xjmzRSS9qpqG9mZqEtvb+9Csq5GaPbt2xd0Rjayc+6dxGAIOEWAqSan3OWVsQgYr9wZrzE+Oz8eCUpDAAJJEDBi5q2XfyBvnv7LRZEZ8maSIEwdZZlF8CoCk0a3RcCkQZU6IQABJVApb4YVTfSNJAn4PIYhYOr0FJ+dn+RLQl0QgEDzBEx0xiQBI2aaZ1r2GnwewxAwCJiyv9+0HwKFJICYKaRbnDMKAeOcy5Iz2GfnJ0eJmiAAgTQJIGbSpOt33T6PYURgiMD4\/fbSOgh4RsAcPBmeZmJptmeOTqg5CJiEQLpYjc\/Od9Ef2AwBCFwgoGLmH196ZMlqJjbNo5cYAj6PYURgiMDwpkMAAo4T0NVMuiT7jR98Q86f+7ugNSYBGDHjuHObNB8B0yRAl2\/32fku+wXbIQCBygQq7QDccl3\/LwXNR8BWMgI+j2FEYIjAlOx1prkQKA8BI2bMydkmKkO+THn6AAKmPL5e0lKfnV9it9J0CJSOwOvf\/HM5\/8Yx8mVK5nmfxzAiMERgSvY601wIlJvAwpLs7x4RWfajhVwZjjHws18gYPz0a6RW+ez8SAAoBAEIeEug2v4yJP7643KfxzAiMERg\/HlTaQkEINAQAV3FpFNMb57+iiy7\/JWgDk38JSrTEM5C3YSAKZQ7sjXGZ+dnS5KnQQACLhColvhLVMYF7y210ecxjAgMERg330qshgAEUiVQLSrDCqZUsSdeOQImcaTuVOiz893xApZCAAJ5EtCozOvf2ssKpjyd0OCzfR7DiMAQgWnwteA2CECgbAQ0KvPGyS\/La38zKr\/yvsvY7deBDoCAccBJaZnos\/PTYka9EICA\/wR+8t8\/K2+89IhcdPXZoLE6tWR2\/PW\/9e600OcxjAgMERh33kQshQAECkfATC9pZGbZZRcRlSmYhxAwBXNIlub47PwsOfIsCEDAbwI6vTT33cOixxYwvVQcX\/s8hhGBIQJTnDcNSyAAAecJmNVL9vSSTi2xeikf1yJg8uFeiKf67PxCAMYICEDAWwJv\/OAZ+cn\/vGshT0Y3xrvshrFgmokrGwI+j2FEYIjAZPMW8RQIQKC0BDQq89On\/5OcP\/d3C3kyJipTWigZNRwBkxHoZh9z5MgR2bZtW1DN7t27paenp2qVxqlaYM2aNbJ\/\/35pbW1dUt5n5zfLm\/shAAEIxCEQ5Mk8Px4sxSbhNw65xsv6PIZ5E4E5ceKEDA4Oyp49ewJPm\/+vWrVqiee1bF9fnzzwwAPS1dUlKnwmJydlZGREWlpaFpX32fmNvxLcCQEIQKBxAiZPRsUMCb+Nc4xyp89jmDcCJixCRkdHpaOjo2IURstOT0\/L0NBQXf\/77Py6jacABCAAgRQJGCGjK5cuWbOcJdgpsPZ5DPNGwKhg0cuIkvD3pl\/Mzc3J8PCwdHd315xiMuWN8\/v7+6WzszP4cXt7e\/DFBQEIQAACzRNAyDTP0K5hZmZG9Euv2dlZGRgYkImJiWDGwafLKwFjR1yqRVmMgFm3bp08\/PDDcvz48Ug5MLbTVczoFxcEIAABCCRHACGTDMvx8XHRL\/tCwCTDNpVawlNG9QTMyZMnFxJ39d5Tp07VzIEZGxuTtrY2IjCpeI9KIQABCFwggJBprjfYEZipqalAzCBgmmOa6t3NTCHZCcDhpF+f5w9TdQiVQwACEGiSAEKmSYAi4vMY5s0UUjjiUiuJN\/w7FTC7du2SvXv3LllK7bPzm381qAECEIBA+gQQMo0z9nkM80bAxFlGrQ5VEWP2fqmW8KtdxmfnN\/5KcCcEIACB7AmokHnlix8PNsRj1VI0\/j6PYd4IGHVltY3sTOJub2\/vQha2vZHdhg0bKua\/IGCivSCUggAEIJAlgbnnnpRXD35clr3rFYRMHfAImCx7ZsGe5bPzC4YacyAAAQjEInDmL3bIT\/7HZ+WSf7F8YUM8jihYjNDnMcyrCEysnh+xsM\/Oj4iAYhCAAAQKTcAImUu7r5SL\/tklwYZ477x+TN5xpV\/7njTiBJ\/HMARMicNvjbwM3AMBCECgiAQ0P+a1rx+U1\/5mRFTI6FlLF7\/nI6IRmTKffo2AKWJvzcgmn52fEUIeAwEIQCAzAipkfvxQn8jbvxPkx+ilIqbluk9mZkORHuTzGEYEhghMkd41bIEABCCQCAFN9P3h9t8O8mPKvGIJAZNId3KzEp+d76ZHsBoCEIBAdAImP+ZX\/vllC0KmTPkxPo9hRGCIwET\/JKAkBCAAAQcJmGmlf\/zBM0F+jEn0vbz7kPf5MQgYBztsUib77PykGFEPBCAAARcIvPbkF+Xlh\/oCAWMSfX3Pj\/F5DCMCQwTGhc8dbIQABCCQCAGzWunMl7aXIj8GAZNIt3GzEp+d76ZHsBoCEIBA8wTsaSV7IzzfppV8HsOIwBCBaf6TgBogAAEIOEpAk3w1GuPrtBICxtGOmYTZPjs\/CT7UAQEIQMB1AiYac+65J72bVvJ5DCMCQwTG9c8e7IcABCCQCAETjdFdfC\/719fLsstfCVYpuTythIBJpGu4WYnPznfTI1gNAQhAID0CdjTmst9ZL7+y6mdy\/uczzu7m6\/MYRgSGCEx6nwTUDAEIQMBRAgu5MVd1yIrN\/17ePP2XTh4SiYBxtAMmYbbPzk+CD3VAAAIQ8JWAHY25ou8\/y9t+9Tl569VnnYrG+DyGEYEhAuPrZw\/tggAEIJAIARONedctH5PLPrRG5p4fdyYag4BJpAu4WYnPznfTI1gNAQhAIHsC5nDIi67qEI3G\/NM\/\/lUQjbn4PR+Rd94wlr1BEZ\/o8xhGBIYITMTXgGIQgAAEyk3AnlJacdt2uex31slrk72FjsYgYErcZ312fondStMhAAEINEzATCld8hu3yLsH\/5u88dKXg2mli67okstvOtRwvWnc6PMYRgSGCEwa7wx1QgACEPCagD2ldNWdB4KdfIsYjUHAeN0NazfOZ+eX2K00HQIQgEDTBOwpJRUx7\/ytD8nPvj1QqNwYn8cwIjBEYJp+iakAAhCAQJkJ2KuUVMi8+cqzhYnGIGAc6ZlHjhyRbdu2Bdbu3r1benp66lp+4sQJGRwclD179siqVauWlPfZ+XXhUAACEIAABCIRMFNKmhezcvvXg917TTSm5bp+abnuk5HqSbqQz2OYNxEYW4hoB6glSkwHmZubk+HhYTl27JgcOHAAAZP0m0N9EIAABEpEwEwpvfXjadFITMtv3CJzz39uYd+YPM5UQsA40AE1+jI5OSkjIyPS0tIio6Oj0tHRUTMKo47VcnoRgXHAyZgIAQhAoOAE7LyYq7d\/PRAxdjTmXd2H5B1XdmXWCgRMZqgbf5ARIkNDQ0El4e\/DNZ8+fVp27Nghvb29QVkETOPsuRMCEIAABBYTePmhPnntyS+K7hez4qP3B7800ZgsN79DwDjQM8MRF43ITE9PixE04Sbo7\/W68cYbI+XA9Pf3S2dnZ3BPe3t78MUFAQhAAAIQqEZABYwKGT2CQKeU9MoiwXdmZkb0S6\/Z2VkZGBiQiYkJ6erKLvKTRa\/wJgcmjoDRfJmDBw\/KvffeGzg5ShKv7QwVM\/rFBQEIQAACEKhFIJzcq2XTTvAdHx8X\/bIvBEyB+2mcKSQtu3bt2kCNRl2FNDY2Jm1tbURgCtwHMA0CEIBAEQloXszJO64Rs0LJ2GimlJLewdeOwExNTQViBgFTxJ7xS5vCU0bVkng192XLli1y\/PjxJa2p5GCf5w8L7E5MgwAEIOAVAXuF0nu\/8P2FtqU9peTzGObNFFIjy6i1B0WNwPioXr36dKAxEIAABApOwBYxK3d8XfRk67SnlBAwBe8UxrxqG9mZ\/V50xVE4iQkB44hzMRMCEICAJwRObf9tsfeKSXNKCQHjSadppBk+O78RHtwDAQhAAALNE6gmYpKeUvJ5DPNmCqn57lS5Bp+dnxYz6oUABCAAgfoEVMSce+5JMRvemTuS3PjO5zEMAVOnj\/ns\/PqvFyUgAAEIQCBNAtUiMfpMs0qpmbOUfB7DEDAImDTfTeqGAAQgAIE6BGqJGDOl1OhSawRMibufz84vsVtpOgQgAIFCEaglYnRK6aeTvYG977x+LNZZSj6PYURgiMAU6iXGGAhAAAJlJWBEjL1PjGHRaF4MAqasvUlEfHZ+id1K0yEAAQgUkkAtEaMGx82L8XkMIwJDBKaQLzFGQQACECgrgXoi5o2Xviyvf3tAouTFIGDK2ouIwJTY8zQdAhCAQH4EVMTotXL71ysaoVNKZ\/\/6X8myS9tl+YeermooAiY\/H+b+ZJ+dnztcDIAABCAAgYoEqh0AaRc2eTHn52aqJvf6PIYxhcQUEh8fEIAABCBQQAJzzz0pP9z+27Litu2y4qP3V7Xwp9\/olbdefVbe1X1oyQolBEwBHZuVST47PyuGPAcCEIAABBojYERMeLfecG2aE6O5Me+8YUwufs9HFn7t8xhGBIYITGNvFXdBAAIQgEAmBM78xQ4586XtS44cCD+80golBEwmLirmQ3x2fjGJYxUEIAABCIQJ1FuZZMqbnXs1CqPRGJ\/HMCIwRGD4pIAABCAAAQcI1FuZFBYxusz6\/769XzZv3iwTExPS1dXlQCujm4iAQcBE7y2UhAAEIACB3AhETepVA+3jB278g1kETG5ey\/HBPoffcsTKoyEAAQhAoAECUfNhjIj52bcH5OSJZ+WyG8akY82F5N4GHl24W4jAEIEpXKfEIAhAAAIQqE4gaj6M1qB\/hD\/\/1V657t8dYgqpbJ2KCEzZPE57IQABCBSfwMk7rpGLfq2j6k69pgU+j2FEYIjAFP9NxUIIQAACEFhEIGo+DAKmxB3HZ+eX2K00HQIQgIDzBKLkw\/g8hhGBIQLj\/EtMAyAAAQiUlUC9pdUImLL2DE6jLrHnaToEIACB4hMwU0lX3XlA3nXLx5YYjIApvg9Ts9Bn56cGjYohAAEIQCAzAmYq6b1f+L5cdFXHouf6PIZ5NYV05MgR2bZtW+C83bt3S09PT8UOdPr0admyZYscP348+P2GDRtkZGREWlpaSqVeM3u7eBAEIAABCKRKoNqqJARMqtiTqfzEiRMyODgoe\/bsCSo0\/1+1atWiB8zNzcnw8LB0d3cHAsd8v3LlShkaGkLAJOMOaoEABCAAgQwJvPbkF+Xlh\/qWHPiIgMnQCY0+SqMvk5OTC5GU0dFR6ejoqBqFsZ8Tvtf+nc\/Ob5Q190EAAhCAQPEIVNrgzucxzJspJBUsepkoSvj7Wl0NAVO8FxGLIAABCEAgHoFKCb0ImHgMcykdjrioKJmenq44LWQbaPJhNm3aVDFaY5zf398vnZ2dwa3t7e3BFxcEIAABCECgSAR0Gkmnk358998GZs3OzsrAwACHORbJSWFbGhEwJv9F66qXxGs\/T8WMfnFBAAIQgAAEikbgH257m\/zZqeXyZz9csWDaxMQEZyEVzVHGnrhTSFHEi9ZtIjBjY2PS1tZGBKaoHQC7IAABCEAgIGBHYaampmR8fJwITJH7RnjKqFYSb72VR3Y7fZ4\/LLI\/sQ0CEIAABBonoFGYFbdtl+ffu042b96MgGkcZfp3Rl1GrZaouDl16lTVaSMETPr+4gkQgAAEIJAeATsKg4BJj3NiNVfbyM5EXHp7e2X16tWLNrEzD1+zZo3s379fWltbF9lDBCYx91ARBCAAAQhkSECjMP\/wgf8gf\/CnjxOByZB7YR7lkoCZmZmRRx55RG699VYnVklhb7rdHL7wNQToC+XsCyYK86Fj1yBg0u0CxazdJQHjkq3qbexNt8\/DF76GAH2hnH3B7Atz9wtXy9C+R1mFlG43KF7tLr34LtmKgEm\/r9Mf0mXsEl+XbOWzIdl+q2ck\/dXfvyI37HoCAZMs2uLXVmkju6JabTYssjfdK6qtahf2pusd+MLXEKAvlLcv\/Or3\/pdc\/JU\/lmX3PC0dN9ycLoiMa\/fmKIG0uOncse5iqGvpuSAAAQhAAAKuEfjrD3xf3nXLx+SqOw+4ZnpNexEwEdypIka\/uCAAAQhAAAKuEXj3xW95F31RHyBgXOuJ2AsBCEAAAhCAAAKGPgABCEAAAhCAgHsEiMC45zMshgAEIAABCJSeAAKm9F0AABCAAAQgAAH3CCBg3PMZFkMAAhCAAARKTwABU\/ouAAAIQAACEICAewQQMDV8pqdW79u3LygxMTFRmF0M9eTtvr6+4ETtDRs2RD5V+9ChQ5HKJt2N49hrM1+5cqUcOHBAVq1albRJNeuLY699gGi1A0HTND6OrcYOc7hpd3e39PT0pGnekrrj2Guz1YqKzvf06dOLDorN4zMjKt8wW+Oo3bt3Z9onotqr9tllXfhsMJug5tV3M32xc3oYAqYKeO18OpjqCdUvvPDCwv\/Dp1Vn7Td78Nm4caMMDw9LvYHIvEhRxU6SbYpjr36oTk5OLogs\/f7w4cMVTwlP0ka7rjj22n1E+4X2FxWVIyMj0tLSkpaJC\/XGsdU2xgxeWQ9Wce1Vnh0dHZkOqI32BdM2HViHhoaCwXZwcFD27NmTmQCPy9dua7gvp955RSSOvUYcKtuuri4p+meDEVsPPPBAYK\/yzesPyCx8mdczEDBVyOuHp176wpgXrbe3N\/coTPiDsd6Loe04evSo3HLLLfLaa69lNrgarHHttd2RxyDQjL1ZDwKN2KoDwd133y1nz56VTZs2ZSoO4thbhHcujr1adteuXbJ3717J64+cOPba71lYHGQ1GMWxN1y26J8N4T\/GtD\/v3LlTbr\/99swEbVZ+zPM5CJgK9MMh9jxD7mHzwoNkvUFTf2\/+YrGjG1l1urj25i1gmrHXFr1Z8G3EVrXxgx\/8oHzlK1+pG7lLug1x7M1rUK0Vlaj1roUHrKTZRakvDt9wRK7onw2VIjBZ2xyHbyUBo9HyIvwRHKUvuVIGAVNDwNidLe9wtjEzHHGJ+pdfXh+wjdqr7c16Skaf2Yi9Zoou63n5uLZqXzl48KD80R\/9kezYsSMXAWOH0Wv1XTvfwfT9rHNK4vDV92t6ejowNa+8uTj2GqZ5CsW49po\/JDWivHXr1iA6nuUVx97w6d\/m+6ynbbPkk8ezEDAImFT7XZyXPvxX4ec\/\/\/nMk3gbtdeIn0996lOZ2RzHVjuE3d7eHil3KumOEcdeLWuzDH+ftG2V6otjr8krMiKr6PbafxCZXL+sp77i8K2UU5K13XHsVb52ovTmzZuDKfx6+YpZ9GufnoGAqSFgTGdzeQrJNC\/PCIz9QVNvysu8+HmIFyNC4tprGGfdT+KEtLXsU089tSinK+sP0zj2hl\/LrNnG7QvVpgyyZNwIXxM5yjqaURa+RYh0+SRYwm1BwFTxrj1lVISEQmNmOOwe\/qugWmfNS8DEtTeP1QU2s7j22vdmHY6PY6u9PN22OctQfBx7qwmYLHMI4tgbfg\/z+MyIY6\/yzUMUNvquFUEgxuUbbmvWq9J8Fi6mbQiYKl62\/5pxfRm1iWpknfQW\/pCst+w7j7B7rb\/0o9hrR2uyFl9xlqHa7cxr4Ipjr2tL6sPiNUqkMekBJg5ffbZZkXbPPffksjImjr2VppCynK6N+1lmix3dUkETeM0S+6T9Xub6EDA1vO\/iRnbVQsJ5RWAUb63Nqmx7q0UJsk7ejGqvEYbbtm0LelHWSbxx2BZBwMS1184hyINtXHvtjexcsDePpcjhj9s475pJhM3rXYvbH+z+m8ceXGUQNgiYMniZNkIAAhCAAAQ8I4CA8cyhNAcCEIAABCBQBgIImDJ4mTZCAAIQgAAEPCOAgPHMoTQHAhCAAAQgUAYCCJgyeJk2QgACEIAABDwjgIDxzKE0BwIQgAAEIFAGAgiYMniZNkIAAhCAAAQ8I4CA8cyhNAcCLhLQPVT0gMn7779fsj6Tx0Ve2AwBCIggYOgFEIBA7gTss5pyNwYDIAABJwggYJxwE0ZCwG8Cugvz2rVrpaury++G0joIQCAxAgiYxFBSEQQgECZQ7dwl+6wgPStm586dcvvttwdn8ph7jh49GlSX17b8eBMCECg2AQRMsf2DdRBwnkClQy414qLX0NBQcFbWwYMH5d577w1+pgff6TUyMiIqbrI+JNN54DQAAiUhgIApiaNpJgTyIhA+qTn8vQoUvXp6egIxMzg4KHv27MnlhOS8GPFcCEAgPgEETHxm3AEBCMQkYEdc7OkjXXFk57+EfxfzMRSHAARKRAABUyJn01QI5EXAFiYPP\/ywdHR0BBGX8PJpBExeHuK5EHCPAALGPZ9hMQScI2CmjTZu3CiPPfbYwhRRePk0U0jOuRaDIZAbAQRMbuh5MATKRUBzXbZt2yYbNmxYlKCrFDQao5dZgaQrjzTBVy9ETbn6Ca2FQFQCCJiopCgHAQg0RUCFSF9fn9x1112BYFGxYi+fNpWzjLopzNwMgdIQQMCUxtU0FAIQgAAEIOAPAQSMP76kJRCAAAQgAIHSEEDAlMbVNBQCEIAABCDgD4H\/D+ToBTjm9saYAAAAAElFTkSuQmCC","height":337,"width":560}}
%---
