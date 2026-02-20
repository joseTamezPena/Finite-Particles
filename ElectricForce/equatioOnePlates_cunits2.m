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
figure(1)
plotrange = [0,0.99];
fplot(totfInf,plotrange) %[output:95a4f799]
hold on %[output:95a4f799]
relEq = sqrt(1-v^2) %[output:50ce538e]
relEqApx = taylor(relEq,v,Order=taylorOrder) %[output:18f3f6fe]
fplot(relEq,plotrange) %[output:95a4f799]
fplot(relEqApx,plotrange) %[output:95a4f799]
xlabel('v/c')  %[output:95a4f799]
ylabel('ratio')  %[output:95a4f799]
legend({'Force Carriers',"Lorentz Mass","Lorentz Apx"},'Location','northeast') %[output:95a4f799]
hold off %[output:95a4f799]



hold off %[output:95a4f799]
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
%[output:95a4f799]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQ9sX9d1329VuzYdV5Vpp4tEJqAbS25WBLKdBeSYAHKGtjI6SAXqJBQ1rIYgZB7iuAxiUyTtxrKySCJFCwgRp4BqqKq2jpLq2GuiboNtbLFrVzWXKok2eK4tp+VkksliW1TrOpRnRxnOUy51+fje7\/fe7\/d7f+55nwcQEsn77jv3c+773S\/PPffen\/vpT3\/6U8MFAQhAAAIQgAAEPCLwcwgYj7yFqRCAAAQgAAEIBAQQMHQECEAAAhCAAAS8I4CA8c5lGAwBCEAAAhCAAAKGPgABCEAAAhCAgHcEEDDeuQyDIQABCEAAAhBAwNAHIAABCEAAAhDwjgACxjuXYTAEIAABCEAAAggY+gAEIAABCEAAAt4RQMB45zIMhgAEIAABCEAAAUMfgAAEIAABCEDAOwIIGO9chsEQgAAEIAABCCBg6AMQgAAEIAABCHhHAAHjncswGAIQgAAEIAABBAx9AAIQgAAEIAAB7wggYLxzGQZDAAIQgAAEIICAoQ9AAAIQgAAEIOAdAQSMdy7DYAhAAAIQgAAEEDD0AQhAAAIQgAAEvCOAgPHOZRgMAQhAAAIQgAAChj4AAQhAAAIQgIB3BBAw3rkMgyEAAQhAAAIQQMDQByAAAQhAAAIQ8I4AAsY7l2EwBCAAAQhAAAIImAR9YGZmxsgXFwQgAAEIQMA3Ap2dnUa+tF0ImDoeFeEyODhopqamtPme9kAAAhCAgHICF666zvzDb46Z\/\/mFtepEDAKmTud9\/vnnzdatW834+Ljp6OgodVcXkTUxMeGFrQISe7PtTvCFryVAX6hmX3h1\/rzZ++ybZu7lU+bP7v1109PTky2InGtHwCQUMJOTk6V3vhVbPtgq2LE327cdvvC1BOgL1ewLY0\/8nRn\/T39tVj45ZHwZF9J4CgGDgEnTX1palg\/VluJcVhl84YuAybYPlJnvc6+cM587+qL519f\/vfnaFz+LgMmnK5TrKT4NApKv89hjj5nbb7\/di7lO7M22r8MXvpYAfaFafeHM2fNm8x981\/R\/9H1mwy\/93yANgghMtn2glLX7JGBKCRCjIAABCEAgNwIiXu468qL5QPuV5mv9H\/Juqj4NKKaQFE0hpXE8ZSEAAQhAQB8ByXs58u0fmu\/9\/j8PGqf5j3AEDAJG3xtMiyAAAQhUkICdOnp4y4fMx29YhYCpYB9Y0mTN6rXqvqX9EIAABLQQcPNehjZev9gszWMYERgiMFreX9oBAQhAoLIEJO9FRMzxu26uzB\/hCBgETGVfeBoOAd8JcMyJ7x5s3H73eADJeZHcl29+9uYgede9iMA0ztj7OzU733vn0AAIVJgAx5xU2PnGmO7u7mDXdTkqQJZMu3kvCJhq941KzB\/iYghAwF8CPh1z4i\/lclpuj4YY+9ofmd0nrwj2e3HzXhAw5fRb7lYRgckdOQ+EAAQSEOCzKQEkpUWs79s3P2BWr12\/LO8FAaPU8WmbxYdEWmKUhwAE8iDAZ1MelMv5DFfAPP7AlmV5LwiYcvotd6v4kMgdOQ+EAAQSEPD9s+ns2bNm+\/bt5tSpU0tau2nTJjM6Omra2toSUGhdEcvT1rh3717T19fX1AMWFhbM8PCw6e\/vb+lhwGl8n6ZsU40t4GZWIdWBrtn5BfQ3HgkBCLSIgO+fTVbADA0NtXRwbwSvsLz33nvNoUOHzNq1a40VHr29vU2LmEbsqXdPGt+nKVvvuWX7PQIGAVO2Pok9EIBAAgK+D0z1BIwVEcePHw9o2MMIT58+bXbs2LFI6ODBg8H\/bTRn\/fr1Rn7W3t5upOy2bdvM3NyccX\/u4o0TK3Lvnj17zP79+4O6jh07ZkZGRoJb16xZE4gdWcosEZZz586ZZ599NvjZ448\/vux7G4GJs0d8OTY2Zl577TXzkY98xDzwwAPmS1\/6krFtD0eD0vg+TdkE3a5URRAwCJhSdUiMgQAEkhHwfWCqJ2BkQBfhIdNJMs1kIyRCR0TJQw89FERuwtM0IjSmp6fNZz7zmUDU2AiPW587PWUF0b59+4LoS9RlBYYVRlKXXL\/3e78XCBgRNPIca0v4exEw69ati7XHbZ\/YYNsgdQqnXbt2mZ07dwZCSq40vk9TNlnPK08pBAwCpjy9EUsgAIHEBMIDk+zCKl9lvWSDNXeTtagcGBslEYHh5o64UZJbbrkliMBYwRGOlNj2h0VHXLm4n9fiaAWGFTB2qikczXHFldQnwseNDtkIz8svv7zkd66AiRNUW7duXYxK1bIVAVPWNyIHuzQ7Pwd8PAICEMiIQPizSXZiHXtiOqOnNV\/t0MauJXuV1IrARP1OBv+uri4jAsad2gkLFVfAyCDvXnbqx420JInAhKezpM4777xzMQJjp4jC0aCwgImz54033jBHjhxZkrws7T1w4EBgvp0+C7ct\/PNmxU7zXs63BiIwRGDy7XE8DQIQaAkBLRGYqCTeOCEgkY6wgKkVgQmLgijwcTkw7tTNU089ZU6cOLEoMMIRmKQCJs4e8WXc78SOe+65x9x3332LU1xp\/rBOU7YlHTPHShAwCJgcuxuPggAEWkXA94GpmRwYNwITrkfEhYgNmWb63Oc+t5gDIz8\/evTo4hSO64e4VUg2l8XWKfk4Ingkt0byb+wUUhIBE86Bce2RKSRXwNhokyzjJgcm\/o1BwCBgWvV5Sj0QgECOBLQLmFqrkFwBI8jjVve4P4+aPgqLGHeKR6aIJDokl5uvI3k6kiD8xBNPLK4WSiJgRPDE2ROOwITzg5hCin6xVAoYV73GfZ64mxbFLa+Te33\/kMjx85RHQQACORLgsylH2CV7VBrfpylbsmbWNUedgLGJT7V2UbQq2C7Dc8OD4d0fNTu\/bu+gAAQgUFoCfDaV1jWZG5bG92nKZm54ix+gRsDYkNuqVasCRLfddlvsDor1lqhFhRWTZHu32DdUBwEIQCCWgOaBCbfXJpDG92nK+sZdlYCZn58PNhSS\/QPitoBOu0W0df7AwIDp7u4O\/Cu7L8oXFwQgAIGiCGgemIpi6stz6\/l+ZmbGyJdcs7OzZnBwMNGeMb6039qpRsDYBtUTKPb3GzduNI888kiww2OSHBjXsSJm5IsLAhCAQFEE6g1iRdnFc7MnUM\/3ExMTRr7cS+MsQmUFzJkzZxaX08VtMS3Otx1lfHzcdHR0EIHJ\/t3kCRCAQAIC9QaxBFUUWqTeMuosjUu7+679w\/fkyZOLBz5a++wZSa04vTppm+v53o3ATE1NBWIGAZOUboHlkkZg3CmmWjsx1usoBTaVR0MAAhUm4Ptnk48CRg5t\/K3f+q3F\/EoZb3bv3m1eeOEFs2XLltxOrk7j+zRlfXudKheBEQeFl1nXUuOane9bZ8VeCEDgEgHfP5vqCRh3q4tNmzYFu+DK5Z7+LFEFSQGQn0WdWi37xaxcuXLJ7+yGcjZ9QOr+8pe\/vKRrhaMp9g\/jm266yXz\/+983999\/v5EVqzJ2HD58OLj3wx\/+cCBgwnu4uPvJuMcDJPl5XH9P4\/s0ZX17vyopYMSh7qFa9mRRu2mR60TNzvets2IvBCBQDQHjbnVhBYos0Aif\/mz\/IJV\/5fNbPq\/Dp1bffffdgbBwUwVkiiW8GZ7UEbVtv\/zczZ189NFHF7f1l+kjueT0azmnafPmzUsWkbjRffe8I7lHIjd33HGHift53MnYcm+acSlNWd\/er0oImPC5Gm4HkP9bdR\/eAyZtR\/HN+dgLAQj4SyA8ML372rR550flPczx8l\/uMpe9t2sReK0IjLTN3Vrf\/tH58MMPB6dQ2xSAcB3uZ\/2111675NRqt84oAVPLHrfeZ555xmzYsCGI\/FgR8vjjjwcCRoSSe7mCKOrARjvGJDmzqdE\/rBEw\/r7jTVuu2flNw6ECCECgMALhz6b5P91l5h99sDB76j34mk89aK759M5EAia8uaid5pepHhEwduv+8HSNrVymgKJOrbZCISxgkuZOynPlEhHzO7\/zO4tRHFnRagWMO\/UlZd0jDGzCr\/zcnaaK+3kc0zTjUpqy9XxYtt+ri8C0GrBm57eaFfVBAAL5EahyBMYVMOGTmq0HwrmNcREYibzX2jtM6nMjMJJDs2vXruBk6LfeeiuYurJ5lb\/xG78RHPRoT9iuNyVl22FtjpotiOpRacalNGXz672teRICpg5Hzc5vTReiFghAoAgCvn821ZqyqZcD4w78bg6je59MIbl5LnECRqInckXlQMYJC5uMa5cmxwkYiax89atfDZZef+c73wlyZeQ5dvWS5MDE\/ZwcmPpvFQIGAVO\/l1ACAhAoHQEtAkZWA7mXFQW1ViG5AiZ8arWdmqkVgZF7JFJy4cIFs3r1avPkk08uscFdIRSOwMip0iJMjh49umQvMTuF5E4HffGLXzTf+973gimv8Gopa2ec\/bU6XBrfpylbuk5exyAEDALGtz6LvRCAQMqVKADTRSCNKElT1jdKCBgEjG99FnshAAEETKX7QBpRkqasb1ARMAgY3\/os9kIAAgiYSveBNKIkTVnfoCJgEDC+9VnshQAEEDCV7gNpREmasr5BRcAgYHzrs9gLAQggYCrdB9KIkjRlfYOKgEHA+NZnsRcCEEDAVLoPpBElacr6BhUBg4Dxrc9iLwQgoEDA1DvMMUsn1zrAt95zRRDY85Zq7dVSr55mfp9GlKQp24xNRdyLgEHAFNHveCYEINAkAd8HJh8FjN2ATlxnT59u0o0N3Z7G92nKNmRMgTchYBAwBXY\/Hg0BCDRKwPeBqZ6AqbWR3blz58yzzz5rZNO78AZxdiM8G2VZuXKlOX78eIBZfidHAcgmdrKBntwrh\/nKGUvu5Z5T5P5c6jx8+LD51Kc+Zb7yla+Y\/fv3m\/b29sWjBm666aZgc7u5ubnFs47sLr1y0GPSowLq9Yk0vk9Ttt5zy\/Z7BAwCpmx9EnsgAIEEBHwfmJo5SkAOSLRb\/7tHCbjTO4Jw27Zt5u677w5OiZZyIixGR0dN1GnUUj7u7CLrDtllVy5bn5xKLTvzWmFy5syZQMDIydM7duwIDp7s7OxcPLVafi4HQdY6tiCB600a36cpm+TZZSqDgEHAlKk\/YgsEIJCQQHhguvDjGfOTH88kvDv\/Yj9\/VadZcVXn4oNrCRj33CI5bFG+FwHy8MMPB6Kgt7c3EBHhOtwIh5yFZEWE5Kq4dUYJmHoRIff8IlufFSNRp1m7kReJ3EiUZ8WKFYtRm2Y8kEaUpCnbjE1F3IuAQcAU0e94JgQg0CSB8MC08NJXzMJLE03Wmt3tbTcOmLYbP59IwEik48SJE0G0RASMnQ4SESACxj2N2k4HuZbLFNAtt9yS6DBHdwrICqMoCvagSIni2EumoCTiYk+0Dh8yac9HkvKuoGmWchpRkqZss3blfT8CBgGTd5\/jeRCAQAsIVDkC4wqYe+65x9x3330mvCKo1mGObgTGio9a4sUKEPnXnf6xomTz5s1meHh4MTIUjsiIr\/78z\/888LqcQN3s6qU0oiRN2RZ0y1yrQMAgYHLtcDwMAhBoDQHfB6ZmcmDCkQ4rLGyU5KGHHjIyhbRnz57FKZu4KaRHHnkkcEitvJQ4W22k6IEHHjBf+tKXgnpsjo2dvhI7du3aZXbu3BnkxkgS8P333x9EbRq90vg+TdlG7SnqPgQMAqaovsdzIQCBJgj4PjBZUSCrgdzLriKqtQrJFTA22mFXGtkVRLUiMHKPTD1duHDBrF692jz55JNLbLjzzjuXCBqbgyPTRTLlZC8rmHbv3m0ef\/xxI8nFBw4cCH7trpBy7RXRMz093VQibxrfpynbRHcs5FYEDAKmkI7HQyEAgeYIaB6YmiOT\/92tWh6d1PI0vk9TNunzy1IOAYOAKUtfxA4IQCAFAc0DUwoMpSiKgCnGDQgYBEwxPY+nQgACTRFAwDSFz+ub0\/g+TVnfoCBgEDC+9VnshQAEnLOQBgYGTHd3N0wqRGB2dtYMDg4GeTaykV6tCwHjWcdIs95ekrDczY7CTdXsfM\/cirkQgIBDQJYCyyA2NTUFlwoSENE6Pj4e7PSLgFHSAUS8SBZ43FkWbjPtvOXJkyfNoUOHItfmI2CUdAyaAQGFBETEyBdXawm8On\/e3HXkb8zQxi7zsQ+uam3lLapNhEs98SKP0jyGqZlCskvyVq262Nluu+22YKvpesrUnqMhuztGbS6k2fkteo+oBgIQgIAaAmfOnjeb\/+C7pv+j7zNDG6\/3vl2axzBVAmZ+fj5Yh+\/uiBjX+0TwyOZCsj5fRAwCxvv3lAZAAAIQaIrAc6+cM587+qIa8UIEpqnukP\/NUYdqRVlhTxWV8zKS5MC4iXJJQ3f5t54nQgACEIBAIwTGnvg7c+TbPwyiLhJ98flypxbTJPz61mY1ERgLPomAkcRdu52zODqJgHEdK2JGvrggAAEIQMB\/Ala8fPOzN5sPtF\/pfYMmJiaMfLlXkhVLvjW8kgJGpow2bNgQLD9LugpJsr07OjoC\/xKB8a2bYy8EIACB5QQk3+XIt38QRF60iBdppRuBkVVqImYQMB68AfUiMHHnb0jTohysOQHKA3diIgQgAIFMCGgVL2FYmsewSkZgXAcnjcBoVK+ZfCpQKQQgAIGSE7ArjWSJ9Nf6P1Rya5szDwHTHL9c746KwNQ6pwIBk6t7eBgEIACBQgloWyZdDyYCph4hxb\/X7HzFbqNpEIAABJYRsMm6WvZ4SeJizWOYuimkJA5NU0az89NwoCwEIAABXwlUJd8lyj+axzAETJ03UrPzff0wwm4IQAACSQmIeLnryItGjgfQtNIoafs1j2EIGARM0veAchCAAAS8IlC1fBciMF51z+yN1axes6fHEyAAAQgUQ6CK+S4ImGL6WmmfioAprWswDAIQgMAyAm6+i4ZjAZp1seYxjCkkppCafT+4HwIQgEApCFQ5WTfOAQiYUnTNYozQ7PxiiPJUCEAAAq0nQL5LNFPNYxgRGCIwrf8koUYIQAACORFwoy5V2t8lKV4ETFJSCstpdr5Cd9EkCECgQgREvEiy7l9+\/1wll0gncbXmMYwIDBGYJO8AZSAAAQiUigBTRsncgYBJxkllKc3OV+kwGgUBCKgmwJRROvdqHsOIwBCBSfc2UBoCEIBAQQSYMkoPHgGTnpmaOzQ7X42TaAgEIKCewHOvnDOfO\/qiIVE3nas1j2FEYIjApHsbKA0BCEAgRwJMGTUHGwHTHD+v79bsfK8dg\/EQgIB6AjZRVxpaxYMYW+FgzWMYERgiMK14R6gDAhCAQMsIEHVpGUqDgGkdS+9q0ux875yBwRCAgHoCRF1a62LNYxgRGCIwrX1bqA0CEIBAAwSIujQALcEtCJgEkLQW0ex8rT6jXRCAgF8ERLzcdeRF8+r8efPwlg+Zj9+wyq8GlNhazWMYERgiMCV+9TANAhDQTMBGXcaemDZDG7vM0MbrNTe3kLYhYArBXo6HanZ+OQhjBQQgUEUCbq4LUZfseoDmMYwIDBGY7N4caoYABCAQIkCuS75dAgGTL+9SPU2z80sFGmMgAAH1BI58+4fB6dFysa9LPu7WPIapjMCMjY2Zrq4u09fXF9lDzp49a7Zv325OnToV\/H7Tpk1mdHTUtLW1LSuv2fn5vD48BQIQqDoB9wwjjgLItzdoHsPUCRgRLwcOHDB79+6NFDALCwtmeHjY9Pb2Br+3369Zs8YMDQ0hYPJ9t3gaBCCgmICbpGuFywfar1Tc4vI1DQFTPp8ss8hGVVaturj87rbbbouNwIRvPnbsmDlx4kRkFEaz8z1wKyZCAAKeEnCTdIm6FOdEzWOYmgiMCJj5+XkjkRQ3wpKk2yBgklCiDAQgAIH6BEjSrc8ozxIImDxpN\/ms8BRRveps5GbLli2RERvr\/IGBAdPd3R1U19nZGXxxQQACEIDARQLudJFME5GkW1zPmJmZMfIl1+zsrBkcHDSTk5Omp6enOKMyeLKaCIxlk0bA2LJyb70kXpe9iBn54oIABCAAAWPc1UVMFxXfIyYmJox8uRcCpni\/1LUgqYBJIl7kYTYCMz4+bjo6OojA1PUABSAAgaoQcI8AQLiUx+tuBGZqaioQMwiY8vgn1pIkAqbeyiO3cs3zhx64ExMhAIESEmB1UQmdEmOS5jGsklNIstR6bm4udtoIAePPy4mlEIBAfgTCeS4cAZAf+0afhIBplFwB90VFYOzP+vv7zbp165ZsYmdNXL9+vTl48KBpb29fYrVm5xfgHh4JAQh4SCAsXJgu8seJmscwdRGYVncrzc5vNSvqgwAEdBFwl0RLyxAu\/vlX8xiGgKnTHzU7379XEYshAIG8CLCyKC\/S2T5H8xiGgEHAZPv2UDsEIOAVgedeORccuPjq\/HkiLl55LtpYBIwCJzbaBM3Ob5QJ90EAAvoIhJdE9390teHcIv\/9rHkMIwJDBMb\/N5QWQAACDROQqaLJ\/\/GDIOLysQ+uMkMbr0e4NEyzfDciYMrnk9ws0uz83CDyIAhAoFQEJNryl98\/h3AplVeyMUbzGEYEhghMNm8NtUIAAqUjYIWL5LjIJauKmCoqnZtaahACpqU4\/apMs\/P98gTWQgACjRJAuDRKzv\/7NI9hRGCIwPj\/htICCEAgkgD7uNAxEDAV7gOanV9ht9J0CKgmwM65qt2bqnGaxzAiMERgUr0MFIYABMpLwK4okgRdWQItK4okz4WrugQQMNX1vdHs\/Aq7laZDQA2BqBVFkpj78RtWqWkjDWmcgOYxjAgMEZjG3wzuhAAECiMQNU3EiqLC3FHaByNgSuua7A3T7Pzs6fEECECglQTcaIudJuKAxVYS1leX5jGMCAwRGH1vLC2CgDIC4dVEsmMu00TKnJxRcxAwGYH1oVrNzveBPzZCoKoE4qItTBNVtUc01m7NYxgRGCIwjb0V3AUBCGRCgGhLJlgrWykCprKuN6xCqrDvaToE8iJAtCUv0tV7DgKmej5fbLFm51fYrTQdAoUTsKLluVfmjezfIvu2kNtSuFvUGaB5DGMKiSkkdS8sDYJAmQmIcJHDFF3R8vEbrmHDuTI7zWPbEDAeO69Z0zU7v1k23A8BCCQj4Oa1yP8l2sLy52TsKNUcAc1jGBEYIjDNvR3cDQEIRBKIEy0f++A17JJLn8mNAAImN9Tle5Bm55ePNhZBwG8CiBa\/\/afRes1jWKkiMMeOHTMjIyNL+tDevXtNX19fYf1Ks\/MLg8qDIaCEQDgRV5plp4eItChxsufN0DyGlUbAiHg5evSoOXjwoGlvbw+6zNmzZ8327dvNli1bMhMxY2NjpqurK7Z+zc73\/L3EfAgUQiC85BnRUogbeGhCAprHsFIIGCtUhoaGTE9PzxK3CHwRGa6wSei3usWk3gMHDphaUR7Nzq8LiAIQgEBAwE4NPffKOWPPIJIlz6weooOUnYDmMaySAsYKplWrLh43f9tttxGBKftbiH0QyJFA3NQQ+7Tk6AQe1RICCJiWYKxdSZ5TSCJg5ufnzZo1a8zw8LDp7e2tK2AGBgZMd3d30IjOzs7giwsCENBBwBUs9v+Sy\/L+a64MVgyRz6LDz1VpxczMjJEvuWZnZ83g4KCZnJxcNsPhO49SRGAsxLyTeBcWFhIJmHU\/6m\/Kzyuuqi12VrQt\/\/3Ph+5ZcVXHMhvcesN12PvrPbuphnEzBDwlICLlolCZN3ZaSJrCbrieOhSzlxCYmJgw8uVeCBhlnSSpgPmTP7zbfPL22830+avNc98\/t4yCfOh9\/IOrzMduuDglZa8LP76ogGtdF348u+zXP4m478LCpbqS1Bv1zDjBc0nsXBJJtqwVRgiiep7k92UmUE+wkMtSZu9hW1oCbgRmamoqEDMImLQUS14+qYDZunXrEufbEPOZswtL\/nqzf8HJv7LLplx5hJ7DgsYVQHHCxwqnqLJJBFJYDMWJIFcAEQ0q+QuhyLyoKSH7ftrkWzs9pKjZNAUCywiQA5NBp3BXHq1bty5YLn3q1KnIJ61fvz6TVUiNCpgoI92\/8OT3bljafnC6H5h5CJtm3eYKGSt0rCCyvwsLIfl9PQHkRnfCwkd+h+hp1nPVut++e6\/OnzdyMKIVL2HBYv+oqBYdWlt1AggYpT2glQImDpEbrbn0\/\/OLxd1EwbwiNnm5MyyAwuJGxI8rjOKET1jw2HwgxE5enizXc2pFV+SPhGBKl8MRy+U0rCmMAAImY\/RF7QOTh4BJKmzcvxrtX44XBc2q4AP5A+1tiysiMnZH4dWLkIkSNq7gefeN5yPtjBM7l117cX8hifgwlVW4ixMZIO+EXLLvSlRkxb4fIlaYDkqElEIVJICAydjpRQmYJM3K2\/n1pqKsuAlPR10UOVcmaZKqMjZqI4LHjfBYsRM3peUKncuvu7g83kZ0EDn5d5F600CIlfx9whN1EMh7DMuTWqHLqKOWTUc1\/s477zSyS28RV1mc7\/41KsnDl4TO0lVRdkrKCpoqRW5q9Y84ofPO61PBbVERHStyLr+2x9hpKxvJufy6pTtGF9E3fXxmLaHiinM7DURkxUcvY3OZCJRlDMuCSaECxjaoVgQmi0anqdMH57uDQj1xI213BwUEzqXeUEvkREVyoqI4CJyL2+7LJVM\/coWnf6xQcaMq9v9VjCKm+TyiLATSEvBhDEvbJlu+FAKmUePzuM9354cjNxcHlIsDix1gLEc7eNhESPu9rJiSS3YkrfplRc47bzy\/uNoqLooTF8HRMEWVRjS7ibXhvlX1\/kT7IZA1Ad\/HsFp8EDB1eo9m50vTowRO3PSUK3TCA5GN5FQ1F8eyiRM4cREcWTJuc3AkelOGqSnbJ+RfWZrsRvUufn9pFZ2Nprj9gahe1kMS9UMgOQHNY1hpBMzp06fNtm3bzNzc3DLPZLUPTJIuoNn5SdpfT+REDWj1Ijp2gLODX5WmDewKK5tzI9GbvMRNWJhcjMbNB+66JFaWipOwQLkYibsYkSM\/JekbRDkIFEcUJagUAAAeRElEQVRA8xhWCgHjLmfevHlzcD5Rf3+\/sRvcSQJvT08xSZOand\/KVyrqr3Y7MNrf1RI7dqC0A2OUuBHhE\/69e18r21NEXUnEjV0p5UZt\/uoHV5vZ89ctRkpc0ZlGYIbFSdWjaUX0AZ4JgVYT0DyGlULAhJN4x8bGTFdXV3BCtMA\/cuSIGR0dNW1tFwewPC\/Nzs+TY\/hZcYInLHrk+3rCx63bjeZIhMBe4ShPVNTHCqR6XNx6w2XF1lqXTMe4V3g6Jvy9rU+WiXevesl0Xvm66bjy9eBf+d69Zs5fZ6bO3RiIGfm\/5Npcdm33kn2ErOCrUtSrnj\/5PQQ0E9A8hpVSwMjy6unp6WDptMAXQXPw4EHT3t6eez\/T7PzcYbbgge4A70Z2LgqfeHHg3hcWGWHR0AIza1YRFg9hQeT+3v7fFVe2vPxOhIyIG3dKKrwkXKI2dim45NloSCLO2kfUDwEtBDSPYaUQMNJRRKTIFRYtTz31lDlx4gQRGC1vk2ftaETclCG6YZOJ33716wFxybVxhY1dIXXF+28Pfl+WBGLPugfmQqD0BBAwObgovK2\/CJoDBw6YNWvWmEOHDpm1a9fmYMXyR2h2fiFAeWihBKLybMLCxq6MupRvU0z+WaGgeDgElBDQPIaVIgJjxYsk7haVrBvXVzU7X8n7STNaQMAVNnIMg43c2KpFzBCtaQFoqoBAzgQ0j2GlEDDsxJtzj+ZxEEhAIMk0lJtbU4Y9bBI0iyIQqBQBBEwO7pbE3SJzXYjA5OBkHqGCgAgbuxNxVG4NokaFm2mEEgIImIwdaSMwp06dinwSG9ll7ACqh0CTBOqJGqafmgTM7RBokAACpkFwGm7T7HwN\/qEN5SWQRNRczK35ZHkbgWUQ8JyA5jGsFDkwZe4fmp1fZu7YppOAiBqbILzw4jFjVvwwaKi7wzBLunX6nlYVQ0DzGIaAqdOnNDu\/mNeJp0LgEgF39VNUPo1MPRGlocdAoHECmscwBAwCpvE3gzshkAEBN0rjiprwjsKsesoAPlWqI4CAUefS5A3S7PzkFCgJgeIIxEVpmHYqzic82R8CmscwIjBEYPx5E7EUAj8jYKM0cdNO5NHQVSBwkQACpsI9QbPzK+xWmq6MwJLk4JcmFltndxBG0ChzOM1JTEDzGEYEhghM4heBghDwhQCCxhdPYWfWBBAwWRNuUf2ym+\/IyEhQ2969e01fX19szW7ZTZs2xZ52rdn5LcJONRAoPYEoQePm0MheNPaE7tI3BgMhkIKA5jFMTQTm9OnTZseOHWbfvn2Ba+3\/o06xFofKadcHDx40bW1tZnh4ODj1emhoaFm30Oz8FO8ARSGgikBUDo17YGXbjZ9X1V4aU10CmscwNQImfJaSCJSurq7IKEy4bK1zmDQ7v7qvNC2HwFICtQQN+TP0Fp8JaB7D1AgYESxy2ShK+Hu3A0ZFYHp7eyPFjmbn+\/xSYjsEsiLgLtt++9XHjD2V20ZomG7Kijz1ZkFA8ximSsC4EReJqkxPT0dOC0knkSmnbdu2mbm5OTM5OWl6enoi+451\/sDAgOnu7g7KdHZ2Bl9cEICAfgI2OnPhx7OLxyCwukm\/331u4czMjJEvuWZnZ83g4GDNcc7XtlZSwIi4OXr0aJAD097eHuTDuNGbcLRm69atS\/wrYka+uCAAgeoRsIImfJaTHHtAdKZ6\/aGMLZ6YmDDy5V61\/lAvYxuS2KRKwLgiJE6ULCwsBEm77pSRmwAcTvq1EZjx8XHT0dFBBCZJr6IMBCpCwIqZd1\/7P+ads38WtJroTEWcX+JmuhGYqampQMwgYErssPCUUVwSb6MCRqPzS+xOTIOAlwSsoLG5M6xs8tKNqowmB8YDd6ZZRh01hSS5MKOjo8Gy6qgpJASMB50AEyFQIgK1VjaxTLtEjlJuCgLGEwfHbWRnoy79\/f2LyboSoTlw4EDQMjay88TBmAkBTwnUmmoib8ZTp3piNgLGE0dlYaZm52fBizohAIH6BN55\/Xnz7hvPm4WfndvEEu36zCjRGAHNY5iaJN7GXFv\/Ls3Or996SkAAAlkTWFzVFBIzbKCXNflq1K95DEPA1OnDmp1fjdeXVkLAHwJxScCIGX98WDZLNY9hCBgETNneN+yBAASMCXYAfvvVr5u3Tu43K66+bHF5NjkzdI80BBAwaWgpK6vZ+cpcRXMgoJYAYkatazNvmOYxjAgMEZjMXyAeAAEItI5AXM4MS7Nbx1hTTQgYTd5M2RbNzk+JguIQgEDJCCBmSuaQEpqjeQwjAkMEpoSvHCZBAAJpCfz9f\/myWbHyjcUjDSTx94oPXDyfiau6BBAw1fW90ez8CruVpkNALYHwYZPsMaPW1YkapnkMIwJDBCbRS0AhCEDAPwJxy7LJl\/HPl41ajIBplJyC+zQ7X4F7aAIEIJCQQFS+TNuNA0wxJeTnazHNYxgRGCIwvr6X2A0BCDRA4N3Xps3Ci0fN22ceMytWvs7+Mg0w9OkWBIxP3mqxrZqd32JUVAcBCHhGIG5\/GaaYPHNkDXM1j2FEYIjA6HlTaQkEINAwgbe+\/Sfmwtsng1VMJP42jLF0NyJgSueS\/AzS7Pz8KPIkCEDAFwJEZXzxVDI7NY9hRGCIwCR7CygFAQhUjoDsLfP2q4+Zy1afIyrjqfcRMJ46rhVma3Z+K\/hQBwQgoJ9AVFSGFUx++F3zGEYEhgiMH28hVkIAAoUTkBVMkisTjsqQ9Fu4a2INQMCU1zeZW6bZ+ZnD4wEQgIBaAhKVees7+xePLrARGUkA5ioPAc1jGBEYIjDledOwBAIQ8I6ARGVk5dJbJ\/ebFVdfFuTKML1UHjciYMrji9wt0ez83GHyQAhAQC0BO7208NKE+YUPXr2Y9Mv0UrEu1zyGEYEhAlPs28XTIQABdQSiVi8hZIpxMwKmGO6leKpm55cCMEZAAAJqCQTHFrw0Yd4+8\/Vgeok8mfxdrXkMIwJDBCb\/N4onQgAClSIQzpO57Noec\/XN48E0E1e2BBAw2fJtWe3Hjh0zIyMjQX179+41fX19sXVbp0qB9evXm4MHD5r29vZl5TU7v2XgqQgCEIBAAgJReTIk\/CYA10QRzWOYmgjM6dOnzY4dO8y+ffsCV9v\/r127dpnrpey2bdvMQw89ZHp6eowInxMnTpjR0VHT1ta2pLxm5zfxTnArBCAAgYYJWCEjK5euXL+KhN+GSda\/UfMYpkbAhEXI2NiY6erqiozCSNnp6WkzNDRU1\/uanV+38RSAAAQgkCEBhEyGcH9WteYxTI2AEcEilxUl4e9tN1lYWDDDw8Omt7e35hSTLW+dPzAwYLq7u4Mfd3Z2Bl9cEIAABCDQGgKycomITGtYzszMGPmSa3Z21gwODprJyclgxkHTpUrAuBGXuCiLFTAbN240jzzyiDl16lSiHBjX6SJm5IsLAhCAAARaR4CITGtYTkxMGPlyLwRMa9hmUkt4yqiegDlz5sxi4q7cOzc3VzMHZnx83HR0dBCBycR7VAoBCEDgEgGETHO9wY3ATE1NBWIGAdMc00zvbmYKyU0ADif9ap4\/zNQhVA4BCECgSQIImSYBGmM0j2FqppDCEZdaSbzh34mA2bNnj9m\/f\/+ypdSand\/8q0ENEIAABLInIELmzW8dNgsvfYVVSylxax7D1AiYNMuoxaEiYuzeL3EJv9JPNDs\/5XtAcQhAAAKFEoiKyLCPTG2XaB7D1AgYcWHcRnY2cbe\/v38xC9vdyG7Tpk2R+S8ImEI\/q3g4BCAAgUgCImRe\/+PPmAvn\/3oxIvOem8bN5dfpWmXTCvcjYFpB0dM6NDvfU5dgNgQgAIGAwMILT5s3Dn\/GrPjF1xEyMX1C8ximKgKTxTut2flZ8KJOCEAAAnkTmP\/TXebv\/+uXzS\/8ytWBkOGspUse0DyGIWDqvGmanZ\/3hwzPgwAEIJAlAStkruq9zlz2T64MTr9uu\/HzWT6y9HVrHsMQMAiY0r+AGAgBCEAgKQG7YunNvxg1gZD55S5zxftvr6yQQcAk7TkKy2l2vkJ30SQIQAACAQERMnM7P2EuW31uMT+miiuWNI9hRGCIwPBxBwEIQEAtgXB+zIqrOk2VViwhYNR27foN0+z8+q2nBAQgAAH\/CdhpJTfR94r3f9K85+Zx\/xtX4T\/CicBU2Pnq31waCAEIQMAhYKeVLizMXMyPqUCir+Y\/whEwCBg+4CAAAQhUioBMK80\/+mAgYK7+zfcZzdNKCJhKde2ljdXs\/Aq7laZDAAIVJyDRmB99bZt590fTi4m+sn\/Myo8dUUVG8xhGBIYIjKqXlcZAAAIQSEPgzaf\/2Lz2tW3mF66\/waz87X9m3n3jeVX7xyBg0vQGZWU1O1+Zq2gOBCAAgYYIuNGYKz\/8q+ayjr9RM62keQwjAkMEpqEXnpsgAAEIaCNgozFX\/tqt5pr+3zYLL00Y31crIWC09dIU7dHs\/BQYKAoBCECgEgTcaMx7un\/dXH7DPxpZteTrbr6axzAiMERgKvGhRCMhAAEIpCHgRmOu3X6Xeeu7g14eEomASeN1ZWU1O1+Zq2gOBCAAgZYScKMx1277d+Yn\/+8\/e5fkq3kMIwJDBKalLzyVQQACENBGwEZjrvnUg+bqf7HRvHmi35skXwSMtt6Yoj2anZ8CA0UhAAEIVJqAjcYIhPft+A\/m7Ve\/HiT5ygGRbTd+vrRsNI9hRGCIwJT2xcMwCEAAAmUjIHvGLLzwtHnvXYeCnXzf+t5gYGJZD4hEwJStB+Voj2bn54iRR0EAAhBQQ8CdUrrm0zvNwktfKW00RvMYRgSGCIyaDxUaAgEIQCAvAu6U0poHv2Uu\/HjG\/ON3B4Ml12WKxiBg8uoRJXyOZueXEDcmQQACEPCKgDul1PZrt5YuGqN5DCMCQwTGqw8LjIUABCBQNgLhKaV3Xn9+MTdmZe+RYMVSURcCpijyKZ977NgxMzIyEty1d+9e09fXV7eG06dPmx07dph9+\/aZtWvXLiuv2fl14VAAAhCAAAQSEZAppbmdnzAShZEEX7nKkBujeQxTE4FxhYh0nFqixPbGhYUFMzw8bE6ePGkOHTqEgEn0mlIIAhCAAASiCITzYqSMRGNk35jLru0xKz92JHdwCJjckad\/oERfTpw4YUZHR01bW5sZGxszXV1dNaMw4lgpJxcRmPTMuQMCEIAABJYTCOfFFJngi4DxoIdaITI0NBRYG\/4+3ISzZ8+aXbt2mf7+\/qAsAsYDJ2MiBCAAAU8I2LyY1Q9+K5hWkquIKSUEjAcdJhxxkYjM9PS0sYIm3AT5vVy33HJLohyYgYEB093dHdzT2dkZfHFBAAIQgAAE4gjIhnc\/ePATRo4gkP1i5JJozD\/IUQRtnZlNKc3MzBj5kmt2dtYMDg6ayclJ09PTo8pZanJg0ggYyZc5fPiwuf\/++wMnJ0nidb0uYka+uCAAAQhAAAK1CEQl90r5f\/jL\/sz2jJmYmDDy5V4ImBL30zRTSFJ2w4YNgRpNugppfHzcdHR0EIEpcR\/ANAhAAAJlJBCV3Ct2ynlKb313sOXnKbkRmKmpqUDMIGDK2DN+ZlN4yiguiVdyX7Zv325OnTq1rDVRDtY8f1hid2IaBCAAAXUE5h78RNAm2bnXXjbBV77PYpWS5jFMzRRSI8uopcMkjcBoVK\/qPh1oEAQgAIGSE4gSMWKyJPi+\/epjLT+GAAFT8g5hzYvbyM7u9yIrjsJJTAgYT5yLmRCAAASUELDLrD\/wB3+3pEV2z5i2GwdM242fb0lrETAtwehnJZqd76dHsBoCEICA\/wTiREyrp5Q0j2FqppCy6s6anZ8VM+qFAAQgAIH6BGSvmPMvPLN49IB7hyT3vvPG801PKWkewxAwdfqYZufXf70oAQEIQAACWRIQETP\/p7tMeDpJnmlXKb3n5nFzxfs\/2ZAZmscwBAwCpqGXgpsgAAEIQKA1BETAiJCJEjE2L0YEjAiZtBcCJi0xReU1O1+Rm2gKBCAAAa8JSE7MO69NL1libRvUTF6M5jGMCAwRGK9feoyHAAQgoIVALREjbWwkLwYBo6V3NNAOzc5vAAe3QAACEIBAhgTqiRibF\/OLvUfM5dfVP9tI8xhGBIYITIavIlVDAAIQgEBaAnGb3dl60uwXg4BJS19Rec3OV+QmmgIBCEBAFYF6IiZpXozmMYwIDBEYVS89jYEABCCghYCImLZ\/equ55tM7Y5tkT7Ve9evPRpZBwGjpDQ20Q7PzG8DBLRCAAAQgkBMBOcX6zGevN6sf\/JZp+7VbY59aK7lX8xhGBIYITE6vIo+BAAQgAIG0BBZeeNr84MFP1BUxkty78NLEsp17ETBpiSsqr9n5itxEUyAAAQioJSAiRlYnrdn1LXPZe7ti22mTe90VSprHMCIwRGDUvvQ0DAIQgIAWArV263XbGF6hhIDR0gMaaIdm5zeAg1sgAAEIQKAgApLUe\/l7uyIPf3RNsiuUfv6qTvO\/3r7dbN261UxOTpqenvr7xhTUtIYeSwSGCExDHYebIAABCEAgXwJJk3rFKitiZmdnzL8cnkXA5OuqcjyNCEw5\/IAVEIAABCBgTNKkXstKllmfOf28+cdfPUIEpmodCAFTNY\/TXghAAALlJiD5MAv\/++nIgx\/DlssY9vS\/7ze3\/i4CptxezcA6BEwGUKkSAhCAAASaIpA0H0bzGEYOTJ0upNn5Tb093AwBCEAAAoURSJoPo3kMQ8AgYAp7AXkwBCAAAQg0TuDNp\/842B\/mVx79aWwlCJjG+Xp\/p2bne+8cGgABCECg4gTqHfqoeQwjAkMEpuKvP82HAAQg4C+BelNJCBh\/fdu05Zqd3zQcKoAABCAAgcIJ1JpK0jyGqYrAHDt2zIyMjASdae\/evaavry+yY509e9Zs377dnDp1Kvj9pk2bzOjoqGlra1tWXrPzC3\/rMAACEIAABFpCIG5VkuYxTI2AOX36tNmxY4fZt29f0Bns\/9euXbukcywsLJjh4WHT29sbCBz7\/Zo1a8zQ0BACpiWvEpVAAAIQgECeBOI2uEPA5OmFBp8l0ZcTJ04sRlLGxsZMV1dXbBTGfUz4Xvd3mp3fIGpugwAEIACBEhKISujVPIapicCIYJHLRlHC39fqawiYEr6JmAQBCEAAAqkIRCX0ImBSISymcDjiIqJkeno6clrItdDmw2zZsiUyWmOdPzAwYLq7u4NbOzs7gy8uCEAAAhCAQJkIyDEDktQ797tHArNmZ2fN4OAghzmWyUlhWxoRMDb\/Reqql8TrPk\/EjHxxQQACEIAABMpG4G8\/9XNm3\/R7zZNvXL1o2uTkJIc5ls1R1p60U0hJxIvUbSMw4+PjpqOjgwhMWTsAdkEAAhCAQEBAojAv\/scvmxUjz5qpqSkzMTFBBKbMfSM8ZVQribfeyiO3nZrnD8vsT2yDAAQgAIHGCUgU5r13HTIvXPmrZuvWrQiYxlFmf2fSZdRiiYibubm52GkjBEz2\/uIJEIAABCCQHQGJwsw\/+qD50T1\/hYDJDnPrao7byM5GXPr7+826deuWbGJnn75+\/Xpz8OBB097evsQgIjCt8w81QQACEIBAfgQkCvO3H\/m35t\/84RNEYPLDXp4n+SRgZmZmzGOPPWZuv\/12L1ZJYW+2\/Ry+8LUE6AvV7AtyUvX0d58zff\/9JwiYbLtAOWv3ScD4ZKt4G3uz7fPwha8lQF+oZl+wu\/Pe8\/JqM3TgcVYhZdsNyle7Ty++T7YiYLLv6\/SHbBn7xNcnW\/lsaG2\/ld15v\/EXf21u3vMUAqa1aMtfW9RGdmW12m5Y5G66V1ZbxS7szdY78IWvJUBfqG5f+KVX\/pu54hu\/b1bc96zpuvnj2YLIuXY1RwlkxU3mjmUXQ1lLzwUBCEAAAhDwjcCffPhV86F\/9fvmmk\/v9M30mvYiYBK4U0SMfHFBAAIQgAAEfCOg9fgbBIxvPRF7IQABCEAAAhAwCBg6AQQgAAEIQAAC3hFAwHjnMgyGAAQgAAEIQAABQx+AAAQgAAEIQMA7AggY71yGwRCAAAQgAAEIIGDoAxCAAAQgAAEIeEcAAVPDZXJq9YEDB4ISk5OTpdnFUE7e3rZtW3Ci9qZNmxKfqn3kyJFEZVvdi9PY6zJfs2aNOXTokFm7dm2rTapZXxp73QNE4w4EzdL4NLZaO+zhpr29vaavry9L85bVncZel61UVHa+Z8+eXXJQbBGfGUn5htlaR+3duzfXPpHUXrHPLevDZ4PdBLWovpvri13QwxAwMeCl88lgKidUv\/zyy4v\/D59Wnbff3MFn8+bNZnh42NQbiOyLlFTstLJNaeyVD9UTJ04siiz5\/ujRo5GnhLfSRreuNPa6fUT6hfQXEZWjo6Omra0tKxMX601jq2uMHbzyHqzS2is8u7q6ch1QG+0Ltm0ysA4NDQWD7Y4dO8y+fftyE+Bp+bptDfflzDuvMSaNvVYcCtuenh5T9s8GK7YeeuihwF7hW9QfkHn4sqhnIGBiyMuHp1zywtgXrb+\/v\/AoTPiDsd6LIe04fvy4ufXWW82bb76Z2+Bqsaa113VHEYNAM\/bmPQg0YqsMBPfcc485d+6c2bJlS67iII29ZXjn0tgrZffs2WP2799vivojJ4297nsWFgd5DUZp7A2XLftnQ\/iPMenPu3fvNnfccUdugjYvPxb5HARMBP1wiL3IkHvYvPAgWW\/QlN\/bv1jc6EZenS6tvUULmGbsdUVvHnwbsVVs\/OhHP2q+8Y1v1I3ctboNaewtalCtFZWo9a6FB6xWs0tSXxq+4Yhc2T8boiIweduchm+UgJFoeRn+CE7Sl3wpg4CpIWDczlZ0ONuaGY64JP3Lr6gP2EbtlfbmPSUjz2zEXjtFl\/e8fFpbpa8cPnzYfOELXzC7du0qRMC4YfRafdfNd7B9P++ckjR85f2anp4OTC0qby6NvZZpkUIxrb32D0mJKN95551BdDzPK4294dO\/7fd5T9vmyaeIZyFgEDCZ9rs0L334r8KvfvWruSfxNmqvFT\/33ntvbjansdUNYcu5KElyp1rdMdLYK2VdluHvW21bVH1p7LV5RVZkld1e9w8im+uX99RXGr5ROSV5253GXuHrJkpv3bo1mMKvl6+YR7\/W9AwETA0BYzubz1NItnlFRmDcD5p6U172xS9CvFgRktZeyzjvfpImpC1ln3nmmSU5XXl\/mKaxN\/xa5s02bV+ImzLIk3EjfG3kKO9oRlX4liHSpUmwhNuCgInxrjtlVIaEQmtmOOwe\/qsgrrMWJWDS2lvE6gKXWVp73XvzDsensdVdnu7anGcoPo29cQImzxyCNPaG38MiPjPS2Ct8ixCFjb5rZRCIafmG25r3qjTNwsW2DQET42X3rxnfl1HbqEbeSW\/hD8l6y76LCLvX+ks\/ib1utCZv8ZVmGarbzqIGrjT2+rakPixek0QaWz3ApOErz7Yr0u67775CVsaksTdqCinP6dq0n2Wu2JEtFWTK1i6xb7Xfq1wfAqaG933cyC4uJFxUBEbw1tqsyrU3LkqQd\/JmUnutMBwZGQl6Ud5JvGnYlkHApLXXzSEogm1ae92N7Hywt4ilyOGP2zTvmk2ELepdS9sf3P5bxB5cVRA2CJgqeJk2QgACEIAABJQRQMAocyjNgQAEIAABCFSBAAKmCl6mjRCAAAQgAAFlBBAwyhxKcyAAAQhAAAJVIICAqYKXaSMEIAABCEBAGQEEjDKH0hwIQAACEIBAFQggYKrgZdoIAQhAAAIQUEYAAaPMoTQHAj4SkD1U5IDJnTt3mrzP5PGRFzZDAALGIGDoBRCAQOEE3LOaCjcGAyAAAS8IIGC8cBNGQkA3AdmFecOGDaanp0d3Q2kdBCDQMgIImJahpCIIQCBMIO7cJfesIDkrZvfu3eaOO+4IzuSx9xw\/fjyorqht+fEmBCBQbgIImHL7B+sg4D2BqEMuJeIi19DQUHBW1uHDh839998f\/EwOvpNrdHTUiLjJ+5BM74HTAAhUhAACpiKOppkQKIpA+KTm8PciUOTq6+sLxMyOHTvMvn37CjkhuShGPBcCEEhPAAGTnhl3QAACKQm4ERd3+khWHLn5L+HfpXwMxSEAgQoRQMBUyNk0FQJFEXCFySOPPGK6urqCiEt4+TQCpigP8VwI+EcAAeOfz7AYAt4RsNNGmzdvNt\/85jcXp4jCy6eZQvLOtRgMgcIIIGAKQ8+DIVAtApLrMjIyYjZt2rQkQVcoSDRGLrsCSVYeSYKvXIiaavUTWguBpAQQMElJUQ4CEGiKgAiRbdu2mbvvvjsQLCJW3OXTtnKWUTeFmZshUBkCCJjKuJqGQgACEIAABPQQQMDo8SUtgQAEIAABCFSGAAKmMq6moRCAAAQgAAE9BP4\/xTXAKVQ3C6AAAAAASUVORK5CYII=","height":337,"width":560}}
%---
