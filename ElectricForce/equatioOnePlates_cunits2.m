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
taylorOrder = 9;

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


%psin = simplify(sqrt(dot(r_vec_l1,r_vec_l1))/dot(vec_l1,r_vec_l1))
psin = simplify(dot(r_vec_l1,r_vec_l1)/dot(vec_l1,r_vec_l1)) %[output:1862f7b5]

%%
%[text] ## The Relative Net Force of the Two Points
% Force of the two point in plate and moving particle
dforce = (psin*vec_l1 + psin*vec_l2 )/(2*R^2); % net Force in the x direction first points.

dforce = simplify(expand(dforce),'Steps',100) % Two charges %[output:3d0726da]
dforce_x = dforce(1) %[output:69e50e72]
%dforce_x = simplify(taylor(expand(subs(dforce_x/2,v,v+v_a) + dforce_x/2),v_a,order=2))
%%
%[text] ## Check at v = 0, r =0
df_v_zero = simplify(subs(dforce_x,[v,r],[0,0]),'Steps',100) % Due to the two point charges %[output:4bf04a0c]
df_v_one = simplify(subs(dforce_x,[v,r],[v,0]),'Steps',100) % Due to the two point charges %[output:7b7123d6]
%%
%[text] ## Integrate all the angles
%dforce_xT = expand(taylor(dforce_x,v,Order=taylorOrder));
dforce_xT = dforce_x;

dcircle = int(dforce_xT,a,0,pi); % From zero to pi i.e., the full circle
%%
dcircle = simplify(dcircle,steps=100); % Net force of the two full circles with charge density
df_v_zero = simplify(subs(dcircle,[v,r],[v,0]),'Steps',100) % Net force of the full zero radius circle (point) %[output:3b43bbb2]
%[text] ## 
%%
%[text] ## Integrate for all values of r
syms R_o real
assumeAlso (R_o>0)
dtotf = r*dcircle/pi  %[output:34b522cd]
dtotf = simplify(expand(taylor(dtotf,v,Order=taylorOrder))); % Taylor
%dtotf = simplify(subs(dtotf,v,0)) %static
totInt = int(dtotf,r,0,R_o); % Add all the circles
%%
totf = simplify(totInt,'Criterion','preferReal','Steps',100) % The force on a charge due a charge surface density %[output:4739c8cb]

%%
%[text] ## Set the Plates Radius to Infinity
sympref('PolynomialDisplayStyle', 'ascend');
totfInf = simplify(limit(totf,R_o,Inf)) % Net force by plate charge density %[output:38ff51e6]
%totfInf = subs(totfInf,v,-v)

%%
figure(1)
plotrange = [0,0.99];
fplot(totfInf,plotrange) %[output:0d212725]
hold on %[output:0d212725]
relEq = sqrt(1-v^2) %[output:079a3f68]
relEqApx = taylor(relEq,v,Order=taylorOrder) %[output:50ce538e]
fplot(relEq,plotrange) %[output:0d212725]
fplot(relEqApx,plotrange) %[output:0d212725]
xlabel('v/c')  %[output:0d212725]
ylabel('ratio')  %[output:0d212725]
legend({'Force Carriers',"Lorentz Mass","Lorentz Apx"},'Location','northeast') %[output:0d212725]
hold off %[output:0d212725]



hold off %[output:0d212725]
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
%[output:1862f7b5]
%   data: {"dataType":"symbolic","outputData":{"name":"psin","value":"\\frac{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 -2\\,d\\,v\\,\\sqrt{r^2 +d^2 }}{d^2 +r^2 -d\\,v\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:3d0726da]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n-\\frac{d^4 \\,{\\left(v-v^3 \\right)}-d^3 \\,{\\left(\\sigma_1 -\\sigma_1 \\,v^2 \\right)}+r^2 \\,{\\left(d^2 \\,{\\left(v-v^3 \\right)}-d\\,{\\left(\\sigma_1 +\\sigma_1 \\,v^2 \\right)}\\right)}}{{{\\left(r^2 +d^2 \\right)}}^2 \\,{\\left(r^2 +d^2 -d^2 \\,v^2 \\right)}} & 0 & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\sqrt{r^2 +d^2 }\n\\end{array}"}}
%---
%[output:69e50e72]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce_x","value":"\\begin{array}{l}\n-\\frac{d^4 \\,{\\left(v-v^3 \\right)}-d^3 \\,{\\left(\\sigma_1 -\\sigma_1 \\,v^2 \\right)}+r^2 \\,{\\left(d^2 \\,{\\left(v-v^3 \\right)}-d\\,{\\left(\\sigma_1 +\\sigma_1 \\,v^2 \\right)}\\right)}}{{{\\left(r^2 +d^2 \\right)}}^2 \\,{\\left(r^2 +d^2 -d^2 \\,v^2 \\right)}}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\sqrt{r^2 +d^2 }\n\\end{array}"}}
%---
%[output:4bf04a0c]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\frac{1}{d^2 }"}}
%---
%[output:7b7123d6]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_one","value":"-\\frac{-1+v}{d^2 }"}}
%---
%[output:3b43bbb2]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"-\\frac{\\pi \\,{\\left(-1+v\\right)}}{d^2 }"}}
%---
%[output:34b522cd]
%   data: {"dataType":"symbolic","outputData":{"name":"dtotf","value":"\\begin{array}{l}\n\\frac{r\\,{\\left(\\pi \\,d\\,r^2 \\,{\\left(\\sigma_1 -d\\,v+v^2 \\,\\sigma_1 +d\\,v^3 \\right)}+\\pi \\,d^3 \\,{\\left(-1+v^2 \\right)}\\,{\\left(d\\,v-\\sigma_1 \\right)}\\right)}}{\\pi \\,{{\\left(r^2 +d^2 \\right)}}^2 \\,{\\left(r^2 +d^2 -d^2 \\,v^2 \\right)}}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\sqrt{r^2 +d^2 }\n\\end{array}"}}
%---
%[output:4739c8cb]
%   data: {"dataType":"symbolic","outputData":{"name":"totf","value":"\\begin{array}{l}\n\\frac{2\\,v^2 }{3}-\\frac{v}{2}+\\frac{v^3 }{4}+\\frac{2\\,v^4 }{15}+\\frac{v^5 }{12}+\\frac{2\\,v^6 }{35}+\\frac{v^7 }{24}+\\frac{2\\,v^8 }{63}-\\frac{d}{\\sqrt{\\sigma_4 }}-\\frac{d^2 \\,v^3 }{2\\,\\sigma_4 }+\\frac{d^4 \\,v^3 }{4\\,{\\sigma_4 }^2 }+\\frac{d^3 \\,v^2 }{\\sigma_3 }-\\frac{d^4 \\,v^5 }{4\\,{\\sigma_4 }^2 }-\\frac{d^3 \\,v^4 }{\\sigma_3 }+\\frac{d^6 \\,v^5 }{6\\,{\\sigma_4 }^3 }+\\frac{d^5 \\,v^4 }{\\sigma_2 }-\\frac{d^6 \\,v^7 }{6\\,{\\sigma_4 }^3 }-\\frac{d^5 \\,v^6 }{\\sigma_2 }+\\frac{d^8 \\,v^7 }{8\\,{\\sigma_4 }^4 }+\\frac{d^7 \\,v^6 }{\\sigma_1 }-\\frac{d^7 \\,v^8 }{\\sigma_1 }+\\frac{d^9 \\,v^8 }{9\\,{\\sigma_4 }^{9\/2} }+\\frac{d^2 \\,v}{2\\,\\sigma_4 }-\\frac{d\\,v^2 }{\\sqrt{\\sigma_4 }}+1\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =7\\,{\\sigma_4 }^{7\/2} \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =5\\,{\\sigma_4 }^{5\/2} \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =3\\,{\\sigma_4 }^{3\/2} \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 ={R_o }^2 +d^2 \n\\end{array}"}}
%---
%[output:38ff51e6]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInf","value":"1-\\frac{v}{2}+\\frac{2\\,v^2 }{3}+\\frac{v^3 }{4}+\\frac{2\\,v^4 }{15}+\\frac{v^5 }{12}+\\frac{2\\,v^6 }{35}+\\frac{v^7 }{24}+\\frac{2\\,v^8 }{63}"}}
%---
%[output:079a3f68]
%   data: {"dataType":"symbolic","outputData":{"name":"relEq","value":"\\sqrt{1-v^2 }"}}
%---
%[output:50ce538e]
%   data: {"dataType":"symbolic","outputData":{"name":"relEqApx","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}-\\frac{v^6 }{16}-\\frac{5\\,v^8 }{128}"}}
%---
%[output:0d212725]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQGMHdV5709cU1hCiFlIG7xL3tJgkyjqM5CXt1sHCSKlMq+qXamErtdPKvHzi6wXQjcqrHcXt9hOY3vXi6WsIJEsZLmWqrVdCmqg6hPw+oCSbNkSN1gVosTQbM2ukxfAdkrJmkKdp282Zzk7vnNn5t65M3O++Y10Bd49c+Y7v+\/MPf\/9znfO+cDPf\/7znxsuCEAAAhCAAAQg4BGBDyBgPPIWpkIAAhCAAAQgEBBAwNARIAABCEAAAhDwjgACxjuXYTAEIAABCEAAAggY+gAEIAABCEAAAt4RQMB45zIMhgAEIAABCEAAAUMfgAAEIAABCEDAOwIIGO9chsEQgAAEIAABCCBg6AMQgAAEIAABCHhHAAHjncswGAIQgAAEIAABBAx9AAIQgAAEIAAB7wggYLxzGQZDAAIQgAAEIICAoQ9AAAIQgAAEIOAdAQSMdy7DYAhAAAIQgAAEEDD0AQhAAAIQgAAEvCOAgPHOZRgMAQhAAAIQgAAChj4AAQhAAAIQgIB3BBAw3rkMgyEAAQhAAAIQQMDQByAAAQhAAAIQ8I4AAsY7l2EwBCAAAQhAAAIIGPoABCAAAQhAAALeEUDAeOcyDIYABCAAAQhAAAFDH4AABCAAAQhAwDsCCBjvXIbBEIAABCAAAQggYBL0gZmZGSMfLghAAAIQgIBvBN674hPmxmuW+WZ2rL0ImBhEIlwGBgbM1NRULEwKQAACEIAABMpE4L0rrjX\/duMW8+iXr1cnYhAwMT3tueeeMxs2bDBjY2Omo6OjTP3yPFtEZI2Pj3thqxiPva3tTvCFryVAX6hmX3jt9Fnzv779uln6xsvmL+\/+vOnp6WktiJxrR8AkFDATExOld74VWz7YKtixt7VvO3zhawnQF6rXF06cOmvuOPSS+bt\/fMVc+sSg8WVcSOMpBAwCJk1\/ybQsX6qZ4jyvMvjCFwHT2j5QZr6jj\/\/QHHr+x2brp98xg3f8DwRMPl2hXE\/xaRCQfJ2HH37Y3Hrrraazs7NcIGtYg72tdRF84WsJ0Beq1RdEuIiAeWD9J83SN\/4pSIMgAtPaPlDK2n0SMKUEiFEQgAAEIJAbAZk6uu7rf2cG13SZwTVXezdVnwYUU0iKppDSOJ6yEIAABCCgi4CIl3Xf+r757MeXmW\/2fTJonOY\/whEwCBhdbzCtgQAEIFBBAjZpV1YevfBHv7FAAAFTwc5gm6zZ+RV2K02HAAQgoIqATdqV\/V4+1n4RAkaVdxtsDAKmQXDcBgEIQAACuRBwk3bDO+5qHsOYQmIKKZcXjIdAAALZE+CYk+yZ+lKjrDSVTzhpN2w\/AsYXj7bATs3ObwEuqoQABHIiwDEnOYEu6WO6u7vNncN\/Yr786OuLknYRMCV1WBFmIWCKoM4zIQCBOAI+HXMS1xZ+n46APRpCzjhavnLVoqRdBEw6lqpLI2BUu5fGQcBbAnw3eeu6pg23vm9fd6955N71i5J2ETBN49VTAV8SenxJSyCgiQDfTZq8ma4taXyfpmw6K4ovTRJvjA80O7\/47ocFEIBAowR8\/246deqU2bRpkzl27NgiBGvXrjUjIyOmra2tUTQN3Wd52pt3795tent7G6rL3jQ3N2eGhoZMX19fpocBp\/F9mrJNNbaAmxEwCJgCuh2PhAAEmiXg+8BkBczg4GCmg3sjXIXl3XffbQ4cOGBWrFhhrPBYvXp10yKmEXvi7knj+zRl455btt+rFDCjo6Omq6urbsdz1faqVavM\/v37TXt7+3n+0ez8snVG7IEABJIT8P27KU7AWBHx2GOPBVDsYYTHjx83W7ZsWQAl391y2WiO+30uZTdu3GhOnjxpor7no8SK3Ltr1y6zd+\/eYGw4cuSIGR4eDp61fPnyQOzIMmaJsJw5c8Y8++yzwc8eeeSR8\/5tIzBR9ogvZdx6\/fXXzac\/\/Wlz7733mq997WvGtj0cDUrj+zRlk\/e+cpRUJ2CkE+zbt8\/UC\/\/ZTnTfffcFyl865uTkZM2wpWbnl6MLYgUEINAIAd+\/m+IEjHyXi\/CQ6SSZZrIREmElosR+f4enaeT7fHp62nzpS18KRI2N8Lj1udNTVhDt2bMniL7UuqzAsH\/oSl1y\/cEf\/EEgYETQyHOsLeF\/i4BZuXJlpD1u+8QG2wapUzjt2LHDbNu2beGP7DS+T1O2kX5Y5D1qBIx9GZYtWxbwvOWWWyIjMG7niIOv2flxbef3EIBAeQmEv5tkQzP5lPWS7e3dLe5r5cDYKIkIDDd3xI2S3HDDDUEExgqOcKTEtj8sOqLKRf28Hkc7hlgBY6eawtEcV1xJfSJ8rAhyn\/uDH\/xg0e\/ixqg041KasmXtO1F2qRIwp0+fDpSwdPyoucu0c5uane9bZ8VeCEDgfQLh7yY5C2f08enSIhpc02UG11y9YF+9CEyt39nUABEw7tROWKi4AmbDhg2LeNipHzfSkiQCE57Okko3b968EIGxU0ThaFBYwETZ8+abb5pDhw4tmgWwswnyLDt9Fm5b+OdR0SN5bpKype08EYapETC2fXECxf5+zZo15sEHHwxCk0lyYPr7+43sfCiX3cLZN2djLwQgoIeAlghMrSTeKCEgf5iGBUy9CExYFNTyftSY4U7dPPnkk4vSDMIRmKQCJsoe8WXU78SOu+66y9xzzz0LU1xxf1i7R0zMzs6agYEBBIwPr35SAXPixImFUF7U3Ki0N7y0Tn4mYkY+XBCAAASKIhA3iBVlV9LnNpMD40ZgwvXYnEaZZvrKV76ykAMjPz98+HDNBRtRq5BsLoubJyljjOTWSP6knUJKImDCOTCuPTKF5AoYdyFKIzkw4+PjRj7uRQQmac8ssFxSAeNOMdULIdovibGxMdPR0UEEpkDf8mgIQOB9AtoFTL1VSK6AESJRq3vcn9eaPnL7U\/iPVZkikuiQXG6+jkTsJUH48ccfX1gtlETAiOCJsiccgQnnB6WdQnIjMPbYAQSMB98ecQJGmhBeZl0vicv3LwkPXIaJEIBAAwT4bmoAmpJb0vg+TVnf8FQuB0YcFLUkzqrtWqpco3r1rbNiLwQgoCcCgy8bJ5BGlKQp27hFxdxZCQFTaztnN1xYb+tqzc4vpsvxVAhAIAsCfDdlQdHPOtL4Pk1Z32ioEzBZO0Cz87NmRX0QgEB+BPhuyo912Z6UxvdpypatnXH2IGBiCGl2flzn4PcQgEB5CfDdVF7ftNqyNL5PU7bVdmddPwIGAZN1n6I+CEAgBwK+D0xxy6hbiTDt7rs2DeHo0aMLBz5a++wZSVmcXp20zWl8n6Zs0ueXpRwCBgFTlr6IHRCAQAoCvg9MPgoYObTxt37rtxaOqRFhs3PnTvPiiy+a9evX53ZydRrfpymbovuVoigCBgFTio6IERCAQDoCvg9McQKm1kILIeSe\/iyrQ2VfFvlZrVOrZb+YSy+9dNHv7IZydhd2WcTx9a9\/fRH8cDTFRmCuu+468+qrr5qtW7caOa9JIjkHDx4M7v31X\/\/1QMCE93Bx95NxjwdI8vOoHpHG92nKpuuBxZdGwCBgiu+FWAABCKQm4PvAVE\/A2A3f5MRpK1BkI7rw6c8CzZ4MLdtgCJPwqdV33nlnICzcHddlo7fwZnhSV61t++Xn7hE0Dz300MK2\/jJ9JJecft3V1WXWrVu36Cw+d5NU97wjuUciN7fffruJ+nnUydhybxrfpymbuhMWfAMCBgFTcBfk8RCAQCMEwgPTe69Pm3d\/Ut7DHC\/4lS6z9CNdC02tJ2Ckbe7W+nbvrgceeCA4hdrupB6uw90y4\/LLL190arVbZy0BU88et95nnnnG3HTTTYGwsiLkkUceCQSMCCX3cgVRrQMbrRhJcmaTW28aUZKmbCP9sMh7EDAImCL7H8+GAAQaJBAemE7\/+Q5z+qHtDdbW+tsuu227uez3tiUSMO7ZQ3aqRiImMtUjAsZu3R+errGVyxRQrVOrrVAIC5i4HdzDp0qLiPnd3\/3dhSiOHAxsBYw79SX2uEcY2IRf+bk7TRX18yivpBElacq2vhdk+wQEDAIm2x5FbRCAQC4EqhyBcQVM+KRmCz+80igqAiMCSXJo3PPxwg50BYzk0OzYsSM4Gfrtt98Ozkuyx9P85m\/+ZnDQoz1hO25KyrbDPq\/Wpqu1OlMaUZKmbC4dN8OHIGAQMBl2J6qCAATyIuD7wNRMDow78Ls5MG7ujEwhuXkuUQJGoidy1TpKJkpY2GRce8RMlICRyMr9998fLL3+h3\/4hyBXRp5jVy9JDkzUz8mBiX+TEDAImPheQgkIQKB0BLQIGFkN5F5WFNRbheQKmPCp1XZqpl4ERu6RSMm5c+fMlVdeaZ544olFNrgrhOQX4ciICJPDhw+b\/fv3m\/b29kUHBLvTQX\/8x39sXnjhhWDKK7xaytoZZX+9DpfG92nKlq6TxxiEgEHA+NZnsRcCEEi5EgVgugikESVpyvpGCQGDgPGtz2IvBCCAgKl0H0gjStKU9Q0qAgYB41ufxV4IQAABU+k+kEaUpCnrG1QEDALGtz6LvRCAAAKm0n0gjShJU9Y3qAgYBIxvfRZ7IQABBEyl+0AaUZKmrG9QETAIGN\/6LPZCAAIImEr3gTSiJE1Z36AiYBAwvvVZ7IUABBQImLjDHFvp5PAS6zTPEkFgz1uqt1dLmjrTlk0jStKUTWtH0eURMAiYovsgz4cABBog4PvA5KOAsRvQibvs6dMNuK7pW9L4Pk3Zpg3LuQIEDAIm5y7H4yAAgSwI+D4wxQmYehvZnTlzxjz77LNGNr0LbxBnN8KzUZZLL73UPPbYYwFy+Z0cBSCb2MkGenLv2rVrgzOW3Ms9p8j9udR58OBBc9ttt5lvfOMbZu\/evcFGdnYzuuuuuy7Y3O7kyZMLZx3ZXXrloMekRwXE9Y80vk9TNu65Zfs9AgYBU7Y+iT0QgEACAr4PTM0cJSAHJNqt\/92jBNzpHUG4ceNGc+eddwanREs5ERYjIyOm1mnUUj7q7CLrDtllVy5bn5xK3dPTsyBMTpw4EQgYOXl6y5YtwcGTnZ2dC6dWy8\/lIMh6xxYkcL1J4\/s0ZZM8u0xlEDAImDL1R2yBAAQSEggPTOd+NmP+42czCe\/Ov9gvXdxpllzcufDgegLGPbdIDluUf4sAeeCBBwJRYA9eDNfhRjjkLCQrIiRXxa2zloCJiwi55xfZ+qwYqXWatRt5kciNRHmWLFmyELVpxgNpREmass3YVMS9CBgETBH9jmdCAAJNEggPTHMvf8PMvTzeZK2tu73t2n7Tdu1XEwkYiXRMTk4G0RIRMHY6SESACBj3NGo7HeRaLlNAN9xwQ6LDHN0poHonUtuDIiWKYy+ZgpKIiz3ROnzIZFdXVxCtkcsVNM1STiNK0pRt1q6870fAIGDy7nM8DwIQyIBAlSMwroC56667zD333GPCK4LqHeboRmCs+KgnXqwAkf+60z9WlKxbt84MDQ0tRIbCERnx1V\/91V8FXpcTqJtdvZRGlKQpm0G3zLUKBAwCJtcOx8MgAIFsCPg+MDWTAxOOdFhhYaMk9913n5EppF27di1M2URNIT344IOBQ+rlpUTZaiNF9957r\/na174W1GNzbOz0ldixY8cOs23btiA3RpKAt27dGkRtGr3S+D5N2UbtKeo+lQImTahOOrw7Txp2hGbnF9XpeC4EINA8Ad+\/m6wokNVA7mVXEdVbheQKGBvtsCuN7AqiehEYuUemns6dO2euvPJK88QTTyyyYfPmzYsEjc3BkekimXKylxVMO3fuNI888oiR5OJ9+\/YFv3ZXSLn2iuiZnp5uKpE3je\/TlG2+V+ZbgzoBI+JFOlDUMjgXr+34R48eNQcOHKgZ1tPs\/Hy7Gk+DAASyJMB3U5Y0m6srq+XRSa1I4\/s0ZZM+vyzl1AgYq+aXLVsWsL3lllsWkqeiYFtVLb+XxLBa85KanV+WTogdEIBAegJ8N6Vn1qo7EDCtIlu\/XlUC5vTp00EIz02mimq+CB6Zl5TQnkRt4gRMf3+\/6e7uDqqTdf3y4YIABCBQFAEETFHki39unO8lSVk+cs3OzpqBgYFgSkv2rNF0qREw1im11uPXcpjdkEiW2iXJgXHrEDEjHy4IQAACRRGwg5j7x1VRtvDcfAnEiZLx8XEjH\/dCwOTro4aelkTA2O2gJRNcVGoSATM2NmY6OjqIwDTkFW6CAASyJiDfXfKX9dTUVNZVU58HBGRGQMalWrMBbgRG+oeIGQSMB05NImBkyshuAc0qJA+ciokQgEBNAu5ABaLGCbx2+qy549A\/mY+1X2QeWP+JxivK8c6kqQxx0005mpz5oyo3hRS1dE\/I1lKomp2feW+iQghAAAKeETj0\/I\/NHYdeMn2f+aj5Zt8nPbM+3lzNY1jlBEzY3URg4l8ASkAAAhDQSGD08R+a0cenzeCaLjO45mqNTUx18KNvACohYOotcUPA+NZlsRcCEIBAcwROnDprDj3\/IyPRF4m8aBUvQokITHN9xeu7NTvfa8dgPAQgAIEGCIh4Wfet7wd3inARAaP50jyGqYvAZN0RNTs\/a1bUBwEIQKDMBFzx8uiXrw+SdrVfmscwBExM79XsfO0vLu2DAAQgYAlUId+llrc1j2EIGAQM33AQgAAE1BKoUr4LAkZtN26sYZrVa2NEuAsCEICAHwSqlu+CgPGjX+ZmJQImN9Q8CAIQgEBmBKqY74KAyaz76KgIAaPDj7QCAhCoDoGq5rsgYKrTxxO1FAGTCBOFIAABCBROoOr5LgiYwrtguQxAwJTLH1gDAQhAoBYBpoxq9wvNYxirkGK+CzQ7n69BCEAAAhoI2Cmjz358WXCeURX2d0nqN81jGAIGAZP0PaAcBCAAgVIRYMoo3h0ImHhGaktodr5ap9EwCEBAPQGmjJK5WPMYRgSGCEyyt4BSEIAABEpAwEZdtJ8inRVqBExWJD2sR7PzPXQHJkMAAhUmIOLljkMvmddOn1V\/inRWbtY8hhGBIQKT1XtCPRCAAARaRuDQ8z8OxIsk6FblIMYsYCJgsqDoaR2ane+pSzAbAhCoEAGmjJpztuYxjAgMEZjm3g7uhgAEINAiAt955Yz5yuGXgtr7PvNRM7jm6hY9SW+1CBi9vo1tmWbnxzaeAhCAAAQKIEDUJTvomscwIjBEYLJ7U6gJAhCAQJME3OXRRF2ahGmMQcA0z9DbGjQ731unYDgEIKCOgBt1YUfd7NyreQwjAkMEJrs3hZogAAEINECAqEsD0BLegoBJCEpjMc3O1+gv2gQBCPhDgKhL632leQwjAkMEpvVvEE+AAAQgECIg+7rIIYxykevSuu6BgGkd29LXrNn5pYePgRCAgDoCdjfd7756xgyu6TJ9n7mS06Nb6GXNY5jKCMzo6Kjp6uoyvb29NbvFqVOnzKZNm8yxY8eC369du9aMjIyYtra288prdn4L3xmqhgAEILCIgDtdJLvpPrD+k+bGa5ZBqcUENI9h6gSMiJd9+\/aZ3bt31xQwc3NzZmhoyKxevTr4vf338uXLzeDgIAKmxS8T1UMAAtUjQJJucT5HwBTHPvGTbVRl2bJ5RX\/LLbdERmDClR45csRMTk7WjMJodn5iuBSEAAQg0ACBcNSFM4wagNjkLZrHMDURGBEwp0+fNhJJcSMsSXyPgElCiTIQgAAEkhMgSTc5q1aWRMC0km7GdYeniOKqt5Gb9evX14zYWOf39\/eb7u7uoLrOzs7gwwUBCEAAAosJuEm6dnWR5Lxw5UdgZmbGyEeu2dlZMzAwYCYmJkxPT09+RuTwJDURGMsqjYCxZeXeuCRe1xciZuTDBQEIQAAC8wRI0i1PTxgfHzfycS8ETHn8E2lJUgGTRLzIQ2wEZmxszHR0dATPJQLjQUfARAhAIBcCYeHCni65YK\/7EDcCMzU1FYgZBEzxfom1IImAiVt55D5E8\/xhLEwKQAACEKhDgDyX8ncPzWNYJaeQZKn1yZMnI6eNEDDlfymxEAIQKI4AeS7FsU\/7ZARMWmIFlq8VgbE\/6+vrMytXrly0iZ01ddWqVWb\/\/v2mvb19kfWanV+gm3g0BCDgIQHyXPxzmuYxTF0EJuvupdn5WbOiPghAQCcB8lz89avmMQwBE9MvNTvf31cSyyEAgTwIIFzyoNzaZ2gewxAwCJjWvj3UDgEIeEcA4eKdyyINRsDo8WXqlmh2fmoY3AABCKgmYIWLrC6SS5ZEc1q03y7XPIYRgSEC4\/fbifUQgEDTBES4fPfVM2b08R8GG9Kxg27TSEtTAQKmNK7I3xDNzs+fJk+EAATKRKCWcJGIy43XzB+Ky+U\/Ac1jGBEYIjD+v6G0AAIQSEXAnSqyEReESyqE3hRGwHjjquwN1ez87GlRIwQgUGYC4RyXz358WZDjQsSlzF5rzjbNYxgRGCIwzb0d3A0BCJSeQHhVkQiXwTVXG06JLr3rmjYQAdM0Qn8r0Ox8f72C5RCAQBICCJcklHSX0TyGEYEhAqP77aV1EKggAREusqJIlkNLlIXl0BXsBL9oMgKmur43mp1fYbfSdAioI2BXFE38\/Y+CJdFWuMhUEVd1CWgew4jAEIGp7ptNyyGggEB4RZEIFxEtEnXhggACpsJ9QLPzK+xWmg4B7wnUym9hRZH3bs28AZrHMCIwRGAyf2GoEAIQaA2BqGkitvtvDW8NtSJgNHixwTZodn6DSLgNAhDImUCtaSK73X\/OpvA4zwhoHsOIwBCB8ex1xFwIVIeArCL6ziunF1YTsfFcdXyfVUsRMFmR9LAezc730B2YDAH1BIi2qHdxrg3UPIYRgSECk+vLxMMgAIHzCdTKbSHaQk\/JggACJguKntah2fmeugSzIaCGwHdeOWMOPf+jYIpILvZuUePa0jRE8xhGBIYITGleNAyBQBUIRE0Rffbjl3GoYhU6QM5tRMDkDLxMj9Ps\/DJxxhYIaCbAFJFm75a7bZrHMCIwRGDK\/fZhHQQ8JVBLtFx12UVmw3+9kl1yPfWpj2YjYHz0WkY2a3Z+RoioBgIQ+AUBRAtdoWwENI9hRGCIwJTtfcMeCHhFICxaxHjOI\/LKhaqNRcAodu\/o6Kjp6uoyvb29NVup2fmK3UrTINAyAiJY5JLVQ7KKyJ78LNNDN16zzLCtf8vQU3EDBDSPYZWOwIh42bdvn9m9ezcCpoEXg1sgUBUCNspid8W1URbZq+XGay4jp6UqHcHDdiJgPHRaPZNPnTplNm3aZJYtWxYUu+WWWxAwynxMcyDQLIGoqSE2mGuWLPfnSQABkxPtI0eOmOHh4UVPqxcdadQsETCnT582y5cvN0NDQ2b16tWxAqa\/v990d3cHj+zs7Aw+XBCAgB4CbpTF\/r+NtMjBiezTosfX2lsyMzNj5CPX7OysGRgYMBMTE6anp0dV00szhSTi5fDhw2b\/\/v2mvb09gGwjJevXr48UGM14Y25uLpGAeWuyb9FjRLx0dCQTML90cXS5JRd3RJq\/JHTfkrbF9bj1hss2w4R7IVAVAiJS5oXK6YVcFitYmBqqSi\/Q2c7x8XEjH\/dCwLTI11aoDA4OnqcQJfwluSqusMnKjKQC5qX\/e6d574prIx8rKw5u\/PhlC7\/\/j5\/NK99617m588ucS3BfXL31hI8VPa5wsuVdgfR+uWQiLc4mfg+BMhCwguW102fNxN\/\/KEi+tYLF7s9iE3HLYC82QKBRAm4EZmpqKhAzCJhGacbcV3YBs2HDhsD5y1deV\/MvNts8ETLyl9v8f7PZFryWqAkLJCuGwmXP\/Wy2pqiKKh\/lJlcUidCJEkJWBMnviQq16GWh2sQE6kVYrFDJ6j1NbBQFIZAzAXJgcgBe5ikkK2BqzR+6X5KCyS6rdEWN\/L\/MoctVxi9MV\/hYcVRL5FhB5JaJixrVEj\/hCJAVPhdcoWt+NofXhkc4BGzeyolTc+dNCSFY6CpVJYCAycnzeSXx2uYknUKqJ2Ci0ER9mbrCxg1Xl1HYpHF7WASFBZCIn6TCx53Wmo\/mvJ8rtPTyeZFDlCeNd\/SVdaeDZGmzTAfZ\/VkkAopg0edzWtQYAQRMY9xKf1crBUytxrsbYNWK1sjP3C\/fskZssnSsFT4ibkT0uELIip7wz93nR4kd+bmd7mI6K0uP5V+XfW9EpIhYCa8Qmn9P5vdjIYclf\/\/wxHITQMCU2z8ttS4P5ycVNvaLWkTOx9rbKvllLQLHRnLee\/O5Bd+\/+8ZU8P9RYidK6EhEh2hOS1+hVJWHIyuuWLECH7GSCimFK04gjzGsKMSFLaN2E3dXrlwZbCx37NixmhxWrVrVklVISaAX6Xz3L0+Z138\/32Z+9YS9akVt5HeyrXlVLzeyY4WOO43lih\/LyBU5F1wxv+ePjeSQn5NtT3L7ttQcjqwgVrLlTW3VJVDkGNZq6oUJmFY3LKv6y+r8WsnDsjzUDgxhcTMftalu5KZWf6g1fRU3bbV4aqoDgRPzoiUVKlINkZWsvrWoBwLvEyjrGJaFj0ohYIpaRp0EoI\/OryVuRLwcev7H5zXZRm\/CAsf+OwkjzWWsyHn3zecW8nPsdFVUFEdyb2wExyYda47gWJEi\/xURbaMp0i\/sXis2ohIWKvQzzW8PbSsDAR\/HsKTcEDAxpLQ5PzwtNR++PxMMPOHoTXjQkcFm\/mfz+TcMPuZ9URMSOLVycXyO3ri5KfN9Zj6ZNi7qZxNr6StJv5IpB4FsCWgbw1w6hQqYWsuma7lu8+bNRnbpLeLS7PxaPN2\/pmWLdbncwaueyAmLGoTOvMCRpGMbrYmK3ryfazOfe5NXcnE4ejLv7\/fzraIEipSz\/rYiRX5W5byrIr6feCYE4ghoHsMKFTAWfL0ppDjntPr3mp3fKLt6IkfqdKcNws+wUZxaERwreGzkx5Zt1M6y31dL3NSL3Mi0VJKkYtc\/wkBEiI2aWEFaz0+Wu+SkzIuS+WMyiLqVvUdhHwTOJ6B5DCuFgClzp9Ps\/FZzD\/\/OQTH2AAAfIElEQVR1L3\/Z2ysuquPa5goZGUStwHH\/O\/\/\/81NbUfeWXRC50S0RN9Ozr5h333guECBXnnvBdC97eVHbrJiZOnNtUGbqzCfMQ\/\/8nyLdavOdXIEozBAnrX4TqB8CxRHQPIYhYGL6lWbnF\/dK1X5yrchBWPTInbacjSzUmtaKa1uUmAkLoHr11Kqjni3WXrfOtLb\/xvJ\/M50XvWm6l\/1TINZE1ISTiWfOXmFmz14eTEOJjb+24majOYk4ztf8HgJVJqB5DCuNgDl+\/LjZuHGjOXny5Hl9rar7wFT5pWu07a4gcP\/fFQ+uKIoTE2kFRjgqVKsdYeFjoyC2rCui3k+cXhxZCtdrV0vNvXTYLLlkqZFcm7CwcfNskkxFNeoD7oMABMpDAAHTYl+4W\/qvW7fODA0Nmb6+PmM3uJME3loHKbbYrKB6zc7Pgx\/PKJaAFTbvvPYXgSFRwubCq24Nfi9RG6I1xfqMp0MgSwKax7BSRGDCSbyjo6Omq6vL9Pb2BgLi0KFDZmRkxLS1zc\/X53lpdn6eHHlWuQi4CcR28z43YuNGaxA15fId1kAgDQHNY1gpBYwsr56eng6WTgt8ETT79+837e3tafyWSVnNzs8EEJWoIiDCJipaE56CuvCqL6hqO42BgEYCmsewUggY6TQiUuQKi5Ynn3zSTE5OEoHR+GbRJi8IiKixOxGHp6CI1HjhQoysMAEETA7Od\/NgZOpIBM2+ffvM8uXLzYEDB8yKFStysOL8R2h2fiFAeagKAogaFW6kERUgoHkMK0UExooXSdwtKlk3qh9rdn4F3l2amCOBOFFzweU9ZsnFHSQK5+gTHgUBzWNYKQQMO\/HykkFAJ4G4nBpZ\/STTUOTT6PQ\/rSqeAAImBx9I4m6RuS5EYHJwMo+oPAF39VOtfBqiNJXvIgDImAACJmOg4epsBObYsWM1n8RGdi12ANVDoEACImpqbcBHgnCBTuHRagggYNS4Mn1DNDs\/PQ3ugEDrCcRFaWTaib1pWu8HnqCDgOYxrBQ5MGXuJpqdX2bu2AYBl4DNpak17WR3EW679qtAgwAEQgQ0j2EImJjurtn5vOkQ8JWAmxw89\/L4QjPmE4Lnj0VA0PjqXezOkoDmMQwBg4DJ8l2hLggUQiBO0DDlVIhbeGgJCCBgSuCEokzQ7PyimPJcCLSaAIKm1YSp3xcCmscwVREYWYo9PDwc9Kvdu3cHh0FGXW7ZtWvXRh5VoNn5vryA2AmBZglE5dCYcx81bZ\/sJSm4WcDcX1oCmscwNQLm+PHjZsuWLWbPnj1BR7L\/X+sIAveASDnhemhoKDiyQM5hCl+anV\/aNw7DINBiArUEjeTP2H1oyJ9psQOoPjcCmscwNQImvBGenKXU1dVVMwoTLltvEz3Nzs\/tDeJBECgxAXfZ9ttH95ollywNrLUJweTPlNh5mBZLQPMYpkbAuKdZi0fD\/3a9XCsCs3r16ppiR7PzY3s+BSBQQQI2OjP34jPGLH1hQcwsaes0F1zRzeqmCvYJn5useQxTJWDciItEVaanp2tOC0lnlCmnjRs3mpMnT5qJiYnIQySt8\/v7+013d3fQjzs7O4MPFwQgoJ\/AgqB56YgxS35MdEa\/y71v4czMjJGPXLOzs2ZgYKDuOOdrgyspYETcHD582Ozfv9+0t7fHRms2bNiwyL8iZuTDBQEIVIuAFTPv\/Mt3zbmz31sQM+TOVKsflL214+PjRj7uVe8P9bK3J8o+VQJGGmkTcaOmkObm5oKkXXfKyE0ADif92gjM2NiY6ejoIALja0\/Hbgi0gMA7\/\/KdQMjUyp2RE7Ylj4YLAnkTcCMwU1NTgZhBwOTthRTPC08ZRSXxNipgNDo\/BV6KQgACMQSiVjZxdhNdp0gC5MAUST\/hs9Mso641hSS5MCMjI0aWVbuXZucnREsxCEAgJYGoRGDETEqQFG+agOYxTM0Ukng5aiM7G3Xp6+tbSNaVCM2+ffuCzsFGdk2\/I1QAAQhEEHjv9Wnz9vN\/Zn7545cYe26TXaI9\/98vwA4CLSOAgGkZ2vJXrNn55aePhRDQR8BGZ2zejIgYlmjr83NZWqR5DFMVgWlFh9Hs\/Fbwok4IQCA5gbCYkTvnN85jv5nkFClZj4DmMQwBE9P3NTuf1x4CECgPAcRMeXyhyRLNYxgCBgGj6V2lLRBQQaDWNJMkAJMzo8K9uTYCAZMr7nI9TLPzy0UaayAAgVoEosQMZzTRX5IQ0DyGEYEhApPkHaAMBCBQAgKymuncO0fNu6f+MrDGrmZi07wSOKekJiBgSuqYPMzS7Pw8+PEMCEAgewJ2afZ\/nD266MBJmWZCzGTP2+caNY9hRGCIwPj8bmI7BCpPQMTM3EuHzTsnHjZLLn1jUWSm7dqvVp5P1QEgYCrcAzQ7v8JupekQUEnA5su89bejZumvXhS0se3afpJ\/VXo7WaM0j2FEYIjAJHsLKAUBCHhFwObLvHPiL8ySS5aSL+OV97IzFgGTHUvvatLsfO+cgcEQgEBqAjZf5p3XHjZLrzzDFFNqgn7foHkMIwJDBMbvtxPrIQCBxASilmST+JsYoXcFETDeuSw7gzU7PztK1AQBCPhEgKiMT95qzlbNYxgRGCIwzb0d3A0BCHhNwK5iksMl5cRsuSTxl6iM125dMB4Bo8OPDbVCs\/MbAsJNEICAWgI\/\/euvG3JldLlX8xhGBIYIjK63ldZAAAJNE5CojKxekuXYEpVhx9+mkRZWAQKmMPTFP1iz84uniwUQgECZCdhcGZleumjVssBUK2bYJK\/MnnvfNs1jGBEYIjB+vIVYCQEIFErgraf\/1My9PG6WfOiNRfvKIGQKdUvswxEwsYj0FtDsfL1eo2UQgECrCNSKypD02yrazdereQwjAkMEpvk3hBogAIHKEYiaXrJipnJAStpgBExJHZOHWZqdnwc\/ngEBCOgncPrPd5h3T\/8l00sldLXmMYwIDBGYEr5ymAQBCPhIYO7Fp82bB78UCBmb9Mv0UrGeRMAUy7\/Qp2t2fqFgeTgEIKCWgEwv\/fSvdwZLsREyxbpZ8xhGBIYITLFvF0+HAATUEiBPpnjXImCK90FhFmh2fmFQeTAEIFApAiJk3nrqoJl7+RvzG+NdspT9ZHLqAZrHMFURmCNHjpjh4eGgW+zevdv09vZGdhHrVCmwatUqs3\/\/ftPe3n5eec3Oz+n94TEQgAAEAgJWyLz1tyPB1NLSX70IIdPivqF5DFMjYI4fP262bNli9uzZE3QH+\/8rVqw4r3tI2Y0bN5r77rvP9PT0GBE+k5OTZmRkxLS1tS0qr9n5LX5vqB4CEIBATQJWyPz0f3\/dXPSfly06roCN8bLtNJrHMDUCJixCRkdHTVdXV80ojJSdnp42g4ODsT1Fs\/NjG08BCEAAAi0k4AqZX\/61S4KoDEcVZAtc8ximRsCIYJHLipLwv22XmJubM0NDQ2b16tV1p5hseev8\/v5+093dHfy4s7Mz+HBBAAIQgEDzBBAyzTN0a5iZmTHykWt2dtYMDAyYiYmJYMZB06VKwLgRl6goixUwa9asMQ8++KA5duxYohwY1+kiZuTDBQEIQAAC2RFAyGTDcnx83MjHvRAw2bBtSS3hKaM4AXPixImFxF259+TJk3VzYMbGxkxHRwcRmJZ4j0ohAAEIvE+glpCR38qmeOTIxPcUNwIzNTUViBkETDy3wko0M4XkJgCHk341zx8W5iweDAEIQCABASIyCSDFFNE8hqmZQgpHXOol8YZ\/JwJm165dZu\/evectpdbs\/OZfDWqAAAQg0HoCImTkvKW3v\/dnhmTfdLw1j2FqBEyaZdTiUBExdu+XqIRf6SaanZ\/uNaA0BCAAgWIJiJD5yTc3mn\/\/l+8gZBK6QvMYpkbAiC+jNrKzibt9fX0LWdjuRnZr166tmf+CgEn4hlAMAhCAQI4EJBpz+qHtwY6+RGTqg0fA5Ngxy\/Yozc4vG2vsgQAEIJCGgCtk3A3x7AnYaerSWlbzGKYqAtOKDqjZ+a3gRZ0QgAAE8iRgE31tRObi1VcsHFHwwevGzAVX6Nr7JC1bzWMYAiamN2h2ftoXgfIQgAAEykrA5secffFp88tXX2Mu\/Z3\/Yt578zmz9PIec8n1Y8EOv1W8NI9hCBgETBXfadoMAQgoJTD34tPmR9s\/Z5Z+pCs4Y0mOJzj3s5nK7iGDgFHa0ZM0S7Pzk7SfMhCAAAR8I+BOK4mQueTz1xmz9IWgGVXbDE\/zGEYEhgiMb99N2AsBCEAgEQF3WunDv\/0\/zcU9nzJzL49X6sBIBEyirqKzkGbn6\/QYrYIABCCwmIBdrSTRmA\/\/9qZgaskKGe2JvprHMCIwRGD4roMABCCgnoDdzfetp\/\/UXPSpm83lG\/\/EvDPzzSDR98KrvmA+eP2YSgYIGJVuTdYozc5PRoBSEIAABPQQkCTf17+50Yiguey27ebi3\/iUefv7A2rzYzSPYURgiMDo+WaiJRCAAAQSErDTSjYac+7s91ROKyFgEnYIjcU0O1+jv2gTBCAAgaQEJApzctvnFqIxkh\/zb98fULV\/jOYxjAgMEZik7zrlIAABCKgk4EZjlm9\/yrz7xnPm7RcGVOwfg4BR2WWTNUqz85MRoBQEIAAB\/QRsboy09CN3HDBtn7rZzL38De+nlTSPYURgiMDo\/2aihRCAAAQSEghHY2QXX1ly\/c5rf+HlsQQImISO11hMs\/M1+os2QQACEGiWgHscgY3G+DqtpHkMIwJDBKbZd537IQABCKgj4O7i+6GbvxhMK8nl27QSAkZd10zeIM3OT06BkhCAAASqScCNxizf8VRwSKRMK9nVSmU\/W0nzGEYEhghMNb+VaDUEIACBhARsNOa9n0wvJPj6Eo1BwCR0ssZimp2v0V+0CQIQgECrCMgOvnIUgezge9nvbQseE07yvfSzh1r1+Ibq1TyGEYEhAtPQS8FNEIAABKpIQASMCBnZwVf2jLGXJPm+NdlXupOuETBV7KW\/aLNm51fYrTQdAhCAQMME7A6+UoFdpWSjMbLcWpZdL728x5QhGqN5DCMCQwSm4ZeYGyEAAQhUmcDJ7Z8zZ198etGUkvCwS67l\/y+86lbTdu1XC8OEgCkMffEP1uz84uliAQQgAAG\/CYQ3vnNbY5dcSzTmkuvHgumlvC\/NYxgRGCIweb9PPA8CEICAKgLuUuuPfeuHi9rmRmM+eN2YueCKnlzbjoDJFXfjDzty5IgZHh4OKti9e7fp7e2Nrez48eNmy5YtZs+ePWbFihXnldfs\/Fg4FIAABCAAgUQEopZa25vdaEyeuTGaxzA1ERhXiEiHqSdKFjrU3JwZGhoyR48eNQcOHEDAJHpNKQQBCEAAAlEEbF7MldufCg6EdC93pVJe0RgEjAd9VaIvk5OTZmRkxLS1tZnR0VHT1dVVNwojjpVychGB8cDJmAgBCEDAAwI2L8bdL8aanfcuvggYDzqMFSKDg4OBteF\/h5tw6tQps2PHDtPX1xeURcB44GRMhAAEIOAJAbtfjHuOkmu6LLd++\/sDQWJvK6MxCBgPOkw44iIRmenpaWMFTbgJ8nu5brjhhkQ5MP39\/aa7uzu4p7OzM\/hwQQACEIAABKII2OTe8KZ3rY7GzMzMGPnINTs7awYGBszExITp6ck3gbjVPUNNDkwaASP5MgcPHjRbt24NnJwkidd1hIgZ+XBBAAIQgAAE6hGwm94t\/ZWuRTv3uvdkneA7Pj5u5ONeCJgS99M0U0hS9qabbgrUaNJVSGNjY6ajo4MITIn7AKZBAAIQKCMBd4VSeJm1G43518m+4J\/NTim5EZipqalAzCBgytgzfmFTeMooKolXcl82bdpkjh07dl5rajlY8\/xhid2JaRCAAARUEUgiYqTBNhpz4VVfMB+8fqxpBprHMDVTSI0so5aekTQCo1G9Nv1mUAEEIAABCKQiIMus3\/vJ9KIzlMIVZLncGgGTyj3FFY7ayG7uF\/u9yIqjcBITAqY4f\/FkCEAAAlUkkETEuMutP7T6UMM7+CJgqtjDftFmzc6vsFtpOgQgAIFCCSQRMVlMKWkew9RMIbWqJ2p2fquYUS8EIAABCMQTqLdrr3u3O6V06epDqQ6F1DyGIWBi+phm58e\/XpSAAAQgAIFWEkgaiWl0SknzGIaAQcC08t2kbghAAAIQiCFgRUzUEmv39rSrlBAwFe5+mp1fYbfSdAhAAAKlIpBGxKSZUtI8hhGBIQJTqpcYYyAAAQhUlYCIGLmWb38qFoGdUjo3N1N34zsETCxKvQU0O1+v12gZBCAAAf8I2M3ukooYKScHQsrBkG3X9pu2a796XqM1j2FEYIjA+PeWYzEEIAABpQTs2Ultn7o52OwuyVUvLwYBk4Sg0jKana\/UZTQLAhCAgNcERMSc+PLV5kM3fzGxiHHzYpZ9\/tmF9msew4jAEIHx+kXHeAhAAAIaCcy9+LT50fbPBQJGhEySq1ZeDAImCTmlZTQ7X6nLaBYEIAABFQRO\/\/kOc\/qh7ebK7U8ZmVJKev3rd\/vMe28+FxwG+f3ZTrNhwwZOo04KT1M5BIwmb9IWCEAAAn4RSLpbb7hVNi\/m9Q\/0mP\/WP4WA8cvt2ViLgMmGI7VAAAIQgEBjBNLsEeM+webFHH3FGDkQMnyYcWPWlOcucmBifIGAKU9nxRIIQAACVSRgVyYt\/ZWuRHvEuIwmn3rYHP+bu83nfv+Q+djKHlX4EDAIGFUdmsZAAAIQ0EjAJvVedtt2c9nvbUvcRM1\/hCNgEDCJXwQKQgACEIBAcQQaSepFwBTnr8KfrNn5hcPFAAhAAAIQSEUgbT6M5jGMCAwRmFQvD4UhAAEIQKA4AnaTu4s+dXOifBgETHG+KvzJmp1fOFwMgAAEIACB1ARsPkyS\/WE0j2FEYIjApH55uAECEIAABIolkDQfBgFTrJ8Kfbpm5xcKlodDAAIQgEBTBCQfRq7l25+KrEfzGEYEhghMUy8QN0MAAhCAQDEEkpyXhIApxjeleKpm55cCMEZAAAIQgEDDBOxU0se+9UOz9CNd59WjeQwjAkMEpuEXhxshAAEIQKB4Aie+fLWJ2qUXAVO8fwqzQLPzC4PKgyEAAQhAIDMC9VYlaR7DVEVgjhw5YoaHh4NOsXv3btPb21uzg5w6dcps2rTJHDt2LPj92rVrzcjIiGlra6tU+C2zt4eKIAABCECgUAJRG9whYAp1S7KHHz9+3GzZssXs2bMnuMH+\/4oVKxZVMDc3Z4aGhszq1asDgWP\/vXz5cjM4OIiASYabUhCAAAQgUCICdoO78FlJCJgSOSnKFIm+TE5OLkRSRkdHTVdXV2QUxq0nfK\/7O83O98CtmAgBCEAAAgkJvP7Njeatp\/\/UuAm9mscwNVNIIljkslGU8L\/r+T+JgOnv7zfd3d1BNZ2dncGHCwIQgAAEIFAmAv982wfMT6\/5vHnnd\/4kMGt2dtYMDAyYiYkJ09PTUyZTm7ZFlYBxIy4iSqanp2tOC7nUbD7M+vXra0ZrrHp17xExIx8uCEAAAhCAQJkISARGIjH\/\/R+vMv\/v35cumIaAKZOXQraEp4ySCBib\/yJVxSXxjo2NmY6ODiIwJe4DmAYBCEAAAsbIsuo3P3xNEIWZmpoy4+PjRGDK3DHSTiElES\/SXs3zh2X2J7ZBAAIQgEBjBGwURnJhvvfqj82GDRsQMI2hzOeucMSlXhJv3Moj12IETD7+4ykQgAAEIJAdAcmF+dDNXzSvfnozAiY7rK2pKekyanm6iJuTJ09GThshYFrjI2qFAAQgAIF8CNgjBk7+\/mHz+\/3DRGDywd74U6I2srMRl76+PrNy5cpFm9jZp61atcrs37\/ftLe3LzKACEzj\/uBOCEAAAhAojoBEYV5s+4Tp\/847CJji3FDck30SMDMzM+bhhx82t956qxfLvLG3tf0avvC1BOgL1ewLdl+Yzx+9GgHT2i5Qztp9EjA+2Srext7W9nn4wtcSoC9Usy\/YM5L2TH\/EfPH+R9kHprXdoHy1+\/Ti+2QrAqb1fZ3+0FrGPvH1yVa+G7Ltt3JG0rf\/9nvm+l1PImCyRVv+2uyL7+7EW1ar7Y6LPtgqDLG3tT0JvvC1BOgL1e0LH37lb8yF3\/4js+SeZ03X9Te2FkTOtavZibdV3GTuWLZhls2AuCAAAQhAAAK+Efg\/n\/5hsKT6I3cc8M30uvYiYBK4U0SMfLggAAEIQAACvhH46IXvqYu+iA8QML71ROyFAAQgAAEIQAABQx+AAAQgAAEIQMA\/AkRg\/PMZFkMAAhCAAAQqTwABU\/kuAAAIQAACEICAfwQQMP75DIshAAEIQAAClSeAgKl8FwAABCAAAQhAwD8CCJg6PpNTq\/ft2xeUmJiYKM0uhnLy9saNG4MTtdeuXZv4VO1Dhw4lKpt1N05jr8t8+fLl5sCBA2bFihVZm1S3vjT2ugeIRh0I2krj09hq7bCHm65evdr09va20rzz6k5jr8tWKio731OnTi06KLaI74ykfMNsraN2796da59Iaq\/Y55b14bvBboJaVN\/N9cUu6GEImAjw0vlkMJUTqn\/wgx8s\/H\/4tOq8\/eYOPuvWrTNDQ0MmbiCyL1JSsZNlm9LYK1+qk5OTCyJL\/n348OGap4RnaaNbVxp73T4i\/UL6i4jKkZER09bW1ioTF+pNY6trjB288h6s0torPLu6unIdUBvtC7ZtMrAODg4Gg+2WLVvMnj17chPgafm6bQ335ZZ3XmNMGnutOBS2PT09puzfDVZs3XfffYG9wreoPyDz8GVRz0DARJCXL0+55IWxL1pfX1\/hUZjwF2PciyHteOyxx8zNN99s3nrrrdwGV4s1rb2uO4oYBJqxN+9BoBFbZSC46667zJkzZ8z69etzFQdp7C3DO5fGXim7a9cus3fvXlPUHzlp7HXfs7A4yGswSmNvuGzZvxvCf4xJf965c6e5\/fbbcxO0efmxyOcgYGrQD4fYiwy5h80LD5Jxg6b83v7F4kY38up0ae0tWsA0Y68revPg24itYuNnPvMZ8+1vfzs2cpd1G9LYW9SgWi8qUe9dCw9YWbNLUl8avuGIXNm\/G2pFYPK2OQ3fWgJGouVl+CM4SV\/ypQwCpo6AcTtb0eFsa2Y44pL0L7+ivmAbtVfam\/eUjDyzEXvtFF3e8\/JpbZW+cvDgQfOHf\/iHZseOHYUIGDeMXq\/vuvkOtu\/nnVOShq+8X9PT04GpReXNpbHXMi1SKKa11\/4hKRHlzZs3B9HxPK809oZP\/7b\/znvaNk8+RTwLAYOAaWm\/S\/PSh\/8qvP\/++3NP4m3UXit+7r777txsTmOrG8Lu7OxMlDuVdcdIY6+UdVmG\/521bbXqS2OvzSuyIqvs9rp\/ENlcv7ynvtLwrZVTkrfdaewVvm6i9IYNG4Ip\/Lh8xTz6taZnIGDqCBjb2XyeQrLNKzIC437RxE152Re\/CPFiRUhaey3jvPtJmpC2lH3mmWcW5XTl\/WWaxt7wa5k327R9IWrKIE\/GjfC1kaO8oxlV4VuGSJcmwRJuCwImwrvulFEZEgqtmeGwe\/ivgqjOWpSASWtvEasLXGZp7XXvzTscn8ZWd3m6a3Oeofg09kYJmDxzCNLYG34Pi\/jOSGOv8C1CFDb6rpVBIKblG25r3qvSNAsX2zYETISX3b9mfF9GbaMaeSe9hb8k45Z9FxF2r\/eXfhJ73WhN3uIrzTJUt51FDVxp7PVtSX1YvCaJNGY9wKThK8+2K9LuueeeQlbGpLG31hRSntO1ab\/LXLEjWypIAq9dYp+136tcHwKmjvd93MguKiRcVARG8NbbrMq1NypKkHfyZlJ7rTAcHh4OelHeSbxp2JZBwKS1180hKIJtWnvdjex8sLeIpcjhr9s075pNhC3qXUvbH9z+W8QeXFUQNgiYKniZNkIAAhCAAASUEUDAKHMozYEABCAAAQhUgQACpgpepo0QgAAEIAABZQQQMMocSnMgAAEIQAACVSCAgKmCl2kjBCAAAQhAQBkBBIwyh9IcCEAAAhCAQBUIIGCq4GXaCAEIQAACEFBGAAGjzKE0BwI+EpA9VOSAyW3btpm8z+TxkRc2QwACxiBg6AUQgEDhBNyzmgo3BgMgAAEvCCBgvHATRkJANwHZhfmmm24yPT09uhtK6yAAgcwIIGAyQ0lFEIBAmEDUuUvuWUFyVszOnTvN7bffHpzJY+957LHHguqK2pYfb0IAAuUmgIApt3+wDgLeE6h1yKVEXOQaHBwMzso6ePCg2bp1a\/AzOfhOrpGRESPiJu9DMr0HTgMgUBECCJiKOJpmQqAoAuGTmsP\/FoEiV29vbyBmtmzZYvbs2VPICclFMeK5EIBAegIImPTMuAMCEEhJwI24uNNHsuLIzX8J\/y7lYygOAQhUiAACpkLOpqkQKIqAK0wefPBB09XVFURcwsunETBFeYjnQsA\/AggY\/3yGxRDwjoCdNlq3bp159NFHF6aIwsunmULyzrUYDIHCCCBgCkPPgyFQLQKS6zI8PGzWrl27KEFXKEg0Ri67AklWHkmCr1yImmr1E1oLgaQEEDBJSVEOAhBoioAIkY0bN5o777wzECwiVtzl07ZyllE3hZmbIVAZAgiYyriahkIAAhCAAAT0EEDA6PElLYEABCAAAQhUhgACpjKupqEQgAAEIAABPQT+P\/sfAjj9AaPUAAAAAElFTkSuQmCC","height":337,"width":560}}
%---
