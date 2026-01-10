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


%psin = simplify(sqrt(1 - dot(vp,vp))*sqrt(dot(r_vec_l1,r_vec_l1)));
psin = simplify(sqrt(1 - dot(vp,vp)));


%%
%[text] ## The Relative Net Force of the Two Points
% Force of the two point in plate and moving particle
dforce = (psin*vec_l1 + psin*vec_l2 )/(2*R^2); % net Force in the x direction first points.

dforce = simplify(expand(dforce),'Steps',100) % Two charges %[output:16ef0998]
dforce_x = dforce(1) %[output:1862f7b5]
%dforce_x = simplify(taylor(expand(subs(dforce_x/2,v,v+v_a) + dforce_x/2),v_a,order=2))
%%
%[text] ## Check at v = 0, r =0
df_v_zero = simplify(subs(dforce_x,[v,r],[0,0]),'Steps',100) % Due to the two point charges %[output:3d0726da]
df_v_one = simplify(subs(dforce_x,[v,r],[v,0]),'Steps',100) % Due to the two point charges %[output:69e50e72]
%%
%[text] ## Integrate all the angles
%dforce_xT = expand(taylor(dforce_x,v,Order=taylorOrder));
dforce_xT = dforce_x;

dcircle = int(dforce_xT,a,0,pi); % From zero to pi i.e., the full circle
%%
dcircle = simplify(dcircle,steps=100); % Net force of the two full circles with charge density
df_v_zero = simplify(subs(dcircle,[v,r],[v,0]),'Steps',100) % Net force of the full zero radius circle (point) %[output:4bf04a0c]
%[text] ## 
%%
%[text] ## Integrate for all values of r
syms R_o real
assumeAlso (R_o>0)
dtotf = r*dcircle/pi  %[output:7b7123d6]
dtotf = simplify(expand(taylor(dtotf,v,Order=taylorOrder))); % Taylor
%dtotf = simplify(subs(dtotf,v,0)) %static
totInt = int(dtotf,r,0,R_o); % Add all the circles
%%
totf = simplify(totInt,'Criterion','preferReal','Steps',100) % The force on a charge due a charge surface density %[output:3b43bbb2]

%%
%[text] ## Set the Plates Radius to Infinity
sympref('PolynomialDisplayStyle', 'ascend');
totfInf = simplify(limit(totf,R_o,Inf)) % Net force by plate charge density %[output:34b522cd]
%totfInf = subs(totfInf,v,-v)

%%
figure(1)
plotrange = [0,0.99];
fplot(totfInf,plotrange) %[output:3800c3f5]
hold on %[output:3800c3f5]
relEq = sqrt(1-v^2) %[output:4739c8cb]
relEqApx = taylor(relEq,v,Order=taylorOrder) %[output:38ff51e6]
fplot(relEq,plotrange) %[output:3800c3f5]
fplot(relEqApx,plotrange) %[output:3800c3f5]
xlabel('v/c')  %[output:3800c3f5]
ylabel('ratio')  %[output:3800c3f5]
legend({'Force Carriers',"Lorentz Mass","Lorentz Apx"},'Location','northeast') %[output:3800c3f5]
hold off %[output:3800c3f5]



hold off %[output:3800c3f5]
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
%[output:16ef0998]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce","value":"\\left(\\begin{array}{ccc}\n\\frac{d\\,\\sqrt{1-v^2 }}{{{\\left(r^2 +d^2 \\right)}}^{3\/2} } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:1862f7b5]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce_x","value":"\\frac{d\\,\\sqrt{1-v^2 }}{{{\\left(r^2 +d^2 \\right)}}^{3\/2} }"}}
%---
%[output:3d0726da]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\frac{1}{d^2 }"}}
%---
%[output:69e50e72]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_one","value":"\\frac{\\sqrt{1-v^2 }}{d^2 }"}}
%---
%[output:4bf04a0c]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\frac{\\pi \\,\\sqrt{1-v^2 }}{d^2 }"}}
%---
%[output:7b7123d6]
%   data: {"dataType":"symbolic","outputData":{"name":"dtotf","value":"\\frac{d\\,r\\,\\sqrt{1-v^2 }}{{{\\left(r^2 +d^2 \\right)}}^{3\/2} }"}}
%---
%[output:3b43bbb2]
%   data: {"dataType":"symbolic","outputData":{"name":"totf","value":"\\frac{{\\left(\\frac{d}{\\sqrt{{R_o }^2 +d^2 }}-1\\right)}\\,{\\left(-128+64\\,v^2 +16\\,v^4 +8\\,v^6 +5\\,v^8 \\right)}}{128}"}}
%---
%[output:34b522cd]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInf","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}-\\frac{v^6 }{16}-\\frac{5\\,v^8 }{128}"}}
%---
%[output:4739c8cb]
%   data: {"dataType":"symbolic","outputData":{"name":"relEq","value":"\\sqrt{1-v^2 }"}}
%---
%[output:38ff51e6]
%   data: {"dataType":"symbolic","outputData":{"name":"relEqApx","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}-\\frac{v^6 }{16}-\\frac{5\\,v^8 }{128}"}}
%---
%[output:3800c3f5]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQtsXcd55z+rcmz6odi0ldYknaW7lV0gyMp2tkuC2cIy4NRtUapobYeSdlFb0GbVTeIyrcyXjEaSN5YoUsKasJ1CawiKioKS\/EJrFYvYXsBJnXJNbJxE2BqtoyZhbZLpxrGkxlEYrV158R1l2OHhfZy5957HzP0dgJBIzpn55vfN3Pnzm9dF77\/\/\/vvCAwEIQAACEIAABDwicBECxiNvYSoEIAABCEAAAhEBBAwNAQIQgAAEIAAB7wggYLxzGQZDAAIQgAAEIICAoQ1AAAIQgAAEIOAdAQSMdy7DYAhAAAIQgAAEEDC0AQhAAAIQgAAEvCOAgPHOZRgMAQhAAAIQgAAChjYAAQhAAAIQgIB3BBAw3rkMgyEAAQhAAAIQQMDQBiAAAQhAAAIQ8I4AAsY7l2EwBCAAAQhAAAIIGNoABCAAAQhAAALeEUDAeOcyDIYABCAAAQhAAAFDG4AABCAAAQhAwDsCCBjvXIbBEIAABCAAAQggYGgDEIAABCAAAQh4RwAB453LMBgCEIAABCAAAQQMbQACEIAABCAAAe8IIGC8cxkGQwACEIAABCDQ9AJm79690tnZKX19fbQGCEAAAhCAAAQ8IdDUAkbFy4EDB2TPnj0IGE8aLGZCAAIQgAAElEBTCphTp07Jli1b5Kqrropawa\/\/+q8jYOgPEIAABCAAAY8INK2AOX36tLS1tcnw8LD09PQgYDxqtJgKAQhAAAIQaEoBY9y+sLCQSMD86K83Rq\/83GUdZVvMisva5fxP5qLf6\/9LPSus91e0LM3r+6cuvDH\/s39pmhCAAAQgAIFGEOjo6BD9Cu1BwFSJwMzOzsrfPvmrif1+XevypG0lfpYkw1f\/XsTkFxc45vsLoucisX9v\/y5JOaSBAAQgAIFwCXR1dcn4+HhwIgYBU0XAvPLKK7Jp06bI+e3tpSMrrs0+Lmhs0XN+YTbK7vxPLvx74f8XIjv\/\/LOflUpTygYTzVFBo9EfO+pz8bXd0SsXft6YerlymJ6elomJiYaydbXBJT32utByTwtfd2ZJ34BtUlK1pSsyX2Pb5OSkdHdf+NwP5UHAJBQwRXa+ETtTLz0TCYLRXf3S1vr+MvGjwscWRnYjNtNbKnJ0qsxMgxnhc+FnjQ1BGnFYZLY2I+xN92MPvunxhW16bDXnIvMtsm31egUBE4CAMY3ApaEaIaNRHVvYaLRHf\/be26+UbFtGxFx8TfeiyFl5zQVVb6I6SRuli71J80wzHfamSbfYg0CpmvvUHnyyteiCwLe24JvvXT5lEDABCRhdr\/PMM8\/IXXfd1bC5TiN03n37lcXojRE5pSI6diTn4mu7oraoAqdUBCcNe10av2ta7HUl5pYevm68XFLD1oWWe9oi80XAuPszmDdCdn6jnKQix47kVIrimCmpauKmUbaRDwQgAIFmJhDyGNbUEZgkjTpk5yepfyPSGIFjpqXe\/eF0ySmqUuLGdVqqEfaSBwQgAIFQCIQ8hiFgqrTSkJ1fhA6aRNwgbIrgKWyAAAR8JBDyGIaAQcAUsk+6Cpt\/ETlhbRMspHMwCgIQ8IYAAsYbVzXe0JCd33ha6edoFhWfe\/PpqLBS01EqZi65\/q7o97qAmGmo9P1CCRCAQDEJhDyGEYEhAlPMXudolR2x0UXERuCYbFTUmK3fiBpHuCQvLAHd\/aJfPM1HIOn1AAiY5msbizUO2fmhu7VatAZRE3oLCLt+KlwGBgZET1rlaT4CSa8HCHkMIwJDBKbper4Km3JTUPb00yXX393w04ebDjYVTo1AGtecpGYsGTeUgMv1AAiYhqL3K7OQne+XJ9Kz1o7UxNfU2DugmHpKzwfk7E6AzyZ3ZqG84eJ7l7S+8SECQwTGtzabib32mppSokYXCV+I1tydiT0UAoE4gZAHJrxdmYCL713S+sYdAYOA8a3N5mavmXoqJWjMAmGmnXJzT9MV7PvAdOrUKdmyZYucOHFiie96e3tldHRUWlpaMvWp4WkK3bNnj\/T19dVlw8LCggwPD8vGjRsbehO0i+9d0tZV2RxeRsAgYHJodmEUaa+lOffmM4t3RdmLg1tu+lwYlaUWhSPg+8BkBMzQ0FBDB\/daHKUsH3jgATl06JCsWbNGjPDo6empW8TUYk+1d1x875K2WrlF+z0CBgFTtDbprT0qaMyll3FBY86lQdB4697CGe77wFRNwBgRcfz48Yj95ORkJHROnjwpg4ODi\/44ePBg9H8TzVm7dq3oz1pbW6O0mzdvlvn5ebF\/bjuznFjRd3fv3i379++P8jp27JiMjIxEr7a1tUViR7cya4TlzJkz8vLLL0c\/e\/bZZ5d9byIw5exRX+7du1feeust+djHPiaf\/\/zn5aGHHhJT93g0yMX3LmkL18irGISAQcD41ma9sRdB442rvDTU94GpmoDRAV2Fh04n6TSTiZCos1SU7Nu3LxI08WkaFRozMzPyqU99KhI1JsJj52dPTxlBNDY2FkVfSj1GYBhhpHnp8wd\/8AeRgFFBo+UYW+Lfq4C58cYby9pj109tMHXQPJXTrl27ZMeOHZGQ0sfF9y5pfesICBgEjG9t1lt7bUGz8PrEYj3M1m12OXnr2lwMjw9Mb5z6qehXUZ8Pt14q+mWeUmtgTJREBYa9dsSOktx6661RBMYIjnikxOQfFx3l0pX7eSWORmAYAWOmmuLRHFtcaX4qfOzokInwfPvb317yO1vAlBNUmzZtWoxKVbIVAVPUHpGBXSE7PwN8FFGBQLlFwReETFe0w0nFDQ8Ekgxie5\/\/nux9fqawsIbu7JShO29YJmBKrYEpFZ3Rwb+zs1NUwNhTO3GhYgsYHeTtx0z92JGWJBGY+HSW5rl169bFCIyZIopHg+ICppw9b7\/9thw5cmTJ4mWt74EDByLzzfRZvG7xnydpJ4VtIDUYRgSmCjQETA2tildqIlBK0BCdqQllU7wUSgSmlIApJwQ00hEXMJUiMHFRUKphlFsDY0\/dvPjiizI1NbUoMOIRmKQCppw96styv1M7tm3bJtu3b1+c4nIZl1zS+tZxEDAIGN\/abFPYW266yURnWAzcFM2gYiV9H5jqWQNjR2Di+ai4ULGh00yf\/exnF9fA6M+PHj26OIVjwy23C8msZTF56nocFTy6tkbX35gppCQCJr4GxrZHp5BsAWOiTbqNmzUw5bsBAgYBw0jgAQGiMx44KWMTQxcwlXYh2QJGsZfb3WP\/vNT0UVzE2FM8OkWk0SF97PU6uk5HFwg\/\/\/zzi7uFkggYs4PK7Iqy7YlHYOLrg5hCKt25EDAImIw\/dimuXgL2+TNmMbCZauJ04Hrp+vO+7wLGH9LFs9TF9y5pi1fTyhYhYBAwvrVZ7I0RMILG3tnEVFP4zSTkgSl879VXQxffu6Stz6rs30bAIGCyb3WUmBoBI2ZKHaTHupnUsOeSccgDUy5APSrUxfcuaT1CEJkalICxT0qsdo+FvUXNnuuMOzBk5\/vWWLHXjUCldTOIGTeWRUzNZ1MRvZKNTS6+d0mbjfWNKyUYAWPv5Vc89kFHcVzx1ej2im8ETOMaFzkVhwBipji+aJQlIQ9MjWIUaj4uvndJ6xuvYASMvc1NT3EsJ0pK7flX8XP48GF58MEHl92AGrLzfWus2NsYAoiZxnDMOxffP5uqbaNOk6\/r6btm3Hj11VcXL3w09pnIf7WofyPr4+J7l7SNtDGLvIIRMOZuCrPtLf69gVlOwMS35Zn0xvn9\/f3S1dUV\/Vgv8NIvHgj4TgAx468HfR+YfBQwemnjb\/7mby7eUK3jycMPPyyvvfaabNiwIbObq6v5fnZ2VvRLn7m5ORkYGEh07YBvvSEoAaPHTOvBP\/pUuktCxY02APtiLr3101ylbjvRNBT7Zypm9IsHAiERKCVm2M1UXA9XG8SKa\/kFy6oJGPuzt7e3NzoFVx\/79mc9H0XPZdGflbq1Wv8wXbVq1ZLfmQPl9AJFfVfz\/sIXvrAEVzyaYv7wvfnmm+U73\/nOYrTeRO\/15Y9+9KPR+BM\/w8VeY1lu7WXSNZnxP6zLXSUwMTEh+mU\/Sa4dKHqbidvXlAImfkDSH\/\/xH8tzzz23eDlYKQEzPj4u7e3t0a+IwPjWzLHXlUB8NxNXGrgSTD99yALGHECnN04bgaIHv8Vvf1bKdrTdPlFXf6eHxt1\/\/\/2RsLBvo9boRKmoe6lj+zUfM2bceeed8tRTTy0e669\/KOujt1\/rH9Dr16+PxJS53NFem2nfd6TvaOTm3nvvlXI\/L3cztr5bzfd2BGZ6ejoSMwiY9PtkzSUknUIqVYA2BvuW0FICJkTn1wybF5uKQPycGSNmuGwy32YQH8Tee2tG3v1BcS9zvPhDnbJydecitEoRGK2bfbS++Yx+7LHHoj80jUCI52HfoXTNNdcs2cxh51lKwFSyx873q1\/9qtx2222RsDIi5Nlnn40EjJkBMJW0BVGpCxuNGElyZ1Ot41I1sZNvK66v9GAiMPEpo0o7i+LI4guAa20o9bmCtyFQfALlxAzbsrP3XXxgOv3kLjn91M7sDUlY4tX37JSrP7kjkYCJfyabRbc61aMCxhzdH5+uMZnrFFCpW6uNUIgLmHIXOpr84rdKq4j53d\/93cUozhNPPLEoYOLLDuwrA8od9eFyBEiSCEyzjGHBCBjXbdQa8tMFv3aoUu+qiD8hq9eEnzskg0BJAguvPyLn\/uGv5fxPvx79XiMyK6\/tiv7lSZ9AM0dgbAETv6nZkI\/vNCoXgdFdq\/a0TynP2QJG19Ds2rUruhn67Nmz0Thi\/mD+xCc+EV30aG7YrjYlZepRSiiVGo9MOpdxySVt+q22sSUEI2AUSzkVW+5qdrPoq9L2t5Cd39imRG7NSqDcehmmmNJtEb5\/NlWasqm2BsYe+O3lA\/Z7OoVkr3MpJ2A0eqKP2cFaTcCosDCLbs3SgnICRsekRx99NNog8o1vfCNaK6PlmN1Lugam3M\/rWQNDBCbdvudN7r5\/SHgDGkODIMAUU3Zu9P2zqdz0jxEFlXYh2QImvinD\/EFaKQKj72ik5Pz583LdddfJCy+8sMRx8dPZ438EVzoM1f5DWjeIfOtb34qmvOK7pYyd5eyv1JJcfO+SNrvW25iSgorANAbJ0lxCdn4avMgTAkqAqEz67YDPpvQZF7UEF9+7pC1qfcvZhYCp4rGQne9bY8VePwkQlUnHb3w2pcPVh1xdfO+S1oe62zYiYBAwvrVZ7PWYQLTw981noggN27Hrc2TIA1N9ZMJ\/28X3Lml9I4eAQcD41maxNwACpaIyLTf1s4PJwbchD0wOGJoyqYvvXdL6BhMBg4Dxrc1ib0AEosPXTv35sqgM58pUd3LIA1P12jd3Chffu6T1jSoCBgHjW5vF3kAJaFTm7Df2R4JGHxOR0akmnuUEQh6Y8HdlAi6+d0nrG3cEDALGtzaLvYETYHopmYN9H5iqXeaYjEJtqeJbrF1yUe4PPPBAyct\/XfKpJ62L713S1mNTHu8iYBAwebQ7yoRAIgKlFv0yvXQBne8Dk48CxhxAp\/zN7dOJGnKDE7n43iVtg81MPTsEDAIm9UZGARCol4AdlWH3UnMImEoH2Z05c0Zefvnl6Ibl+AFx5iA8E2VZtWqVmFPX9Xd6FYAeYnfixIno3d7eXtE7luyn3Onsmufhw4flnnvukUceeUT2798vra2ti7dV33zzzXLw4EGZn58Xk4d9L1\/8QLxa+4WLKHFJW6s9eb2HgEHA5NX2KBcCzgTi00vNvE7G94GpnqsE9IJEc\/S\/fZWAPb2jjWvz5s1y\/\/33R7dEazoVFqOjo1LqNmpNX+7uItNQ9ZRdfUx+eiu1Xi1ghMkbb7wRCRi9eXpwcDC6eLKjo2Px1mr9uV4EWenagiSdwsX3LmmTlF2kNAgYBEyR2iO2QCARAdbJLJ9CUib\/\/JPZRPzySPRzl3VEZ\/+Yp5KA0UHX3Bytly3q9ypAHnvssUgU9PT0RCIinocd4dC7kIyI0HuF7DxLCZhqU1r2\/UUmPyNGSt1mbUdeNHKjUZ4VK1YsRm3q8YGLKHFJW49NebyLgEHA5NHuKBMCDSPQrOtk4gOTclh4faJhXBudkUbL7PVLlQSDRjqmpqaiaIkKGDMdpCJABYx9G7WZDrLt1embW2+9NdFljvYUkBFGpepuLorUKI55dApKIy7mRuv4JZOdnZ2R0NLHFjT1snURJS5p67Ur6\/cRMAiYrNsc5UEgFQIagdAB\/NybTy+e8hvygt\/4wNRMERhbwGzbtk22b98u8dubK13maEdgjPioJF6MANF\/7ekfI0rWr18vw8PDi5GheERGffWXf\/mXUbvXG6gr3TSdpHO4iBKXtEnKLlIaBAwCpkjtEVsgUDcB+zwZs+A3RCHj+8BUzxqYeKTDCAsTJdm3b5\/oFNLu3bsXp2zKTSE98cQTUZurtC6lnK0mUvT5z39eHnrooSgfs8bGTF+pHbt27ZIdO3ZEa2N0EfCDDz4YRW1qfVx875K2Vnvyeg8Bg4DJq+1RLgRSJVBq51JIQsb3gcmIAt0NZD9mF1GlXUi2gDHRDrPTyOz+qRSB0Xd06un8+fNy3XXXyQsvvLDEhq1bty4RNGYNjk4X6ZSTeYxgevjhh+XZZ58VXVx84MCB6Nf2DinbXhU9MzMzdS3kdfG9S9pUO2QKmSNgEDApNCuyhEBxCJTauRSCkAl5YCpO60lmSaO2Rycrze0MoJDbCQIGAZO0z5AOAl4TCG0LdsgDk28NDQGTj8cQMAiYfFoepUIgJwKhCBkETE4NqADFuvjeJW0BquZkAgIGAePUYEgMgVAImJuwzdbj+DbfotfTDEz9\/f3S1dVVdHOxr4EE5ubmZGBgIFpnowfpVXoQMA0E71tWITvfN19gLwTSIODrGhndCqyD2PT0dBpYyLPgBFS0jo+PRyf9ImAK7qy8zEPA5EWeciGQLQEfhYyKGP3iaT4CKlyqiRelEvIYFtQUkm5PGxkZiVpyucu4TDO30+plXubEx3g3CNn5zdflqTEEqhPwUchUrxUpmpVAyGNYMAJG9+Obg4O0odp3YJQSJXqCon0EtH05mJ0+ZOc3a4em3hBIQiD0c2SSMCCN\/wRCHsOCETDxuzMq3TsRTxv\/HgHjf6elBhBoFAGETKNIkk8eBBAweVB3LNO+Ul1fjX8fFyXxCEy5ezBKrfRPOvfoWAWSQwACBSagQuZHX\/tDOf\/Tr0d3LemupUuuv7vAFmNasxKw10a57FjyjVcwEZh4xKXacc32zaKVtqLZx1kb5+q2Rf3igQAEmo+ACpkzL9wjsuIfIyFz+c3jcvG1lbeyNh8lapwngYmJCdEv+0my5TpPm2spuykFjIqbo0ePRmtg9F6LatGaTZs2RdvV2tvbI8ZEYGpparwDgbAILPztUTk397iooFl5Tbdccct4JGh4IJA3ATsCo9vsVcwgYPL2SoXyk04hxa851yztBcDxa85Dnj8ssDsxDQLeEFh4\/RE5++p+WXHFymhaKYR7lryBj6FVCYQ8hgUTgYlPGZVbxIuAqdreSQABCDgSsE\/11SjMJdffhZBxZEjydAggYNLh2tBcXbZRl5pCmp+fL3kWTMjOb6gDyAwCEIimk85+Y7+8e+rPWehLeygEgZDHsGAiMNpSyh1kV+qmUI3QHDhwIGpgHGRXiH6GERAIhoAKmR9\/c0Dee\/uVSMis6jnC+phgvOtXRRAwfvmrodaG7PyGgiIzCEBgGYF3\/mpUzs0+zvoY2kZuBEIew4KKwKTRQkJ2fhq8yBMCEFhOQBf66q3XrI+hdWRNIOQxDAFTpTWF7PysOxLlQaCZCcTXx3B+TDO3huzqHvIYhoBBwGTXkygJAhCQc\/\/wNTn7zYHoIDw9P2bVx49ABQKpEUDApIa2+BmH7Pzi08dCCIRLQNfHvHvmwkYCzo8J18951yzkMYwIDBGYvPsX5UOgaQno+TE\/evkPRVZ+i2sJmrYVpFtxBEy6fAude8jOLzR4jINAExFgWqmJnJ1xVUMew4jAEIHJuDtRHAQgUI7A2Vf3RfcrsVuJNtIoAgiYRpH0MJ+Qne+hOzAZAsETiG67\/vL9i9NKHIIXvMtTrWDIYxgRGCIwqXYeMocABGojcPZ\/\/5mc\/daArPz5S1nkWxtC3hIRBEwTN4OQnd\/EbqXqEPCCAIt8vXBToY0MeQwjAkMEptCdD+MgAAGJzo75pxfuIRpDY3AmgIBxRhbOCyE7PxwvURMIhE9AozHn3nh6cZEvJ\/mG7\/NG1DDkMYwIDBGYRvQR8oAABDIiEEVjvny\/rLzuDCf5ZsTc52IQMD57r07bQ3Z+nWh4HQIQyJGAuel65Yc65ZLr75KWmz6XozUUXVQCIY9hRGCIwBS132EXBCBQhYC9yJd7lWgupQggYJq4XYTs\/CZ2K1WHQFAE7GgMa2OCcm3dlQl5DCMCQwSm7g5CBhCAQP4EWBuTvw+KaAECpoheycimkJ2fEUKKgQAEMiRgbrnW6wiIxmQIvqBFhTyGEYEhAlPQbodZEIBArQQ0GvP2n\/bKpWuvkkuuv1suv2W81qx4z3MCCBjPHViP+SE7vx4uvAsBCBSbgC7wPfX0f5GV7X8XXQ5JNKbY\/krLupDHMCIwRGDS6jfkCwEIFIAAdyoVwAk5moCAyRG+S9HHjh2TkZGR6JU9e\/ZIX1\/fstcXFhZkeHhYjh8\/vuR3bW1tcujQIVmzZs2Sn4fsfBe2pIUABPwloNEYveF6xQf\/hmiMv26syfKQx7BgIjAnT56UwcFBGRsbi5xs\/h8XJKVawN69e6MfDw0NLft1yM6vqTfwEgQg4C2Bf\/ofX4iuIuCGa29d6Gx4yGNYMAJGoy9TU1MyOjoqLS0toqKks7OzZBTGbgHqXE178OBBaW1tRcA4dw9egAAEfCKw8NpX5PSzvxMt8OXwO588V5utCJjauGX6VjyKUimqYgwz00k9PT1lhY5xfn9\/v3R1dUWvdnR0RF88EIAABHwkYBb46pSSXkXAAl8fvVje5tnZWdEvfebm5mRgYEAmJyelu7s7qIoGE4GJR1w0IjMzM1NyWsh4sFr0RdMZAWN7XcWMfvFAAAIQ8JkAU0o+e6+87RMTE6Jf9oOAKbCvaxEwSaI0RsCMj49Le3s7EZgCtwFMgwAE3AloNOaHf9orH1jz42iB76qeI9G\/PP4SsCMw09PTkZhBwBTYn65TSGb6aOPGjRXDaiHPHxbYnZgGAQhkTMCOxlzZc0Quvjas6YaMcRamuJDHsGCmkOJTRtUW8dq7lirtVArZ+YXpYRgCAQgUgoC9wJcTfAvhkrqNCHkMC0bAuG6jVqceOXJkcddSuVYSsvPr7hlkAAEIBEdAp5R+8Ce\/IZfe\/FOmlALwbshjWDACRttZuYPsSk0XxbddI2AC6KlUAQIQaBiBfxz\/hKy46m\/kAzf8EruUGkY1+4wQMNkzL0yJITu\/MJAxBAIQKCSB00\/ukoXXH4nOjGm5qV9abvpcIe3EqPIEQh7DgorApNGIQ3Z+GrzIEwIQCIuArot5609+Q674tV+IppSuuuPlsCoYeG1CHsMQMFUab8jOD7zfUj0IQKBBBKKt1l\/6FFNKDeKZZTYhj2EIGARMln2JsiAAAY8J6JTSubnH5AP\/+gqmlDzxIwLGE0elYWbIzk+DF3lCAAJhE3jnK1+Sf\/ry\/XJZz7XcpeSBq0Mew4jAEIHxoAtiIgQgUCQCOqU0O7hGVv1OB1uti+SYErYgYAruoDTNC9n5aXIjbwhAIGwCrIvxw78hj2FEYIjA+NELsRICECgkgfmdt4veas26mEK6Z\/FCYu5CKqZ\/UrUqZPWaKjgyhwAEmoaALu49+419rIspoMdDHsOIwBCBKWCXwyQIQMA3AvHzYrjVuhgeRMAUww+5WBGy83MBSqEQgECwBMziXt2hxBUExXBzyGMYERgiMMXoZVgBAQgEQSC6DPLxzdG6mEs\/+svReTF6szVPPgQQMPlwL0SpITu\/EIAxAgIQCJKALu49\/9OvR+tiVMBcfst4kPUseqVCHsOIwBCBKXr\/wz4IQMBTArq4952\/Go3uUVp5Tbes+vgRT2vir9kIGH99V7flITu\/bjhkAAEIQKAKAT259+3Dn1o89I7LILNtMiGPYURgiMBk25soDQIQaDoCukPp\/+77RDSdxOLebN2PgMmWd6FKC9n5hQKNMRCAQNAEVMR8f+ftF6aTfv5SubLniFx8bXfQdS5C5UIew4jAEIEpQh\/DBghAoAkI6A6lNz59g1z6b66SS9dexY3WGfgcAZMB5KIWEbLzi8ocuyAAgXAJmG3W7FDKxschj2FEYIjAZNOLKAUCEIDAzwgYEfPeD19hh1LKrQIBkzLgImcfsvOLzB3bIACB8AnoWTH\/7x++Jld84hdk5Yc6hR1Kjfd5yGMYERgiMI3vMeQIAQhAICGBtx7fLGe\/\/mfsUErIyzUZAsaVWE7pjx07JiMjI1Hpe\/bskb6+vrKWGKdqgrVr18rBgweltbV1WfqQnZ+TmygWAhCAwBICeuDd6ad2RtNJbLNubOMIeQwrVATGFiDGhdWEiEl38uRJGRwclLGxsehH5v9r1qxZ1ho07ebNm2Xfvn3S3d0tWu7U1JSMjo5KS0vLkvQhO7+x3YTcIAABCNROQA+802gM26xrZ1jqzZDHsMIIGBURR48eXRIJOXXqlGzZskU2bNhQMZqiTouLkL1790pnZ2fJ9zTtzMyMDA0NVW0pITu\/auVJAAEIQCBDAkbErPrtfysrVv2Qs2IawD7kMawQAsYIFRUUGhGxH4WvYqTcFI9Jq2n0MaIk\/r1Jt7CwIMPDw9LT01NVFOk7xvn9\/f3S1dUVZdPR0RF98UAAAhCAQGMJLIqY9b8V3Witl0Bym7Ub49nZWdEvfebm5mRgYEAmJyeXja9uuRYvdVACxo64lIuyGAFz5513yhNPPCEnTpxItAbGdp2KGf2zYG9+AAAe+ElEQVTigQAEIACBxhMwIubqvt+X9z\/wZW6zdkQ8MTEh+mU\/CBhHiC7J651Cik8ZVRMwb7zxxmJUR9+dn5+vuAZmfHxc2tvbicC4OJW0EIAABGokYK4eMCKG26yTg7QjMNPT05GYQcAk51dTynoW8dYzhWQvAI4v+g15\/rAmJ\/ESBCAAgYwIIGLqBx3yGFaIKaT6XXRhEa+9MLfSIt7471TA7N69W\/bv379sK3XIzm8Ed\/KAAAQgkCYBI2Iu+3d3yAdu+nshEuNGO+QxLBgB47KNOr4wuNyCX20mITvfrRuQGgIQgEA+BGwRc1nPtZERqz5+JB9jPCs15DEsNwFj7zy68cYbo+3SuqC21FPpoDk7fbmD7MzC3Y0bNy6uwrYPsuvt7S25\/gUB41lPxVwIQCBYAkbEfPC3\/pNc\/Es\/RsQk9DQCJiGoEJOF7PwQ\/UWdIACBcAkYEXPN5v8qF33wNTm\/MMv9SVXcHfIYllsExmbeiHNg0uqyITs\/LWbkCwEIQCAtAogYN7Ihj2EImCZWr27dgNQQgAAEikHAiJjrdr4k7545EEViLr95XC6+dulBqMWwNl8rEDAp8S+1bbpUUVu3bk107H8aZobs\/DR4kScEIACBLAiYw+4QMZVphzyGFT4Ck0VHqFRGyM7Pmy3lQwACEKiHACKmOr2Qx7BCCJjqLsgvRcjOz48qJUMAAhBoDIHTT+6S00\/tFCIxpXmGPIYhYKr0oZCd35iPD3KBAAQgkC+Btx7fLBqNQcQs90PIY1hhBIweRLd58+boTqL4k\/QcmDS6UMjOT4MXeUIAAhDIgwAihghMHu1OzEFzPT09sn79ehkeHhY9dM4ccDc0NJTbNeAImFyaBIVCAAIQcCYwv\/N2ee8HM7L6M4fYnfQzeiGPYYWIwMTPgbHvKlL4R44cKXtSrnMLd3whZOc7oiA5BCAAgcITMCLmw1\/8nvzorzc2\/RbrkMewQgoY+2LG+L1FWfeekJ2fNUvKgwAEIJAFARUx+rTtfKnpRUzIY1ghBIw2NPtCRVu0vPjiizI1NUUEJoteTxkQgAAEAiDw3lszMr\/jdln5oc6mFzEImAwatL0Opq+vLxI0Bw4ckLa2Njl06JCsWbMmAyuWFxGy83MBSqEQgAAEMiCgIuaNT98gl35k3RIRc9UdL2dQenGKCHkMK0QEptRt0UVxf8jOLwpj7IAABCCQBgFz5cDV9+yUqz+5Y3E6qZlETMhjWCEETKXLHNNo1C55hux8Fw6khQAEIOAjAfvepJaPrGs6ERPyGFYIAaOdQhfu5rnWpVzHDNn5Pn4YYTMEIAABVwL2ab3NJmJCHsMKIWBMBObEiRMl2yUH2bl2V9JDAAIQgIBNwD7orplEDAKmiftByM5vYrdSdQhAoAkJ2GfEaPX1nBh9Vn38SLA0Qh7DChGBKXLLCdn5ReaObRCAAATSIGCfEXP+J7Py428OBC1iQh7DEDBVekjIzk\/jw4E8IQABCBSZQPyMGBUxZ\/7nr8rKa7qDjMSEPIYhYBAwRf6swTYIQAACDScQPyPm3R++Iu9MbZRLrr9bLr9lvOHl5ZkhAiZP+jmXHbLzc0ZL8RCAAARyIxDfXm1ETMtN\/dJy0+dys6vRBYc8hgUVgdGt2CMjI5H\/9+zZI3qib6nHHJx3\/PjxxV9v3bpV9Nbr+BOy8xvdUcgPAhCAgE8E4turjYjRKIxGY0J4Qh7DghEwJ0+elMHBQRkbG4vanPl\/qSsIdNv2tm3bZPv27VWvKAjZ+SF0TuoAAQhAoB4C8Z1JRsRc2XNELr62u56sC\/FuyGNYMAImfhCe3qXU2dlZMgqjYmf37t2yf\/9+aW1trdjIQnZ+IXoXRkAAAhDImYC9M0lNOffm03L2mwMSgogJeQwLRsDYt1lrA4x\/b\/cP+7brpAKmv79furq6omw6OjqiLx4IQAACEPCfgFnUa+5M0hotvP6ILLw+4aWImZ2dFf3SZ25uTgYGBmRyclK6u\/2PKNmtLSgBY0dcNCIzMzNTcl2LvVZGYVQ66deoVxuaihn94oEABCAAgTAIxBf1aq00CqPRGN8iMRMTE6Jf9oOAKXA7jU8ZVRIwmnZ+fl5GR0elpaUlitbY38ejNZs2bZLx8XFpb28nAlPgNoBpEIAABOohEF\/Uq3npab3nF2ZlVc8RWXGZH5F3OwIzPT0diRkETD0tI+V3XaaQ4qbYC4Dji35Dnj9M2SVkDwEIQMA7AvFFvbaIueqOl72rT8hjWDBTSPGIS6VFvKUETLlFvSE737ueiMEQgAAEMiDwxqdvkJUf6pS2nS9Fpfl85UDIY1gwAibpNmpzBkxPT0+0Q8l839bWxjkwGXwwUAQEIACBohMw62HsRb2+XjmAgCl6a\/uZfeUOsjMiZePGjdEq7PhBdr29vYvrYeJVDdn5nrgVMyEAAQhkTqDUol4fT+sNeQwLJgKTVusO2flpMSNfCEAAAiEQ0PUwP33tK3Ldzpek5SProir5dtBdyGMYAqZKLwvZ+SF8wFAHCEAAAmkSKLWo16czYkIewxAwCJg0+z55QwACEPCagDnk7sp198nqzxxarItur37v7VcKf0YMAsbr5lef8SE7vz4yvA0BCECgOQiUWg+jNTdnxBR5e3XIYxgRGCIwzfEJRC0hAAEI1EGg1FSSETH676qPH6kj9\/ReRcCkx7bwOYfs\/MLDx0AIQAACBSIQPx9GTSv69uqQxzAiMERgCvTxgCkQgAAEikug3FRSkbdXI2CK255Styxk56cOjwIgAAEIBEag1H1JWsWibq8OeQwjAkMEJrCPF6oDAQhAIF0Cuh5GH3PVgCmtiNurETDptoVC5x6y8wsNHuMgAAEIFJSA2VptXzVgTC3azqSQxzAiMERgCvoRgVkQgAAEikvgna98Sd56fPOSU3ptEaP\/L8LOJARMcdtQ6paF7PzU4VEABCAAgYAJlNtaXaSdSSGPYURgiMAE\/PFC1SAAAQikR8BMJV36kXXL1sMUZWcSAiY9\/xc+55CdX3j4GAgBCECg4ATKba1Ws4uwMynkMYwIDBGYgn88YB4EIACBYhMwt1Z\/+Ivfk5WrO5cYe\/abA3LuzadzuzMJAVPstpOqdSE7P1VwZA4BCECgiQiUOqXXVD\/PnUkhj2FEYIjANNFHDFWFAAQgkA6BSlNJeS7qRcCk428vcg3Z+V44ACMhAAEIeELATCX94lPvL7M4r0W9IY9hRGCIwHjy0YCZEIAABIpP4Lv3XCRXrrtPVn\/m0DJj8zipFwFT\/DaTmoUhOz81aGQMAQhAoEkJVDrgTpFkvR4m5DGMCAwRmCb9mKHaEIAABNIhUO6AO1Oaihh9sjipFwGTjo8bnuuxY8dkZGQkynfPnj3S19dXtYyTJ0\/K4OCgjI2NyZo1a5alD9n5VeGQAAIQgAAEnAlUuitJM8tyPUzIY1gwERhbiGgDqSRKTGtcWFiQ4eFhefXVV+XQoUMIGOduygsQgAAEIFCKwOknd8npp3aWvCvJFjFX9hyRi6\/tTg0iAiY1tI3LWKMvU1NTMjo6Ki0tLbJ3717p7OysGIVRx2o6fYjANM4X5AQBCEAAAiI6laRP286XSuLQqaT33n5FrrrjZVlxWUcqyBAwqWBtbKZGiAwNDUUZx7+Pl3bq1CnZtWuXbNy4MUqLgGmsP8gNAhCAQLMTMGfD6I4k3ZlU6kl7PQwCxoNWGI+4aERmZmZGjKCJV0F\/r8+tt96aaA1Mf3+\/dHV1Re90dHREXzwQgAAEIACBSgTeenyz6M6kUmfD6HvmkLtLrr9bLr9lvCEwZ2dnRb\/0mZubk4GBAZmcnJTu7vSmqhpiuGMmwayBcREwul7m8OHD8uCDD0ZOTrKI1+aqYka\/eCAAAQhAAALVCOg1Ay0fWVfybBh9t9GXPk5MTIh+2Q8CppqXcvy9yxSSpr3tttsiNZp0F9L4+Li0t7cTgcnRxxQNAQhAwEcC1c6G0To18pA7OwIzPT0diRkETIFbTnzKqNwiXl37smXLFjlx4sSy2pRycMjzhwV2J6ZBAAIQCIpAtQW9Wtk0DrkLeQwLZgqplm3U2mCSRmBCVK9BfTpQGQhAAAIFJlDpskdjdhqXPiJgCtwobNPKHWRnznvRHUfxRUwIGE+ci5kQgAAEPCdQ7YRerV6jD7lDwHjeaOoxP2Tn18OFdyEAAQhAwI1AtRN6TW66Hubcm8\/I5TeP133IXchjWDBTSG7NKHnqkJ2fnAIpIQABCECgEQTMCb0f\/uL3ZOXqzrJZNmo9TMhjGAKmSosM2fmN6IzkAQEIQAACbgSqbavW3Bq1HibkMQwBg4Bx63mkhgAEIACBuggk2VatBTTifBgETF2u8vvlkJ3vt2ewHgIQgIC\/BJJsq9ba1XtfUshjGBEYIjD+fgJgOQQgAAFPCSTZVm2qVs99SQgYTxtII8wO2fmN4EMeEIAABCBQG4Ek26rtqaSWm\/ql5abPORUW8hhGBIYIjFNnIDEEIAABCDSGQJLbqk1JtV41gIBpjK+8zCVk53vpEIyGAAQgEBCBpFEYrXItW6tDHsOIwBCBCeijgKpAAAIQ8IuASxSmlq3VCBi\/2kNDrQ3Z+Q0FRWYQgAAEIFATgbce3yy6tfoXn3q\/6vuuW6tDHsOIwBCBqdphSAABCEAAAukRMFcMrP7MIbly3X1VCzr7zQE59+bTctUdL8uKyzoqpkfAVMUZboKQnR+u16gZBCAAAb8IuERhtGZJt1aHPIYRgSEC41cvx1oIQAACARJwjcIkvbUaARNgY0lapZCdn5QB6SAAAQhAIH0CrlGYJFurQx7DiMAQgUm\/V1ICBCAAAQhUJeAahUkylYSAqYo93AQhOz9cr1EzCEAAAn4ScI3CmK3V5U7pDXkMIwJDBMbPXo7VEIAABAIkYKIwV9+zU67+5I5ENaw0lYSASYQwzEQhOz9Mj1ErCEAAAn4TcI3CVJpKCnkMIwJDBMbvno71EIAABAIjUMtamHJTSQiYwBqHS3VCdr4LB9JCAAIQgEB2BFzuSDJWlZpKCnkMIwJDBCa7HklJEIAABCCQiIC5I+m6nS9Jy0fWJXrHTCWdX5iNTunVBwGTGF2+CY8dOyYjIyOREXv27JG+vr6SBi0sLMjw8LAcP348+v3WrVtlaGioZNqQnZ+vtygdAhCAAAQqEdAojD5tO19KDCo+lRTyGBZMBObkyZMyODgoY2NjkaPN\/9esWbPM8Xv37o1+pqLl1KlTsmXLFtmwYUNJwROy8xP3CBJCAAIQgEDmBPSCR13Q6xqF0XuS9L6kK3uOyKt\/L7Jp0yaZnJyU7u7uzOuQZoHBCBiNvkxNTcno6Ki0tLSIipTOzs6yURgbqi1o4rARMGk2P\/KGAAQgAIFKBN749A3RFJJe9Ojy6F1JOpX0d1eMI2BcwOWRNi5CKokS2z4TgdFoTCl1ioDJw5uUCQEIQAACSsBEYT78xe\/JytWdiaGYqaQTp7tly65pIjCJyeWQMB5x0YjMzMxM2bUtaqK+c+DAAent7V2M3JSLwPT390tXV1f0646OjuiLBwIQgAAEIJA2ge\/ec5Fcue6+xFGY2dlZ0a\/V70\/L1W89Ilsfu0i2PcQUUtp+qjn\/WgSMKUzfnZ+fLyliTATGNkzFjH7xQAACEIAABNIm4Hqw3cTEhOiXPgc++75c1yry418+whqYtB1Va\/61TiFpefYC4PiiXyNgxsfHpb29nQhMrQ7iPQhAAAIQqJmARmF0HYxGYqo9JgKj6d44OS13fPCRaEHvxdeyiLcau1x+H58yclnEqyJF0x88eFBaW1uX2M8amFzcSaEQgAAEIGARcI3CmFd1DNv22U2y\/zGmkArboGrdRm3OhGlrayu5XgYBU1iXYxgEIACBpiFQ65bqkMewYLZRaysud5CdESkbN26M5gDjB9klWcQb4h76pun5VBQCEIBAAARq2VKNgAnA8bVWIWTn18qE9yAAAQhAIHsCtWypDnkMCyoCk0ZzCtn5afAiTwhAAAIQSI+Ay2JetSLkMQwBU6Wdhez89LoYOUMAAhCAQBoEdDGvXvSoB9sleUIewxAwCJgkfYA0EIAABCBQAAKui3kRMAVwWl4mhOz8vJhSLgQgAAEI1E7AZTFvyGMYERgiMLX3It6EAAQgAIHMCZx+cpecfmqn\/OJT71ctGwFTFVG4CUJ2frheo2YQgAAEwiaQdDFvyGMYERgiMGH3cmoHAQhAIEAC8ztvj2rVtvOlirVDwATo\/KRVCtn5SRmQDgIQgAAEikUg6WLekMcwIjBEYIrVK7EGAhCAAAQSEUiymBcBkwhlmIlCdn6YHqNWEIAABJqDQJILHkMew4jAEIFpjp5OLSEAAQgERiDJNBICJjCnu1QnZOe7cCAtBCAAAQgUj0C1aaSQxzAiMERgitcjsQgCEIAABBIRqDaNhIBJhDHMRCE7P0yPUSsIQAACzUNA70X6\/s7bZfVnDsmV6+5bVvGQxzAiMERgmqenU1MIQAACARKoNI2EgAnQ4UmrFLLzkzIgHQQgAAEIFJdApWmkkMcwIjBEYIrbK7EMAhCAAASqEtBpJBUxV39yx7JpJARMVXzhJgjZ+eF6jZpBAAIQaC4C5aaRQh7DiMAQgWmuXk5tIQABCARIoNw0EgImQGcnrVLIzk\/KgHQQgAAEIFBsAuUOtQt5DCMCQwSm2L0S6yAAAQhAIBGBUtNICJhE6PJPdOzYMRkZGYkM2bNnj\/T19ZU06tSpU7JlyxY5ceJE9Pve3l4ZHR2VlpaWZelDdn7+HsMCCEAAAhBoFIFS00ghj2HBRGBOnjwpg4ODMjY2FrUF8\/81a9YsaRsLCwsyPDwsPT09kcAx37e1tcnQ0BACplE9iXwgAAEIQCBTAmYa6cNf\/J6sXN0ZlY2AydQFtRWm0ZepqanFSMrevXuls7OzbBTGLiX+rv27kJ1fG2neggAEIACBohL47j0XLTmVN+QxLJgIjAoWfUwUJf59pcaGgClqV8QuCEAAAhBwIaDTSPro1QJEYFzI5Zg2HnFRUTIzM1NyWsg206yH2bBhQ8lojVGv\/f390tXVFb3a0dERffFAAAIQgAAEikTATCPN\/95RWbn6X8nc3JwMDAzI5OSkdHd3F8nUum0JKgJjTxklETBm\/YtSrLaI1yatYka\/eCAAAQhAAAJFIvDeWzOiu5HGZlbLC29fsWgaAqZIXorZ4jqFlES82OG38fFxaW9vJwJT4DaAaRCAAAQgIJGAefuDvyTnfvu\/yvT0tExMTBCBKXLDiEdcKi3irbbzyK5nyAugiuxPbIMABCAAgdoI2NupQx7DgplCSrqNWpuDipv5+fmy00YImNo6DW9BAAIQgED+BOzt1F\/\/zj\/Kpk2biMDk75bKFpQ7yM5EXDZu3Cg33njjkkPsTI5r166VgwcPSmtr65JCQlavRfcn9kEAAhCAQG0EzHbq1y79ZQRMbQj9fwsB478PqQEEIACBZiNgrhX4zse2ImCazfmmvgiYZvU89YYABCDgLwGzDuYH2\/4XAsZfN9ZnOQKmPn68DQEIQAAC2ROwz4P5vf4R1sBk74L8S0TA5O8DLIAABCAAAXcCug7mux\/7ffnP\/\/15BIw7Pv\/fQMD470NqAAEIQKAZCeg6mBM\/vlT6v3YOAdOMDQAB04xep84QgAAE\/Cdg1sHc8eoNCBj\/3eleAwSMOzPegAAEIACB\/AmYdTD\/4f9cL\/\/tS09yF1L+LsnWAgRMtrwpDQIQgAAEGkNg4bWvyPd33h7di3Tfo88hYBqD1Z9cEDD++ApLIQABCEBgKQFdyPv821fILbtfRMA0W+PwScDMzs7KM888I3fddZd0dHQU3lXYm66L4AtfQ4C20LxtQdfBPP3M0wiYdJtAMXP3ScD4ZKt6G3vTbfPwha8hQFto3rZw+sldcvqpnaIH2nV3d6cLIuPcg7nMMS1uPnV8n2xFwKTVYv8lX9pDuox94uuTrXw2NLbdmoW8K7a\/LJ23\/PvGZp5zbgiYKg4wHb+\/v1+6urpydlfl4ufm5mRgYEB8sFVrgr3pNif4wtcQoC00b1tQ399y9G5Z\/ZlDcuW6+9IFkXHuCJgqwHXuWEXB9PR0xq6hOAhAAAIQgED9BP7so2\/Kml\/7j5GICelBwCTwpooY\/eKBAAQgAAEI+EbgFy55L7jpI\/UBAsa3loi9EIAABCAAAQggYGgDEIAABCAAAQj4R4AIjH8+w2IIQAACEIBA0xNAwDR9EwAABCAAAQhAwD8CCBj\/fIbFEIAABCAAgaYngIBp+iYAAAhAAAIQgIB\/BBAwFXy2d+9eOXDgQJRicnKyMMcwnzx5UjZv3izz8\/PS29sro6Oj0tLSUrH16YF8R44cSZS20c3YxV6beVtbmxw6dEjWrFnTaJMq5udi77Fjx2RkZCTKb+3atXLw4EFpbW3NzF4XW41RCwsLMjw8LD09PdLX15eZrVqQi702Wx\/4njp1SrZs2SInTpzI7TMjKd84W9MI9uzZk2mbSGpvvO348NlgDkHNq+1m2rFzKgwBUwa8Nj4dTHVA+va3v734\/ywHp1Km2YPP+vXrEw1EpiMlFTuNbIsu9uqH6tTU1KLI0u+PHj2aqShwsdduI9outL2oqEwiKBvB2MVWuzwzeGU9WLnaqzw7OzszHVBtTi72mrQ6sA4NDUVCbXBwUMbGxjIT4C72xttfvC03on1Wy8PFXiMOla3e51P0zwYjzPbt2xfZm+cfkNX84PPvETBlvKcfnvpohzEdbePGjblHYeIfjNU6htbj+PHjsm7dOnnnnXcyG1wNVld7bXfkMQjUY2\/Wg0AttupAsG3bNjlz5oxs2LAhU3HgYm8R+pyLvZp29+7dsn\/\/\/kwjcJX6S7XPBvNuXBxkNaC58rUFYdE\/G+J\/jGl7fvjhh+Xee+\/NTNBm5cc8y0HAlKAfD7HnGXKv9pdStUFTf2\/+YrGjG1k1urh91ezNW8DUY68terPgW4utauOv\/MqvyF\/8xV9kPoXkYm9eg6rtNxd74wNWFv6v97PBvJ+X7S58S0Vgsv48c7G3lIDRadsi\/BGcR9tMq0wETAUBYze2vMPZxsz4X1VJ\/\/LL80PKXnuT1F6tb9ZTMlpmLXzNFF3W8\/Kutir7w4cPyx\/90R\/Jrl27chEwSduCvTbCtP2s16G58NX+NTMzE5ma17o5F3vzjr7U0tfMH5IaUd66dWsUHc\/yceEbv\/3bfJ\/1tG2WfPIoCwGDgEm13bl0etsQHRAeffTRzBfx1mqv+UB+4IEHMrPZxVY7hN3R0ZFo7VSjG4aLvZrWZhn\/vtG2lcrPxV6zrsiIrKLba\/9BZNb6Zb2+z4VvqTUlWdvtYq\/ytRdKb9q0KZrCz2PhfBZ9Ja8yEDAVBIxpbD5PIfkUJrZtzUO8GBFifzC6THll3U5cQtqa9qtf\/eqSNV1Zf5i62BvvllmzdW0L5aYMsmRcC18TOco6mtEsfIsQ6cpLXGRRLgKmDGV7yqgICwqNmfEpmPhfBeUaTV5TSK725rG7wGbmaq\/9btbrNlxstben2zZnGYp3sbecgMlyDYGLvfF+mMdnhou9yjcPUVhrXyuCQHTlG69r1rvSshAQeZeBgCnjAfuvGd+3UZtwZtaL3uIfktW2fecRdq\/0l34Se+1oTdbiy2Ubql3PvAYuF3t921IfF68ukbtGDQIufLVMsyNt+\/btueyMcbG31BRSltO1rp9lttjRM7p0Aa\/ZYt8of5OPcBt1pUbg40F25ULCeUVglG+lw6pse8tFCbJevJnUXiMMzUF2WS\/idWFbBAHjaq+9hiAPtq722gfZ+WBvHluR45+3Ln3NLITVPHzga7ffPM7gagaBQwSmGbxMHSEAAQhAAAKBEUDABOZQqgMBCEAAAhBoBgIImGbwMnWEAAQgAAEIBEYAAROYQ6kOBCAAAQhAoBkIIGCawcvUEQIQgAAEIBAYAQRMYA6lOhCAAAQgAIFmIICAaQYvU0cIQAACEIBAYAQQMIE5lOpAwEcCeoaKXjC5Y8cOyfpOHh95YTMEIMBBdrQBCECgAATsu5oKYA4mQAACHhAgAuOBkzARAqET0FOYb7vtNunu7g69qtQPAhBoEAEETINAkg0EILCcQLl7l+y7gvSumIcffljuvffe6E4e887x48ejDPM6Nh5\/QgACxSaAgCm2f7AOAt4TKHXJpUZc9BkaGoruyjp8+LA8+OCD0c\/04jt9RkdHRcVN1pdkeg+cCkCgSQggYJrE0VQTAnkRiN\/UHP9eBYo+fX19kZgZHByUsbGxXG5IzosR5UIAAu4EEDDuzHgDAhBwJGBHXOzpI91xZK9\/if\/OsRiSQwACTUQAAdNEzqaqEMiLgC1MnnjiCens7IwiLvHt0wiYvDxEuRDwjwACxj+fYTEEvCNgpo3Wr18vzz333OIUUXz7NFNI3rkWgyGQGwEETG7oKRgCzUVA17qMjIxIb2\/vkgW6SkGjMfqYHUi680gX+OqDqGmudkJtIZCUAAImKSnSQQACdRFQIbJ582a5\/\/77I8GiYsXePm0yZxt1XZh5GQJNQwAB0zSupqIQgAAEIACBcAggYMLxJTWBAAQgAAEINA0BBEzTuJqKQgACEIAABMIh8P8BRix3KR\/kvQoAAAAASUVORK5CYII=","height":0,"width":0}}
%---
