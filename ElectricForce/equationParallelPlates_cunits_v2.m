%[text] # Two Parallel Charged Plates and a Point Particle
% The variables
syms a d v r real positive;
assume(d>0)
assume(r>=0)
assume(v>=0)
assumeAlso (v < 1)
assume(a>=0)
taylorOrder = 12;

vp = [0,v,0];
%%
%[text] ## The vectors
%[text] d= distance from the charge to the plate
%[text] r=radius of a circle in the plate
%[text] The speed of light is 1 (c=1)
R = (d^2+r^2)^(1/2) % total distance from a point in the plate to the charge %[output:0a21a822]

% Get two points per circle. 
vec_l1 = [d,r*cos(a),r*sin(a)]/R  % velocity vector of light from a point 1 in left positive plate  %[output:4613f0bd]
vec_l2 = [d,-r*cos(a),-r*sin(a)]/R  % velocity vector of light from a point 2 in left positive plate  %[output:161358a7]
vec_r1 = [-d,r*cos(a),r*sin(a)]/R  % velocity vector of light vector from a point 1 right negative plate %[output:554fe653]
vec_r2 = [-d,-r*cos(a),-r*sin(a)]/R  % velocity vector of light vector from a point 2 right negative plate %[output:86e39d00]

r_vec_l1 = (vec_l1 - vp) % left plate relative velocity of particle (v in light speed fraction =v/c) %[output:9f9767cb]
r_vec_r1 = (vec_r1 - vp) % right plate relative velocity of particle %[output:34f9ec89]
r_vec_l2 = (vec_l2 - vp) % left plate relative velocity of particle (v in light speed fraction =v/c) %[output:6900aafa]
r_vec_r2 = (vec_r2 - vp) % right plate relative velocity of particle %[output:123f7d6a]
%%
%[text] ## The Cosine 

%dprd_l1 = simplify(simplifyFraction(dot(vec_l1,r_vec_l1)^2/dot(r_vec_l1,r_vec_l1))) % Cosine^2 left plate
%dprd_l2 = simplify(simplifyFraction(dot(vec_l2,r_vec_l2)^2/dot(r_vec_l2,r_vec_l2))) % Cosine^2 left plate
%dprd_r1 = simplify(simplifyFraction(dot(vec_r1,r_vec_r1)^2/dot(r_vec_r1,r_vec_r1))) % Cosine^2 right plate
%dprd_r2 = simplify(simplifyFraction(dot(vec_r2,r_vec_r2)^2/dot(r_vec_r2,r_vec_r2))) % Cosine^2 right plate

dprd_l1 = simplify(simplifyFraction(dot(vec_l1,r_vec_l1)/dot(r_vec_l1,r_vec_l1))) % Cosine left plate %[output:05796605]
dprd_l2 = simplify(simplifyFraction(dot(vec_l2,r_vec_l2)/dot(r_vec_l2,r_vec_l2))) % Cosine left plate %[output:34d747be]
dprd_r1 = simplify(simplifyFraction(dot(vec_r1,r_vec_r1)/dot(r_vec_r1,r_vec_r1))) % Cosine right plate %[output:90cd3a58]
dprd_r2 = simplify(simplifyFraction(dot(vec_r2,r_vec_r2)/dot(r_vec_r2,r_vec_r2))) % Cosine right plate %[output:61660ac0]


%%
%[text] ## The Net Force of the Four Points

dforce = ((dprd_l1/R^2)*vec_l1 - (dprd_r1/R^2)*vec_r1); % net Force in the x direction first points.
dforce = dforce + ((dprd_l2/R^2)*vec_l2 - (dprd_r2/R^2)*vec_r2); % net Force in the x direction second points.
dforce = simplify(expand(dforce),'Steps',100) % Two charges %[output:6d204237]
dforce_x = dforce(1) %[output:9cd5c9e8]
%%
%[text] ## Check at v = 0, r =0
df_v_zero = simplify(subs(dforce,[v,r],[0,0]),'Steps',100) % Due to the four point charges %[output:769f5ec5]
df_v_one = simplify(subs(dforce,[v,r],[0.01,0]),'Steps',100) % Due to the four point charges %[output:45e059fe]
%%
%[text] ## Integrate all the angles
dforce_xT = expand(taylor(dforce_x,v,Order=taylorOrder));
dcircle = int(dforce_xT,a,0,pi); % From zero to pi i.e., the full circle
%%
dcircle = simplify(dcircle,steps=100); % Net force of the two full circles with charge density
df_v_zero = simplify(subs(dcircle,[v,r],[0,0]),'Steps',100) % Net force of the full zero radius circle (point) %[output:926989fe]
%[text] ## 
%%
%[text] ## Integrate for all values of r

dtotf = simplify(expand(r*dcircle),steps=100)/pi; % Add all the circles 
%dtotf = taylor(dtotf,v,Order=taylorOrder);
totf = int(dtotf,r,0,r);
totf = simplify(totf,steps=100); % The force on a charge due a charge surface density

%%
%[text] ## Set the Plates Radius to Infinity

totfInf = simplify(limit(totf,r,Inf)) % Net force by plate charge density %[output:733de68c]
totfInfV0 = subs(totfInf,v=0) %[output:84edc3b9]

totfInfRel = simplify(totfInf/totfInfV0) %[output:3830c242]




%[text] ## 
%%
%[text] ## Plots
figure(1)
plotrange = [0,0.99];
fplot(totfInfRel,plotrange) %[output:649ddfc4]
hold on %[output:649ddfc4]
relEq = sqrt(1-v^2) %[output:988a2589]
relEqApx = taylor(relEq,v,Order=taylorOrder) %[output:4a5e14ba]
fplot(relEq,plotrange) %[output:649ddfc4]
fplot(relEqApx,plotrange) %[output:649ddfc4]
xlabel('v/c')  %[output:649ddfc4]
ylabel('ratio')  %[output:649ddfc4]
legend({'Force Carriers',"Lorentz Mass","Lorentz Apx"},'Location','northeast') %[output:649ddfc4]
hold off %[output:649ddfc4]

%%
%[text] ## 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":22}
%---
%[output:0a21a822]
%   data: {"dataType":"symbolic","outputData":{"name":"R","value":"\\sqrt{r^2 +d^2 }"}}
%---
%[output:4613f0bd]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_l1","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:161358a7]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_l2","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:554fe653]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_r1","value":"\\left(\\begin{array}{ccc}\n-\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:86e39d00]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_r2","value":"\\left(\\begin{array}{ccc}\n-\\frac{d}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:9f9767cb]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l1","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }}-v & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:34f9ec89]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_r1","value":"\\left(\\begin{array}{ccc}\n-\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }}-v & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:6900aafa]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l2","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & -v-\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:123f7d6a]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_r2","value":"\\left(\\begin{array}{ccc}\n-\\frac{d}{\\sqrt{r^2 +d^2 }} & -v-\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:05796605]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_l1","value":"\\frac{d^2 +r^2 -r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }}{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 -2\\,r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:34d747be]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_l2","value":"\\frac{d^2 +r^2 +r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }}{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 +2\\,r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:90cd3a58]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_r1","value":"\\frac{d^2 +r^2 -r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }}{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 -2\\,r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:61660ac0]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_r2","value":"\\frac{d^2 +r^2 +r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }}{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 +2\\,r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:6d204237]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce","value":"\\left(\\begin{array}{ccc}\n\\frac{4\\,d\\,r^2 +4\\,d^3 -4\\,\\cos \\left(2\\,a\\right)\\,d\\,r^2 \\,v^2 +4\\,d^3 \\,v^2 }{{\\left({\\left(1-2\\,\\cos \\left(2\\,a\\right)\\,v^2 +v^4 \\right)}\\,r^2 +{\\left(1+2\\,v^2 +v^4 \\right)}\\,d^2 \\right)}\\,{{\\left(r^2 +d^2 \\right)}}^{3\/2} } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:9cd5c9e8]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce_x","value":"\\frac{4\\,d\\,r^2 +4\\,d^3 -4\\,\\cos \\left(2\\,a\\right)\\,d\\,r^2 \\,v^2 +4\\,d^3 \\,v^2 }{{\\left({\\left(1-2\\,\\cos \\left(2\\,a\\right)\\,v^2 +v^4 \\right)}\\,r^2 +{\\left(1+2\\,v^2 +v^4 \\right)}\\,d^2 \\right)}\\,{{\\left(r^2 +d^2 \\right)}}^{3\/2} }"}}
%---
%[output:769f5ec5]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\left(\\begin{array}{ccc}\n\\frac{4}{d^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:45e059fe]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_one","value":"\\left(\\begin{array}{ccc}\n\\frac{40000}{10001\\,d^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:926989fe]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\frac{4\\,\\pi }{d^2 }"}}
%---
%[output:733de68c]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInf","value":"4-\\frac{4\\,v^2 }{3}-\\frac{4\\,v^4 }{15}-\\frac{4\\,v^6 }{35}-\\frac{4\\,v^8 }{63}-\\frac{4\\,v^{10} }{99}"}}
%---
%[output:84edc3b9]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInfV0","value":"4"}}
%---
%[output:3830c242]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInfRel","value":"1-\\frac{v^2 }{3}-\\frac{v^4 }{15}-\\frac{v^6 }{35}-\\frac{v^8 }{63}-\\frac{v^{10} }{99}"}}
%---
%[output:988a2589]
%   data: {"dataType":"symbolic","outputData":{"name":"relEq","value":"\\sqrt{1-v^2 }"}}
%---
%[output:4a5e14ba]
%   data: {"dataType":"symbolic","outputData":{"name":"relEqApx","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}-\\frac{v^6 }{16}-\\frac{5\\,v^8 }{128}-\\frac{7\\,v^{10} }{256}"}}
%---
%[output:649ddfc4]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQ1wXtV555+4JqAEXCOgAUuwosHQTCbDR5aRKsjg7JaSdsfutkBkebclXm\/WuyFEaUBfsAGbDbY+8DTaQFov63U8zVo2X5PgzDaQziQpQUVDCXjabEIcghYkwQaQTagjWCjsPNcccXR13\/e99\/24955zf3fmHVt67z33Ob\/n3Hv+es5zznnP22+\/\/bZwQAACEIAABCAAAYcIvAcB45C3MBUCEIAABCAAgYAAAoaGAAEIQAACEICAcwQQMM65DIMhAAEIQAACEEDA0AYgAAEIQAACEHCOAALGOZdhMAQgAAEIQAACCBjaAAQgAAEIQAACzhFAwDjnMgyGAAQgAAEIQAABQxuAAAQgAAEIQMA5AggY51yGwRCAAAQgAAEIIGBoAxCAAAQgAAEIOEcAAeOcyzAYAhCAAAQgAAEEDG0AAhCAAAQgAAHnCCBgnHMZBkMAAhCAAAQggIChDUAAAhCAAAQg4BwBBIxzLsNgCEAAAhCAAAQQMLQBCEAAAhCAAAScI4CAcc5lGAwBCEAAAhCAAAKGNgABCEAAAhCAgHMEEDDOuQyDIQABCEAAAhAovIAZHh6WtrY26erqojVAAAIQgAAEIOAIgUILGBUvO3fulO3btyNgHGmwmAkBCEAAAhBQAoUUMHNzc7Jp0yZZuXJl0Ao+8YlPIGB4HiAAAQhAAAIOESisgDl8+LCsWrVKBgYGpLOzEwHjUKPFVAhAAAIQgEAhBYxx+\/z8fCwB88tHuoNLfu19rSVbzPRrp0jrCS8H3+v\/zXFW8wkL\/59+7dSF\/8+8dqqEvzvr5BPEPkdPXvarl4JrzL80WQhAAAIQgEASAq2traIf3w4ETIUIzPT0tPz47o\/F9vsZzUtPXRXxuzgFTh45T1reEUUz74giI3Cen3u3hNk5ERVEepj\/24Jn2a+OCat3xdDin+PYwjkQgAAEIOAmgfb2dhkdHfVOxCBgKgiYRx99VDZs2BA4v6WlpS6tNyxobNHz1vx0cI+3fnXs32P\/nwn+\/ad3fhd1TpRhRuyo+NH\/G5Gj56o40kN\/b4shFTsmMnTmyceiR\/qz+V34u1qATE5OytjYWF3Z1mJPpWuxtxKh2r6Hb238yl0N28axDd6nOX6XGdv27t0rHR0djQWRcukImJgCJs\/ON2Jn4rv3BYJgaGuPrGp+e4n4UeFjCyO7rZUSO8eEz1IBZISN\/hstdJoWfn\/pOceSpcOHEYd5ZmvbjL2NfTvBt3F8Yds4tlpynvnm2bZavYKA8UDAmEaQpKEaIaNRHVvYaLRHf\/fmy49Gti0jdDSCYyI6zy+7QJ6dm5flp3TIc4dfk2fnXou81o7e6P812nPffffK1b+3JvjYkZ5aG3Yjrk\/CtxH3T1om9iYllux8l\/i6ZGveBYFrf4y55vskTyECxiMBo\/k69913n1x55ZV1G+s0QueNlx9diN4YkRMV0Vn2TqKzipvlp7QHYkYFjv687H0twc9G4Dzy9JFYQscIn0s+eHJwfqmITpKGX825jeBbjR1xr8HeuKSqO88lvi7Zqt7A3urapGviqtZaFlrAxIHns3qNU\/8456jIsSM55aI4KnCWNbXKcae2B0Vr1OaYuGkNhM2xCM58RaFjR3OMoDmr+diwVd4jOXGYcg4EIACBehDwuQ9DwFRoIT47vx4PR5wyjMAxw1JvvDQZOUQVJW6OO\/XdpDMTvTEiR+\/9g58dKTlspULGFjRG4GQVwYnDinMgAAEI1JOAz30YAgYBU89nJXFZccRNJWFjbmoPTWkUxwgc\/TdquOqSDx5LLjaCRoeoiN4kdiEXQAACOSaAgMmxcxptms\/ObzS7WspPKmzeFTmlpwmWEjhR4sZEb2xxQ+SmFo9yLQQgkAUBn\/swIjBEYLJ4pqq+p0kqfv25e4MyooajVMwcf+aVwfeaY2MPQ5W7sRmieuTpwwvRm7C4eTeheOU70ZpjeTeIm6pdyoUQgEADCSBgGgg370X77Py8s09inx2x0SRiI3BMGSpqjjulI5gJlUTUmOtV3KiYsROMo8SNLWZ0SAphk8SLnJuUgM7W0Q9H8QjE3R7A5z6MCAwRGG+f\/ErRmlpFjYKLG7VB2HjbzDKrmAqX3t7eYBVYjuIRiLs9AAKmeG1jocY+O7+oblVhU2oIyh5+Ov7Mq4Lp3dUcJt9m\/LHnyw5HGWHDDKlqKBf7mkZsc1Jsou7UPsn2AD73YURgiMC489Q2yFI7UhPOqbFnQFUz9BQ2Oa6w0RlSmm\/DMFSDnO5BsT53TB64p6FVSOL7JOc21OgGFI6AQcA0oFm5X6SdUxMlajRJ+Fi05qq6VNYeinr3\/++uVGzPikLU1AW584X43DE575wGVyCJ75Oc22Cz6148AgYBU\/dG5WuBZugpStCYBOFahp2iuKmYKTUMZYsaHYLqvvh0X9FTrwgCrndMc3NzsmnTJjl48OCi2q1du1aGhoakqakpVb8bnuam27dvl66urppsmJ+fl4GBAenu7q7rTtBJfJ\/k3Joqm8HFCBgETAbNzo9b2rk0rz9338JeUXZycNN5n697Ze0ZUboSsT0bikhN3XHntkDXOyYjYPr7++vauVfjMGV5ww03yO7du2X16tVihEdnZ2fNIqYaeypdk8T3Sc6tdN+8fY+AQcDkrU06a48KGrPpZVjQmHVpGiFoFBiixtlmU7XhrndMlQSMEREHDhwIGO3duzcQOocOHZK+vr4Fbrt27Qr+b6I5559\/vujvmpubg3M3btwos7OzYv\/ehl5KrOi127Ztkx07dgRl7d+\/XwYHB4NLV61aFYgdncqsEZYjR47Iww8\/HPzu\/vvvX\/KzicCUskd9OTw8LC+++KJ89KMflZtvvlluvfVWMXUPR4OS+D7JuVU3xowuRMAgYDJqev7fNktBY0RNueEnM+RETo2bbdH1jqmSgNEOXYWHDifpMJOJkKi3VJTcfvvtgaAJD9Oo0JiampJPf\/rTgagxER67PHt4ygiikZGRIPoSdRiBYYSRlqXH5z73uUDAqKDR+xhbwj+rgDn33HNL2mPXT20wddAyldPWrVvllltuCYSUHkl8n+Rc154EBAwCxrU266y9tqCZf2psoR5m6nY9ZjlVgmPn1Aw\/OLVwOkNPlcjl7\/twx2SSv\/Nn6TGLwvuMReXAmCiJCgw7d8SOklx00UVBBMYIjnCkxNQ\/LDpKnVfq9+U4GoFhBIwZagpHc2xxpeWp8LGjQybC89Of\/nTRd7aAKSWoNmzYsBCVKmcrAiavT0QKdvns\/BTwcYsyBEolBR8TMu3BDKdq16GJC96e\/RSVT2OiNP1XnB23SM5LiUD43TT84DNii9KUzIh9m\/4r2sRuR+UiMFHfaeff1tYmKmDsoZ2wULEFjHby9mGGfuxIS5wITHg4S8vcvHnzQgTGDBGFo0FhAVPKnpdfflnGx8cXJS9rfXfu3BmYb4bPwnUL\/75WsRPbmTk5kQhMBUcgYHLSUgtgRpSgSTM6YxCXi9KwPk1+GqIvEZioJN5SQkAjHWEBUy4CExYFUd4rlQNjD9185zvfkYmJiQWBEY7AxBUwpexRX5b6Tu24\/vrr5cYbb1wY4krSLyU5Nz+tO54lCBgETLyWwlmpEig13GSiM41KBo6qZKkEYYadUm0SS27mesdUSw6MHYEJl6PiQsWGDjN99rOfXciB0d\/v27dvYQjHBlpqFpLJZTFlaj6OCh7NrdH8GzOEFEfAhHNgbHt0CMkWMCbapNO4yYEp\/ZwhYBAw2b6FuXssAnmJzqixpYadEDSxXFm3k3wXMOVmIdkCRoGWmt1j\/z5q+CgsYuwhHh0i0uiQHna+jubpaILwgw8+uDBbKI6AMTOozKwo255wBCacH8QQUvRjg4BBwNTthUpB6RCw158xycBmqKmeqwPHrU2cPBpmOsWlGf881wVM\/JpyZphAEt8nOdc10ggYBIxrbRZ7QwSMoLFnNmUx1GSbZfJo7MRgO0LTffEZwawUjuoJ+NwxVU+lGFcm8X2Sc12jh4BBwLjWZrG3DAEjZqIW0kszbyYsZo7l0RyWsKAxScHMckrerH3umJLTKNYVSXyf5FzXKHolYOyVEivtY2FPUbPHOmsJ1bnmfOz1m0C5vJmsxIwSt4ecxh97IfhZD43I6LRthpvitUufO6Z4BIp7VhLfJznXNaLeCBh7Lr86wV7oKOyUcDa6nfGNgHGtCWNvHAJ5FTNG0JgVg8OL6yFoSnvX544pTpsu8jlJfJ\/kXNeYeiNg7GluuopjKVESNedfxc+ePXvkpptuWrIDqs\/Od62xYm99CORZzNiCJjzcxKJ6i\/3v+rup0jTq+rT26FKSrr5r+o3HH398YcNHU7KJ\/FeK+tezPkl8n+TcetqYRlneCBizN4WZ9hb+2cAsJWDC0\/LM+cb5PT090t7eHvxaN\/DSDwcEXCfggpgx+TMMNyFg6vW8VStgdNPG3\/\/931\/YoVr7k9tuu01+9KMfyfr161PbubqSKJmenhb96DEzMyO9vb2xth2oF9+0yvFKwOgy07rwjx7l9pJQcaMNwN6YS3f9NFup2\/BNQ7F\/p2JGPxwQ8IlAlJjJejZTmG+p2U0anTmruSnIoSnKUakTyzuHShEY+927du3aYBVcPezdn3V9FF2XRX8XtWu1\/mG6YsWKRd+ZBeV0A0W9Vsv+0pe+tAhXOJpi\/vC94IIL5Omnn16I1pvovV78kY98JOh\/wmu42DmWpXIv4+Zkhv+wLrWVwNjYmOjHPuJsO5D3NhO2r5ACJrxA0he\/+EV54IEHFjYHixIwo6Oj0tLSEnxFBMa1Zo69SQmEZzNlsaVBJZvtZOCo3BnfZzb5LGDMAnS647QRKLrwW3j3Z20jdrTdXlFXv9NF46677rpAWNi7UWt0IirqHrVsv5Zj+owrrrhC7rnnnoVl\/fUPZT1092v9A3rdunWBmDKbO9q5mfZ+R3qNRm6uueYaKfX7Ujtj67WVfG9HYCYnJwMxg4Cp9EbJ8Pu4Q0hRJmpjsHcJjRIwPjo\/Q3dxa4cIhNeZMWImjc0mk2AqF53xcd2ZcCf25otT8sYv3t1hPAm7NM497jfaZPlpbQu3KheBCa9Ma97Rd9xxR\/CHphEI4TLsPZROOeWURZM57DKjBEw5e+xyv\/\/978tll10WCCsjQu6\/\/\/5AwJgRAFNJWxBFbdhoxEicPZuq7ZcqiZ00fN+oe3gTgQkPGZWbWRSGGU4ArrahNMpJlAuBvBAoJWaynJYdxcaImWP\/vhCc4ts07XDHdPjurXL4ni15aSpL7Dj56i1y8idviSVgwu9kk7OiQz0qYMzS\/eHhGlO4DgFF7VpthEJYwJTa0NGUF95VWkXMH\/3RHy1Ece66664FARNOO7C3DCi11EeSJUDiRGCK0od5I2CSTqPWkJ8m\/NqhSt2rInz4rF5z+6bDMCcIzD\/1ZXnjpUl58+VHA3ubzuuRLLYyiAPr2Iymw2KGmoyY0WtdHWoqcgTGFjDhnZpNewgn6paKwOisVXvYJ6o92QJGc2i2bt0a7Ax99OjRoB8xfzBffvnlwUaPZoftSkNSph5RQimqPzLnJemXkpwb51nK0zneCBiFWkrFltqa3SR9lZv+5rPz89QQscVdAqXyZfI2xGQIm+iMPatJVwS+9JyVTokZ199N5YZsKuXA2B2\/nT5gX6dDSHaeSykBo9GTQMi+s3FjJQGjwsIk3ZrUglICRvukr3zlK8EEkR\/+8IdBrozex8xe0hyYUr+vJQeGCIy779O6Wu76S6KuMCgMAhUIuDLEVE7MmOhM3iMzrr+bSg3\/GFFQbhaSLWDCkzLMH6TlIjB6jUZK3nrrLTnjjDPkoYceWtSyw6uzh\/8ILrcYqv2HtE4QefLJJ4Mhr\/BsKWNnKfvLPWpJfJ\/kXNdecF5FYBoB32fnN4IXZUJACbgWlVGboyIzeRYzvJuK+6wl8X2Sc10jioCp4DGfne9aY8VeNwm4FpWxxUzUasB5iczwbnLzeaiH1Ul8n+TcetiWZhkIGARMmu2NexWcgCb+mp2y8zodO+yivObM+NwxFfwxqVj9JL5Pcm7FG+fsBAQMAiZnTRJzikAgKiqjs5g08TfPR5SY6b+iLZNVgH3umPLcBvJgWxLfJzk3D3VLYgMCBgGTpL1wLgTqSiBYfG3uG0uiMnlbVyaq0mExk\/Y6Mz53THVtZB4WlsT3Sc51DRUCBgHjWpvFXk8JaFRm\/qkxef25e4MamoiMDjXl\/TBiJrzOTCNXAPa5Y8q7v7O2L4nvk5ybdb2S3h8Bg4BJ2mY4HwINJeDq8JKBUkrM1Dv51\/WOqdJmjo1sZEl3o7ZtUe433HBD5Oa\/jbQ5bMOGDRti7W\/kejspxxQBg4BJ65njPhBITCAq6deF4SVT0eEHn3lnevbi7QzqIWZc75hcFDBmATr1r9l9OnGjrsMFSXyf5Nw6mJZqEQgYBEyqDY6bQaAaAnZUxsxecknIlMqXqWWIyfWOqZKAKbeQ3ZEjR+Thhx8OIhDhBeLMQngmyrJixQoxq67rd7oVgC5id\/DgweDatWvXiu6xZB+lVmfXMvfs2SNXX321fPnLX5YdO3ZIc3Pzwm7VF1xwgezatUtmZ2fFlGHvyxdeEK+aZ0GvSeL7JOdWa09W1yFgEDBZtT3uC4HEBMLDSy7lyZjK1muIyfWOqZatBHSDRLP0v72VgD28o7w3btwo1113XbBLtJ6nwmJoaEiidqPW80vtXWR8p6vs6mHK012pdWsBI0yeffbZQMDoztN9fX3BxpOtra0Lu1br73UjyHLbFsR5KJL4Psm5ce6dp3MQMAiYPLVHbIFALAKu58mYSupCeeOPPb9ox2wdXuq++PSKHMIdkzL5p19NV7wuqxN+7X2twWaf5ignYLRuZudo3WxRf1YBcscddwSioLOzMxAR4TLsCIfuhWREhO4rZJcZJWAqRYTs\/YtMeUaMRO1mbUdeNHKjUZ5ly5YtRG1q8UMSUZLk3FpsyuJaBAwCJot2xz0hUDcCrufJKIhqhpjCHZNy0FlceT00WmYP+5UTDBrpmJiYCKIlKmDMcJCKABUw9m7UZjjIrrcO31x00UWxNnO0h4CMMIpiaDaK1CiOOXQISiMuZkfr8CaTbW1tgdDSwxY0tfooiShJcm6tdqV9PQIGAZN2m+N+EGgIAXsatot5MgZK1BBTVFSmyBEYW8Bcf\/31cuONN0p49+ZymznaERgjPsqJFyNA9F97+MeIknXr1snAwMBCZCgckVFffetb3wpcrDtQl9tpOs7DkUSUJDk3zr3zdA4CBgGTp\/aILRComYBvQmb8sReCCI1ZKM8k\/rreMdWSAxOOdBhhYaIkt99+u+gQ0rZt2xaGbEoNId11111BmyuXl1LKVhMpuvnmm+XWW28NyjE5Nmb4Su3YunWr3HLLLUFujCYB33TTTUHUptojie+TnFutPVldh4BBwGTV9rgvBBpKwPWZSzacqKjMH5\/9itz5xc\/EWgukoaCrLNyIAp0NZB9mFlG5WUi2gDHRDjPTyMz+KReB0Wt06Omtt96SM844Qx566KFFNmzevHmRoDE5ODpcpENO5jCC6bbbbpP7779fNLl4586dwdf2DCnbXhU9U1NTNSXyJhElSc6t0pWZXYaAQcBk1vi4MQTSIOCjkNGozOxPD8qJPxhxVsCk4fu07lGv6dFx7U0iSpKcG\/f+eTkPAYOAyUtbxA4INJRA1BRsl9aSCcO556++J\/3X\/jsETENbTbzCETDxONX7LAQMAqbebYryIJBrAr4IGZ\/\/ss51A8qBcUl8n+TcHFQtkQkIGARMogbDyRDwhYDrQsZ0TD09PdLe3u6LW6hHDAIzMzPS29sbK\/qGgIkB1NdTfHa+rz6jXhBIQsBVIaNTgbUTm5ycTFJdzvWEgIrW0dHRYKXfcofPfRgRGCIwnjzOVAMCtRFwMdlXRYx+yh2PPK2r\/R6biq3HJR9cKf1XtNUGi6szJ6DCpZJ4USMRMJm7Kp4BOj1tcHAwOLnUZlymJPtc3czLrPgYvpPPzo9HlbMgUCwCLgqZOB5SAaO7Y6uYMWvK1GNX7Dj35pzsCPjch3kTgdH5+GbhIG0q9h4YUaJEV1C0l4C2Nwezz\/fZ+dk9UtwZAvkn4LOQ0f2Xhh+cCpyg0ZhadsXOvyeLbaHPfZg3Aia8d0a5fSfC54Z\/RsAU+4Gn9hCwCYSFjNkB23VK4cXxjg0tnS2XnrPS9aphv0UAAeNAc7C3VFdzwz+HRUk4AlNqH4yoTP+4Y48OYMNECEAgJgEVMv\/4RK+8+fKjwa7K779gVI47tSPm1fk+zQwtmS0L4u6Ine9aFdc6OzcqyYwl14h5E4EJR1wqLdds7yxqlq6Ocp69nLX5Xqct6ocDAhAoHgEVMkceulpk2Quy\/JQOOfHC0UDQ+HD84GdH5LP7frxo7yXyZNzz7NjYmOjHPsr1c+7V8JjFhRQwKm727dsX5MDovhaVojUbNmwIpqu1tLQE0IjAuNrcsRsC9SMw\/+N98vrMnaKCxgwr+SJkSPitXzvJoiQ7AqPT7FXMIGCy8ETMe8YdQgpvc67F2wnA4W3OfR4\/jImW0yAAgTIE5p\/6ssw\/deyvXRUyLm9PEK4mQsb9pu9zH+ZNBCY8ZFQqiRcB4\/4DSQ0gkDcCb744JW\/MfSMQMhqFOf7MK70TMmbmElOw89b6ytuDgHHAX0mmUUcNIc3OzkauBeOz8x1wKyZCwCkCOpx09Ic7AjHjW6KvOiI8c0mnYJMjk+8m6nMf5k0ERptQqYXsonYK1QjNzp07g5bHQnb5fgCxDgKuEbBnLPmW6IuQcas1ImDc8lddrfXZ+XUFRWEQgMASAq\/+zZC8Pn2nLDtxuXf5MQgZNxq8z32YVxGYRjQnn53fCF6UCQEILCVgEn19zI9ByOS7xfvchyFgKrQ9n52f78cO6yDgF4FwfsyKznFv1o8xniJHJn9t1uc+DAGDgMnfE4dFEPCYwOv\/5wcyf2hwYf0Yn6ZdRwkZZi1l25gRMNnyz\/TuPjs\/U7DcHAIFJ+D7sFJ4aAkhk02D97kPIwJDBCabp4q7QgACQRTmlz\/4U3nrtb\/zblsC27320BJCJt2Gj4BJl3eu7uaz83MFGmMgUGACb7z0qBx9stfrYSUiMtk0cJ\/7MCIwRGCyeaq4KwQgsISAPazk027X4YqGtyhg9+vGPQwImMaxzX3JPjs\/9\/AxEAIFJBBeBG\/FJePeUlAhc+34j+WRp4+IDi3dsf5Dcuk5K72tbxYV87kPIwJDBCaLZ4p7QgACFQjobtc6W8nXtWPs6v\/gZ0fks\/t+HGxV0H3x6cH2BCpoOGongICpnaGzJfjsfGedguEQKAgBjcYc+fZ1Isuf9DrJ17hz+MFnZPjBqeBH9lmqTyP3uQ8jAkMEpj5PCaVAAAINI6BJvq88dLW3WxKEwRkhw4yl2psUAqZ2hs6W4LPznXUKhkOggATefHEq2OV6\/qkxL3e6Drs0PPWa\/JjqGr3PfRgRGCIw1T0VXAUBCGRCQFfyPfpEr8iyF7zcIDJKyJhE30s+uFLu7P4Q+TEJWh4CJgEs30712fm++Yr6QKBIBIoy5dr41E70JT8mfkv3uQ8jAkMEJv6TwJkQgECuCNhTrpvO6xEf91UiP6a2JoeAqY2f01f77HynHYPxEIDAAoGiRWPIj4nf+H3uw4jAEIGJ\/yRwJgQgkFsCRcuNUUfYC+ExrBTdNBEwuX1kG2+Yz85vPD3uAAEIpE3Ajsas6BwPZiz5fjDturSHfe7DiMAQgfH93Ub9IFA4AhqN0XVjln\/ghELMVDLRmPHHng8WwmNbgnebPAKmcI9\/MZxfYLdSdQh4T0DXjfnlw3+6sIqvz3sq2c5kWGlx00bAeP+oFzP8VmC3UnUIFIbA0ce+HuyptPw32sTnHa7DDrWHlYq82zUCxpFHff\/+\/TI4OBhYu337dunq6lpi+fz8vAwMDMiBAwcWfbdq1SrZvXu3rF69etHvfXa+I27FTAhAoEYCdjSmKNOtzbCSCpnxx16Qoi6C53Mf5k0OzKFDh6Svr09GRkaCR938PyxIot4Dw8PDwa\/7+\/uXfO2z82t8J3I5BCDgGIFX\/2ZI3jiysxBbEdiuKfIieD73Yd4IGI2+TExMyNDQkDQ1NYmKkra2tsgojN2w1bl67q5du6S5uRkB49gLGXMhAIFkBIqY4GsI2cNKRdlbCQGT7PnI5OxwFKVcVMUYaIaTOjs7Swod4\/yenh5pb28PLm1tbQ0+HBCAAARcJFDUBF8zrGT2VvJ17Zjp6WnRjx4zMzPS29sre\/fulY6ODheba0mbvYnAhCMuGpGZmpqKHBYyNCpFX\/Q8I2Bsgipm9MMBAQhAwGUCRR1SUp\/5HI0ZGxsT\/dgHAibHT2o1AiZOlMYImNHRUWlpaSECk+M2gGkQgEByArqf0pGHri7M7tY2IV+nXNsRmMnJyUDMIGCSPxupXZF0CMkMH3V3d5cNq\/k8fpiac7gRBCCQewJHH79dXp+5U5af0iFFWTPGOEVnKemwko8L4Pnch3kzhBQeMqqUxGvPWio3U8ln5+f+jYqBEIBAqgQ0wffVR7oLt2aMQrajMd0Xny53dn8oVfaNupnPfZg3AibpNGp16vj4+MKspVKNx2fnN+qBoVwIQMBdAkGC7yPdhRxSUq\/5Fo3xuQ\/zRsBowyu1kF3UcFF42jUCxt0XLpZDAAL1JzB3\/x8WbhsCQ9GnaAwCpv7PhjMl+ux8Z5yAoRCAQCYEijxLSYHrAnjrvvqE07kxPvdhXkVgGvGE++z8RvCiTAhAwC8C9sJ3779wVI4\/8yq\/KlihNq7PVPK5D0PAVGi8Pju\/UG8hKgsBCFRNwF74TgWMCpmiHa7mxvjchyFgEDBFew9RXwhAoEoC9pDSis7xYE+lIh0uRmMQMEVqoaG6+uz8AruVqkMAAlUSOPrY1+Xok73y3rPPkfdfMCrHnerX8vTUqq5iAAAgAElEQVRxsJhVfHWH6wPXXhjnkszO8bkPIwJDBCazB4sbQwACbhLQIaW5b\/y2LP\/ACcFwUtHyYtRrGo3RBF898rwxJALGzWesLlb77Py6AKIQCECgkATIiznmdhONyevidz73YURgiMAU8uVLpSEAgfoQ+MWf\/54sb\/lJIbcgMATzPN0aAVOfdu5kKT4730mHYDQEIJA7Aq\/8ry\/JG0d2FnILAuOMvCb4+tyHEYEhApO7lyEGQQAC7hGw14s5qXO8kMm99pBSXhJ8ETDuPUt1s9hn59cNEgVBAAIQEBHNi\/nFX1woJ5y\/UprO65Gm8z5fSC55GlLyuQ8jAkMEppAvGCoNAQg0hoCKmCPfvk6W\/fo\/FDovJi9DSgiYxrRzJ0r12flOOAAjIQABJwmQ3HvMbVmvGeNzH0YEhgiMky9HjIYABPJPQJN7X5+5M1j0rogr9xoPZTmkhIDJ\/3PSMAt9dn7DoFEwBCAAgXcIsHLvMRD2kNIDn7lQLj1nZSptxOc+jAgMEZhUHiJuAgEIFJfA\/I++J69OdAcr9xZ5hpI9pNR\/RZv0X3F2wxsFAqbhiPN7A5+dn1\/qWAYBCPhGQEXM4fv\/MJihVNTtB4xPzc7WaUy19rkPIwJDBMa39yT1gQAEckpAZyi9sOMj8r7OU4P9k1TIFPVIay8lBExRW5iI+Oz8AruVqkMAAhkRCDaCvPc\/FX77gbTyYnzuw4jAEIHJ6DXGbSEAgaISUBHz0tc+Le8972eFXivG+L+RG0IiYIr6lBGBKbDnqToEINBoAi+MXi7LVv5DMM36\/ReMFnb7AeVsplrXOy8GAdPoVpzj8n12fo6xYxoEIFAQAi\/euVHefu+3ETHWVOvnDr8md6z\/UF2mWvvch3k1hLR\/\/34ZHBwMHvvt27dLV1dXyVeAcaqecP7558uuXbukubl5yfk+O78g70eqCQEI5JzA4bu3yhuv\/AXTrN\/x09o7n5BHnj4i9Vgvxuc+LFcCxhYg5nmrJETMeYcOHZK+vj4ZGRkJfmX+v3r16iWPrp67ceNGuf3226Wjo0P0vhMTEzI0NCRNTU2LzvfZ+Tl\/p2EeBCBQIAKvfu9rcvSJXnnvB08s\/Fox6naTF1PrejE+92G5ETAqIvbt27coEjI3NyebNm2S9evXl42mqLPDImR4eFja2toir9Nzp6ampL+\/v+LrwWfnV6w8J0AAAhBIkQAiZjHseuTF+NyH5ULAGKGigkIjIvah8FWMlBriMefqOXoYURL+2Zw3Pz8vAwMD0tnZWVEU6TXG+T09PdLe3h4U09raGnw4IAABCECgvgRUxLzy7euCtWKKvuCdkq1GxExPT4t+9JiZmZHe3l7Zu3fvkv61vp5LvzSvBIwdcSkVZTEC5oorrpC77rpLDh48GCsHxnaNihn9cEAAAhCAQP0JIGIWMzX7KMVN7h0bGxP92AcCpv7tdKHEWoeQwkNGlQTMs88+uxDV0WtnZ2fL5sCMjo5KS0sLEZgGtgGKhgAEIGAI6NYDL\/7578mJv3t64VftNUziJvfaEZjJyclAzCBgGvxs1ZLEW8sQkp0AHE769Xn8sMHupHgIQAACNRFAxCzFd+34j0X3Urqz+0PSffHpFfn63IflYgipogdinBCOuJRL4g1\/pwJm27ZtsmPHjiVTqX12fgysnAIBCEAgUwK2iFl+SoesuGQ8U3vycHOzGaQKGBUy5Q6f+zBvBEySadThxOBSCb\/aKHx2fh4eRGyAAAQgUImAipj\/e\/vlsuIPW9l64B1YJrm3kojxuQ\/LTMDYM4\/OPffcYLq0JtRGHeUWmrPPL7WQnUnc7e7uXsjCtheyW7t2bWT+CwKm0muF7yEAAQikQwARs5SzJvde8KW\/lXLbDyBg0mmfubyLz87PJXCMggAEIFCCACImWsRoXozOUHryP\/\/2khN87sMyi8DYlOuxDkyjnnifnd8oZpQLAQhAoFEEEDHlRUx4DyWf+zAETIWnzGfnN+oFQ7kQgAAEGkkAERNNV6dZh9eK8bkPy1TARE2bjnLL5s2bYy3734gHxmfnN4IXZUIAAhBIg4Audvfynk+T2BuCHV4rxuc+LFMBY7iXG0JK40Eodw+fnZ81W+4PAQhAoBYCiJhoematGN3NevlLP5ENGzawkF0tDc3VaxEwrnoOuyEAgSIQOHz3Vnnlr75EJCbkbLOb9R+f\/Yoc+LMvIGCK8DCE64iAKaLXqTMEIOASgRfv3ChH\/+7riJiQ08yCdyf85Jty\/83r2cyxUY1aF6LbuHFjsCdR+Ii7DkwjbEPANIIqZUIAAhCoLwFETDTP\/d\/\/B9EhpW\/ccLlces7K+kLPuLRc5MCYheY6Oztl3bp1MjAwILronFngrr+\/PzPliIDJuIVyewhAAAIxCcxu+bi89atped+lb7Ji7zvMtA9b\/+8\/J\/v++3\/NrB+N6b7Ep+VCwISTeO29ihT++Ph4yZVyE9c44QUImITAOB0CEIBAhgRUxLz50qPsYm0JGJJ4G9ggwwLG3pgxvG9RA82ILBoBkzZx7gcBCECgNgKImHf5+dyH5SICo6jtDRVt0fKd73xHJiYmiMDU9jxzNQQgAIHCEHjzxSmZveXjIsteCCIx779wVI4\/86rC1N+uKAImBbfbeTBdXV2BoNm5c6esWrVKdu\/eLatXr07BiqW38Nn5mQDlphCAAARSIKAi5tnPnC3LP3BCoUWMz31YLiIwUbtFp9C+Y93CZ+fHAsBJEIAABBwloFsOPL\/l43Liv\/iELG\/5iZzUOS7HndrhaG2qM9vnPiwXAoaVeKtrmFwFAQhAAALlCRgRc8qmz8o\/\/b9vFU7EIGBSeEI0cTfLXJdSVfTZ+Sm4lVtAAAIQyJyArtZ7+J4tsrL7X4ssf7JQIsbnPixXEZiDBw9GNnQWssv8+ccACEAAAk4T0IXudO+kFX\/wz2XZipcKI2IQME4329qM99n5tZHhaghAAAJuEQimV\/9iSk759FXy+nP3FkLE+NyH5SICk+dHwGfn55k7tkEAAhBoBAEVMXqcePnp8tb8tLz\/glGvE3t97sMQMBWeEJ+d34iXA2VCAAIQyDMBs0bM8t9oWxAxK3\/n4TybXJNtPvdhCBgETE0PBxdDAAIQcI2AWSPmhA+vCUSMHisuGXetGrHsRcDEwuTnST4730+PUSsIQAAClQmY6dUf6PtLef35L3q7+aPPfZhXERidij04OBi03O3bt4uu6Bt1mIXzDhw4sPD15s2bRXe9Dh8+O7\/yI84ZEIAABPwlYKZXn3bddnnjyE5pOq9Hms77vFcV9rkP80bAHDp0SPr6+mRkZCRofOb\/UVsQ6MJ5119\/vdx4440Vtyjw2flePaVUBgIQgEAVBMzMpDNuHZdXJ7q92zfJ5z7MGwETXghP91Jqa2uLjMKo2Nm2bZvs2LFDmpubyzZ5n51fxbPOJRCAAAS8I2BmJp2y6Vo5+kSvV9Orfe7DvBEw9m7W+nSFf7afOHu367gCpqenR9rb24NiWltbgw8HBCAAAQi4T8Ak9Z589RY54fxfl\/mnxpwWMdPT06IfPWZmZqS3t1f27t0rHR1+7QPllYCxIy4akZmamorMa7FzZdTB5Vb6NerVfkRVzOiHAwIQgAAE\/CBgknrP2PLdIB\/G5TVixsbGRD\/2gYDJcTsNDxmVEzB67uzsrAwNDUlTU1MQrbF\/DkdrNmzYIKOjo9LS0kIEJsdtANMgAAEI1ELAJPXaIsbFNWLsCMzk5GQgZhAwtbSMBl+bZAgpbIqdABxO+vV5\/LDBLqF4CEAAAs4RMEm9Z331GfnlI92B\/S6vEeNzH+bNEFI44lIuiTdKwJRK6vXZ+c69WTAYAhCAQAoEnv3M2aIr9Z7e95dy5K8\/5vQaMT73Yd4ImLjTqM0aMJ2dncEMJfPzqlWrWAcmhRcDt4AABCCQdwImH0aTek\/8F1cE06tdXSMGAZP31vaOfaUWsjMipbu7O8jCDi9kt3bt2oV8mHBVfXa+I27FTAhAAAKpE7CTepd\/4IRAxJzUOe7cxo8+92HeRGAa1bp9dn6jmFEuBCAAAR8IaD7Maz\/6nmhSryx\/0snp1T73YQiYCk+Zz8734QVDHSAAAQg0kkA4qVenV7s0M8nnPgwBg4Bp5LNP2RCAAAScJmAWuTtpzafktGt3OzczCQHjdPOrzXifnV8bGa6GAAQgUAwCdj7M8Wef49TMJJ\/7MCIwRGCK8QailhCAAARqIGAPJb3x0qPObPyIgKnB6a5f6rPzXfcN9kMAAhBIk4BZH2bVlu\/K\/FNfdiKp1+c+jAgMEZg0n3\/uBQEIQMBZAvZQUtOH1wT5MG++\/GiQ1Lvsffnc4BcB42xzq91wn51fOx1KgAAEIFAsAvZ+SUbE5Hlmks99GBEYIjDFevtQWwhAAAI1EtB8GD10KOmtX03nOqkXAVOjs12+3Gfnu+wXbIcABCCQFQEztVq3Gjj5k7eISerN43YDPvdhRGCIwGT1DuC+EIAABJwl8Or3viYv3rkxWKVXh5Jef+5eOfpEb+62G0DAONvEajfcZ+fXTocSIAABCBSXgD21WiloUq\/mw6zoHM9NUq\/PfRgRGCIwxX37UHMIQAACNRAwQ0knfHhNkA9jRIz+u+KS8RpKrt+lCJj6sXSuJJ+d75wzMBgCEIBAzgiEp1abpN685MP43IcRgSECk7PXAeZAAAIQcIuA2bX6rK8+I8tPa1tI6j2pc1yOO7Uj08ogYDLFn+3NfXZ+tmS5OwQgAAF\/CNir9GqtNKFXE3uzFjE+92FEYIjA+PMGoSYQgAAEMiIQHkpSM0xSr67Um9WBgMmKfA7u67Pzc4AXEyAAAQh4Q8AMJf3mPW8HdcrDInc+92FEYIjAePPyoCIQgAAEsibw86vfIyet+ZScdu3uwBSzyF1WQ0kImKxbRIb399n5GWLl1hCAAAS8JBBe4E4rmeXO1T73YURgiMB4+RKhUhCAAASyIhBe4E7t0HwYPdJeHwYBk1UrSHjf\/fv3y+DgYHDV9u3bpaurq2IJhw4dkr6+PhkZGZHVq1cvOd9n51eEwwkQgAAEIJCYQHivJC0gq\/VhfO7DvInA2EJEG0s5UWJa4\/z8vAwMDMjjjz8uu3fvRsAkfky5AAIQgAAEoggcvnurHL5ny8JeSXpOFvslIWAcaJ8afZmYmJChoSFpamqS4eFhaWtrKxuFUcfqeXoQgXHAyZgIAQhAwCECOpSkh9lmQP+f9tRqBIwDDcYIkf7+\/sDa8M\/hKszNzcnWrVulu7s7OBcB44CTMRECEICAQwTM2jA6I0lnJpnjyF9\/TJY1taaSD4OAcaDBhCMuGpGZmpoSI2jCVdDv9bjoooti5cD09PRIe3t7cE1ra2vw4YAABCAAAQiUI\/DinRtFZyaZtWH03EZPrZ6enhb96DEzMyO9vb2yd+9e6ejIdluDercUb3JgkggYzZfZs2eP3HTTTYGT4yTx2uBVzOiHAwIQgAAEIFCJgG4z0PThNQtrw+j5jZxaPTY2JvqxDwRMJS9l+H2SISQ997LLLgvUaNxZSKOjo9LS0kIEJkMfc2sIQAACLhKIWhtG69GoqdV2BGZycjIQMwiYHLec8JBRqSRezX3ZtGmTHDx4cEltohzs8\/hhjt2JaRCAAAS8IhCV0JvG1Gqf+zBvhpCqmUatT0fcCIyP6tWrtwOVgQAEIJBjAlGbPaq5jZ5ajYDJcaOwTSu1kJ1Z70VnHIWTmBAwjjgXMyEAAQg4TiBqhV4zlPTW\/LQ0YtdqBIzjjaYW8312fi1cuBYCEIAABJIRiFqh15TQqKnVPvdh3gwhJWtG8c\/22fnxKXAmBCAAAQjUg4BZofesrz4jy09rWyiyUVOrfe7DEDAVWqTPzq\/Hw0gZEIAABCCQjEDUtGotoRFTq33uwxAwCJhkTx5nQwACEIBATQRKTavWQus9tRoBU5Or3L7YZ+e77RmshwAEIOAugahp1Vqbek+t9rkPIwJDBMbdNwCWQwACEHCUQKlp1Vqdek6tRsA42kDqYbbPzq8HH8qAAAQgAIHqCJSaVm2GkuoxtdrnPowIDBGY6p48roIABCAAgZoIlNqtup5DSQiYmlzk9sU+O99tz2A9BCAAAfcJlIvC1GNWks99GBEYIjDuvwGoAQQgAAFHCZSLwpihJP13xSXjVdUQAVMVNj8u8tn5fniIWkAAAhBwm8CLd24UnVr9m\/e8vaQiZoG7pvN6pOm8zyeuqM99GBEYIjCJHwgugAAEIACB+hEwWwycdu1uOWnNp5YUXMtQEgKmfn5yriSfne+cMzAYAhCAgKcEykVhahlK8rkPIwJDBMbT1wHVggAEIOAOgUpRGDOU9P4LR+X4M6+KXTEETGxU\/p3os\/P98xY1ggAEIOAugUpRGDOUtPJ3HpZl72uNVVGf+zAiMERgYj0EnAQBCEAAAo0lUCkKU81QEgKmsT7Ldek+Oz\/X4DEOAhCAQAEJVIrCmKGkkzrH5bhTOyoS8rkPIwJDBKbiA8AJEIAABCCQDgEThTn56i1y8idvibzp0Sd6g\/2S4gwlIWDS8Vsu7+Kz83MJHKMgAAEIFJxApSiM4jny1x+TZU2tFRe487kPIwJDBKbgrwqqDwEIQCBfBOLkwsQdSkLA5Mu3qVrjs\/NTBcnNIAABCEAgNoFyeySZQn75SLe8+fKj0rzumZLl+tyHEYEhAhP7geJECEAAAhBIh4DZI+mMLd+Vpg+vKXnTuQfODtaF0fVhog4ETDr+qvku+\/fvl8HBwaCc7du3S1dXV2SZ8\/PzMjAwIAcOHAi+37x5s\/T39xfO+TUDpwAIQAACEGgYAY3C6LFqy3dL3kOTeTWpt9SsJARMw9xTv4IPHTokfX19MjIyEhRq\/r969eolNxkeHg5+p6Jlbm5ONm3aJOvXr48UPD47v370KQkCEIAABOpNQDd41ITeSlEYHUrSI2rHap\/7MG+GkDT6MjExIUNDQ9LU1CQqUtra2kpGYeyGZguacAP02fn1ftgoDwIQgAAE6kvg2c+cHQwh6UaPpY5y2wz43Id5I2DCIqScKLEbgYnAaDSmo2PpokA+O7++jxmlQQACEIBAvQmYKMxZX31Glp\/WVrJ4szZMOKHX5z7MKwFjR1w0IjM1NVUyt0VbgYqcnTt3ytq1axciN6UiMD09PdLe3h583draGnw4IAABCEAAAo0m8POr3yMnrflU2SiM2mASeg+f1iPT09OBWTMzM9Lb2yt79+6N\/CO90bY3svxCCxgDVoXM7OxspIgx6tV2gooZ\/XBAAAIQgAAEGk0gzsJ2aoNJ6P0fj3XIV\/\/n5CKzEDCN9lIN5Vc7hKS3tBOAw0m\/RsCMjo5KS0sLEZgafMSlEIAABCBQHQGNwmgejEZiyh2a0DszMy2vnHVsWvXk5KSMjY0RgakOezpXhYeMkiTxqkjR83ft2iXNzc2LDPZ5\/DAdz3AXCEAAAhColUDcKEw4odfnPsybIaRqp1GbNWFWrVoVmS\/js\/NrfaC4HgIQgAAE0iEQd0q1WmMn9Prch3kjYNRppRayMyKlu7s7SGIKL2QXJ4nXx\/HDdB477gIBCEAAAvUgEGdKtbmPSej9+9evlA0bNjCEVA8HuFaGz+rVNV9gLwQgAIEiE4g7pVoZmYTe\/\/1rn5c\/uY4cmEK2GwRMId1OpSEAAQjkkkDcZF413iT0\/quBGSIwufRmg41CwDQYMMVDAAIQgEBsAprMqxs96sJ2lQ6T0Ltl73tkw+dZB6YSL+++R8B451IqBAEIQMBZAkmSebWSJqH3H39rXM46d+lq886CEBGvkngb4QgETCOoUiYEIAABCFRLIEkyr\/ZhJ\/6kW1b\/y1E5\/syrqr1lLq9DwFRwCwIml+0WoyAAAQgUlsDhu7fK4Xu2yG\/e83ZFBj73YQgYBEzFB4ATIAABCEAgXwTiJvMiYPLlt1St8dn5qYLkZhCAAAQgUDcCs1s+HpS1ast3y5bpcx9GBIYITN0eKAqCAAQgAIF0CMRN5kXApOOPXN7FZ+fnEjhGQQACEIBALAJxknl97sOIwBCBifWgcBIEIAABCOSLQJwNHhEw+fJZqtb47PxUQXIzCEAAAhCoK4E4w0g+92FEYIjA1PWBojAIQAACEEiPQKVhJARMer7I3Z18dn7uYGMQBCAAAQgkIlBpGMnnPowIDBGYRA8LJ0MAAhCAQH4I6L5Iz2\/5uJx27W45ac2nlhiGgMmPr1K3xGfnpw6TG0IAAhCAQN0JlBtG8rkPIwJDBKbuDxMFQgACEIBAegTKDSMhYNLzQ+7u5LPzcwcbgyAAAQhAIDEBHUZSEXPyJ29ZMozkcx9GBIYITOKHhQsgAAEIQCBfBEoNIyFg8uWnVK3x2fmpguRmEIAABCDQMAKlhpF87sOIwBCBadgDRcEQgAAEIJAOgVKL2iFg0uGfy7v47PxcAscoCEAAAhCoikDUMJLPfZhXEZj9+\/fL4OBg4Pjt27dLV1dXZCOYm5uTTZs2ycGDB4Pv165dK0NDQ9LU1LTkfJ+dX9UTwkUQgAAEIJBLAlHDSD73Yd4ImEOHDklfX5+MjIwEDcv8f\/Xq1Ysa2vz8vAwMDEhnZ2cgcMzPq1atkv7+fgRMLh9LjIIABCAAgUoEzDDSWV99Rpaf1hacjoCpRC0H32v0ZWJiYiGSMjw8LG1tbSWjMLbJ4Wvt73x2fg7chgkQgAAEIFBHAj+\/+j2LVuX1uQ\/zJgKjgkUPE0UJ\/1yufSBg6vj0UBQEIAABCGRGQIeR9NCtBYjAZOaGZDcOR1xUlExNTUUOC9klm3yY9evXR0ZrjHrt6emR9vb24NLW1tbgwwEBCEAAAhDIEwEzjDT7J\/tk+Wn\/TGZmZqS3t1f27t0rHR0deTK1Zlu8isDYQ0ZxBIzJf1GKlZJ4bdIqZvTDAQEIQAACEMgTgTdfnBKdjTQydZo89PKJC6YhYPLkpZAtSYeQ4ogXO\/w2OjoqLS0tRGBy3AYwDQIQgAAEJBAwL\/\/6OfL6H\/wXmZyclLGxMSIweW4Y4YhLuSTeSjOP7Hr6nACVZ39iGwQgAAEIVEfAnk7tcx\/mzRBS3GnU2hxU3MzOzpYcNkLAVPfQcBUEIAABCGRPwJ5O\/XdPvyAbNmwgApO9W8pbUGohOxNx6e7ulnPPPXfRInamxPPPP1927dolzc3Ni27is3rNuz+xDwIQgAAEqiNgplP\/6ITfQsBUh9D9qxAw7vuQGkAAAhAoGgGzrcDTH92MgCma8019ETBF9Tz1hgAEIOAuAZMH84vr\/xYB464ba7McAVMbP66GAAQgAIH0CdjrwfxJzyA5MOm7IPs7ImCy9wEWQAACEIBAcgKaB\/Pzj\/5H+Q\/\/7UEETHJ87l+BgHHfh9QAAhCAQBEJaB7MwX88QXp+8DoCpogNAAFTRK9TZwhAAALuEzB5ML\/z+NkIGPfdmbwGCJjkzLgCAhCAAASyJ2DyYP7N358pf\/a1u9kLKXuXpGsBAiZd3twNAhCAAATqQ2D+R9+T57d8PNgX6VNfeQABUx+s7pSCgHHHV1gKAQhAAAKLCWgi74MvnygXbvsOAqZojcMlATM9PS333XefXHnlldLa2pp7V2FvY10EX\/gaArSF4rYFzYO59757ETCNbQL5LN0lAeOSrept7G1sm4cvfA0B2kJx28Lhu7fK4Xu2iC5o19HR0VgQKZfuzWaOjeLm0oPvkq0ImEa12HfLpT00lrFLfF2ylXdDfdutSeRdduPD0nbhpfUtPOPSEDAVHGAe\/J6eHmlvb8\/YXeVvPzMzI729veKCrVoT7G1sc4IvfA0B2kJx24L6\/sJ9V8lp1+6Wk9Z8qrEgUi4dAVMBuI4dqyiYnJxM2TXcDgIQgAAEIFA7ga9\/5DlZ\/bv\/NhAxPh0ImBjeVBGjHw4IQAACEICAawROP\/5N74aP1AcIGNdaIvZCAAIQgAAEIICAoQ1AAAIQgAAEIOAeASIw7vkMiyEAAQhAAAKFJ4CAKXwTAAAEIAABCEDAPQIIGPd8hsUQgAAEIACBwhNAwBS+CQAAAhCAAAQg4B4BBEwZnw0PD8vOnTuDM\/bu3ZubZZgPHTokGzdulNnZWVm7dq0MDQ1JU1NT2danC\/KNj4\/HOrfezTiJvTbzVatWye7du2X16tX1NqlseUns3b9\/vwwODgblnX\/++bJr1y5pbm5Ozd4kthqj5ufnZWBgQDo7O6Wrqys1W\/VGSey12brAd25uTjZt2iQHDx7M7J0Rl2+YrWkE27dvT7VNxLU33HZceDeYRVCzarupPtgZ3QwBUwK8Nj7tTLVD+ulPf7rw\/zQ7pyjT7M5n3bp1sToi8yDFFTv1bItJ7NWX6sTExILI0p\/37duXqihIYq\/dRrRdaHtRURlHUNaDcRJb7fuZzivtziqpvcqzra0t1Q7V5pTEXnOudqz9\/f2BUOvr65ORkZHUBHgSe8PtL9yW69E+K5WRxF4jDpWt7ueT93eDEWa33357YG+Wf0BW8oPL3yNgSnhPX5566ANjHrTu7u7MozDhF2OlB0PrceDAAVmzZo28+uqrqXWuBmtSe213ZNEJ1GJv2p1ANbZqR3D99dfLkSNHZP369amKgyT25uGZS2Kvnrtt2zbZsWNHqhG4cs9LpXeDuTYsDtLq0JLytQVh3t8N4T\/GtD3fdtttcs0116QmaNPyY5b3QcBE0A+H2LMMuVf6S6lSp6nfm79Y7OhGWo0ubF8le7MWMLXYa4veNPhWY6vaePHFF8s3v\/nN1IeQktibVadq+y2JveEOKw3\/1\/puMNdnZXsSvlERmLTfZ0nsjRIwOmybhz+Cs2ibjbonAqaMgLEbW9bhbGNm+K+quH\/5ZfmSsnNv4tqr9U17SEbvWQ1fM0SX9rh8UluV\/Z49e+QLX\/iCbN26NRMBE7ct2LkRpu2nnYeWhK8+X1NTU4GpWeXNJbE36+hLNc+a+UNSI8qbN28OouNpHkn4hnf\/Nj+nPWybJp8s7oWAQRBogJMAAAQQSURBVMA0tN0leehtQ7RD+MpXvpJ6Em+19poX8g033JCazUlstUPYra2tsXKn6t0wktir59oswz\/X27ao8pLYa\/KKjMjKu732H0Qm1y\/t\/L4kfKNyStK2O4m9ytdOlN6wYUMwhJ9F4nwaz0pW90DAlBEwprG5PITkUpjYtjUL8WJEiP1iTDLklXY7SRLS1nO\/\/\/3vL8rpSvtlmsTe8GOZNtukbaHUkEGajKvhayJHaUczisI3D5GurMRFGvdFwJSgbA8Z5SGh0JgZHoIJ\/1VQqtFkNYSU1N4sZhfYzJLaa1+bdt5GElvt6em2zWmG4pPYW0rApJlDkMTe8HOYxTsjib3KNwtRWO2zlgeBmJRvuK5pz0pLQ0BkfQ8ETAkP2H\/NuD6N2oQz0056C78kK037ziLsXu4v\/Tj22tGatMVXkmmodj2z6riS2OvalPqweE0SuatXJ5CEr97TzEi78cYbM5kZk8TeqCGkNIdrk77LbLGja3RpAq+ZYl8vf1OOsBt1uUbg4kJ2pULCWUVglG+5xapse0tFCdJO3oxrrxGGZiG7tJN4k7DNg4BJaq+dQ5AF26T22gvZuWBvFlORw+\/bJM+aSYTVMlzga7ffLNbgKoLAIQJTBC9TRwhAAAIQgIBnBBAwnjmU6kAAAhCAAASKQAABUwQvU0cIQAACEICAZwQQMJ45lOpAAAIQgAAEikAAAVMEL1NHCEAAAhCAgGcEEDCeOZTqQAACEIAABIpAAAFTBC9TRwhAAAIQgIBnBBAwnjmU6kDARQK6hopuMHnLLbdI2nvyuMgLmyEAARayow1AAAI5IGDv1ZQDczABAhBwgAARGAechIkQ8J2ArsJ82WWXSUdHh+9VpX4QgECdCCBg6gSSYiAAgaUESu27ZO8VpHvF3HbbbXLNNdcEe\/KYaw4cOBAUmNWy8fgTAhDINwEETL79g3UQcJ5A1CaXGnHRo7+\/P9gra8+ePXLTTTcFv9ON7\/QYGhoSFTdpb5LpPHAqAIGCEEDAFMTRVBMCWREI79Qc\/lkFih5dXV2BmOnr65ORkZFMdkjOihH3hQAEkhNAwCRnxhUQgEBCAnbExR4+0hlHdv5L+LuEt+F0CECgQAQQMAVyNlWFQFYEbGFy1113SVtbWxBxCU+fRsBk5SHuCwH3CCBg3PMZFkPAOQJm2GjdunXywAMPLAwRhadPM4TknGsxGAKZEUDAZIaeG0OgWAQ012VwcFDWrl27KEFXKWg0Rg8zA0lnHmmCrx6ImmK1E2oLgbgEEDBxSXEeBCBQEwEVIhs3bpTrrrsuECwqVuzp06ZwplHXhJmLIVAYAgiYwriaikIAAhCAAAT8IYCA8ceX1AQCEIAABCBQGAIImMK4mopCAAIQgAAE\/CHw\/wFasBt0oSTm9gAAAABJRU5ErkJggg==","height":337,"width":560}}
%---
