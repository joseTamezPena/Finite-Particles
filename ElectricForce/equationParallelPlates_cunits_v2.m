%[text] # Two Parallel Charged Plates and a Point Particle
clear;
% The variables
syms a d v r real positive;
assume(d>0)
assume(r>=0)
assume(v>=0)
assumeAlso (v < 1)
assume(a>=0)
taylorOrder = 12;

vp = [-v,0,0];
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
%[text] ## The Cosine Abnormality 

%dprd_l1 = simplify(simplifyFraction(dot(vec_l1,r_vec_l1)^2/dot(r_vec_l1,r_vec_l1))) % Cosine^2 left plate
%dprd_l2 = simplify(simplifyFraction(dot(vec_l2,r_vec_l2)^2/dot(r_vec_l2,r_vec_l2))) % Cosine^2 left plate
%dprd_r1 = simplify(simplifyFraction(dot(vec_r1,r_vec_r1)^2/dot(r_vec_r1,r_vec_r1))) % Cosine^2 right plate
%dprd_r2 = simplify(simplifyFraction(dot(vec_r2,r_vec_r2)^2/dot(r_vec_r2,r_vec_r2))) % Cosine^2 right plate

dprd_l1 = simplify(simplifyFraction(dot(r_vec_l1,r_vec_l1)/dot(vec_l1,r_vec_l1))) % Cosine left plate %[output:05796605]
dprd_l2 = simplify(simplifyFraction(dot(r_vec_l2,r_vec_l2)/dot(vec_l2,r_vec_l2))) % Cosine left plate %[output:34d747be]
dprd_r1 = simplify(simplifyFraction(dot(r_vec_r1,r_vec_r1)/dot(vec_r1,r_vec_r1))) % Cosine right plate %[output:90cd3a58]
dprd_r2 = simplify(simplifyFraction(dot(r_vec_r2,r_vec_r2)/dot(vec_r2,r_vec_r2))) % Cosine right plate %[output:61660ac0]


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
fplot(totfInfRel,plotrange) %[output:4124e227]
hold on %[output:4124e227]
relEq = sqrt(1-v^2) %[output:4a5e14ba]
relEqApx = taylor(relEq,v,Order=taylorOrder) %[output:8cd61783]
fplot(relEq,plotrange) %[output:4124e227]
fplot(relEqApx,plotrange) %[output:4124e227]
xlabel('v/c')  %[output:4124e227]
ylabel('ratio')  %[output:4124e227]
legend({'Force Carriers',"Lorentz Mass","Lorentz Apx"},'Location','northeast') %[output:4124e227]
hold off %[output:4124e227]

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
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l1","value":"\\left(\\begin{array}{ccc}\nv+\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:34f9ec89]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_r1","value":"\\left(\\begin{array}{ccc}\nv-\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:6900aafa]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l2","value":"\\left(\\begin{array}{ccc}\nv+\\frac{d}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:123f7d6a]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_r2","value":"\\left(\\begin{array}{ccc}\nv-\\frac{d}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:05796605]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_l1","value":"\\frac{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 +2\\,d\\,v\\,\\sqrt{r^2 +d^2 }}{d^2 +r^2 +d\\,v\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:34d747be]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_l2","value":"\\frac{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 +2\\,d\\,v\\,\\sqrt{r^2 +d^2 }}{d^2 +r^2 +d\\,v\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:90cd3a58]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_r1","value":"\\frac{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 -2\\,d\\,v\\,\\sqrt{r^2 +d^2 }}{d^2 +r^2 -d\\,v\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:61660ac0]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_r2","value":"\\frac{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 -2\\,d\\,v\\,\\sqrt{r^2 +d^2 }}{d^2 +r^2 -d\\,v\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:6d204237]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce","value":"\\left(\\begin{array}{ccc}\n\\frac{4\\,d}{{{\\left(r^2 +d^2 \\right)}}^{3\/2} }-\\frac{4\\,d\\,r^2 \\,v^2 }{{\\left(-r^2 +{\\left(-1+v^2 \\right)}\\,d^2 \\right)}\\,{{\\left(r^2 +d^2 \\right)}}^{3\/2} } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:9cd5c9e8]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce_x","value":"\\frac{4\\,d}{{{\\left(r^2 +d^2 \\right)}}^{3\/2} }-\\frac{4\\,d\\,r^2 \\,v^2 }{{\\left(-r^2 +{\\left(-1+v^2 \\right)}\\,d^2 \\right)}\\,{{\\left(r^2 +d^2 \\right)}}^{3\/2} }"}}
%---
%[output:769f5ec5]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\left(\\begin{array}{ccc}\n\\frac{4}{d^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:45e059fe]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_one","value":"\\left(\\begin{array}{ccc}\n\\frac{4}{d^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:926989fe]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\frac{4\\,\\pi }{d^2 }"}}
%---
%[output:733de68c]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInf","value":"4+\\frac{8\\,v^2 }{3}+\\frac{8\\,v^4 }{15}+\\frac{8\\,v^6 }{35}+\\frac{8\\,v^8 }{63}+\\frac{8\\,v^{10} }{99}"}}
%---
%[output:84edc3b9]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInfV0","value":"4"}}
%---
%[output:3830c242]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInfRel","value":"1+\\frac{2\\,v^2 }{3}+\\frac{2\\,v^4 }{15}+\\frac{2\\,v^6 }{35}+\\frac{2\\,v^8 }{63}+\\frac{2\\,v^{10} }{99}"}}
%---
%[output:4a5e14ba]
%   data: {"dataType":"symbolic","outputData":{"name":"relEq","value":"\\sqrt{1-v^2 }"}}
%---
%[output:8cd61783]
%   data: {"dataType":"symbolic","outputData":{"name":"relEqApx","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}-\\frac{v^6 }{16}-\\frac{5\\,v^8 }{128}-\\frac{7\\,v^{10} }{256}"}}
%---
%[output:4124e227]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQ9sHcd95yeOXIup48q000YibTBtZF8R9OQ4cMljAtgBGlgXlLq7OC5FAgeHpxYGmrgMavOPrMayUlsUSQsIEaeAYAiq0IKS6thIotwdHN9dkjrmmUmdmMAZra34wiqk4qtjUqnrUK4d5fBbeajhct\/b3ff2z8zsZ4EHm+Ts7G8+v9k3X\/3mNzPv+OUvf\/lLxQUBCEAAAhCAAAQcIvAOBIxD3sJUCEAAAhCAAAQCAggYOgIEIAABCEAAAs4RQMA45zIMhgAEIAABCEAAAUMfgAAEIAABCEDAOQIIGOdchsEQgAAEIAABCCBg6AMQgAAEIAABCDhHAAHjnMswGAIQgAAEIAABBAx9AAIQgAAEIAAB5wggYJxzGQZDAAIQgAAEIICAoQ9AAAIQgAAEIOAcAQSMcy7DYAhAAAIQgAAEEDD0AQhAAAIQgAAEnCOAgHHOZRgMAQhAAAIQgAAChj4AAQhAAAIQgIBzBBAwzrkMgyEAAQhAAAIQQMDQByAAAQhAAAIQcI4AAsY5l2EwBCAAAQhAAAIIGPoABCAAAQhAAALOEUDAOOcyDIYABCAAAQhAAAFDH4AABCAAAQhAwDkCCBjnXIbBEIAABCAAAQggYOgDEIAABCAAAQg4RwAB45zLMBgCEIAABCAAAQRMgj6wsLCg5MMFAQhAAAIQcI3Aj97Zofpueq9rZsfai4CJQSTCZWhoSM3OzsbCpAAEIAABCEDAJgLn\/s0O9a\/Xflh9f3+PurZ1o02mNW0LAiYG4TPPPKP6+\/vV5OSkamtraxp4nhWIyJqamnLCVuGAvXn2BvjmS9ctvrxr+fYGG\/n+ePmcmv7uy+qpn7xT\/crpp9Xj9+1UXV1d+YIouHYETEIBMz09bb3ztdhywVbBjr35vu3wha8mQF+oVl84vXROHfveT9Sx772s\/vP7fqa+9Lk\/Vq6MC2k8hYBBwKTpL5mW5Us1U5zrKoMvfBEw+fYBG\/mKePn0sb9XEoH52h9\/UJ158blgFgEBU0xfsOopLg0Ckq\/z2GOPqdtuu021t7dbxTHKGOzN10Xwha8mQF+oRl8Q8bLjL34QNFbEi+S8uDSGpfUSERiPIjBpnU95CEAAAhDwg4ApXp77s3+32igEjB\/+bagVPju\/ISDcBAEIQAACVhEQ8XLDA\/87WCr9pb7fXmObz2MYERgiMFa9iBgDAQhAAALJCUiiruS8jNzaoUZufd+6GxEwyVlaUXJ8fFx1dHSo3t7emvacOHFC7d69O\/h7T0+POnDggGppaamU861wFkZAAAIQgEBDBMaf+JEaf2K+pniRShEwDaEt5yYRL4cOHVJjY2M1BYw4VModPnw4EC2jo6Nqy5YtamRkBAFTjtt4KgQgAAEIJCRgLpOWaaOoyIuuCgGTEGqZxZaWltSuXbvUpk2bAjO2b99eU8BI9GVmZmY16hL+2WyHz84v0188GwIQaJ4Ax5w0z9DFGiTq8o9v\/Koa+k83xR4R4PMY5k0OjAiY5eXlIJIiEZXu7u5UEZha5X12vosvLjZDAAIXCHDMSbV7QmdnZ7DretyWGT6PYd4IGN2VV1ZWYgWMlD116pQaGBhQZ86cqbvBj3b+4OCgkg4jl3SYuE5T7VeL1kMAAnkTcOmYk7xZVK1+fXRBrc3pzMjc4uJicJ4fG9k50EuSCBiZMjp+\/HiQA9Pa2hrkw8hVLwfGbLqIGflwQQACECiLgM\/\/si6LqSvPjfO9nIknH\/NCwDjg3TgBE\/V3icYMDw+riYkJtXXr1jWtjPpXDhEYBzoCJkLAcwJxg5jnza908+J8b0Zg4qI1LoOs3BRSowLGR\/XqcsfFdghUnUDcIGY7H73wYm5ubo2p9ba1yLNNmqd+Rr2VrEnt0ONNX19fpocBp\/F9mrJJ22VLucoJGAEfNYUkuTBRe8H47HxbOiF2QAAC6Qm4\/t2kBYxM3Xd1daUHkOEdwvKee+5RR44cCaLwcZH8DB\/dUFVpfJ+mbEPGlHhTJQRMlArW+8UIezayK7EH8mgIQKAhAq4PTHECRn9vnzx5MuCjo+B6yl9Dk1xGuWQbDYnmbNu2bTW\/0VysYf7eBF5LrMi9+\/fvVwcPHgxyJc3NT2W1q4gdSSeQVa9nz55VTz31VPC7xx9\/fN3POgJTyx69N9krr7yiPvShD6n77rtPff7zn1e67eFoUBrfpynbUEcs8SbvBEzWLH12ftasqA8CECiOgOvfTXECRv6RqSPjIkx0hEQIywrShx56KIjchP+BKkJjfn5e\/dEf\/VEganSEx6zP3HW9Xg6k9qa5+am58ONP\/uRP1myEqm3RG6Oatl133XU17THbJxEg3QaxXTjt27dP7d27NxBScqXxfZqyxfXebJ6EgInh6LPzs+lC1AIBCJRBIPzdJLuzysfW69rWjUo++orKgdFREr1Duo5cmFGSG2+8cc2ii3CkpJboqFWu1u\/rcdQCQwsYvY9YOJpjChipT+8AL0LEfO6LL7645m+mgImyI824lKasrX2nll0IGASMa30WeyEAgYh\/hetzcWyFEz5ssF4EJupv+ow7ETDm1E44OmIKmP7+\/jU49NSPudo0SQQmPJ0lld55551KC5iw0Ir6We6pZc+rr76qjh07tiYP00xzCC8iSSNK0pS1te8gYBr0jM\/ObxAJt0EAAhYQ8CUCE5XEG54WCkdgTAFTLwITFgVRbquVA2NO3Tz55JPrjp+Raaq0AqaWPeLLWn8TO+6++2517733rm7zkWZcSlPWgm6dygQiMERgUnUYCkMAAnYQcH1gaiYHxhQw4Xr02Xayt9dnPvOZ1RyY8OpT04u1ViHpXBbzvDwRPJJbI\/k3aQRMOAfGtEemkEwBo6NNvb295MDUed0QMAgYO76NsQICEEhFwHcBU28VkilgBFqt1T3m76Omj8IixpzikSkivTu7ma8jeTqSIPzEE0+srhZKMoUkgqeWPeEITDg\/iCmk6FcDAYOASfWlSWEIQMAOAq4LGDsoumlFGt+nKesaDQQMAsa1Pou9EIBAyqW0APOLQBpRkqasa5QQMAgY1\/os9kIAAgiYSveBNKIkTVnXoCJgEDCu9VnshQAEEDCV7gNpREmasq5BRcAgYFzrs9gLAQggYCrdB9KIkjRlXYOKgEHAuNZnsRcCEPBAwMQto87TyWl339Urop599tnVAx+1ffqMpCxOr07a5jSiJE3ZpM+3pRwCBgFjS1\/EDghAIAUB1wcmFwWMHNr48Y9\/XMn+LHKJsHnwwQfV888\/r3bu3Ln6+xRubKhoGt+nKduQMSXehIBBwJTY\/Xg0BCDQKAHXB6Y4AaPbJ3x6enqCbfblMk9\/lv1RZF8W+V3UqdWyX8wVV1yx5m96Qzl9crXU\/cADD6xxQziaoiMwN9xwg3rppZfUnj17lJzXJJGco0ePBvf+zu\/8TiBgwnu4mPvJmMcDJPl9rb6RxvdpyjbaF8u6DwGDgCmr7\/FcCECgCQKuD0z1BIze8E1OnNYCRTaiC5\/+LPhEFMglm86ZO+rK7+TU6rvuuisQFuZp1AsLC2vOU9JuiNq2X0daRCTdeuut6tFHH13d1l+mj+SSYwU6OjrUjh07AjGlD3c0z1kyzzuSeyRyc8cdd6havzfPawp3kzS+T1O2ie5Yyq0IGARMKR2Ph0IAAs0RCA9Mb70yr978p\/nmKs3x7kt\/vUNteE\/H6hPqCZjwzrTyswiQhx9+WE1MTKwKhHAd5hlKV1111ZpTq806owRMPXvMer\/97W+rm2++ORBWWoQ8\/vjjgYDRU0tRgijqwEYpF25rEhekESVpyiZ5tk1lEDAIGJv6I7ZAAAIJCYQHpuW\/2aeWH70\/4d3FF7vy9vvVlX+wN5GAMc8e0lM1Mh0kUz0iYPTW\/eHpGl25TAFFnVqtzxsKC5haBzrq+kwBI78TEfOJT3xiNYrzyCOPrAoYc+pLyppHGOiEX\/m9OU1V6\/e1vJRGlKQpW3yvaO6JCBgETHM9iLshAIFSCFQ5AmMKmPBJzdoZ4ZVGtSIwIpDMaZ8oZ5oCRnJo9u3bF5wM\/frrrwdTV\/rwxY997GPBQY\/6hO24KSndjiihJGcnIWDqv1oIGARMKV++PBQCEGiOgOv\/sm4mB8Yc+M0cGDN3RqaQzEMfawkYiZ7IpQ9ujBMwIix0Mq4+ZLGWgJHIyhe\/+MVg6fX3v\/\/9IFdGnqNXL0kOTK3fkwMT\/34gYBAw8b2EEhCAgHUEfBEwshrIvLQoqLcKyRQw4VOr9dRMvQiM3CORkvPnz6vNmzerb3zjG2tsMFcIyR\/MCIwIGBEmx48fV4cPH1atra2rERjJgTGngz73uc+p5557LpjyCq+W0nbWsr9eh0vj+zRlrevkMQYhYBAwrvVZ7IUABDzYyA4nNk4gjShJU7Zxi8q500sBo8N54YxwE7Gp7kUZayUddoPPzi+ny\/FUCEAgCwJ8N2VB0c060vg+TVnXaHgnYPTcZL1tnc15Uh0OnJmZCTZKkoSuKKGjw5quORh7IQABPwn4PDD56bHsWpXG92nKZmdhMTV5I2B0QtimTZsCctu3b6+5rbPMUepkqjjMPjs\/ru38HQIQsJcA3032+iZvy9L4Pk3ZvO3Oun6vBMzy8nKw5r7ekri49f5MIWXdxagPAhDIg4DPA1MevHyqM43v05R1jZE3AkaDjxMo+u+yJbQsn9PnYcTlwAwODqrOzs7gMe3t7cGHCwIQgEBZBHwemMpi6spz43wvG\/XJR67FxUU1NDSkfEyDqKyAOX369JolcGfOnKmbA2N2bBEz8uGCAAQgUBaBuEGsLLuSPjfuMMek9TRSLrzEOk0dwv2ee+4J9napt1dLmjrTlo3z\/dTUlJKPeSFg0lIuoXzSCIw+bEtMNA\/cCndI3VEmJydVW1sbEZgSfMojIQCB9QTiBjHbmbkoYPQGdMJWnz5dBuc435sRmNnZ2UDMIGDK8FTKZ8YJGKkuvMy6nhqP6ygpzaM4BCAAgUwIuP7dFCdg6m1kd\/bsWfXUU08Fg3J4gzg9UOvv9SuuuEKdPHkyYC5\/k6MAZBM7nT7Q09MTnLFkXrVWsUqdR48eVbfffrv6whe+oA4ePBhsZKfHnRtuuCGI7EtEX9dhjjfhDfEa7QhpfJ+mbKP2lHVf5aaQBLQ4VDqVuYui\/D5qK2mfnV9Wp+O5EIBA8wRc\/25q5igBWayhv6\/NowTM6R0hPDAwoO66665gRaqU06kCUadRS\/laZxdpb8kKVrl0fXIqtWzFoYWJTk2Qk6eHh4eDgyclX1KfWi2\/l4Mg6x1bkKRnpPF9mrJJnm1TmUoImCjVG6Xuw3vAaLHT39\/vZfjNpo6ILRCAQDoC4YHp\/M8X1C9+fiFx08brne9qV5e86+Lih3oCRtqmT46W72X9j86HH344EAU6BSBch\/ldL2chaREhqQFmnVECJi4iZJ5fpOvTYiQq8m9GXiRyI1GeSy65ZDVq04yP0oiSNGWbsamMe70TMFlD9Nn5WbOiPghAoDgC4e+mlRe+oFZeWJu4WZw18U9quX5QtVz\/2dWC9QSDRDrMzUX1dJCIABEw5mnUejrItECmb2688cZEhzmaU0BmbmS4RXoDVIni6Evv4q5PtA4fMtnR0bG6H1mSHeLjKV4okWZcSlM26fNtKYeAifGEz863pRNiBwQgkJ5AlSMwpoC5++671b333rtuRVC9wxzNCIwWH\/XEi3jHnKrS3tKiZMeOHWv2HwtHZMRXX\/\/614Pb5ATqZlcvpRmX0pRN3wvLvQMBg4AptwfydAhAoCECrg9MzeTAhCMdAlDySsxjYmQKaf\/+\/atTNrWmkGQ\/MH1\/LUfUslVHiu677z71+c9\/PrhdjqQRgaSnr8SOffv2qb179yrJgZEk4D179qw7tiZNJ0jj+zRl09hgQ1kEDALGhn6IDRCAQEoCrg9MWhTIaiDz0quI6q1CMgWMjnbolUZ69U+9CIzcI1NP58+fV5s3b1bf+MY31thw5513rkm0DS\/80IW1YJIk3ccffzzYCf7QoUPBn80VUqa9aY6yqdUl0vg+TdmUXbD04ggYBEzpnRADIACB9AR8HpjS0yj3jqyWRydtRRrfpymb9Pm2lEPAIGBs6YvYAQEIpCDg88CUAoMVRREw5bgBAYOAKafn8VQIQKApAgiYpvA5fXMa36cp6xoUBAwCxrU+i70QgICxlNY8aBYw1SCQ5oBGBEw1+kRkK312foXdStMh4DwBWekipwzLWTdc1SPQ2dmp5Iw+2em33uXzGEYEhghM9d58WgwBTwiYh\/Z50qTVZvx4+Zya\/u7L6umXzqprWzeqkVs71DVXbvStmQ23R4RLnHiRyhEwDSN2\/0afne++d2gBBCDgG4HTS+fUse\/9RI0\/MR8Il76b3qtGbn2fb80srD0+j2FEYIjAFPYi8SAIQAACtQho4XLsey8HRRAu2fQVBEw2HJ2sxWfnO+kQjIYABLwjIKJl\/IkfBe368G9tCiIuEn3hap6Az2MYERgiMM2\/IdQAAQhAoAECIlymv\/uTIM9FR1wQLg2ArHMLAiZbnk7V5rPznXIExkIAAt4QkOmiTx\/7+9UE3Yd3\/rb6yPs3edM+mxri8xhGBIYIjE3vGrZAAAIeE0C4FO9cBEzxzK15os\/OtwYyhkAAAl4TEOEiOS4yZXRhSfT7gikjrvwJ+DyGEYEhApP\/G8QTIACBShIICxdWFhXfDRAwxTO35ok+O98ayBgCAQh4RYC9XOxxp89jGBEYIjD2vGlYAgEIOE2AiIt97kPA2OeTwizy2fmFQeRBEICA1wQQLva61+cxjAgMERh73zwsgwAErCaAcLHaPYFxCBj7fbTGwvHxcdXR0aF6e3tjLT916pQaHh5WExMTauvWrevK++z8WDgUgAAEIBBBAOHiTrfweQzzLgIj4uXQoUNqbGwsVsCsrKyo0dFR9eyzz6ojR44gYNx5J7EUAhAogYC5cy4HLZbggAYeiYBpAFrRtywtLaldu3apTZsu7Oa4ffv2WAEjjhXBIxcRmKI9xvMgAAEXCEi0Rbb611v+s4+LC167aCMCxgF\/iYBZXl5WW7ZsCaIq3d3ddQWMlN+3b5\/q6+sLREycgBkcHFSdnZ0Bifb29uDDBQEIQMBXAubp0PL\/IlzY8t8Nby8sLCj5yLW4uKiGhobU9PS06urqcqMBCa30bgpJTwvFCZgTJ04EiG688cZEOTAmTxEz8uGCAAQg4BuB8B4ucjp0302bOavIIUdPTU0p+ZgXAsYBByYRMJK4e\/ToUbVnz55ApSZJ4p2cnFRtbW1EYBzoA5gIAQikJxAlXGTLf06HTs+y7DvMCMzs7GwgZhAwZXslwfOTCBiZMrr55puDcBqrkBJApQgEIOAtgagVRRJxQbj44XJyYBzyY5yA0cm+c3Nz61oVpVB9dr5DbsVUCEAgQwJRibmcU5QhYIuq8nkMq2wOjO5fRGAsetMwBQIQyJVAVGIuJ0Pnirz0yhEwpbsguQFRERj9O1lxFM7CRsAkZ0tJCEDATQIk5rrptyysRsBkQdHROnx2vqMuwWwIQCABgVrTROS3JIDnURGfxzDvppCy7nc+Oz9rVtQHAQiUT4BpovJ9YJMFPo9hCJiYnuaz8216ybAFAhBojkB4NRH7tzTH05e7fR7DEDAIGF\/eU9oBgcoRYDVR5VyeusEImNTI\/LnBZ+f74yVaAoFqETCTcqXlEm3p\/93NSpZCc0HAJODzGEYEhggMbzsEIOAAgVrRlg\/\/1pVs8++A\/8oyEQFTFnkLnuuz8y3AiwkQgEAMATO3hWgL3SUtAZ\/HMCIwRGDSvg+UhwAEcibAEuicAVeoegRMhZwdbqrPzq+wW2k6BKwjIKLFXAJNtMU6FzlpkM9jGBEYIjBOvpQYDQFfCETt2yLJuGw454uHy20HAqZc\/qU+3WfnlwqWh0OgwgSipojYt6XCHSLHpvs8hhGBIQKT46tD1RCAgCYQNUV0betGxWGK9JE8CSBg8qRred0+O99y9JgHAS8IhPdsEdHCFJEXrnWiET6PYURgiMA48RJiJARcIqBFy3d+eFY9\/dJZJaKFKSKXPOiPrQgYf3yZuiU+Oz81DG6AAARqEqiV1\/KR91\/JDrn0m9II+DyGEYEhAlPai8WDIeA6gbBokfawrb\/rXvXLfgSMX\/5M1RqfnZ8KBIUhAIGAQJRoIa+FzmErAZ\/HMCIwRGBsfe+wCwLWENCi5Ts\/XFbHvvdyYBeixRr3YEgdAgiYCncPn51fYbfSdAjEEqgnWjhAMRYfBSwh4PMYRgSGCIwlrxlmQKB8AvWmhxAt5fsHC9ITQMCkZ+bNHT473xsn0RAINEEA0dIEPG61noDPY5iXEZjx8XHV0dGhent7IzvX0tKS2rVrl5qbmwv+3tPTow4cOKBaWlrWlffZ+da\/eRgIgZwIRO3TIo9ig7mcgFNtaQR8HsO8EzAiXg4dOqTGxsYiBczKyooaHR1V3d3dwd\/1z1u2bFEjIyMImNJeMx4MgfwI6G38n35pWZmby11z5UbV\/7ub2aclP\/TUXDIBBEzJDkjyeB1V2bRpU1B8+\/btNSMw4fpOnDihZmZmIqMwPjs\/CVfKQMBVArWScGWfFjaXc9Wr2J2WgM9jmDcRGBEwy8vLSiIpZoQlibMRMEkoUQYC9hOolc\/CNv72+w4L8yGAgMmHay61hqeI4h6iIzc7d+6MjNho5w8ODqrOzs6guvb29uDDBQEIlEvAPOFZixexSO\/Rwsqhcv3D08shsLCwoOQj1+LiohoaGlLT09Oqq6urHINyeqo3ERjNJ42A0WXl3rgkXpO\/iBn5cEEAAsUTqBdlYWqoeH\/wRPsITE1NKfmYFwLGPj+tsyipgEkiXqRyHYGZnJxUbW1tRGAc6AOY6BcBoix++ZPW5E\/AjMDMzs4GYgYBkz\/3pp+QRMDErTwyjfB5\/rBp2FQAgZwIEGXJCSzVVo6Az2NYJaeQZKn1mTNnak4bIWAq947T4JIJ1FoxJMucP\/L+TYpclpIdxOOdJYCAcch1UREY\/bu+vj513XXXrdnETjdt27Zt6vDhw6q1tXVNa312vkNuxVTPCJiCJZx8y4ohz5xNc0ol4PMY5l0EJuue4rPzs2ZFfRCoRSBOsEjyrQgXWT3EBQEIZEfA5zEMARPTT3x2fnavCDVBYC2BqJ1vpYQIFL2RnJ4egh0EIJAfAZ\/HMAQMAia\/N4eaK0MgSYQFwVKZ7kBDLSKAgLHIGUWb4rPzi2bJ8\/whgGDxx5e0xG8CPo9hRGCIwPj99tK6TAggWDLBSCUQKJwAAqZw5PY80Gfn20MZS2wiUCt\/RWwkh8UmT2ELBOIJ+DyGEYEhAhP\/BlDCawL1oiuStyKihS36ve4CNM5jAggYj50b1zSfnR\/Xdv7uHwEtVk4vrajv\/PCsevqls6uNJLrin79pEQR8HsOIwBCB4Q33lICeCvrx8jn1nR8uq\/CGcURXPHU8zYKAQQABU+Hu4LPzK+xW75oeJ1akwey\/4p3baRAEYgn4PIYRgSECE\/sCUMAuAvXEiljKVJBd\/sIaCJRJAAFTJv2Sn+2z80tGy+MTEAjnrMh0kPxOCxUiKwkgUgQCFSbg8xhmVQTmxIkTavfu3Wu62tjYmOrt7S2t+\/ns\/NKg8uB1BJJEVcycFXa1pRNBAAJJCPg8hlkjYES8HD9+fM2J0EtLS8HJ0Tt37ixNxPjs\/CSdnzLZEtDRE1n9oxNrzaiKjqwgVrLlTm0QqCoBn8cwKwSMFiojIyOqq6trTT8T+OPj42uETZEd0WfnF8mxas8KCxVpv7kKCKFStR5BeyFQDgGfxzAETEyf8tn55bxOfj01auonSURFEm3lwwUBCEAgTwI+j2FWCBhxHlNIeXZh6m6GgBlNkXrCe6roukWQhKd+ECrNkOdeCECgWQIImGYJJryfJN6EoCiWOYFaIiUqmiIPlz1V5JIt9kmozdwdVAgBCGREAAGTEUgXq\/HZ+S76oxmbzekeM5JSS6QQTWmGNvdCAAI2EPB5DLNmCskGR0fZ4LPzbWXeqF1JBYrUr\/NPwpEUpnwapc99EICAjQR8HsNKEzDmyqPrrrsuWC49NzcX6f9t27axCsnGN6NAm\/QUj\/xXn+0jj9c\/679rk7RAMaMo8jemewp0Go+CAARKJ4CAKd0F5Rngs\/PLo3rxyWFhoqd26okTHUHR4kR+llwULVCIotjgWWyAAARsIODzGFZaBMZ0bJn7wMgeMx0dHTU3yvPZ+Xm9XFGi5PTSyuoW+LWiJlqYmEIEcZKXl6gXAhCoAgGfx7BKCxgRL4cOHVL1jivw2flJXl5zakYLjwvRkbWCRH4XToY16w9P6Wixcm1rSzCtc\/Fn9kZJ4hfKQAACEEhCwOcxrFQBE7VsOsohd955p5JderO6dMRn06YLS2G3b9\/ubQQmLEC00NAi5MJ\/LxwOqP8rQsT8OYq7KUi0+JD\/mlM5iJKseiz1QAACEGiMAAKmMW6J76o3hZS4khQF5XnLy8tqy5YtanR0VHV3d8cKmMHBQdXZ2Rk8pb29Pfg0eoUTTmuJBS0k9HMk6hFV1qwvifjQ9Zk7wZpRkIvCoyUoGvW3RtvOfRCAAAQgkC+BhYUFJR+5FhcX1dDQkJqenl53VE++VuRfe6kRmPybV\/8JKysriQTM5f\/Qd6EjnLt6VcDE2a7Fhllu4e37w\/fqeqPqDN+zeO6qdcU62t4f\/O6Sd10QVReSWC+ID1OAXPwb0zRx\/uPvEIAABFwlMDU1peRjXggYV71Zw+6kAuaFr\/epT972yYZa\/4ufX1DBta7zK9F\/Px9zXxJjtKDRZS9puRg1eufbYueSd7WtVmWW12Uvlms84pTEVspAAAIQgEA2BMwIzOzsbCBmEDDZsI2s5dSpU2pgYECdOXNm3d\/z2gcmqYDp7++3wvlRosYUSKYYCpc9\/\/PFVa5R9yQVTPVFzgUxpMuYIigspnLsSlQNAQgAofpeAAAeUElEQVRAAAJvEyAHJueuYAqJHTt2BNM6fX19Sm9wJwm8XV1dmVvhmoDJHECNCrWYiRM6Iop0GS2e4oRQWACFI0Hyd4RPUZ7mORCAgO8EEDA5ezicxGvuzSLwjx07pg4cOKBaWi7mdWRhEgImC4rRdZhCRkROWOCExU9S4SPiRkSPnvpC8OTnQ2qGAATcJ4CAydmHYQEjy6vn5+eDpdMCXwTN4cOHVWtra6aWIGAyxZlJZeHojwgf\/bukosecwqoldi69OvuIXiYAqAQCEIBAhgQQMBnCrFWViBS5wqLlySefVDMzM7lEYJI0y2fnJ2m\/7WXqCZ43fzobmG+KoHB7dAQnSuhc+B3Jy7b3AeyDAARqE\/B5DLNmGXU4GqJ3yZW9Wo4cOaK2bt1aSh\/12fmlAC3xoabYiYrs1BI6ZkTn0qsv7AW0Vvggckp0K4+GAATqEPB5DLNCwGjxIom7eSTrNtO7fXZ+M1x8vjdK6OjpqziRc+lVXav5ORuu6no7XweB43N\/oW0QsJmAz2OYFQKm6J1403Q2n52fhgNl1xIwRc5brz5zYarq7VVZ+mfzDh3FQeDQkyAAgSIJ+DyGWSFgxJmSuFtmrkutDuWz84t8iar4LBE5sgIrqcCRFVbhKSqSjavYc2gzBLIj4PMYZoWA0RGYubm5SK\/ltZFdki7is\/OTtJ8y+REICxxJOo6aogonGuufETf5+YaaIeALAZ\/HMCsEjM0dxWfn28y96rZFiRthEp6euihmLiQXk3dT9Z5D+yGwloDPYxgCJqa3++x8XnQ3CaSN3OhpKRE3RG3c9DlWQ6BRAj6PYQgYBEyj7wX3WUggStwkidogbCx0JiZBIAMCCJgMILpahc\/Od9Un2J2eAMImPTPugIAPBHwew4jAEIHx4R2lDQ0S0MvB3\/jxl4MaJJG4XsSGBOIGQXMbBEoigIApCbwNj\/XZ+TbwxQY7CSSN2Fx2zW1BA8ivsdOPWAUBn8cwIjBEYHjDIZCYgClsZOO+N199ZvWwTalEIjR6sz5ETWKsFIRAbgQQMLmhtb9in51vP30sdIFA3DSUudQbUeOCR7HRJwI+j2FEYIjA+PSu0haLCIiw0RGacG6NKWrk\/y+75pMWWY4pEPCHAALGH1+mbonPzk8Ngxsg0CSBpKKGSE2ToLkdAm8T8HkMIwJDBIYXHQKlEogTNeTUlOoeHu44AQSM4w5sxnyfnd8MF+6FQJ4ERNTUWtp9Ycrpwuqnlus\/m6cZ1A0B5wn4PIYRgSEC4\/wLSgP8J2AmCtfLp2Hqyf++QAvTEUDApOPlVWmfne+Vo2hM5QjETT1JlAZBU7luQYNDBHwew4jAEIHhhYeAFwTiojRMO3nhZhqRkgACJiWwsoqfOHFC7d69O3j82NiY6u3trWmKWbanp0cdOHBAtbS0rCvvs\/PL8hPPhUBRBMxcmpUXplYfa264J0u45WcuCPhIwOcxzJsIzKlTp9Tw8LCamJgI+qD+\/61bt0aKkvHxcXX48OFAtIyOjqotW7aokZERBIyPbzBtgsDbBExB88aPH1vdRZjEYLqIrwQQMA54ViIqMzMzq5EUESgdHR2RUZhw2fDPZnN9dr4DbsVECORKwMyjQdDkiprKSyLg8xjmTQRGBItcOooS\/jksSsIRmO7u7kix47PzS3qfeCwErCWAoLHWNRjWIAGfxzCvBIwZcZGoyvz8fOS0kPQDmXIaGBhQZ86cUdPT06qrqyuye2jnDw4Oqs7OzqBMe3t78OGCAAT8JlAvh4ZVTn773uXWLSwsKPnItbi4qIaGhuqOc662tZICRsTN8ePHgxyY1tZWFRet6e\/vX+NfETPy4YIABKpFQAuaqL1oEDTV6gs2t3ZqakrJx7zq\/UPd5rbUs80rASMNjZtCWllZCZJ2zSkjMwE4nPSrIzCTk5Oqra2NCIyrPR27IZADgShBQ0JwDqCpMhUBMwIzOzsbiBkETCqExRYOTxnVSuJtVMD46PxiPcTTIOA3AREzv\/j5gnrr1WdUeMk20Rm\/fW9z68iBsdk7b9uWZhl11BSS5MJE7QXjs\/MdcCsmQsBZAkRnnHWdV4b7PIZ5M4UkPa7WRnY66tLX17earCsRmkOHDgUdlY3svHpfaQwErCPw1ivz6pfveDmIzrz+7EF1yeUbAhv1dBOb6VnnMm8MQsB448r0DfHZ+elpcAcEIJAFgXrRGc5vyoIwdWgCPo9hXkVg8uiyPjs\/D17UCQEIpCOgxcxbr\/yjenPpK2siM\/JDy\/WfTVchpSFgEPB5DEPAxHR1n53PWw4BCNhHQAsac6rpQlSmEzFjn7ust8jnMQwBg4Cx\/gXEQAhUlYAWM1HHHBCZqWqvSNduBEw6Xl6V9tn5XjmKxkDAcwL18mZIAvbc+U00z+cxjAgMEZgmXg1uhQAEyiCAmCmDupvPRMC46bdMrPbZ+ZkAohIIQKBUAqvTTP\/4tDp\/7u8CW1ieXapLrHq4z2MYERgiMFa9bBgDAQg0TqBeZIacmca5unwnAsZl7zVpu8\/ObxINt0MAAhYT0GJm5e9PKHXJy2siM4gZix2XsWk+j2FEYIjAZPy6UB0EIGAbgSgxI0uzL7v2NiUJwFz+EkDA+Ovb2Jb57PzYxlMAAhDwjkB4nxmdL8MOwN65OmiQz2MYERgiMH6+tbQKAhCIJbAamXlh6kLZ8+9VLb\/dG0RlRNhwuU8AAeO+Dxtugc\/ObxgKN0IAAl4RkMMmX\/\/eX6tfnHtWqQ3PBW3TkRnyZdx2tc9jGBEYIjBuv51YDwEIZEpAR2Ve+9txteE3Nq4KGaaYMsVcWGUImMJQ2\/cgn51vH20sggAEbCIQnmIiKmOTd5LZ4vMYRgSGCEyyt4BSEIBAZQnoKSY5k2nD5rOrU0wt1w+yisnyXoGAsdxBeZrns\/Pz5EbdEICAnwRqrWIiV8ZOf\/s8hhGBIQJj51uHVRCAgNUEiMpY7Z5V4xAwbvgpFyt9dn4uwKgUAhCoHAGiMva63OcxjAgMERh73zwsgwAEnCIQFZXReTLsK1OOKxEw5XC34qk+O98KwBgBAQh4SUCiMq9\/\/6B6c+krQftYwVSOm30ew4jAEIEp563iqRCAQCUISFTmjdNfVisvTKlLLt+AkCnY6wiYgoE3+rgTJ06o3bt3B7ePjY2p3t7emlVpp0qBbdu2qcOHD6vW1tZ15X12fqOcuQ8CEIBAWgJ6eun1Zw+qjds2BbczvZSWYvryPo9h3kRgTp06pYaHh9XExETgYf3\/W7duXedxKTswMKAeeugh1dXVpUT4zMzMqAMHDqiWlpY15X12fvpXgTsgAAEINE\/gZ\/\/tAcWeMs1zTFKDz2OYNwImLELGx8dVR0dHZBRGys7Pz6uRkZFY\/\/vs\/NjGUwACEIBAjgQkKiNTSzLFxPRSPqB9HsO8ETAiWOTSoiT8s+4aKysranR0VHV3d9edYtLltfMHBwdVZ2dn8Ov29vbgwwUBCEAAAs0T0NNLImZ+5bcuJ0+mSaQLCwtKPnItLi6qoaEhNT09Hcw4+HR5JWDMiEutKIsWMLfeeqt65JFH1NzcXKIcGNPpImbkwwUBCEAAAtkRIE8mG5ZTU1NKPuaFgMmGbS61hKeM4gTM6dOnVxN35d4zZ87UzYGZnJxUbW1tRGBy8R6VQgACELhIIErIXHbNJ4OkX\/aTie8pZgRmdnY2EDMImHhupZVoZgrJTAAOJ\/36PH9YmrN4MAQgAIEEBMJCRsTLJS3t6vIPTiJkEvCTIj6PYd5MIYUjLvWSeMN\/EwGzf\/9+dfDgwXVLqX12fsL+TzEIQAACpRNY\/pt9auWFL1zIkXl7P5lfvWFSXXq1X3kdWYP2eQzzRsCkWUYtDhURo\/d+qZXw67t6zfpFoT4IQAACeRNAyKQjjIBJx6u00rU2stOJu319fatZ2OZGdj09PZH5LwiY0lzJgyEAAQjUJCBTS69986h67W8PqHd1X70akdEb44HuIgEETIV7g8\/Or7BbaToEIOABgVpChqklBIwH3bv5JiBgmmdIDRCAAATyJKCFzM\/++wNBRGbDb2wMknwRMiTx5tnvrK8bAWO9izAQAhCAQEAgKiKz4aquSq9a8nkM8yaJN6\/312fn58WMeiEAAQiUSUCEzD99aUCdP\/d3auO\/3RTkyFRVyPg8hiFgYt4yn51f5hcMz4YABCCQN4GV57+lXvnSgNqw+WxlT8D2eQxDwCBg8v4OoX4IQAACpRKQpdfLj94fRGM2btsU2CIrllqu\/2ypdhXxcARMEZQtfYbPzrcUOWZBAAIQyJyAzo8xhYwk+l52zW1eCxmfxzAiMERgMv+ioEIIQAACthIQISMRmdf\/7q\/Vr\/zm5UFERoSMr3vIIGBs7YkF2OWz8wvAxyMgAAEIWElAhMyZvR9V51cWgqml4IgCD5de+zyGEYEhAmPllwtGQQACECiCgM6PkZVKl3\/svd6tWELAFNGLLH2Gz863FDlmQQACECiUgF52fe75b3mX6OvzGEYEhghMoV8UPAwCEICArQRk2fVP7v+o2vCejtWl164n+iJgbO1tBdjls\/MLwMcjIAABCDhFwFytJNNKV+36jHpz6SvO5sf4PIYRgSEC49SXC8ZCAAIQKIKATvKV\/\/7a7\/+huvT9\/6LeevUZ5\/aPQcAU0VssfYbPzrcUOWZBAAIQsIaATvKVaaUr+\/\/jajTGlf1jfB7DiMAQgbHmiwJDIAABCNhIwEzylWjMu7o+oFZemHJiWgkBY2OPKsgmn51fEEIeAwEIQMALAqtHEnzgFvXe4b9S\/zzTp87\/fMHqaSWfxzAiMERgvPhioREQgAAEiiCgc2PkWe\/59BGlNjwXRGPktGvZzffSq7uKMCPxMxAwiVH5V9Bn5\/vnLVoEAQhAoBgCOhrz7ls+pa4a+HP1xo+\/vCpkrvjwsWKMSPAUn8cwIjBEYBK8AhSBAAQgAIEwAXPfGInGbPiNjer154aCYrYk+SJgKtxvfXZ+hd1K0yEAAQhkQsBM8L3y9vvVlX+wV6288IXVJN8ruo8Fyb5lXT6PYV5FYE6cOKF2794d9JOxsTHV29sb22dOnTqlhoeH1cTEhNq6deu68j47PxYOBSAAAQhAIBGB1771l+qVLw2ojR+4RW25\/5tBcq8NSb4+j2HeCBhTiEhvqydKdG9cWVlRo6Oj6tlnn1VHjhxBwCR6TSkEAQhAAAJRBMJTSi0fuKX0aAwCxoG+KtGXmZkZdeDAAdXS0qLGx8dVR0dH3SiMOFbKyUUExgEnYyIEIAABywmYU0qb7\/+mEhEj0Zh\/+cFQKTv5ImAs7zBinhYiIyMjgbXhn8NNWFpaUvv27VN9fX1BWQSMA07GRAhAAAKOENCrlHRejJitc2NkyXVRK5UQMA50mHDERSIy8\/PzSguacBPk73LdeOONiXJgBgcHVWdnZ3BPe3t78OGCAAQgAAEI1CKgp5R0XoyU07kx8v+yb8xl13wyc4ALCwtKPnItLi6qoaEhNT09rbq67NqjptmGe5MDk0bASL7M0aNH1Z49ewInJ0niNUGLmJEPFwQgAAEIQKAeAT2lJGWuvH1vMKWUdzRmampKyce8EDAW99M0U0hS9uabbw7UaNJVSJOTk6qtrY0IjMV9ANMgAAEI2ErgzP0fVeee\/5bSeTFi55s\/fUa9NtOX+ZlKZgRmdnY2EDMIGFt7hlIqPGVUK4lXcl927dql5ubm1rUmysE+zx9a7E5MgwAEIOAdAVlmLcutZdM72cFXrrwTfH0ew7yZQmpkGbV0nqQRGB\/Vq3ffDjQIAhCAgOUEopJ7xWQ5iuD1HwxlHo1BwFjeIbR5tTay0\/u9yIqjcBITAsYR52ImBCAAAU8I6OReicIEB0K+feURjUHAeNJpGmmGz85vhAf3QAACEIBA8wSiVijpWrNcbu3zGObNFFLz3Sm6Bp+dnxcz6oUABCAAgXgCskLpzN6Pqg2\/3hEcP2Be5nLrX71hUl16dWNLoH0ewxAwMX3MZ+fHv16UgAAEIACBPAmYy6zDIkae+89P9zW1g6\/PYxgCBgGT57tJ3RCAAAQgEEMgTsToBN9GdvBFwFS4+\/ns\/Aq7laZDAAIQsIpAnIhpdErJ5zGMCAwRGKteYoyBAAQgUFUCWsS89U\/z6tq\/+NE6DI2sUkLAVLU3KaV8dn6F3UrTIQABCFhJIE7EiNFpVin5PIYRgSECY+VLjFEQgAAEqkxAjh6QKyqxV36fdEoJAVPhXuSz8yvsVpoOAQhAwGoC9ZZYa8PNKaV3dx+LXGrt8xhGBIYIjNUvMcZBAAIQqCoBETGn\/\/h9auMHbqkZiTGnlFquH1Qt1392DS4ETFV7DzkwFfY8TYcABCBQPoFaxw6ELau11BoBU74PS7PAZ+eXBpUHQwACEIBAYgJaxGy+\/5uq5QO31LzvzZ8+o16b6QsOhNz0e08F5Xwew5hCYgop8UtEQQhAAAIQKIdAUhGj82LOrywoOYLg2R8q1d\/fr6anp9cdZlxOS7J7KgIGAZNdb6ImCEAAAhDIjcDy3+xTy4\/er+IiMWKAPoLg\/7zxSfWpkccQMLl5xeKKfQ6\/WYwd0yAAAQhAIIKALK+utdFduLjeL+bkd5W6\/vePEYGpWo9CwFTN47QXAhCAgN0E4vaIMa2fn\/uy+pcfDAUCptETrW2lwRQSU0i29k3sggAEIACBCAI6H+bK2+9XV\/7B3rqMfP5HOAIGAcMXBAQgAAEIOEYgaVIvAsYxx2Zprs\/Oz5ITdUEAAhCAQLEEZCrp3PPfUr\/56C9rPtjnMYwIDBGYYt84ngYBCEAAApkRkJ16N\/x6R82dehEwmaF2ryKfne+eN7AYAhCAAARMAnoq6T2fPqLefcun1sHxeQwjAkMEhm8DCEAAAhBwmIDsD\/Pat\/5Sbdn3TbXhPR1rWoKAccSxJ06cULt37w6sHRsbU729vZGWLy0tqV27dqm5ubng7z09PerAgQOqpaWlUurVEbdiJgQgAAEIxBCotbQaAeNA1zl16pQaHh5WExMTgbX6\/7du3brG+pWVFTU6Oqq6u7sDgaN\/3rJlixoZGUHAOOBrTIQABCAAgbUEaq1KQsA40FMk+jIzM7MaSRkfH1cdHR01ozBmk8L3mn\/z2fkOuBUTIQABCEAgIYGoVUk+j2He5MCIYJFLR1HCP9fzfxIBMzg4qDo7O4Nq2tvbgw8XBCAAAQhAwCYC\/\/f2d6ifvf\/31Bv\/4c8DsxYXF9XQ0BBnIdnkpLAt4YiLiJL5+fnIaSHzXp0Ps3PnzshojVav5j0iZuTDBQEIQAACELCJgCTzvvKlAXX3i5vV3GsbV03jNGqbvBSypREBo\/NfpKq4JN7JyUnV1tZGBMbiPoBpEIAABCCglOwN8\/\/+dYP6xR\/+lZqdnVVTU1NEYGzuGGmnkJKIF2mvz\/OHNvsT2yAAAQhAoDECOgqz+f5vBlGY\/v5+BExjKIu5KzxlVC+JN27lkWkxAqYY\/\/EUCEAAAhDIjoBeVn16+xgCJjus+dSUdBm1PF3EzZkzZ2pOGyFg8vERtUIAAhCAQDEE9LLq0\/\/+gPovDxwiAlMM9safUmsjOx1x6evrU9ddd92aTez007Zt26YOHz6sWltb1xjgUgRmYWFBPfbYY+q2225zYpUU9jbe15PcCd8klBov4xJfl2wVj2Bv4\/3SvFOiMIsLC6r3f\/0CAZMNUrdqcUnAuGSr9ALszfddgC98NQH6QjX7grkiaeTQ46qrqytfEAXX7s0+MHlxc+nFd8lWBExePfZivfSHfBm7xNclW\/luyLbfyoqk\/\/oPP1Uf3P8kAiZbtPbXpl98cyM7W63WGxa5YKswxN58exJ84asJ0Beq2xd+7Yf\/U1321T9Tl9z7lOr44EfyBVFw7URgYoDLXKzsYihr6bkgAAEIQAACrhH4Hx\/6kXr3LZ9S7\/n0EddMr2svAiaBO0XEyIcLAhCAAAQg4BqB9172lnfRF\/EBAsa1noi9EIAABCAAAQggYOgDEIAABCAAAQi4R4AIjHs+w2IIQAACEIBA5QkgYCrfBQAAAQhAAAIQcI8AAsY9n2ExBCAAAQhAoPIEEDCV7wIAgAAEIAABCLhHAAFTx2dy6OOhQ4eCEtPT09bsYigHVw4MDAQHUvb09CQ+lPLYsWOJymbdjdPYazLfsmWLOnLkiNq6dWvWJtWtL4295vlbtc7TytP4NLZqO\/TZYN3d3aq3tzdP89bVncZek61UZDvfpaWlNeeslfGdkZRvmK121NjYWKF9Iqm9Yp9Z1oXvBr0Jall9t9AXu6SHIWBqgJfOJ4OpHPD44osvrv5\/+LDHov1mDj47duxQo6OjKm4g0i9SUrGTZZvS2CtfqjMzM6siS34+fvx45CGbWdpo1pXGXrOPSL9Ic8p5FvansdV8nh68ih6s0torPDs6OgodUBvtC7ptMrCOjIwEg+3w8LCamJgoTICn5Wu2NdyXs+ifcXWksVeLQ2Er5\/nY\/t2gxdZDDz0U2Ct8y\/oHZJwfXP47AqaG9+TLUy55YczTrMs+DCv8xRj3Ykg7Tp48qW655Rb12muvFR6BSWuv6Y4yBoFm7C16EGjEVhkI7r77bnX27Fm1c+fOQsVBGntteOfS2Ctl9+\/frw4ePLjuRPuiBog09po2hcWBjfaG22b7d0P4H2PSnx988EF1xx13FCZoi\/Jjmc9BwETQD4fYywy5h80LD5Jxg6b8Xf+LxYxuFNXp0tpbtoBpxl5T9BbBtxFbxcabbrpJffWrX42N3GXdhjT2ljWo1otK1HvXwgNW1uyS1JeGbzgiZ\/t3Q1QEpmib0\/CNEjASLe\/r67MmFSFJn7K9DAKmjoAxO1vZ4WxtZjjikvRffmV9wTZqr7S36CkZeWYj9uopuqLn5dPaKn3l6NGj6k\/\/9E\/Vvn37ShEwZhi9Xt818x103y86pyQNX3m\/5ufnA1PLyptLY69mWqZQTGuv\/oekRJTvvPPOIDpe5JXG3vDp3\/rnoqdti+RTxrMQMAiYXPtdmpc+\/K\/CL37xi4Un8TZqrxY\/99xzT2E2p7HVDGG3t7cnyp3KumOksVfKmizDP2dtW1R9aezVeUVaZNlur\/kPIp3rV3R+Xxq+UTklRdudxl7hayZK9\/f3B1P4cfmKRfRrn56BgKkjYHRnc3kKSTevzAiM+UUTN+WlX\/wyxIsWIWnt1YyL7idpQtpS9tvf\/vaanK6iv0zT2Bt+LYtmm7Yv1JoyKJJxI3x15KjoaEZV+NoQ6fJJsITbgoCp4V1zysiGhEJtZjjsHv5XQa3OWpaASWtvGasLTGZp7TXvLTocn8ZWc3m6aXORofg09tYSMEXmEKSxN\/welvGdkcZe4VuGKGz0XbNBIKblG25r0avSfBYuum0ImBpeNv814\/oyah3VKDrpLfwlGbfsu4ywe71\/6Sex14zWFC2+0ixDNdtZ1sCVxl7XltSHxWuSSGPWA0wavvJsvSLt3nvvLWVlTBp7o6aQipyuTftdZoqdlpaWYMpWL7HP2u9Vrg8BU8f7Lm5kVyskXFYERvDW26zKtLdWlKDo5M2k9mphuHv37qAXFZ3Em4atDQImrb1mDkEZbNPaa25k54K9ZSxFDn\/dpnnXdCJsWe9a2v5g9t8y9uCqgrBBwFTBy7QRAhCAAAQg4BkBBIxnDqU5EIAABCAAgSoQQMBUwcu0EQIQgAAEIOAZAQSMZw6lORCAAAQgAIEqEEDAVMHLtBECEIAABCDgGQEEjGcOpTkQgAAEIACBKhBAwFTBy7QRAhCAAAQg4BkBBIxnDqU5EHCRgOyhIgdM7t27VxV9Jo+LvLAZAhBQCgFDL4AABEonYJ7VVLoxGAABCDhBAAHjhJswEgJ+E5BdmG+++WbV1dXld0NpHQQgkBkBBExmKKkIAhAIE6h17pJ5VpCcFfPggw+qO+64IziTR99z8uTJoLqytuXHmxCAgN0EEDB2+wfrIOA8gahDLiXiItfIyEhwVtbRo0fVnj17gt\/JwXdyHThwQIm4KfqQTOeB0wAIVIQAAqYijqaZECiLQPik5vDPIlDk6u3tDcTM8PCwmpiYKOWE5LIY8VwIQCA9AQRMembcAQEIpCRgRlzM6SNZcWTmv4T\/lvIxFIcABCpEAAFTIWfTVAiURcAUJo888ojq6OgIIi7h5dMImLI8xHMh4B4BBIx7PsNiCDhHQE8b7dixQ33ta19bnSIKL59mCsk512IwBEojgIApDT0PhkC1CEiuy+7du1VPT8+aBF2hINEYufQKJFl5JAm+ciFqqtVPaC0EkhJAwCQlRTkIQKApAiJEBgYG1F133RUIFhEr5vJpXTnLqJvCzM0QqAwBBExlXE1DIQABCEAAAv4QQMD440taAgEIQAACEKgMAQRMZVxNQyEAAQhAAAL+EPj\/YJlLOLKy92MAAAAASUVORK5CYII=","height":0,"width":0}}
%---
