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
R = (d^2+r^2)^(1/2) % total distance from a point in the plate to the charge %[output:42ee9aee]

% Get two points per circle. 
vec_l1 = [d,r*cos(a),r*sin(a)]/R  % velocity vector of light from a point 1 in left positive plate  %[output:54318feb]
vec_l2 = [d,-r*cos(a),-r*sin(a)]/R  % velocity vector of light from a point 2 in left positive plate  %[output:1ddbdec5]
vec_r1 = [-d,r*cos(a),r*sin(a)]/R  % velocity vector of light vector from a point 1 right negative plate %[output:47dd6bef]
vec_r2 = [-d,-r*cos(a),-r*sin(a)]/R  % velocity vector of light vector from a point 2 right negative plate %[output:5f663134]

r_vec_l1 = (vec_l1 - vp) % left plate relative velocity of particle (v in light speed fraction =v/c) %[output:1e1bc6c0]
r_vec_r1 = (vec_r1 - vp) % right plate relative velocity of particle %[output:273974d1]
r_vec_l2 = (vec_l2 - vp) % left plate relative velocity of particle (v in light speed fraction =v/c) %[output:89b2099e]
r_vec_r2 = (vec_r2 - vp) % right plate relative velocity of particle %[output:16cfefbf]
%%
%[text] ## The Cosine 

dprd_l1 = simplify(simplifyFraction(dot(vec_l1,r_vec_l1)^2/dot(r_vec_l1,r_vec_l1))) % Cosine^2 left plate %[output:13cd5908]
dprd_l2 = simplify(simplifyFraction(dot(vec_l2,r_vec_l2)^2/dot(r_vec_l2,r_vec_l2))) % Cosine^2 left plate %[output:9d0148b3]
dprd_r1 = simplify(simplifyFraction(dot(vec_r1,r_vec_r1)^2/dot(r_vec_r1,r_vec_r1))) % Cosine^2 right plate %[output:49ca2830]
dprd_r2 = simplify(simplifyFraction(dot(vec_r2,r_vec_r2)^2/dot(r_vec_r2,r_vec_r2))) % Cosine^2 right plate %[output:891b6e8e]

%dprd_l1 = simplify(simplifyFraction(dot(vec_l1,r_vec_l1)/dot(r_vec_l1,r_vec_l1))) % Cosine^2 left plate
%dprd_l2 = simplify(simplifyFraction(dot(vec_l2,r_vec_l2)/dot(r_vec_l2,r_vec_l2))) % Cosine^2 left plate
%dprd_r1 = simplify(simplifyFraction(dot(vec_r1,r_vec_r1)/dot(r_vec_r1,r_vec_r1))) % Cosine^2 right plate
%dprd_r2 = simplify(simplifyFraction(dot(vec_r2,r_vec_r2)/dot(r_vec_r2,r_vec_r2))) % Cosine^2 right plate


%%
%[text] ## The Net Force of the Four Points

dforce = ((dprd_l1/R^2)*vec_l1 - (dprd_r1/R^2)*vec_r1); % net Force in the x direction first points.
dforce = dforce + ((dprd_l2/R^2)*vec_l2 - (dprd_r2/R^2)*vec_r2); % net Force in the x direction second points.
dforce = simplify(expand(dforce),'Steps',100) % Two charges %[output:54936e39]
dforce_x = dforce(1) %[output:8d6c3959]
%%
%[text] ## Check at v = 0, r =0
df_v_zero = simplify(subs(dforce,[v,r],[0,0]),'Steps',100) % Due to the four point charges %[output:3bb373ae]
df_v_one = simplify(subs(dforce,[v,r],[0.01,0]),'Steps',100) % Due to the four point charges %[output:5ade3090]
%%
%[text] ## Integrate all the angles
dforce_xT = expand(taylor(dforce_x,v,Order=taylorOrder));
dcircle = int(dforce_xT,a,0,pi); % From zero to pi i.e., the full circle
%%
dcircle = simplify(dcircle,steps=100); % Net force of the two full circles with charge density
df_v_zero = simplify(subs(dcircle,[v,r],[0,0]),'Steps',100) % Net force of the full zero radius circle (point) %[output:604093cb]
%[text] ## 
%%
%[text] ## Integrate for all values of r

dtotf = simplify(expand(r*dcircle),steps=100)/pi; % Add all the circles 
%dtotf = taylor(dtotf,v,Order=taylorOrder);
totf = int(dtotf,r,0,r);
totf = simplify(totf,steps=100); % The force on a charge due a charge surface density

%%
%[text] ## Set the Plates Radius to Infinity

totfInf = simplify(limit(totf,r,Inf)) % Net force by plate charge density %[output:698d8451]
totfInfV0 = subs(totfInf,v=0) %[output:3fd25b58]

totfInfRel = simplify(totfInf/totfInfV0) %[output:9018417d]




%[text] ## 
%%
%[text] ## Plots
figure(1)
plotrange = [0,0.99];
fplot(totfInfRel,plotrange) %[output:49b3241f]
hold on %[output:49b3241f]
relEq = sqrt(1-v^2) %[output:7f1affa0]
relEqApx = taylor(relEq,v,Order=taylorOrder) %[output:1390390a]
fplot(relEq,plotrange) %[output:49b3241f]
fplot(relEqApx,plotrange) %[output:49b3241f]
xlabel('v/c')  %[output:49b3241f]
ylabel('ratio')  %[output:49b3241f]
legend({'Force Carriers',"Lorentz Mass","Lorentz Apx"},'Location','northeast') %[output:49b3241f]
hold off %[output:49b3241f]

%%
%[text] ## 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":22}
%---
%[output:42ee9aee]
%   data: {"dataType":"symbolic","outputData":{"name":"R","value":"\\sqrt{r^2 +d^2 }"}}
%---
%[output:54318feb]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_l1","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:1ddbdec5]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_l2","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:47dd6bef]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_r1","value":"\\left(\\begin{array}{ccc}\n-\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:5f663134]
%   data: {"dataType":"symbolic","outputData":{"name":"vec_r2","value":"\\left(\\begin{array}{ccc}\n-\\frac{d}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:1e1bc6c0]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l1","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }}-v & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:273974d1]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_r1","value":"\\left(\\begin{array}{ccc}\n-\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }}-v & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:89b2099e]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l2","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }} & -v-\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:16cfefbf]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_r2","value":"\\left(\\begin{array}{ccc}\n-\\frac{d}{\\sqrt{r^2 +d^2 }} & -v-\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:13cd5908]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_l1","value":"\\frac{{{\\left(d^2 +r^2 -r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }\\right)}}^2 }{{\\left(r^2 +d^2 \\right)}\\,{\\left(r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 -2\\,r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }\\right)}}"}}
%---
%[output:9d0148b3]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_l2","value":"\\frac{{{\\left(d^2 +r^2 +r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }\\right)}}^2 }{{\\left(r^2 +d^2 \\right)}\\,{\\left(r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 +2\\,r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }\\right)}}"}}
%---
%[output:49ca2830]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_r1","value":"\\frac{{{\\left(d^2 +r^2 -r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }\\right)}}^2 }{{\\left(r^2 +d^2 \\right)}\\,{\\left(r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 -2\\,r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }\\right)}}"}}
%---
%[output:891b6e8e]
%   data: {"dataType":"symbolic","outputData":{"name":"dprd_r2","value":"\\frac{{{\\left(d^2 +r^2 +r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }\\right)}}^2 }{{\\left(r^2 +d^2 \\right)}\\,{\\left(r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 +2\\,r\\,v\\,\\cos \\left(a\\right)\\,\\sqrt{r^2 +d^2 }\\right)}}"}}
%---
%[output:54936e39]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce","value":"\\left(\\begin{array}{ccc}\n\\frac{4\\,d\\,r^2 +4\\,d^3 +4\\,d^3 \\,v^2 +4\\,d\\,r^2 \\,v^2 -12\\,d\\,r^2 \\,v^2 \\,{\\cos \\left(a\\right)}^2 +4\\,d\\,r^2 \\,v^4 \\,{\\cos \\left(a\\right)}^2 }{{\\left({\\left(2\\,v^2 -4\\,v^2 \\,{\\cos \\left(a\\right)}^2 +v^4 +1\\right)}\\,r^2 +{\\left(1+2\\,v^2 +v^4 \\right)}\\,d^2 \\right)}\\,{{\\left(r^2 +d^2 \\right)}}^{3\/2} } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:8d6c3959]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce_x","value":"\\frac{4\\,d\\,r^2 +4\\,d^3 +4\\,d^3 \\,v^2 +4\\,d\\,r^2 \\,v^2 -12\\,d\\,r^2 \\,v^2 \\,{\\cos \\left(a\\right)}^2 +4\\,d\\,r^2 \\,v^4 \\,{\\cos \\left(a\\right)}^2 }{{\\left({\\left(2\\,v^2 -4\\,v^2 \\,{\\cos \\left(a\\right)}^2 +v^4 +1\\right)}\\,r^2 +{\\left(1+2\\,v^2 +v^4 \\right)}\\,d^2 \\right)}\\,{{\\left(r^2 +d^2 \\right)}}^{3\/2} }"}}
%---
%[output:3bb373ae]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\left(\\begin{array}{ccc}\n\\frac{4}{d^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:5ade3090]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_one","value":"\\left(\\begin{array}{ccc}\n\\frac{40000}{10001\\,d^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:604093cb]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\frac{4\\,\\pi }{d^2 }"}}
%---
%[output:698d8451]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInf","value":"4-\\frac{8\\,v^2 }{3}+\\frac{8\\,v^4 }{15}+\\frac{8\\,v^6 }{105}+\\frac{8\\,v^8 }{315}+\\frac{8\\,v^{10} }{693}"}}
%---
%[output:3fd25b58]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInfV0","value":"4"}}
%---
%[output:9018417d]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInfRel","value":"1-\\frac{2\\,v^2 }{3}+\\frac{2\\,v^4 }{15}+\\frac{2\\,v^6 }{105}+\\frac{2\\,v^8 }{315}+\\frac{2\\,v^{10} }{693}"}}
%---
%[output:7f1affa0]
%   data: {"dataType":"symbolic","outputData":{"name":"relEq","value":"\\sqrt{1-v^2 }"}}
%---
%[output:1390390a]
%   data: {"dataType":"symbolic","outputData":{"name":"relEqApx","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}-\\frac{v^6 }{16}-\\frac{5\\,v^8 }{128}-\\frac{7\\,v^{10} }{256}"}}
%---
%[output:49b3241f]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQ+QXVWd539kgtAImRBgNEnDNAOBsViHPy7bbYMFzujEmZ1kdgaw05mtwWzKzY7ItCP0PyglcSXpTpNduwScLJWNKd10AkIpcXcFrVIGydKLCKnRRYxID3QHxkASxdgwYNj63Xia0zf3vXfve\/fde8+5n1vVlXT3uef+zud37jvf\/p3fOee4N998803hggAEIAABCEAAAg4ROA4B45C3MBUCEIAABCAAgYAAAoaOAAEIQAACEICAcwQQMM65DIMhAAEIQAACEEDA0AcgAAEIQAACEHCOAALGOZdhMAQgAAEIQAACCBj6AAQgAAEIQAACzhFAwDjnMgyGAAQgAAEIQAABQx+AAAQgAAEIQMA5AggY51yGwRCAAAQgAAEIIGDoAxCAAAQgAAEIOEcAAeOcyzAYAhCAAAQgAAEEDH0AAhCAAAQgAAHnCCBgnHMZBkMAAhCAAAQggIChD0AAAhCAAAQg4BwBBIxzLsNgCEAAAhCAAAQQMPQBCEAAAhCAAAScI4CAcc5lGAwBCEAAAhCAAAKGPgABCEAAAhCAgHMEEDDOuQyDIQABCEAAAhAovYAZHh6WtrY26erqojdAAAIQgAAEIOAIgVILGBUvmzdvlg0bNiBgHOmwmAkBCEAAAhBQAqUUMAcOHJDVq1fL\/Pnzg17woQ99CAHD+wABCEAAAhBwiEBpBczBgwdl0aJFMjAwIJ2dnQgYhzotpkIAAhCAAARKKWCM26enp2MJmF880h3c8lsntVbsMXNOWixHfjUV\/F7\/H3XNse6f0zK7rhcOHL1j32\/+pWtCAAIQgAAE0iDQ2toq+uXbhYCpEYGZnJyUp+5+Xyy\/HznpdFl84ktBWf2\/uVp\/87NYlViFHtv3Vn1hgWO+Pyp6jhP79\/bvkj6T8hCAAAQg4BeB9vZ2GRkZ8U7EIGBqCJhHH31UVq5cGTh\/8eLoyEqcrh4WNOb75w6+KipwnjvwalCNCqDnD7wa\/Gt+dmR6Mvi\/EUKLT3w5KFtLGJlojgoajf7YUZ\/jT+8I6jj68\/rbFaftlcqMj4\/L6Ohow2wbsSHJvdibhFbysvBNzizuHbCNS6q+ckXma2zbvn27dHQc\/dz35ULAxBQwRXC+ETTa+fT\/zx98VZ47MC0qaFT0PLf3UdGI0TtaT7fEztGIkIodLVdJ9JjpLRU5OlVmpsGM8Dn6s3RDkEYcFoFtnBcae+NQqr8MfOtnV+tO2NYi1Njvi8y3yLY1Rr2kq5AMtDg5MC45P2yrihwjeh555mDQbPOziamfRAobjfyoyGmf\/3Rk3zIi5vjTOmZEztzTjqp6E9WJ2yldYqttwt64nq2vHHzr4xbnLtjGoVR\/mSLzLbJt9RM\/eicRGIciMLWcrdGXe++9V6666qrYc51G0Jhojvn+kWcOBY8zOT0qaN6awjoqcqIiOnYk5\/jT24M6VOBERXDqsbcWg2b+HnubSVeC6GHS\/ttci6rX7pK9Ltmq1LE3vZ6NgEmPpXM1+ez8Ws54K3pzKJiq0uu7PzkkRtyY+03U5r2LDsuRX03Kexf9Ut678LC88fKjxzzCTEnVEje1bOP3EIAABCBQm4DPY1ipIzC1Xe\/etEGcNqVRJixuVNgcjeIcTUbW66wFJwb\/dl3wa\/n1ryalff6P5PJz58vrL43HFjdJp6XSaBt1QAACEPCFAALGF0\/W0Q6fnV8Hjli3qIjRKI1GbaIiNipsLjtnfiBwNGrz3oW\/nBE0UeImKmqDsInlCgpBAAIlJ+DzGEYEpkbn9tn5Wb7XJjIz9tgLwWPDwkbFzJmnnhhEaM5a0CJdF7wRRG3MNFQtYfOWyPFrmWCWPuJZEICAfwR8HsMQMAiYXN\/YatEaW9Rcds6pgbjRHBu9Xnv+K8G\/lYTNCWdeFfxeE4iJ1uTqYh4OAQjkSAABkyP8vB\/ts\/PzZlvp+bVEjZl+MqJG61FhYyI2eqSDETjmGRqhMUu\/ETVF9Tx2JSWgq3X0i6t8BOIeD+DzGEYEhgiME29+LVHTfek7xRY0RtRUi9YgapxwPUZWIKDCpbe3V3SnVa7yEYh7PAACpnx9Y6bFPjvfZbe+tV\/Nwch8GhU0evUvPfuYZmq0ptIUlIoaM\/10wplXp777sMvMsb1YBNI65qRYrcKaOASSHA\/g8xhGBIYITJz3xYkyKmpMkvDYYy\/OLOm2Vz11X7pwZnm3aZSdVxPOqbFXQDH15EQ3KI2RPg9MpXFinQ1N4vskZes0J7fbEDAImNw6X7MfXGnayU4OjhI0apedUxMlajRKczRac3Wzm0H9EIgk4PPAhMurE0ji+yRlXeOOgEHAuNZn67a30rSTCppqU052pEannqIEjUkQZtqpbvdwY0ICrg9MBw4ckNWrV8uePXtmtXzZsmUyNDQkLS0tCYk0VtzwNLVs2LBBurq6GqrUnLfX3d2d6knQSXyfpGxDjc3hZgQMAiaHbleMR9qCJjzlpIJG96MxwibKYjuX5rXn751Z4m0nB7ec\/4liNBYrvCPg+sBkBEx\/f3+qg3s9jlaWN954o2zdulWWLFkicQ76rec5ad2TxPdJyqZlX1b1IGAQMFn1tcI\/x86hGX5gYsZeE6EJr3IKN0gFzesvPxoImbCgMYnBCJrCdwNnDHR9YKolYIyI2LVrV+CT7du3B0Jn79690tfXN+OnLVu2BP830ZwLL7xQ9GcLFiwIyq5atUr27dsn9s9tJ1cSK3rv+vXrZdOmTUFdO3fulMHBweDWRYsWBWJHlzIPDAzIoUOH5OGHHw5+dt999x3zvYnAVLJHfTk8PCz79++X97znPfLpT39aPvOZz4hpezgalMT3Sco60\/l\/YygCBgHjWp\/NzF4jaOxdg+2E4KgVTrZxCJrMXFXKB7k+MNUSMDqgq\/DQ6SSdZjIREnW2ipLbbrstEDThaRoVGhMTE\/LRj340EDUmwmPXZ09PGUG0cePGIPoSdRmBYYSR1qXX3\/7t3wYCRgWNPsfYEv5eBcx5551X0R67fWqDaYPWqZzWrVsnt9xySyCk9Eri+yRlXXuREDAIGNf6bG72GkETNd1UKzqjRtuCZvrp0Zl2mKXbrHLKzbVOPjg8MJkp0aI2RsW\/OeBVbYzKgTFREhUYKgxM5MKOklxyySVBBMYIjnCkxLQ\/LDoqlav082ocjcAwAqazszPIlwlHc2xxpfWp8LGjQybC8+Mf\/3jW72wBU0lQrVy5ciYqVc1WBExR34gM7PLZ+Rng8\/YRlaIzcZKBDRSTQxNOCj4qZNqDFU4qbrggEGcQG37gWbGnPotGrX9p26x9mapFYKJ+p4N\/W1ubqICxp3bCQsUWMDrI25eZ+rEjLXEiMOHpLK1zzZo1MxGYsNCK+l7vqWTPyy+\/LGNjY7OSl7W9mzdvDsw302fhtoV\/HqefFK1fNGIPEZga9BAwjXSv8txbLTpTa6qpmqAhOlOePpS0pb5EYKKSeMPTQuEIjC1gqkVgwqIginGlHBh76uab3\/ym7N69e0ZghCMwcQVMJXvUl5V+p3bccMMNctNNN81McSUZl5KUTdoH8y6PgEHA5N0HvXt+tehMnKkmBVJpuslEZ0gG9q7bJG6Q6wNTIzkwtoAJ16PiQsWGTjN9\/OMfn8mB0Z\/v2LFjZgrHBl5pFZLJZTF1aj6OCh7NrdH8GzOFFEfAhHNgbHt0CskWMCbapNNS5MBUfjUQMAiYxB+c3BCfgBEzR\/99Mbgx7qom+ylR001EZ+L7wceSvguYaquQbAGjvq20usf+edT0UVjE2FM8OkWk0SG97HwdzdPRBOEHHnhgZrVQHAFjVlCZVVG2PeEITDg\/iCmk6DcYAYOA8fGzvbBtOrqi6eBMroIRM7X2nIkSM\/ozkwxsxAy7AxfW9akb5rqASR1IiSpM4vskZV1DiIBBwLjWZ72x10RnTOKlfcRB3LwZM92kOwTbK5uYavKmm1RsiM8Dk\/\/ea6yFSXyfpGxjVmV\/NwIGAZN9r+OJxxCISgK+7Jz5cvm58yNP1K6E0Ew1RW2kR96MXx3P54HJL0+l35okvk9SNn1Lm1ujVwLG3imx1jkW9hI1e64zjNtn5ze3a1F7vQQqiZmV\/2Zh1aMNws+rljeDmKnXO8W5j8+m4vgia0uS+D5J2azb0ejzvBEw9lp+hWJvdBSGFM5GtzO+ETCNdinuT5NAWMzUkzOj9iBm0vRKMeryeWAqBuHiWpHE90nKFrfF0ZZ5I2DsZW66i2MlURK15l\/Fz7Zt2+Tmm28+5gRUn53vWmctu72ImbL3gNntd\/2zqdYy6mZ6O+nuu2bcePzxx2cOfDT2mch\/rah\/mu1J4vskZdO0MYu6vBEw5mwKs+wt\/L2BWUnAhJflmfLG+T09PdLe3h78WA\/w0i8uCORFICoBWHcB7r504azt2mvZR2SmFqHi\/t71gclFAaOHNv7pn\/5pcGyAXjqe3HrrrfLDH\/5QVqxYMfPzZveaWr6fnJwU\/dJrampKent7Yx070Gy7067fKwGj20ybjlXtLAkVN9oB7IO59NRPc5S6Ddl0FPtnKmb0iwsCRSDQTDHDaqYieDjahlqDWHEtP2pZLQFjf\/YuW7Ys2AVXL\/v0Z90fRfdl0Z9FnVqtf5jOmzdv1u\/MhnJ6gKLeq3V\/9rOfnYUrHE0xf\/hedNFF8swzz8xE6030Xm9+97vfHYw\/4T1c7BzLSrmXcXMyjZG1fD86Oir6ZV9xjh0oep8J21dKARPeIOlTn\/qU3H\/\/\/TOHg0UJmJGREVm8eHHwKyIwrnXz8thbaZ+ZJMuylVZ4NROb5hWvD9UaxIpn8WyLqgkYswGdnjhtBIpu\/BY+\/VlrtKPt9o66+jvdNO76668PhIV9GrVGJ6Ki7lHb9ms9ZsxYunSp3HPPPTPb+usfynrp6df6B\/Ty5csDMWUOd7RzM+3zjvQejdxce+21UunnlU7G1ntr+d6OwIyPjwdiBgFT4Dci7hRSVBO0M9inhEYJGB+dX2B3YlqDBCodZ6BCxhw4GfcRRsyEN83jsMm4BJtTLjyIvbF\/Ql7\/2URzHpZCrcf\/TpvMPaNtpqZqAkbbZm+tbz6jb7\/99uAPTSMQwnXYZyiddtppsxZz2HVGCZhq9tj1PvTQQ3LFFVcEwsqIkPvuuy8QMGYGwDTSFkRRBzYaMRLnzKZ6x6VaYicF1+ZWhTcRmPCUUbWVRWHa4QTgejtKbl7kwRCoQqBS8m\/SfBl9RCUxw7Ls7LtgeGA6ePc6OXjP2uwNifnEU69ZK6d++JZYAib8mWySbnWqRwWM2bo\/PF1jKtcpoKhTq41QCAuYSgc6mvpsAaM\/UxHzl3\/5lzNRnLvuumtGwITTDuwjAypt9ZFkC5A4EZiyjGHeCJiky6g15KcJv3aoUs+qCF8+q9eYnzsU84hApXyZpFNMimT66c\/J6y+NyxsvPxoQajm\/RzjKILvOUuYIjC1gwic1Gw+EVxpVisDoqlV72ifKg7aA0RyadevWBSdDHz58OBhHzB\/MH\/zgB4ODHs0J27WmpEw7ooRS1HhkyiUZl5KUza73pvMkbwSM4qikYisdzW6Svqotf\/PZ+el0IWpxlYAeLvndnxycOWRSd\/5NulmeHZUxu\/+afBmmmJrbM1z\/bGokB8Ye+O30AfsPUp1CsvNcKgkYjZ7oZVaw1hIwKixM0q1JLagkYHRM+vznPx8sEPn+978f5Mroc8zqJc2BqfTzRnJgiMA0991zpnbXPyScAY2huRFgiik39A092PXPpkrTP0YUVFuFZAuY8KIM8wdptQiM3qORkiNHjsjChQvlwQcfnOWL8O7s4T+Cq22Gav8hrQtEnnzyyWDKK7xaythZyf5qnSOJ75OUbahD5nCzVxGYZvDz2fnN4EWdbhNIa4qp0iomojLp9Q8+m9Jj6VpNSXyfpKxrHBAwNTzms\/Nd66zYmx0BojLZsa73SXw21UvO\/fuS+D5JWdfIIGAQMK71WezNmEBaURk1WxN\/yZVJx4E+D0zpEPK3liS+T1LWNWIIGASMa30We3MkMPzAs3JU0LwYHFlQz\/EFan7UcmxdxaRTTFzxCPg8MMUjUN5SSXyfpKxrRBEwCBjX+iz2FoBAWlGZYPO1A189JirDvjK1nezzwFS79eUukcT3Scq6RhUBg4Bxrc9ib8EIaFRGIzIqahqNyuhuv689\/5WghSYio8uyuY4l4PPAhL+rE0ji+yRlXeOOgEHAuNZnsbegBKKiMmkdXcD0kn8CptZhjs3s5uEl1kmepYLgxhtvjDz8N0k9jZRNIkqSlG3EpjzuRcAgYPLodzzTcwJRUZl6d\/sNJ\/0yvXS087g+MLkoYMwGdMrfnD6dx6ucxPdJyubRlkaeiYBBwDTSf7gXAlUJhKMy\/UvbpNEzmMxOv2UXMq4PTLUETLWN7A4dOiQPP\/xwcMJyeIM4sxGeibLMmzdPzK7r+js9CkA3sduzZ09w77Jly0TPWLKvSruza53btm2Ta665Rj73uc\/Jpk2bZMGCBTOnVV900UWyZcsW2bdvn5g67HP5whvi1fvxkcT3ScrWa09e9yFgEDB59T2eWyICaSX9hlcvlTlPxvWBqZGjBPSARLP1v32UgD29o6\/XqlWr5Prrrw9OidZyKiyGhoYk6jRqLV\/p7CLzquouu3qZ+vRUaj1awAiT5557LhAwevJ0X19fcPBka2vrzKnV+nM9CLLasQVxPhaS+D5J2TjPLlIZBAwCpkj9EVtKQEATfrf\/3xfkkWcOzST9Jp1eYhn2sVNIyuTXv5osbA\/6rZNag8M+zVVNwOiga06O1sMW9XsVILfffnsgCjo7OwMREa7DjnDoWUhGROi5QnadUQKmVkTIPr\/I1GfESNRp1nbkRSM3GuWZM2fOTNSmEUclESVJyjZiUx73ImAQMHn0O54Jgd\/sJ\/OCDD8wEdDQPWVUyOhKpiRX1OZ4ZZheCg9MykFXcRX10miZ7ZdqgkEjHbt37w6iJSpgzHSQigAVMPZp1GY6yG63Tt9ccsklsQ5ztKeAjDCKYmgOitQojrl0CkojLuZE6\/Ahk21tbYHQ0ssWNI36KIkoSVK2Ubuyvh8Bg4DJus\/xPAjMIlDp2IJ6ojJmGXYZ8mTCA1OZIjC2gLnhhhvkpptukvDpzdUOc7QjMEZ8VBMvRoDov\/b0jxEly5cvl4GBgZnIUDgio776+te\/HvR7PYG62knTcT4ekoiSJGXjPLtIZRAwCJgi9UdsKTkBnV4yu\/2aPWUQMtGdwvWBqZEcmHCkwwgLEyW57bbbRKeQ1q9fPzNlU2kK6a677goAV8tLqWSriRR9+tOfls985jNBPSbHxkxfqR3r1q2TW265JciN0STgm2++OYja1Hsl8X2SsvXak9d9CBgETF59j+dCoCKBNFYv2XkyPkZkXB+YjCjQ1UD2ZVYRVVuFZAsYE+0wK43M6p9qERi9R6eejhw5IgsXLpQHH3xwlg1r1qyZJWhMDo5OF+mUk7mMYLr11lvlvvvuE00u3rx5c\/Bre4WUba+KnomJiYYSeZP4PklZ1z6SEDAIGNf6LPaWiEAam+P5KmR8Hphc6+JpLY+O2+4kvk9SNu7zi1IOAYOAKUpfxA4IVCXQ6OZ4UUuwXU729Xlgcu1VQMDk4zEEDAImn57HUyFQJwGNylw39lTdy7B9ETIImDo7kAe3JfF9krKuoUHAIGBc67PYC4GAgD29VE\/Cr+tCxgxMPT090t7eTq8oEYGpqSnp7e0N8mx0I71qFwKmRB0j3FSfnV9it9J0jwg0mvDrqpDRpcA6iI2Pj3vkTZoSl4CK1pGRkWCnXwRMXGolK4eAKZnDaa6zBKISfm9f8S65\/Nz5sdrkYrKvihj94iofARUutcSLUvF5DPNqCkmXpw0ODgY9udJhXKab22X1MC+z4yMRmPJ9ENBi\/wiEE351Lxnd6TfO5aKQidMuypSTAALGAb\/renyzcZCaa5+BESVKdAdFewto+3Awu7zPznfArZgIgYYIfPcnh4KN8eo5dwkh0xB6bi4IAZ\/HMG8iMOGzM6qdOxEuG\/4eAVOQNw8zIJASAZ1eMlGZpAm\/YSFjTsBOyTSqgUBTCSBgmoo3ncrtI9W1xvD3YVESjsBUOgcjKtM\/7txjOi2jFghAIC0CjaxcUiHzyyd65Y2XHw1OVX77RSNy\/OnVV4CkZTf1QCAJATs3KsmKpSTPKEJZbyIw4YhLre2a7ZNFqy1Fs7ezNg7TZYv6xQUBCLhJoJGVSypkDj14jcicF2XuaR1y8sUjgaDhgkBRCIyOjop+2VecJddFsT+uHaUUMCpuduzYEeTA6LkWtaI1K1euDJarLV68OOBKBCZu96IcBIpNoBEhM\/3UDnlt6g5RQWOmlRAyxfZ3WayzIzC6zF7FDAKmwN6PO4UUPuZcm2QnAIePOfd5\/rDA7sQ0CGRKICxkLjtnvtzR\/S7RfJla1\/TTn5Ppp4\/+tatCxuXjCWq1ld+7R8DnMcybCEx4yqhSEi8Cxr0XEIshkCWB8BLsOHvJvLF\/Ql4\/8NVAyPh48nWW\/HlWugQQMOnybEptSZZRR00h7du3L3IvGJ+d3xRHUCkEPCFQj5DR6aTD398UiBkSfT3pCI43w+cxzJsIjPaxShvZRZ0UqhGazZs3B12Tjewcf0MxHwJNJKB7yXx8x1PB2Us6pRRnUzx7xRKJvk10DlXXJICAqYnI3wI+O99fr9EyCKRPoB4h88o\/DMlrk3fInJPnkh+TvkuoMQYBn8cwryIwMXyZuIjPzk8MgxsgAAEJ7+4bJyJjEn3Jj6EDZU3A5zEMAVOjN\/ns\/KxfJJ4HAZ8I6JTSdWNPxT6mIJwfM69zjP1jfOoQBW2Lz2MYAgYBU9DXDrMg4AaBpELmtX\/6rkzvHZzZP4Zl12742VUrETCuei4Fu312fgp4qAICEPgNgaTnLTGtRNfJgoDPYxgRGCIwWbxDPAMCpSGQRMjotNIvvvt3cuTV73EsQWl6SLYNRcBky7tQT\/PZ+YUCjTEQ8IxAEiHz+kuPyuEne5lW8qwPFKE5Po9hRGCIwBThHcMGCHhLIImQsaeVOO3a2y6RacMQMJniLtbDfHZ+sUhjDQT8JhAWMpWWX4c3wZt32ZjfYGhdUwn4PIYRgSEC09SXh8ohAIHZBMKrlioJGT3tWlcrsXcMPagRAgiYRug5fq\/PznfcNZgPAacJxBEyGo059I3rReY+SZKv097Oz3ifxzAiMERg8nuzeDIEIBCcsWRviBcVkdEk358\/eA1HEtBfEhNAwCRG5s8NPjvfHy\/REgi4TyAsZG5f8S65\/Nz5Mw17Y\/9EcMr19NOjnHTtvrsza4HPYxgRGCIwmb1IPAgCEKhNIHxoZFjI6E6+h5\/oFZnzIgdE1sZZ+hIImBJ3AZ+dX2K30nQIFJ6ALWQuO2e+3NH9LjlrwYkzdrPkuvAuLISBPo9hRGCIwBTiJcMICEAgmsDYYy\/K8APPBrkyYSFjL7luOb9HOFeJXhQmgIApcZ\/w2fklditNh4BzBFTEqJhRIdO\/tE26L104E5EhGuOcOzMz2OcxjAgMEZjMXiQeBAEINE5AhczwAxNBRbaQITemcbY+1oCA8dGrMdvks\/NjIqAYBCBQQAJhIaPLr\/WyozHzOseCFUtc5SXg8xhGBIYITHnfbFoOAQ8IRAkZjcbovjFz33EiK5U88HEjTUDANELP8Xt9dr7jrsF8CEDgNwQ0L2bssReCqSVdqdR96TvlhkuOk188\/Hczu\/hyplI5u4vPYxgRGCIw5XyraTUEPCQQJWQ+vuCR4Eylub\/TJpxw7aHTSzyGeSVgdu7cKYODg4E7N2zYIF1dXce4dnp6WgYGBmTXrl2zfrdo0SLZunWrLFmyZNbPfVav5XuVaTEEykEgLGSuaXtN\/uak0SAaw3LrcvQB00qfxzBvBMzevXulr69PNm7cGPjN\/D8sSKK67vDwcPDj\/v7+Y37ts\/PL9RrTWgiUj4AKGbP8WqeWPnbiTrn6vK9xFEGJuoLPY5g3AkajL7t375ahoSFpaWkRFSVtbW2RURi776pzteyWLVtkwYIFCJgSvdg0FQJlIWCfs3Tp3B\/Il5asJ8G3JM5HwDjg6HAUpVpUxTTHTCd1dnZWFDrG+T09PdLe3h7c2traGnxxQQACEHCJgAqZ5Xc+Ia\/\/bEKGztwsl53\/TzL3tA4hwdclL9a2dXJyUvRLr6mpKent7ZXt27dLR0dH7ZsdKuFNBCYccdGIzMTEROS0kPFPreiLljMCxvapihn94oIABCDgIgFzztKfTW+TT7Y\/wJSSi06sYvPo6Kjol30hYArs5HoETJwojREwIyMjsnjxYiIwBe4DmAYBCCQjoPkxX37oe\/Klc9fLWacfIsE3Gb7ClrYjMOPj44GYQcAU1l0S5LHoZRJxa4kTM33U3d1dNazm8\/xhgd2JaRCAQIYEVMgcfnyTfOJf\/U+mlDLknsWjfB7DvJlCCk8Z1UritVctVVup5LPzs3h5eAYEIOAGAc2P+ckT35I\/+Oe\/CfaMefKUm+XK9g+5YTxWViTg8xjmjYBJuoxanTo2NjazaqmS9312Pu88BCAAgTCBN\/ZPyC8e6RaZ86J89Rddsv+MNWLOWYKWewR8HsO8ETDarSptZBc1XRRedo2Ace\/FxGIIQKB5BA7c9xfBxnfjh86XweduCUSMHlHA5RYBBIxb\/krVWp+dnyooKoMABLwj8Mo\/DMnrhzbL5KunS\/+P\/oO8MOciuX3Fu+Tyc+d711ZfG+TzGOZVBKYZHdBn5zeDF3VCAAJ+EbBPtu770Wq578XL5LJz5ssd3e8KDo7kKjYBn8cwBEyNvuez84v92mEdBCBQFAJBXsxvTrZ+7F\/+SPp\/tFo06VenlHRqCSFTFE8da4fPYxgCBgFT3DcPyyAAgUIRMFNKc05qlS+8skFGHn49sK9\/aRuJvoXy1FvGIGAK6pgszPLZ+Vnw4xkQgIBfBA4\/9mU5\/GSvvO3sc+XtF42Tulm7AAAgAElEQVTIf3n8HTL8wEQQhTERGb9a7HZrfB7DiMAQgXH77cR6CEAgcwI6pXTgq+8NDoR8+8Uj8s9v\/zMZe+wFhEzmnqj9QARMbUbelvDZ+d46jYZBAAJNJ2DnxZxw5tWBkNG8GIRM09EneoDPYxgRGCIwiV4GCkMAAhCwCfzsC38icxf\/aNYRBCpk9HiCscdeDKaWWHqdX59BwOTHPvcn++z83OFiAAQg4AWBn\/+vzwb7xegRBJoXc\/zpHUG7VMhcN\/aUPPLMIZZe5+Rpn8cwIjBEYHJ6rXgsBCDgEwF7v5hTOsdmRIy28bs\/OSQf3\/EUS69zcDgCJgfoRXmkz84vCmPsgAAE\/CCgeTE\/+\/uL5cQL50vL+T3Scv4nZjVMp5V0xZJeLL3Oxuc+j2FEYIjAZPMW8RQIQKAUBFTEHPrG9TLnt38wKy\/GbrwRMiy9bn6XQMA0n3Fhn+Cz8wsLHcMgAAHnCUQl99qNYsVSNi72eQwjAkMEJpu3iKdAAAKlI6DJva9N3RFsejevc0x0B9\/wFRYynHqdbjdBwKTL06nafHa+U47AWAhAwEkC4Z17zQqlKCFjViyx9Do9V\/s8hhGBIQKT3ptCTRCAAAQiCEz\/8Dvyyu7uYOfe8AqlakKGU68b704ImMYZOluDz8531ikYDgEIOEdARczB+\/4iWKGku\/bq7r3VLnvpta5Y6r50Iade1+F1n8cwIjBEYOp4JbgFAhCAQHICukLpxU3vlpM6Tw8EjAqZWhdLr2sRqv57BExj\/Jy+22fnO+0YjIcABJwkEBwE+ZW\/Oeb4gVqNYel1LULRv\/d5DCMCQwSmvreCuyAAAQjUSUBFzEtf\/Ki87fyfVNwrJqpqll4nB46ASc7Mmzt8dr43TqIhEICAkwReHPmgzJn\/g2CZtX2GUq3GsPS6FqG3fu\/zGEYEhghM\/DeBkhCAAARSJrD\/jlXy5tu+kVjEqBn2YZEsvWYKKeWumW11O3fulMHBweChGzZskK6urooGGFWqBS688ELZsmWLLFiw4JjyPqvXbL3D0yAAAQhEEzh49zp5\/ed\/H2uZdVQNKmSW3\/lEIGhYej2bkM9jWKEiMLYAMS6oJURMub1790pfX59s3Lgx+JH5\/5IlS47p71p21apVctttt0lHR4foc3fv3i1DQ0PS0tIyq7zPzufDFAIQgEBRCLzynS\/K4Sd65W3nnFxzr5hKNrP0+lgyPo9hhREwKiJ27NgxKxJy4MABWb16taxYsaJqNEVdFhYhw8PD0tbWFnmflp2YmJD+\/v6a767Pzq\/ZeApAAAIQyJBAGiJGzWXp9VtO83kMK4SAMUJFBYVGROxL4asYqTTFY8pqGb2MKAl\/b8pNT0\/LwMCAdHZ21hRFeo9xfk9Pj7S3twfVtLa2Bl9cEIAABCCQLgEVMT\/\/xvXBXjFxNryr9vSyCpnJyUnRL72mpqakt7dXtm\/ffsz4mq7nsq\/NKwFjR1wqRVmMgFm6dKncddddsmfPnlg5MLZrVMzoFxcEIAABCKRPIE0RU8al16Ojo6Jf9oWASb+fztTY6BRSeMqoloB57rnnZqI6eu++ffuq5sCMjIzI4sWLicA0sQ9QNQQgAAFDQI8e2P+FP5GT\/\/idsXftrUavTELGjsCMj48HYgYB0+R3q5Ek3kamkOwE4HDSr8\/zh012J9VDAAIQaIhA2iJGjVEho1NLY4+9GJyt1L\/0bOm+9J0N2Vnkm30ewwoxhZSG88MRl2pJvOHfqYBZv369bNq06Zil1D47Pw3u1AEBCECgmQRsETP3tA6Zd9lYKo8ryx4yPo9h3giYJMuow4nBlRJ+9S3x2fmpfApQCQQgAIEmE1AR88+3fVDm\/UVroqMH4phl7yGjEZn7P3axV6de+zyG5SZg7JVH5513XrBcWhNqo65qG83Z5SttZGcSd7u7u2eysO2N7JYtWxaZ\/4KAifP6UwYCEIBA8wk0U8So9fYeMj5thoeAaX7fLOwTfHZ+YaFjGAQgAIEIAs0WMfpIkx+jkZn+pW3SfelCpyMyPo9huUVg7L6Zxj4wzXrbfXZ+s5hRLwQgAIFmEchCxBghM\/zARNAMFTKa7Ovi5fMYhoCp0SN9dr6LLyM2QwACEMhKxNhLr10VMj6PYbkKmKhl01Gv5po1a2Jt+9+M19pn5zeDF3VCAAIQyIKAbnb38raPNiWxN2y\/y3vI+DyG5SpgTCepNoWUxYtQ7Rk+Oz9vtjwfAhCAQCMEshQxaqeLQsbnMawQAqaRDtzse312frPZUT8EIACBZhM4ePc6+fn\/\/mwmkRjTFpeEjM9jGAKmxtvls\/Ob\/cFC\/RCAAASyILD\/jlVy+HtfzlTEmIiMvauv7uhbtGRfn8ewwggY3Yhu1apVwZlE4SvuPjDNeFF8dn4zeFEnBCAAgTwI5CViii5kfB7DCiFgzEZznZ2dsnz5chkYGBDddM5scNff35\/bMeA+Oz+PDxmeCQEIQKBZBPatfb8c+dWknHT5G6nv2BvH5vDxBEU4Z8nnMawQAiacxGufVaTwx8bGKu6UG6dTNVLGZ+c3woV7IQABCBSRgIqYN156NLVTrOtpY5GEjM9jWCEFjH0wY\/jcono6UyP3+Oz8RrhwLwQgAIGiEiiCiDFTS9eNPSWPPHMot5OvfR7DCiFg1NH2gYq2aPnmN78pu3fvJgJT1E8K7IIABCBQMAJv7J+Qfbe8X2TOi0Ek5u0Xj8gJZ16dm5V5RmQQMBm43c6D6erqCgTN5s2bZdGiRbJ161ZZsmRJBlYc+wifnZ8LUB4KAQhAIAMCKmKe+9jZMvcdJxZCxOQVkfF5DCtEBCbqtOgM+nesR\/js\/FgAKAQBCEDAUQJ65MALa98vJ\/\/hh2Tu4h\/JKZ1jcvzpHbm3JsuIjM9jWCEEDDvx5v4+YQAEIAABLwkYEXPa6o\/Lr\/\/l64URMVlFZBAwGXRrTdzNM9elUhN9dn4GbuUREIAABHInoLv1Hrxnrczv\/ncic58slIhptpDxeQwrVARmz549kR2djexyf\/8xAAIQgIDTBHSjOz07ad6f\/2uZM++lwomYZgkZBIzT3bYx4312fmNkuBsCEICAWwSC5dU\/m5DTPnq1vPb8VwopYqKETCNHFPg8hhUiAlPkV8Bn5xeZO7ZBAAIQaAYBFTF6nfzBd8qR6Ul5+0UjhUjsjWqrJvs2etaSz2MYAqbGG+Kz85vx4UCdEIAABIpMwOwRM\/d32mZEzPwPPFxkk6URIePzGIaAQcAU+sXFOAhAAAJpEzB7xJx4wZWBiNFr3mVjaT8m9fpUyIw99oIMPzAR7OwbZ2oJAZO6G9yp0Gfnu+MFLIUABCCQLgGzvPodfV+S1174VC6HP9bboiRCxucxzKsIjC7FHhwcDPrEhg0bRHf0jbrMxnm7du2a+fWaNWtET70OXz47v96Xh\/sgAAEI+EDALK8+4\/oN8vqhzdJyfo+0nP8JZ5pmCxk1un9pm3RfujCIzpjL5zHMGwGzd+9e6evrk40bNwZ+M\/+POoJAN8674YYb5Kabbqp5RIHPznfmLcVQCEAAAk0iYFYmLfzMmLyyuzv3c5PqaWZYyFx2znzpX3q2XH7ufPF5DPNGwIQ3wtOzlNra2iKjMCp21q9fL5s2bZIFCxZU7S8+O7+eF4V7IAABCPhGwKxMOm31dXL4id7CLq+Ow92sWlJRo5GYv7rgePnCDR+W7du3S0dH\/scoxGlD3DLeCBj7NGttfPh7G4h92nVcAdPT0yPt7e1BNa2trcEXFwQgAAEIuE\/AJPWees1aOfHC35bpp0edFjGTk5Oy86EfyIaHX5EjJ50uc371ktz\/sYul8w\/Odd9ZVgu8EjB2xEUjMhMTE5F5LXaujLKottOvicDYXlcxo19cEIAABCDgBwGT1Ltw7beDfJii7xFTjfro6Kjol14qYF79\/eXy1Rs\/QASmqF01PGVUTcBo2X379snQ0JC0tLQE0Rr7+3C0ZuXKlTIyMiKLFy8mAlPUDoBdEIAABBokYJJ6bRFT9D1iopqsERj90mt8fDwQM0whNdg5mnl7kimksB12AnA46ZccmGZ6jbohAAEIFIuASeo9685n5RePdAfGubBHTCWKPo9h3kwhhSMu1ZJ4owRMpaRen51frI8NrIEABCBQDALPfexs0Z1639n3JTn0rfc5tUdMmKDPY5g3AibuMmqzB0xnZ2ewQsl8v2jRIvaBKcZnB1ZAAAIQyJWAyYfRpN6T\/3BpsLzatT1iDEAETK5dKf7DK21kZ0RKd3d3kMQU3shu2bJlM\/kwZVKv8clSEgIQgEC5CNhJvXPfcWIgYk7pHCvswY9MIZWrf8Zqrc\/qNRYACkEAAhAoKQHNh3n1h98RTeqVuU86ubza5zHMmymkZr1fPju\/WcyoFwIQgIAvBMJJvbq82qWVST6PYQiYGm+Zz8735QOGdkAAAhBoFgGzyd0pV35Ezrhuq3Mrk3wewxAwCJhmvffUCwEIQMALAnY+zAlnn+vUyiQEjBddsL5G+Oz8+ohwFwQgAIHyEbCnkl5\/6VFnDn70eQwjAkMEpnyfRLQYAhCAQB0EzP4wi9Z+W6af\/pwTSb0ImDoc7cstPjvfFx\/RDghAAAJZELCnklouuDLIh3nj5UeDpN45JxXzgF+fxzAiMERgsnjveQYEIAABLwjY5yUZEVPklUkIGC+6XX2N8Nn59RHhLghAAALlJqD5MHrpVNKRX00WOqnX5zGMCAwRmHJ\/EtF6CEAAAgkJmKXVetTAqR++RUxSbxGPG0DAJHSuT8V9dr5PfqItEIAABLIk8Mp3vij771gV7NKrU0mvPf8VOfxEb+GOG\/B5DCMCQwQmy3eeZ0EAAhDwhoC9tFobpUm9mg8zr3OsMEm9CBhvulvyhvjs\/OQ0uAMCEIAABAwBM5V04gVXBvkwRsTov\/MuGysEKJ\/HMCIwRGAK8ZJhBAQgAAEXCYSXVpuk3qLkwyBgXOxVKdnss\/NTQkQ1EIAABEpNwJxafdadz8rcM9pmknpP6RyT40\/vyJWNz2MYERgiMLm+XDwcAhCAgA8E7F16tT2a0KuJvXmLGASMD72rzjb47Pw6kXAbBCAAAQiECISnkvTXJqlXd+rN6\/J5DCMCQwQmr\/eK50IAAhDwioCZSvq9e94M2lWETe4QMF51sWSN8dn5yUhQGgIQgAAEahH46TXHySlXfkTOuG5rUNRscpfXVJLPYxgRGCIwtd5Hfg8BCEAAAjEJhDe409vyPLkaARPTcT4W89n5PvqLNkEAAhDIm0B4gzu1R\/Nh9Mp6fxifxzCvIjA7d+6UwcHBoJNs2LBBurq6avbjvXv3Sl9fn2zcuFGWLFlyTHmfnV8TDgUgAAEIQCAxgfBZSVpBXvvD+DyGeSNgbCGinaWaKDG9cXp6WgYGBuTxxx+XrVu3ImASv6bcAAEIQAACUQQO3r1ODt6zduasJC2Tx3lJCBgH+qdGX3bv3i1DQ0PS0tIiw8PD0tbWVjUKo47VcnoRgXHAyZgIAQhAwCECOpWklzlmQP+f9dJqBIwDHcYIkf7+\/sDa8PfhJhw4cEDWrVsn3d3dQVkEjANOxkQIQAACDhEwe8PoiiRdmWSuQ996n8xpac0kHwYB40CHCUdcNCIzMTEhRtCEm6C\/1+uSSy6JlQPT09Mj7e3twT2tra3BFxcEIAABCECgGoH9d6wSXZlk9obRss1eWj05OSn6pdfU1JT09vbK9u3bpaMj32MN0u4p3uTAJBEwmi+zbds2ufnmmwMnx0nitcGrmNEvLghAAAIQgEAtAnrMQMsFV87sDaPlm7m0enR0VPTLvhAwtbyU4++TTCFp2SuuuCJQo3FXIY2MjMjixYuJwOToYx4NAQhAwEUCUXvDaDuatbTajsCMj48HYgYBU+CeE54yqpTEq7kvq1evlj179hzTmigH+zx\/WGB3YhoEIAABrwhEJfRmsbTa5zHMmymkepZR69sRNwLjo3r16tOBxkAAAhAoMIGowx7V3GYvrUbAFLhT2KZV2sjO7PeiK47CSUwIGEeci5kQgAAEHCcQtUOvmUo6Mj0pzTi1GgHjeKdpxHyfnd8IF+6FAAQgAIFkBKJ26DU1NGtptc9jmDdTSMm6UfzSPjs\/PgVKQgACEIBAGgTMDr1n3fmszD2jbabKZi2t9nkMQ8DU6JE+Oz+Nl5E6IAABCEAgGYGoZdVaQzOWVvs8hiFgEDDJ3jxKQwACEIBAQwQqLavWStNeWo2AachVbt\/ss\/Pd9gzWQwACEHCXQNSyam1N2kurfR7DiMAQgXH3EwDLIQABCDhKoNKyam1OmkurETCOdpA0zPbZ+WnwoQ4IQAACEKiPQKVl1WYqKY2l1T6PYURgiMDU9+ZxFwQgAAEINESg0mnVaU4lIWAacpHbN\/vsfLc9g\/UQgAAE3CdQLQqTxqokn8cwIjBEYNz\/BKAFEIAABBwlUC0KY6aS9N95l43V1UIETF3Y\/LjJZ+f74SFaAQEIQMBtAvvvWCW6tPr37nnzmIaYDe5azu+RlvM\/kbihPo9hRGCIwCR+IbgBAhCAAATSI2COGDjjuq1yypUfOabiRqaSEDDp+cm5mnx2vnPOwGAIQAACnhKoFoVpZCrJ5zGMCAwRGE8\/DmgWBCAAAXcI1IrCmKmkt188IieceXXshiFgYqPyr6DPzvfPW7QIAhCAgLsEakVhzFTS\/A88LHNOao3VUJ\/HMCIwRGBivQQUggAEIACB5hKoFYWpZyoJAdNcnxW6dp+dX2jwGAcBCECghARqRWHMVNIpnWNy\/OkdNQn5PIYRgSECU\/MFoAAEIAABCGRDwERhTr1mrZz64VsiH3r4id7gvKQ4U0kImGz8Vsin+Oz8QgLHKAhAAAIlJ1ArCqN4Dn3rfTKnpbXmBnc+j2FEYIjAlPyjguZDAAIQKBaBOLkwcaeSEDDF8m2m1vjs\/ExB8jAIQAACEIhNoNoZSaaSXzzSLW+8\/KgsWP5sxXp9HsOIwBCBif1CURACEIAABLIhYM5IWrj229JywZUVH3rg\/rODfWF0f5ioCwGTjb8afsrOnTtlcHAwqGfDhg3S1dUVWef09LQMDAzIrl27gt+vWbNG+vv7S+f8hoFTAQQgAAEINI2ARmH0WrT22xWfocm8mtRbaVUSAqZp7kmv4r1790pfX59s3LgxqNT8f8mSJcc8ZHh4OPiZipYDBw7I6tWrZcWKFZGCx2fnp0efmiAAAQhAIG0CesCjJvTWisLoVJJeUSdW+zyGeTOFpNGX3bt3y9DQkLS0tIiKlLa2topRGLuj2YIm3AF9dn7aLxv1QQACEIBAugSe+9jZwRSSHvRY6ap2zIDPY5g3AiYsQqqJErsTmAiMRmM6Oo7dFMhn56f7mlEbBCAAAQikTcBEYc6681mZe0ZbxerN3jDhhF6fxzCvBIwdcdGIzMTERMXcFu0FKnI2b94sy5Ytm4ncVIrA9PT0SHt7e\/Dr1tbW4IsLAhCAAAQg0GwCP73mODnlyo9UjcKoDSah9+AZPTI5ORmYNTU1Jb29vbJ9+\/bIP9KbbXsz6y+1gDFgVcjs27cvUsQY9Wo7QcWMfnFBAAIQgAAEmk0gzsZ2aoNJ6P3vj3XInf9jfJZZCJhme6mB+uudQtJH2gnA4aRfI2BGRkZk8eLFRGAa8BG3QgACEIBAfQQ0CqN5MBqJqXZpQu\/U1KT8\/Kyjy6rHx8dldHSUCEx92LO5KzxllCSJV0WKlt+yZYssWLBglsE+zx9m4xmeAgEIQAACjRKIG4UJJ\/T6PIZ5M4VU7zJqsyfMokWLIvNlfHZ+oy8U90MAAhCAQDYE4i6pVmvshF6fxzBvBIw6rdJGdkakdHd3B0lM4Y3s4iTx+jh\/mM1rx1MgAAEIQCANAnGWVJvnmITef3ztKlm5ciVTSGk4wLU6fFavrvkCeyEAAQiUmUDcJdXKyCT0\/r\/f+oT89fXkwJSy3yBgSul2Gg0BCECgkATiJvOq8Sah998OTBGBKaQ3m2wUAqbJgKkeAhCAAARiE9BkXj3oUTe2q3WZhN6124+TlZ9gH5havLz7PQLGO5fSIAhAAALOEkiSzKuNNAm9v\/z9MTnrvGN3m3cWhIh4lcTbDEcgYJpBlTohAAEIQKBeAkmSeXUMO\/lH3bLkj0bkhDOvrveRhbwPAVPDLQiYQvZbjIIABCBQWgIH714nB+9ZK793z5s1Gfg8hiFgEDA1XwAKQAACEIBAsQjETeZFwBTLb5la47PzMwXJwyAAAQhAIDUC+9a+P6hr0dpvV63T5zGMCAwRmNReKCqCAAQgAIFsCMRN5kXAZOOPQj7FZ+cXEjhGQQACEIBALAJxknl9HsOIwBCBifWiUAgCEIAABIpFIM4BjwiYYvksU2t8dn6mIHkYBCAAAQikSiDONJLPYxgRGCIwqb5QVAYBCEAAAtkRqDWNhIDJzheFe5LPzi8cbAyCAAQgAIFEBGpNI\/k8hhGBIQKT6GWhMAQgAAEIFIeAnov0wtr3yxnXbZVTrvzIMYYhYIrjq8wt8dn5mcPkgRCAAAQgkDqBatNIPo9hRGCIwKT+MlEhBCAAAQhkR6DaNBICJjs\/FO5JPju\/cLAxCAIQgAAEEhPQaSQVMad++JZjppF8HsOIwBCBSfyycAMEIAABCBSLQKVpJARMsfyUqTU+Oz9TkDwMAhCAAASaRqDSNJLPYxgRGCIwTXuhqBgCEIAABLIhUGlTOwRMNvwL+RSfnV9I4BgFAQhAAAJ1EYiaRvJ5DPMqArNz504ZHBwMHL9hwwbp6uqK7AQHDhyQ1atXy549e4LfL1u2TIaGhqSlpeWY8j47v643hJsgAAEIQKCQBKKmkXwew7wRMHv37pW+vj7ZuHFj0LHM\/5csWTKro01PT8vAwIB0dnYGAsd8v2jRIunv70fAFPK1xCgIQAACEKhFwEwjnXXnszL3jLagOAKmFrUC\/F6jL7t3756JpAwPD0tbW1vFKIxtcvhe+3c+O78AbsMECEAAAhBIkcBPrzlu1q68Po9h3kRgVLDoZaIo4e+r9Q8ETIpvD1VBAAIQgEBuBHQaSS89WoAITG5uSPbgcMRFRcnExETktJBds8mHWbFiRWS0xqjXnp4eaW9vD25tbW0NvrggAAEIQAACRSJgppH2\/fUOmXvG78rU1JT09vbK9u3bpaOjo0imNmyLVxEYe8oojoAx+S9KsVYSr01axYx+cUEAAhCAAASKROCN\/ROiq5E2TpwhD7588oxpCJgieSlkS9IppDjixQ6\/jYyMyOLFi4nAFLgPYBoEIAABCEggYF7+7XPltT\/\/zzI+Pi6jo6NEYIrcMcIRl2pJvLVWHtnt9DkBqsj+xDYIQAACEKiPgL2c2ucxzJsppLjLqLU7qLjZt29fxWkjBEx9Lw13QQACEIBA\/gTs5dTfe+ZFWblyJRGY\/N1S3YJKG9mZiEt3d7ecd955szaxMzVeeOGFsmXLFlmwYMGsh\/isXovuT+yDAAQgAIH6CJjl1D888fcRMPUhdP8uBIz7PqQFEIAABMpGwBwr8Mx71iBgyuZ8014ETFk9T7shAAEIuEvA5MH87Ib\/g4Bx142NWY6AaYwfd0MAAhCAQPYE7P1g\/rpnkByY7F2Q\/xMRMPn7AAsgAAEIQCA5Ac2D+el7\/pP8x\/\/2AAImOT7370DAuO9DWgABCECgjAQ0D2bPL0+Unu++hoApYwdAwJTR67QZAhCAgPsETB7MBx4\/GwHjvjuTtwABk5wZd0AAAhCAQP4ETB7MX\/3jmfJfv3g3ZyHl75JsLUDAZMubp0EAAhCAQDoEpn\/4HXlh7fuDc5E+8vn7ETDpYHWnFgSMO77CUghAAAIQmE1AE3kfePlkuXj9NxEwZescLgmYyclJuffee+Wqq66S1tbWwrsKe5vrIvjC1xCgL5S3L2gezFfu\/QoCprldoJi1uyRgXLJVvY29ze3z8IWvIUBfKG9fOHj3Ojl4z1rRDe06OjqaCyLj2r05zLFZ3Fx68V2yFQHTrB77Vr30h+YydomvS7by2ZBuvzWJvHNueljaLr483cpzrg0BU8MB5sXv6emR9vb2nN1V\/fFTU1PS29srLtiqLcHe5nYn+MLXEKAvlLcvqO8v3nG1nHHdVjnlyo80F0TGtSNgagDXuWMVBePj4xm7hsdBAAIQgAAEGifw5Xc\/L0v++N8HIsanCwETw5sqYvSLCwIQgAAEIOAagXee8IZ300fqAwSMaz0ReyEAAQhAAAIQQMDQByAAAQhAAAIQcI8AERj3fIbFEIAABCAAgdITQMCUvgsAAAIQgAAEIOAeAQSMez7DYghAAAIQgEDpCSBgSt8FAAABCEAAAhBwjwACporPhoeHZfPmzUGJ7du3F2Yb5r1798qqVatk3759smzZMhkaGpKWlpaqvU835BsbG4tVNu1unMRem\/miRYtk69atsmTJkrRNqlpfEnt37twpg4ODQX0XXnihbNmyRRYsWJCZvUlsNUZNT0\/LwMCAdHZ2SldXV2a26oOS2GuzdYHvgQMHZPXq1bJnz57cPjPi8g2zNZ1gw4YNmfaJuPaG+44Lnw1mE9S8+m6mL3ZOD0PAVACvnU8HUx2QfvzjH8\/8P8vBKco0e\/BZvnx5rIHIvEhxxU6afTGJvfqhunv37hmRpd\/v2LEjU1GQxF67j2i\/0P6iojKOoEyDcRJb7eeZwSvrwSqpvcqzra0t0wHV5pTEXlNWB9b+\/v5AqPX19cnGjRszE+BJ7A33v3BfTqN\/1qojib1GHCpbPc+n6J8NRpjddtttgb15\/gFZyw8u\/x4BU8F7+uGpl74w5kXr7u7OPQoT\/mCs9WJoO3bt2iVXXnmlvPLKK5kNrgZrUnttd+QxCDRib9aDQD226kBwww03yKFDh2TFihWZioMk9hbhnUtir5Zdv369bNq0KdMIXLX3pdZng7k3LA6yGtCS8rUFYdE\/G8J\/jGl\/vvXWW+Xaa6\/NTNBm5cc8n4OAiaAfDhIsB5MAAATvSURBVLHnGXKv9ZdSrUFTf2\/+YrGjG1l1urB9tezNW8A0Yq8terPgW4+tauOll14qX\/va1zKfQkpib16Dqu23JPaGB6ws\/N\/oZ4O5Py\/bk\/CNisBk\/XmWxN4oAaPTtkX4IziPvtmsZyJgqggYu7PlHc42Zob\/qor7l1+eH1J27k1ce7W9WU\/J6DPr4Wum6LKel09qq7Lftm2bfPKTn5R169blImDi9gU7N8L0\/azz0JLw1fdrYmIiMDWvvLkk9uYdfannXTN\/SGpEec2aNUF0PMsrCd\/w6d\/m+6ynbbPkk8ezEDAImKb2uyQvvW2IDgif\/\/znM0\/irdde84F84403ZmZzElvtEHZra2us3Km0O0YSe7WszTL8fdq2RdWXxF6TV2REVtHttf8gMrl+Wef3JeEblVOStd1J7FW+dqL0ypUrgyn8PBLns3hX8noGAqaKgDGdzeUpJJfCxLateYgXI0LsD8YkU15Z95MkIW0t+9BDD83K6cr6wzSJveHXMmu2SftCpSmDLBnXw9dEjrKOZpSFbxEiXXmJiyyei4CpQNmeMipCQqExMzwFE\/6roFKnyWsKKam9eawusJkltde+N+u8jSS22svTbZuzDMUnsbeSgMkyhyCJveH3MI\/PjCT2Kt88RGG971oRBGJSvuG2Zr0qLQsBkfczEDAVPGD\/NeP6MmoTzsw66S38IVlr2XceYfdqf+nHsdeO1mQtvpIsQ7XbmdfAlcRe15bUh8VrkshdWoNAEr76TLMi7aabbsplZUwSe6OmkLKcrk36WWaLHd2jSxN4zRL7tPxNPcJp1NU6gYsb2VUKCecVgVG+1Tarsu2tFCXIOnkzrr1GGJqN7LJO4k3CtggCJqm9dg5BHmyT2mtvZOeCvXksRQ5\/3iZ510wirNbhAl+7\/+axB1cZBA4RmDJ4mTZCAAIQgAAEPCOAgPHMoTQHAhCAAAQgUAYCCJgyeJk2QgACEIAABDwjgIDxzKE0BwIQgAAEIFAGAgiYMniZNkIAAhCAAAQ8I4CA8cyhNAcCEIAABCBQBgIImDJ4mTZCAAIQgAAEPCOAgPHMoTQHAi4S0D1U9IDJW265RbI+k8dFXtgMAQiwkR19AAIQKAAB+6ymApiDCRCAgAMEiMA44CRMhIDvBHQX5iuuuEI6Ojp8byrtgwAEUiKAgEkJJNVAAALHEqh07pJ9VpCeFXPrrbfKtddeG5zJY+7ZtWtXUGFe28bjTwhAoNgEEDDF9g\/WQcB5AlGHXGrERa\/+\/v7grKxt27bJzTffHPxMD77Ta2hoSFTcZH1IpvPAaQAESkIAAVMSR9NMCORFIHxSc\/h7FSh6dXV1BWKmr69PNm7cmMsJyXkx4rkQgEByAgiY5My4AwIQSEjAjrjY00e64sjOfwn\/LuFjKA4BCJSIAAKmRM6mqRDIi4AtTO666y5pa2sLIi7h5dMImLw8xHMh4B4BBIx7PsNiCDhHwEwbLV++XO6\/\/\/6ZKaLw8mmmkJxzLQZDIDcCCJjc0PNgCJSLgOa6DA4OyrJly2Yl6CoFjcboZVYg6cojTfDVC1FTrn5CayEQlwACJi4pykEAAg0RUCGyatUquf766wPBomLFXj5tKmcZdUOYuRkCpSGAgCmNq2koBCAAAQhAwB8CCBh\/fElLIAABCEAAAqUhgIApjatpKAQgAAEIQMAfAv8f7\/FpdAApUjkAAAAASUVORK5CYII=","height":0,"width":0}}
%---
