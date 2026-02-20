%[text] # Two Parallel Charged Plates and a Point Particle
clear;
% The variables
syms a d v r real positive;
assumeAlso(d>0)
assumeAlso(r>=0)
assumeAlso(v>=0)
assumeAlso (v < 1)
assumeAlso(a>=0)
taylorOrder = 11;

vp = [v,0,0];
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
%[text] ## The Abnormality 

psinl = simplify(dot(r_vec_l1,r_vec_l1)/dot(r_vec_l1,vec_l1)) %[output:93ff2ce7]
psinr = simplify(dot(r_vec_r1,r_vec_r1)/dot(r_vec_r1,vec_r1)) %[output:201b8515]





%%
%[text] ## The Net Force of the Four Points

dforce = (psinl*vec_l1 - psinr*vec_r1)/(2*R^2); % net Force in the x direction first points.
dforce = dforce + (psinl*vec_l2 - psinr*vec_r2)/(2*R^2); % net Force in the x direction second points.
dforce = simplify(expand(dforce),'Steps',100) % Two charges %[output:2e9c3bb9]
dforce_x = dforce(1) %[output:37eff396]
%%
%[text] ## Check at v = 0, r =0
df_v_zero = simplify(subs(dforce,[v,r],[0,0]),'Steps',100) % Due to the four point charges %[output:8c92e749]
df_v_one = simplify(subs(dforce,[v,r],[0.9,0]),'Steps',100) % Due to the four point charges %[output:15dcd664]
%%
%[text] ## Integrate all the angles
%dforce_xT = expand(taylor(dforce_x,v,Order=taylorOrder));
dforce_xT = dforce_x;
dcircle = simplify(int(dforce_xT,a,0,pi),'Criterion','preferReal','Steps',100) % From zero to pi i.e., the full circle %[output:38239ed5]
%%
dcircle = simplify(dcircle,steps=100); % Net force of the two full circles with charge density
df_v_zero = simplify(subs(dcircle,[v,r],[0.5,0]),'Steps',100) % Net force of the full zero radius circle (point) %[output:51b3b2c0]
%[text] ## 
%%
%[text] ## Integrate for all values of r

dtotf = simplify(expand(r*dcircle),'Criterion','preferReal','Steps',100)/pi; % Add all the circles 
dtotf = taylor(dtotf,v,Order=taylorOrder);
totf = int(dtotf,r,0,r);
totf = simplify(totf,'Criterion','preferReal','Steps',100); % The force on a charge due a charge surface density

%%
%[text] ## Set the Plates Radius to Infinity

totfInf = simplify(limit(totf,r,Inf),'Criterion','preferReal','Steps',100) % Net force by plate charge density %[output:15368a71]
totfInfV0 = limit(totfInf,v,0) %[output:71384b41]

totfInfRel = simplify(totfInf/totfInfV0) %[output:86b47f7e]

%totfInfRel = taylor(totfInfRel,v,order=taylorOrder)



%[text] ## 
%%
%[text] ## Plots
figure(1)
plotrange = [0,0.99];
fplot(totfInfRel,plotrange) %[output:46ef278e]
hold on %[output:46ef278e]
relEq = sqrt(1-v^2) %[output:5d154b59]
relEqApx = taylor(relEq,v,Order=taylorOrder) %[output:63cfc48a]
fplot(relEq,plotrange) %[output:46ef278e]
fplot(relEqApx,plotrange) %[output:46ef278e]
xlabel('v/c')  %[output:46ef278e]
ylabel('ratio')  %[output:46ef278e]
legend({'Force Carriers',"Lorentz Mass","Lorentz Apx"},'Location','northeast') %[output:46ef278e]
hold off %[output:46ef278e]

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
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l1","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }}-v & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:34f9ec89]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_r1","value":"\\left(\\begin{array}{ccc}\n-v-\\frac{d}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & \\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:6900aafa]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_l2","value":"\\left(\\begin{array}{ccc}\n\\frac{d}{\\sqrt{r^2 +d^2 }}-v & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:123f7d6a]
%   data: {"dataType":"symbolic","outputData":{"name":"r_vec_r2","value":"\\left(\\begin{array}{ccc}\n-v-\\frac{d}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\cos \\left(a\\right)}{\\sqrt{r^2 +d^2 }} & -\\frac{r\\,\\sin \\left(a\\right)}{\\sqrt{r^2 +d^2 }}\n\\end{array}\\right)"}}
%---
%[output:93ff2ce7]
%   data: {"dataType":"symbolic","outputData":{"name":"psinl","value":"\\frac{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 -2\\,d\\,v\\,\\sqrt{r^2 +d^2 }}{d^2 +r^2 -d\\,v\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:201b8515]
%   data: {"dataType":"symbolic","outputData":{"name":"psinr","value":"\\frac{r^2 \\,v^2 +d^2 +r^2 +d^2 \\,v^2 +2\\,d\\,v\\,\\sqrt{r^2 +d^2 }}{d^2 +r^2 +d\\,v\\,\\sqrt{r^2 +d^2 }}"}}
%---
%[output:2e9c3bb9]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce","value":"\\left(\\begin{array}{ccc}\n\\frac{2\\,d}{{{\\left(r^2 +d^2 \\right)}}^{3\/2} }-\\frac{2\\,d\\,r^2 \\,v^2 }{{\\left(-r^2 +{\\left(-1+v^2 \\right)}\\,d^2 \\right)}\\,{{\\left(r^2 +d^2 \\right)}}^{3\/2} } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:37eff396]
%   data: {"dataType":"symbolic","outputData":{"name":"dforce_x","value":"\\frac{2\\,d}{{{\\left(r^2 +d^2 \\right)}}^{3\/2} }-\\frac{2\\,d\\,r^2 \\,v^2 }{{\\left(-r^2 +{\\left(-1+v^2 \\right)}\\,d^2 \\right)}\\,{{\\left(r^2 +d^2 \\right)}}^{3\/2} }"}}
%---
%[output:8c92e749]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\left(\\begin{array}{ccc}\n\\frac{2}{d^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:15dcd664]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_one","value":"\\left(\\begin{array}{ccc}\n\\frac{2}{d^2 } & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:38239ed5]
%   data: {"dataType":"symbolic","outputData":{"name":"dcircle","value":"\\frac{2\\,\\pi \\,d}{{{\\left(r^2 +d^2 \\right)}}^{3\/2} }-\\frac{2\\,\\pi \\,d\\,r^2 \\,v^2 }{{\\left(-r^2 +{\\left(-1+v^2 \\right)}\\,d^2 \\right)}\\,{{\\left(r^2 +d^2 \\right)}}^{3\/2} }"}}
%---
%[output:51b3b2c0]
%   data: {"dataType":"symbolic","outputData":{"name":"df_v_zero","value":"\\frac{2\\,\\pi }{d^2 }"}}
%---
%[output:15368a71]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInf","value":"2+\\frac{4\\,v^2 }{3}+\\frac{4\\,v^4 }{15}+\\frac{4\\,v^6 }{35}+\\frac{4\\,v^8 }{63}+\\frac{4\\,v^{10} }{99}"}}
%---
%[output:71384b41]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInfV0","value":"2"}}
%---
%[output:86b47f7e]
%   data: {"dataType":"symbolic","outputData":{"name":"totfInfRel","value":"1+\\frac{2\\,v^2 }{3}+\\frac{2\\,v^4 }{15}+\\frac{2\\,v^6 }{35}+\\frac{2\\,v^8 }{63}+\\frac{2\\,v^{10} }{99}"}}
%---
%[output:5d154b59]
%   data: {"dataType":"symbolic","outputData":{"name":"relEq","value":"\\sqrt{1-v^2 }"}}
%---
%[output:63cfc48a]
%   data: {"dataType":"symbolic","outputData":{"name":"relEqApx","value":"1-\\frac{v^2 }{2}-\\frac{v^4 }{8}-\\frac{v^6 }{16}-\\frac{5\\,v^8 }{128}-\\frac{7\\,v^{10} }{256}"}}
%---
%[output:46ef278e]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQ9sXcd1p09VOTEdxyvLTltLTEBvTbvdIpXjICXLBLC7aGBtAWqBOClFLlBXFQoBSVxma\/OPrDay0lgUSQsbNm4BrVdQBaSUlMRGEgW76xiLJLXNNeuosXY3SG3FW0YmFTeOSXUdhXbsKotz1aGHl\/e9++5799\/M\/S7wYIucO\/fMd+a++fHMmZmf+9nPfvYz4YIABCAAAQhAAAIOEfg5BIxD3sJUCEAAAhCAAAQCAggYOgIEIAABCEAAAs4RQMA45zIMhgAEIAABCEAAAUMfgAAEIAABCEDAOQIIGOdchsEQgAAEIAABCCBg6AMQgAAEIAABCDhHAAHjnMswGAIQgAAEIAABBAx9AAIQgAAEIAAB5wggYJxzGQZDAAIQgAAEIICAoQ9AAAIQgAAEIOAcAQSMcy7DYAhAAAIQgAAEEDD0AQhAAAIQgAAEnCOAgHHOZRgMAQhAAAIQgAAChj4AAQhAAAIQgIBzBBAwzrkMgyEAAQhAAAIQQMDQByAAAQhAAAIQcI4AAsY5l2EwBCAAAQhAAAIIGPoABCAAAQhAAALOEUDAOOcyDIYABCAAAQhAAAFDH4AABCAAAQhAwDkCCBjnXIbBEIAABCAAAQggYOgDEIAABCAAAQg4RwAB45zLMBgCEIAABCAAAQRMA31gfn5e9MMFAQhAAAIQcI3AxSuulZ5fv8E1s2PtRcDEIFLhMjQ0JLOzs7EwKQABCEAAAhAoE4Gfvuv98uqvbJNn\/uQ3pb29vUymtWwLAiYG4VNPPSUDAwMyOTkpmzdvbhl4lhWoyJqamnLCVuWAvVn2BvhmS9ctvrxr2faGsvIdf3ROHv\/Bz8uVT0zK8f\/y59Ld3Z0tiJxrR8A0KGCmp6dL73wjtlywVbFjb7ZvO3zhawjQF6rXF8Yf\/Qc59vSLsue9r8nIx\/5AXBkXkngKAYOASdJfUi3Ll2qqONdUBl\/4ImCy7QNl5Ht28VU59vQPAvHylY++R84990wwi4CAyacvlOopLg0Cmq\/z8MMPyx133OHEXCf2ZtvV4QtfQ4C+UI2+oOLlY8e+GzT25MfeE\/zXpTEsqZeIwHgUgUnqfMpDAAIQgIAfBFS8bPvLb0v\/+35JRm6\/fqVRCBg\/\/NtUK3x2flNAuAkCEIAABEpFoJZ4IQJTKjflbwwCJn\/mPBECEIAABBojYKaNPnDDhlWRF3O3z2OYl1NI4+Pj0tHRIX19fTV7wIkTJ2T37t3B73t7e+XAgQPS1ta2przPzm\/s9aAUBCAAAQiUkYAm6upqo\/C0kW2rz2OYdwJGxcuhQ4dkbGyspoBRh2q5w4cPB6JldHRUNm3aJCMjIwiYMr6l2AQBCEAAAqsImGXSD27\/VdHoS60LAeNAx1lcXJSdO3fKhg2XHLl169aaAkajLzMzMytRl\/C\/q6JeHXArJkIAAnUIcMxJNbvHk8+fl7\/+zuvypXs+KO\/aeHldCAgYB\/qICpilpaUgkqIRlZ6enkQRmFrlfXa+A27FRAhAoAYBjjmpdtfo6uoKdl2POx7A5zHMuymk5eXlWAGj3f7MmTOyY8cOOXfuXN0NfozzBwcHRTuMXtph4jpNtV8tWg8BCGRNwKVjTrJmUbX6zdEFtTansyNzCwsLwXl+bGTnQC9pRMDolNHx48eDHJiNGzcG+TB61cuBsZuuYkY\/XBCAAASKIuDzX9ZFMXXluXG+1zPx9GNfCBgHvBsnYKJ+r9GY4eFhmZiYkM7OzlWtjPorhwiMAx0BEyHgOYG4Qczz5le6eXG+tyMwcdEal0FWbgqpWQHjo3p1ueNiOwSqTiBuECs7H7Pw4vTp06tMrbetRZZtMjzNM+qtZG3UDjPe9Pf3p3oYcBLfJynbaLvKUq5yAkbBR00haS5M1F4wPju\/LJ0QOyAAgeQEXP9uMgJGp+67u7uTA0jxDmV5zz33yJEjR4IofFwkP8VHN1VVEt8nKduUMQXeVAkBE6WCzX4xyp6N7ArsgTwaAhBoioDrA1OcgDHf2ydPngz4mCi4mfI30DSXUS\/dRkOjOVu2bFnJb7QXa9g\/t4HXEit67\/79++XgwYNBrqS9+amudlWxo+kEuur1\/Pnz8vjjjwc\/e+SRR9b820Rgatlj9iZ76aWX5L3vfa988pOflE996lNi2h6OBiXxfZKyTXXEAm\/yTsCkzdJn56fNivogAIH8CLj+3RQnYPSPTBMZV2FiIiRKWFeQPvDAA0HkJvwHqgqNubk5+cM\/\/MNA1JgIj12fvet6vRxI401781N74ccf\/dEfrdoI1dhiNka1bbvxxhtr2mO3TyNApg1qu3Lat2+f7N27NxBSeiXxfZKy+fXedJ6EgInh6LPz0+lC1AIBCBRBIPzdpGfi6Kesl264Zm+6FpUDY6IkZod0E7mwoyS33HLLqkUX4UhJLdFRq1ytn9fjaASGETBmH7FwNMcWMFqf2QFehYj93Oeee27V72wBE2VHknEpSdmy9p1adiFgEDCu9VnshQAEIv4K163lxx+dKy2bkds7Vh02WC8CE\/U7c8adChh7aiccHbEFzMDAwCoeZurHXm3aSAQmPJ2lle7atUuMgAkLrah\/6z217Hn55Zfl2LFjq\/Iw7TSH8CKSJKIkSdnSdp4ahiFgEDCu9VnshQAEIgSMqxGYqCTe8LRQOAJjC5h6EZiwKIjqOLVyYOypm8cee2zN8TM6TZVUwNSyR0VGrd+pHXfffbfce++9K9t8JBElScq69mIhYBAwrvVZ7IUABBLmQZQRWCs5MLaACddjzrbTvb0+\/vGPr+TAhFef2kxqrUIyuSz2eXkqeDS3RvNvkgiYcA6MbY9OIdkCxkSb+vr6yIGp03kRMAiYMn63YRMEIOD5d1OcgKm3CskWMIqp1uoe++dR00dhEWNP8egUkdmd3c7X0TwdTRB+9NFHV1YLNTKFpIKnlj3hCEw4P4gppOiXAQHj+ZcEowAEIOAnAZ+nBvz0WHqtSuL7JGXTszCfmhAwCJh8ehpPgQAEUiXg88CUKigPK0vi+yRlXUOFgEHAuNZnsRcCEPAgBwYnNk8giShJUrZ5i4q5EwGDgCmm5\/FUCECgJQI+D0wtganAzUl8n6Ssa+gQMAgY1\/os9kIAAkRgKt0HkoiSJGVdg4qAQcC41mexFwIQ8EDAxK1CytLJSXffNSuiTp06tXLgo7HPnJGUxunVjbY5iShJUrbR55elHAIGAVOWvogdEIBAAgKuD0wuChg9tPF3fud3RPdn0UuFzf333y\/f+c53ZPv27Ss\/T+DGpoom8X2Ssk0ZU+BNCBgETIHdj0dDAALNEnB9YIoTMKZ9yqe3tzfYZl8v+\/Rn3R9F92XRn0WdWq37xVx11VWrfmc2lDMnV2vdn\/70p1e5IRxNMRGYm2++WZ5\/\/nnZs2eP6HlNGsk5evRocO+73\/3uQMCE93Cx95Oxjwdo5Oe1+kYS3ycp22xfLOo+BAwCpqi+x3MhAIEWCLg+MNUTMGbDNz1x2ggU3YgufPqz4lNRoJduOmfvqKs\/01Or77rrrkBY2KdRz8\/PrzpPybghatt+E2lRkXT77bfLF77whZVt\/XX6SC89VqCjo0O2bdsWiClzuKN9zpJ93pHeo5GbO++8U2r93D6vKdxNkvg+SdkWumMhtyJgEDCFdDweCgEItEYgPDC98dKcvP7D8h7meNkvdMj6d3SsNLqegAnvTKv\/VgHy4IMPysTExIpACNdhn6F0zTXXrDq12q4zSsDUs8eu95vf\/KbceuutgbAyIuSRRx4JBIyZWooSRFEHNmq5cFsb6RVJREmSso08u0xlEDAImDL1R2yBAAQaJBAemJY+v0+WvnBfg3fnX+zqj9wnV\/\/u3oYEjH32kJmq0ekgnepRAWO27g9P15jKdQoo6tRqc95QWMDUOtDR1GcLGP2ZipgPfehDK1Gchx56aEXA2FNfWtY+wsAk\/OrP7WmqWj+v5aUkoiRJ2fx7RWtPRMAgYFrrQdwNAQgUQqDKERhbwIRPajbOCK80qhWBUYFkT\/tEOdMWMJpDs2\/fvuBk6AsXLgRTV+bwxQ9+8IPBQY\/mhO24KSnTjiihpGcnIWDqv1oIGARMIV++PBQCEGiNgOt\/WbeSA2MP\/HYOjJ07o1NI9qGPtQSMRk\/0Mgc3xgkYFRYmGdccslhLwGhk5bOf\/Wyw9Prv\/u7vglwZfY5ZvaQ5MLV+Tg5M\/PuBgEHAxPcSSkAAAqUj4IuA0dVA9mVEQb1VSLaACZ9abaZm6kVg9B6NlFy8eFGuu+46+drXvrbKBnuFkP7CjsCogFFhcvz4cTl8+LBs3LhxJQKjOTD2dNCf\/umfyjPPPBNMeYVXSxk7a9lfr8Ml8X2SsqXr5DEGIWAQMK71WeyFAAQ82MgOJzZPIIkoSVK2eYuKudNLAWPCeeGMcBuxre5VGRslHXaDz84vpsvxVAhAIA0CfDelQdHNOpL4PklZ12h4J2DM3GS9bZ3teVITDpyZmQk2StKEriihY8KarjkYeyEAAT8J+Dww+emx9FqVxPdJyqZnYT41eSNgTELYhg0bAnJbt26tua2zzlGaZKo4zD47P67t\/B4CECgvAb6byuubrC1L4vskZbO2O+36vRIwS0tLwZr7ekvi4tb7M4WUdhejPghAIAsCPg9MWfDyqc4kvk9S1jVG3ggYAz5OoJjf65bQunzOnIcRlwMzODgoXV1dwWPa29uDDxcEIACBogj4PDAVxdSV58b5Xjfq049eCwsLMjQ0JD6mQVRWwJw9e3bVErhz587VzYGxO7aKGf1wQQACECiKQNwgVpRdjT437jDHRutpplx4iXWSOpT7PffcE+ztUm+vliR1Ji0b5\/upqSnRj30hYJJSLqB8oxEYc9iWmmgfuBXukKajTE5OyubNm4nAFOBTHgkBCKwlEDeIlZ2ZiwLGbECnbM3p00VwjvO9HYGZnZ0NxAwCpghPJXxmnIDR6sLLrOup8biOktA8ikMAAhBIhYDr301xAqbeRnbnz5+Xxx9\/PBiUwxvEmYHafK9fddVVcvLkyYC5\/k6PAtBN7Ez6QG9vb3DGkn3VWsWqdR49elQ+8pGPyGc+8xk5ePBgsJGdGXduvvnmILKvEX1Thz3ehDfEa7YjJPF9krLN2lPUfZWbQlLQ6lDtVPYuivrzqK2kfXZ+UZ2O50IAAq0TcP27qZWjBHSxhvm+to8SsKd3lPCOHTvkrrvuClakajmTKhB1GrWWr3V2kfGWrmDVy9Snp1LrVhxGmJjUBD15enh4ODh4UvMlzanV+nM9CLLesQWN9Iwkvk9StpFnl6lMJQRMlOqNUvfhPWCM2BkYGPAy\/FamjogtEIBAMgLhgeniT+bln39yKXGzjNfPX9Eu6654c\/FDPQGjbTMnR+v3svmj88EHHwxEgUkBCNdhf9frWUhGRGhqgF1nlICJiwjZ5xeZ+owYiYr825EXjdxolGfdunUrUZtWfJRElCQp24pNRdzrnYBJG6LPzk+bFfVBAAL5EQh\/Ny0\/+xlZfnZ14mZ+1sQ\/qe2mQWm76RMrBesJBo102JuLmukgFQEqYOzTqM10kG2BTt\/ccsstDR3maE8B2bmR4RaZDVA1imMus4u7OdE6fMhkR0fHyn5kjewQH0\/xUokk41KSso0+vyzlEDAxnvDZ+WXphNgBAQgkJ1DlCIwtYO6++265995716wIqneYox2BMeKjnnhR79hTVcZbRpRs27Zt1f5j4YiM+uqrX\/1qcJueQN3q6qUk41KSssl7YbF3IGAQMMX2QJ4OAQg0RcD1gamVHJhwpEMBal6JfUyMTiHt379\/Zcqm1hSS7gdm7q\/liFq2mkjRJz\/5SfnUpz4V3K5H0qhAMtNXase+fftk7969ojkwmgS8Z8+eNcfWJOkESXyfpGwSG8pQFgGDgClDP8QGCEAgIQHXByYjCnQ1kH2ZVUT1ViHZAsZEO8xKI7P6p14ERu\/RqaeLFy\/KddddJ1\/72tdW2bBr165VibbhhR+msBFMmqT7yCOPBDvBHzp0KPi1vULKtjfJUTa1ukQS3ycpm7ALFl4cAYOAKbwTYgAEIJCcgM8DU3Iaxd6R1vLoRluRxPdJyjb6\/LKUQ8AgYMrSF7EDAhBIQMDngSkBhlIURcAU4wYEDAKmmJ7HUyEAgZYIIGBawuf0zUl8n6Ssa1AQMAgY1\/os9kIAAtZSWvugWcBUg0CSAxoRMNXoE5Gt9Nn5FXYrTYeA8wR0pYueMqxn3XBVj0BXV5foGX2602+9y+cxjAgMEZjqvfm0GAKeELAP7fOkSSvNeGHpVZn+2xflyefPy7s2Xi4Pbv8V35rYUntUuMSJF30AAqYlzG7f7LPz3fYM1kMAAj4SOLv4qhx7+gcy\/uhcIFxGbr9e+t\/3Sz42NZc2+TyGEYEhApPLS8RDIAABCMQROPb0izL+6D8ExVS0qHjhao0AAqY1fk7f7bPznXYMxkMAAt4Q0KjLx459V3TaCOGSrlt9HsOIwBCBSfdtoTYIQAACDRKwhcv7f3lDEHHRaSOu9AggYNJj6VxNPjvfOWdgMAQg4AUBFS46VaRTRibignDJxrU+j2FEYIjAZPPWUCsEIACBEAE7QVeFS\/\/7rpMP3LABThkSQMBkCLfsVfvs\/LKzxz4IQMAPAuGVRQ9u\/1WES06u9XkMIwJDBCan14jHQAACVSOAcCne4wiY4n1QmAU+O78wqDwYAhDwmgB7uZTHvT6PYURgiMCU503DEghAwGkCCJfyuQ8BUz6f5GaRz87PDSIPggAEvCaAcCmve30ew4jAEIEp75uHZRCAQKkJIFxK7Z7AOARM+X20ysLx8XHp6OiQvr6+WMvPnDkjw8PDMjExIZ2dnWvK++z8WDgUgAAEIBBBAOHiTrfweQzzLgKj4uXQoUMyNjYWK2CWl5dldHRUTp06JUeOHEHAuPNOYikEIFAAAXsDuksnRLMcugA3JHokAiYRrmIKLy4uys6dO2XDhkubIm3dujVWwKhjVfDoRQSmGL\/xVAhAoPwEnvje+WDn3CefPy9my382oCu\/35hCcsNHogJmaWlJNm3aFERVenp66goYLb9v3z7p7+8PREycgBkcHJSurq6ARnt7e\/DhggAEIOArAY22qGCZ\/tsfBIcsqnBh51w3vD0\/Py\/60WthYUGGhoZkenpauru73WhAg1Z6N4VkpoXiBMyJEycCRLfccktDOTA2TxUz+uGCAAQg4BsBk9+i5xTpZbb856widzw9NTUl+rEvBIwD\/mtEwGji7tGjR2XPnj2BSm0kiXdyclI2b95MBMaBPoCJEIBAcgIIl+TMynqHHYGZnZ0NxAwCpqzesuxqRMDolNGtt94ahNNYheSAUzERAhDIjEB4RZE5HTqzB1JxrgRI4s0Vd2sPixMwJtn39OnTax4UpVB9dn5rpLkbAhBwmYBOEWl+i+a56PTQyO3XB9NFXH4R8HkMq2wOjOmiRGD8ellpDQQgUJsAibnV6x0IGId8HhWBMT\/TFUfhLGwEjEPOxVQIQKApAlHTRLqiiMTcpnA6dRMCxil3pWusz85PlxS1QQACZSMQ3niO\/JayeSh7e3wew7ybQkq7O\/js\/LRZUR8EIFA8AXuayGw8N\/Ab15HfUrxrCrHA5zEMARPTpXx2fiFvEw+FAAQyIRA1TfT+X75a2DE3E9zOVOrzGIaAQcA48yJiKAQgsJpAONqiOS1ME9FLbAIImAr3B5+dX2G30nQIOE0gvOkc2\/w77c5Mjfd5DCMCQwQm05eHyiEAgXQI1Iq2sJooHb6+1oKA8dWzDbTLZ+c30HyKQAACBRMI57YQbSnYIY493ucxjAgMERjHXkfMhYD\/BMht8d\/HebUQAZMX6RI+x2fnlxA3JkGg0gSItlTa\/Zk03ucxjAgMEZhMXhoqhQAEGiNgRMsT3zu\/ci4RK4kaY0epeAIImHhG3pbw2fneOo2GQaDkBKKmiMhtKbnTHDXP5zGMCAwRGEdfS8yGgHsEoqaIPnDD1eyS654rnbEYAeOMq9I31Gfnp0+LGiEAgTABpojoE0US8HkMIwJDBKbId4tnQ8BLArVWEbG1v5fuLnWjEDCldk+2xvns\/GzJUTsEqkWAvJZq+duV1vo8hhGBIQLjynuInRAoHQEjWp743pIce\/pF0bOINBmXvJbSuaqyBiFgKut6EZ+dX2G30nQINE0gSrS88+rLZeA3riMZt2mq3JgVAZ\/HMCIwRGCyem+oFwLeEEC0eOPKyjUEAVM5l7\/ZYJ+dX2G30nQIxBIIixa9QaeI2GQuFh0FSkTA5zGMCAwRmBK9apgCgWIJ1BMtnPpcrG94enMEEDDNcfPiLp+d74WDaAQEWiRQb3pIE3I16sIFAVcJ+DyGeRmBGR8fl46ODunr64vsc4uLi7Jz5045ffp08Pve3l45cOCAtLW1rSnvs\/NdfSGxGwKtEqi15JnVQ62S5f6yEfB5DPNOwKh4OXTokIyNjUUKmOXlZRkdHZWenp7g9+bfmzZtkpGREQRM2d4+7IFASgSidsRlyXNKcKmmtAQQMKV1zZuGmajKhg0bgh9u3bq1ZgQm3JwTJ07IzMxMZBTGZ+c74FZMhEDTBFSwGNHy5PPng\/83+7RoPssHbrj0XcEFAZ8J+DyGeROBUQGztLQkGkmxIyyNdEwETCOUKAOB8hOol4TLNv7l9x8Wpk8AAZM+08xqDE8RxT3IRG62b98eGbExzh8cHJSurq6guvb29uDDBQEIFEtABYtex57+gTzxvfOikRaNspiN5UjCLdY\/PL0YAvPz86IfvRYWFmRoaEimp6elu7u7GIMyeqo3ERjDJ4mAMWX13rgkXpu\/ihn9cEEAAvkTqBVlIZ8lf1\/wxHISmJqaEv3YFwKmnL5aZVWjAqYR8aIVmwjM5OSkbN68mQiMA30AE\/0iYHJZnnx+aU2URfNYmBryy9+0pnUCdgRmdnY2EDMImNa5Zl5DIwImbuWRbaTP84eZO4MHQKBJAkRZmgTHbRAIEfB5DKvkFJIutT537lzNaSMEDN8BEMiXgC1YzP\/bJztrTgurhvL1CU\/zgwACxiE\/RkVgzM\/6+\/vlxhtvXLWJnWnali1b5PDhw7Jx48ZVrfXZ+Q65FVM9IxA1LaRNNGcNMS3kmcNpTmEEfB7DvIvApN1LfHZ+2qyoDwK1CNQTLCTf0m8gkB0Bn8cwBExMv\/HZ+dm9MtRcdQKNCBamhareS2h\/HgR8HsMQMAiYPN4hnuE5AQSL5w6mec4SQMA467rWDffZ+a3ToYaqEjCJtmcXl1eWNpscFjMlRISlqr2DdpeJgM9jGBEYIjBletewpaQEolYJIVhK6izMgoBFAAFT4e7gs\/Mr7FaaXoeA2Z5ft+V\/4ntLcuzpF4PSukJIL5Ju6T4QcIeAz2MYERgiMO68iViaCYF60RUzDcSy5kzQUykEMieAgMkccXkf4LPzy0sdy7IiYJJtX1h6dVV0hemgrIhTLwSKJeDzGEYEhghMsW8XT8+MQHgqyERajFghupIZeiqGQGkIIGBK44r8DfHZ+fnT5IlZEYgTK\/pcVgdlRZ96IVBeAj6PYURgiMCU983DskgC4WmgcGQFsULHgQAEDAEETIX7gs\/Or7BbnWk6YsUZV2EoBEpJwOcxrFQRmBMnTsju3btXdYKxsTHp6+srrGP47PzCoPLgNQTqTQFpYV3CrDkr+t8P3HB18P+czkxHggAE4gj4PIaVRsCoeDl+\/PiqE6EXFxeDk6O3b99emIjx2flxHZ\/fZ0PAnvLRfVbsfxuxovkqKlZ0+bL+1+zBko1F1AoBCPhKwOcxrBQCxgiVkZER6e7uXtWPFP74+PgqYZNnR\/PZ+XlyrNqz7IiKtr2WUCGqUrWeQXshkC8Bn8cwBExMX\/LZ+fm+Rn4+DaHip19pFQR8IeDzGFYKAaMdhSkkX14XP9sRlUyrm8EZAWOmfoio+Ol\/WgUBVwkgYHLyHEm8OYHmMWsIGCGi\/zW71IajK0ak6H81R0UvEmrpTBCAQJkJIGDK7J2MbfPZ+RmjK131dhRFjTN5KVGRlCiRQjJt6VyKQRCAQAwBn8ew0kwhlbUX+uz8sjJv1q5a+Shan56sbF\/2ycp2JAWR0ix97oMABMpIwOcxrDABY688uvHGG4Pl0qdPn470\/5YtW1iFVMY3I0ebwlM8ZxeXg\/wTM+Vj56KoWVECRX\/O\/ik5Oo1HQQAChRNAwBTuguIM8Nn5xVF988lhYaK\/0akdvWqJEyNQTMKsiaAgUMrgUWyAAATKRMDnMaywCIzt4CL3gdE9Zjo6OmpulOez87N6yWxRos+4lGNyKWLSiDAxQsREUTRR1v4Zm7pl5TnqhQAEfCPg8xhWaQGj4uXQoUNS77gCn53fyItqT82YiMglEbJakLwpVC6JlPBlREetqMmlKMqllT1cEIAABCCQDgGfx7BCBUzUsukol+3atUt0l960LhPx2bDh0oC5detWbyMwYQFihIZhGSVE4sSI\/j5KkFxKgG0LqlahYsoRMUmr51IPBCAAgWQEEDDJeCUuXW8KKXFlDdygz1taWpJNmzbJ6Oio9PT0xAqYwcFB6erqCmpvb28PPs1e4YTTN8XE6uiFTr3Yl4oNvcL32\/8299R6hl2fLSxswfGm8ECMNOtj7oMABCBQFIH5+XnRj14LCwsyNDQk09PTa47qKcq+tJ5baAQmrUY0W8\/y8nJDAubKv++\/1BFevXZFwDTyTCM4tOz8v9wbdZ+pN+p34fsWXr1mVbF1V7TLu\/4l2hH8\/0YT+bgkPvQy4uRNYXKpDBcEIAABCPhHYGpqSvRjXwgYz\/zcqIB59qv98uE7PtxU6\/\/5J5dUcL3r4vLaMhcbuC+uXhU09rWu7dK\/f976+borNq8UMeVNObtsuK64Z\/N7CEAAAhAohoAdgZmdnQ3EDAImQ1+cOXNGduzYIefOnVvzlKz2gWlUwAwMDJTC+VGixhZIthAKl734k4UVrlE3ThL3AAAeX0lEQVT3NCqYbCETFkRGDIWFkAomBFCGLw9VQwACEKhBgByYjLuGLSS2bdsWTOv09\/eL2eBOE3i7u7tTt8I1AZM6gBoV2mLGiB0jjuzfqSgKi6E4IRQtblZHgWxhhPDJy+s8BwIQ8JEAAiZjr4aTeO29WRT+sWPH5MCBA9LW9mZeRxomIWDSoBhdR1gEhQWQiQjZAqme+LGFz6WIziXRoz9H8GTnR2qGAATcJoCAydh\/YQGjy6vn5uaCpdMKXwXN4cOHZePGjalagoBJFWcqlUUJH\/MzO+KjgqiW4IkTO0xppeIqKoEABBwggIDJwUkqUvQKi5bHHntMZmZmMonANNIsn53fSPtdKGOEjEZzbGFjBE+c2NEITlRUB6HjgvexEQIQqEfA5zGsNMuow9EQs0uu7tVy5MgR6ezsLKSX+uz8QoAW+NBmhY6J6Fx2TfeaqavLrk0\/N6tARDwaAhDwjIDPY1gpBIwRL5q4m0Wybiv90Wfnt8LF53trCZ3XfzS7KsJjMzC5OHYkZ\/013f8S2Wl+00OfOdM2CEAgewI+j2GlEDB578SbpMv47PwkHCi7moAtct54+angl2bKyvy7EYFDBIeeBQEIZEnA5zGsFAJGnaeJu0XmutTqQD47P8uXhrpV0MwHy8xNDk4jAueyay8dV6HRG70QOPQkCECgFQI+j2GlEDAmAnP69OlIP2W1kV0jncJn5zfSfspkR8AIHBOx0SkqvcIRHDM9hbjJzhfUDAFfCfg8hpVCwJS54\/js\/DJzr7ptZorq9ZefCiI5taI3UeKGqE3Vew\/th8CbBHwewxAwMT3dZ+fzkrtJoBFxE1459abQYdWUm17Hagg0R8DnMQwBg4Bp7q3grlISiJqWqjclhbAppRsxCgKpEUDApIbSvYp8dr573sDiZgmEozaab4OwaZYm90HAHQI+j2FEYIjAuPMmYmnqBBoVNmYTP10dRY5N6m6gQghkRgABkxna8lfss\/PLTx8LiyJghM1rL3wxMCEcsTE5Nm995x3B7xE2RXmK50KgPgGfxzAiMERgeP8h0DCBuBwbFTZEaxrGSUEIZE4AAZM54vI+wGfnl5c6lrlGQIVNvWiNHpip+9gQqXHNs9jrOgGfxzAiMERgXH8\/sb+kBOxoTdQUFKKmpI7DLK8IIGC8cmeyxvjs\/GQkKA2B1gnEiRqmn1pnTA0QsAn4PIYRgSECw9sOgUIJxIkaEoULdQ8Pd5wAAsZxB7Zivs\/Ob4UL90IgSwIqaswxCvb0k736iXyaLD1A3b4Q8HkMIwJDBMaX95R2eE7AThRefnZqpbWsfPLc8TSvJQIImJbwuX2zz8532zNYX3UCtaae7AMuidJUvZfQfp\/HMCIwRGB4wyHgDQETpdHTu82ybm2cihrNpUHQeONqGtIgAQRMg6CKLnbixAnZvXt3YMbY2Jj09fXVNMku29vbKwcOHJC2trY15X12ftH+4vkQyJqAPe104dRBWXfl+uCRRtDo\/7fd9ImszaB+CBRGwOcxzJsIzJkzZ2R4eFgmJiaCjmL+v7OzM1KUjI+Py+HDhwPRMjo6Kps2bZKRkREETGGvGQ+GQPYE6iUHm9VOCJrs\/cAT8iOAgMmPddNP0ojKzMzMSiRFBUpHR0dkFCZcNvxv2wifnd80bG6EgCcE6uXRIGg8cXLFm+HzGOZNBEYFi14mihL+d1iUhCMwPT09kWLHZ+dX\/L2m+RBYQyBO0JBDQ6dxjYDPY5hXAsaOuGhUZW5uLnJaSDugTjnt2LFDzp07J9PT09Ld3R3ZL43zBwcHpaurKyjT3t4efLggAAG\/CdiC5rUXHhZzSjdJwX773fXWzc\/Pi370WlhYkKGhobrjnKvtraSAUXFz\/PjxIAdm48aNEhetGRgYWOVfFTP64YIABKpFIGovGnvZ9lvf+eEgQZgLAkUSmJqaEv3YV70\/1Iu0tZVneyVgFETcFNLy8nKQtGtPGdkJwOGkXxOBmZyclM2bNxOBaaW3cS8EPCNgBE14t2CWbHvmaMeaY0dgZmdnAzGDgCmxE8NTRrWSeJsVMD46v8TuxDQIOEnACJqo6SaiM0661HmjyYFxwIVJllFHTSFpLkzUXjA+O98Bt2IiBJwlYMTMGy99X15f\/FLQDnJnnHWns4b7PIZ5M4WkvavWRnYm6tLf37+SrKsRmkOHDgWdko3snH03MRwCzhAwgmb5uydE1r0Y2H1pVVMXm+k540X3DEXAuOez1Cz22fmpQaIiCEAgEQFyZxLhonALBHwew7yKwLTg45q3+uz8LHhRJwQgkIzAGy\/Nyc9+7kV54+WnxJyyzVRTMoaUrk3A5zEMARPT8312Pi89BCBQPgLhRGDObSqfj1yyyOcxDAGDgHHpXcRWCFSKQNRUE3kzleoCLTcWAdMyQncr8Nn57noFyyFQPQIrkZnvPykXX\/1WAMBEZziAsnr9odEW+zyGEYEhAtPoe0A5CECgJARUzFz41udE1v3jmiXaiJmSOKkkZiBgSuKIIszw2flF8OSZEIBAugQ0CXj5u8cRM+li9aY2n8cwIjBEYLx5UWkIBKpOwIgZPdpA1j\/DNFPVO4SIIGAq3Al8dn6F3UrTIeA9ASNm\/vnVUys5M5oA\/NZ33SF6rAFXNQj4PIYRgSECU423mFZCoMIETALwhaf\/WtZd9SOONKhQX0DAVMjZ4ab67PwKu5WmQ6CyBMJHGpiVTBw26WeX8HkMIwJDBMbPt5ZWQQACsQRWIjOnDsq6K9ezLDuWmHsFEDDu+Sw1i312fmqQqAgCEHCewGvff0J++sLDq5Zlt900SL6M4571eQwjAkMExvHXE\/MhAIE0CWjy72tnvyh6aradL8MUU5qU86sLAZMf69I9yWfnlw42BkEAAqUiwBRTqdzRlDE+j2FEYIjANPVScBMEIFAtAv\/0Xz8tuiTb7C+j00uXzmXqrhYIx1qLgHHMYWma67Pz0+REXRCAQDUI6BTT64tfklf+ZlzW\/+LlJP6W3O0+j2FEYIjAlPz1wzwIQKCsBKKiMuTKlMtbCJhy+SNXa3x2fq4geRgEIOAtAROVucBy7NL52OcxjAgMEZjSvXAYBAEIuElAhcyFpz8nr73wsKy\/7vzK9BJRmeL8iYApjn3hT\/bZ+YXDxQAIQMBbAlFRGfaVyd\/dPo9hRGCIwOT\/RvFECECgMgRMVGb52Sl5yy9fSdJvzp5HwOQMvNnHnThxQnbv3h3cPjY2Jn19fTWrMk7VAlu2bJHDhw\/Lxo0b15T32fnNcuY+CEAAAs0Q0KRfppeaIdf8PT6PYd5EYM6cOSPDw8MyMTEReNr8f2dn5xrPa9kdO3bIAw88IN3d3aLCZ2ZmRg4cOCBtbW2ryvvs\/OZfCe6EAAQg0DwBjcpoREZ3\/DVnML3t5kn2lGkeaewf69PT08F459PljYAJi5Dx8XHp6OiIjMJo2bm5ORkZGYn1JQImFhEFIAABCDRFIGp6iTyZplAiYNLFlm9tKlj0MqIk\/G9jzfLysoyOjkpPT0\/dKSZT3giYwcFB6erqCn7c3t4efLggAAEIQKB1AkbI6DLsy7dsIE+mRaTz8\/OiH70WFhZkaGhIiMC0CDXL28MRl1pRFiNgbr\/9dnnooYfk9OnTDeXA2LarmNEPFwQgAAEIpEcAIZMOy6mpKdGPfSFg0mGbSS1JBczZs2dXEnf13nPnztXNgZmcnJTNmzcTgcnEe1QKAQhAYDUBTfg1ERn9jZlaWncF0e+4vmJHYGZnZwMxg4CJo1bg71uZQrITgMNJv+TAFOhUHg0BCFSewNLn98nrS1+SdW\/\/UZDwi5BJ1iV8HsO8SuK1E3PrJfGGf6cCZv\/+\/XLw4ME1S6l9dn6y14DSEIAABIojoEJm+e+Py1tu+DFCJoEbfB7DvBEwSZZRq0NVxJi9X2ol\/Gof8dn5Cd4BikIAAhAoBQEVMhe+9Tl5S+ePg9Ow11\/TLVe+ZzJI\/OVaS8DnMcwbAaNuq7WRnUnc7e\/vX1kHb29k19vbG5n\/goDh6wACEIBAOQmokPmn\/\/ZpuaLnWoRMHRchYMrZf3Oxymfn5wKQh0AAAhDIkIARMpf\/+obgqAIiMqth+zyGeRWByeId8dn5WfCiTghAAAJFEEDIRFP3eQxDwMS8aT47v4gvGZ4JAQhAICsCuo\/MK18\/umpqqeqrlnwewxAwCJisvkuoFwIQgEAhBGoJmbabPlGIPUU+FAFTJP2Cn+2z8wtGy+MhAAEIZEpAhcwP\/2KHvPGjp4JkX7OPTJWEjM9jGBEYIjCZfoFQOQQgAIGiCSx\/5xvy0l\/skHVX\/Ug02Xf9L3TIW995h1RByCBgiu59BT7fZ+cXiJVHQwACEMidgCb6Ln3hvkDEVOXQSJ\/HMCIwRGBy\/xLhgRCAAASKImDyY8JC5m03T8pl13YXZVZmz0XAZIa2\/BX77Pzy08dCCEAAAtkQMPkxP\/3+E\/KWf31lEJHxcQ8Zn8cwIjBEYLL5dqBWCEAAAg4QeOUbf3UpP+bK9Su7+urSa1\/yYxAwDnTCrEz02flZMaNeCEAAAi4RsKeVrvy3W4NzlvTyIdHX5zGMCAwRGJe+Z7AVAhCAQGYEVMic2\/tbQf1X\/vbNIuufCQ6JdDk\/BgGTWXcpf8U+O7\/89LEQAhCAQP4EzGqlt1x\/g1yz82Oy\/OyUs\/kxPo9hRGCIwOT\/7cATIQABCJScwMomeD+ck7d1\/bZcdsOP5Y2XnxLX8mMQMCXvaFma57Pzs+RG3RCAAAR8IGCSfC\/\/tdvkHR8fkwvPDAXNcmVayecxjAgMERgfvmNoAwQgAIHMCNjRmHd87EiQG2Omla56\/7HMnptGxQiYNCg6WofPznfUJZgNAQhAoBACJhrz9tt+X67Z8Wfy428PlX5ayecxjAgMEZhCvgh4KAQgAAEXCYSjMet\/8XJ5Zaa\/tKuVEDAu9rKUbPbZ+SkhohoIQAAClSNgVipd\/ZH75Orf3SvLz34mmFYqW5Kvz2MYERgiMJX74qHBEIAABNIgYO8bs2nf12Xd29YH00oXl+dLk+SLgEnD047W4bPzHXUJZkMAAhAoFQETjbnuvq9L26\/dthKN0bOVik7y9XkM8yoCc+LECdm9e3fQscfGxqSvry+2k585c0aGh4dlYmJCOjs715T32fmxcCgAAQhAAAINEVj+zjfkB\/f9lpgppYs\/mV+JxhR5JIHPY5g3AsYWItrb6okS0xuXl5dldHRUTp06JUeOHEHANPSaUggCEIAABKIImARf\/d2m+74eFDG5MUVFYxAwDvRVjb7MzMzIgQMHpK2tTcbHx6Wjo6NuFEYdq+X0IgLjgJMxEQIQgIADBPR0a43I6J4xOqVkR2Py3gAPAeNAhzFCZGRkJLA2\/O9wExYXF2Xfvn3S398flEXAOOBkTIQABCDgCAGzZ4yZUlKzX3vhi3Lh20Py1nd+WN72nslcWoKAyQVzaw8JR1w0IjM3NydG0IRr19\/rdcsttzSUAzM4OChdXV3BPe3t7cGHCwIQgAAEIFCLgE4pnf3o9aLHEJgppTyiMfPz86IfvRYWFmRoaEimp6elu7vbK2d5kwOTRMBovszRo0dlz549gZMbSeK1va5iRj9cEIAABCAAgXoEovJitHyW+8ZMTU2JfuwLAVPifppkCknL3nrrrYEabXQV0uTkpGzevJkITIn7AKZBAAIQKCuBcF6M2plVNMaOwMzOzgZiBgFT1p4hIuEpo1pJvJr7snPnTjl9+vSa1kQ52Of5wxK7E9MgAAEIeEdA82J0zxiT3GsamGU0xucxzJsppGaWUWvnaTQC46N69e7bgQZBAAIQKDkBs1+M2fTOmPv6j56SC88Mybq2drnyPZPB2UppXAiYNCjmUEetjezMfi+64iicxISAycExPAICEIAABFYIhDe9M7\/QKSVdqaRnKukqJV2t1OqFgGmVoMP3++x8h92C6RCAAAScJmCSey97R0cwpWRfGo3RE67TWG7t8xjmzRRSVj3ZZ+dnxYx6IQABCEAgnkA9EWMn+F7Vc6zpKSWfxzAETEwf89n58a8XJSAAAQhAIGsC5+77reARZq8Y+3kmwfftPcfksmuT7+Pi8xiGgEHAZP1uUj8EIAABCMQQqCdiWplSQsBUuOv57PwKu5WmQwACECgdgXoiptkpJZ\/HMCIwRGBK9xJjEAQgAIGqEqgnYpRJ0iklBExVe5KI+Oz8CruVpkMAAhAoLQEVMW3\/5ja5+nf3RtpoppTabhqUtps+UbcdPo9hRGCIwJT2JcYwCEAAAlUlEBeJMVNKyueq9x+riQkBU9UeRASmwp6n6RCAAASKJaAiJmqfGNuqC98ektdffkredvNk5ColBEyxPiz06T47v1CwPBwCEIAABOoSMPvE1JtO0gp0914VMlFLrX0ew5hCYgqJrxAIQAACECgpARUxZz96fbBb79tv+\/2aVtbKi0HAlNSxeZjls\/Pz4MczIAABCECgNQJGxIQPgAzXGpUX4\/MYRgSGCExrbxZ3QwACEIBA5gRqnWId9eD\/92S\/XFyeFz2C4G\/\/17wMDAzI9PT0msOMMzc64wcgYBAwGXcxqocABCAAgTQIqIh56S92yLv+8h9iq9P9Yl574WH5P699WH7vrikETCwxDwv4HH7z0F00CQIQgIDXBFTAvP7SXOS5SeGGm7yY\/\/zfRW77vWNEYLzuGRGNQ8BUzeO0FwIQgEC5CcTtEWNbP\/P1h+VXXrlHnvsFBEy5vZqBdQiYDKBSJQQgAAEINE2g0aRefYDPYxg5MDFdyGfnN\/32cCMEIAABCBRKoNGkXp\/HMAQMAqbQl5CHQwACEIBAcwSWPr9PXvnGX9VN6kXANMfWi7t8dr4XDqIREIAABCpMIO64AZ\/HMCIwRGAq\/OrTdAhAAAJuE4jLh0HAuO3flqz32fktgeFmCEAAAhAoBQGdRtLppKj9YXwew7yKwJw4cUJ2794ddKixsTHp6+uL7FyLi4uyc+dOOX36dPD73t5eOXDggLS1ta0p77PzS\/HmYQQEIAABCLRMoNZUks9jmDcC5syZMzI8PCwTExNBRzD\/39nZuapjLC8vy+joqPT09AQCx\/x706ZNMjIygoBp+TWiAghAAAIQyJtArakkBEzenmjieRp9mZmZWYmkjI+PS0dHR80ojP2I8L3273x2fhOYuQUCEIAABEpKIGpVks9jmDcRGBUsepkoSvjf9fpbIwJmcHBQurq6gmra29uDDxcEIAABCECgTATOfvR6+afO35YLv7kzMGthYUGGhoY4C6lMTgrbEo64qCiZm5uLnBay7zX5MNu3b4+M1hj1at+jYkY\/XBCAAAQgAIEyEdCEXj0v6T\/873fKP\/50\/YppnEZdJi+FbGlGwJj8F60qLol3cnJSNm\/eTASmxH0A0yAAAQhAQEQTev\/xtfXy2r\/\/M5mdnZWpKU6jLnW\/SDqF1Ih40Qb7PH9YaodiHAQgAAEINEXAPmbg9CuXy8DAAFNITZHM6abwlFG9JN64lUe2yQiYnBzIYyAAAQhAIDUCZln18+\/dhYBJjWpGFTW6jFofr+Lm3LlzNaeNEDAZOYlqIQABCEAgFwImCnP23x2QP\/j0ISIwuVBv4SG1NrIzEZf+\/n658cYbV21iZx63ZcsWOXz4sGzcuHGVBS5FYObn5+Xhhx+WO+64w4lVUtjbQmdv4Fb4NgCphSIu8XXJVnUJ9rbQMa1bNQpz6v++KINPvIaASQepW7W4JGBcslV7AfZm+y7AF76GAH2hmn3BXpH0n\/7q89Ld3Z0tiJxr92YfmKy4ufTiu2QrAiarHvtmvfSHbBm7xNclW\/luSLffahTmy3\/zLXnP\/scQMOmiLX9t5sW3N7Irq9VmwyIXbFWG2JttT4IvfA0B+kJ1+8K\/+t7\/kLd++U\/kh3f\/TwRMtt2gfLXrXKzuYqhr6bkgAAEIQAACrhH43LtfkHf\/x4fk7bf9vmum17WXKaQG3KkiRj9cEIAABCAAAdcI+Hr8DQLGtZ6IvRCAAAQgAAEICAKGTgABCEAAAhCAgHMEEDDOuQyDIQABCEAAAhBAwNAHIAABCEAAAhBwjgACxjmXYTAEIAABCEAAAggY+gAEIAABCEAAAs4RQMDUcZke+njo0KGgxPT0dGk2AdKDK3fs2BEcSNnb29vwoZTHjh1rqGzavTiJvTbzTZs2yZEjR6SzszNtk+rWl8Re+\/ytWudpZWl8EluNHeZssJ6eHunr68vSvDV1J7HXZqsVlZ3v4uLiqnPWivjOaJRvmK1x1NjYWK59olF71T67rAvfDWYT1KL6bq4vdkEPQ8DUAK+dTwdTPeDxueeeW\/n\/8GGPefvNHny2bdsmo6OjEjcQmRepUbGTZpuS2KtfqjMzMysiS\/99\/PjxyEM207TRriuJvXYf0X6R5JTzNOxPYqv9PDN45T1YJbVXeXZ0dOQ6oDbbF0zbdGAdGRkJBtvh4WGZmJjITYAn5Wu3NdyX0+ifcXUksdeIQ2Wr5\/mU\/bvBiK0HHnggsFf5FvUHZJwfXP49AqaG9\/TLUy99YezTrIs+DCv8xRj3Ymg7Tp48Kbfddpu88soruUdgktpru6OIQaAVe\/MeBJqxVQeCu+++W86fPy\/bt2\/PVRwksbcM71wSe7Xs\/v375eDBg2tOtM9rgEhir21TWByU0d5w28r+3RD+Y0z78\/333y933nlnboI2Lz8W+RwETAT9cIi9yJB72LzwIBk3aOrvzV8sdnQjr06X1N6iBUwr9tqiNw++zdiqNr7vfe+TL3\/5y7GRu7TbkMTeogbVelGJeu9aeMBKm10j9SXhG47Ilf27ISoCk7fNSfhGCRiNlvf395cmFaGRPlX2MgiYOgLG7mxFh7ONmeGIS6N\/+RX1BdusvdrevKdk9JnN2Gum6PKel09qq\/aVo0ePyh\/\/8R\/Lvn37ChEwdhi9Xt+18x1M3887pyQJX32\/5ubmAlOLyptLYq9hWqRQTGqv+UNSI8q7du0KouN5XknsDZ\/+bf6d97RtnnyKeBYCBgGTab9L8tKH\/yr87Gc\/m3sSb7P2GvFzzz335GZzElvtELaei9JI7lTaHSOJvVrWZhn+d9q2RdWXxF6TV2REVtnttf8gMrl+eef3JeEblVOSt91J7FW+dqL0wMBAMIUfl6+YR7\/26RkImDoCxnQ2l6eQTPOKjMDYXzRxU17mxS9CvBgRktRewzjvfpIkpK1lv\/nNb67K6cr7yzSJveHXMm+2SftCrSmDPBk3w9dEjvKOZlSFbxkiXT4JlnBbEDA1vGtPGZUhodCYGQ67h\/8qqNVZixIwSe0tYnWBzSypvfa9eYfjk9hqL0+3bc4zFJ\/E3loCJs8cgiT2ht\/DIr4zktirfIsQhc2+a2UQiEn5htua96o0n4WLaRsCpoaX7b9mXF9GbaIaeSe9hb8k45Z9FxF2r\/eXfiP22tGavMVXkmWodjuLGriS2OvakvqweG0k0pj2AJOErz7brEi79957C1kZk8TeqCmkPKdrk36X2WKnra0tmLI1S+zT9nuV60PA1PG+ixvZ1QoJFxWBUbz1Nquy7a0VJcg7ebNRe40w3L17d9CL8k7iTcK2DAImqb12DkERbJPaa29k54K9RSxFDn\/dJnnXTCJsUe9a0v5g998i9uCqgrBBwFTBy7QRAhCAAAQg4BkBBIxnDqU5EIAABCAAgSoQQMBUwcu0EQIQgAAEIOAZAQSMZw6lORCAAAQgAIEqEEDAVMHLtBECEIAABCDgGQEEjGcOpTkQgAAEIACBKhBAwFTBy7QRAhCAAAQg4BkBBIxnDqU5EHCRgO6hogdM7t27V\/I+k8dFXtgMAQiIIGDoBRCAQOEE7LOaCjcGAyAAAScIIGCccBNGQsBvAroL86233ird3d1+N5TWQQACqRFAwKSGkoogAIEwgVrnLtlnBelZMffff7\/ceeedwZk85p6TJ08G1RW1LT\/ehAAEyk0AAVNu\/2AdBJwnEHXIpUZc9BoZGQnOyjp69Kjs2bMn+JkefKfXgQMHRMVN3odkOg+cBkCgIgQQMBVxNM2EQFEEwic1h\/+tAkWvvr6+QMwMDw\/LxMREISckF8WI50IAAskJIGCSM+MOCEAgIQE74mJPH+mKIzv\/Jfy7hI+hOAQgUCECCJgKOZumQqAoArYweeihh6SjoyOIuISXTyNgivIQz4WAewQQMO75DIsh4BwBM220bds2+cpXvrIyRRRePs0UknOuxWAIFEYAAVMYeh4MgWoR0FyX3bt3S29v76oEXaWg0Ri9zAokXXmkCb56IWqq1U9oLQQaJYCAaZQU5SAAgZYIqBDZsWOH3HXXXYFgUbFiL582lbOMuiXM3AyByhBAwFTG1TQUAhCAAAQg4A8BBIw\/vqQlEIAABCAAgcoQQMBUxtU0FAIQgAAEIOAPgf8Pm\/D9KUv8MQ8AAAAASUVORK5CYII=","height":337,"width":560}}
%---
