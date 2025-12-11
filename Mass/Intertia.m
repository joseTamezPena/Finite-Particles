%[text] # The mass of a finite particle
%[text] 
%[text] This script will show that a finite particle oposes to be acelerated i.e., it has a mass.
%[text] The first finite size particle model is very simple:
%[text] - Two point elements with charge Q/2 on the horizontal axis separated by a distance of 2\*r\_o 
%[text] - The two elements are bounded by a nuclear force
%[text] - The particle is acelerated with constant acceleration \

% Solving for the arrival time of a particle at constant aceleration
% a   : Aceleration
% r_o : Particle radious
% c   : Speed of light
% K   : eletrostatic constant
% t   : time to reach a target
% Q   : Particle charge
syms a Q r_o c K t real;
assumeAlso(a > 0)
assumeAlso(r_o > 0)
assumeAlso(c > 0)
assumeAlso(K > 0)
assumeAlso(t > 0)
assumeAlso(Q > 0)

%%
%[text] ## Solving for the arrival time
%[text] I'll solve for the forces that the particle experiment at time equal zero. 
eq1 = -a/2*t^2-c*t+2*r_o==0 % right corpuscle moving towards the left element %[output:00446ff5]
t_1=solve(eq1,t) % The time it takes to reach the left element  %[output:58cbf678] %[output:7fe5b240]
t_1 = t_1(1) %[output:5fd8352d]
eq2 = -a/2*t^2+c*t-2*r_o==0 % The equation of the corpuscle traveling from left to right element %[output:3198a421]
t_2= solve(eq2,t) % The time it takes to reach the right element %[output:602e0f15] %[output:1088dd07]
t_2=t_2(2) % I only take the shortest time %[output:3d8c3eda]
vpa(subs(t_1,[a,r_o,c],[1,1.0e-15,1]))-vpa(subs(t_2,[a,r_o,c],[1,1.0e-15,1])) %[output:936c4610]
D_a=simplify(1/(t_2^2)-1/(t_1^2)) % The force is due to differences in arrival time %[output:208c73d9]
D_a= taylor(D_a,a,Order=4) % Lets take an aproximation %[output:437b6ce5]
F=K*Q^2/(4*c^2)*D_a % The reaction force experimented by the acelerated particle  %[output:636e6739]
%[text] Hence if you acelerate a electically charged finite-sized particle, it will opose the aceleration. In other words particles have intertia and is origen is the fact that electric charges have dimensions. 
%%
%[text] ## Energy and Mass
U= -K*Q^2/4*int(1/t^2,t,Inf,2*r_o) % The energy required to bound the two elements of the particle %[output:6eb36bdc]
syms U
F_net=expand(simplify(subs(F,Q^2,8*U*r_o/K))) % The reaction force in units of the energy %[output:89bc52ad]
%%
%[text] where $r\_o$ is very small and c is very large. 
%[text] The first term of the total reactive force is:
F_net = taylor(F_net,Order=2) % The first term %[output:5d82f2e5]
%[text] 
%[text] The sum of all the forces is equal zero, and $F\_{net}+ma=0$, then $m=\\frac{2U}{c^2}$ or the $U \\approx mc^2$  the famus Einstein equation. 
%%
%[text] ## The expanded model of the particle
%[text] 
%[text] Let use a more complex model of the particle this time lets use 6 elements to model the finite-sized paticle. Each element is in the corner of a regular octahedron:
%[text]{"align":"center"} ![](text:image:5e5f).
%[text] The particle resides in the origen, and each one of the six elements is at a distance $r\_o$ form the origen. The constant aceleration is in the positve $x$ direction. The analysis involve estimating the time of arrival of each corpuscle to the other elements.
eq3 = (c*t)^2 - r_o^2 - (r_o-1/2*a*t^2)^2 == 0 % corpuscle moving from off x-axis left %[output:5f1045af]
t_dL= solve(eq3,t) % The time it takes for corpuscles to reach diagonal left element  %[output:6efb3ffb] %[output:4e10cb92]
t_dL = t_dL(2) %[output:2621b0ae]
eq4 = (c*t)^2 - r_o^2 - (r_o + 1/2*a*t^2)^2 == 0 % diagnoal corpuscle to right element %[output:6832ac3d]
t_dR= solve(eq4,t) % Diagonal Time to reach the right element %[output:063e43c3] %[output:2bdd1101]
t_dR = t_dR(2) % The shortest time %[output:06871cbc]

eq5= (1/2*a*t^2)^2-(c*t)^2 + (2*r_o)^2==0 % The equation from off x-axis to the opposite element %[output:06c8b3db]
t_5= solve(eq5,t) %[output:9fd4cdc5] %[output:0ca207f2]
t_5=t_5(2) % The shortest time %[output:84461151]
eq6= (1/2*a*t^2)^2-(c*t)^2 + 2*r_o^2 == 0 % Off axis to the next off axis %[output:89560e0f]
t_6= solve(eq6,t) %[output:6333bfee] %[output:6b6c9327]
t_6=t_6(2) % The shortest time %[output:68b0dc63]


%%
%[text] Now, I'll compute the reactive forces for each one of the element to element interactions
syms d v real positive
cos_L = cos(pi/4);
cos_R = cos(pi/4);
cos_5 = (1/2*a*t_5^2)/sqrt((2*r_o)^2+(1/2*a*t_5^2)^2);
cos_6 = (1/2*a*t_6^2)/sqrt(2*r_o^2+(1/2*a*t_6^2)^2);

F_1 = -1/(c*t_1)^2 - 4/(c*t_dL)^2*cos_L; % Left element
F_2 = +1/(c*t_2)^2 + 4/(c*t_dR)^2*cos_R; % Right element
F_3 = 1/(c*t_5)^2*cos_5; % The two off axis elements

F_3 = F_3 + 2/(c*t_6)^2*cos_6; % The  two sides of the off-axis elements
F_3 = F_3 - 1/(c*t_dL)^2*cos_L; % add the right element to the off-axis
F_3 = F_3 + 1/(c*t_dR)^2*cos_R; % add the left element to off-axis

F_netp = simplify(F_1 + F_2 + 4*F_3,steps=100); % One left, one right plus the four off-axis
F_net = K*(Q/6)^2*simplify(taylor(F_netp,a,Order=2)) %[output:9b789c89]
%%
% The Binding energy of the regular octahedron
U_1 = K*(Q/6)^2/(2*r_o) % Right to left %[output:4ec02bd5]
U_2 = K*(Q/6)^2/(sqrt(2)*r_o) % off diagonal  %[output:3143e79d]
U_net = simplify(3*U_1 + 12*U_2) %[output:61b52991]
syms U_net Q2
qeU = U_net - (K*Q2*(4*sqrt(sym(2)) + 1))/(24*r_o) == 0 %[output:9f04f30d]
simplify(solve (qeU,Q2)) %[output:2736dd13]
F_net %[output:8adbf34a]
F_u=simplify(subs(F_net,Q^2,(24*U_net*r_o)/(K*(4*sqrt(sym(2)) + 1)))) %[output:6a5836a1]
vpa(subs(F_u,[a,c,r_o],[a,1,1])) %[output:08177f46]
%[text] The inertial mass of the regular octahedron is lower that that of the two element particle.
%%
%[text] ## 
%[text] ## Add the central element to the regular octahedron
eq1p = -a/2*t^2-c*t+r_o==0 % right corpuscle moving towards the left element %[output:9e63200e]
t_1p=solve(eq1p,t) % The time it takes to reach the left element  %[output:406e4844] %[output:28fade7c]
t_1p = t_1p(1) %[output:6aba67dc]
eq2p = -a/2*t^2+c*t-r_o==0 % The equation of the corpuscle traveling from left to right element %[output:9439c2ae]
t_2p= solve(eq2p,t) % The time it takes to reach the right element %[output:401ad15f] %[output:5019f9e6]
t_2p=t_2p(2) % I only take the shortest time %[output:0bf4bcc9]
eq_c = (1/2*a*t^2)^2-(c*t)^2 + r_o^2==0 % Points off axis %[output:7b2c02ec]
t_c= solve(eq_c,t) %[output:073e3cd3] %[output:47a4f2ae]
t_c=t_c(2) % The shortest time %[output:8df05cfd]
cosp = (1/2*a*t_c^2)/sqrt(r_o^2+(1/2*a*t_c^2)^2);

F_cp = -1/(c*t_1p)^2 + 1/(c*t_2p)^2 + 4/(c*t_c)^2*cosp;
F_1p = F_1 - 1/(c*t_1p)^2; % Left element
F_2p = F_2 + 1/(c*t_2p)^2; % Righ element
F_3p = F_3 + 1/(c*t_c)^2*cosp; %Off axis elements
F_netcp = simplify(F_cp + F_1p + F_2p + 4*F_3p,steps=100); % One left, one right plus the four off-axis
F_netc = K*(Q/7)^2*simplify(taylor(F_netcp,a,Order=2)) %[output:464dd6d8]
%%
% The Binding energy of the regular octahedron with a central element
U_p = K*(Q/7)^2/(r_o) % Elements pairs separed by r_o %[output:8e44930c]
U_1 = K*(Q/7)^2/(2*r_o) % Elements  pairs separed by the sqrt(2)r_o %[output:8df842b5]
U_2 = K*(Q/7)^2/(sqrt(2)*r_o) % Elements pairs separated by 2r_o %[output:87d4e21a]
U_net = simplify(3*U_1 + 12*U_2 + 6*U_p) % The 21 pairs of element to element interactions. %[output:83cdff46]
syms U_net Q2
qeU = U_net - (3*K*Q2*(4*sqrt(sym(2)) + 5))/(98*r_o) == 0 %[output:0e3ea66d]
simplify(solve (qeU,Q2)) %[output:73c1cb3e]
F_netc %[output:9591f1ac]
F_uc=simplify(subs(F_netc,Q^2,(98*U_net*r_o)/(3*K*(4*sqrt(sym(2)) + 5)))) %[output:8b762e00]
vpa(subs(F_uc,[a,c,r_o],[a,1,1])) %[output:96f77fb4]
%[text] The inertial mass of the regular octahedron with a central element is even lower. 
%[text] The implications is that it may be possible to find a charge distribution that will explain the oserved mass of charged particles.
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[text:image:5e5f]
%   data: {"align":"baseline","height":141,"src":"data:image\/jpeg;base64,\/9j\/4AAQSkZJRgABAQAAAQABAAD\/2wCEAAkGBwgHBgkIBwgKCgkLDRYPDQwMDRsUFRAWIB0iIiAdHx8kKDQsJCYxJx8fLT0tMTU3Ojo6Iys\/RD84QzQ5OjcBCgoKDQwNGg8PGjclHyU3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3N\/\/AABEIAI0AjQMBIgACEQEDEQH\/xAAbAAACAwEBAQAAAAAAAAAAAAACAwAFBgQBB\/\/EAD4QAAEDAgQCBwYEBAUFAAAAAAEAAgMEEQUhMUESUQYTImFxgbEUIzJikaEzQ1JyJEJTspLBwtLiFTRzgtH\/xAAVAQEBAAAAAAAAAAAAAAAAAAAAAf\/EABQRAQAAAAAAAAAAAAAAAAAAAAD\/2gAMAwEAAhEDEQA\/APrQRApYKIKg0YKWCiBQMBRXSwUQKgYCiBSwVx1VTJK91LRuAePxZrXEXcObu7a9zsCDJMXoocSjw+SW08jSR2Twg5WaXaBxvkN1YA3CwuNQxDE2xNj93FEGWPaN3cTnXJO\/YN87q2wPGTG5lFXyEg2bFUO35NcefI7756hpgUQKWD3ogUBqLwL1BFFFEEUUUQUgKIFLRgqoYCvQgCIFFGCiBQArmfI+pkdT0ry1rTaWcfy\/K35u\/bxyUBTzyTyOpqR3CW5TTf0+4c3emp5FsMMcETY4m8LB535knc8zujhhjgibFEwNY3QBFkMzoEGPxBwkrauXUCp1NsuGzDr4FLexr2ljwHNIsQRcEIWky0XWZh0jC\/W2bru9SjuHAOabgi4Kot8Dxp0DmUeISExkhsU7zmOTXH0Pkc8zqQV88e0OBDgCCLEHdW+B42aQtpK55NPkI5nG\/V\/K48uR+qg1wRBLBRAoCUUUQRRRRBQgoggCIFVBgogUu6U3jrnuigc5lO08MswNi4\/pYfV22gzvwlES+skfBTuLImHhlnbrfdrO\/mdvHTuiijgibFCwMjaLNaNkUUbIY2xxMaxjRZrWiwARFQAuPGHmPCqtzTY9U5rT3kWH3IXaqvpC7+Cjh\/rTsH+E8f8AoQUtg21hYN0C5479U25uQLak6eK6XJDBbjGeTjrffPfxQCQluG1rhOKW5UWeB43\/ANP4aWtdek0ZKT+D3H5O\/bw02A0uvnDgrPAsbOGltNVuJob2a4\/kf8P7fDSDbAoktrg4BzSCCLgjQoggJRRRBnwV6EAK8poXYk48JLaMGzpAbGU7tb3cz5DmKPYmPxB5ZG4spWm0krci8jVrT6u8hnmLVkbIo2xxsaxjRZrWiwA5L33VPEG9iKNgsBk0ABD1zTk0OceTWnv302O\/qoCXiHikcRaKw+ZwHpdDwSnJ0gb+xv8A9vvdUGe5UnSCRvtNJESBwh8mZ37LR\/eVcGBjr8Zc4HUFxtvtpufsqWua04lNwNAbHExuQt2jxE\/YtUFYZAfgY53la\/1sk3eJZCYsjY3BF9\/roPqux4tdc4HvZDbly5fX6opfGwm17HkcivHAWumuF0gxNbky7ANA3TbbTZVAkIHBeu6wfpd9kDntHxAt8Rlvvpt6ILPAsadhZEFQS6hOm5h8Pl7ttsltmPbIxr2Oa5jgC1zTcEcwvmhzFx6qywPGX4U8RS8T6Fxza0XMJ5tHLmPMclBvQV6lQysmjbJE9r43gOa5puHA6EJl0GOxhtYaQmhY2VzTd8Ljw9c3dodsfKx03QYT01waqDY6s1FA9h4eCbsxgg6XbkALW7VuSsAVm+k\/R4VvFXULB7WB7xgy64f7vXQ7KjeU4ppImzUpifG4dmSMggi1tR3AfRMIXxKhqamhlMtFPLTSXzMbi03Gzhv4ELVYb05xCENbiEEVUwfzt92\/6fCfsoPoS8VNh3SnCK\/haKn2eV2QjqfdknkDoT4Eq5tkqPCqS3WSVUv6535\/tsz\/AEq7Lg0FztALlVNNEW4fAXfG6MOf+45n7kqCtmbYrkbYyS2IyIuLjLIKyqGrgFzJLnezhv3D6IoHBLIRzyRwi80jIxze4D1VfLi+Ht+GpbJ\/4gXj6jJEdTglnkq2XHos+qp5XHm8ho9SfsuOXG6p\/wCHFDH3m7z\/AJKi4dG29w2x5tyS3MeB2X3AH84\/zWflrq2T4qp4B2YA37gX+645W9YQZi6Q3y4yXeqDaYN0qiwScQzzNkpHu7UUbg9zHE6sGpz1bbO+Wev0gHJYjoP0U9h4MSxOMe1awQn8kcz83p9VtlBnQUV0C9BVRQdI8A9qLq2hYPadZIx+b3j5vVZNgvmNO9fTQqDH8E64urKFnvtZYm\/md4+b18UVlgy678Pr67D\/APsaqWFv9O\/Ez\/Cch5LmjAcARoU5jVBfydL5zh9THW0gcTE4CWndYi4tfhdy118kyr6Ygjho8OJaMgZpg23k0H1VC1gNwQCDzXJCwsYYjrESw+Wn2IQWNTj+JTX4Xww\/siufq4n0VW+pqpi\/rqqof2s\/eFo0GwsEwhKDfivu7vQJ6pgNwxtzqbZnzUITXBCQqEkICE4hA6wFyQBzKBLgACTkAt30M6LdQY8TxOP33xQQOH4fzOH6uQ28dF9D+jOceJ4nFmDxU8Dxpye4c+Q211tbcXUBhEClgokGfUUUVR6CiQIgUFHjmDlznVlEy8hzliH8\/wAw+b18VSx2c0Oabg6ELcAqkxfCzxOrKRlyTeaJurvmb38xv46xVQ1t1zVcfV1LH7St4f8A2GY+xP0VjA0PAcw8TSLgjQqV9MZKN5YCXx9tttbjbzFx5oKshJY3J2VruOy6MnAFpuCLhLYOyf3O5czyQLLUBCeQgIQJItcrT9FOjomdHiOIxnqx2oIXD4uTnDlyHnyQdGcB9rLK+uZ\/Dg3hid+af1H5eQ38LX2181QwFEClgolAwFElgoroKJRRRVEUCiiAgUSWEQQVmIUJhe6qpmFzXdqaJoue9zRz5jfXXXyBrXsa9hDmuza4aEK2GRXDNAaWR00DS6BxvLG0XLTu9o9R5jPIxWYlp\/Zp5ae1hG7s\/tOY+xt5JEY7HmeXM8lfdIaZvFTVsXC5kg6oubnfVzT4fF9QqVguwHM3vvdFLIVtgGCCvIqqxv8ACA3Ywj8bv\/b\/AHeGswXCDiDhNUN\/g2nIf1jy\/b6+GuvbYAACwGgGyIYEQKAFENVQaIFAvQUDAUQKXdEFBTKKKKoiiiiCKBRRAYK98EARoK\/E6YtoKmOJnFC8cfABnG8HiDm91xm3vuN70mA0TcXijqSSaEgHiz993C+duZ8hvbWAlE0BoAaAANgiiYA1oa0ANAsABoEYSwjCIMFECgXoQMBRAoAvWopgRXSxqiQf\/9k=","width":141}
%---
%[output:00446ff5]
%   data: {"dataType":"symbolic","outputData":{"name":"eq1","value":"2\\,r_o -c\\,t-\\frac{a\\,t^2 }{2}=0"}}
%---
%[output:58cbf678]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:7fe5b240]
%   data: {"dataType":"symbolic","outputData":{"name":"t_1","value":"-\\frac{c-\\sqrt{c^2 +4\\,a\\,r_o }}{a}"}}
%---
%[output:5fd8352d]
%   data: {"dataType":"symbolic","outputData":{"name":"t_1","value":"-\\frac{c-\\sqrt{c^2 +4\\,a\\,r_o }}{a}"}}
%---
%[output:3198a421]
%   data: {"dataType":"symbolic","outputData":{"name":"eq2","value":"-2\\,r_o +c\\,t-\\frac{a\\,t^2 }{2}=0"}}
%---
%[output:602e0f15]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:1088dd07]
%   data: {"dataType":"symbolic","outputData":{"name":"t_2","value":"\\left(\\begin{array}{c}\n\\frac{c+\\sqrt{c^2 -4\\,a\\,r_o }}{a}\\\\\n\\frac{c-\\sqrt{c^2 -4\\,a\\,r_o }}{a}\n\\end{array}\\right)"}}
%---
%[output:3d8c3eda]
%   data: {"dataType":"symbolic","outputData":{"name":"t_2","value":"\\frac{c-\\sqrt{c^2 -4\\,a\\,r_o }}{a}"}}
%---
%[output:936c4610]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-0.0000000000000000000000000000040000000001946799263974006163732"}}
%---
%[output:208c73d9]
%   data: {"dataType":"symbolic","outputData":{"name":"D_a","value":"\\frac{a^2 }{{{\\left(c-\\sqrt{c^2 -4\\,a\\,r_o }\\right)}}^2 }-\\frac{a^2 }{{{\\left(c-\\sqrt{c^2 +4\\,a\\,r_o }\\right)}}^2 }"}}
%---
%[output:437b6ce5]
%   data: {"dataType":"symbolic","outputData":{"name":"D_a","value":"-\\frac{a}{r_o }-\\frac{a^3 \\,r_o }{c^4 }"}}
%---
%[output:636e6739]
%   data: {"dataType":"symbolic","outputData":{"name":"F","value":"-\\frac{K\\,Q^2 \\,{\\left(\\frac{a}{r_o }+\\frac{a^3 \\,r_o }{c^4 }\\right)}}{4\\,c^2 }"}}
%---
%[output:6eb36bdc]
%   data: {"dataType":"symbolic","outputData":{"name":"U","value":"\\frac{K\\,Q^2 }{8\\,r_o }"}}
%---
%[output:89bc52ad]
%   data: {"dataType":"symbolic","outputData":{"name":"F_net","value":"-\\frac{2\\,U\\,a}{c^2 }-\\frac{2\\,U\\,a^3 \\,{r_o }^2 }{c^6 }"}}
%---
%[output:5d82f2e5]
%   data: {"dataType":"symbolic","outputData":{"name":"F_net","value":"-\\frac{2\\,U\\,a}{c^2 }"}}
%---
%[output:5f1045af]
%   data: {"dataType":"symbolic","outputData":{"name":"eq3","value":"c^2 \\,t^2 -{r_o }^2 -{{\\left(r_o -\\frac{a\\,t^2 }{2}\\right)}}^2 =0"}}
%---
%[output:6efb3ffb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:4e10cb92]
%   data: {"dataType":"symbolic","outputData":{"name":"t_dL","value":"\\left(\\begin{array}{c}\n\\sqrt{\\frac{2\\,{\\left(a\\,r_o +\\sqrt{c^4 +2\\,a\\,c^2 \\,r_o -a^2 \\,{r_o }^2 }+c^2 \\right)}}{a^2 }}\\\\\n\\sqrt{\\frac{2\\,{\\left(a\\,r_o -\\sqrt{c^4 +2\\,a\\,c^2 \\,r_o -a^2 \\,{r_o }^2 }+c^2 \\right)}}{a^2 }}\n\\end{array}\\right)"}}
%---
%[output:2621b0ae]
%   data: {"dataType":"symbolic","outputData":{"name":"t_dL","value":"\\sqrt{\\frac{2\\,{\\left(a\\,r_o -\\sqrt{c^4 +2\\,a\\,c^2 \\,r_o -a^2 \\,{r_o }^2 }+c^2 \\right)}}{a^2 }}"}}
%---
%[output:6832ac3d]
%   data: {"dataType":"symbolic","outputData":{"name":"eq4","value":"c^2 \\,t^2 -{r_o }^2 -{{\\left(r_o +\\frac{a\\,t^2 }{2}\\right)}}^2 =0"}}
%---
%[output:063e43c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:2bdd1101]
%   data: {"dataType":"symbolic","outputData":{"name":"t_dR","value":"\\left(\\begin{array}{c}\n\\sqrt{\\frac{2\\,{\\left(\\sqrt{c^4 -2\\,a\\,c^2 \\,r_o -a^2 \\,{r_o }^2 }-a\\,r_o +c^2 \\right)}}{a^2 }}\\\\\n\\sqrt{-\\frac{2\\,{\\left(a\\,r_o +\\sqrt{c^4 -2\\,a\\,c^2 \\,r_o -a^2 \\,{r_o }^2 }-c^2 \\right)}}{a^2 }}\n\\end{array}\\right)"}}
%---
%[output:06871cbc]
%   data: {"dataType":"symbolic","outputData":{"name":"t_dR","value":"\\sqrt{-\\frac{2\\,{\\left(a\\,r_o +\\sqrt{c^4 -2\\,a\\,c^2 \\,r_o -a^2 \\,{r_o }^2 }-c^2 \\right)}}{a^2 }}"}}
%---
%[output:06c8b3db]
%   data: {"dataType":"symbolic","outputData":{"name":"eq5","value":"4\\,{r_o }^2 -c^2 \\,t^2 +\\frac{a^2 \\,t^4 }{4}=0"}}
%---
%[output:9fd4cdc5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:0ca207f2]
%   data: {"dataType":"symbolic","outputData":{"name":"t_5","value":"\\left(\\begin{array}{c}\n\\sqrt{\\frac{2\\,{\\left(c^2 +\\sqrt{-{\\left(c^2 +2\\,a\\,r_o \\right)}\\,{\\left(-c^2 +2\\,a\\,r_o \\right)}}\\right)}}{a^2 }}\\\\\n\\sqrt{\\frac{2\\,{\\left(c^2 -\\sqrt{-{\\left(c^2 +2\\,a\\,r_o \\right)}\\,{\\left(-c^2 +2\\,a\\,r_o \\right)}}\\right)}}{a^2 }}\n\\end{array}\\right)"}}
%---
%[output:84461151]
%   data: {"dataType":"symbolic","outputData":{"name":"t_5","value":"\\sqrt{\\frac{2\\,{\\left(c^2 -\\sqrt{-{\\left(c^2 +2\\,a\\,r_o \\right)}\\,{\\left(-c^2 +2\\,a\\,r_o \\right)}}\\right)}}{a^2 }}"}}
%---
%[output:89560e0f]
%   data: {"dataType":"symbolic","outputData":{"name":"eq6","value":"2\\,{r_o }^2 -c^2 \\,t^2 +\\frac{a^2 \\,t^4 }{4}=0"}}
%---
%[output:6333bfee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:6b6c9327]
%   data: {"dataType":"symbolic","outputData":{"name":"t_6","value":"\\left(\\begin{array}{c}\n\\sqrt{\\frac{2\\,{\\left(\\sqrt{c^4 -2\\,a^2 \\,{r_o }^2 }+c^2 \\right)}}{a^2 }}\\\\\n\\sqrt{-\\frac{2\\,{\\left(\\sqrt{c^4 -2\\,a^2 \\,{r_o }^2 }-c^2 \\right)}}{a^2 }}\n\\end{array}\\right)"}}
%---
%[output:68b0dc63]
%   data: {"dataType":"symbolic","outputData":{"name":"t_6","value":"\\sqrt{-\\frac{2\\,{\\left(\\sqrt{c^4 -2\\,a^2 \\,{r_o }^2 }-c^2 \\right)}}{a^2 }}"}}
%---
%[output:9b789c89]
%   data: {"dataType":"symbolic","outputData":{"name":"F_net","value":"-\\frac{\\sqrt{2}\\,K\\,Q^2 \\,a}{18\\,c^2 \\,r_o }"}}
%---
%[output:4ec02bd5]
%   data: {"dataType":"symbolic","outputData":{"name":"U_1","value":"\\frac{K\\,Q^2 }{72\\,r_o }"}}
%---
%[output:3143e79d]
%   data: {"dataType":"symbolic","outputData":{"name":"U_2","value":"\\frac{\\sqrt{2}\\,K\\,Q^2 }{72\\,r_o }"}}
%---
%[output:61b52991]
%   data: {"dataType":"symbolic","outputData":{"name":"U_net","value":"\\frac{K\\,Q^2 \\,{\\left(4\\,\\sqrt{2}+1\\right)}}{24\\,r_o }"}}
%---
%[output:9f04f30d]
%   data: {"dataType":"symbolic","outputData":{"name":"qeU","value":"U_{\\textrm{net}} -\\frac{K\\,Q_2 \\,{\\left(4\\,\\sqrt{2}+1\\right)}}{24\\,r_o }=0"}}
%---
%[output:2736dd13]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\frac{24\\,U_{\\textrm{net}} \\,r_o }{K\\,{\\left(4\\,\\sqrt{2}+1\\right)}}"}}
%---
%[output:8adbf34a]
%   data: {"dataType":"symbolic","outputData":{"name":"F_net","value":"-\\frac{\\sqrt{2}\\,K\\,Q^2 \\,a}{18\\,c^2 \\,r_o }"}}
%---
%[output:6a5836a1]
%   data: {"dataType":"symbolic","outputData":{"name":"F_u","value":"\\frac{U_{\\textrm{net}} \\,a\\,{\\left(\\frac{4\\,\\sqrt{2}}{93}-\\frac{32}{93}\\right)}}{c^2 }"}}
%---
%[output:08177f46]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-0.28325963172588838499777682906625\\,U_{\\textrm{net}} \\,a"}}
%---
%[output:9e63200e]
%   data: {"dataType":"symbolic","outputData":{"name":"eq1p","value":"r_o -c\\,t-\\frac{a\\,t^2 }{2}=0"}}
%---
%[output:406e4844]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:28fade7c]
%   data: {"dataType":"symbolic","outputData":{"name":"t_1p","value":"-\\frac{c-\\sqrt{c^2 +2\\,a\\,r_o }}{a}"}}
%---
%[output:6aba67dc]
%   data: {"dataType":"symbolic","outputData":{"name":"t_1p","value":"-\\frac{c-\\sqrt{c^2 +2\\,a\\,r_o }}{a}"}}
%---
%[output:9439c2ae]
%   data: {"dataType":"symbolic","outputData":{"name":"eq2p","value":"-r_o +c\\,t-\\frac{a\\,t^2 }{2}=0"}}
%---
%[output:401ad15f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:5019f9e6]
%   data: {"dataType":"symbolic","outputData":{"name":"t_2p","value":"\\left(\\begin{array}{c}\n\\frac{c+\\sqrt{c^2 -2\\,a\\,r_o }}{a}\\\\\n\\frac{c-\\sqrt{c^2 -2\\,a\\,r_o }}{a}\n\\end{array}\\right)"}}
%---
%[output:0bf4bcc9]
%   data: {"dataType":"symbolic","outputData":{"name":"t_2p","value":"\\frac{c-\\sqrt{c^2 -2\\,a\\,r_o }}{a}"}}
%---
%[output:7b2c02ec]
%   data: {"dataType":"symbolic","outputData":{"name":"eq_c","value":"{r_o }^2 -c^2 \\,t^2 +\\frac{a^2 \\,t^4 }{4}=0"}}
%---
%[output:073e3cd3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Solutions are only valid under certain conditions. To include parameters and conditions in the solution, specify the 'ReturnConditions' value as 'true'."}}
%---
%[output:47a4f2ae]
%   data: {"dataType":"symbolic","outputData":{"name":"t_c","value":"\\left(\\begin{array}{c}\n\\sqrt{\\frac{2\\,{\\left(c^2 +\\sqrt{-{\\left(c^2 +a\\,r_o \\right)}\\,{\\left(-c^2 +a\\,r_o \\right)}}\\right)}}{a^2 }}\\\\\n\\sqrt{\\frac{2\\,{\\left(c^2 -\\sqrt{-{\\left(c^2 +a\\,r_o \\right)}\\,{\\left(-c^2 +a\\,r_o \\right)}}\\right)}}{a^2 }}\n\\end{array}\\right)"}}
%---
%[output:8df05cfd]
%   data: {"dataType":"symbolic","outputData":{"name":"t_c","value":"\\sqrt{\\frac{2\\,{\\left(c^2 -\\sqrt{-{\\left(c^2 +a\\,r_o \\right)}\\,{\\left(-c^2 +a\\,r_o \\right)}}\\right)}}{a^2 }}"}}
%---
%[output:464dd6d8]
%   data: {"dataType":"symbolic","outputData":{"name":"F_netc","value":"-\\frac{2\\,\\sqrt{2}\\,K\\,Q^2 \\,a}{49\\,c^2 \\,r_o }"}}
%---
%[output:8e44930c]
%   data: {"dataType":"symbolic","outputData":{"name":"U_p","value":"\\frac{K\\,Q^2 }{49\\,r_o }"}}
%---
%[output:8df842b5]
%   data: {"dataType":"symbolic","outputData":{"name":"U_1","value":"\\frac{K\\,Q^2 }{98\\,r_o }"}}
%---
%[output:87d4e21a]
%   data: {"dataType":"symbolic","outputData":{"name":"U_2","value":"\\frac{\\sqrt{2}\\,K\\,Q^2 }{98\\,r_o }"}}
%---
%[output:83cdff46]
%   data: {"dataType":"symbolic","outputData":{"name":"U_net","value":"\\frac{3\\,K\\,Q^2 \\,{\\left(4\\,\\sqrt{2}+5\\right)}}{98\\,r_o }"}}
%---
%[output:0e3ea66d]
%   data: {"dataType":"symbolic","outputData":{"name":"qeU","value":"U_{\\textrm{net}} -\\frac{3\\,K\\,Q_2 \\,{\\left(4\\,\\sqrt{2}+5\\right)}}{98\\,r_o }=0"}}
%---
%[output:73c1cb3e]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\frac{98\\,U_{\\textrm{net}} \\,r_o }{3\\,K\\,{\\left(4\\,\\sqrt{2}+5\\right)}}"}}
%---
%[output:9591f1ac]
%   data: {"dataType":"symbolic","outputData":{"name":"F_netc","value":"-\\frac{2\\,\\sqrt{2}\\,K\\,Q^2 \\,a}{49\\,c^2 \\,r_o }"}}
%---
%[output:8b762e00]
%   data: {"dataType":"symbolic","outputData":{"name":"F_uc","value":"\\frac{U_{\\textrm{net}} \\,a\\,{\\left(\\frac{20\\,\\sqrt{2}}{21}-\\frac{32}{21}\\right)}}{c^2 }"}}
%---
%[output:96f77fb4]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-0.17693946440657614399839169122886\\,U_{\\textrm{net}} \\,a"}}
%---
