A = [-0.703029,30.092,25.56759,33.066178,24.99735];
B = [108.4773,6.832514,6.09613,-11.363417,55.18696];
C = [-42.52157,6.793435,4.054656,11.432816,-33.69137];
D = [5.862788,-2.53448,-2.671301,-2.772874,7.948387];
E = [0.678565,0.082139,0.131021,-0.158558,-0.136638];

T = 700+273; % K
t = T/1000; % 1/K
Cp = A + B.*t + C.*t^2 + D.*t^3 + E./t^2 % J/molK

dHf = [-74.4,-241.83,-110.53,0,-393.91]