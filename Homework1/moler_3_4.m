function moler_3_4
% moler_3_4 interpola il profilo di una mano con i dati acquisiti tramite
% ginput, utilizzando le funzioni splinetx e pchiptx

load x_moler_3_4.mat x;
load y_moler_3_4.mat y;
n = length(x);
s = (1:n)';
t = (1:.05:n)';
clf reset

u = splinetx(s,x,t);
v = splinetx(s,y,t);
nexttile
plot(x,y,'.',u,v,'-');
title('Plot1 splinetx ');

e = pchiptx(s,x,t);
g = pchiptx(s,y,t);
nexttile
plot(x,y,'.',e,g,'-');
title('Plot2 pchiptx');

end

%il grafico con spline fornisce un risultato migliore di pchip, infatti la
%figura 3.11 di Moler è interpolata con spline