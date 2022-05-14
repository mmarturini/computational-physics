function moler_3_5
%moler_3_5 esegue tre plot di una mano dopo essere passati in coordinate
%polari, con origine nel centro inferiore del palmo.

%caricamento file contenenti i punti della mano
load x_moler_3_4.mat x;
load y_moler_3_4.mat y;

%passaggio in coordinate polari, xo e yo sono i punti medi dei due punti
%all'estremità inferiore del palmo
xo = 0.4766;
yo = 0.0894;
x = x - xo;
y = y - yo;
theta = atan2(y,x);
r = sqrt(x.^2 + y.^2);

clf reset
nexttile 
plot(theta,r);
title('Plot 1');

%ordinamento degli angoli theta in ordine crescente e
%riordinamento degli r per far corrispondere di nuovo ad ogni r il suo
%corrispondente angolo theta
[thetaf, ind]= sort(theta);
r = r(ind);
t = (thetaf(1):.005:thetaf(end))';
p = pchiptx(thetaf,r,t);
s = splinetx(thetaf,r,t);

nexttile
plot(thetaf,r,'o',t,[p s],'-');
title('Plot 2');
axis([0 3.2 0 1]);

nexttile
plot(x,y,'o',p.*cos(t),p.*sin(t),'-',s.*cos(t),s.*sin(t),'-');
title('Plot 3');
axis([-0.3 0.3 0 1]);

end

%Il plot 1 rappresenta semplicemente r in funzione di theta, di conseguenza
%i punti vengono collegati da segmenti. La figura risulta ribaltata poichè
%theta ordina i punti della mano in senso antiorario
%Il plot 2 risulta anch'esso ribaltato con la differenza di utilizzare
%splinetx e pchiptx per l'interpolazione, dopo aver riordinato theta.
%Infine il plot 3 riutilizza coseni e seni per le coordinate polari per 
%graficare il profilo della mano e dunque la figura non risulta ribaltata.
%I grafici ottenuti in moler_3_4 sono migliori rispetto a quelli
%ottenuti in questo programma,questo è dovuto al fatto che molti dati
%risultano avere un theta simile rendendo più difficile
%l'interpolazione,inoltre il plot è anche influenzato dalla scelta del
%punto (xo,yo) che non è un punto della mano ma è scelto
%approssimativamente come origine per il sistema di coordinate polari.