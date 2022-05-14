function moler_4_13
% moler_4_13 contiene la funzione y = gammalninv(x)
% e richiede di inserire manualmente il valore di x 
% per cui calcolarla (poichè non sono richiesti
% input espliciti). Infine stampa a schermo il risultato

x = input('Inserire valore di x: ');
y = gammalninv(x);
disp(['y = ' num2str(y)]);

function y = gammalninv(x)
% y = gammalninv(x) calcola l'inversa del logaritmo della 
% funzione gamma, ossia y tale che gammaln(y) = x

% ymin è l'estremo sinistro dell'intervallo di monotonia di
% gammaln(y). xmin e ymin definiscono i valori minimi di x e y
% per cui y = gammalninv(x) è definita in modo biunivoco
ymin = fminsearch(@(y)gammaln(y),1.5);
xmin = gammaln(ymin);

if x < xmin
    error('x deve essere nell intervallo di monotonia (x > -0.12)');
end

f = @(y) gammaln(y)-x;
y = fzerotx(f,[ymin,x+8]);
% l'intervallo [ymin, x+8] permette di individuare sempre uno 
% zero di f, poichè gammaln(y)-y non cambia segno prima di 7.5

end

end
