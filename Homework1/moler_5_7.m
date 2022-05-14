function moler_5_7
%moler_5_7 dopo aver generato 11 punti (contenuti nel vettore t), 
%valuta in questi la funzione erf, i cui valori nel vettore y.
%Viene eseguito un fit con polinomi che vanno dal grado 1 al 10,l'errore
%massimo tra i dati e l'n-esimo polinomio viene inserito nel vettore err1,
%dunque  err1(i) contiene l'errore massimo tra i dati e il polinomio con
%grado i-esimo.In modo analogo viene creato un array err2 contenente gli
%errori massimi tra i dati e i polinomi con esponenti solo dispari.
%Infine la variabile res contiene l'errore massimo tra i punti e la
%funzione definita in functionmoler_5_7.m

t = (0:0.1:1)';
y = erf(t);
err1 = zeros(10,1);

%fit con polinomi da grado 1 a 10
for n=1:10
     p = polyfit(t,y,n); 
     err1(n)=norm(polyval(p,t)-y,inf);
end
fprintf('err1 = \n');
disp(err1);

%fit con polinomi dispari di grado da 1 a 9
C = [];
err2= zeros(5,1);
for j= 1:2:9 
    for i = 1:2:j
        C = [C t.^i];
    end
    coeff = pinv(C)*y;
    z = C*coeff;
    err2(j) = norm(z-y,inf);
end
fprintf('err2 = \n');
disp(err2(1:2:9));

%inizializzaione coefficienti di functionmoler_5_7
lambda=[1 2 3 4 5];
%errore massimo tra i dati e la funzione
fprintf('res = %d\n',err3_moler_5_7(lambda,t,y));

end

%L'errore massimo tra i dati e i polinomi diminuisce all'aumentare del
%grado massimo sia nel primo caso che nel caso con polinomi dispari.
%L'errore massimo ottenuto invece con la funzione fornita risuta maggiore 
% rispetto a quelli ottenuti dai polinomi,bisogna tenere però conto del
% fatto che gli 11 punti generati si trovano tra 0 e 1 e dunque la proprietà
% di essere esponenzialmente soppressa non risulta significativa nel fit.
%La funzione sarà dunque migliore ad approssimare erf,in tutto il suo
%dominio, rispetto ai polinomi che divergono per grandi t ma per i primi
%punti questi risultano approssimare meglio la funzione erf.



