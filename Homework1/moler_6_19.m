function moler_6_19
% moler_6_19 calcola l'integrale di f in seguito definita, con precisione 
% tra 1.e-3 e 1.e-13 con il metodo delle somme parziali descritto. 
% a) un calcolo eseguito con gli algoritmi di quadratura non 
%    converge dato il numero infinito di oscillazioni di f intorno a 0
% Infine si indaga l'uso dell'accelerazione di Aitken

f = @(x)(1./x)*cos(log(x)/x);
% errore di default
tol = 1.e-7;
x(1) = 1;
T(1) = 0;
I(1) = 0;
k = 2; err = 1;

while err > tol
    fz = @(x)log(x)+(k-1.5)*pi*x;
    x(k) = fzerotx(fz,[eps,x(k-1)]);
    T(k) = quadtx(f,x(k),x(k-1));
    I(k) = I(k-1)+(T(k)+T(k-1))/2;
    err = abs(I(k)-I(k-1));
    k = k+1;    
end

fprintf('n_count = %i\nresult = %.10d\nerror = %d\n', k-1,I(end),err);

% si valutano i primi 10 elementi dell'accelerazione di Aitken e si
% confrontano con i corrispondenti valori di T
Aitkend = zeros(12,1);
for k = 2:11
    Aitkend(k) = T(k+1)-((T(k+1)-T(k))^2)/(T(k+1)-2*T(k)+T(k-1));
end

% b) la successione delle accelerazioni di Aitken converge più 
%    velocemente a 0 rispetto alle T; questo vale per ogni successione
%    iterativa
fprintf('Confronto primi elementi T(k) e Aitkend(k):\n');
disp([T(2:11)' Aitkend(2:11)]);



