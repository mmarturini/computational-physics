function g = myinverse(f,ab,N)
% g = inverse(f,ab,N) returns the griddedInterpolant object g that 
% approximately inverts the function f, input as a function handle,
% over the monotonicity interval ab=[ab(1) ab(2)], using N reference 
% points for reverse interpolation. In other words, g(f(x)) ~= x for 
% any x in ab and f(g(x)) ~= x for any x in f(ab).
% g = inverse(f,ab) defaults to g = inverse(f,ab,1e4).

% default N
if nargin < 3 || isempty(N)
    N = 1.e4;
end

% aggiunta per input f char
if ischar(f)
    f = vectorize(f);
    if contains(f,'(')
        f = str2func(['@(x)' f]);
    else
        f = str2func(f);
    end
end

% inizializzazione
a = ab(1);
b = ab(2);
y = (a:(b-a)/(N-1):b);

% controllo dominio funzione: se f non è definita o non è reale
% in un punto xo a partire da b, considera solo intervallo [xo,b]
for i=N:-1:1
    if ~isfinite(f(y(i))) || ~isreal(f(y(i)))
        a = y(i+1);
        y = (a:(b-a)/(N-1):b);
        break
    end
end    

% controllo monotonia funzione: se f non è monotona
% in un punto xo a partire da b, considera solo intervallo [xo,b]
for i=N:-1:3
    diff1 = f(y(i))-f(y(i-1));
    diff2 = f(y(i-1))-f(y(i-2));
    if sign(diff1) ~= sign(diff2)
        a = y(i-1);
        y = (a:(b-a)/(N-1):b);
        break
    end
end

% ordina y al contrario per funzioni decrescenti
% x = f(y) deve essere un vettore crescente
if f(b) < f(a)
    y = (b:(a-b)/(N-1):a);
end

x = f(y);
% per avere un controllo grafico del risultato si esegue il plot della
% funzione inversa nell'intervallo considerato
plot(x,y);
g = griddedInterpolant(x,y); 
