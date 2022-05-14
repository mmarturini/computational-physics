function moler_2_15
%%
%(a) al crescere della dimensione della matrice di golub
%il numero di condizionamento aumenta esponenzialmente, come si evince dal
%grafico in scala logaritmica

clf reset
for n=1:10
    sum = 0;
    for i=1:100
        a = condest(golub(n));
        sum = sum+a;
    end
    mean = sum/100;
    semilogy(n,mean,'o');
    hold on;
end

%%
%(b)
%con il diagonal pivoting, si nota che gli elementi sulla diagonale di U
%sono tutti 1 e gli elementi di U in generale sono numeri interi, ma non
%piccoli, perciò non è evidente che golub(n) sia mal condizionata

%la porzione di codice per rilevare l'andamento è stata lasciata come
%commento
%{
for n=2:8
    lugui(golub(n))
end
%}

%%
%(c)
%si osserva che fino a circa n=10, il determinante viene 1.
%dalla matrice di dimensione 10 in poi non da risultati stabili
%questo perchè det usa la fattorizzazione LU per il calcolo del determinante
%che è sensibile a errori di floating point, che aumentano con n perchè
%la matrice è mal condizionata

clf reset
for n=1:20
    sum = 0;
    for i=1:100
        a = det(golub(n));
        sum = sum+a;
    end
    mean = sum/100;
    semilogy(n,mean,'o')
    hold on;
end