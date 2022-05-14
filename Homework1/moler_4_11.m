function moler_4_11
% moler_4_11 trova n massimo per cui n! viene calcolato 
% con le funzioni prod e gamma:
% a) esattamente in doppia precisione, confrontando il risultato in doppia
%   precisione con quello simbolico;
% b) senza overflow;

diff = 0; n = 1;
while(diff < eps) 
    vd = prod(1:n);
    vs = prod(1:sym(n),'native');
    diff = vd-vs;
    n = n+1;
end
disp(['n_max(prod) = ' num2str(n-2)]);

%%
diff = 0; n = 0;
while(diff < eps) 
    vd = gamma(n+1);
    vs = gamma(sym(n+1));
    diff = vd-vs;
    n = n+1;
end
disp(['n_max(gamma) = ' num2str(n-2)]);

%%
v = 0; n = 1;
while(v < realmax)
    v = prod(1:n);
    n = n+1;
end
disp(['n_max(prod) = ' num2str(n-2)]);

%%
v = 0; n = 0;
while(v < realmax)
    v = gamma(n+1);
    n = n+1;
end
disp(['n_max(gamma) = ' num2str(n-2)]);
    