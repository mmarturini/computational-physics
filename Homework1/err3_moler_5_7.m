function res = err3_moler_5_7(lambda,t,y)
%   [res] = functionmoler_5_7(lambda,t,y) produce un fit di una funzione,
%   lineare nei coefficienti lambda(j), con variabile indipendente t 
%   e variabile dipendente y. La funzione produce come output il valore
%   massimo della norma del residuo.

m = length(t);
n = length(lambda);

%creazione design matrix
X = zeros(m,n);
X(:,1) = lambda(1);
X(:,2) = lambda(2)*exp(-t.^2);
X(:,3) = lambda(3)*(1/(1+t))*exp(-t.^2);
X(:,4) = lambda(4)*((1/(1+t)).^2)*exp(-t.^2);
X(:,5) = lambda(5)*((1/((1+t))).^3)*exp(-t.^2);

coeff= pinv(X)*y;
z = X*coeff;
res = norm(z-y,inf);
return

