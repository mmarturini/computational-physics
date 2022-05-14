function [M,cp,kh2] = mymultidiagL(p)
% [M,cp,kh2] = mymultidiagL(p) returns the matrix M relevant for the 
% determination of the multidiagonal 1D Laplacian, the corresponding
% coefficients cp (in symbolic form) and the symfun kh2(k) which
% characterizes the spectrum.

%costruzione M
M(1,1:p+1) = 2;
M(2:p+1,:)=(0:p).^(2*(1:p)');  

%costruzione cp
a = zeros(p+1,1);
a(2) = -1;
cp = (M\sym(a))';

%spettro
syms k ;
kh2 = symfun(2*cp*cos((0:p)'*k),k);

end