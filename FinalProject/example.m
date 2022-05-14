function [T,T0,k] = example(V0,b,factor)
% try first example(1,2);

if nargin < 3
    factor = 6;
end

V = @(x) V0*(x.^2<b^2);
%L = laplacian1D(1024,2,'DBC');
[T,k] = Tcoeff(V,1024,0.1);
%[T,k] = Tcoeff(V,1024,0.1);
%[T,k] = Tcoeff(V,8*2048,0.01,[]); % the best, up to now
%[T,k] = Tcoeff(V,4*2048,0,04,[]);
T0 = sqbarrierT_exact(k*b,2*V0*b^2);

I = k < factor*sqrt(2*V0);
k = k(I);
T = T(I);
T0 = T0(I);

clf
plot(k,T,'.-')
hold on
plot(k,T0,'.-')
hold off

