function f = icosft(cf)
% f = icosft(cf) 
% inverse cosine-fourier transform with doubling;
% acts columnwise unless f is a row vector.

[N,nc] = size(cf);
flag = isrow(cf);
if flag
    cf = cf.'; N = nc; nc = 1 ;
end

g = zeros(2*N,nc);
g(1:N,:) = (exp(1i*pi*(0:N-1)'/N/2)*ones(1,nc)).*cf*sqrt(2*N);
g(1,:) = sqrt(2)*g(1,:);
g(2*N:-1:N+2,:) = (exp(-1i*pi*(1:N-1)'/N)*ones(1,nc)).*g(2:N,:);    
g = ifft(g);
f = g(1:N,:);
if isreal(cf)
    f = real(f);
end

if flag
    f = f.';
end