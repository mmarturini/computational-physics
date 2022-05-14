function cf = cosft(f)
% cf = cosft(f) 
% cosine-fourier transform with doubling;
% acts columnwise unless f is a row vector.

[N,nc] = size(f);
flag = isrow(f);
if flag
    f = f.'; N = nc; nc = 1 ;
end

g = zeros(2*N,nc);
g(1:N,:) = f;
g(N+1:end,:) = flipud(f);    
g = fft(g);
cf = (exp(-1i*pi*(0:N-1)'/N/2)*ones(1,nc)).*g(1:N,:)/sqrt(2*N);
if isreal(f)
    cf = real(cf);
end
cf(1,:) = cf(1,:)/sqrt(2);

if flag
    cf = cf.';
end