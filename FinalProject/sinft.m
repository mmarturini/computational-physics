function sf = sinft(f)
% sf = sinft(f) 
% sine-fourier transform with doubling;
% acts columnwise unless f is a row vector.

[N,nc] = size(f);
flag = isrow(f);
if flag
    f = f.'; N = nc; nc = 1 ;
end

g = zeros(2*N+2,nc);
g(2:N+1,:) = f;
g(N+3:end,:) = - flipud(f);    
g = fft(g);
sf = 1i*g(2:N+1,:)/sqrt(2*(N+1));
if isreal(f)
    sf = real(sf);
end

if flag
    sf = sf.';
end