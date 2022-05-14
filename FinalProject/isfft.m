function f = isfft(tf)
% inverse symmetric fast-fourier
% f = isfft(tf)
% acts columnwise unless f is a row vector

[N,nc] = size(tf);
flag = isrow(tf);
if flag
    tf = tf.'; N = nc; nc = 1 ;
end

n = floor(N/2);
tf = tf([n+1:N,1:n],:);

n = floor((N-1)/2);
f = ifft(exp(-1i*pi*(1-1/N)*[0:n,n+1-N:-1]')*ones(1,nc) .* tf);

if flag
    f = f.';
end
