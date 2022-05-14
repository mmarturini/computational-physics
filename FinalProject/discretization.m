N=1e3;
k=2*pi*(-N/2:N/2-1)'/N;
plot(k,k.^2,'.-')
xlim([0 3.5])
hold
kh2=@(k)(2*sin(k/2)).^2;
plot(k,kh2(k),'.-')

[m,cp,kh2]=mymultidiagL(2);
plot(k,double(kh2(k)),'.-')
%taylor(kh2,'order',8)	%per verificare che la prima correz avviene a k^6 come da teoria

[m,cp,kh2]=mymultidiagL(5);
plot(k,double(kh2(k)),'.-')
%taylor(kh2,'order',10)

legend('continuo','p = 1','p = 2','p = 5')