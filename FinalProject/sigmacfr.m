function disp2 = sigmacfr(s0)
%sigmacfr accetta come input un vettore s0 contenente diverse sigma iniziali di
%gaussiane, e restituisce la dispersione sulle posizioni di cui plotta il
%grafico in funzione di t.
%Il calcolo è effettuato esclusivamente per un oscillatore armonico con
%funzioni d'onda iniziali gaussiane centrate in x = -3.

N = 3e2;
a = 0.1;
x = -(N-1)/2:(N-1)/2; 
x = a*x';
V = x.^2/2;
dt = 0.01;
ts = 10;
nt = round(ts/dt);
T = dt:dt:nt*dt;
disp2 = zeros(nt,length(s0));
H = myhamiltonian1D(N,a,V,[],'DBC');
U = expm(-1i*H*dt);
for j = 1:length(s0)
    psi = exp(-(x+3).^2/s0(j)^2/2);
    psi = psi/norm(psi);
    for i = 1:nt
        psi = U*psi;
        xmean = real(psi'*diag(x)*psi);
        xmean2 = real(psi'*diag(x.^2)*psi);
        disp2(i,j) = xmean2-xmean^2;
    end
end

figure; clf
plot(T,disp2(:,1))
legend(['$\sigma(0)$ = ' num2str(s0(1))],'Location','northeast','Interpreter',...
    'latex','FontSize',12);
hold on
for j = 2:length(s0)
    plot(T,disp2(:,j),'Displayname',['$\sigma(0)$ = ' num2str(s0(j))]);   
end  
xlabel('t');
ylabel('$\Delta x^2$','Interpreter','latex');





