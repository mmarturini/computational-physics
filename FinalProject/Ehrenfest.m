function diff = Ehrenfest(x0,s0)
%Ehrenfest restituisce la differenza tra i valori di x medio calcolato a
%partire dalla funzione d'onda e la traiettoria classica, e ne fa il grafico
%in funzione di t nel caso dell'oscillatore armonico con una gaussiana come
%funzione d'onda iniziale.
%x0 = posizione di partenza; s0 = deviazione standard

N = 3e2;
a = 0.1;
x = -(N-1)/2:(N-1)/2; 
x = a*x';
V = x.^2/2;
dt = 0.01;
ts = 50;
nt = round(ts/dt);
T = dt:dt:nt*dt;
xmean = zeros(1,nt);
H = myhamiltonian1D(N,a,V,[],'DBC');
U = expm(-1i*H*dt);
psi = exp(-(x-x0).^2/s0^2/2);
psi = psi/norm(psi);
for i = 1:nt
    psi = U*psi;
    xmean(i) = real(psi'*diag(x)*psi);
end
diff = xmean-x0*cos(T);

figure(1); clf
plot(T,diff)


