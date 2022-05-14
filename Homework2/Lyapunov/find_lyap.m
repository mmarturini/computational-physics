function res = find_lyap(theta,par)

% find_lyap calcola i coefficienti di Lyapunov per diversi angoli
% iniziali e per diversi parametri (masse, lunghezze) dell'equazione del doppio
% pendolo, richiamando più volte la funzione swinger_lyapexp. Gli argomenti theta e
% par sono rispettivamente un vettore contenente gli angoli e una matrice n x 4 
% contenente n linee di parametri della forma [m1 m2 L1 L2].
% L'output è una struttura contenente il tempo 'T' e un matrice 'out(j,k)', dove
% ogni elemento è a sua volta una struttura che ha per campi la matrice 'L' dei
% coefficienti, la loro media temporale 'mean', la loro somma a ogni istante'sum' e
% le traiettorie 'path1' e 'path2' per theta1 e theta2.

if nargin < 2
    par = [1 1 1 1];
end

if nargin < 1
    theta = pi/2;
end

s = swinger_lyapexp;
s.outflag = 0;
t0 = s.tstart;
dt = s.tstep;
tf = s.tstop;

res.theta = theta;
res.par = par;
res.T = (t0+dt:dt:tf);

%j indicizza il numero di angoli k il numero dei parametri 
%(k=1 vuol dire [1 1 1 1] k=2 [1 1 1 1; 2 2 2 2;] ...) 
for j = 1:length(theta)
    s.ystart = [theta(j) theta(j) 0 0];
       for k = 1:size(par,1)        %ritorna la lunghezza della dimensione 1 di par, cioè numero di colonne
           s.odeparams = par(k,:);  %k-esima riga di par
           o = swinger_lyapexp(s);    
           res.out(j,k).L = o.L;
           res.out(j,k).path1 = o.path1;
           res.out(j,k).path2 = o.path2;
           res.out(j,k).sum = o.sum;
           res.out(j,k).mean = o.mean;
       end
end

%l'output di ritorno è la struttura res divise in campi
%res.theta = theta;
%res.par = par;
%res.T;
%res.out; out è una matrice di dim length(theta)Xsize(par,1)
%di cui ogni elemento è rappresentato da una struttura con i campi visti
%sopra nel for

%esempio di chiamata s=find_lyap([pi/2 pi/4],[1 1 1 1]);

