function out = swinger_lyapexp(setup)
% computes Lyapunov spectra of a dynamical system, using QR factorization.
% Notice that QR factorization differs from Gram-Schmidt only by
% some sign and standard roundoff.
%
% K. Ramasubramanian, M.S. Sriram, "A comparative study of computation 
% of Lyapunov spectra with different algorithms", Physica D 139 (2000) 
% 72-86.

% Oltre ai coefficineti di Lyapunov, vengono salvati nella struttura di output
% anche la loro somma, media temporale, le traiettorie di theta1 e theta2 e il tempo.

%esempi di chiamata: 
%s=swinger_lyapexp per avere in s la strutt di def che puoi modificare
%s=swinger_lyapexp([]) fa partire il programma con strutt di def

narginchk(0,1) %narginchk(minArgs,maxArgs) throw error nel caso non sono risp quei lim
                %posso passare o zero argomenti o 1 argomento

%% prep1

default.odefun      = @swinger_plusJ;
default.odesolver   = @ode113;      %[t,y] = ode113(odefun,tspan,y0) ode113(@(t,y)y,[0 5],1)
                                    %[t,y] = ode113(odefun,tspan,y0,options)
                                    %opts = odeset('reltol',1.e-4,'abstol',1.e-6,'outputfcn',@odephas2);
default.odeparams   = [1 1 1 1];
default.abstol      = 1e-14;
default.reltol      = 1e-12;
default.tstart      = 0;
default.tstep       = 0.002;
default.tstop       = 15;
default.ystart      = [pi/2 pi/2 0 0];
default.outflag     = 1;        %se voglio o no i valori dei coeff lyap stampati
default.outintrvl   = 5;

if nargin == 0  %out = swinger_lyapexp restituisce struttura out
    out = default;
    return;
end

if isempty(setup)   %out = swinger_lyapexp([]) fa partire il progr con valori di default
    setup = struct([]); %struct crea una struttura
end

% fill from 'setup' all new values, keeping the same ordering of 'default'
% ridefinisce default da setup e poi la chiamo s
fnames = fieldnames(setup);
for j = 1:length(fnames)
    fname = fnames{j};
    default.(fname) = setup.(fname);
end
s = default;

%% prep2

% number of degrees of freedom
n = length(s.ystart);

%  number of steps
niter = round((s.tstop-s.tstart)/s.tstep); %round arrotonda all'intero piu vicino

% initialization
t = s.tstart;
y = eye(n);     %mat identita nxn
y = [s.ystart'; y(:)]; %y(:) crea un vettore colonna con tutte le colonne di y impilate sotto la prima in ordine
                        %y rappresenta le condizioni iniziali del problema
                        %dell'ode
T = zeros(niter,1);
L = zeros(niter,n);
uL = zeros(1,n);
path1 = zeros(niter,1);
path2 = zeros(niter,1);
opt = odeset('Abstol',s.abstol,'Reltol',s.reltol);

%% main

for j = 1:niter
    T(j) = t;
    [~,y] = s.odesolver(s.odefun,[t t+s.tstep],y,opt,s.odeparams);  %la tilde annula il ritorno del tempo
    t = t+s.tstep;      %la y è la sol dell'eq diff quindi y(1) è theta1 y(2) theta2 diversi da y1 e y2
    y1 = y(end,1:n)';   %def y1 e y2 per calcolare i coeff di lyap
                        %16 colonne e indeterminate righe, prendo l'ultima
                        %riga, i primi 4 elementi
    y2 = reshape(y(end,n+1:n*(n+1))',n,n);%prendo i rim 16 e la metto in matrice
    [Q,R] = qr(y2);     
    for q = 1:n
        uL(q) = uL(q) + log(abs(R(q,q)));
    end
    L(j,:) = uL/(t-s.tstart);
    if s.outflag && mod(j,s.outintrvl) == 0
        fprintf('% 4.2f\t\t',t);
        fprintf('%+.8f\t',L(j,:));
        fprintf('\r');
    end
    y = [y1; Q(:)]; %trasforma Q in colonna e rimette y per il for
    path1(j) = y(1);
    path2(j) = y(2);
end

out.sum = sum(L,2);
out.mean = mean(L);
out.path1 = path1;
out.path2 = path2;
out.T = T;          %vettore dei tempi
out.L = L;          %matrice per cui ogni riga coeff di lyap ai diversi tempi
out.setup = s;      %campo della struttura di output che ha al suo interno la strutt di partenza
