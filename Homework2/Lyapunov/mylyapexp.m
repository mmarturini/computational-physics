function out = mylyapexp(setup)
% computes Lyapunov spectra of a dynamical system, using QR factorization.
% Notice that QR factorization differs from Gram-Schmidt only by
% some sign and standard roundoff.
%
% K. Ramasubramanian, M.S. Sriram, "A comparative study of computation 
% of Lyapunov spectra with different algorithms", Physica D 139 (2000) 
% 72-86.

narginchk(0,1)

%% prep1

default.odefun = @lorenz_plusJ;
default.odesolver = @ode113;
default.odeparams = [8/3 28 10];
default.abstol = 1e-14;
default.reltol = 1e-12;
default.tstart = 0;
default.tstep = 0.5;
default.tstop = 200;
default.ystart = [27 9 13];
default.outflag = 1;
default.outintrvl = 5;

if nargin == 0
    out = default;
    return;
end

if isempty(setup)
    setup = struct([]);
end

% fill from 'setup' all new values, keeping the same ordering of 'default'
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
niter = round((s.tstop-s.tstart)/s.tstep);

% initialization
t = s.tstart;
y = eye(n);
y = [s.ystart'; y(:)];
T = zeros(niter,1);
L = zeros(niter,n);
uL = zeros(1,n);
opt = odeset('Abstol',s.abstol,'Reltol',s.reltol);

%% main

for j = 1:niter
    T(j) = t;
    [~,y] = s.odesolver(s.odefun,[t t+s.tstep],y,opt,s.odeparams);
    t = t+s.tstep;
    y1 = y(end,1:n)';
    y2 = reshape(y(end,n+1:n*(n+1))',n,n);
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
    y = [y1; Q(:)];
end

out.T = T; 
out.L = L;
out.setup = s;
