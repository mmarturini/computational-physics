function out = randomwalk2(in)
% n-dimensional random walk
% in = randomwalk returns the default setup as a struct.
% out = randomwalk(in) returns the structure out with field 
% 'final_positions' and field 'in' containing the corresponding setup.
% out = randomwalk([]) runs with default settings;


% TODO implement variable step update algorithm?

narginchk(0,1)

%% set defaults (2-dim case)

default.number_of_particles = 1e4;
default.boundary_conditions = [0 0 0 0]; %s.boundary_conditions=[1 1 1 1];
default.sigma = 1;
default.initial_positions = @(N)0.01*randn(N,2); %matrice N*2 di random numbers da normale standard *0.01
                                                %[zeros(N,1) ones(N,1)]
default.drift_field = @(t,x) zeros(size(x));    %s.drift_field=@(t,x)-4*x; forza che attrae lungo x
default.random_jumps = @(t,x,s)normrnd(0,s,size(x)); %r = normrnd(mu,sigma) numero random da gauss con quei param
                      %s.random_jumps=@(t,x)2*rand(size(x))-1;
                      %s.random_jumps=@(t,x)0.1*tan(pi*(rand(size(x))-1/2)); lorentz-Cacuhy
default.time_step = 1e-3; % for fixed step algorithms
default.step_tolerance = 1e-6; % for variable step algorithms
default.time_span = 0.5;
default.step_algorithm = 'euler';
default.space_is_phase_space = 0;
default.plot.viewbox = [-1 1 -1 1];
default.plot.figure_number = 5;
default.plot.figure_position = [850 260 500 500];
default.plot.frame_skips = 0;
default.hist.varpair = 'xr';  % '' 'xy' 'xz' 'xr' ....
default.hist.nbin = 30;
default.hist.maxheight = [5 5];
default.parallel_procs = 1;
default.save_walks = true;    

% E' stato aggiunto il campo 'sigma' per scegliere l'ampiezza dei salti
% gaussiani.
% Posto il campo 'save_walks' = true, vengono salvati i cammini dei walkers in
% una matrice. Se questi sono ottenuti dall'opzione grafica, l'allocazione è dinamica
% e l'esecuzione può risultare lenta; per ovviare il problema, si possono utilizzare
% processi in parallelo, per cui la matrice è di dimensione nota e viene quindi
% preallocata

%per tre dimensioni
%s.boundary_conditions=[0 0 0 0 0 0];
%s.plot.viewbox=[-1 1 -1 1 -1 1];
%s.initial_positions=@(N)0.01*randn(N,3);

%% input handling

if nargin == 0
    out = default;
    return;
end
if isempty(in)
    in = struct([]);
end

% fill from 'in' all new values, keeping the same ordering of default
in = recsetup(in,default);

%% prep

% fill short-named variables and perform some consistency check
vbox = in.plot.viewbox;
dim = length(vbox)/2;   %di default [1 -1 1 -1] significa quadrato 2d quindi due dimensioni (x,y)
bc = zeros(1,2*dim); % default free boundary conditions
bc(1:length(in.boundary_conditions)) = in.boundary_conditions;

N = in.number_of_particles;
s = in.sigma;
nproc = max(1,in.parallel_procs);% numero proc in parallelo
Nperproc = floor(N/nproc); %Y = floor(X) rounds each element of X to the nearest integer less than or equal to that element
                            %numero di part per processo in parallelo
if N ~= nproc*Nperproc
    warning('number of particles not divisible by number of processes')
    N = nproc*Nperproc;
end
if isnumeric(in.initial_positions)
    x = in.initial_positions;   %in_pos è matrice N*2 prima colonna indica le x la seconda le y
                                %quindi un elemento indica posizione in in x e y della particella
else
    x = in.initial_positions(N);
end
if in.save_walks == true
    positions = x;
end

B = in.drift_field;
rj = in.random_jumps;
dt = in.time_step;
nt = round(in.time_span/dt);
if dim ~= 2 && dim ~= 3
    skip = inf; % do not plot in dimensions other that 2 or 3
else
    skip = in.plot.frame_skips;
end

%preallocazione matrice dei cammini dei random walkers nel caso si
%utilizzino processi in parallelo
if nproc > 1 && isfinite(nt) %nproc è massimo tra 1 e in.parallel_procs, nt sarebbe numero di step temporali
    x = reshape(x,Nperproc,nproc,dim); %x è un array di dimensioni Nperproc*nproc*dim
    x = permute(x,[1 3 2]); %A=rand(3,4,2) B=permute(A,[3 2 1]) B diventa una mat 2*4*3 
    if in.save_walks == true
        positions = zeros(Nperproc,dim,nproc,nt); %matrice di zeri Nperproc*dim*nproc*nt
        positions(:,:,:,1) = x;
    end
end

update = @update0;
if in.space_is_phase_space
    if mod(dim,2) ~= 0
        error('a classical phase space must have even dimensionality')
    end
    update = @update1;
end

onestep = in.step_algorithm;

%% main

% handle graphics
if skip < nt && nproc == 1
    hf = figure(in.plot.figure_number);
    clf(in.plot.figure_number,'reset');
    set(hf,'numbertitle','off','name','Random walk')
    set(hf,'position',in.plot.figure_position);
    ha = axes('outerposition',[0 .02 1 1]);
    if dim == 2
        h = plot(ha,x(:,1),x(:,2),'.');
    else
        h = plot3(ha,x(:,1),x(:,2),x(:,3),'.');
    end
    box('on'); axis(vbox); axis square; shg
    if ~isempty(in.hist.varpair)
        nbin = in.hist.nbin;
        set(hf,'position',in.plot.figure_position + [0 -250 0 250]);
        ha.OuterPosition = [0.1 .3 0.8 0.8];
        hb = axes('outerposition',[0 0.05 0.55 0.375]);
        hc = axes('outerposition',[0.45 0.05 0.55 0.38]);
        [xb,binb] = histplot(in.hist.varpair(1),x);
        [xc,binc] = histplot(in.hist.varpair(2),x);
        hxb = histogram(hb,xb,binb,'Normalization','pdf');
        hxc = histogram(hc,xc,binc,'Normalization','pdf');
        hc.YTickLabel = [];
        hb.YLim = [0 in.hist.maxheight(1)];
        hc.YLim = [0 in.hist.maxheight(2)];
    end
    ttl = title(sprintf('t = %-8.1f',0));
    run = uicontrol('style','toggle','string','start','value',0,...
        'units','normalized','position',[.375 .01 .1 .05]);
    stop = uicontrol('style','toggle','string','stop','value',0,...
                'units','normalized','position',[.55 .01 .1 .05]);
else
    run.Value = 1;
    stop.Value = 0;
end

if nproc > 1 && isfinite(nt)
    parfor jproc = 1:nproc     %parfor esegue ciclo for in parallelo
        for jt = 2:nt
            t = (jt-1)*dt;
            x(:,:,jproc) = update(t,x(:,:,jproc)); %#ok<PFBNS> %ignora il warning?
            if in.save_walks == true %#ok<PFBNS>
                positions(:,:,jproc,jt) = x(:,:,jproc);
            end
        end
    end
else
    %no parallelo
    jt = 0;
    while ~stop.Value && jt < nt
        if ~run.Value
            pause(0.2)
            if strcmp(run.String,'pause')
                set(run,'string','resume','callback','run.Value=1;')
            end
        else
            jt = jt+1;
            t = jt*dt;
            x = update(t,x);
            if in.save_walks == true
                positions = cat(2,positions,x); %concatena matrice x alla mat positions a destra(2=righe)
            end                                 %salva in positions i cammini in maniera dinamica
            if skip < nt
                if mod(jt+1,skip+1) == 0
                    h.XData = x(:,1);
                    h.YData = x(:,2);
                    if dim == 3
                        h.ZData = x(:,3);
                    end
                    if ~isempty(in.hist.varpair)
                        hxb.Data = histplot(in.hist.varpair(1),x);
                        hxc.Data = histplot(in.hist.varpair(2),x);
                    end
                    ttl.String = sprintf('t = %-5.2f',t);
                    drawnow
                end
                set(run,'string','pause','callback','run.Value=0;')
            end
        end
    end
end

%in parallelo
if nproc > 1
    if in.save_walks == true
        positions = permute(positions,[1 3 2 4]);
        positions = reshape(positions,N,dim*nt); 
    else
        %non sto salvando i cammini sto solo salvando le posiz finali
        x1=x; %x1 è il risultato di x = update(t,x);
        x = permute(x,[1 3 2]);
        x = reshape(x,N,dim); %questi due comandi mi permettono di salvare in x le posizioni finali dei random walker
    end                         %nel caso 2D matrice N*2 in cui prima colonna le x seconda le y
end

%salva nel campo positions la matrice delle posizioni se save è 1
%altrimenti solo le posizioni finali
if in.save_walks == true
    out.positions = positions;
else
    out.positions = x; 
end
out.in = in;

if skip < nt && nproc == 1
    set(stop,'string','close','value',0,'callback','close(gcf)')
end

%% nested funcs

    function [u,bin] = histplot(w,x)
        if w == 'x'
            u = x(:,1);
            bin = linspace(1.25*vbox(1),1.25*vbox(2),nbin);
        elseif w == 'y'
            u = x(:,2);
            bin = linspace(1.25*vbox(3),1.25*vbox(4),nbin);
        elseif w == 'z' && dim == 3
            u = x(:,1);
            bin = linspace(1.25*vbox(5),1.25*vbox(6),nbin);
        elseif w == 'r'
            u = sqrt(sum(x.^2,2));
            dbox = sqrt(sum(vbox(2:2:end).^2));
            bin = linspace(0,1.25*dbox,nbin);
        end
    end

    function x = update0(t,x)
        xnew = feval(onestep,B,t,x,dt);
        x = xnew + sqrt(dt)*rj(t,x,s);
        %x = xnew + sqrt(dt)*(rj(t,x) + rj(t+dt,xnew))/2;
        for jdim = 1:dim
            jo = 2*jdim-1;
            je = 2*jdim;
            if bc(jo)
                I = x(:,jdim) < vbox(jo);
                x(I,jdim) = -x(I,jdim) + 2*vbox(jo);
            end
            if bc(je)
                I = x(:,jdim) > vbox(je);
                x(I,jdim) = -x(I,jdim) + 2*vbox(je);
            end
        end
    end

    function x = update1(t,x)
        xnew = feval(onestep,B,t,x,dt);
        x = xnew + sqrt(dt)*rj(t,x,s);
        %x = xnew + sqrt(dt)*(rj(t,x) + rj(t+dt,xnew))/2;
        for jdim = 1:dim/2
            jo = 2*jdim-1;
            je = 2*jdim;
            if bc(jo)
                I = x(:,jdim) < vbox(jo);
                x(I,jdim) = -x(I,jdim) + 2*vbox(jo);
                x(I,jdim+dim/2) = -x(I,jdim+dim/2);
            end
            if bc(je)
                I = x(:,jdim) > vbox(je);
                x(I,jdim) = -x(I,jdim) + 2*vbox(je);
                x(I,jdim+dim/2) = -x(I,jdim+dim/2);
            end
        end
    end



end

%% aux functions

function yp = euler(F,t,y,h) %#ok<DEFNU>
    yp = y + h*F(t,y);
end

function yp = rk4(F,t,y,h) %#ok<DEFNU>
    s1 = F(t,y);
    s2 = F(t+h/2,y+h/2*s1);
    s3 = F(t+h/2,y+h/2*s2);
    s4 = F(t+h,y+h*s3);
    yp = y + h*(s1 + 2*s2 + 2*s3 + s4)/6;
end

function s = recsetup(a,b)
    s = b;
    for fname = fieldnames(a)'
        if isstruct(a.(fname{1}))
            s.(fname{1}) = recsetup(a.(fname{1}),b.(fname{1}));
        else
            s.(fname{1}) = a.(fname{1});
        end
    end
end