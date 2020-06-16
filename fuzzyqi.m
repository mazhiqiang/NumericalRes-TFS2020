function varargout = fuzzyqi(cfg)

% Original author: Lucian Busoniu
% Modified by: Zhian Kuang

c = [3,2.4,0.7,0.95];
if nargin < 1, cfg = struct(); end;

Local_cfg.run = 0;                        % run learning
Local_cfg.resume = 0;                     % resume learning
Local_cfg.replay = 0;                     % replay learned policy
Local_cfg.approx = 0;                     % special mode: export Q-function approximator from result of fuzzyQI
Local_cfg.problem = '';                   % what problem to solve
Local_cfg.loadapprox = '';                % load approximator data from file
Local_cfg.datadir = '';                   % data dir
Local_cfg.datafile = 'fzqidata';          % save data to file

Local_cfg.gamma = 0.98;                   % discount factor
Local_cfg.eps = .01;                      % threshold for convergence
Local_cfg.maxiter = 1000;                 % max number of iterations for Q-iteration
Local_cfg.term = 'zero';                  % how to handle terminal states ('ignore' or 'zero' the next-state Q-values/BFs)
Local_cfg.serial = 0;                     % serial (asynchronous) or parallel (synchronous) Q-iteration
Local_cfg.storeact = 1;                   % whether to store activation (membership) vectors (-1 for auto)

Local_cfg.returnprecision = [];           % replay s.t. return can be evaluated with this precision
Local_cfg.interph = 0;                    % use interpolated (averaged) policy to replay
Local_cfg.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
Local_cfg.tend = 30;                      % end time for replay
Local_cfg.userew = [];                    % use different reward function

Local_cfg.truncatehist = [];              % whether history should be truncated at this time value (in seconds)
                                    % before outputting and using in plots
Local_cfg.plottarget = 'screen';          % 'screen', '', 'latex', or 'beamer'. If 'screen' figures will not be closed
Local_cfg.savetheta = 0;                  % save param history in stats
Local_cfg.savedir = '';                   % save figure in this directory
Local_cfg.savefig = '';                   % save figure under this name

Local_cfg.visualize = 0;                  % visualization level (0 = none, 1 = iteration-level)
Local_cfg.viscfg = struct;                % visualization config options
Local_cfg.verb = 5;                       % verbosity: the higher, the more detailed the messages displayed
Local_cfg.noplot = 0;                     % whether to suppress figure plots
Local_cfg.silent = 0;                     % suppress all output
Local_cfg.initstepdisp = .1;              % feedback every 10% of MDP init
Local_cfg.iterdisp = 10;                  % feedback after every 10 iterations
Local_cfg.itervis = 10;                   % visualize once every 10 iterations
Local_cfg.itersave = 25;                  % save after each 25 iterations

ELocal_cfg.fuzzy_params = {};             % extra parameters for problem calling in 'fuzzy' mode
ELocal_cfg.model_params = {};             % extra parameters for problem calling in 'model' mode

KEEPFIELDS = {'problem', 'gamma', 'storeact', 'serial'};

if ischar(cfg),
    cfg = str2cfg(cfg, [fieldnames(Local_cfg); fieldnames(ELocal_cfg)]);
end;
cfg = checkparams(cfg, ELocal_cfg);

try     
    if ~isempty(cfg.problem), cfg = checkparams(cfg, feval(cfg.problem, 'fuzzy', cfg.fuzzy_params{:})); end;
catch
end;

cfg = checkparams(cfg, Local_cfg);

if cfg.silent, cfg.verb = -Inf; cfg.noplot = 1; end;

cfg.noinit = (cfg.resume || cfg.replay || cfg.approx) && exist([cfg.datafile '.mat'], 'file');

if cfg.noinit,        % load data file, making sure that cfg is not overwritten
    cfg1 = cfg; kf = KEEPFIELDS;
    load(cfg.datafile);

    cfg = copyfields(cfg, cfg1, kf);
    KEEPFIELDS = kf; clear cfg1 kf;
    dispx(['Data loaded from [' cfg.datafile '].'], cfg.verb, 1);
end;

if ~cfg.noinit,
    model = feval(cfg.problem, 'model', cfg.model_params{:});
end;


if cfg.loadapprox,
    load(cfg.loadapprox, 'X', 'U', 'DIMS', 'XMFS', 'MDP');
    dispx(['Initialization (approximator and sample data) loaded from [' cfg.loadapprox '].'], cfg.verb, 1);
end;



if ~cfg.noinit && isempty(cfg.loadapprox),

    checkparams(cfg, [], {'xgrids', 'ugrids'});
    X = cfg.xgrids; U = cfg.ugrids;

    DIMS.p = length(X); DIMS.q = length(U);     % # of states and outputs
    DIMS.dimx = []; DIMS.dimu = [];             % dimensions of grids
    XMFS = {};
    for p = 1:DIMS.p,
        XMFS{p} = gen_mfs(X{p});
        DIMS.dimx(end+1) = length(X{p});
    end;
    for q = 1:DIMS.q,
        DIMS.dimu(end+1) = length(U{q});
    end;

    DIMS.N = prod(DIMS.dimx);
    DIMS.M = prod(DIMS.dimu);
end;


% -----------------------------------------------
% Q-iteration
if cfg.run || cfg.resume,

    if ~cfg.resume && isempty(cfg.loadapprox),

        % estimate the storage size required for membership degrees, assuming most next states
        % activate 2^p membership functions
        actstorage = DIMS.N * DIMS.M * 2^DIMS.p;
        if cfg.storeact < 0,        % auto

            cfg.storeact = actstorage < 1e8;
        end;
        % override storeact when the variable size would be too large
        if cfg.storeact && actstorage > 1e8, 
            cfg.storeact = 0; 
            dispx(['Disabling (overriding) membership storage, size too large:' num2str(actstorage)], cfg.verb, -1);
        end;
        if cfg.storeact,	
            dispx(['Computing MDP and membership data for ' num2str(DIMS.N*DIMS.M) ' (x,u) pairs...'], cfg.verb, 0);
        else
            dispx(['Computing MDP data for ' num2str(DIMS.N*DIMS.M) ' (x,u) pairs...'], cfg.verb, 0);
        end;

        t = cputime;
        % init MDP structures
        if cfg.storeact,  MDP.F = sparse(DIMS.N * DIMS.M, DIMS.N);
        else              MDP.F = zeros(DIMS.N, DIMS.M, DIMS.p);
        end;
        MDP.R = zeros(DIMS.N, DIMS.M);
        MDP.T = zeros(DIMS.N, DIMS.M);

        % iterate over (xi, uj) \in (X0, U0) and populate MDP structures
        prog = cfg.initstepdisp;           
        tab = dec2base(0:(2^DIMS.p-1), 2) - 47; roll =0 * DIMS.dimx;
        Xflat = flat(X); Uflat = flat(U);
        for i = 1:DIMS.N,
            for j = 1:DIMS.M,
                [xplus MDP.R(i, j) MDP.T(i, j)] = feval(model.fun, model, Xflat(:, i), Uflat(:, j),c);
                if cfg.storeact,
                    if ~MDP.T(i, j) || cfg.term(1) == 'i',
                        % compute & store membership degrees directly
                        [ind, mu] = mdegs_p(xplus, X, roll, DIMS.dimx, DIMS.p, tab);
                        MDP.F(i+(j-1)*DIMS.N, ind) = mu;
                    % else: next state is terminal, and handle terminal states = zero
                    % mdegs remain zero, i.e., Q-values of next state always zero
                    end;
                else    
                    MDP.F(i, j, :) = xplus;
                end;
            end;            % FOR over action discretization
            if any(any(MDP.T)) && cfg.term(1) ~= 'i' && ~cfg.storeact,
                error('Handling terminal states without storeact=1 is not implemented');
            elseif any(any(MDP.T)) && cfg.term(1) == 'i',
                dispx('WARNING! Terminal states encountered & will be ignored.', cfg.verb, -1);
            end;

            if i/DIMS.N > prog,     % progress feedback
                dispx([num2str(prog * 100) '% completed...'], cfg.verb, 2);
                prog = prog + cfg.initstepdisp;
            end;
        end;                % FOR over state samples
        
        % record how much time MDP data computation took (disregarding that progress display is
        % also counted here)
        qistats.tinit = cputime - t;

        dispx('done.', cfg.verb, 2);
        save(cfg.datafile);
        dispx(['Initialization data saved to [' cfg.datafile ']'], cfg.verb, 1);
    end;                    % IF need to initialize MDP
    
    % -----------------------------------------------
    % Do Q-iteration
    dispx('Performing fuzzy Q-value iteration...', cfg.verb, 0);
    
    % Q-iteration parameters on the config: cfg.gamma, cfg.eps, cfg.maxiter
    % init param vector and iteration index when not resuming
    if ~cfg.resume,
        qistats.delta = [];
        qistats.t = 0;
        theta = zeros(DIMS.N, DIMS.M);
        if cfg.savetheta, 
            qistats.theta = cell(cfg.maxiter+1, 1);     % also allow for theta_0
            qistats.theta{1} = theta;                   % save theta_0 on the stats
        end;
        k = 1;
    end;

    % init visualization config if needed
    if cfg.visualize,
        vcfg = cfg.viscfg;
        vcfg.gview = [];
        vcfg.fuzzyqiter = 1;
        % visualize initial state of the algorithm
        vcfg.ell = 0;
        [figh vcfg.gview] = feval(model.visualizefun, vcfg);
    end;        
    
    tab = dec2base(0:(2^DIMS.p-1), 2) - 47; roll =0 * DIMS.dimx;
    delta = inf;
    t = cputime;
    while k <= cfg.maxiter && delta > cfg.eps,       % main loop
        
        thetaold = theta;     
        % compute policy optimal in Q

        % loop over (xi, uj) samples
        if cfg.storeact,        % can use the activation matrix MDP.F
            if cfg.serial,      % iterate over samples (need to since we have to reuse theta at each sample)
                for i = 1:DIMS.N,
                    for j = 1:DIMS.M,
                        theta(i, j) = MDP.R(i, j) + cfg.gamma * max(MDP.F(i+(j-1)*DIMS.N, :) * theta, [], 2);
                    end;    % FOR j
                end;        % FOR i
            else                % can get away with a matrix operation
                % This is most likely the most efficient version
                theta = MDP.R + cfg.gamma * reshape(max(MDP.F * theta, [], 2), DIMS.N, DIMS.M);
            end;            % IF serial
        else                % only the states have been stored
            for i = 1:DIMS.N,
                for j = 1:DIMS.M,
                    [ind, mu] = mdegs_p(MDP.F(i, j, :), X, roll, DIMS.dimx, DIMS.p, tab);
                    if cfg.serial,
                        theta(i, j) = MDP.R(i, j) + cfg.gamma * max(mu' * theta(ind, :));
                    else
                        theta(i, j) = MDP.R(i, j) + cfg.gamma * max(mu' * thetaold(ind, :));
                    end;    % IF serial
                end;% FOR j
            end;    % FOR i
        end;        % IF cfg.storeact
        
        % compute infinity-norm of function (not matrix)
        delta = max(max(abs(theta - thetaold)));        

        % update stats
        qistats.t = qistats.t + (cputime - t);
        qistats.delta(end + 1) = delta;
        if cfg.savetheta, qistats.theta{k+1} = theta; end;     % save theta on stats
        
        % visual feedback of algorithm progress
        if ~mod(k, cfg.iterdisp), 
            dispx(['k=' num2str(k) ' iteration completed, delta=' num2str(delta)], cfg.verb, 2);
        end;
        % visualization
        if cfg.visualize && ~mod(k, cfg.itervis),
            vcfg.ell = k;
            [figh vcfg.gview] = feval(model.visualizefun, vcfg);
        end;
        
        % data backup
        if ~mod(k, cfg.itersave),
            save(cfg.datafile);
            dispx(['Intermediary data at k=' num2str(k) ' saved to [' cfg.datafile '].'], cfg.verb, 1);
        end;
        t = cputime;
        k = k + 1;
    end;        % while not converged and allowed more iterations

    if delta < cfg.eps,	dispx('Convergence detected. Algorithm stopped.', cfg.verb, 0);
    else                dispx(['maxiter=' num2str(cfg.maxiter) ' exhausted. Algorithm stopped'], cfg.verb, 0);
    end;
    
    % finalize visualizer
    if cfg.visualize,
         % make sure last iteration is visualized
        if vcfg.ell < k - 1,    vcfg.ell = k - 1;       % show last iteration
        else                    vcfg.fuzzyqiter = 0;    % don't show anything
        end;
        vcfg.finalize = 1;
        [figh vcfg.gview] = feval(model.visualizefun, vcfg);
    end;
    
    % output optimal param and Q-iteration statistics
    varargout = {theta, qistats};
end;


% -----------------------------------------------
% Replay
if cfg.replay,

    % use different reward function
    if ~isempty(cfg.userew),
        model = feval(cfg.problem, 'changerew', model, cfg.userew);
    end;
    
    % compute locally optimal policy (local optimal discrete action for each basis
    % function/tile)
    if cfg.interph
        [Qstar ui] = max(theta, [], 2); clear Qstar;
        ui = lin2ndi(ui, DIMS.dimu);
        hstar = zeros(DIMS.N, DIMS.q);
        for q = 1:DIMS.q,
            hstar(:, q) = U{q}(ui(:, q));
        end;
    end;

    % initial state
    if ~isempty(cfg.x0),      % specified initial state
        x0 = cfg.x0(:);
    else                      % zeros
        x0 = zeros(p, 1);
    end;
    dispx(['Controlling from x0=' num2str(reshape(x0, 1, [])) ], cfg.verb, 0);

    % history
    if ~isempty(cfg.returnprecision),
        K = ceil(log(cfg.returnprecision * (1-cfg.gamma) / model.maxr) / log(cfg.gamma));
        cfg.tend = K * model.Ts;
    end;
    t = 0 : model.Ts : cfg.tend;
    Ns = length(t)-1;       % number of samples / time instances at which control is applied
    x = zeros(DIMS.p, length(t)); x(:, 1) = x0;
    u = zeros(DIMS.q, length(t)); u(:, end) = NaN;
    r = zeros(1, length(t)); r(1) = NaN;

    % steps loop
    tab = dec2base(0:(2^DIMS.p-1), 2) - 47; roll =0 * DIMS.dimx;
    for k = 1:Ns,
        % compute mdegs of current state
        [ind, mu] = mdegs_p(x(:, k), X, roll, DIMS.dimx, DIMS.p, tab);
        
        % compute optimal action
        if cfg.interph,     % either interpolated
            u(:, k) = (mu' * hstar(ind, :))';
        else                % or crisp (ties broken randomly)
            Qa = mu' * theta(ind, :);
            ui = find(Qa == max(Qa)); ui = ui(ceil(rand * length(ui)));
            ui = lin2ndi(ui, DIMS.dimu);
            for q = 1:DIMS.q,
                u(q, k) = U{q}(ui(:, q));
            end;
        end;
        
        % apply to system
%         [x(:, k+1) r(k+1) term] = feval(model.fun, model, x(:, k), u(:, k));
        [x(:, k+1) r(k+1) term] = feval(model.funapplied, model, x(:, k), u(:, k),c);
        if term, Ns = k; u(:, k+1) = NaN; break; end;      % entered terminal state
    end;
    save('Qiteration.mat','u');
    R = discreturn(cfg, r, Ns, term);
    
    % plot history & optionally save figures
    hist.t = t(1:Ns+1); hist.x = x(:, 1:Ns+1); hist.u = u(:, 1:Ns+1); hist.r = r(1:Ns+1); hist.R = R(1:Ns+1);
    if ~isempty(cfg.truncatehist), hist = truncatehist(hist, cfg.truncatehist); end;

    if ~cfg.noplot,
        if isfield(model, 'plotfun'),
            figh = feval(model.plotfun, hist);
            if isempty(figh), figh = plothistory(hist); end;
        else
            figh = plothistory(hist);
        end;
        setfigprop(cfg);
        saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
    end;
    
    % output history and possibly figure handles
    fig = ~cfg.noplot && (strcmp(cfg.plottarget, 'screen') || isempty(cfg.plottarget));
    if fig,         varargout = {hist, figh};
    else            varargout = {hist, []};
    end; 
end;        % IF replay

if cfg.approx,
    acfg.xgrids = X;
    acfg.ugrids = U;
    approx = triang(model, acfg);
    varargout = {approx, theta};
end;

% -----------------------------------------------
% Backup data

if cfg.run || cfg.resume, 
    % save into the same directory as the problem
    if isempty(cfg.datadir),
        datadir = fileparts(which(cfg.problem));
    else
        datadir = cfg.datadir; 
        if (datadir(end) == '\'), datadir = datadir(1:end-1); end;
    end;
    % if default file name is used, generate a unique file name
    if strcmp(cfg.datafile, Local_cfg.datafile),
        c = fix(clock);
        for i = 1:length(c),
            cfg.datafile = [cfg.datafile '-' num2str(c(i), '%02d')];
        end;
        delete([Local_cfg.datafile '.mat']);
    elseif ~strcmp(cfg.savedir, pwd)
        delete([cfg.datafile '.mat']);   % need to save in a different directory, so delete anyway
    end;    % otherwise just overwrite
    cfg.datafile = [datadir '\' cfg.datafile];
    save(cfg.datafile);
    dispx(['Fuzzy Q-iteration finished. Data was saved to [' cfg.datafile '].'], cfg.verb, 1);
end;


% END fuzzyqi() RETURNING varargout =================================================================
