function out = deployment_problem(what, varargin)

% gamma = 0.98;

switch what
        
    case 'model'      
        % Note that hooks exist for different model types
        % However, for now, the model type has no effect 
        if nargin >= 2 && ~isempty(varargin{1})
            cfg = varargin{1};
        else
            cfg = struct;
        end;

        MODEL.type = 'default';
        MODEL.wrap = 1;                          % wrapping enabled (1 -- default) or disabled (0)
        MODEL.orientation = 'horiz';             % 'horiz' or 'vert', arm orientation
        % reward/goal
        MODEL.rewtype = 'lqr';
        MODEL.Q = diag([1.5 0.3 1 1]);
        MODEL.R = zeros(2, 2);
        % integration

        MODEL.Ts = .01;
        MODEL.odemethod = 'euler';
        MODEL.odesolver = 'euler';
        MODEL.odesteps = 5;              

        model = parseconfig(cfg, MODEL);
        
        model.det = 1;      % always deterministic
        model.p = 4;        % always 4 states
        model.q = 2;        % only 1 actions
        % set bounds in standard format
        model.maxx = [0.98; 1.2; 1.2; 1.2]; 
        model.maxu = [3;0]; 
        model.fun = @deployment_mdpstandalone;
        model.funapplied = @deployment_SMC;
        model.plotfun = @deployment_plot;
        
        % compute max reward
        switch model.rewtype
            case 'lqr'
                model.maxr = model.maxx' * model.Q * model.maxx + model.maxu' * model.R * model.maxu;
        end;
        
        out = model;
   
end;


