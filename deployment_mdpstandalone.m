function [xplus, rplus, terminal] = deployment_mdpstandalone(m, x, u,c)

u = max(-m.maxu, min(m.maxu, u));

    dt = m.Ts / m.odesteps; xplus = x;
    for i = 1:m.odesteps
        xplus = xplus + dt * deployment_transstandalone(x, u);
    end;

if m.wrap  
    xplus(1) =  max(-0.98, min(0.98, xplus(1)));
    xplus(2) =  max(-1.2, min(1.2, xplus(2)));
    xplus(3) = max(-1.2, min(1.2, xplus(3)));
    xplus(4) = max(-1.2, min(1.2, xplus(4)));
else        % no wrapping, just bound everything
    xplus = max(-m.maxx, min(m.maxx, xplus));
end;

% Compute reward
if m.rewtype(1) == 'l'           % LQR (only supported type)
    c1 = c(1);
    c2 = c(2);
    c3 = c(3);
    beta1 = c(4);

f1 = -2*x(2)/(x(1)+1)*(1+x(4))-3*cos(x(3))*sin(x(3))-c3*sign(x(1))*abs(x(1))^beta1/(x(1)+1); 
    s = c1*x(3)+c2*x(4)+f1;
    rplus = -1*(s)^2; 

end;

terminal = 0;       % task is continuing



function xdot = deployment_transstandalone(x, u)

xdot =zeros(4,1);
dx2 = (1+x(1))*((x(4)+1)^2-1+3*cos(x(3))*cos(x(3)))-u(1)-3;
dx4 = -2*x(2)*(1+x(4))/(x(1)+1)-3*cos(x(3))*sin(x(3));

xdot(1) = x(2);
xdot(2) = dx2;
xdot(3) = x(4);
xdot(4) = dx4;
