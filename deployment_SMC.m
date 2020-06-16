function [xplus, rplus, terminal] =deployment_SMC(m, x, u,c)
% Discrete-time dynamics of the deployment
%  [XPLUS, RPLUS, TERMINAL] = RARM_MDP(M, X, U)
% Standalone, new version (does not require the m/rarm directory).
%


% limit torque
    c1 = c(1);
    c2 = c(2);
    c3 = c(3);
    beta1 = c(4);
    
u = max(-m.maxu, min(m.maxu, u));
dt = m.Ts / m.odesteps; xplus = x;
for i = 1:m.odesteps
    xplus = xplus + dt * deployment_trans_SMC(x, u,c1,c2,c3,beta1);
end;


% Normalized state
if m.wrap
    xplus(1) =  max(-0.98, min(0.98, xplus(1)));
    xplus(2) =  max(-1.5, min(1.5, xplus(2)));
    xplus(3) = max(-1.5, min(1.5, xplus(3)));
    xplus(4) = max(-1.5, min(1.5, xplus(4)));
else        % no wrapping, just bound everything
    xplus = max(-m.maxx, min(m.maxx, xplus));
end;

% Compute reward

    f1 = -2*x(2)/(x(1)+1)*(1+x(4))-3*cos(x(3))*sin(x(3))-c3*sign(x(1))*abs(x(1))^beta1/(x(1)+1); 
     s = c1*x(3)+c2*x(4)+f1;
    rplus = -1*(s)^2; 

terminal = 0;       % task is continuing



function xdot = deployment_trans_SMC(x, u,c1,c2,c3,beta1)



xdot =zeros(4,1);
dx2 = (1+x(1))*((x(4)+1)^2-1+3*cos(x(3))*cos(x(3)))-u(1)-3;
dx4 = -2*x(2)*(1+x(4))/(x(1)+1)-3*cos(x(3))*sin(x(3));

xdot(1) = x(2);
xdot(2) = dx2;
xdot(3) = x(4);
xdot(4) = dx4;
 
f1 = -2*x(2)/(x(1)+1)*(1+x(4))-3*cos(x(3))*sin(x(3))-c3*sign(x(1))*abs(x(1))^beta1/(x(1)+1); 
fu = (1+x(1))*((x(4)+1)^2-1+3*cos(x(3))*cos(x(3)));
    f2 = (x(1)+1)*((1+x(4))^2-1+3*cos(x(3))*cos(x(3)));
    p13 = 3*sin(x(3))*sin(x(3))-3*cos(x(3))*cos(x(3));
    p14 = -2*x(2)/(x(1)+1);
    p11 = 2*x(2)/(x(1)+1)/(x(1)+1)*(1+x(4))-beta1*c3*sign(x(1))*abs(x(1))^(beta1-1)/(x(1)+1)+beta1*c3*sign(x(1))*abs(x(1))^(beta1-1)/(x(1)+1)^2;
    p12 = -2/(x(1)+1)*(1+x(4));
    
    s = c1*x(3)+c2*x(4)+f1;
    
if s^2<0.1&&x(4)>-1

    lambda = 40;   
    
    
     ueq = -1/p12*(c1*x(4)+c2*f1+p11*x(2)+p12*fu+p13*x(4)+p14*f2);
     usw = -1/p12*(lambda*s);

    u_smc = -(ueq+usw);
    if u_smc > 0 && u_smc <6
         xdot =zeros(4,1);
    dx2_smc = (1+x(1))*((x(4)+1)^2-1+3*cos(x(3))*cos(x(3)))-u_smc;
    dx4_smc = -2*x(2)*(1+x(4))/(x(1)+1)-3*cos(x(3))*sin(x(3));

    xdot(1) = x(2);
    xdot(2) = dx2_smc;
    xdot(3) = x(4);
    xdot(4) = dx4_smc;
    end

end





