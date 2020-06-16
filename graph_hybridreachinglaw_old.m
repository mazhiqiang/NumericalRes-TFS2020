clc;clear;close all
load Qiteration.mat;

Omega = 1e-3;
h = 10;
dure = 1.01;
num = dure/h/Omega;
odesteps = 5;
dt = h*Omega/odesteps;
x = zeros(4,num);
% u = ans.u;   
uhybrid = zeros(1,num);
s2 = zeros(1,num);
x(:,1)=[-0.98;0.69;0.116;0];

c = [3,2.4,0.7,0.95];
c1 = c(1);
c2 = c(2);
c3 = c(3);
beta1 = c(4);

for k = 1:num-1
    xplus = x(:,k);
    for i = 1:odesteps
        [xdot1,uplus,s] = rarm_trans_SMC(x(:,k), u(:,k),c1,c2,c3,beta1);
        xplus = xplus + dt *xdot1;
    end
    x(:,k+1) = xplus;
    u_hybrid(k) = uplus(1);
    s2(k) = s;
end
t = (1:num)*h*Omega;

s2_max = max(abs(s2));
x(1,:) = (x(1,:)+1)*3500;
figure;
area(t,s2.^2,'facecolor',[0 0.8 0.8],'edgecolor','w');
alpha(0.6);
label_true = 1;
for i = 2:num-1
    hold on;
        if s2(i)^2<0.05&&u_hybrid(i)<6&&u_hybrid(i)>0
        plot(t(i-1:i),x(1,i-1:i),'Color','b','linewidth',1,'linestyle','-');
        if label_true == 1
            hold on;
            plot(t(i),x(1,i),'bo');
            hold on;
            label_true = 0;
        end
        

        else
            if label_true ==0
                plot(t(i),x(1,i),'rs');
                label_true= 1;
            end
        
        plot(t(i-1:i),x(1,i-1:i),'Color',[1-sqrt(abs(s2(i))/s2_max),sqrt(abs(s2(i))/s2_max),0],'linewidth',1,'linestyle','-');
    end
end

text(0.7,700,'Switching points','Fontsize',12,'FontName','Times New Roman')    
text(0.02,50,'\rho(x)=s^2(k)','Fontsize',12,'FontName','Times New Roman')
color_num = num;
mymap = zeros(color_num,3);
mymap(:,1) = 1;
for i = 2:num-1
    mymap(i-1,1) =  1-sqrt(abs(s2(i))/s2_max);
    mymap(i-1,2) =  sqrt(abs(s2(i))/s2_max);
    
    if s2(i)^2<0.05&&u_hybrid(i)<6&&u_hybrid(i)>0
        mymap(i-1,1) =  0;
        mymap(i-1,2) =  0;
         mymap(i-1,3) =  1;
    end
end
mymap(i,1) =  0;
        mymap(i,2) =  0;
         mymap(i,3) =  1;

mymap(i+1,1) =  0;
        mymap(i+1,2) =  0;
         mymap(i+1,3) =  1;
colormap(mymap);

xlim([0.01,1]);
ylim([0,1000]);
grid on;
colorbar('YTickLabel',{'0 rad','0.1 rad','0.2 rad',...
    '0.3 rad','0.4 rad','0.5rad','0.6 rad','0.7 rad','0.8 rad','0.9 rad','1 rad'},'location','southoutside')
xlabel(' True anomaly (rad)');
ylabel('Deployed Length (m)');

set(gca,'Fontsize',12,'FontName','Times New Roman','box','on','position',[0.15,0.28,0.80,0.7]);

function [xdot,u_new,s]= rarm_trans_SMC(x, u,c1,c2,c3,beta1)

u_new = u+3;
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
    
if s^2<0.05&&x(4)>-1

    lambda = 40;   
        
    ueq = -1/p12*(c1*x(4)+c2*f1+p11*x(2)+p12*fu+p13*x(4)+p14*f2);
    usw = -1/p12*(lambda*s);

    u_smc = -(ueq+usw);
    if u_smc > 0 && u_smc <6 && x(3)>-1
         xdot =zeros(4,1);
    dx2_smc = (1+x(1))*((x(4)+1)^2-1+3*cos(x(3))*cos(x(3)))-u_smc;
    dx4_smc = -2*x(2)*(1+x(4))/(x(1)+1)-3*cos(x(3))*sin(x(3));
    xdot(1) = x(2);
    xdot(2) = dx2_smc;
    xdot(3) = x(4);
    xdot(4) = dx4_smc;
    u_new = u_smc;
    end

end

end





