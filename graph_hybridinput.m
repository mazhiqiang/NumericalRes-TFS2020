clc;clear;close all
load Qiteration.mat;

Omega = 1e-3;
h = 10;
dure = 20;
num = dure/h/Omega;
odesteps = 5;
dt = h*Omega/odesteps;
x = zeros(4,num);

uhybrid = zeros(1,num);
s2 = zeros(1,num);
x(:,1)=[-0.98;0.69;0.116;0];

c = [3,2.4,0.7,0.95];
c1 = c(1);
c2 = c(2);
c3 = c(3);
beta1 = c(4);

for k = 1:num
    xplus = x(:,k);
    for i = 1:odesteps
        [xdot1,uplus,s] = rarm_trans_SMC(x(:,k), u(:,k),c1,c2,c3,beta1);
        xplus = xplus + dt *xdot1;
    end
    x(:,k+1) = xplus;
    u_hybrid(k) = uplus(1)*(1.17e-3)^2*3500*10;
    s2(k) = s;
end
t = (1:num)*h*Omega;

s2_max = max(abs(s2));


label_true = 1;
figure('position',[100 100 800 450]);
h1p = [0.08,0.105,0.9,0.88];
h2p =[0.35,0.65,0.6,0.3]; 
h1=axes('Position',h1p);

% h1=axes('box','on');

h1pxmax = 20;
h1pxmin = 0;
h1pymax = 0.3;
h1pymin = -0.01;

a = [0.1;1];
b = [-0.005;0.3];%Target location

c = [(h2p(1)-h1p(1))/h1p(3)*(h1pxmax-h1pxmin)+h1pxmin;
    (h2p(1)-h1p(1)+h2p(3))/h1p(3)*(h1pxmax-h1pxmin)+h1pxmin];
d = [(h2p(2)-h1p(2))/h1p(4)*(h1pymax-h1pymin)+h1pymin;
    (h2p(2)-h1p(2)+h2p(4))/h1p(4)*(h1pymax-h1pymin)+h1pymin];%Enlarged location

[sortsetA,indexA] = sort([a;c]);
[sortsetB,indexB] = sort([b;d]);
dislocation1 = [3,1,4,2];
dislocation2 = [1,3,2,4];
location_flag1 = [1;2;3;4];%At right side
location_flag2 = [3;4;1;2];%At left side
if indexA == location_flag1 
    pointA1 = [a(2),b(1)];
    pointA2 = [a(2),b(2)];
    pointB1 = [c(1),d(1)];
    pointB2 = [c(1),d(2)];
elseif indexA == location_flag2
    pointA1 = [a(1),b(1)];
    pointA2 = [a(1),b(2)];
    pointB1 = [c(2),d(1)];
    pointB2 = [c(2),d(2)];
else
    if b(2)<d(1)||b(2)==d(1)
        pointA1 = [a(1),b(2)];
        pointA2 = [a(2),b(2)];
        pointB1 = [c(1),d(1)];
        pointB2 = [c(2),d(1)]; 
    elseif b(1)>d(2)||b(1)==d(2)
        pointA1 = [a(1),b(1)];
        pointA2 = [a(2),b(1)];
        pointB1 = [c(1),d(2)];
        pointB2 = [c(2),d(2)]; 
    else
        warning('Parameters error');
    end
end
flag_num = 0;
num_k = 1;
for i = 1:num-1
    hold on;
        if s2(i)^2<0.05&&u_hybrid(i)<0.2874&&u_hybrid(i)>0
         h1h = plot(h1,t(i:i+1),u_hybrid(1,i:i+1),'Color','b','linewidth',1,'linestyle','-');
        if label_true == 1
            hold on;
            plot(h1,t(i),u_hybrid(1,i),'bo');
            flag_num(num_k) = i;
            num_k = num_k+1;
            hold on;
            label_true = 0;
        end
        

        else
            if label_true ==0
                plot(h1,t(i),u_hybrid(1,i),'rs');
                label_true= 1;
            end
        
        plot(h1,t(i:i+1),u_hybrid(1,i:i+1),'r--');
    end
end
hold on;
plot(h1,[a(1) a(1)],[b(1) b(2)],'k-.',...
    [a(2) a(2)],[b(1) b(2)],'k-.',...
    [a(1) a(2)],[b(1) b(1)],'k-.',...
    [a(1) a(2)],[b(2) b(2)],'k-.',...
    [pointA1(1) pointB1(1)],[pointA1(2) pointB1(2)],'k-.',...
    [pointA2(1) pointB2(1)],[pointA2(2) pointB2(2)],'k-.');

ylim(h1,[h1pymin,h1pymax]);
grid(h1);
xlabel(h1,'True anomaly (rad)','Fontsize',12,'FontName','Times New Roman');
ylabel(h1,'Input (N)','Fontsize',12,'FontName','Times New Roman');
h1h = plot(h1,t(flag_num(1)-1),u_hybrid(1,flag_num(1)-1),'r--',t(flag_num(1)),u_hybrid(1,flag_num(1)),'bo',t(flag_num(1)+1),u_hybrid(1,flag_num(1)+1),'b-');

flag_num = 0;
num_k = 1;
 h2 = axes('Position',h2p);
for i = 1:num-1
    hold on;
        if s2(i)^2<0.05&&u_hybrid(i)<0.2874&&u_hybrid(i)>0
         plot(h2,t(i:i+1),u_hybrid(1,i:i+1),'Color','b','linewidth',1,'linestyle','-');
        if label_true == 1
            hold on;
            plot(h2,t(i),u_hybrid(1,i),'bo');
            flag_num(num_k) = i;
            num_k = num_k+1;
            hold on;
            label_true = 0;
        end
        

        else
            if label_true ==0
                plot(h2,t(i),u_hybrid(1,i),'rs');
                label_true= 1;
            end
        
        plot(h2,t(i:i+1),u_hybrid(1,i:i+1),'r--');
    end
end

xlim(h2,a);
ylim(h2,b);
grid(h2);
set(h1,'Fontsize',12,'FontName','Times New Roman','box','on');
% legend(hl([1 2 3 4]),'Q-iteration reaching law','Switching point to SMC','SMC reaching law','Switching point to Q-iteration','location','Southeast');
legend(h1h([1 2 3]),'Q-iteration reaching law','Switching point to SMC','SMC reaching law','location','Southeast');



function [xdot,u_new,s]= rarm_trans_SMC(x, u,c1,c2,c3,beta1)

u_new = u+3;
xdot =zeros(4,1);
dx2 = (1+x(1))*((x(4)+1)^2-1+3*cos(x(3))*cos(x(3)))-u(1)-3;
dx4 = -2*x(2)*(1+x(4))/(x(1)+1)-3*cos(x(3))*sin(x(3));

u(1) = u(1);
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





