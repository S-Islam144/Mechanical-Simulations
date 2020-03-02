clc;
clear all;
close all;

%Sand_Crushing numerical solution coded by "Kh. Saidul Islam"

dt = 0.01;
tend = 50;
t = 0:dt:tend;
s = [0.4 0.3 0.3 0.2 0];
b = [0 0 0 0; ...
    0.6 0 0 0; ...
    0.2 0.4 0 0; ...
    0.1 0.4 0.3 0; ...
    0.1 0.2 0.7 1.0];
m1 = zeros(1,length(t));
m2 = zeros(1,length(t));
m3 = zeros(1,length(t));
m4 = zeros(1,length(t));
m5 = zeros(1,length(t));

m1(1) = 100;
m2(1) = 0;
m3(1) = 0;
m4(1) = 0;
m5(1) = 0;

f_ty_01 = @(t,m1)-s(1)*m1;
f_ty_02 = @(t,m1,m2)-s(2)*m2+s(1)*m1*b(2,1);
f_ty_03 = @(t,m1,m2,m3)-s(3)*m3+s(2)*m2*b(3,2)+s(1)*m1*b(3,1);
f_ty_04 = @(t,m1,m2,m3,m4)-s(4)*m4+s(3)*m3*b(4,3)+s(2)*m2*b(4,2)+s(1)*m1*b(4,1);
f_ty_05 = @(t,m1,m2,m3,m4,m5)-s(5)*m5+s(4)*m4*b(5,4)+s(3)*m3*b(5,3)+s(2)*m2*b(5,2)+s(1)*m1*b(5,1);

for ii = 1:length(t)-1
    % for finding the time required for 99.9% of m5
    %if m5(ii)>=99.9
        %re_time = t(ii-1);
        %break
    %end
    
    a1=f_ty_01(t(ii),m1(ii));
    a2=f_ty_01(t(ii)+(dt/2),m1(ii)+(dt*a1/2));
    a3=f_ty_01(t(ii)+(dt/2),m1(ii)+(dt*a2/2));
    a4=f_ty_01(t(ii)+(dt),m1(ii)+dt*a3);
    m1(ii+1)=m1(ii)+(dt/6)*(a1+(2*(a2+a3))+a4);
    
    b1=f_ty_02(t(ii),m1(ii),m2(ii));
    b2=f_ty_02(t(ii)+(dt/2),m1(ii)+(dt*b1/2),m2(ii)+(dt*b1/2));
    b3=f_ty_02(t(ii)+(dt/2),m1(ii)+(dt*b2/2),m2(ii)+(dt*b2/2));
    b4=f_ty_02(t(ii)+(dt),m1(ii)+dt*b3,m2(ii)+dt*b3);
    m2(ii+1)=m2(ii)+(dt/6)*(b1+(2*(b2+b3))+b4);
    
    c1=f_ty_03(t(ii),m1(ii),m2(ii),m3(ii));
    c2=f_ty_03(t(ii)+(dt/2),m1(ii)+(dt*c1/2),m2(ii)+(dt*c1/2),m3(ii)+(dt*c1/2));
    c3=f_ty_03(t(ii)+(dt/2),m1(ii)+(dt*c2/2),m2(ii)+(dt*c2/2),m3(ii)+(dt*c2/2));
    c4=f_ty_03(t(ii)+(dt),m1(ii)+dt*c3,m2(ii)+dt*c3,m3(ii)+dt*c3);
    m3(ii+1)=m3(ii)+(dt/6)*(c1+(2*(c2+c3))+c4);
    
    d1=f_ty_04(t(ii),m1(ii),m2(ii),m3(ii),m4(ii));
    d2=f_ty_04(t(ii)+(dt/2),m1(ii)+(dt*d1/2),m2(ii)+(dt*d1/2),m3(ii)+(dt*d1/2),m4(ii)+(dt*d1/2));
    d3=f_ty_04(t(ii)+(dt/2),m1(ii)+(dt*d2/2),m2(ii)+(dt*d2/2),m3(ii)+(dt*d2/2),m4(ii)+(dt*d2/2));
    d4=f_ty_04(t(ii)+(dt),m1(ii)+dt*d3,m2(ii)+ dt*d3,m3(ii)+dt*d3,m4(ii)+dt*d3);
    m4(ii+1)=m4(ii)+(dt/6)*(d1+(2*(d2+d3))+d4);
    
    e1=f_ty_05(t(ii),m1(ii),m2(ii),m3(ii),m4(ii),m5(ii));
    e2=f_ty_05(t(ii)+(dt/2),m1(ii)+(dt*e1/2),m2(ii)+(dt*e1/2),m3(ii)+(dt*e1/2),m4(ii)+(dt*e1/2),m5(ii)+(dt*e1/2));
    e3=f_ty_05(t(ii)+(dt/2),m1(ii)+(dt*e2/2),m2(ii)+(dt*e2/2),m3(ii)+(dt*e2/2),m4(ii)+(dt*e2/2),m5(ii)+(dt*e2/2));
    e4=f_ty_05(t(ii)+(dt),m1(ii)+dt*e3,m2(ii)+dt*e3,m3(ii)+dt*e3,m4(ii)+dt*e3,m5(ii)+dt*e3);
    m5(ii+1)=m5(ii)+(dt/6)*(e1+(2*(e2+e3))+e4);
end

plot(t,m1,t,m2,t,m3,t,m5,t,m5)
plot(t,m1,'LineWidth',2)
grid on
grid minor
hold on
plot(t,m2,'LineWidth',2)
grid on
grid minor
hold on
plot(t,m3,'LineWidth',2)
grid on
grid minor
hold on
plot(t,m4,'LineWidth',2)
grid on
grid minor
hold on
grid on
grid minor
plot(t,m5,'LineWidth',2)
hold off
title('Discontinous Milling Process')
xlabel('Time(t)','FontWeight','bold','FontSize',14)
ylabel('Sand Grain Size(m)','FontWeight','bold','FontSize',14)
legend('m1','m2','m3','m4','m5')
    
