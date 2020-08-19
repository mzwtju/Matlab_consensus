clc
clear all
close all
%%model
% symbol k11 k12 k21 k22
N = 500;
B1= [1;0];
B2 = [0;1];
K1= [-1,-0.8];
% K2 = [k21,k22];
dt = 0.1;
[u1,u2,u3,u4,u5]=deal([0,0]);
[v1,v2,v3,v4,v5]=deal([0,0]);
x_dot = [0;0];
x1_0=[4.75,0.19];
x2_0 = [-2.34,4.38];
x3_0 = [-2.62,-4.7];
x1=[4.75,0.19];
x2 = [-2.34,4.38];
x3 = [-2.62,-4.7];
v1 = [0,-0.04];
v2=[-0.06,-0.01];
v3 = [-0.01,0.03];
x_test = [0,0];
x_test2 = [0,0];
delta_matrix = eye(5);
hx = [0,0];
hv = [0,0];
A_matrix = [0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1;1,0,0,0,0]; 
L_matrix = delta_matrix - A_matrix;
% A_matrix = [0,1,0;0,0,1;1,0,0];
% L_matrix = delta_matrix - A_matrix;
theta1(1:2,1:2) = [x1;v1];
theta2(1:2,1:2) = [x2;v2];
theta3(1:2,1:2) = [x3;v3];
[hiX,hiY,hivX,hivX_dot,hivY,hivY_dot,g] = deal(zeros(3,1257));
h = [hx;hv];
w = 0.2;
r = 5;
t = 0:0.1:40*pi;
tspan = [0,40*pi];
T=length(t)-1;
for i = 1:3
        g(i,:) = sign(sin(w*t/2+pi*(i-1)/3));
        hiX(i,:) = r*cos((w*t+2*pi*(i-1)/3));
        hiY(i,:) = r*sin((w*t+2*pi*(i-1)/3));
        hivX(i,:) = -w*r*sin((w*t+2*pi*(i-1)/3));
        hivY(i,:) = w*r*cos((w*t+2*pi*(i-1)/3));
        hivX_dot(i,:) = -w*w*r*cos((w*t+2*pi*(i-1)/3));
        hivY_dot(i,:) = -w*w*r*sin((w*t+2*pi*(i-1)/3));
%         hiX(i,:) = r*cos((w*t+2*pi*(i-1)/3)).*g(i,:);
%         hiY(i,:) = r*sin((w*t+2*pi*(i-1)/3)).*g(i,:);
%         hivX(i,:) = -w*r*sin((w*t+2*pi*(i-1)/3)).*g(i,:);
%         hivY(i,:) = w*r*cos((w*t+2*pi*(i-1)/3)).*g(i,:);
%         hivX_dot(i,:) = -w*w*r*cos((w*t+2*pi*(i-1)/3)).*g(i,:);
%         hivY_dot(i,:) = -w*w*r*sin((w*t+2*pi*(i-1)/3)).*g(i,:);
end
i = 0;

h1 = [hiX(1,:);hiY(1,:);hivX(1,:);hivY(1,:)];
h2 = [hiX(2,:);hiY(2,:);hivX(2,:);hivY(2,:)];
h3 = [hiX(3,:);hiY(3,:);hivX(3,:);hivY(3,:)];
% h4 = [hiX(4,:);hiY(4,:);hivX(4,:);hivY(4,:)];
% h5 = [hiX(5,:);hiY(5,:);hivX(5,:);hivY(5,:)];
h1v_dot = [hivX_dot(1,:)',hivY_dot(1,:)'];
h2v_dot = [hivX_dot(2,:)',hivY_dot(2,:)'];
h3v_dot = [hivX_dot(3,:)',hivY_dot(3,:)'];
% h4v_dot = [hivX_dot(4,:)',hivY_dot(4,:)'];
% h5v_dot = [hivX_dot(5,:)',hivY_dot(5,:)'];
%% Riccati
% A = (B2*K1+B1*B2');
% B = B2;
% Q = eye(2);
% R = [1];
% [P_conj,gain,vectors] = idare(A,B,Q,R,[],[]);
% P = conj(P_conj);
% lambda = eig(L_matrix);
% lambda_real = sort(real(lambda)); 
% K2 = inv(lambda_real(2))*B2'*P_conj;
K2= [0.2761, 0.5141];
%%

for i = 1:T

  h1_ = [h1(1,i),h1(2,i);h1(3,i),h1(4,i)];
  h2_ = [h2(1,i),h2(2,i);h2(3,i),h2(4,i)];
  h3_ = [h3(1,i),h3(2,i);h3(3,i),h3(4,i)];
  u1(i,:) =K1*(theta1-h1_)+K2*((theta3-h3_)-(theta1-h1_))+h1v_dot(i,:);
%        ...
%        ;
  u2(i,:) =K1*(theta2-h2_)+K2*((theta1-h1_)-(theta2-h2_))+h2v_dot(i,:);
%        +...
%        ;
  u3(i,:) =K1*(theta3-h3_)+K2*((theta2-h2_)-(theta3-h3_))+h3v_dot(i,:);
%        ...
%        ;
    if i == 1
      v1(i,:) =  u1(i,:).*dt;
      v2(i,:) =  u2(i,:).*dt;
      v3(i,:) =  u3(i,:).*dt ;
   else
      v1(i,:) = v1(i-1,:) + u1(i,:).*dt;
      v2(i,:) = v2(i-1,:) + u2(i,:).*dt;
      v3(i,:) = v3(i-1,:) + u3(i,:).*dt ;   
    end
   
   if i == 1
      x1(i,:) = x1_0 + v1(i,:).*dt;
      x2(i,:) = x2_0 + v2(i,:).*dt;
      x3(i,:) = x3_0 + v3(i,:).*dt ;
       x_test(i,1) =hivX(1,i)*dt;
       x_test(i,2) =hivY(1,i)*dt;
       x_test2(i,1) =hivX(2,i)*dt;
       x_test2(i,2) =hivY(2,i)*dt
   else

      x1(i,:) = x1(i-1,:) + v1(i,:).*dt;
      x2(i,:) = x2(i-1,:) + v2(i,:).*dt;
      x3(i,:) = x3(i-1,:) + v3(i,:).*dt ;   
       x_test(i,1) =x_test(i-1,1)+hivX(1,i)*dt;
       x_test(i,2) =x_test(i-1,2)+hivY(1,i)*dt;
       x_test2(i,1) =x_test2(i-1,1)+hivX(2,i)*dt;
       x_test2(i,2) =x_test2(i-1,2)+hivY(2,i)*dt;
   end
  theta1(1:2,1:2) = [x1(i,:);v1(i,:)];
  theta2(1:2,1:2) = [x2(i,:);v2(i,:)];
  theta3(1:2,1:2) = [x3(i,:);v3(i,:)];
end
[t_2,v1]=ode45(@(u1,v1) u1,tspan,0);
%% plot
figure(1)
plot(x1(:,1),x1(:,2),'r')
hold on 
plot(x2(:,1),x2(:,2),'b')
plot(x3(:,1),x3(:,2),'k')
plot(x1_0(1,1),x1_0(1,2),'ro');
plot(x2_0(1,1),x2_0(1,2),'bo');
plot(x3_0(1,1),x3_0(1,2),'ko');
% figure(2)
% % plot(1:T,v1(:,1))
% plot(x_test(:,1),x_test(:,2))

