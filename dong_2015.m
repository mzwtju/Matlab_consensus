clc
clear all
close all
%%model
% symbol k11 k12 k21 k22
N = 500;
B1= [1;0];
B2 = [0;1];
K1= [-2,-1.2];
% K2 = [k21,k22];
dt = 0.1;
[u1,u2,u3,u4,u5]=deal([0,0]);
[v1,v2,v3,v4,v5]=deal([0,0]);
x_dot = [0;0];
[x1,x2,x3,x4,x5]=deal( [0,0]);
x1_0 = [-0,-5.7]
x2_0 = [6.7,1.7];
x3_0 = [0.4,6.9];
x4_0 = [-6.5,2.5];
x5_0 = [-4.4,-5.4];
% x1_0 = [4.16,-1.3]
% x2_0 = [0.78,4.113];
% x3_0 = [-0.2,5.08];
% x4_0 = [-3.8,-1.56];
% x5_0 = [-0.7,-4.9];
v1_0 = [0.03,1.4];
v2_0=[-1.2,-0.4];
v3_0 = [-0.8,-1.2];
v4_0 = [0.8,-1.2];
v5_0 = [1.4,0.4];
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
theta4(1:2,1:2) = [x4;v4];
theta5(1:2,1:2) = [x5;v5];

h = [hx;hv];
w = 0.214;
r = 7;
t = 0:0.1:10*pi;
T=length(t)-1;
for i = 1:5         
        g(i,:) = sign(sin(w*t/2+pi*(i-1)/5));
        hiX(i,:) = r*cos((w*t+2*pi*(i-1)/5)-1).*g(i,:);
        hiY(i,:) = r*sin((w*t+2*pi*(i-1)/5)-1);
        hivX(i,:) = -w*r*sin((w*t+2*pi*(i-1)/5)).*g(i,:);
        hivY(i,:) = w*r*cos((w*t+2*pi*(i-1)/5));
        hivX_dot(i,:) = -w*w*r*cos((w*t+2*pi*(i-1)/5)).*g(i,:);
        hivY_dot(i,:) = -w*w*r*sin((w*t+2*pi*(i-1)/5));
%                 hiX(i,:) = r*cos((w*t+2*pi*(i-1)/5)-1);
%         hiY(i,:) = r*sin((w*t+2*pi*(i-1)/5)-1);
%         hivX(i,:) = -w*r*sin((w*t+2*pi*(i-1)/5));
%         hivY(i,:) = w*r*cos((w*t+2*pi*(i-1)/5));
%         hivX_dot(i,:) = -w*w*r*cos((w*t+2*pi*(i-1)/5));
%         hivY_dot(i,:) = -w*w*r*sin((w*t+2*pi*(i-1)/5));
end
i = 0;

h1 = [hiX(1,:);hiY(1,:);hivX(1,:);hivY(1,:)];
h2 = [hiX(2,:);hiY(2,:);hivX(2,:);hivY(2,:)];
h3 = [hiX(3,:);hiY(3,:);hivX(3,:);hivY(3,:)];
h4 = [hiX(4,:);hiY(4,:);hivX(4,:);hivY(4,:)];
h5 = [hiX(5,:);hiY(5,:);hivX(5,:);hivY(5,:)];
h1v_dot = [hivX_dot(1,:)',hivY_dot(1,:)'];
h2v_dot = [hivX_dot(2,:)',hivY_dot(2,:)'];
h3v_dot = [hivX_dot(3,:)',hivY_dot(3,:)'];
h4v_dot = [hivX_dot(4,:)',hivY_dot(4,:)'];
h5v_dot = [hivX_dot(5,:)',hivY_dot(5,:)'];
%% Riccati
A = (B2*K1+B1*B2');
B = B2;
Q = eye(2);
R = [1];
[P_conj,gain,vectors] = idare(A,B,Q,R,[],[]);
P = conj(P_conj);
lambda = eig(L_matrix);
lambda_real = sort(real(lambda)); 
K2 = inv(lambda_real(2))*B2'*P_conj;
K2 = [0.3416,0.7330] ;
%%
theta1(1:2,1:2) = [x1_0;v1_0];
theta2(1:2,1:2) = [x2_0;v2_0];
theta3(1:2,1:2) = [x3_0;v3_0];
theta4(1:2,1:2) = [x4_0;v4_0];
theta5(1:2,1:2) = [x5_0;v5_0];

for i = 1:T
  h1_ = [h1(1,i),h1(2,i);h1(3,i),h1(4,i)];
  h2_ = [h2(1,i),h2(2,i);h2(3,i),h2(4,i)];
  h3_ = [h3(1,i),h3(2,i);h3(3,i),h3(4,i)];
  h4_ = [h4(1,i),h4(2,i);h4(3,i),h4(4,i)];
  h5_ = [h5(1,i),h5(2,i);h5(3,i),h5(4,i)];
  u1(i,:) =K1*(theta1-h1_)+K2*((theta5-h5_)-(theta1-h1_))+h1v_dot(i,:);
       
  u2(i,:) =K1*(theta2-h2_)+K2*((theta1-h1_)-(theta2-h2_))+h2v_dot(i,:);
       
  u3(i,:) =K1*(theta3-h3_)+K2*((theta2-h2_)-(theta3-h3_))+h3v_dot(i,:);
%   ;
  u4(i,:) =K1*(theta4-h4_)+K2*((theta3-h3_)-(theta4-h4_))+h4v_dot(i,:);
%   ;
  u5(i,:) =K1*(theta5-h5_)+K2*((theta4-h4_)-(theta5-h5_))+h5v_dot(i,:);
%   ;

    if i == 1
      v1(i,:) =  u1(i,:).*dt;
      v2(i,:) =  u2(i,:).*dt;
      v3(i,:) =  u3(i,:).*dt ;
      v4(i,:) =  u4(i,:).*dt;
      v5(i,:) =  u5(i,:).*dt ;
   else
      v1(i,:) = v1(i-1,:) + u1(i,:).*dt;
      v2(i,:) = v2(i-1,:) + u2(i,:).*dt;
      v3(i,:) = v3(i-1,:) + u3(i,:).*dt ;   
      v4(i,:) = v4(i-1,:) + u4(i,:).*dt;
      v5(i,:) = v5(i-1,:) + u5(i,:).*dt ;   
    end
   
   if i == 1
      x1(i,:) = x1_0 +v1(i,:).*dt;
      x2(i,:) = x2_0 +v2(i,:).*dt;
      x3(i,:) = x3_0 +v3(i,:).*dt ;
      x4(i,:) = x4_0 +v4(i,:).*dt;
      x5(i,:) = x5_0 +v5(i,:).*dt ;
       x_test(i,1) =hivX(1,i)*dt;
       x_test(i,2) =hivY(1,i)*dt;
       x_test2(i,1) =hivX(2,i)*dt;
       x_test2(i,2) =hivY(2,i)*dt
   else

      x1(i,:) = x1(i-1,:) + v1(i,:).*dt;
      x2(i,:) = x2(i-1,:) + v2(i,:).*dt;
      x3(i,:) = x3(i-1,:) + v3(i,:).*dt ;   
      x4(i,:) = x4(i-1,:) + v4(i,:).*dt;
      x5(i,:) = x5(i-1,:) + v5(i,:).*dt ;   
       x_test(i,1) =x_test(i-1,1)+hivX(1,i)*dt;
       x_test(i,2) =x_test(i-1,2)+hivY(1,i)*dt;
       x_test2(i,1) =x_test2(i-1,1)+hivX(2,i)*dt;
       x_test2(i,2) =x_test2(i-1,2)+hivY(2,i)*dt;
   end
theta1(1:2,1:2) = [x1(i,:);v1(i,:)];
theta2(1:2,1:2) = [x2(i,:);v2(i,:)];
theta3(1:2,1:2) = [x3(i,:);v3(i,:)];
theta4(1:2,1:2) = [x4(i,:);v4(i,:)];
theta5(1:2,1:2) = [x5(i,:);v5(i,:)];
end
figure(1)
plot(x1(:,1),x1(:,2),'r')
hold on 
plot(x2(:,1),x2(:,2),'b')
plot(x3(:,1),x3(:,2),'k')
plot(x4(:,1),x4(:,2),'g')
plot(x5(:,1),x5(:,2),'c')
plot(x1_0(1,1),x1_0(1,2),'ro');
plot(x2_0(1,1),x2_0(1,2),'bo');
plot(x3_0(1,1),x3_0(1,2),'ko');
plot(x4_0(1,1),x4_0(1,2),'go');
plot(x5_0(1,1),x5_0(1,2),'co');
% plot(hiX(1,:),hiY(1,:));
figure(2)
plot(x_test2(:,1),x_test2(:,2))


