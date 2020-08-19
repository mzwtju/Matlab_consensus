close all
% clear
% clc
% 
% %% Scale
% N = 8;
% dimx = 3;
% 
% %% Derived Data
% Bt = [0 1 0;
%       0 0 1];
% Bb = [eye(2), zeros(2,1)];
% alp = 1;
% tau0 = 10;
% 
% %% Chosen Scheme
% % sgm = @(t) 1;
% sgm = @(t) 0 ...
%            + 1*(t >=      0 && t < 1*tau0 || t >= 4*tau0 && t < 5*tau0) ...
%            + 2*(t >= 1*tau0 && t < 3*tau0 || t >= 9*tau0) ...
%            + 3*(t >= 3*tau0 && t < 4*tau0 || t >= 8*tau0 && t < 9*tau0) ...
%            + 4*(t >= 5*tau0 && t < 8*tau0);
% % 
% 
% %% Randomly Generated Data
% x0 = rand(dimx,N) - 0.5; %-0.5 0.5
% for nodeidx = 1:N
%     x0(:,nodeidx) = nodeidx * x0(:,nodeidx);
% end
% 
% %% Simulation
% fps = 12;
% tF = 30;
% tUnitStep = 1/fps;
% tSteps = tF/tUnitStep; %360
% taxis = linspace(0,tF,tSteps+1);
% 
% x_valts = [reshape(x0, [], 1), NaN(N*dimx, tSteps)]; %生成24*361矩阵
% % K1 = zeros(2,3);
% % % To assign eigenvalues of A+BK1+BK2
% % K2 = [3.9988  -4.9905  -3.0122
% %      -7.0005  -4.9783  -5.9995];
% % A = [4 -2 2
% %      1  3 5
% %      2  7 4];
%  B = [0 0;eye(2)];
% % C = [1 0 0 ;0 1 0 ;0 0 1 ];
% % D = zeros(3,2)
% % [z1,p1,~] = ss2zp(A+B*K1+B*K2,B,C,D,1)
% % [z2,p2,~] = ss2zp(A+B*K1+B*K2,B,C,D,2)
% % % [z3,p3,~] = ss2zp(A+B*K1+B*K2,B,C,D,3)
% % ii = 1;
% % t = 1;
% % v = @(t) [v(t);
% %         - 6*cos(2*t + (pi*ii)/4) - 42*sin(2*t + (pi*ii)/4) - 5*sin(16*t + (pi*ii)/4);
% %  16*cos(16*t + (pi*ii)/4) - 24*cos(2*t + (pi*ii)/4) - 108*sin(2*t + (pi*ii)/4) - 4*sin(16*t + (pi*ii)/4)];
% % % sqr = @(x)[ x.^2;x*3];
% % % sqr(2)
% K2 = [3.9988  -4.9905  -3.0122
%      -7.0005  -4.9783  -5.9995];
% K3 = [-13.9520  8.5232  -2.45
%         7.0585 -2.4500   5.0634];
% N = 8;
% for ii = 0:N-1
%     h = @(t) [h(t); 6*sin(ii*t); 12*sin(ii*t); 6*cos(2*t)+sin(16*t)];
% end
%  x_dot = @(t)(kron(eye(8),eye(3)) * h(t))
%  N = 3
% syms t h(t)
% h = @(t) [];
% ii = 0;
% % h = [6*sin(ii*t); 12*sin(ii*t); 6*cos(2*t)+sin(16*t)];
% for ii = 0:N
%     h = [h(t); 6*sin(ii*t); 12*sin(ii*t); 6*cos(2*t)+sin(16*t)];
% end
% h
% xdot = @(t) ...
%     (kron(alp*G.L, B*K3) - kron(eye(N), B*K2)) * h(t) 
% test_fun = kron(eye(N),eye(3))*h(t)
% tidx = 2;
% alp = 1;
% G.L = eye(8);
%     tjust = taxis(tidx-1);
%     tnow  = taxis(tidx);
%    x_dot = @(t,x)(kron(alp*G.L, B*K3) - kron(eye(N), B*K2)) * h(t)
% [~, x_odetrj] = ode45( x_dot, [tjust, tnow], x_valts(:,tidx-1) )
N = 8
h = @(t) []; % 初始化为空
for ii = 0:N-1
       h = @(t) [h(t) ];
end
tmp = h(1) % 看一下 h 的结果
% x_dot = @(t)(kron(eye(8),eye(3)) * h(t))