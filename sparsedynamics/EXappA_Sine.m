% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data:
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

clear all, close all, clc
figpath = '../figures/';
addpath('./utils');

% generate Data
polyorder = 5;
usesine = 0;

x0list = [-1.25:.25:-.25 .25:.25:1.25]; % need more data to ID
xall = [];
dxall = [];
for k=1:length(x0list)
    x0 = x0list(k)
    
    % Integrate
    dt = 0.001;
    tspan=[dt:dt:5];
    N = length(tspan);
%     [t,x] = ode45(@(t,x)-x + (1/6)*x.^3 - (1/120)*x.^5,tspan,x0); 
    [t,x] = ode45(@(t,x)-sin(x),tspan,x0);
    % fourth order central difference
    for i=3:length(x)-3
        dx(i-2) = (1/(12*dt))*(-x(i+2)+8*x(i+1)-8*x(i-1)+x(i-2));
    end
    
    % concatenate
    xall = [xall; x(3:end-3)];
    dxall = [dxall; dx'];
end
dx = dxall;
x = xall;


%% pool Data  (i.e., build library of nonlinear time series)
clear Xi

% Case 1
polyorder = 5;
usesine = 0;
lambda = 0.1;

% % Case 2
% polyorder = 0;
% usesine = 1;
% lambda = 10;

% % Case 3
% polyorder = 5;
% usesine = 1;
% lambda = 10;

Theta = poolData(x,1,polyorder,usesine);
% normalize columns of Theta
for k=1:size(Theta,2)
    normTheta(k) = norm(Theta(:,k));
    Theta(:,k) = Theta(:,k)/normTheta(k);
end
m = size(Theta,2);
Xi = sparsifyDynamics(Theta,dx,lambda,1)  

% reverse of normalization from L66-69
Theta = poolData(x,1,polyorder,usesine);
for k=1:length(Xi)
    Xi(k,:) = Xi(k,:)/normTheta(k);
end

%% Simulate prediction, compare with true model
x0list = [-1.25:.25:-.25 .25:.25:1.25];%.1 -2 pi/2 1 -1 1.1];
x = [];
dx = [];
xall2 = [];
dxall2 = [];
for k=1:length(x0list)
    x0 = x0list(k)
    % x0=[2];  % Initial condition
    
    % Integrate
    dt = 0.001;
    tspan=[dt:dt:5];
    N = length(tspan);
    [t,x]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0);  % approximate    
    
    
    for i=3:length(x)-3
        dx(i-2) = (1/(12*dt))*(-x(i+2)+8*x(i+1)-8*x(i-1)+x(i-2));
    end
    
    % concatenate
    xall2 = [xall2; x(3:end-3)];
    dxall2 = [dxall2; dx'];
end

%% plotting
figure
subplot(2,1,1)
plot(dt*(1:length(xall)),xall,'k','LineWidth',1.5)
hold on, grid on
plot(dt*(1:length(xall2)),xall2,'r--','LineWidth',1.5)
subplot(2,1,2)
plot(dt*(1:length(xall)),dxall,'k','LineWidth',1.5)
hold on, grid on
plot(dt*(1:length(xall2)),dxall2,'r--','LineWidth',1.5)
xlabel('Time')