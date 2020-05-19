% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

clear all, close all, clc

figpath = '../figures/';
addpath('./utils');
addpath('./data');

dt = 0.02;
r=2;

polyorder = 5;
usesine = 0;

%% load data from first run and compute derivative
load PODcoefficients
x = [alpha(1:5000,1:r) alphaS(1:5000,1)];
M = size(x,2);
% compute Derivative
eps = 0;
for i=3:length(x)-3
    for k=1:M
        dx(i-2,k) = (1/(12*dt))*(-x(i+2,k)+8*x(i+1,k)-8*x(i-1,k)+x(i-2,k));
    end
end  

%% load data from second run and compute derivative
load PODcoefficients_run1
x1 = [alpha(1:3000,1:r) alphaS(1:3000,1)];

M = size(x1,2);
% compute Derivative
eps = 0;
for i=3:length(x1)-3
    for k=1:M
        dx1(i-2,k) = (1/(12*dt))*(-x1(i+2,k)+8*x1(i+1,k)-8*x1(i-1,k)+x1(i-2,k));
    end
end  

%% concatenate
x = [x(3:end-3,:); x1(3:end-3,:)];
dx = [dx; dx1];

%% pool Data
Theta = poolData(x,M,polyorder,usesine);
m = size(Theta,2);
dx = dx(1:end,:);

lambda = 0.000001;
Xi = sparsifyDynamics(Theta,dx,lambda,3)
% note that there are constant order terms... this is because fixed points
% are not at (0,0,0) for this data
poolDataLIST({'x','y','z'},Xi,M,polyorder,usesine);


%% - second figure: initial portion
figure(1)

tspan = [0:.02:100];
options = odeset('RelTol',1e-8,'AbsTol',1e-8*ones(1,r+1));
x0 = x(1,:);%
[tD,xD]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);
%
subplot(1,2,1)
color_line3(x(1:5000-6,1),x(1:5000-6,2),x(1:5000-6,end),.02*(1:length(x(1:5000-6,end))),'LineWidth',1.5);
view(27,16)
grid on
xlabel('x')
ylabel('y')
zlabel('z')
axis([-200 200 -200 200 -160 20])
subplot(1,2,2)
dtD = [0; diff(tD)];
color_line3(xD(:,1),xD(:,2),xD(:,end),.02*(1:length(xD(:,1))),'LineWidth',1.5);
view(27,16)
grid on
xlabel('x')
ylabel('y')
zlabel('z')
axis([-200 200 -200 200 -160 20])

%% - second figure: Lorenz for t=20
figure(2)

tspan = [0:.02:95];
options = odeset('RelTol',1e-8,'AbsTol',1e-8*ones(1,r+1));
x0 = x(5001,:);%
[tD,xD]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);
%
subplot(1,2,1)
color_line3(x(5001:end,1),x(5001:end,2),x(5001:end,end),.02*(1:length(x(5001:end,end))),'LineWidth',1.5);
view(27,16)
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)
axis([-200 200 -200 200 -160 20])
subplot(1,2,2)
dtD = [0; diff(tD)];
color_line3(xD(:,1),xD(:,2),xD(:,end),.02*(1:length(xD(:,end))),'LineWidth',1.5);
view(27,16)
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)
set(gca,'FontSize',13)
axis([-200 200 -200 200 -160 20])