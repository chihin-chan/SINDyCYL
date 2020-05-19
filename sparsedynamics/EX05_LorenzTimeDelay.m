% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

clear all, close all, clc
figpath = '../figures/';
addpath('./utils');

%% generate Data
polyorder = 5;
usesine = 0;

sigma = 10;  % Lorenz's parameters (chaotic)
beta = 8/3;
rho = 28;

n = 3;
x0=[-8; 8; 27];  % Initial condition

% Integrate
dt = 0.001;
tspan=[dt:dt:100];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,options);
%% make a movie of the attractor when i increase shift-stack number
clear V, clear dV
stackmax = 10;  % changing the number of shift-stacked rows
r = 3;  % number of modes to keep
X = zeros(stackmax,size(x,1)-stackmax);
for k=1:r
    X(k,:) = x(k:end-stackmax-1+k,1);
end
[U,S,V] = svd(X,'econ');

% compute derivative using fourth order central difference
% use TVRegDiff if more error 
for i=3:length(V)-3
    for k=1:r
        dV(i-2,k) = (1/(12*dt))*(-V(i+2,k)+8*V(i+1,k)-8*V(i-1,k)+V(i-2,k));
    end
end  

%% concatenate
x = V(3:end-3,1:r);
dx = dV;

%% pool Data  (i.e., build library of nonlinear time series)
clear Xi
polyorder = 3;
Theta = poolData(x,r,polyorder,usesine);
% normalize columns of Theta (required in new time-delay coords)
for k=1:size(Theta,2)
    normTheta(k) = norm(Theta(:,k));
    Theta(:,k) = Theta(:,k)/normTheta(k);
end 
m = size(Theta,2);
% compute Sparse regression: sequential least squares
% requires different lambda parameters for each column
Xi(:,1) = sparsifyDynamics(Theta,dx(:,1),.01,1);
Xi(:,2) = sparsifyDynamics(Theta,dx(:,2),.2,1);
Xi(:,3) = sparsifyDynamics(Theta,dx(:,3),2,1); 
Xi
Theta = poolData(x,r,polyorder,usesine);
for k=1:length(Xi)
    Xi(k,:) = Xi(k,:)/normTheta(k);
end
%%
figure
subplot(3,1,1)
plot(t(1:length(Theta)),Theta*Xi(:,1),'k')
xlim([40 60])
hold on
plot(t(1:length(Theta)),dx(:,1),'r--')
ylabel('u')
subplot(3,1,2)
plot(t(1:length(Theta)),Theta*Xi(:,2),'k')
xlim([40 60])
hold on
plot(t(1:length(Theta)),dx(:,2),'r--')
ylabel('v')
subplot(3,1,3)
plot(t(1:length(Theta)),Theta*Xi(:,3),'k')
xlim([40 60])
hold on
plot(t(1:length(Theta)),dx(:,3),'r--')
ylabel('w')
xlabel('Time')
l1=legend('$\Theta(V)\Xi$','$\dot{V}$')
set(l1,'interpreter','latex')
%%
tspan = [0 10.];
x0 = x(50,:);
[tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);  % approximate
figure
plot3(xB(:,1),xB(:,2),xB(:,3),'k');
hold on
axis([-.01 .01 -.02 .01 -.015 .015])
xlabel('u'), ylabel('v'), zlabel('w')
set(gcf,'Position',[100 100 400 300])

figure
subplot(1,3,1)
plot(xB(:,1),xB(:,2),'k');
axis([-.01 .01 -.01 .01])
hold on
xlabel('u'), ylabel('v')
subplot(1,3,2)
plot(xB(:,1),xB(:,3),'k');
hold on
axis([-.01 .01 -.01 .015])
xlabel('u'), ylabel('w')
subplot(1,3,3)
plot(xB(:,2),xB(:,3),'k'); 
hold on
axis([-.01 .01 -.01 .015])
xlabel('v'), ylabel('w')
set(gcf,'Position',[100 100 700 240])