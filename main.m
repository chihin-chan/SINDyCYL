clc
close all
clear 

nx = 513;
ny = 257;
dt = 0.25*20/nx;
snaps = 160;

% main program for SVD
[Y, x, y] = ASVEDE_extract(snaps,nx,ny);

% SVD
[U,S,V] = svd(Y, 'econ');
pcolor(reshape(U(:,1),[nx, ny])'); shading interp; 

weight1 = U(:,2)'*Y;
weight2 = U(:,3)'*Y;
weight3 = U(:,4)'*Y;

% Plotting Mode Shapes and the dynamics of POD Weights
tspan = dt:dt:(dt*snaps);

subplot(4,2,1)
pcolor(reshape(Y(:,i),[nx ny])'); shading interp;
title('Full Flow')
caxis([-0.05 0.05]); colormap('jet')
subplot(4,2,2)
pcolor(reshape(U(:,2),[nx ny])'); shading interp;
title('Mode 1')
caxis([-0.005 0.005]); colormap('jet')
subplot(4,2,3)
pcolor(reshape(U(:,3),[nx ny])'); shading interp;
title('Mode 2')
caxis([-0.005 0.005]); colormap('jet')
subplot(4,2,4)
pcolor(reshape(U(:,4),[nx ny])'); shading interp;
title('Mode 3')
caxis([-0.005 0.005]); colormap('jet')


% Pooling POD, weights together

x = [weight1' weight2' weight3'];
dx = [];
M = 3;

% Computing derivatives using finite difference
for i=3:length(x)-3
    for k=1:M
        dx(i-2,k) = (1/(12*dt))*(-x(i+2,k)+8*x(i+1,k)-8*x(i-1,k)+x(i-2,k));
    end
end

% Preparing for SINDy
x = x(3:end-3,:);
polyorder = 5;
usesine = 0;
Theta = poolData(x,M,polyorder,usesine);

% SINDy-ing
lambda = 0.00002;
Xi = sparsifyDynamics(Theta,dx,lambda,3);
% note that there are constant order terms... this is because fixed points
% are not at (0,0,0) for this data
poolDataLIST({'x','y','z'},Xi,M,polyorder,usesine);

% ODE to Solve
options = odeset('RelTol',1e-8,'AbsTol',1e-8*ones(1,M));
x0 = x(1,:);%
[tD,xD]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);

figure;
subplot(1,2,1)
color_line3(xD(:,1), xD(:,2), -xD(:,3), linspace(0,10,length(xD(:,1)))); view(27,16)
title("Reconstructed from SINDy");
grid on
subplot(1,2,2)
title("Dynamics from DNS");
color_line3(x(:,1), x(:,2), -x(:,3), linspace(0,10,length(x(:,1)))); view(27,16)
grid on


% Fancy Post Processing
v = VideoWriter('SINDyCYL')
v.FrameRate = 15;
open(v)
fig = figure();


for i=3:length(tspan)-6
    subplot(3,4,[1,4])
    pcolor(reshape(Y(:,i),[nx ny])'); shading interp;
    titletext = "Full Flow, Time-Step: " + tspan(i);
    title(titletext)
    caxis([-0.1 0.1]); colormap('jet')
    pbaspect([7 4 1]);
    set(gca, 'XTick',[], 'YTick', []);
    

    subplot(3,4,[5,6])
    plot3(xD(1:i-2,1), xD(1:i-2,2), -xD(1:i-2,3), '-o'); view(27,16)
    axis([min(xD(:,1)) max(xD(:,1)) min(xD(:,2)) max(xD(:,2)) -max(xD(:,3)) -min(xD(:,3))])
    title("Reconstructed from SINDy");
    xlabel('u_1');ylabel('u_2');zlabel('u_3');
    grid on
    subplot(3,4,[7,8])
    plot(x(1:i,1), x(1:i,2));
    xlabel('u_1');ylabel('u_2');
    title("Top View");
    axis([min(xD(:,1)) max(xD(:,1)) min(xD(:,2)) max(xD(:,2))]);
    grid on

    ytop = 129 +50;
    ybot = 129-50;
    xfront = 100;
    xback = 250;

    subplot(3,4,9)
    pcolor(reshape(U(:,1),[nx ny])'); shading interp;
    title('Mode 0 (u_0)')
    caxis([-0.005 0.005]); colormap('jet')
    xlim([xfront xback]), ylim([ybot ytop]);
    pbaspect([1 1 1]);
    set(gca, 'XTick',[], 'YTick', []);
    
    subplot(3,4,10)
    pcolor(reshape(U(:,2),[nx ny])'); shading interp;
    title('Mode 1 (u_1)')
    caxis([-0.005 0.005]); colormap('jet')
    xlim([xfront xback]), ylim([ybot ytop]);
    pbaspect([1 1 1]);
    set(gca, 'XTick',[], 'YTick', []);
    
    subplot(3,4,11)
    pcolor(reshape(U(:,3),[nx ny])'); shading interp;
    title('Mode 2 (u_2)')
    caxis([-0.005 0.005]); colormap('jet')
    xlim([xfront xback]), ylim([ybot ytop]);
    pbaspect([1 1 1]);
    set(gca, 'XTick',[], 'YTick', []);
    
    subplot(3,4,12)
    pcolor(reshape(U(:,4),[nx ny])'); shading interp;
    title('Mode 3 (u_3)')
    caxis([-0.005 0.005]); colormap('jet')
    xlim([xfront xback]), ylim([ybot ytop]);
    pbaspect([1 1 1]);
    set(gca, 'XTick',[], 'YTick', []);


    drawnow;

    x0 = 600;
    y0 = 600;
    width= 800;
    height = 900;
    set(gcf,'position',[x0,y0,width,height])

    frame = getframe(fig);
    writeVideo(v,frame);
end
close(v)

