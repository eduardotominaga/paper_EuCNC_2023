%% Reset
clearvars
close all
clc

%% Parameters
N=1e3;          % Number of network realizations
K_total=1e3;    % Total number of MTDs on each realization
l=250;          % Length of the side of the square area [m]
rng(1);         % Set the random number generator
sigma=25;       % Intensity of the alarm event
nBins=25;       % Number of bins for the histograms
range=0:1:l;    % Range of distances

%% Determining the position of the alarm event:
da_x=l/4;               % x-position of the alarm event [m]
da_y=l/2;               % x-position of the alarm event [m]

%% Obtaining the empirical PDF:

% ATPF:
gaussianATPF = @(d,sigma) exp(-(d^2)/(2*sigma^2));

K_act=0;
dx_act=[];
dy_act=[];
for n=1:N
    
    % Determining the positions of all MTDs:
    dx_MTD=l*rand(1,K_total);       % x-position of the MTDs [m]
    dy_MTD=l*rand(1,K_total);       % y-position of the MTDs [m]
    
    % Determining which MTDs are active:
    for k=1:K_total
        d=sqrt((dx_MTD(k)-da_x)^2+(dy_MTD(k)-da_y)^2);
        if rand<gaussianATPF(d,sigma)
            K_act=K_act+1;
            dx_act=[dx_act dx_MTD(k)];
            dy_act=[dy_act dy_MTD(k)];
        end
    end
end

%% Generating the theoretical truncated Gaussian PDFs:

% Theoretical PDF for the x-axis:
PDx=makedist('Normal',da_x,sigma);
PDx=truncate(PDx,0,l);

% Theoretical PDF for the x-axis:
PDy=makedist('Normal',da_y,sigma);
PDy=truncate(PDy,0,l);

%% Colors for the plot:
blue = [57 106 177]./255;
red = [204 37 41]./255;
black = [83 81 84]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;
yellow = [240 180 50]./255;
colorMap=[black;red;green;blue;yellow;purple;brown];

%% Plotting the PDFs

% Parameters for the plot:
fontSize=18;
lineWidth=2.5;

fig1=figure(1);
    fig1.Position=[100 300 550 400];
         
    % Plot the empirical PDF using histogram:
    histogram(dx_act,nBins,'FaceColor',blue,'Normalization','pdf'); 
    hold on
        
    % Plot the empirical PDF using ksdensity:
%     [f,x1]=ksdensity(dx_act);
%     p1=plot(x1,f,'LineWidth',lineWidth);    
%     hold on
    
    % Plot the theoretical PDF:
    plot(range,pdf(PDx,range),'Color',red,'LineWidth',lineWidth)
    xline(da_x,'--','LineWidth',lineWidth)

    ylim([0 0.025])
    grid on
    leg=legend('Empirical PDF','Theoretical PDF');
    set(leg,'Interpreter','latex','FontSize',fontSize,'Location','Northeast');
    xlabel('$d_{k,x}$','Interpreter','latex','FontSize',fontSize)
    ylabel('PDF','Interpreter','latex','FontSize',fontSize)
    set(gca,'TickLabelInterpreter','latex','FontSize',fontSize)
    saveas(fig1,'PDF_x.png')
    saveas(fig1,'PDF_x.eps','epsc')

    
fig2=figure(2);
    fig2.Position=[700 300 550 400];

    % Plot the empirical PDF using histogram:
    histogram(dy_act,nBins,'FaceColor',blue,'Normalization','pdf'); 
    hold on

    % Plot the theoretical PDF:
    plot(range,pdf(PDy,range),'Color',red,'LineWidth',lineWidth)
    xline(da_y,'--','LineWidth',lineWidth)
    grid on
    
    ylim([0 0.025])
    leg=legend('Empirical PDF','Theoretical PDF');
    set(leg,'Interpreter','latex','FontSize',fontSize,'Location','Northeast');
    xlabel('$d_{k,y}$','Interpreter','latex','FontSize',fontSize)
    set(gca,'TickLabelInterpreter','latex','FontSize',fontSize)
    ylabel('PDF','Interpreter','latex','FontSize',fontSize)
    saveas(fig2,'PDF_y.png')
    saveas(fig2,'PDF_y.eps','epsc')

fig3=figure(3);
    fig3.Position=[1300 300 550 500];
    dscatter(dx_act',dy_act');
    colormap('jet')
%     hold on
%     scatter(da_x,da_y,100,'x','LineWidth',lineWidth)
    box on
    grid on
    xlabel('$x$ [m]','Interpreter','latex','FontSize',fontSize)
    ylabel('$y$ [m]','Interpreter','latex','FontSize',fontSize)
    xlim([0 l])
    ylim([0 l])
    set(gca,'TickLabelInterpreter','latex','FontSize',fontSize)
    saveas(fig3,'heatmapActiveMTDs.png')
    saveas(fig3,'heatmapActiveMTDs.eps','epsc')

    
%% Plot the 3D Gaussian distribution
% mu=[da_x da_y];             % Mean vector
% Sigma=[sigma^2 0;0 sigma^2];    % Covariance matrix
% 
% % Create a grid of evenly spaced points in 2D:
% x1=0:5:l;
% x2=0:5:l;
% [X1,X2]=meshgrid(x1,x2);
% X=[X1(:) X2(:)];
% 
% % Compute the 2d Gaussian PDF:
% y=mvnpdf(X,mu,Sigma);
% y=reshape(y,length(x2),length(x1));
% y=truncate(y,[0 0],[l l]);
% 
% % Plot the PDF values:
% fig4=figure(4);
%     surf(x1,x2,y,'FaceAlpha',0.5)
%     hold
%     histogram2(dx_act',dy_act','Normalization','pdf')
    
    
