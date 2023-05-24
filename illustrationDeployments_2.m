%% Reset
clearvars
close all
clc

%% Parameters
d=250;          % Length of the side of the square area [m]
fontSize=18;    % Font size for the plots
markerSize=100; % Marker size for the plots

%% Define my preferred color map for plots:
blue = [57 106 177]./255;
red = [204 37 41]./255;
black = [83 81 84]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;
yellow = [240 180 50]./255;
myColorMap=[black;red;green;blue;yellow;purple;brown];

%% Illustration: centralized mMIMO
% Determining the position of the BS:
dx_BS=d/2;      % x-position of the BS [m]
dy_BS=d/2;      % x-position of the BS [m]

% Plot the spatial distribution of BS and MTD:
fig1=figure(1);
    fig1.Position=[500 200 1000 475];
    hold on
    p1=scatter(dx_BS,dy_BS,2*markerSize,'MarkerEdgeColor',blue,'MarkerFaceColor',blue,...
        'DisplayName','BS in the CmMIMO deployment');
    grid on
    xlim([0 d])
    ylim([0 d])
    box on
    xlabel('$d_X$ [m]','Interpreter','latex','FontSize',fontSize)
    ylabel('$d_Y$ [m]','Interpreter','latex','FontSize',fontSize)
    set(gca,'TickLabelInterpreter','latex','FontSize',fontSize)
%     leg=legend('BS');
%     set(leg,'Interpreter','latex','FontSize',fontSize,'Location','Southeast')
%     saveas(fig1,'CmMIMO.png')
%     saveas(fig1,'CmMIMO.eps','epsc')

%% Illustration: Distributed mMIMO, grid deployment

Q=16;           % Number of APs

% Determining the position of the APs:
range_AP=d/(2*sqrt(Q)):d/sqrt(Q):d-d/(2*sqrt(Q));     
dx_AP=zeros(1,Q);                           % x-position of the APs [m]
dy_AP=zeros(1,Q);                           % y-position of the APs [m]
for idx=1:sqrt(Q)
   dx_AP(1+sqrt(Q)*(idx-1):sqrt(Q)*idx)=range_AP(idx); 
end
for idx=1:sqrt(Q)
   dy_AP(1+sqrt(Q)*(idx-1):sqrt(Q)*idx)=range_AP; 
end

% Plot the spatial distribution of APs and MTD:
% fig2=figure(2);
%     fig2.Position=[-750 500 480 500];
    p2=scatter(dx_AP,dy_AP,markerSize,'MarkerEdgeColor',yellow,'MarkerFaceColor',yellow,...
        'DisplayName','AP in the DmMIMO/grid deployment');
    grid on
    xlim([0 d])
    ylim([0 d])
    box on
    xlabel('$d_X$ [m]','Interpreter','latex','FontSize',fontSize)
    ylabel('$d_Y$ [m]','Interpreter','latex','FontSize',fontSize)
    set(gca,'TickLabelInterpreter','latex','FontSize',fontSize)
%     leg=legend('AP');
%     set(leg,'Interpreter','latex','FontSize',fontSize,'Location','Southeast')
%     saveas(fig2,'DmMIMO_grid.png')
%     saveas(fig2,'DmMIMO_grid.eps','epsc')
    
%% Illustration: Distributed mMIMO, Radio Stripes

Q=64;   % Number of APs in the radio stripe
% Define the positions of the APs:
range_AP=d/(2*(Q/4)):d/(Q/4):d-d/(2*(Q/4));   
dx_AP=zeros(1,Q);                   % x-position of the APs [m]
dy_AP=zeros(1,Q);                   % y-position of the APs [m]
for idx=1:4         % A square area has four sides.
    switch idx
        case 1
            dy_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=range_AP;
        case 2
            dx_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=range_AP;
            dy_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=d*ones(1,Q/4);            
        case 3
            dx_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=d*ones(1,Q/4);
            dy_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=flip(range_AP);
        case 4
            dx_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=flip(range_AP);
    end
end

% Plot the geographic distribution of APs and MTD:
% fig3=figure(3);
%     fig3.Position=[-750 500 480 500];
    p3=scatter(dx_AP,dy_AP,markerSize,'MarkerEdgeColor',green,'MarkerFaceColor',green,...
        'DisplayName','AP in the DmMIMO/linear deployment');
    grid on
    xlim([0 d])
    ylim([0 d])
    box on
    xlabel('$d_X$ [m]','Interpreter','latex','FontSize',fontSize)
    ylabel('$d_Y$ [m]','Interpreter','latex','FontSize',fontSize)
    set(gca,'TickLabelInterpreter','latex','FontSize',fontSize)
%     leg=legend('AP');
%     set(leg,'Interpreter','latex','FontSize',fontSize,'Location','Southeast')
    
%% Inserting the legend and saving the figure:
%     leg=legend([p1 p2 p3],'BS in the CmMIMO deployment','AP in the DmMIMO/grid deployment',...
%         'AP in the DmMIMO/linear deployment');
    leg=legend;
    % -0.0786    0.1072    0.9828    0.1406
    set(leg,'Interpreter','latex','FontSize',fontSize,'Location','Eastoutside')    
%     leg.Position=[0.05    0.1072    0.9828    0.1406];
    saveas(fig1,'deployments.png')
    saveas(fig1,'deployments.eps','epsc')
    