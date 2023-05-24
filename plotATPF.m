%% Reset
clearvars
close all
clc

%% Parameters
nu=25:25:125;
l=250;
d=0:1:l;

%% Gaussian Alarm Triggering Probabibilty Function:
gaussianATPF = @(d,sigma) exp(-(d.^2)/(2*sigma^2));

%% Colors for the plot:
blue = [57 106 177]./255;
red = [204 37 41]./255;
black = [83 81 84]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;
yellow = [240 180 50]./255;
colorMap=[black;red;green;blue;yellow;purple;brown];

%% Plot the curve
fig1=figure(1);
    fig1.Position=[100 300 550 400];
    hold on
    for s=1:length(nu)
        plot(d,gaussianATPF(d,nu(s)),'Color',colorMap(s,:),...
            'DisplayName',strcat('$\nu=\;$',num2str(nu(s))),...
            'LineWidth',2.5)
    end
    grid on
    leg=legend;
    set(leg,'Interpreter','latex','FontSize',18,'Location','Northeast');
    xlabel('$d_{k,a}$ [m]','Interpreter','latex','FontSize',18)
    ylabel('MTD Activation probability','Interpreter','latex','FontSize',18)
    set(gca,'TickLabelInterpreter','latex','FontSize',18)
    saveas(fig1,'gaussianATPF.png')
    saveas(fig1,'gaussianATPF.eps','epsc')