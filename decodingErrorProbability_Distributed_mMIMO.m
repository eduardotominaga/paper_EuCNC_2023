%% Description
% Uplink of a distributed mMIMO setup.
% Industrial path loss model validated by 3GPP.
% Compute the decoding error probability versus number of antenna elements.

%% Reset
clearvars;
close all;
clc

%% Parameters                     
K=16;                       % Number of MTDs
S=4;                        % Number of antennas at each AP
h_AP=6;                     % Heigth of the AP [m]
h_MTD=1.5;                  % Height of the MTDs [m]
l=250;                      % Length of the side of the square area [m]
R=1;                        % Target data rate [bits/s/Hz]
fc=3.5;                     % Carrier frequency [GHz]
rho_dBm=20;                 % Tranmit power of the MTD [dBm]
rho=10^((rho_dBm-30)/10);   % Tranmit power of the MTD [W]
B=20e6;                     % Bandwidth [Hz]
NF_dB=7;                    % Noise Figure [dB]
NF=10^(NF_dB/10);           % Noise Figure [linear scale]
N=1e2;                      % Number of network realizations
MC=1e4;                     % Number of channel realizations
decoding="MMSE";            % Select MRC, ZF or MMSE
rng(1);                     % Set the random number generator

% Select the deployment of APs:
% deployment="Grid";          
deployment="RadioStripes";

% Determine the number of APs:
if deployment=="Grid"
    Q=[4 9 16 25];
elseif deployment=="RadioStripes"
    Q=4:4:24;
end
M=Q*S;  % Total number of receive antennas 

% Define the spatial distribution of MTDs:
% distributionMTDs="Uniform";
distributionMTDs="Gaussian";

% Compute the noise variance:
sigma2_dB=-204+10*log10(B)+NF;
sigma2=10^(sigma2_dB/10);

%% Decoding error probability computation
epsilon=zeros(1,length(M));     % Vector of decoding error probabilities
       
% Loop 1: compute the outage probability for different numbers of
% active MTC devices:
parfor m=1:length(M)
    str=['Progress of the simulation: ',num2str(m),'/',num2str(length(M))];
    disp(str) 
    
    gamma=zeros(N,MC);       % Matrix of instantaneous SINR values    
    % Loop 2: average over several network realizations
    for n=1:N               

        % Generate the wireless channel matrix:
        G=generateChannelMatrixDistributed(M(m),K,S,Q(m),sigma2,l,...
            h_AP,h_MTD,fc,deployment,MC,distributionMTDs);

        % If the number of users is lower or equal than the number of antennas,...
        if K<=M(m)
            % Loop 3: compute the instantaneous SINR for each channel
            % realization:
            for mc=1:MC
                % Select the linear decoding method:
                switch decoding
                    case "MRC"      
                        b=G(:,1,mc);
                    case "ZF"
                        b=G(:,:,mc)*inv(G(:,:,mc)'*G(:,:,mc));
                        b=b(:,1);
                    case "MMSE"
                        b=inv(G(:,:,mc)*G(:,:,mc)'+(1/rho)*eye(M(m)))*G(:,:,mc);
                        b=b(:,1);
                    otherwise
                        error("Please select a valid linear decoding method!")
                end
                % Compute the interference from other users:
                I=0;
                for j=2:K
                    I=I+abs(b'*G(:,j,mc))^2;
                end
                % Instantaneous SINR:
                gamma(n,mc)=(rho*abs(b'*G(:,1,mc))^2)/(rho*I+(vecnorm(b)^2));
            end                
        end
    end
    % Compute the decoding error probability:
    epsilon(m)=mean(log2(1+vec(gamma))<R);
end  
    
%% Saving the results:
varName=matlab.lang.makeValidName(strcat(...
    "epsilon_Distributed_",deployment,"_",decoding,"_",distributionMTDs));
assignin('base',varName,epsilon);
fileName=strcat("Results_Distributed_",deployment,"_",...
    decoding,"_",distributionMTDs,".mat");
save(fileName,varName,'K','S','Q','M')

%% Plotting the results
       
fig=figure(1);
    fig.Position=[400 400 500 350];
    hold on
    plot(M,epsilon(1,:))
    set(gca,'yscale','log')
    grid on
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('$K$ (Number of MTC devices)','Interpreter','latex','FontSize',12)
    ylabel('Decoding Error Probability','Interpreter','latex','FontSize',12)
    ylim([1e-6 1])
    leg=legend('show');
    set(leg,'Interpreter','latex','FontSize',12,'Location','Northeast')
%     saveas(fig,'Figure.png')
%     saveas(fig,'Figure.eps','epsc')
     
%% This part of the code terminates all the Matlab processes is the script run on a server:
% if getenv('COMPUTERNAME')~="OY2106111"  % If this is not my personal computer...    
%     exit;                               % Terminate all the Matlab processes
% end
