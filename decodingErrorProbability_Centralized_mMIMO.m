%% Description
% Uplink of a centralized mMIMO setup.
% Industrial large-scale fading model validated by 3GPP
% Compute the decoding error probability versus number of active MTDs.

%% Reset
clearvars;
close all;
clc

%% Parameters
M=16:16:96;                 % Number of receive antennas                      
K=16;                       % Number of MTDs
h_BS=6;                     % Heigth of the BS [m]
h_MTD=1.5;                  % Height of the MTDs [m]
l=250;                      % Length of the side of the square area [m]
R=1;                        % Target data rate [bits/s/Hz]
fc=3.5;                     % Carrier frequency [GHz]
rho_dBm=20;                 % Tranmit power of the MTD [dBm]
rho=10^((rho_dBm-30)/10);    % Tranmit power of the MTD [W]
B=20e6;                     % Bandwidth [Hz]
NF_dB=7;                    % Noise Figure [dB]
NF=10^(NF_dB/10);           % Noise Figure [linear scale]
N=1e2;                      % Number of network realizations
MC=1e3;                     % Number of channel realizations
decoding="MMSE";            % Select MRC, ZF or MMSE
rng(1);                     % Set the random number generator

% Define the spatial distribution of MTDs:
% distributionMTDs="Random";
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
    gamma=zeros(N,MC);      % Matrix of instantaneous SINR values
    % Loop 2: average over several network realizations
    for n=1:N

        % Generate the matrix of wireless channel gains:
        G=generateChannelMatrixCentralized(M(m),K,sigma2,...
            l,h_BS,h_MTD,fc,MC,distributionMTDs);

        % If the number of users is lower or equal than the number of antennas,...
        if K<=M(m)
            % Loop 3: compute the instantaneous SINR for each channel
            % realization:
            for mc=1:MC
                % Select the linear decoding method:
                switch decoding
                    case "MRC"      
                        b=G(:,1,mc);   % Receive beamforming vector
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
    "epsilon_Centralized_",decoding,"_",distributionMTDs));
assignin('base',varName,epsilon);
fileName=strcat("Results_Centralized_",decoding,...
    "_",distributionMTDs,".mat");
save(fileName,varName,'K','M')

%% Plotting the results
        
fig7=figure(1);
    fig7.Position=[400 400 500 350];
    hold on
    plot(M,epsilon(1,:))
    set(gca,'yscale','log')
    grid on
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('$M$','Interpreter','latex','FontSize',12)
    ylabel('Decoding Error Probability','Interpreter','latex','FontSize',12)
    ylim([1e-6 1])
    leg=legend('show');
    set(leg,'Interpreter','latex','FontSize',12,'Location','Northeast')
     
%% This part of the code terminates all the Matlab processes is the script run on a server:
% if getenv('COMPUTERNAME')~="OY2106111"  % If this is not my personal computer...    
%     exit;                               % Terminate all the Matlab processes
% end
