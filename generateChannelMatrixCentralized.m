function G=generateChannelMatrixCentralized(M,K,sigma2,l,h_BS,h_MTD,fc,MC,distributionMTDs)
    % Generate small scale fading samples:
    H=(1/sqrt(2))*(randn(M,K,MC)+1i*randn(M,K,MC));

    % Path-loss model:
    n_NLOS=3.19;                        % Path-loss exponent
    sigma_NLOS_dB=7.56;                 % Variance of the shadowing [dB]
    X_dB=sigma_NLOS_dB*randn(1,K,MC);   % Shadowing samples

    % Determining the position of the BS:
    dx_BS=l/2;              % x-position of the BS [m]
    dy_BS=l/2;              % x-position of the BS [m]
    
    % Determining the positions of the MTDs:
    if distributionMTDs=="Random"
        dx_MTD=l*rand(1,K);         % x-position of the MTDs [m]
        dy_MTD=l*rand(1,K);         % y-position of the MTDs [m]
    elseif distributionMTDs=="Gaussian"
        x_mean=l/4;                 % Mean value in the x-axis 
        y_mean=l/4;                 % Mean value in the y-axis
        stdDev=50;                  % Standard deviation in both axis
        dx_MTD=x_mean+stdDev*randn(1,K);   % x-position of the MTDs [m]
        dy_MTD=y_mean+stdDev*randn(1,K);   % y-position of the MTDs [m]
        % Guarantee that all MTDs are inside the square cell:
        for k=1:K
            if dx_MTD(k)<0
                dx_MTD(k)=0;
            end
            if dx_MTD(k)>l
                dx_MTD(k)=l;
            end
            if dy_MTD(k)<0
                dy_MTD(k)=0;
            end
            if dy_MTD(k)>l
                dy_MTD(k)=l;
            end
        end
    else
        error("Error! Please select a valid spatial distribution for the MTDs (random or Gaussian).")
    end
    

    % Computing the 3D distances between the BS and the MTDs:
    d_3D=zeros(1,K);
    for k=1:K
        d_3D(k)=sqrt((dx_BS-dx_MTD(k)).^2+(dy_BS-dy_MTD(k)).^2+(h_BS-h_MTD).^2);
    end

    % Attenuation due to the distance for all the K MTDs:
    A_dB=zeros(1,K);
    for k=1:K
        A_dB(k)=32.45+20*log10(fc)+10*n_NLOS*log10(d_3D(k));   
    end
    PL_dB=A_dB+X_dB;                    % Path-loss [dB]
    PL=10.^(PL_dB/10);                  % Path-loss [linear scale]
    beta=1./PL;                         % Large scale fading coefficient        

    % Normalize the large scale fading coefficients:
    beta=beta/sigma2;

    % Matrices of wireless channel gains:
    G=sqrt(beta).*H;
end