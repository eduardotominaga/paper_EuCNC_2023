function G=generateChannelMatrixDistributed(M,K,S,Q,sigma2,l,...
    h_AP,h_MTD,fc,deployment,MC,distributionMTDs)
    % Generate the small scale fading samples:
    H=(1/sqrt(2))*(randn(M,K,MC)+1i*randn(M,K,MC));

    % Path-loss model:
    n_NLOS=3.19;                        % Path-loss exponent
    sigma_NLOS_dB=7.56;                 % Variance of the shadowing [dB]
    X_dB=sigma_NLOS_dB*randn(1,K,MC);   % Shadowing samples

    % Determining the positions of the MTDs:
    if distributionMTDs=="Uniform"
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

    % Determining the positions of the APs:
    switch deployment
        case "Grid"
            range_AP=l/(2*sqrt(Q)):l/sqrt(Q):l-l/(2*sqrt(Q));     
            dx_AP=zeros(1,Q);           % x-position of the APs [m]
            dy_AP=zeros(1,Q);           % y-position of the APs [m]
            for q=1:sqrt(Q)
               dx_AP(1+sqrt(Q)*(q-1):sqrt(Q)*q)=range_AP(q); 
            end
            for q=1:sqrt(Q)
               dy_AP(1+sqrt(Q)*(q-1):sqrt(Q)*q)=range_AP; 
            end

        case "RadioStripes"
            range_AP=l/(2*(Q/4)):l/(Q/4):l-l/(2*(Q/4));   
            dx_AP=zeros(1,Q);           % x-position of the APs [m]
            dy_AP=zeros(1,Q);           % y-position of the APs [m]
            for idx=1:4                 % A square area has four sides.
                switch idx
                    case 1
                        dy_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=range_AP;
                    case 2
                        dx_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=range_AP;
                        dy_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=l*ones(1,Q/4);            
                    case 3
                        dx_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=l*ones(1,Q/4);
                        dy_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=flip(range_AP);
                    case 4
                        dx_AP(1+(Q/4)*(idx-1):(Q/4)*idx)=flip(range_AP);
                end
            end
        otherwise
            error("Error! Please select a valid deployment for the APs!")
    end

    % 3D distances between the APs and the MTDs:
    d_3D=zeros(Q,K);
    for k=1:K
        for q=1:Q
            d_3D(q,k)=sqrt((dx_AP(q)-dx_MTD(k))^2+(dy_AP(q)-dy_MTD(k))^2+(h_AP-h_MTD)^2);
        end
    end

    % Attenuation due to the distance:
    A_dB=zeros(Q,K);
    for k=1:K
        for q=1:Q
            A_dB(q,k)=32.45+20*log10(fc)+10*n_NLOS*log10(d_3D(q,k));
        end
    end

    PL_dB=A_dB+X_dB;                                % Path-loss [dB]
    PL=10.^(PL_dB/10);                              % Path-loss [linear scale]
    beta=1./PL;                                     % Attenuation coefficient        

    % Normalize the large scale fading coefficients:
    beta=beta/sigma2;

    % Generate the vectors of wireless channel gains:
    G=zeros(M,K,MC);
    for k=1:K
        for idx=1:MC
            for q=1:Q
                G(1+S*(q-1):S*q,k,idx)=sqrt(beta(q,k,idx))*H(1+S*(q-1):S*q,k,idx);
            end
        end
    end
end
