clear all 
close all

% CODE: Linear Random Coefficient Model Estimation
% GOAL: Variance-Covariance Matrix using Within Dimensions Variability
% Coefficients 
% TYPE OF PANEL: Balanced 
% Use it when the panel is balanced, the errors are white noise and the
% regressors are exogeneous. 
% AUTHOR: Monika Avila M?rquez
% DATE: 10/03/18
% If you use this matlab code for estimation of a linear RCM for Three 
% Dimensional Panel Data please paste the following disclaimer: 
% 


%% Load Data 

filename='Peru_Model1.mat';    
load(filename)

%% FGLS Estimation 

% 2.1.1 Variance-Covariance Estimation: Within Dimensions Variability 

% Auxiliar Matrices for Method 1: Within Variation 

% Fixing i dimension 

    auMF1=(1:N2*T:N1*N2*T);    % To get the initial row of data to keep. Used from line 56 Whitin Code
    auMF2=(N2*T:N2*T:N1*N2*T); % To get the last row of data to keep. Used from line 56 Whitin Code

    % For fixing j

    au5=0:1:(N2-1); % Auxiliar matrix of positions 
    au6=1:1:(N2);   % Auxiliar matrix of positions 
    b1=ones(T,K);   % Matrix of ones(TxK)
    f2=ones(N1*N2*T*N2,K);  % Stores N2 different indicator matrices in one. Stacks N2 matrices
    % This is why it has N2*Sample Size. In each j loop takes the j sample size
    % matrix. 
    au7=1:N1*N2*T:(N1)*N2*T*N2; % Lowest row for the f2 to put the indicator matrix 
    au8=N1*N2*T:N1*N2*T:N2*(N1)*N2*T; % Highest row f 2 to put the indicator matrix 

%% Fixing i 

    % Storage matrices 
    
    estimatesfa=[]; % To store theta estimates 

    for i=1:N1
    
    % Projector and Emulator matrices for each group i   
    
    P=X(auMF1(1,i):auMF2(1,i),:)*(inv(X(auMF1(1,i):auMF2(1,i),:)'...
        *X(auMF1(1,i):auMF2(1,i),:)))*X(auMF1(1,i):auMF2(1,i),:)';  
        % Projection Matrix for each i
        
    M=eye(N2*T)-P;          % Orthogonal Projection Matrix for each i
    
    r=M*Y(auMF1(i):auMF2(i),:); % Residuals for each i              
    
    R=r(:,:)*r(:,:)';       % Outer product of residuals for all i's. 

    % H Matrices

    % H_gamma and H_lambda

    x_1=ones(N2*T,1);                       % Auxiliary variable 
    
    X_1=[x_1 X(auMF1(i):auMF2(i),:)];       % Auxiliary matrix 

    H_gamma1f=ones(N2*T,(K+1)*N2*T);        % Matrix for storage 
    
    H_gamma1=ones(N2*T,(K+1)*(K)*N2*T);     % Matrix with extra space for storage

    H_lambda1f=ones(N2*T,(K+1)*N2*T);       % Matrix for storage 
    
    H_lambda1=ones(N2*T,(K+1)*(K)*N2*T);    % Matrix with extra space for storage

    for k=2:(K+1) % Starts at 2 because we skip the first one that is not necessary

        aa=diag(X_1(:,k));% Put the column k from matrix X in diagonal 

        for kk=2:(K+1)

        l=(kk-1)*N2*T+1;    % Initial column position where to store
        
        u=(kk)*N2*T;        % Last column position where to store

        bb=diag(X_1(:,kk)); % Put the column k' from X in diagonal 

        H_gamma1f([1:end],[l:u])=(aa*kron(I_N2,iota_T))...
            *(kron(I_N2,iota_T')*bb');   % Getting H_re,kk',i
        
        H_gamma11f=H_gamma1f(:,[N2*T+1:end]);% Storage of matrix obtained without 
        % the auxiliar part

        H_lambda1f([1:end],[l:u])=(aa*kron(iota_N2,I_T))...
            *(kron(iota_N2',I_T)*bb');     % Getting H_re,kk',i
        
        H_lambda11f=H_lambda1f(:,[N2*T+1:end]);% Storage of matrix obtained
        % without the auxiliar part
        
        end
        
        ll=(k-1)*(K)*N2*T+1; % Initial column position where to store
        
        uu=(k)*(K)*N2*T;     % Last column position where to store

        H_gamma1([1:end],[ll:uu])=H_gamma11f(:,:);  % Storing all data for the nested loop

        H_lambda1([1:end],[ll:uu])=H_lambda11f(:,:);% Storing all data for the nested loop

    end 

    H_gamma=H_gamma1(:,[N2*T*K+1:end]);   % Matrix with H.kk 

    H_lambda=H_lambda1(:,[N2*T*K+1:end]); % Matrix with H.kk 


    %MH_gammaM and MH_lambdaM

    XX_1=eye(N2*T);                     % Auxiliar matrix 

    H_gammaa=[XX_1 H_gamma];            % Auxiliar matrix 
    
    MH_gamma1M=ones(N2*T,(K*K+1)*N2*T); % Matrix to store final values 

    H_lambdaa=[XX_1 H_lambda];          % Auxiliar matrix 
    
    MH_lambda1M=ones(N2*T,(K*K+1)*N2*T); % Matrix to store final values 

    for k=2:(K*K+1)
        
        MHM1=eye(N2*T); % Matrix to store, inside loop because I want to overwrite it 
        
        l=(k-1)*N2*T+1; % Initial column to store in memory mapping matrix 
        
        u=(k)*N2*T;     % Last column to store in memory mapping matrix 
        
        MHM1=M*H_gammaa([1:end],[l:u])*M; % M_re,kk'iM
        
        MH_gamma1M([1:end],[l:u])=MHM1(:,:); % Storing M_re,kk'iM

        MHM2=eye(N2*T); % Matrix to store, inside loop because I want to overwrite it 
        
        MHM2=M*H_lambdaa([1:end],[l:u])*M;  % M_re,kk'iM
        
        MH_lambda1M([1:end],[l:u])=MHM2(:,:);% Storing M_re,kk'iM

    end

    MH_gammaM=MH_gamma1M([1:end],[N2*T+1:end]);  % Final matrix containing values 
    
    MH_lambdaM=MH_lambda1M([1:end],[N2*T+1:end]);% Final matrix containing values 

    % Matrices  B, C

    B1=ones(N2*T*N2*T,K*K); % Matrix to store values

    C1=ones(N2*T*N2*T,K*K); % Matrix to store values
    
    for j=2:K*K+1
        
        BB=vec(MH_gamma1M([1:end],[(j-1)*N2*T+1:j*N2*T])); % Obtaining vec(MH_kk,reiM)
        
        B1(:,j)=BB(:,:); % Storing values 

        CC=vec(MH_lambda1M([1:end],[(j-1)*N2*T+1:j*N2*T]));% Obtaining vec(MH_kk,reiM)
        C1(:,j)=CC(:,:); % Storing values 

    end

    B=B1(:,[2:end]); % Final matrix B

    C=C1(:,[2:end]); % Final matrix C
    
    %% Z 

    BS=B*D;
    CS=C*D;
    m=vec(M);
    Z=horzcat( BS, CS, m);

    % Errors 

    R1=vec(R);

    estimates_mf=(Z'*Z)\(Z'*R1);
    
    % Estimates fixing alpha
    
    estimatesfa=[estimatesfa,estimates_mf];
    
    end 

    % Average Estimates of gamma, lambda and sigma

    averagest=mean(estimatesfa');

%% Fixing j 


    for j=1:N2
    % Getting the N1*T observations corresponding to j=i
    % Positions 
     % Creating the matrix of positions that we have to extract
        
        a1=0*ones(N1*N2*T,K);
        H_A=[1+T*au5(j):N2*T:(N1)*N2*T+T*au5(j)+1]; % Initial position to put indicator vector
        G=[T*au6(j):N2*T:(N1)*N2*T+T*au6(j)];       % Last position to put indicator matrix

   % Creating a matrix of O's and 1's. 1 for extracting the values 
   % and 0 to delete values that we don't need. The same structure repeats
   % for each individual in dimension i. This is why we repeat N1 times with a loop.  
        
    for l=1:N1 

        a1(H_A(l):G(l),:)=b1(:,:); % a1 has the size of whole matrix.
  

    end
    
     % f2 Stores N2 different indicator matrices in one. Stacks N2 matrices
     % This is why it has N2*Sample Size. In each j loop takes the j sample size
     % matrix. Replaces with the new indicator matrix of positions 
    
    f2([au7(j):au8(j)],:)=a1(:,:);

    %This the matrix of values that we used in iteration 
    
    %Xx=X;%Copying matrix of original data 
    
    %Xx(f2([au7(j):au8(j)],:)~=1)=NaN; % Replacing with NaN values that we don't need  
    
    %aa=0*ones(N1*T,K);%Storage matrix 
    
    %Extracting the values that we will need 
    
    for i=1:K % We do a loop to get all K columns 
        
    Xau=X(:,i); % Extracting column i of the data of explanatory vbles   
    
    Xa(:,i)=Xau(f2([au7(j):au8(j)],i)~=0);% Extracts values that we should keep 
    % for each k signaled by 1 in the indicator matrix created above 
    
    Ya=Y(f2([au7(j):au8(j)],i)~=0);
    
    end
    
    P2=Xa*(pinv(Xa'*Xa))*Xa';           % Projection Matrix 
    
    M2=eye(N1*T)-P2;  % Residual Matrix 
    
    r2=M2*Ya;                           % Residuals 
    
    R2=r2*r2';                         % Outer product of residuals 

    % H Matrices 

    %H_alpha

    x_11=ones(N1*T,1); % Auxiliary variable 
    X_11=[x_11 Xa];    % Auxiliary matrix 

    H_alpha1f1=ones(N1*T,(K+1)*N1*T);     % Matrix for storage 
    H_alpha11=ones(N1*T,(K+1)*(K)*N1*T);  % Matrix with extra space for storage

    H_lambda1f1=ones(N1*T,(K+1)*N1*T);    % Matrix for storage 
    H_lambda11=ones(N1*T,(K+1)*(K)*N1*T); % Matrix with extra space for storage

    for k=2:(K+1) % Starts at 2 because we skip the first one that is not necessary

        aa1=diag(X_11(:,k)); %Put the column in diagonal 

        for kk=2:(K+1)

        l=(kk-1)*N1*T+1; % Position where to store
        u=(kk)*N1*T;

        bb1=diag(X_11(:,kk)); % Put the column in diagonal 

        H_alpha1f1([1:end],[l:u])=(aa1*superkron(I_N1,iota_T))...
            *(superkron(I_N1,iota_T')*bb1'); % Get the matrix H-kk we want;

        H_alpha11f1=H_alpha1f1(:,[N1*T+1:end]); % Gettin only the matrix that we want


        H_lambda1f1([1:end],[l:u])=(aa1*superkron(iota_N1,I_T))...
            *(superkron(iota_N1',I_T)*bb1'); % Storing in the matrix for j loop

        H_lambda11f1=H_lambda1f1(:,[N1*T+1:end]);%Gettin only the matrix that we want

        end
        ll=(k-1)*(K)*N1*T+1; % Position where to store
        uu=(k)*(K)*N1*T;

        H_alpha11([1:end],[ll:uu])=H_alpha11f1(:,:);%Storing all data for the nested loop

        H_lambda11([1:end],[ll:uu])=H_lambda11f1(:,:);%Storing all data for the nested loop

    end 

    H_alpha1=H_alpha11(:,[N1*T*K+1:end]); % Matrix with H.kk 

    H_lambda1=H_lambda11(:,[N1*T*K+1:end]); % Matrix with H.kk 

    %MH_alphaM 

    XX_1=eye(N1*T);
    
    H_alphaa1=[XX_1 H_alpha1];
    
    H_lambdaa1=[XX_1 H_lambda1];
    
    MH_alpha1M1=ones(N1*T,(K*K+1)*N1*T);
    
    MH_lambda1M1=ones(N1*T,(K*K+1)*N1*T);
    
    for k=2:(K*K+1)
        
        l=(k-1)*N1*T+1;
        u=(k)*N1*T;

        MHM12=eye(N1*T);
        MHM12=M2*H_alphaa1([1:end],[l:u])*M2;
        MH_alpha1M1([1:end],[l:u])=MHM12(:,:);

        MHM22=eye(N1*T);
        MHM22=M2*H_lambdaa1([1:end],[l:u])*M2;
        MH_lambda1M1([1:end],[l:u])=MHM22(:,:);
    end

    MH_alphaM1=MH_alpha1M1([1:end],[N1*T+1:end]);
    MH_lambdaM1=MH_lambda1M1([1:end],[N1*T+1:end]);

    % Matrices  A, C

    A2_1=ones(N1*T*N1*T,K*K);
    C2_1=ones(N1*T*N1*T,K*K);
    for i=2:K*K+1

        A2_1(:,i)=vec(MH_alpha1M1([1:end],[(i-1)*N1*T+1:i*N1*T]));

        C2_1(:,i)=vec(MH_lambda1M1([1:end],[(i-1)*N1*T+1:i*N1*T]));
    end

    A2=A2_1(:,[2:end]);

    C2=C2_1(:,[2:end]);
    %% Z 

    AS1=A2*D;
    CS1=C2*D;
    m1=vec(M2);
    Z1=horzcat( AS1, CS1, m1);

    % Errors 

    R2=vec(R2);

    estimates2(:,j)=pinv((Z1'*Z1))*(Z1'*R2);

    end 

    % Estimates of alpha,lamda and sigma

    averagest2=mean(estimates2');

%% Fixing t 
 
    % Auxiliar matrices for obtaining design matrices
    
    au9=0:1:(T-1);
    au10=1:1:(T);

    b3=1*ones(1,K); % Vector of 1's. Indicating values to keep 

    for t=1:T

      f3=0*ones(N1*N2*T,K);

      % Creating indicating matrix with 1's in positions to save 

      au11_W=(t:T:N1*N2*T);  

      for i=1:N1*N2

        f3(au11_W(1,i),:)=b3(:,:); 

      end 

    %Extracting the values that we will need 
     
    for i=1:K % We do a loop to get all K columns 
    Xau1=X(:,i);   
    Xa1(:,i)=Xau1(f3(:,i)~=0);% Column i with values that we keep 
    Ya1=Y(f3(:,i)~=0);

    end

    P3=Xa1*(pinv(Xa1'*Xa1))*Xa1';            % Projector matrix 
    M3=eye(N1*N2)-P3; % Residual matrix 
    r3=M3*Ya1;                              % Residuals 
    R3=r3*r3';                              % Outer product of residuals 


    % H

    x_2=ones(N1*N2,1);% Auxiliary variable 
    X_2=[x_2 Xa1];    % Auxiliary matrix 
    
    H_alpha1f2=ones(N1*N2,(K+1)*N1*N2);   % Matrix for storage 
    H_alpha12=ones(N1*N2,(K+1)*(K)*N1*N2);% Matrix with extra space for storage

    H_gamma1f2=ones(N1*N2,(K+1)*N1*N2);   % Matrix for storage 
    H_gamma12=ones(N1*N2,(K+1)*(K)*N1*N2);% Matrix with extra space for storage

    for k=2:(K+1) % Starts at 2 because we skip the first one that is not necessary

        aa2=diag(X_2(:,k)); % Put the column in diagonal 
        
        for kk=2:(K+1)

        l=(kk-1)*N1*N2+1; % Position where to store
        u=(kk)*N1*N2;

        bb2=diag(X_2(:,kk)); % Put the column in diagonal 

        H_alpha1f2([1:end],[l:u])=(aa2*superkron(I_N1,iota_N2))...
            *(superkron(I_N1,iota_N2')*bb2'); % Storing in the matrix for j loop

        H_alpha11f2=H_alpha1f2(:,[N1*N2+1:end]); % Gettin only the matrix that we want

        H_gamma1f2([1:end],[l:u])=(aa2*superkron(iota_N1,I_N2))...
            *(superkron(iota_N1',I_N2)*bb2'); % Storing in the matrix for j loop
        H_gamma11f2=H_gamma1f2(:,[N1*N2+1:end]); % Gettin only the matrix that we want


        end

        ll=(k-1)*(K)*N1*N2+1; % Position where to store
        uu=(k)*(K)*N1*N2;

        H_alpha12([1:end],[ll:uu])=H_alpha11f2(:,:); % Storing all data for the nested loop
        H_gamma12([1:end],[ll:uu])=H_gamma11f2(:,:); % Storing all data for the nested loop
    end 

    H_alpha2=H_alpha12(:,[N1*N2*K+1:end]); % Matrix with H.kk 
    H_gamma2=H_gamma12(:,[N1*N2*K+1:end]); % Matrix with H.kk 

    %MHM 

    XX_1=eye(N1*N2);
    H_alphaa2=[XX_1 H_alpha2];
    MH_alpha1M2=ones(N1*N2,(K*K+1)*N1*N2);

    H_gammaa2=[XX_1 H_gamma2];
    MH_gamma1M2=ones(N1*N2,(K*K+1)*N1*N2);

    for k=2:(K*K+1)

        l=(k-1)*N1*N2+1;
        u=(k)*N1*N2;

        MHM3=eye(N1*N2);
        MHM3=M3*H_alphaa2([1:end],[l:u])*M3;
        MH_alpha1M2([1:end],[l:u])=MHM3(:,:);

        MHM4=eye(N1*N2);
        MHM4=M3*H_gammaa2([1:end],[l:u])*M3;
        MH_gamma1M2([1:end],[l:u])=MHM4(:,:);
    end

    MH_alphaM2=MH_alpha1M2([1:end],[N1*T+1:end]);
    MH_lambdaM2=MH_gamma1M2([1:end],[N1*N2+1:end]);


    % Matrices  A, B

    %A
    
    A3_1=ones(N1*N2*N1*N2,K*K);
    B3_1=ones(N1*N2*N1*N2,K*K);
    
    for i=2:K*K+1

        A3_1(:,i)=vec(MH_alpha1M2([1:end],[(i-1)*N1*N2+1:i*N1*N2]));


        B3_1(:,i)=vec(MH_gamma1M2([1:end],[(i-1)*N1*N2+1:i*N1*N2]));
    end

    A3=A3_1(:,[2:end]);

    B3=B3_1(:,[2:end]);

    %% Z 

    AS3=A3*D;
    BS3=B3*D;
    m3=vec(M3);
    Z3=horzcat( AS3, BS3, m3);

    % Errors 

    R3=vec(R3);

    estimates3(:,t)=(Z3'*Z3)\(Z3'*R3);

end 

%% Average estimates of variance-covariance elements of alpha, gamma, lambda and epsilon

averagest3=mean(estimates3');

sigmahat=[averagest(1,2*K*(K+1)/2+1); averagest2(1,2*K*(K+1)/2+1); averagest3(1,K*(K+1)/2)];
sigma_h=mean(sigmahat)';

alphahat=[ averagest2(1,[1:K*(K+1)/2]); averagest3(1,[1:K*(K+1)/2])];
alpha_h=mean(alphahat)';

gammahat=[ averagest(1,[1:K*(K+1)/2]); averagest3(1,[K*(K+1)/2+1:2*K*(K+1)/2])];
gamma_h=mean(gammahat)';

lambdahat=[ averagest2(1,[K*(K+1)/2+1:2*K*(K+1)/2]); averagest(1,[K*(K+1)/2+1:2*K*(K+1)/2])];
lambda_h=mean(lambdahat)';

% 2.1.2 Coefficient Estimation 
    
theta_alpha2=alpha_h;
theta_gamma2=gamma_h;
theta_lambda2=lambda_h;
theta_sigma2=sigma_h;


[Omega_theta]=covariance(theta_sigma2,theta_alpha2,theta_gamma2,theta_lambda2,...
              K,D,X1,X2,X3,I_N1,I_N2,I_T,I_N1N2T);

Omega_theta_inv=inv(Omega_theta);   % Inverse of the estimated Omega 

dO=det(Omega_theta); % Determinant of the estimated Omega 

est_coeff=(X'*Omega_theta_inv*X)\(X'*Omega_theta_inv*Y) ; % Beta FGLS 

est_varcov=[alpha_h;gamma_h;lambda_h;sigma_h];

% 2.1.3 Inference 

var_cov_coeff=X'*Omega_theta_inv*X;

% To obtain t_statistics for significance test you need to divide the
% estimated coefficient by the s.e of the coefficient. This s.e. is the
% squared root of the variance of the coefficient that you can find in the
% diagonal of the variance-covariance of the estimated coefficients.
% (var_cov_coeff)
% e.g.
t_B1=est_coeff(2,1)/(sqrt(var_cov_coeff(2,2)))

save('Estimates_FGLS_Within.mat','est_coeff','est_varcov','var_cov_coeff');
