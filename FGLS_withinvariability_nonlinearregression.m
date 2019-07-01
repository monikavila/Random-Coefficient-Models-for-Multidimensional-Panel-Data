clear all 
close all

% CODE: Linear Random Coefficient Model Estimation
% GOAL: Variance-Covariance Matrix using Within Dimensions Variability
% Coefficients 
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
    estimatesfg=[];
    estimatesfl=[];
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

    
    %% Non linear estimation - Gradient Descent 

    % Initializing 

    theta=rand(2*K2+1,1); 
   
    cost=[]
    
    step=0.000000000000000000001;
  
    i=0
    
    for ii=1:100
    % Gradient 

    
    %% IMPORTANT: THIS HAS TO BE ADJUSTED MANUALLY ACCORDING TO THE NUMBER OF VARIABLES 
    %% ADDED TO THE MODEL. (WE ESTIMATE: 2*K*(K+1)/2) 
    % We use a quadratic reparamterization of the diagonal elements. 
    % i.e. When K=7,we estimate 57 variance-covariance elements. The
    % elements that are reparamterized to the quadratic funtion are: 1,8,
    
    vu=[theta(1,1)^2;theta(2,1);theta(3,1);theta(4,1);theta(5,1);theta(6,1);theta(7,1);
        theta(8,1)^2;theta(9,1);theta(10,1);theta(11,1);theta(12,1);theta(13,1);
        theta(14,1)^2;theta(15,1);theta(16,1);theta(17,1);theta(18,1);
        theta(19,1)^2;theta(20,1);theta(21,1);theta(22,1);
        theta(23,1)^2;theta(24,1);theta(25,1);
        theta(26,1)^2;theta(27,1);
        theta(28,1)^2;
        theta(29,1)^2;theta(30,1);theta(31,1);theta(32,1);theta(33,1);theta(34,1);theta(35,1);
        theta(36,1)^2;theta(37,1);theta(38,1);theta(39,1);theta(40,1);theta(41,1);
        theta(42,1)^2;theta(43,1);theta(44,1);theta(45,1);theta(46,1);
        theta(47,1)^2;theta(48,1);theta(49,1);theta(50,1);
        theta(51,1)^2;theta(52,1);theta(53,1);
        theta(54,1)^2;theta(55,1);
        theta(56,1)^2;
        theta(57,1)^2];

    grad_vu=diag([2*theta(1,1);theta(2,1)/abs(theta(2,1));theta(3,1)/abs(theta(3,1));theta(4,1)/abs(theta(4,1));theta(5,1)/abs(theta(5,1));theta(6,1)/abs(theta(6,1));theta(7,1)/abs(theta(7,1));
                  2*theta(8,1);theta(9,1)/abs(theta(9,1));theta(10,1)/abs(theta(10,1));theta(11,1)/abs(theta(11,1));theta(12,1)/abs(theta(12,1));theta(13,1)/abs(theta(13,1));
                  2*theta(14,1);theta(15,1)/abs(theta(15,1));theta(16,1)/abs(theta(16,1));theta(17,1)/abs(theta(17,1));theta(18,1)/abs(theta(18,1));
                  2*theta(19,1);theta(20,1)/abs(theta(20,1));theta(21,1)/abs(theta(21,1));theta(22,1)/abs(theta(22,1));
                  2*theta(23,1);theta(24,1)/abs(theta(24,1));theta(25,1)/abs(theta(25,1));
                  2*theta(26,1);theta(27,1)/abs(theta(27,1));
                  2*theta(28,1);
                  2*theta(29,1);theta(30,1)/abs(theta(30,1));theta(31,1)/abs(theta(31,1));theta(32,1)/abs(theta(32,1));theta(33,1)/abs(theta(33,1));theta(34,1)/abs(theta(34,1));theta(35,1)/abs(theta(35,1));
                  2*theta(36,1);theta(37,1)/abs(theta(37,1));theta(38,1)/abs(theta(38,1));theta(39,1)/abs(theta(39,1));theta(40,1)/abs(theta(40,1));theta(41,1)/abs(theta(41,1));
                  2*theta(42,1);theta(43,1)/abs(theta(43,1));theta(44,1)/abs(theta(44,1));theta(45,1)/abs(theta(45,1));theta(46,1)/abs(theta(46,1));
                  2*theta(47,1);theta(48,1)/abs(theta(48,1));theta(49,1)/abs(theta(49,1));theta(50,1)/abs(theta(50,1));
                  2*theta(51,1);theta(52,1)/abs(theta(52,1));theta(53,1)/abs(theta(53,1));
                  2*theta(54,1);theta(55,1)/abs(theta(55,1));
                  2*theta(56,1);
                  2*theta(57,1);]);    

    gradient=(grad_vu*Z'*(R1-Z*vu))*(-1);

    cost1=R1-Z*vu;

    cost1=cost1'*cost1;

    % Gradient descent 

    theta1=theta-step*gradient;
    vu1=[theta1(1,1)^2;theta1(2,1);theta1(3,1);theta1(4,1);theta1(5,1);theta1(6,1);theta1(7,1);
        theta1(8,1)^2;theta1(9,1);theta1(10,1);theta1(11,1);theta1(12,1);theta1(13,1);
        theta1(14,1)^2;theta1(15,1);theta1(16,1);theta1(17,1);theta1(18,1);
        theta1(19,1)^2;theta1(20,1);theta1(21,1);theta1(22,1);
        theta1(23,1)^2;theta1(24,1);theta1(25,1);
        theta1(26,1)^2;theta1(27,1);
        theta1(28,1)^2;
        theta1(29,1)^2;theta1(30,1);theta1(31,1);theta1(32,1);theta1(33,1);theta1(34,1);theta1(35,1);
        theta1(36,1)^2;theta1(37,1);theta1(38,1);theta1(39,1);theta1(40,1);theta1(41,1);
        theta1(42,1)^2;theta1(43,1);theta1(44,1);theta1(45,1);theta1(46,1);
        theta1(47,1)^2;theta1(48,1);theta1(49,1);theta1(50,1);
        theta1(51,1)^2;theta1(52,1);theta1(53,1);
        theta1(54,1)^2;theta1(55,1);
        theta1(56,1)^2;
        theta1(57,1)^2];

    cost2=R1-Z*vu1;
    cost2=cost2'*cost2;

    if cost2<cost1
    theta=theta1;
    else 
    theta=theta;
    end 
    % Cost function evaluation 
    vu3=[theta(1,1)^2;theta(2,1);theta(3,1);theta(4,1);theta(5,1);theta(6,1);theta(7,1);
        theta(8,1)^2;theta(9,1);theta(10,1);theta(11,1);theta(12,1);theta(13,1);
        theta(14,1)^2;theta(15,1);theta(16,1);theta(17,1);theta(18,1);
        theta(19,1)^2;theta(20,1);theta(21,1);theta(22,1);
        theta(23,1)^2;theta(24,1);theta(25,1);
        theta(26,1)^2;theta(27,1);
        theta(28,1)^2;
        theta(29,1)^2;theta(30,1);theta(31,1);theta(32,1);theta(33,1);theta(34,1);theta(35,1);
        theta(36,1)^2;theta(37,1);theta(38,1);theta(39,1);theta(40,1);theta(41,1);
        theta(42,1)^2;theta(43,1);theta(44,1);theta(45,1);theta(46,1);
        theta(47,1)^2;theta(48,1);theta(49,1);theta(50,1);
        theta(51,1)^2;theta(52,1);theta(53,1);
        theta(54,1)^2;theta(55,1);
        theta(56,1)^2;
        theta(57,1)^2];

    cost3=(R1-Z*vu3)'*(R1-Z*vu3);

    cost=[cost; cost3];
    
    ii=ii+1
    end
    figure(1)
    
    plot(cost)
    
    estimates=theta;
    estimates_mf=[estimates(1,1)^2;estimates(2,1);estimates(3,1);estimates(4,1);estimates(5,1);estimates(6,1);estimates(7,1);
                  estimates(8,1)^2;estimates(9,1);estimates(10,1);estimates(11,1);estimates(12,1);estimates(13,1);
                  estimates(14,1)^2;estimates(15,1);estimates(16,1);estimates(17,1);estimates(18,1);
                  estimates(19,1)^2;estimates(20,1);estimates(21,1);estimates(22,1);
                  estimates(23,1)^2;estimates(24,1);estimates(25,1);
                  estimates(26,1)^2;estimates(27,1);
                  estimates(28,1)^2;
                  estimates(29,1)^2;estimates(30,1);estimates(31,1);estimates(32,1);estimates(33,1);estimates(34,1);estimates(35,1);
                  estimates(36,1)^2;estimates(37,1);estimates(38,1);estimates(39,1);estimates(40,1);estimates(41,1);
                  estimates(42,1)^2;estimates(43,1);estimates(44,1);estimates(45,1);estimates(46,1);
                  estimates(47,1)^2;estimates(48,1);estimates(49,1);estimates(50,1);
                  estimates(51,1)^2;estimates(52,1);estimates(53,1);
                  estimates(54,1)^2;estimates(55,1);
                  estimates(56,1)^2;
                  estimates(57,1)^2];
   
    
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

   
     %% Non linear estimation - Gradient Descent 

    % Initializing 

    theta=rand(2*K2+1,1); 
   
    cost=[]
    step=0.0000000001;
    i=0
    
    for ii=1:100
    % Gradient 

    
    %% IMPORTANT: THIS HAS TO BE ADJUSTED MANUALLY ACCORDING TO THE NUMBER OF VARIABLES 
    %% ADDED TO THE MODEL. (WE ESTIMATE: 2*K*(K+1)/2) 
    % We use a quadratic reparamterization of the diagonal elements. 
    % i.e. When K=7,we estimate 57 variance-covariance elements. The
    % elements that are reparamterized to the quadratic funtion are: 1,8,
    
    vu=[theta(1,1)^2;theta(2,1);theta(3,1);theta(4,1);theta(5,1);theta(6,1);theta(7,1);
        theta(8,1)^2;theta(9,1);theta(10,1);theta(11,1);theta(12,1);theta(13,1);
        theta(14,1)^2;theta(15,1);theta(16,1);theta(17,1);theta(18,1);
        theta(19,1)^2;theta(20,1);theta(21,1);theta(22,1);
        theta(23,1)^2;theta(24,1);theta(25,1);
        theta(26,1)^2;theta(27,1);
        theta(28,1)^2;
        theta(29,1)^2;theta(30,1);theta(31,1);theta(32,1);theta(33,1);theta(34,1);theta(35,1);
        theta(36,1)^2;theta(37,1);theta(38,1);theta(39,1);theta(40,1);theta(41,1);
        theta(42,1)^2;theta(43,1);theta(44,1);theta(45,1);theta(46,1);
        theta(47,1)^2;theta(48,1);theta(49,1);theta(50,1);
        theta(51,1)^2;theta(52,1);theta(53,1);
        theta(54,1)^2;theta(55,1);
        theta(56,1)^2;
        theta(57,1)^2];

    grad_vu=diag([2*theta(1,1);theta(2,1)/abs(theta(2,1));theta(3,1)/abs(theta(3,1));theta(4,1)/abs(theta(4,1));theta(5,1)/abs(theta(5,1));theta(6,1)/abs(theta(6,1));theta(7,1)/abs(theta(7,1));
                  2*theta(8,1);theta(9,1)/abs(theta(9,1));theta(10,1)/abs(theta(10,1));theta(11,1)/abs(theta(11,1));theta(12,1)/abs(theta(12,1));theta(13,1)/abs(theta(13,1));
                  2*theta(14,1);theta(15,1)/abs(theta(15,1));theta(16,1)/abs(theta(16,1));theta(17,1)/abs(theta(17,1));theta(18,1)/abs(theta(18,1));
                  2*theta(19,1);theta(20,1)/abs(theta(20,1));theta(21,1)/abs(theta(21,1));theta(22,1)/abs(theta(22,1));
                  2*theta(23,1);theta(24,1)/abs(theta(24,1));theta(25,1)/abs(theta(25,1));
                  2*theta(26,1);theta(27,1)/abs(theta(27,1));
                  2*theta(28,1);
                  2*theta(29,1);theta(30,1)/abs(theta(30,1));theta(31,1)/abs(theta(31,1));theta(32,1)/abs(theta(32,1));theta(33,1)/abs(theta(33,1));theta(34,1)/abs(theta(34,1));theta(35,1)/abs(theta(35,1));
                  2*theta(36,1);theta(37,1)/abs(theta(37,1));theta(38,1)/abs(theta(38,1));theta(39,1)/abs(theta(40,1));theta(41,1)/abs(theta(41,1));theta(42,1)/abs(theta(42,1));
                  2*theta(42,1);theta(43,1)/abs(theta(43,1));theta(44,1)/abs(theta(44,1));theta(45,1)/abs(theta(45,1));theta(46,1)/abs(theta(46,1));
                  2*theta(47,1);theta(48,1)/abs(theta(48,1));theta(49,1)/abs(theta(49,1));theta(50,1)/abs(theta(50,1));
                  2*theta(51,1);theta(52,1)/abs(theta(52,1));theta(53,1)/abs(theta(53,1));
                  2*theta(54,1);theta(55,1)/abs(theta(55,1));
                  2*theta(56,1);
                  2*theta(57,1);]);    

    gradient=(grad_vu*Z'*(R1-Z*vu))*(-1);

    cost1=R1-Z*vu;

    cost1=cost1'*cost1;

    % Gradient descent 

    theta1=theta-step*gradient;
    vu1=[theta1(1,1)^2;theta1(2,1);theta1(3,1);theta1(4,1);theta1(5,1);theta1(6,1);theta1(7,1);
        theta1(8,1)^2;theta1(9,1);theta1(10,1);theta1(11,1);theta1(12,1);theta1(13,1);
        theta1(14,1)^2;theta1(15,1);theta1(16,1);theta1(17,1);theta1(18,1);
        theta1(19,1)^2;theta1(20,1);theta1(21,1);theta1(22,1);
        theta1(23,1)^2;theta1(24,1);theta1(25,1);
        theta1(26,1)^2;theta1(27,1);
        theta1(28,1)^2;
        theta1(29,1)^2;theta1(30,1);theta1(31,1);theta1(32,1);theta1(33,1);theta1(34,1);theta1(35,1);
        theta1(36,1)^2;theta1(37,1);theta1(38,1);theta1(39,1);theta1(40,1);theta1(41,1);
        theta1(42,1)^2;theta1(43,1);theta1(44,1);theta1(45,1);theta1(46,1);
        theta1(47,1)^2;theta1(48,1);theta1(49,1);theta1(50,1);
        theta1(51,1)^2;theta1(52,1);theta1(53,1);
        theta1(54,1)^2;theta1(55,1);
        theta1(56,1)^2;
        theta1(57,1)^2];

    cost2=R1-Z*vu1;
    cost2=cost2'*cost2;

    if cost2<cost1
    theta=theta1;
    else 
    theta=theta;
    end 
    % Cost function evaluation 
    vu3=[theta(1,1)^2;theta(2,1);theta(3,1);theta(4,1);theta(5,1);theta(6,1);theta(7,1);
        theta(8,1)^2;theta(9,1);theta(10,1);theta(11,1);theta(12,1);theta(13,1);
        theta(14,1)^2;theta(15,1);theta(16,1);theta(17,1);theta(18,1);
        theta(19,1)^2;theta(20,1);theta(21,1);theta(22,1);
        theta(23,1)^2;theta(24,1);theta(25,1);
        theta(26,1)^2;theta(27,1);
        theta(28,1)^2;
        theta(29,1)^2;theta(30,1);theta(31,1);theta(32,1);theta(33,1);theta(34,1);theta(35,1);
        theta(36,1)^2;theta(37,1);theta(38,1);theta(39,1);theta(40,1);theta(41,1);
        theta(42,1)^2;theta(43,1);theta(44,1);theta(45,1);theta(46,1);
        theta(47,1)^2;theta(48,1);theta(49,1);theta(50,1);
        theta(51,1)^2;theta(52,1);theta(53,1);
        theta(54,1)^2;theta(55,1);
        theta(56,1)^2;
        theta(57,1)^2];

    cost3=(R1-Z*vu3)'*(R1-Z*vu3);

    cost=[cost; cost3];
    
    
    
    ii=ii+1
    end
    figure(2)
    
    plot(cost)
    estimates=theta;
    estimates_mf=[estimates(1,1)^2;estimates(2,1);estimates(3,1);estimates(4,1);estimates(5,1);estimates(6,1);estimates(7,1);
                  estimates(8,1)^2;estimates(9,1);estimates(10,1);estimates(11,1);estimates(12,1);estimates(13,1);
                  estimates(14,1)^2;estimates(15,1);estimates(16,1);estimates(17,1);estimates(18,1);
                  estimates(19,1)^2;estimates(20,1);estimates(21,1);estimates(22,1);
                  estimates(23,1)^2;estimates(24,1);estimates(25,1);
                  estimates(26,1)^2;estimates(27,1);
                  estimates(28,1)^2;
                  estimates(29,1)^2;estimates(30,1);estimates(31,1);estimates(32,1);estimates(33,1);estimates(34,1);estimates(35,1);
                  estimates(36,1)^2;estimates(37,1);estimates(38,1);estimates(39,1);estimates(40,1);estimates(41,1);
                  estimates(42,1)^2;estimates(43,1);estimates(44,1);estimates(45,1);estimates(46,1);
                  estimates(47,1)^2;estimates(48,1);estimates(49,1);estimates(50,1);
                  estimates(51,1)^2;estimates(52,1);estimates(53,1);
                  estimates(54,1)^2;estimates(55,1);
                  estimates(56,1)^2;
                  estimates(57,1)^2];
   
   
    
    % Estimates fixing gamma
    
    estimatesfg=[estimatesfg,estimates_mf];

    end 

    % Estimates of alpha,lamda and sigma

    averagest2=mean(estimatesfg');

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

    %% Non linear estimation - Gradient Descent 

    % Initializing 

    theta=rand(2*K2+1,1); 
   
    cost=[]
    step=0.00000000000001;
    i=0
    
    for ii=1:100
    % Gradient 

    
    %% IMPORTANT: THIS HAS TO BE ADJUSTED MANUALLY ACCORDING TO THE NUMBER OF VARIABLES 
    %% ADDED TO THE MODEL. (WE ESTIMATE: 2*K*(K+1)/2) 
    % We use a quadratic reparamterization of the diagonal elements. 
    % i.e. When K=7,we estimate 57 variance-covariance elements. The
    % elements that are reparamterized to the quadratic funtion are: 1,8,
    
    vu=[theta(1,1)^2;theta(2,1);theta(3,1);theta(4,1);theta(5,1);theta(6,1);theta(7,1);
        theta(8,1)^2;theta(9,1);theta(10,1);theta(11,1);theta(12,1);theta(13,1);
        theta(14,1)^2;theta(15,1);theta(16,1);theta(17,1);theta(18,1);
        theta(19,1)^2;theta(20,1);theta(21,1);theta(22,1);
        theta(23,1)^2;theta(24,1);theta(25,1);
        theta(26,1)^2;theta(27,1);
        theta(28,1)^2;
        theta(29,1)^2;theta(30,1);theta(31,1);theta(32,1);theta(33,1);theta(34,1);theta(35,1);
        theta(36,1)^2;theta(37,1);theta(38,1);theta(39,1);theta(40,1);theta(41,1);
        theta(42,1)^2;theta(43,1);theta(44,1);theta(45,1);theta(46,1);
        theta(47,1)^2;theta(48,1);theta(49,1);theta(50,1);
        theta(51,1)^2;theta(52,1);theta(53,1);
        theta(54,1)^2;theta(55,1);
        theta(56,1)^2;
        theta(57,1)^2];

    grad_vu=diag([2*theta(1,1);theta(2,1)/abs(theta(2,1));theta(3,1)/abs(theta(3,1));theta(4,1)/abs(theta(4,1));theta(5,1)/abs(theta(5,1));theta(6,1)/abs(theta(6,1));theta(7,1)/abs(theta(7,1));
                  2*theta(8,1);theta(9,1)/abs(theta(9,1));theta(10,1)/abs(theta(10,1));theta(11,1)/abs(theta(11,1));theta(12,1)/abs(theta(12,1));theta(13,1)/abs(theta(13,1));
                  2*theta(14,1);theta(15,1)/abs(theta(15,1));theta(16,1)/abs(theta(16,1));theta(17,1)/abs(theta(17,1));theta(18,1)/abs(theta(18,1));
                  2*theta(19,1);theta(20,1)/abs(theta(20,1));theta(21,1)/abs(theta(21,1));theta(22,1)/abs(theta(22,1));
                  2*theta(23,1);theta(24,1)/abs(theta(24,1));theta(25,1)/abs(theta(25,1));
                  2*theta(26,1);theta(27,1)/abs(theta(27,1));
                  2*theta(28,1);
                  2*theta(29,1);theta(30,1)/abs(theta(30,1));theta(31,1)/abs(theta(31,1));theta(32,1)/abs(theta(32,1));theta(33,1)/abs(theta(33,1));theta(34,1)/abs(theta(34,1));theta(35,1)/abs(theta(35,1));
                  2*theta(36,1);theta(37,1)/abs(theta(37,1));theta(38,1)/abs(theta(38,1));theta(39,1)/abs(theta(40,1));theta(41,1)/abs(theta(41,1));theta(42,1)/abs(theta(42,1));
                  2*theta(42,1);theta(43,1)/abs(theta(43,1));theta(44,1)/abs(theta(44,1));theta(45,1)/abs(theta(45,1));theta(46,1)/abs(theta(46,1));
                  2*theta(47,1);theta(48,1)/abs(theta(48,1));theta(49,1)/abs(theta(49,1));theta(50,1)/abs(theta(50,1));
                  2*theta(51,1);theta(52,1)/abs(theta(52,1));theta(53,1)/abs(theta(53,1));
                  2*theta(54,1);theta(55,1)/abs(theta(55,1));
                  2*theta(56,1);
                  2*theta(57,1);]);    

    gradient=(grad_vu*Z'*(R1-Z*vu))*(-1);

    cost1=R1-Z*vu;

    cost1=cost1'*cost1;

    % Gradient descent 

    theta1=theta-step*gradient;
    vu1=[theta1(1,1)^2;theta1(2,1);theta1(3,1);theta1(4,1);theta1(5,1);theta1(6,1);theta1(7,1);
        theta1(8,1)^2;theta1(9,1);theta1(10,1);theta1(11,1);theta1(12,1);theta1(13,1);
        theta1(14,1)^2;theta1(15,1);theta1(16,1);theta1(17,1);theta1(18,1);
        theta1(19,1)^2;theta1(20,1);theta1(21,1);theta1(22,1);
        theta1(23,1)^2;theta1(24,1);theta1(25,1);
        theta1(26,1)^2;theta1(27,1);
        theta1(28,1)^2;
        theta1(29,1)^2;theta1(30,1);theta1(31,1);theta1(32,1);theta1(33,1);theta1(34,1);theta1(35,1);
        theta1(36,1)^2;theta1(37,1);theta1(38,1);theta1(39,1);theta1(40,1);theta1(41,1);
        theta1(42,1)^2;theta1(43,1);theta1(44,1);theta1(45,1);theta1(46,1);
        theta1(47,1)^2;theta1(48,1);theta1(49,1);theta1(50,1);
        theta1(51,1)^2;theta1(52,1);theta1(53,1);
        theta1(54,1)^2;theta1(55,1);
        theta1(56,1)^2;
        theta1(57,1)^2];

    cost2=R1-Z*vu1;
    cost2=cost2'*cost2;

    if cost2<cost1
    theta=theta1;
    else 
    theta=theta;
    end 
    % Cost function evaluation 
    vu3=[theta(1,1)^2;theta(2,1);theta(3,1);theta(4,1);theta(5,1);theta(6,1);theta(7,1);
        theta(8,1)^2;theta(9,1);theta(10,1);theta(11,1);theta(12,1);theta(13,1);
        theta(14,1)^2;theta(15,1);theta(16,1);theta(17,1);theta(18,1);
        theta(19,1)^2;theta(20,1);theta(21,1);theta(22,1);
        theta(23,1)^2;theta(24,1);theta(25,1);
        theta(26,1)^2;theta(27,1);
        theta(28,1)^2;
        theta(29,1)^2;theta(30,1);theta(31,1);theta(32,1);theta(33,1);theta(34,1);theta(35,1);
        theta(36,1)^2;theta(37,1);theta(38,1);theta(39,1);theta(40,1);theta(41,1);
        theta(42,1)^2;theta(43,1);theta(44,1);theta(45,1);theta(46,1);
        theta(47,1)^2;theta(48,1);theta(49,1);theta(50,1);
        theta(51,1)^2;theta(52,1);theta(53,1);
        theta(54,1)^2;theta(55,1);
        theta(56,1)^2;
        theta(57,1)^2];

    cost3=(R1-Z*vu3)'*(R1-Z*vu3);

    cost=[cost; cost3];
    
   
    
    ii=ii+1
    end
    figure(3)
    
    plot(cost)
    estimates=theta;
    estimates_mf=[estimates(1,1)^2;estimates(2,1);estimates(3,1);estimates(4,1);estimates(5,1);estimates(6,1);estimates(7,1);
                  estimates(8,1)^2;estimates(9,1);estimates(10,1);estimates(11,1);estimates(12,1);estimates(13,1);
                  estimates(14,1)^2;estimates(15,1);estimates(16,1);estimates(17,1);estimates(18,1);
                  estimates(19,1)^2;estimates(20,1);estimates(21,1);estimates(22,1);
                  estimates(23,1)^2;estimates(24,1);estimates(25,1);
                  estimates(26,1)^2;estimates(27,1);
                  estimates(28,1)^2;
                  estimates(29,1)^2;estimates(30,1);estimates(31,1);estimates(32,1);estimates(33,1);estimates(34,1);estimates(35,1);
                  estimates(36,1)^2;estimates(37,1);estimates(38,1);estimates(39,1);estimates(40,1);estimates(41,1);
                  estimates(42,1)^2;estimates(43,1);estimates(44,1);estimates(45,1);estimates(46,1);
                  estimates(47,1)^2;estimates(48,1);estimates(49,1);estimates(50,1);
                  estimates(51,1)^2;estimates(52,1);estimates(53,1);
                  estimates(54,1)^2;estimates(55,1);
                  estimates(56,1)^2;
                  estimates(57,1)^2];
   
   
    
    % Estimates fixing lambda
    
    estimatesfl=[estimatesfl,estimates_mf];

end 

%% Average estimates of variance-covariance elements of alpha, gamma, lambda and epsilon

averagest3=mean(estimatesfl');

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

coeff_withinNL=(X'*Omega_theta_inv*X)\(X'*Omega_theta_inv*Y) ; % Beta FGLS 

varcov_withinNL=[alpha_h;gamma_h;lambda_h;sigma_h];

% 2.1.3 Inference 

var_cov_coeff=inv(X'*Omega_theta_inv*X);

% To obtain t_statistics for significance test you need to divide the
% estimated coefficient by the s.e of the coefficient. This s.e. is the
% squared root of the variance of the coefficient that you can find in the
% diagonal of the variance-covariance of the estimated coefficients.
% (var_cov_coeff)
% e.g.
t_B2=coeff_withinNL(3,1)/(sqrt(abs(var_cov_coeff(3,3))))

save('Estimates_FGLS_WithinNL.mat','coeff_withinNL','varcov_withinNL');
