clear all 
close all

% CODE: Linear Random Coefficient Model Estimation
% GOAL: Data arrangenment  
% AUTHOR: Monika Avila M?rquez
% DATE: 10/03/18
% If you use this matlab code for estimation of a linear RCM for Three 
% Dimensional Panel Data please paste the following disclaimer: 
% 
%% Set path where you have the functions and codes for model estimation
% For Thinkpad 

%'C:\Users\Monika\Dropbox\1.PHD\1.RESEARCH\2.PROJECTS\1. PARAMETRIC\9.MATLABCODE\OfficialCodestoRelease'
% Path containing codes for model estimation
% addpath('C:\Users\Monika\Dropbox\1.PHD\1.RESEARCH\2.PROJECTS\1.PARAMETRIC\9.MATLABCODE\OfficialCodestoRelease')  
% savepath pathdef.m
% Path containing Auxiliary functions needed for the codes
% addpath('C:\Users\Monika\Dropbox\1.PHD\1.RESEARCH\2.PROJECTS\1.PARAMETRIC\9.MATLABCODE\OfficialCodestoRelease\AuxiliaryFunctions')  
% savepath pathdef.m
% Path containing Data
% addpath('C:\Users\Monika\Dropbox\1.PHD\1.RESEARCH\2.PROJECTS\1.PARAMETRIC\9.MATLABCODE\OfficialCodestoRelease\Data')  
% savepath pathdef.m


%For compueter UNI

% Path containing codes for model estimation
 addpath('C:\Users\avilamar\Dropbox\1.PHD\1.RESEARCH\2.PROJECTS\1.PARAMETRIC\9.MATLABCODE\OfficialCodestoRelease\ModelEstimation')
 savepath pathdef.m
% Path containing Auxiliary functions needed for the codes
 addpath('C:\Users\avilamar\Dropbox\1.PHD\1.RESEARCH\2.PROJECTS\1.PARAMETRIC\9.MATLABCODE\OfficialCodestoRelease\AuxiliaryFunctions')
 savepath pathdef.m
% Path containing Data
 addpath('C:\Users\avilamar\Dropbox\1.PHD\1.RESEARCH\2.PROJECTS\1.PARAMETRIC\9.MATLABCODE\OfficialCodestoRelease\Data')
 savepath pathdef.m

%% Upload Data 

% Set the directory where you have your data 

% For mac users: 

% Name of your path

%path='/Users/moka/Dropbox/1.PHD/1.RESEARCH/2.PROJECTS/PROJECT1PARAMETRIC/9.MATLABCODE/'

%cd(path)

% Upload data 

filename='PeruModel1.csv';
data=csvread(filename);
data=data(1:620,:);

%% Initial conditions of the model 

    K=7;      % Number of variables + the intercep 
    K2=K*(K+1)/2;  % Number of free elements in each variance-covariance matrix
    N1=2;         % Number of individuals in the i(alpha) dimension
    N2=155;        % Number of individuals in the j(gamma) dimension
    T=2;           % Number of individuals in the t(lambda) dimension
    re=3;          % Number of Random Elements 
%% Auxiliar Matrices 

    I_N1=sparse((eye(N1)));     % Identity size i level 
    I_N2=sparse((eye(N2)));     % Identity size j level 
    I_T=sparse((eye(T)));       % Identity size t level 
    I_K=sparse((eye(K)));       % Identity size K level 
    I_N1N2T=sparse(eye(N1*N2*T)); 
    iota_N1=ones(N1,1); % Iota size i level
    iota_N2=ones(N2,1); % Iota size j level
    iota_T=ones(T,1);   % Iota size ti level
    iota_K2=ones(K2,1); % Iota size one off-diagonal plus the diagonal level
    iota_K=ones(K,1);   % Iota size K
    iota_re=ones(re,1); % Iota size k 

    [D,L]=(dupel(K));     % Generation of Duplication and Elimination matrices


%% X matrix 

X=data(:,7:13);

%% Y vector

Y=data(:,6);

%% Auxiliar matrices for Construction of Hau Matrices 

    % When k=k' : H_kk
    % When k=k' : H_kk'+ H_k'k

    au1_1=(0:1:2);              % Goes from 0 to 2, which is the number of random
                                %  components in the model minus 1

    au1_A=kron(au1_1',iota_K2); % We create a matrix that repeats 

    % Matrix of positions 

     % Auxiliars for obtaining matrix with columns kk', for the cases when they
     % are not repeating. 

     au2_1=(0:1:K-1);
     au2_2=(1:1:K);
     au2_3=(K:-1:1);
     au2_4=[0,au2_3(1,[1:K-1])];

    % With this loop we obtain the column of k and k' to have all the possible 
    % non repeating combinations of kk'.

    for j=1:K

      % For column k  
      au2_5([sum(au2_4(1,[1:au2_2(j)]))+1:...
          sum(au2_3(1,[1:au2_2(j)]))],1)=...
           au2_2(j)*ones(au2_3(j),1);

       % For column k'  
      au2_6([sum(au2_4(1,[1:au2_2(j)]))+1:...
          sum(au2_3(1,[1:au2_2(j)]))],1)=...
          au2_2(1,[j:K]);   
    end 

    au2_7=[au2_5,au2_6]; % this is the matrix with the non repeating combinations of kk'
    % The first column is k, and the secon is k'. 
    au2_8=kron(iota_re,au2_7);%The matrix of the combinations of kk' has to be repated re,
    % which corresponds to the random elements in the coefficient
    au2_9_1=ones(K2,2);% Matrix of ones to substract for the final matrix of positions 
    au2_9=kron(iota_re,au2_9_1);;%the columns are two because it corresponds to kk'
    au2_A=au2_8-au2_9; %Final matrix of positions 
    % kk'

    % Auxiliar for initial position
    % Last line 

    au3_1=[re, 0, 0];    % This is the last row and corresponds to epsilon
    au3_2=[au1_A,au2_A]; 

    au3=[au3_2;au3_1];

    %La primera columna multiplica el termino N1N2TKK
    %La segunda columna multiplica el t?rmino N1N2TK
    %La tercera columna multiplica el t?rmino N1N2T

    % Proofs of auxiliar matrices 
    au3_3=N1*N2*T*K*K*au3(:,1);
    au3_4=N1*N2*T*K*au3(:,2);
    au3_5=N1*N2*T*au3(:,3);
    au3_6=[au3_3,au3_4,au3_5];
    au3_7=sum(au3_6,2) + ones(3*K2+1,1);

    % Auxiliar for last position

    % Last line 

    au4_1=[3, 0, 1];
    au4_2=[au1_A,au2_A(:,1),au2_8(:,2)]; 
    au4=[au4_2;au4_1];

    %La primera columna multiplica el termino N1N2TKK
    %La segunda columna multiplica el t?rmino N1N2TK
    %La tercera columna multiplica el t?rmino N1N2T

    % Proofs of auxiliar matrices 
    au4_3=N1*N2*T*K*K*au4(:,1);
    au4_4=N1*N2*T*K*au4(:,2);
    au4_5=N1*N2*T*au4(:,3);
    au4_6=[au4_3,au4_4,au4_5];
    au4_7=sum(au4_6,2);

    % Final Matrix with Initial and Last positions, used in the loop. The
    % initial position is in column one of the matrix and the last position is 
    % in the second column 

    au_5=[au3_7,au4_7];

    % Proof of columns between initial and last position, to check that we're
    % going to take N1*N2*T 

    au_6=au4_7-au3_7;

    % Auxiliar for if condition. For k=k'.

        au10_1=[1,1];

        % Column 1 is k and Column 2 is k' 
        au10_A=[au2_8;au10_1];

    au0_1=(0:1:K-1);
    au0_2=(1:1:K);
    
    % Duplication matrix for vech(RHS)
    
    
    au11=N1*N2*T*0:N1*N2*T*1:N1*N2*T*(3*K*(K+1)/2+1-1);
    au12=N1*N2*T*1:N1*N2*T*1:N1*N2*T*(3*K*(K+1)/2+1);
    
    auH12=eye(N1*N2*T);
 
    row=ones(K2*3+1,1);

    rhs=0*ones(3*K*(K+1)/2+1,1);
    LHS1=0*ones((3*K*(K+1)/2+1)^2,1);

    sigma_alphamA=eye(K);
    sigma_gammamA=eye(K);
    sigma_lambdamA=eye(K);
 
   % [D2,L2]=dupel(3*K*(K+1)/2+1);
    
    % For the rows that we have to replace in matrix RHS after each i in
    % the loop
    
%% Auxiliar matrices for Construction of Anderson Algorithm and Hau

    au15=1:(3*K*(K+1)/2)+1:(3*K*(K+1)/2+1)*((3*K*(K+1)/2+1));
    au16=(3*K*(K+1)/2+1):(3*K*(K+1)/2+1):(3*K*(K+1)/2+1)*((3*K*(K+1)/2+1));


       % To Create the diagonal matrix of X's 

    [X_hat]=sparse(XDiag(N1,N2,T,K,X));

    %The three matrices for the model 

    X1=sparse(X_hat*kron(superkron(I_N1,iota_N2,iota_T),I_K)); % To create the design matrix for element alpha
    X2=sparse(X_hat*kron(superkron(iota_N1,I_N2,iota_T),I_K)); % To create the design matrix for element gamma
    X3=sparse(X_hat*kron(superkron(iota_N1,iota_N2,I_T),I_K));% To create the design matrix for element lambda

% The Projector and Emulator matrices 

    P=X*(inv(X'*X))*X';                      % Projector matrix 
    M=eye(N1*N2*T,N1*N2*T)-X*(inv(X'*X))*X'; % Residual Matrix 
    
    
    r=Y-P*Y;                                 % Residuals 

    R=r*r';                                  % Outer product of residuals 

% H matrices: PROBLEM... MATRIICES ARE TOO BIG. Solution is to cut them into two  
 
 % For setting the Variance Covariance Matrix as a linear
 % combination
 
 % For setting the Variance Covariance Matrix as a linear
 % combination. 

 i=1 % This is for the intercept 
[H_alpha_1,H_gamma_1,H_lambda_1]=HMatrices_BigData060317(i,I_N1,I_N2,I_T,iota_N1,...
                                                 iota_N2,iota_T,N1,N2,T,K,X);

i=2
[H_alpha_2,H_gamma_2,H_lambda_2]=HMatrices_BigData060317(i,I_N1,I_N2,I_T,iota_N1,...
                                                 iota_N2,iota_T,N1,N2,T,K,X);
                                             
                                             i=3
[H_alpha_3,H_gamma_3,H_lambda_3]=HMatrices_BigData060317(i,I_N1,I_N2,I_T,iota_N1,...
                                                 iota_N2,iota_T,N1,N2,T,K,X);
                                             
                                             i=4
[H_alpha_4,H_gamma_4,H_lambda_4]=HMatrices_BigData060317(i,I_N1,I_N2,I_T,iota_N1,...
                                                 iota_N2,iota_T,N1,N2,T,K,X);
                                             
                                             i=5
[H_alpha_5,H_gamma_5,H_lambda_5]=HMatrices_BigData060317(i,I_N1,I_N2,I_T,iota_N1,...
                                                 iota_N2,iota_T,N1,N2,T,K,X);
                                             
                                             i=6
[H_alpha_6,H_gamma_6,H_lambda_6]=HMatrices_BigData060317(i,I_N1,I_N2,I_T,iota_N1,...
                                                 iota_N2,iota_T,N1,N2,T,K,X);

i=7

[H_alpha_7,H_gamma_7,H_lambda_7]=HMatrices_BigData060317(i,I_N1,I_N2,I_T,iota_N1,...
                                                 iota_N2,iota_T,N1,N2,T,K,X);

                                            
save('Peru_Model1.mat','K','K2','N1','N2','T','re','I_N1',...
    'I_T','I_N1N2T','iota_N1','iota_N2','iota_T','I_N2','iota_K','iota_K2',...
    'iota_re','Y','X','X1','X2','X3','M','r', 'L','D','X_hat', 'H_alpha_1',...
    'H_alpha_2','H_alpha_3','H_alpha_4', 'H_alpha_5','H_alpha_6','H_alpha_7', ...
    'H_gamma_1', 'H_gamma_2', 'H_gamma_3', 'H_gamma_4', 'H_gamma_5', 'H_gamma_6','H_gamma_7',...
    'H_lambda_1', 'H_lambda_2', 'H_lambda_3', 'H_lambda_4','H_lambda_5', 'H_lambda_6','H_lambda_7','-v7.3')        

