% CODE: Function to obtain H_re Matrices
%
% GOAL: Obtain H_re=X_hat,re,k*X_hat,re,k'
%
% AUTH: Monika Avila M?rquez
%
% Date: 28-06-16
%-------------------------------------------------


function [H_alphaf,H_gammaf,H_lambdaf]=Hmatrices(I_i,I_j,I_T,iota_i,...
                                                 iota_j,iota_T,N1,N2,T,K,X)
%% to get the H Matrices 

%H_alpha 

x_1=ones(N1*N2*T,1);                      % Auxiliary column that allows the loop to start at 2 
X_1=[x_1 X];                              % Auxiliary matrix that allows the loop to start at 2 
H_alpha1k=ones(N1*N2*T,(K+1)*N1*N2*T);    % Matrix for storage the H matrices  
H_alpha1=ones(N1*N2*T,(K+1)*(K)*N1*N2*T); % Matrix with extra space for storage
H_gamma1k=ones(N1*N2*T,(K+1)*N1*N2*T);    % Matrix for storage the H matrices  
H_gamma1=ones(N1*N2*T,(K+1)*(K)*N1*N2*T); % Matrix with extra space for storage
H_lambda1k=ones(N1*N2*T,(K+1)*N1*N2*T);   % Matrix for storage the H matrices  
H_lambda1=ones(N1*N2*T,(K+1)*(K)*N1*N2*T);% Matrix with extra space for storage

for k=2:(K+1)%starts at 2 because we skip the first one that is not necessary
    
    a=diag(X_1(:,k));%obtain the column that we want to diagonalize
    
    for kk=2:(K+1) 
    
    l=(kk-1)*N1*N2*T+1;% Initial column position where to store
    u=(kk)*N1*N2*T;% last column position where to store
    
    b=diag(X_1(:,kk));%obtain the column that we want to diagonalize
 
    H_alpha1k([1:end],[l:u])=(a*superkron(I_i,iota_j,iota_T))...
                            *(superkron(I_i,iota_j',iota_T')*b');% Get the matrix H-kk we want and Storing in the matrix for j loop
    H_alpha11k=H_alpha1k(:,[N1*N2*T+1:end]); % Gettin only the matrix that we want
    
    H_gamma1k([1:end],[l:u])=(a*superkron(iota_i,I_j,iota_T))...
        *(superkron(iota_i',I_j',iota_T')*b');%Stores H_gamma_ijk
    H_gamma11k=H_gamma1k(:,[N1*N2*T+1:end]);%Stores all H_gamma_ijk
    
    H_lambda1k([1:end],[l:u])=(a*superkron(iota_i,iota_j,I_T))...
        *(superkron(iota_i',iota_j',I_T)*b');%Stores H_lambda_ijk
    H_lambda11k=H_lambda1k(:,[N1*N2*T+1:end]);%Stores all H_lambda_ijk
   
    end
    
    ll=(k-1)*(K)*N1*N2*T+1;% initial position where to store
    uu=(k)*(K)*N1*N2*T;% last position where to store
   
    H_alpha1([1:end],[ll:uu])=H_alpha11k(:,:);%Storing all data for the nested loop
    H_gamma1([1:end],[ll:uu])=H_gamma11k(:,:); % Storage of all H_gamma
    H_lambda1([1:end],[ll:uu])=H_lambda11k(:,:);% Storage of all H_lambda
end 

H_alphaf=H_alpha1(:,[N1*N2*T*K+1:end]); % Matrix with H.kk 
H_gammaf=H_gamma1(:,[N1*N2*T*K+1:end]); % Keeping only the H_gamma, we get rid of the first 
H_lambdaf=H_lambda1(:,[(N1*N2*T*K+1):end]);% Keeping only the H_lambda, we get rid of the first 

end