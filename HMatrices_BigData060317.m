% CODE: Function to obtain H_re Matrices
%
% GOAL: Obtain H_re=X_hat,re,k*X_hat,re,k'
%
% AUTH: Monika Avila M?rquez
%
% Date: 28-06-16
%-------------------------------------------------

function [H_alphaf1,H_gammaf1,H_lambdaf1]=Hmatrices(i,I_N1,I_N2,I_T,iota_N1,...
                                                 iota_N2,iota_T,N1,N2,T,K,X)
%% to get the H Matrices 

%H_alpha 

x_1=sparse(0*ones(N1*N2*T,1));                      % Auxiliary column that allows the loop to start at 2 
X_1=[x_1 X];                              % Auxiliary matrix that allows the loop to start at 2 
H_alpha1K=sparse(0*ones(N1*N2*T,(floor(K/2)+1)*N1*N2*T));    % Matrix for storage the H matrices  
H_gamma1K=sparse(0*ones(N1*N2*T,(floor(K/2)+1)*N1*N2*T));    % Matrix for storage the H matrices  
H_lambda1K=sparse(0*ones(N1*N2*T,(floor(K/2)+1)*N1*N2*T));   % Matrix for storage the H matrices  

% We divide the matrices into two 

% In that case, we do first from 1 to K/2. After the rest. If K is not par,
% then we do from 1 to floor(K/2), and after from floor(K/2)+1 to K 

   
    a=diag(X_1(:,i)); % obtain the column that we want to diagonalize to have X_k
    
    for j=2:K+1
    
    l=(j-1)*N1*N2*T+1;% Initial column position where to store H_kk'
    u=(j)*N1*N2*T;    % Last column position where to store H_kk'
    
    b=diag(X_1(:,j)); % Obtain the column that we want to diagonalize to have (X_k)'
 
    H_alpha1K([1:end],[l:u])=(a*superkron(I_N1,iota_N2,iota_T))...
                            *(superkron(I_N1,iota_N2',iota_T')*b');% Get the matrix H-kk we want and Storing in the matrix for j loop
    H_alphaf1=H_alpha1K(:,[N1*N2*T+1:end]);   % Getting only the matrix that we want
    
    H_gamma1K([1:end],[l:u])=(a*superkron(iota_N1,I_N2,iota_T))...
        *(superkron(iota_N1',I_N2,iota_T')*b'); % Stores all generated H_gamma_ijk
    H_gammaf1=H_gamma1K(:,[N1*N2*T+1:end]);    % Gets rid of first N1*N2*T columns 
    
    H_lambda1K([1:end],[l:u])=(a*superkron(iota_N1,iota_N2,I_T))...
        *(superkron(iota_N1',iota_N2',I_T)*b');  % Stores H_lambda_ijk
    H_lambdaf1=H_lambda1K(:,[N1*N2*T+1:end]);   % Gets rid of first N1*N2*T columns 
   
    end
    
%     H_alpha1=ones(N1*N2*T,K+1)*N1*N2*T);    % Matrix for storage the H matrices  
% H_gamma1=ones(N1*N2*T,(floor(K/2)+1)*N1*N2*T);    % Matrix for storage the H matrices  
% H_lambdaf1=ones(N1*N2*T,(floor(K/2)+1)*N1*N2*T);   % Matrix for storage the H matrices  
% % I have to correct here, because the matrix is growing with the index....
% % !!!!!!!!!!!!!!!
%     for j=(floor(K/2)+1) +1 :K+1
%     
%     l=(j-1)*N1*N2*T+1;% Initial column position where to store H_kk'
%     u=(j)*N1*N2*T;    % Last column position where to store H_kk'
%     
%     b=diag(X_1(:,j)); % Obtain the column that we want to diagonalize to have (X_k)'
%  
%     H_alpha2K([1:end],[l:u])=(a*superkron(I_N1,iota_N2,iota_T))...
%                             *(superkron(I_N1,iota_N2',iota_T')*b');% Get the matrix H-kk we want and Storing in the matrix for j loop
%     H_alphaf2=H_alpha2K(:,[N1*N2*T+1:end]);   % Getting only the matrix that we want
%     
%     H_gamma2K([1:end],[l:u])=(a*superkron(iota_N1,I_N2,iota_T))...
%         *(superkron(iota_N1',I_N2,iota_T')*b'); % Stores all generated H_gamma_ijk
%     H_gammaf2=H_gamma2K(:,[N1*N2*T+1:end]);    % Gets rid of first N1*N2*T columns 
%     
%     H_lambda2K([1:end],[l:u])=(a*superkron(iota_N1,iota_N2,I_T))...
%         *(superkron(iota_N1',iota_N2',I_T)*b');  % Stores H_lambda_ijk
%     H_lambdaf2=H_lambda2K(:,[N1*N2*T+1:end]);   % Gets rid of first N1*N2*T columns 
%    
%     end
    

% We get rid of first K*N1*N2*T columns because they were only auxiliar 


end