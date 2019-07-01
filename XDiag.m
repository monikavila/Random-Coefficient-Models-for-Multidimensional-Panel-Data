function [XDiag]=XDiag(N1,N2,T,K,X)
x=ones(1,K);
X_0=[x;X];
% To Create the diagonal matrix of X's 
    XX=ones(N1*N2*T,(N1*N2*T+1)*K); % Matrix to 
    % store the new matrix that is diagonal 
    
for i=2:(N1*N2*T+1); % We start from second row, 
    % so we don't consider the auxiliar row that 
    % we created 
    l=(i-1)*K+1;     % Get initial position 
    u=i*K;           % Get last position 
    zero=0*ones(1,(N1*N2*T+1)*K); % Storage Matrix, 
    % we don't put outside the loop because we want to overwrite each time.  
    zero([l:u])=X_0(i,:); % We store the row extracted 
    % in the place that we want in the zero vector 
    XX(i,:)=zero; % We replace the row i of 
    % the XX matrix with the new zero vector 
end 

XDiag=XX(2:end,(K+1):end);%Diagonal matrix 

end