% CODE: Function to obtain X_hat=diag(x_ijt')
%
% GOAL: Obtain X_hat=diag(x_ijt')
%
% AUTH: Monika Avila M?rquez
%
% Date: 28-06-16
%-------------------------------------------------

function [HDiag]=HDiag(N1,N2,T,K,H)
    
    HH=ones(N1*N2*T*(3*K*K+1),N1*N2*T*(3*K*K+1)*N1*N2*T); % Intermediate storage Matrix
    H_O=[ones(1,N1*N2*T*(3*K*K+1));H]                     % Storage matrix  
    
    for i=2:(N1*N2*T*(3*K+1)+1); % Start iteration from second row
        
        l=(i-1)*K+1;                      % Get initial position 
        u=i*K;                            % Get last position 
        zero=0*ones(1,N1*N2*T*(3*K*K+1)); % Storage Matrix, 
                                          % Don't put outside the loop 
                                          % because we want to overwrite each time.  
        zero([l:u])=H_0(i,:);             % Store the row extracted 
                                          % in the place that we want in the zero vector 
        KK(i,:)=zero;                     % Replace the row i of 
                                          % the XX matrix with the new zero vector 
    end 

    HDiag=HH(2:end,(N1*N2*T+1):end);      % Diagonal matrix 

end