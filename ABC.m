%% CODE: Function for M(H_{re,kk})M
%
%  GOAL: Obtain M(H_{re,kk})M
%
%  AUTH: Monika Avila Marquez
%
%  Date: 28-06-16
%-----------------------------------------------------------------------

function [Af,Bf,Cf]=ABC(N1,N2,T,K,MH_alpha1M,MH_gamma1M,MH_lambda1M)

A1=ones(N1*N2*T*N1*N2*T,K*K); % Matrix to store all the vec(MH_kkM)
B1=ones(N1*N2*T*N1*N2*T,K*K); % Matrix to store all the vec(MH_kkM)
C1=ones(N1*N2*T*N1*N2*T,K*K); % Matrix to store all the vec(MH_kkM)

for j=2:K*K+1
    
    AA=vec(MH_alpha1M([1:end],[(j-1)*N1*N2*T+1:j*N1*N2*T])); % Matrix to store all the vec(MH_kkM)
    A1(:,j)=AA(:,:);
    
    BB=vec(MH_gamma1M([1:end],[(j-1)*N1*N2*T+1:j*N1*N2*T])); % Matrix to store all the vec(MH_kkM)
    B1(:,j)=BB(:,:);
    
    CC=vec(MH_lambda1M([1:end],[(j-1)*N1*N2*T+1:j*N1*N2*T])); % Matrix to store all the vec(MH_kkM)
    C1(:,j)=CC(:,:);
    
end

Af=A1(:,[2:end]); % Droping the first matrix 

Bf=B1(:,[2:end]); % Droping the first matrix 

Cf=C1(:,[2:end]); % Droping the first matrix 


end 

