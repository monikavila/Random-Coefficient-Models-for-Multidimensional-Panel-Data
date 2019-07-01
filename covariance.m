function [Cov]=covariance(theta_sigma,theta_alpha,theta_gamma,theta_lambda,...
    K,D,X1,X2,X3,I_i,I_j,I_T,I_N1N2T);
 
au0_1=(0:1:K-1);
  au0_2=(1:1:K);
  
 aux=[1,1];
 aux=[aux,0:K-1];
 
theta_alpha2=D*theta_alpha;
theta_gamma2=D*theta_gamma;
theta_lambda2=D*theta_lambda;


for i=1:K
    
   sigma_alpham(:,i)=theta_alpha2([au0_1(i)*K+1:au0_2(i)*K],1);
   sigma_gammam(:,i)=theta_gamma2([au0_1(i)*K+1:au0_2(i)*K],1);
   sigma_lambdam(:,i)=theta_lambda2([au0_1(i)*K+1:au0_2(i)*K],1);
end

%sigma_alpham=reshape(theta_alpha2,[K,K]);
%sigma_gammam=reshape(theta_gamma2,[K,K]);
%sigma_lambdam=reshape(theta_lambda2,[K,K]);


%Omega with the estimates 
Cov=X1*kron(I_i,sigma_alpham)*X1'+X2*kron(I_j,sigma_gammam)*X2'+...
      X3*kron(I_T,sigma_lambdam)*X3'+theta_sigma*I_N1N2T;
end 