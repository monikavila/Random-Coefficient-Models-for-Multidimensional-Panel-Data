clear all 
close all

% CODE: Prediction of Coefficients 
% DATE: 14/06/17
% GOAL: Prediction of Coefficients 
% AUTH: Monika Avila M?rquez

% Name of your path

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


% For WINDOWS

% cd('C:\Users\AVILAMAR\Dropbox\1. PHD\1. RESEARCH\Journal of Econometrics Proposal\3. PERU_Cognitivedevelopment\MATLAB\')

%% Upload Data 

% For mac 

data=load('Peru_Model1.mat');

%% Initial Settings 

N1=2;
N2=155;
T=2;
K=7;

% Je besoin: X, Y, r
X_hat=data.X_hat;
X=data.X;
X1=data.X1;
X2=data.X2;
X3=data.X3;
Y=data.Y;
r=data.r; 
I_N1=data.I_N1;
I_N2=data.I_N2;
I_T=data.I_T;
I_N1N2T=eye(N1*N2*T);
iota_N1=data.iota_N1;
iota_N2=data.iota_N2;
iota_T=data.iota_T;
iota_N1N2=ones(N1*N2,1);
iota_N1T=ones(N1*T,1);
iota_N2T=ones(N2*T,1);
iota_N1N2T=ones(N1*N2*T,1);
I_N1N2=eye(N1*N2);
I_N1T=eye(N1*T);
I_N2T=eye(N2*T);
I_N1N2T=eye(N1*N2*T);
I_K=eye(K);
D=data.D; 

%clearvars data  

% Getting the Variance-Covariance Matrices 

load('Estimates_FGLS_WithinNL');

theta_alpha2=varcov_withinNL(1:28,:);
theta_gamma2=varcov_withinNL(29:56,:);
theta_lambda2=varcov_withinNL(57:84,:);
theta_sigma2=varcov_withinNL(85:85,:);
aux=[1,8,14,19,23,26,28];
size_aux=size(aux);

for i=1:size_aux(1,2) 
theta_alpha2(aux(1,i),1)=(theta_alpha2(aux(1,i),1))^2;
theta_gamma2(i,1)=(theta_gamma2(aux(1,i),1))^2;
theta_lambda2(i,1)=(theta_lambda2(aux(1,i),1))^2;
end 
theta_alpha3=D*theta_alpha2;
theta_gamma3=D*theta_gamma2;
theta_lambda3=D*theta_lambda2;

au0_1=(0:1:K-1);
au0_2=(1:1:K);

for i=1:K
    
    sigma_alpham(:,i)=theta_alpha3([au0_1(i)*K+1:au0_2(i)*K],1);
    sigma_gammam(:,i)=theta_gamma3([au0_1(i)*K+1:au0_2(i)*K],1);
    sigma_lambdam(:,i)=theta_lambda3([au0_1(i)*K+1:au0_2(i)*K],1);
end

a=1

% Sigma 

Sigma=superkron(I_N1,iota_N2T,I_K)*kron(I_N1,sigma_alpham)*superkron(I_N1,iota_N2T',I_K)+superkron(iota_N1,I_N2,iota_T,I_K)*kron(I_N2,sigma_gammam)*superkron(iota_N1',I_N2,iota_T',I_K)...
    +superkron(iota_N1N2,I_T,I_K)*kron(I_T,sigma_lambdam)*superkron(iota_N1N2',I_T,I_K); 
a=2

%% Inverse of Sigma: Laszlo 

inv_sigma_alpham=pinv(sigma_alpham);
inv_A=superkron(I_N1,iota_N2T,I_K)*kron(I_N1,inv_sigma_alpham)*superkron(I_N1,iota_N2T',I_K); 

inv_sigma_gammam=pinv(sigma_gammam);
inv_B2=kron(I_N2,inv_sigma_gammam); 
B1=superkron(iota_N1,I_N2,iota_T,I_K);

inv_sigma_lambdam=pinv(sigma_lambdam);
inv_C2=kron(I_T,inv_sigma_lambdam); 
C1=superkron(iota_N1N2,I_T,I_K);

inv_Q=inv_A-inv_A*B1*(pinv(inv_B2+B1'*inv_A*B1))*B1'*inv_A;

Sigma_inv=inv_Q-inv_Q*C1*(pinv(inv_C2+C1'*inv_Q*C1))*C1'*inv_Q;

a=3
%% Prediction of coefficients 

beta_predicted=pinv(X_hat'*X_hat*(theta_sigma2)^(-2)+Sigma_inv)*((X_hat'*Y*(theta_sigma2)^(-2)+Sigma_inv*kron(iota_N1N2T,coeff_withinNL)))
a=4

save('beta_pr.mat','beta_predicted');
