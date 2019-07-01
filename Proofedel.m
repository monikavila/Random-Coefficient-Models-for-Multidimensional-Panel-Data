% CODE: Proof Duplication and Elimination Matrices 
%
% GOAL: Proof that matrices D and L are correct 
%
% AUTH: Monika Avila M?rquez
%
% BIB.:  Magnus Paper: "The Elimination Matrix
                     % Lemmas and Applications."

%
% Date: 28-06-16
%-------------------------------------------------

clear all;

close all; 

K=6;

[D,L]=dupel(K);

%% Proof 

C=20*eye(K)+5*ones(K,K);

vC=vec(C);

vh=L*vC;

vC2=D*vh;

C2=reshape(vC,K,K);

proof=[vC2,vC];