% CODE: Elimination and Duplication funtion
%
% GOAL: Obtain the Elimination and Duplication funtion    
%
% AUTH: Monika Avila M?rquez
%
% BIB: % Magnus Paper: The Elimination Matrix
        % Lemmas and Applications.
%        
% Date: 28-06-16
%-------------------------------------------------

function [Dup,El]=dupel(K)


K3=K*(K+1)/2;

% Uij and Uii 

    % Positions 

       auj=(1:1:K); % For the positions for u_ij (1 
                    % in a vector of zeros everywhere but in position j)

       w=[]; % To store the positions calculated in each loop

        for i=1:K
            
            for j=1:K
                
                if i>=j
                    
                z(j)=(auj(j)-1)*K+i-0.5*auj(j)*(auj(j)-1); 
                % Position where we have value 1 
                
                else 
                    
                z(j)=0; % Putting zero in the case that i<j
                
                end
                
            end
            w=[w,z]; % Storing the vector of positions created in each loop
            
        end

            %uij 
            
            UU=[];
            
            for i=1:K*K
                
              uu=0*ones(K3,1); %Vector with zeros, we create one in each loop           
                               % so we don't have the one of previous loop
              if w(1,i)~=0
                  
              uu(w(1,i))=1; % We put a one in the position obtained in previous loop 
              
              else 
                  
              uu=0*ones(K3,1); % We keep a vector of zeros is    
                                % the value is zero in the vector of positions
              end 
              
              UU=[UU,uu(:,:)]; % We store the vector created 
              
            end 

% Sum of ij elements 
 
    C2=0*ones(K3,K*K); % Create a matrix to store the sums of each loop

    I_K=eye(K); % Create an identity matrix to get the vectors 

    % e_i and e_j

    au1=(0:1:K);% Auxiliar to get the column positions 

    % from matrix wih vectors u_ij

    for i=1:K
        
        for j=1:K
            
            A2=superkron(UU(:,au1(j)*K+i),...
                I_K(:,j)',I_K(:,i)');       % u_ij(kron)e_j'(kron)e_i'
            
            C2(:,:)=C2(:,:)+A2(:,:);
            
        end
        
        El=0*ones(K3,K*K); % Storage matrix for the 
                           % final Elimination matrix. Inside the loop so we avoid 
                           % duplication of the sums. 
                           
        El(:,:)=C2(:,:)+El(:,:); % ELIMINATION MATRIX
        
    end 

% K

    K3=K*K;
    
    K4=0*ones(K3,K3); % To store the values 
    
    K5=0*ones(K3,K3); % To store the values 
    
    for i=1:K
        
        e1=I_K(:,i); % Vector e_i
        
        for j=1:K
            
            e2=I_K(:,j); % To get the vecto e_j
            
            Ee=e1*e2';
            
            eE=e2*e1';
            
            K1=kron(Ee,eE);
            
            K4(:,:)=K4(:,:)+K1(:,:);
            
        end
        
        K5(:,:)=K5(:,:)+K4(:,:);
        
    end

% N

        II=eye(K3);
        N=0.5*(II+K4);

% D 

        Dup=(N*El')*(inv((El*N*El')));

%% Proof 

        C=20*eye(K)+5*ones(K,K);

        vC=vec(C);

        vh=El*vC;

        vc=Dup*vh;

        proof=[vC,vc];
end 