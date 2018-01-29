%ARMS Lab 2018
%calc_hp.m

function [Hx,Px,Km] = calc_hp(A,B,C,D,np)
%This function calculates prediction matrices for vector x and output
%vector y with prediction horizon np
    
 %Initialization
    
    % number of states
    nx = size(A, 1);
    %number of inputes 
    nu = size(B, 2);
    %number of outputs 
    no=size(C,1);
    
    %zero initialization  
    Px=zeros(np*nx,size(A,2));
    Hx=zeros(np*size(B));
    Km=zeros(no*np,size(C,2));
    S=zeros(size(C));
    
    %start of the main loop 
    for ind1=1:np
        
        % Filling Matrices Px,P,Km recursively 
        Px((1+(ind1-1)*nx):ind1*nx,1:nx)=A^ind1;
        Km(1+no*(ind1-1):no*ind1,:)=S;
        S=S+C*A^(ind1-1);
        
        %Filling Marices Hx, H recurcively  
        for ind2=1:np            
            if(ind1>=ind2)                
                Hx((1+(ind1-1)*nx):ind1*nx,(1+(ind2-1)*nu):ind2*nu)=A^(ind1-ind2)*B;                
            end
        end
    %End of the main loop    
    end
    
%End of function     
end