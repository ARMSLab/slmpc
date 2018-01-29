%ARMS Lab 2018 
%nonline_sim_VSA.m 
function [x, dx] = RK4(x1,u,Ts,k)
% This function calculates next position, speed in small time Ts 
% with initial position x and input u,which is considered to be constant in 
% time Ts. In order to find position RK4 integration method is used 
        
        %Algorithm of RK4
        k1 = k(x1,u); 
        k2 = k(x1 + (Ts/2)*k1,u);         
        k3 = k(x1 + (Ts/2)*k2,u);        
        k4 = k(x1+Ts*k3,u);        
        
        dx=k1;
        
        x = x1 + Ts*(k1 + 2*k2 + 2*k3 + k4)/6;
end
