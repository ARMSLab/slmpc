%ARMS Lab 2018
%discretize.m
function [A,B,K] = discretize(Ac,Bc,Kc,Ts)
    %Applies discretization on continious system of the form
    %dx/dt = Ac*x +Bc*u+Kc and calculates such A,B,K that provides
    %analogical discrete system of form x(k+1)=A*x(k)+B*u(k)+K
    %more info in https://en.wikipedia.org/wiki/Discretization 
    A=expm(Ac*Ts);                      
    K=Ac\(A-eye(size(A)))*Kc;              
    B=Ac\(A-eye(size(A)))*Bc;
end