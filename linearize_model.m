%ARMS Lab 2018
%linearize model.m
function [Ac,Bc,Kc] = linearize_model(x,dx,u,sys)
    %this function is linearization of nonlinear model given by form
    %dx/dt=f(x,u). Approximation is in the form dx/dt =Ac*x+Bc*u+Kc
    %You can find more info about linearization in
    %https://en.wikipedia.org/wiki/Linearization
    %Basically Matrices Ac, Bc are Jacobians of f(x,u) about x and u
    %respectively
    
    %global parameters used to represent constant parameters in dynamical
    %model of the system. 
    g =sys.g;
    l=sys.l;
    b=sys.b;
    %This two matrices, as it was mentioned, Jacobians of f(x,u) about its
    %parameters
    Bc = [0;1];
    Ac = [0 1; (-g/l)*cos(x(1)) -b];
    
    %Kc is offcet created by linear approximation of the system and initial
    %parameters. 
    Kc=dx-Ac*x-Bc*u;
end
