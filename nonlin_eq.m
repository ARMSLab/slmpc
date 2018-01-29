%ARMS Lab 2018
%nonlin_eq.m
function dx = nonlin_eq(x,u,sys)
    %This function represents nonlinear model of simple mathematical
    %pendulum in medium with linear viscocity 
    g = sys.g;
    l = sys.l;
    b= sys.b;
    dx =zeros(2,1);
    dx(1) = x(2);
    dx(2) = (-g/l)*sin(x(1)) -b*x(2)+ u;
end