% ARMS Lab 2018
% Tutorial_Linear.m

close all; clear; clc;
% this is tutorial of LMPC controller for linear system. In this tutorial
% basic parts of MPC controller would be shown.

% Matrices A,B,C,D provides dynamic system of the form  dx/dt =Ax +Bu , y=Cx+Du

A = -[0  -1; 2 1];
B = [0; 1];

rank([B A*B])

C = eye(2);
D = [0;0];
x = [0; 0];

% upper and lower constraints on input u
ui=0;
umax = 40;
umin =-20;

%horizon length of MPC controller
np=20;

%sampling time of the MPC controller 
Ts = 0.001;
%simulation star and end time 
Tinitial = 0;
Tfinal = 1;

%generating reference of the form [1;0] for whole simulation time
ref = 1*[ones(1,Tfinal/Ts +np);zeros(1,Tfinal/Ts +np)];

%model represents the equation that describes system dynamic midel and in
%form of dx/dt =f(x,u)
model = @(x,u) A*x + B*u;

%discretization of the linar system
[Ad,Bd,Kd] = discretize(A,B,[0;0],Ts);

% number of states, number of inputs, number of outputs 
nx = size(A,1);
nu = size(B,2);
no = size(C,1);

%PARAMETERS FOR MPC CONTROLLER

%initialization of vector storing all elements of reference values of
%states for whole horizon length 
rr=zeros(np*nx,1);

% initialization of vector to store output values in for each sampling time
y = zeros(no,Tfinal/Ts);
%vector to store input values derived from MPC controller
uh = zeros(nu,Tfinal/Ts);

%constructing constraints for inputs in the whole horizon in the 
%form of Acon*u <=Bcon
Ac = [1;-1];
Acon = zeros(2*nu*np,nu*np);
for ind1=1:np
    Acon((2*(ind1-1)+1):2*ind1,(1*(ind1-1)+1):ind1)=Ac;
end
Bcon = repmat([umax; -umin],np,1);

% weighting coefficient reflecting the relative importance of states for
% whole horizon 

% ! relative importanse of states shown in vector [1000 1],where 1000 is
% asosiated with state1 and 1 with state2 respectively
Q = diag(repmat([1000 1],1,np));
% weighting coefficient penalizing relative big changes in inputs for
%whole horizon
R = diag(repmat([0.0001],1,np));

%adjusting solver to make 200 iterations at maximum
opts = optimoptions('quadprog', 'MaxIter', 200);

% MAIN SIMULATION LOOP
for t=1:Tfinal/Ts
    
    %writing data from reference to vector rr in appropriate way    
    for ind2 = 1:np
        rr(nx*(ind2-1)+1:ind2*nx,1)=ref(:,t+ind2-1);
    end
    %witing current state to the y
    y(:,t)= C*x+D*ui; 
    %calculating Hessian and gradient of cost function  
    [G,f] = grad_n_hess(R , Q , Ad , Bd , C , D , Kd, rr,np, x);
    
    % using quadprog optimizer to find solution for our cost function in
    % given constraints, solution is given as u = [u(1); u(2); ... ; u(np)]
    u = quadprog(G,f,Acon, Bcon, [], [], [], [], [],opts);
    
    %providing first solution as input to our system
    ui = u(1);
    %storing input
    uh(:,t) = ui;
    %simulation of system by its model to one sampling time ahead with
    %Runge-Kutta 4 order integrator 
    [x, dx] = RK4(x,ui,Ts,model);
end

%plotting data with respect to reference 
tt = Ts:Ts:Tfinal;
subplot(3,1,1)
plot(tt ,y(1,:),'b', tt ,ref(1,1:(Tfinal/Ts)),'r--');
subplot(3,1,2)
plot(tt ,y(2,:),'b',tt,ref(2,1:(Tfinal/Ts)),'r--');
subplot(3,1,3)
plot(tt ,uh);




