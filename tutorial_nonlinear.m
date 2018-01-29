% ARMS Lab 2018

% this is tutorial of LMPC controller for nonlinear system. In this tutorial
% basic parts of Linearized MPC controller would be examined. This script
% shows how to implement controller for nonlinear system provided by
% equation dx/dt = f(x,u) and y=C*x+D*u. For details of derivation please
% refer to 
% Zhakatayev, Altay, et al. "Successive linearization based model 
% predictive control of variable stiffness actuated robots." IEEE AIM 2017
% 
% linearize, discretize and how to find solution
% with simple constraints are shown. We setup out tutorial on 
% simple mathematical pendulum, with equation of motion:
%          dx1/dt = x2
%          dx2/dt = (-g/l)*sin(x1) -b*x2 + u
%
% Model Predictive Control scheme can be formulated as an optimization 
% problem with horizon length np:
%
%   Minimize J = 0.5*(X-rr)'*Q*(X-rr) + 0.5u'*R*u = 0.5*u'G*u + f*u 
%       s.t. Acon*u <= Bcon
%            x(k+1) = f(x(k),u(k))
%            X = [x(1); x(2);...;x(np)];
%            
% In order to solve this problem MATLAB built-in quadprog() function is used.
% Please refer to documentation of quadprog() function for details. 
% In fact, any nonlinear optimization problem solver can be used to come up
% with a solution. For example, qpOASES is suitable for real-time operation
% of robotic systems.
close all; clear; clc;
% global parameters associated with dynamic model of the system 
sys.g = 9.81;
sys.l = 0.1;
sys.b = 0.2;
% this is matrices to represent output as y=C*x+D*u
C = eye(2);  
D = [0;0];
%initial point of states
x = [0.001; 0];
%control input max and min values
ui=0;         % zero cotrol input
umax = 80;
umin =-20;
% IMPORTANT PARAMETERS
np = 40;       % horizon length 
nx = 2;        % number of states 
nu = 1;        % number of inputs
no = size(C,1);% number of outputs
Ts = 0.001;    % step size
Tfinal = 0.5;  % final time
wx = [1000 1]; % relative importance of states
wu = 0.00001;  % penalizing weights of control inputs

% generating simple step reference of the form [0.5 ; 0]  
%!!you should provide additional 'np' points 
ref =0.5*[ones(1,Tfinal/Ts +np);...
          zeros(1,Tfinal/Ts +np)];
% model is anonymous function that represents the equation which describes 
% system dynamic model and in the form of dx/dt =f(x,u).
% In order to make simulation of an arbitrary nonlinear system 
% you should _CREATE_YOUR_OWN_ nonlin_eq(x,u) function that calculates dx/dt 
model = @(x,u) nonlin_eq(x,u,sys); 

% FOR MPC CONTROLLER
%initialization of vectors for reference and otput results 
rr = zeros(np*nx,1);
y  = zeros(no,Tfinal/Ts);
uh = zeros(nu,Tfinal/Ts);
%constraints matrices for inputs in the whole horizon in the form of Acon*u <=Bcon
[Acon,Bcon] = simple_constraints(umax,umin,np,nu);
% weighting coefficient reflecting the relative importance of states and control inputs reduction 
% ! relative importance of states shown in vector [1000 1],where 1000 is
% asociated with state1 and 1 with state2 respectively
% !! In general, values of matrices Q and R (i.e. vectors wx, wu )should be tuned in order to
% achieve satisfactory performance. You can vary values relative to each other
% to observe the difference on the output and influence of such values. 
Q = diag(repmat(wx, 1, np)); 
R = diag(repmat(wu, 1, np));
% Setting the quadprog with 200 iterations at maximum
opts = optimoptions('quadprog', 'MaxIter', 200, 'Display','off');
% MAIN SIMULATION LOOP
for t=1:Tfinal/Ts
    
    rr =  ref_for_hor(rr,ref,t,np,nx);% reference vector for whole horizon  
    y(:,t) = C*x+D*ui;                % evaluating system output 
    [x, dx] = RK4(x,ui,Ts,model);     % i.e. simulate one step forward with Runge-Kutta 4 order integrator
    % linearize_model() function _SHOULD_BE_MODIFIED_ according to your model
    [A, B, K] = linearize_model(x,dx,ui,sys);                           % linearization step
    [Ad,Bd,Kd] = discretize(A,B,K,Ts);                              % discretization step
    [G, f] = grad_n_hess(R, Q, Ad, Bd, C, D, Kd, rr, np, x);        % calculating Hessian and gradient of cost function
    u = quadprog(G, f, Acon, Bcon, [], [], [], [], [],opts);
    ui = u(1:nu);     %providing first solution as input to our system
    uh(:,t) = ui;     %storing input
end
% plotting the results
tt = Ts:Ts:Tfinal;
subplot(3,1,1)
plot(tt ,y(1,:),'b', tt ,ref(1,1:(Tfinal/Ts)),'r--');
title('x_1(t) vs t');
legend('state1','reference');
subplot(3,1,2)
plot(tt ,y(2,:),'b',tt,ref(2,1:(Tfinal/Ts)),'r--');
legend('state2','reference');
title('x_2(t) vs t');
subplot(3,1,3)
plot(tt ,uh);
title('u(t) vs t');

