%ARMS Lab 2018
%grad_n_hess.m
function [G,f] = grad_n_hess(R , Q , A , B , C , D ,K, rr,np, x)
        %This function returns gradient and hessian matrices of the cost
        %function using function calc_hp to obtain horizon matrices.
        [Hx,Px,Km] = calc_hp(A,B,C,D,np);
        K1 = (Km*K)-rr;
        G = Hx'*Q*Hx + R;
        f = Hx'*Q*(Px*x+K1);
        
end