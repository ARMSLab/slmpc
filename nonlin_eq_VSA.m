%ARMS Lab 2018 
%nonline_sim_VSA.m 
function dy = nonlin_eq_VSA(x,u,sys)
% This function calculates next position, speed in small time Ts 
% with initial position x and input u,which is considered to be constant in 
% time Ts. In order to find position RK4 integration method is used 
        
        %this part of the code represents model of the system and must be
        %changed
        
        %BEGINNING

    Bk = sys.Bk;
    u1i = sys.u1i;
    u2i = sys.u2i;
    x0 =sys.x0;
    N1=sys.N1;
    Rl=sys.R1;
    PI=sys.PI;
    beta1 = sys.beta1;
    beta2 =sys.beta2;
    alpha1 = sys.alpha1;
    alpha2 = sys.alpha2;
    Rp = sys.Rp;
    Im = sys.Im;
    bm = sys.bm;
    invM = sys.invM;
    Kk = sys.Kk;
        del1 = x0-Rl*(x(1)+PI/2)+Rp*(x(4)-u1i);
        del2 = x0+Rl*(x(1)+PI/2)-Rp*(x(5)-u2i);
        
        T1 = alpha1*del1*del1 + beta1*del1; %// tendon tension !
        T2 = alpha2*del2*del2 + beta2*del2;
        
        a = ( (T1 - T2)*Rl - Bk(1,1)*x(2)- cos(x(1))*N1);
        b = (u(3)- Bk(2,2)*x(3));
        
        dy01 = x(2);
        dy02 = invM(1,1)*a + invM(1,2)*b;
        dy03 = invM(2,1)*a + invM(2,2)*b; 
        dy04 = x(6) ;
        dy05 = x(7);
        dy06 = (Kk*(u(1)) - (T1)*Rp - bm*x(6))/Im;
        dy07 = (Kk*(u(2)) + (T2)*Rp - bm*x(7))/Im;
        
        dy = [dy01;dy02;dy03;dy04;dy05;dy06;dy07];
        %END
        

end
