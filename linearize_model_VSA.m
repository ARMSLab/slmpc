function [Ac,Bc,Kc] = linearize_model_VSA(x,dx,u,sys)
    
    x0 =sys.x0;
    M =sys.M;
    N1=sys.N1;
    Rl=sys.R1;
    PI=sys.PI;
    beta1 = sys.beta1;
    beta2 =sys.beta2;
    alpha1 = sys.alpha1;
    alpha2 = sys.alpha2;
    u10 = sys.u10;
    u20 =sys.u20;
    Rp = sys.Rp;
    Im = sys.Im;
    bm = sys.bm;
    b1 = sys.b1;
    b2 = sys.b2;
    invM = sys.invM;
    Kk = sys.Kk;
    x1=x(1);
    x4=x(4);
    x5=x(5);
    
    r1 = M(1,1)*M(2,2) - M(2,1)*M(1,2);
    r7 = Rl*(x1 + PI/2);
    r6 = r7 - x0 + Rp*(u10 - x4 + (x0/Rp));
    r5 = r7 + x0 - Rp*(x5 - u20 + (x0/Rp));
    r4 = Rp*beta1 - Rp*alpha1*r6*2;
    r3 = Rp*beta2 + Rp*alpha1*r5*2;
    r2 = N1*sin(x1) - Rl*(Rl*beta1 + Rl*beta2 + Rl*alpha2*r5*2 - Rl*alpha2*r6*2);
    
    a1  = M(2,2)*r2/r1;
    a2  = -M(2,2)*b1/r1;
    a3  =  M(1,2)*b2/r1;
    a4  =  M(2,2)*Rl*r4/r1;
    a5  =  M(2,2)*Rl*r3/r1;
    a6  =  -M(2,1)*r2/r1;
    a7  =  M(2,1)*b1/r1;
    a8  = -M(1,1)*b2/r1;
    a9  = -M(2,1)*Rl*r4/r1;
    a10 = -M(2,1)*Rl*r3/r1;
    a11 = Rp*(Rl*beta1 - Rl*alpha1*r6*2)/Im;
    a12 = -Rp*r4/Im;
    a13 = -bm/Im;
    a14 = Rp*(Rl*beta2 + Rl*alpha2*r5*2)/Im;
    a15 = -Rp*r3/Im;
    a16 = -bm/Im;
    Ac=[0, 1, 0, 0, 0, 0, 0; a1, a2, a3, a4, a5, 0, 0 ; a6, a7, a8, a9, a10, 0, 0; 0, 0, 0, 0, 0, 1, 0; 0, 0, 0, 0, 0, 0, 1; a11, 0, 0, a12, 0, a13, 0; a14, 0, 0, 0, a15, 0, a16];
    
    Bc=zeros(7,3);
    Bc(2,3)=invM(1,2);
    Bc(3,3)=invM(2,2);
    Bc(6,1)=Kk/Im;
    Bc(7,2)=Kk/Im;
    
    Kc=dx-Ac*x-Bc*u;
    %discretization
end