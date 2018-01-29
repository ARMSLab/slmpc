function [Acon,Bcon] = simple_constraints(umax,umin,np,nu)
    Ac = [1;-1];
    Acon = zeros(2*nu*np,nu*np);
    for ind1=1:np
        Acon((2*(ind1-1)+1):2*ind1,(1*(ind1-1)+1):ind1)=Ac;
    end
    Bcon = repmat([umax; -umin],np,1);
end