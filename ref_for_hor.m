function rr = ref_for_hor(rr,ref,t,np,nx)
    for ind1 = 1:np
        rr(nx*(ind1-1)+1:ind1*nx,1)=ref(:,t+ind1-1);
    end
end