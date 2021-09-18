function [ ke_tt, pe2_tt ] = KEKP( p2, m, nx, d, tau,eps,rmax )
% Tensor product for kinetic energy
e_tt = tt_ones(nx);
for k=1:d
    if(k == 1)
        kek_tt{k}=tt_tensor(p2/m(k));
        kep2_tt{k}=tt_tensor(exp(-1i*tau*p2/m(k)));
    else
        kek_tt{k}=e_tt;
        kep2_tt{k}=e_tt;
    end
    for j=2:d
        v_tt=e_tt;
        vp2_tt=e_tt;
        if (j == k)
            v_tt=tt_tensor(p2/m(j));
            vp2_tt=tt_tensor(exp(-1i*tau*p2/m(j)));
        end
        kek_tt{k}=tkron(kek_tt{k},v_tt);
        kep2_tt{k}=tkron(kep2_tt{k},vp2_tt);
    end
end
ke_tt=kek_tt{1};
pe2_tt=kep2_tt{1};
for k=2:d
    ke_tt=ke_tt+kek_tt{k};
    pe2_tt=times(pe2_tt,kep2_tt{k});
end
ke_tt=round(ke_tt,eps,rmax);
pe2_tt=round(pe2_tt,eps,rmax);
end

