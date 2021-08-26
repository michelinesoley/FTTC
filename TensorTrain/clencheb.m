function [G_tt] = clencheb(PE_tt,KE_tt,G_tt,Dp,Dm,dx,nx,dp,np,Npoly,eps,rmax,t)
% Chebyshev polynomials applied to state G_tt
d_tt{1}=round(0*G_tt,eps,rmax);
d_tt{2}=round(0*G_tt,eps,rmax);
for jj = 0:Npoly-1
    d_tt{3}=d_tt{2};
    d_tt{2}=d_tt{1};
    kk = Npoly-1-jj;
    d_tt{1}=round((-1i)^(kk)*besselj(kk,t)*G_tt,eps,rmax);
    d_tt{1}=round(d_tt{1}-d_tt{3},eps,rmax);
    temp_tt=round(times(PE_tt,d_tt{2}),eps,rmax);
    FT_temp_tt=round(times(KE_tt,tt_FT(d_tt{2},dx,nx,1)),eps,rmax);
    temp_tt = round((temp_tt+tt_FT(FT_temp_tt,dp,np,-1))*4/Dm-2*d_tt{2}*Dp/Dm,eps,rmax);
    d_tt{1}=round(d_tt{1}+temp_tt,eps,rmax);
end
G_tt = d_tt{1}-d_tt{3};
end
