function [G_tt] = Trotter(G_tt,dx,nx,KP2_tt,...
    dp,np,eps,rmax,PPhalf_tt)

% Always half potential method
% Apply potential Trotter propagator
G_tt = round(times(G_tt,PPhalf_tt),eps,rmax); % diabatic evolution
% Apply kinetic Trotter propagator
FT_G_tt = tt_FT(G_tt,dx,nx,1);  % FT to momentum space
FT_G_tt = round(times(FT_G_tt,KP2_tt),eps,rmax);
G_tt = tt_FT(FT_G_tt,dp,np,-1); % IFT to coord space
% Apply potential Trotter propagator
G_tt = round(times(G_tt,PPhalf_tt),eps,rmax); % diabatic evolution

end

