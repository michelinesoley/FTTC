function [FT] = tt_FT(tt,dx,nx,ns)
% Fourier transform (ns=1) and inverse Fourier transform (ns=-1) of
% tensor train
FT=tt;
co=core(tt);
d=length(co); % # of cores
nel=0; % number of elements in previous cores
for k=1:d
    sco=size(co{k});
   CC=matricize(co{k},1); % uses h-tucker tensor toolbox
    s = size(CC); % s(1) = n_k, s(2) = r_{k-1}+r_{k}
    s12=s(1)*s(2);
    if(ns ==1)
        CC=fft(CC)*dx/sqrt(2*pi);
    else
        CC=ifft(CC)*dx/sqrt(2*pi)*nx;
    end
    if (k == 1)
        FT.core(1+nel:s12+nel)=reshape(CC,[s12,1]);
    elseif (k < d)
        if(length(sco) == 2)
            FT.core(1+nel:s12+nel)=reshape(CC,[s12,1]);
        else
            CC=reshape(CC,[sco(1),sco(2),sco(3)]);
            FT.core(nel+1:nel+s12)=reshape(permute(CC,[2,1,3]),[s12,1]);
        end
    else
        FT.core(1+nel:s12+nel)=reshape(conj(CC)',[s12,1]);
    end
    nel=nel+s12;
end