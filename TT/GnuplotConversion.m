addpath(genpath('~/Downloads/TT-Toolbox'));
format long

% Read data
load('./resultofstepfinal.mat')

% Save norm
norm=softnorm;
normc=norm;
normdata = [t; norm; normc];

fileID = fopen('norm.dat','w');
fprintf(fileID,'%12.8f %12.8f %12.8f\n',normdata);
fclose(fileID);

% Save autocorrelation function
realxi=real(xi);
realxic=real(xic);
imagxi=imag(xi);
imagxic=imag(xic);
autodata = [t; realxi; imagxi; realxic; imagxic];

fileID = fopen('autocorrelation.dat','w');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f\n',autodata);
fclose(fileID);

numlist={0,100,200,300,400,500,600};
% Save probability densities
for ii=1:length(numlist)
    filename=sprintf('%s%d%s','resultofstep',numlist{ii},'.mat');
    load(filename)
    v1=reshape(abs(full(G_tt(:,:,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2))).^2,[nx nx]);
    v3=reshape(abs(full(Gc_tt(:,:,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2))).^2,[nx nx]);
    newfilename=sprintf('%s%d%s','wave.',ii-1,'.dat');
    fileID = fopen(newfilename,'w');
    for jj=1:(size(x,2))
        xpositions=x;
        ypositions=ones(1,32)*x(jj);
        softslice=transpose(v1(:,jj));
        chebslice=transpose(v3(:,jj));
        densitydata=[xpositions; ypositions; softslice; chebslice];
        fprintf(fileID,'%e %e %e %e\n',densitydata);
        fprintf(fileID,'\n',' ');
    end
    fclose(fileID);
end