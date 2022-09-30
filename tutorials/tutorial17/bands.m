clear;
ef=12.6178; # Fermi energy
v=[ef ef];

fid=fopen('fe_slice_x.dat', 'r');
x=fscanf(fid,'%f\n');
fclose(fid);
dimx=length(x(:));
xmin=x(1);
xmax=x(dimx);

fid=fopen('fe_slice_y.dat', 'r');
y=fscanf(fid,'%f\n');
fclose(fid);
dimy=length(y(:));
ymin=y(1);
ymax=y(dimy);

fid=fopen('fe_slice_bands.dat', 'r');
z=fscanf(fid,'%f\n');
fclose(fid);
zz=reshape(z,[],dimy,dimx);
dim=size(zz)
for n=1:dim(1)
  zz_n(:,:)=zz(n,:,:);
  contour(x,y,zz_n,v,'bk')
  hold on;
end
axis([xmin xmax ymin ymax],"square","nolabel","tic[]")
grid off


print('bands.eps', '-depsc2')

title('Fermi contours of bcc Fe on the (010) plane')
