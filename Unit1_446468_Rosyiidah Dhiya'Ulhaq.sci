// Rosyiidah Dhiya'Ulhaq
// 19/446468/TK/49573

//////////////////////////////////////////////////
// PENGUJIAN 1
/////////////////////////////////////////////////
// konfigurasi variabel terkunci
etha=377; // impedansi ruang bebas
f=3*(10^9); // frekuensi operasi antena
c=3*(10^8); // cepat rambat gelombang di ruang hampa
lamda=c/f; // panjang gelombang
I0=1; // arus maksimum yang masuk ke antena
// panjang keseluruhan antena dipole
//l=(0.5)*lamda;
//l=(2.1)*lamda;
//l=(2.4)*lamda;
//l=(2.7)*lamda;
l=(3)*lamda;
k=2*(%pi)/lamda; //bilangan gelombang (wavenumber)
r=10*lamda; // radius spherical coordinate
theta=[0:0.01:2*(%pi)]; // sudut elevasi
// E0
rms1=%i*etha*I0*exp(-%i*k*r);
rms2=2*(%pi)*r;
rms3=cos((k*l/2)*cos(theta))-cos(k*l/2);
rms4=sin(theta);
E0=abs((rms1/rms2)*(rms3./rms4));
E0dB=10*log10(E0);
// plot pola radiasi antena secara 2D
u=figure(1);
polarplot(theta,E0dB);
set(u,"background",-2);
title('$Plot\ Pola\ Radiasi\ Antena\ secara\ 2D$','fontsize',4)
// azimuth
azimuth=[0:0.01:2*(%pi)];
nE0dB=E0dB-min(E0dB);
nE0dB(1,1)=0;
mtx_az=[];
mtx_th=[];
mtx_et=[];
for i=1:629
    mtx_az(i,:)=azimuth(1,:);
    mtx_th(:,i)=theta(1,:);
    mtx_et(:,i)=nE0dB(1,:);
end;
x=mtx_et.*cos(mtx_th).*cos(mtx_az);
y=mtx_et.*cos(mtx_th).*sin(mtx_az);
z=mtx_et.*sin(mtx_th);
// plot pola radiasi antena secara 3D
G=[0.6 0.7 0.8 0.9 1 1 0.9 0.8 0.7 0.6]';
B=[0.4 0.3 0.2 0.1 0 0 0.1 0.2 0.3 0.4]';
R=ones(B);
cmap=[R,G,B];
f=figure(2);
surf(x,y,z,'thickness',0);
e=gce();
e.cdata_mapping='direct'
e.color_flag=4;
set(f,"color_map",cmap);
set(f,"background",-2);
f.color_map=parulacolormap(32);
colorbar;
title('$Plot\ Pola\ Radiasi\ Antena\ secara\ 3D$','fontsize',4)
xgrid

//////////////////////////////////////////////////
// PENGUJIAN 2
/////////////////////////////////////////////////
// konfigurasi variabel terkunci
etha=377; // impedansi ruang bebas
f=3*(10^9); // frekuensi operasi antena
c=3*(10^8); // cepat rambat gelombang di ruang hampa
lamda=c/f; // panjang gelombang
I0=1; // arus maksimum yang masuk ke antena
// panjang keseluruhan antena dipole
l=(1.5)*lamda;
k=2*(%pi)/lamda; //bilangan gelombang (wavenumber)
r=10*lamda; // radius spherical coordinate
theta=[0:0.01:2*(%pi)]; // sudut elevasi
// nilai n
//n = 10;
//n = 30;
n = 50;
d = lamda/4;
beta = -3(%pi/2);
// E01
rms1=%i*etha*I0*exp(-%i*k*r);
rms2=2*(%pi)*r;
rms3=cos((k*l/2)*cos(theta))-cos(k*l/2);
rms4=sin(theta);
E0=abs((rms1/rms2)*(rms3./rms4));
E0dB=10*log10(E0);
//E0
Psi = (k*d*cos(theta))+beta;
rms21 = 1/n;
rms22 = sin((n/2)*Psi);
rms23 = sin((l/2)*Psi);
Af = abs((rms21)*(rms22./rms23));
AfdB = 10*log10(Af);
E02 = E0dB.*AfdB;
// plot pola radiasi antena secara 2D
v= figure(1);
set(v, "background", -2);
polarplot(theta,E02);
title('$Plot\ Pola\ Radiasi\ Antena\ secara\ 2D$','fontsize',4)
// azimuth
azimuth = [0:0.01:2*(%pi)];
nE02 = E02 - min(E02);
nE02(1,1) = 0;
mtx_az=[];
mtx_th=[];
mtx_et=[];
for i=1:629
    mtx_az(i,:)=azimuth(1,:);
    mtx_th(:,i)=theta(1,:);
    mtx_et(:,i)=nE02(1,:);
end;
x=mtx_et.*cos(mtx_th).*cos(mtx_az);
y=mtx_et.*cos(mtx_th).*sin(mtx_az);
z=mtx_et.*sin(mtx_th);
// plot pola radiasi antena secara 3D
G=[0.6 0.7 0.8 0.9 1 1 0.9 0.8 0.7 0.6]';
B=[0.4 0.3 0.2 0.1 0 0 0.1 0.2 0.3 0.4]';
R=ones(B);
cmap=[R,G,B];
f = figure(2);
surf(x,y,z, 'thickness', 0);
e=gce();
e.cdata_mapping='direct'
e.color_flag=4;
set(f, "color_map", cmap);
set(f, "background", -2);
f.color_map=parulacolormap(64);
colorbar;
title('$Plot\ Pola\ Radiasi\ Antena\ secara\ 3D$','fontsize',4)
xgrid
