% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

%  Copyright (C) Oboue et al., 2022

% Script to plot Figures 23 and 24

%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.

%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/

%  References

%  Oboue et al., 2022
%  Huang, G., M. Bai, Q. Zhao, W. Chen, and Y. Chen, 2021, Erratic noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting, 69, 101-121.
%  Chen, Y., S. Zu, Y. Wang, and X. Chen, 2020, Deblending of simultaneous-source data using a structure-oriented space varying median filter, Geophysical Journal International, 222, 1805�1�723.
%  Zhao, Q., Q. Du, X. Gong, and Y. Chen, 2018, Signal-preserving erratic noise attenuation via iterative robust sparsity-promoting filter, IEEE Transactions on Geoscience and Remote Sensing, 56, 1558-0644.
%-------------------------------------------------------------------------
clc;clear;close all;
%% addpath
addpath('../amf_data/');
addpath('../amfsrcs/');
addpath('../seistr/');
%% %Generate the clean synthetic data including linear events

dcl=amf_levents(200);
dcl=amf_scale(dcl);
[n1,n2]=size(dcl);
%% %Generate the noisy synthetic data having linear events

% Create the mask operator
mask=rand(1,n2);
mask(logical(mask<0.9))=0;
mask(logical(mask>=0.9))=1;

% Adding erratic noise 

err_n=zeros(size(dcl));
for i=1:n1
    randn('state',123456+i);
    err_n(i,:)=0.5*randn(1,n2).*mask;
end

% Adding random noise 

randn('state',201920);
ran_n=0.1*randn(n1,n2);

dnl=dcl+err_n+ran_n;
dinl=dnl;
%%  Denosing using the BP method
% Parameter tuning for the BP method

dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=251;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=7;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag

%
tic
d1l=amf_bandpass(dinl,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%
%% Denosing using the BP+SOSVMF method 
% Parameter tuning：add the key parameters of the SOSVMF method

niter=5;                      % number of nonlinear iterations
liter=20;                     % liter: number of linear iterations (in divn)
order1=2;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=11;                   % rect:  smoothing radius (ndim*1)
rect(2)=5;                    % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=4;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)

%                                                                                                             0,0,dipn,[],n1,n2,ns,2,0.01,n1*n2,n1*n2,type_mf,ifsmooth,d1,[]
tic
d2l=amf_bandpasssosvmf(dinl,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
% 
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method

w=0.0001; % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

% 
tic
d3l = amf_bandpasssosvmffk(dinl,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
%
%% Denoising data using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the key parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.9;              % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

% 
tic
d4l=amf_bandpasssosvmffkcurvelet(dinl,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%
%% Denoising using the AMF method 
%  Parameter tuning: add the key parameters for the local orthogonalization operation

par.dt=dt;         % sampling
par.flo=flo;       % Low frequency in band, default is 0
par.fhi=fhi;       % High frequency in band, default is Nyquist
par.nplo=nplo;     % number of poles for low cutoff
par.nphi=nphi;     % number of poles for high cutoff
par.phase=phase;   % y: minimum phase, n: zero phase
par.verb0=verb0;   % verbosity flag
%
par.niter=niter;                       % number of nonlinear iterations
par.liter=liter;                       % liter: number of linear iterations (in divn)
par.order1=order1;                     % order: accuracy order
par.eps_dv=eps_dv;                     % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=eps_cg;                     % eps_cg: eps for CG    (default: 1)
par.tol_cg=tol_cg;                     % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=11;                        % rect:  smoothing radius (ndim*1)
par.rect(2)=5;                         % "      "        "
par.rect(3)=1;                         % "      "        "
par.verb1=verb1;                       % verbosity flag

par.adj=adj;                           % adjoint flag
par.add=add;                           % adding flag
par.ns=ns;                             % spray radius
par.order2=order2;                     % PWD order
par.eps=eps;                           % regularization (default:0.01);
par.ndn=n1*n2;                         % size of dn (n1*n2)
par.nds=n1*n2;                         % size of ds (n1*n2)
par.type_mf=type_mf;                   % 0 (MF) or 1 (SVMF)
par.ifsmooth=ifsmooth;                 % 1 (if smooth) or 0 (only MF)
% 
par.w=w;                  % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=c1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=c2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=c3;                % Thresholding parameter (alpha)
par.niter1=10;            % Number of iteration

rec = zeros(3, 1);        % 3-D vector denoting smooth radius 
par.rec(1) = 5;
par.rec(2) = 5;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 
%
tic
d_amf=amf(dinl,par);
d5l=d_amf;
toc
%
%% Plot figures

t=[0:n1]*0.004;
x=1:n2;
ngap=5;
d_1l=[d3l,dnl-d3l,zeros(n1,ngap),d4l,d4l-dnl,zeros(n1,ngap),d5l,dnl-d5l];

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1); imagesc(x,t,d_1l);colormap(amf_seis);hold on
text(n2/-15,-0.1,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(n2/6.0,-0.06,'BP+SOSVMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.06,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.19,-0.06,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
annotation(gcf,'arrow',[0.176666666666667 0.155555555555556],...
    [0.785858585858586 0.807070707070707],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.202222222222222 0.232222222222222],...
    [0.764646464646465 0.785858585858586],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.31 0.33],...
    [0.881818181818181 0.866666666666665],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.32 0.32],...
    [0.837373737373738 0.819191919191918],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.31 0.33],...
    [0.881818181818181 0.866666666666665],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.44 0.418888888888889],...
    [0.781818181818181 0.803030303030302],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.466666666666667 0.496666666666667],...
    [0.761616161616162 0.782828282828282],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.58 0.58],...
    [0.845454545454545 0.827272727272725],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.582222222222222 0.594444444444445],...
    [0.87979797979798 0.867676767676765],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.698888888888889 0.677777777777778],...
    [0.786868686868686 0.808080808080807],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.734444444444444 0.764444444444444],...
    [0.759595959595959 0.78080808080808],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.837777777777778 0.837777777777778],...
    [0.847474747474747 0.829292929292927],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.838888888888889 0.856666666666667],...
    [0.878787878787879 0.867676767676766],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'textarrow',[0.261111111111111 0.179999999999999],...
    [0.246474747474747 0.33030303030303],...
    'String',{'Stronger','  signal',' leakage'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'arrow',[0.343333333333333 0.385555555555556],...
    [0.242424242424243 0.326262626262626],'LineWidth',3,'HeadWidth',25,...
    'HeadLength',15);
%
% local similarity maps
rect=[10,10,1];niter=20;eps=0;verb=0;
[simid3l]=amf_localsimi(dinl-d3l,d3l,rect,niter,eps,verb);
[simid4l]=amf_localsimi(dinl-d4l,d4l,rect,niter,eps,verb);
[simid5l]=amf_localsimi(dinl-d5l,d5l,rect,niter,eps,verb);

dsiml=[simid3l,zeros(n1,ngap),simid4l,zeros(n1,ngap),simid5l];

subplot(2,1,2); imagesc(x,t,dsiml);colormap(amf_seis);hold on
text(n2/-15,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(n2/6.0,-0.06,'BP+SOSVMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.06,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.19,-0.06,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 16;%c.Label.FontWeight = bold;
caxis([0,0.5]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');

print(gcf,'-depsc','-r300','fig19-amf.eps');
%% Numerical test using signal-to-noise ratio (S/N)

amf_snr(dcl,dnl) %Noisy
amf_snr(dcl,d1l) %BP
amf_snr(dcl,d2l) %BP+SOSVMF
amf_snr(dcl,d3l) %BP+SOSVMF+FK
amf_snr(dcl,d4l) %BP+SOSVMF+FK+Curvelet
amf_snr(dcl,d5l) %AMF
%%

%% %Generate the clean synthetic data (dc) including linear events

load Curveddata
dc=data;
% [n1,n2]=size(dc);
%% Load the noisy data corrupted by both random and erratic noise

load dncurved 

[n1,n2]=size(dn);
din=dn;
%%  Denosing using the BP method
% Parameter tuning for the BP method

dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=251;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=7;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag

%
tic
d1=amf_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%
%% Denosing using the BP+SOSVMF method 
% Parameter tuning：add the key parameters of the SOSVMF method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=21;                   % rect:  smoothing radius (ndim*1)
rect(2)=7;                    % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=9;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)

%                                                                                                             0,0,dipn,[],n1,n2,ns,2,0.01,n1*n2,n1*n2,type_mf,ifsmooth,d1,[]
tic
d2=amf_bandpasssosvmf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method

w=0.0001; % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

%
tic
d3=amf_bandpasssosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
%
%% Denoising data using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the key parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.4;              % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

% 
tic
d4=amf_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%
%% Denoising using the AMF method 
% % Parameter tuning: add the key parameters for the local orthogonalization operation

%% AMF method
par.dt=dt;                       % sampling
par.flo=flo;                     % Low frequency in band, default is 0
par.fhi=fhi;                     % High frequency in band, default is Nyquist
par.nplo=nplo;                   % number of poles for low cutoff
par.nphi=nphi;                   % number of poles for high cutoff
par.phase=phase;                 % y: minimum phase, n: zero phase
par.verb0=verb0;                 % verbosity flag
%
par.niter=niter;                 % number of nonlinear iterations
par.liter=liter;                 % liter: number of linear iterations (in divn)
par.order1=order1;               % order: accuracy order
par.eps_dv=eps_dv;               % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=eps_cg;               % eps_cg: eps for CG    (default: 1)
par.tol_cg=tol_cg;               % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=11;                  % rect:  smoothing radius (ndim*1)
par.rect(2)=5;                   % "      "        "
par.rect(3)=1;                   % "      "        "
par.verb1=verb1;                 % verbosity flag

par.adj=adj;                     % adjoint flag
par.add=add;                     % adding flag
par.ns=ns;                       % spray radius
par.order2=order2;               % PWD order
par.eps=eps;                     % regularization (default:0.01);
par.ndn=n1*n2;                   % size of dn (n1*n2)
par.nds=n1*n2;                   % size of ds (n1*n2)
par.type_mf=type_mf;             % 0 (MF) or 1 (SVMF)
par.ifsmooth=ifsmooth;           % 1 (if smooth) or 0 (only MF)
% 
par.w=w;                         % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=c1;                       % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=c2;                       % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=c3;                       % Thresholding parameter (alpha)
par.niter1=10;                   % Number of iteration

rec = zeros(3, 1);               % 3-D vector denoting smooth radius 
par.rec(1) = 8;
par.rec(2) = 8;
par.rec(3) = 1;
par.eps1=0;                      % regularization parameter, default 0.0
par.niter2=20;                   % number of CG iterations
par.verb=1;                      % verbosity flag (default: 0) 

%
tic
dout=amf(din,par);
d5=dout;
toc
%
%% Plot figures/curved events

t=[0:n1]*0.004;
x=1:n2;
ngap=5;
d_2=[d3,dn-d3,zeros(n1,ngap),d4,d4-dn,zeros(n1,ngap),d5,dn-d5];

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1); imagesc(x,t,d_2);colormap(amf_seis);
text(n1/-50,-0.1,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(n2/6.0,-0.06,'BP+SOSVMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.06,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.19,-0.06,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
annotation(gcf,'arrow',[0.23 0.231111111111111],...
    [0.841574167507568 0.811301715438951],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25);
annotation(gcf,'arrow',[0.192222222222222 0.186666666666667],...
    [0.675075681130171 0.709384460141271],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25);
annotation(gcf,'arrow',[0.326666666666667 0.326666666666667],...
    [0.84460141271443 0.81029263370333],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25);
annotation(gcf,'arrow',[0.49 0.49],...
    [0.847628657921292 0.813319878910192],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25);
annotation(gcf,'arrow',[0.452222222222222 0.446666666666667],...
    [0.671039354187689 0.705348133198789],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25);
annotation(gcf,'arrow',[0.582222222222222 0.582222222222222],...
    [0.843592330978809 0.80928355196771],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25);
annotation(gcf,'arrow',[0.716666666666667 0.711111111111111],...
    [0.668012108980828 0.702320887991928],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25);
annotation(gcf,'arrow',[0.764444444444444 0.764444444444444],...
    [0.839556004036327 0.805247225025228],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25);
annotation(gcf,'arrow',[0.846666666666667 0.846666666666667],...
    [0.840565085771948 0.806256306760849],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25);
%
% local similarity maps
rect=[10,10,1];niter=20;eps=0;verb=0;
[simid3]=amf_localsimi(dn-d3,d3,rect,niter,eps,verb);
[simid4]=amf_localsimi(dn-d4,d4,rect,niter,eps,verb);
[simid5]=amf_localsimi(dn-d5,d5,rect,niter,eps,verb);

dsim2=[simid3,zeros(n1,ngap),simid4,zeros(n1,ngap),simid5]; 

subplot(2,1,2); imagesc(x,t,dsim2);colormap(amf_seis);hold on
text(n1/-50,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(n2/6.0,-0.06,'BP+SOSVMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.06,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.19,-0.06,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 16;%c.Label.FontWeight = bold;
caxis([0,0.5]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
annotation(gcf,'textarrow',[0.388888888888889 0.333333333333333],...
    [0.395560040363269 0.350151362260343],...
    'String',{'    Stronger','signal leakage'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
 annotation(gcf,'arrow',[0.517777777777778 0.567777777777778],...
    [0.394550958627649 0.34409687184662],'LineWidth',3,'HeadWidth',25);

print(gcf,'-depsc','-r300','fig20-amf.eps');
%% Numerical test using signal-to-noise ratio (S/N) / Curved events

amf_snr(dc,din) %Noisy
amf_snr(dc,d1)  %BP
amf_snr(dc,d2)  %BP+SOSVMF
amf_snr(dc,d3)  %BP+SOSVMF+FK
amf_snr(dc,d4)  %BP+SOSVMF+FK+Curvelet
amf_snr(dc,d5)  % AMF