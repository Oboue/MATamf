% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets
% 
%  Copyright (C) Oboue et al., 2022
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
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
addpath /Users/oboue/Desktop/MATamf/amf_data/
addpath /Users/oboue/Desktop/MATamf/amfsrcs/
addpath /Users/oboue/Desktop/MATamf/seistr/
%% read field reflection seismic data 1
fid=fopen('real2.bin','r');
d2d=fread(fid,[256,256],'float');

dn=amf_scale(d2d,2);
[n1,n2]=size(dn);
din=dn;
%%  Denosing using the BP method
% Parameter tuning for the BP method
%
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
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
rect(1)=50;                   % rect:  smoothing radius (ndim*1)
rect(2)=50;                   % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=8;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)

%
tic
d2=amf_bandpasssosvmf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%

%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method
%
w=0.00;   % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
tic
d3=amf_bandpasssosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
%
%% Denoising data using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the key parameters of the curvelet method
%
c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.5;               % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
% 
tic
d4=amf_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%
%% Denoising using the AMF method 
%  Parameter tuning: add the key parameters for the local orthogonalization operation
%
rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
rec(1) = 10;
rec(2) = 10;
rec(3) = 1;
eps1=0;               % regularization parameter, default 0.0
niter2=20;            % number of CG iterations
verb=1;               % verbosity flag (default: 0) 
%
tic
d_amf=amf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1,rec,eps1,niter2,verb);
d5=d_amf;
toc
%
%% Plot figures

t=[0:n1]*0.004;
x=1:n2;
ngap=10;
ngap0=250;
d_1=[din,zeros(n1,ngap0),d1,din-d1,zeros(n1,ngap),d2,din-d2];
d_2=[d3,din-d3,zeros(n1,ngap),d4,d4-din,zeros(n1,ngap),d5,din-d5];

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1); imagesc(x,t,d_1);colormap(amf_seis);hold on
text(n2/-200,-0.1,'(a)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/9.5,-0.03,'Raw','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.03,'BP','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.17,-0.03,'BP+SOSVMF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
colormap(amf_seis);caxis([-0.4,0.4]);

subplot(2,1,2); imagesc(x,t,d_2);colormap(amf_seis);hold on
text(n2/-200,-0.1,'(b)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/6.0,-0.03,'BP+SOSVMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.03,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.19,-0.03,'AMF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
colormap(amf_seis);caxis([-0.4,0.4]);

annotation(gcf,'rectangle',...
    [0.130769230769231 0.790426908150065 0.774684210526315 0.0685640362225094],...
    'LineWidth',3);
annotation(gcf,'rectangle',...
    [0.13036437246964 0.316946959896517 0.774684210526313 0.0685640362225092],...
    'LineWidth',3);
print(gcf,'-depsc','-r300','fig7.pdf');
print(gcf,'-depsc','-r300','fig7.eps');
%% zoomed sections
inds1=50:100;
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1);amf_imagesc(d_1(inds1,:),25,2,x,t(inds1));
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
colormap(amf_seis);caxis([-0.4,0.4]);

subplot(2,1,2);amf_imagesc(d_2(inds1,:),25,2,x,t(inds1));
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
colormap(amf_seis);caxis([-0.4,0.4]);
annotation(gcf,'textarrow',[0.174218750000007 0.1890625],...
    [0.868301116814964 0.839035769828927],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'      High-amplitude','erratic noise'},...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textbox',...
    [0.151 0.728615863141524 0.06228125 0.0318818040435465],'Color',[1 0 0],...
    'String',{'Random','  noise'},...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textarrow',[0.44921875 0.425000000000001],...
    [0.890357698289269 0.876360808709176],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.50625 0.48359375],...
    [0.891135303265941 0.880248833592535],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.44765625 0.425000000000001],...
    [0.730948678071539 0.720062208398133],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.48125 0.49140625],...
    [0.682737169517885 0.661741835147745],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.560937500000001 0.560937500000001],...
    [0.778382581648522 0.813374805598756],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.603125000000001 0.603125000000002],...
    [0.660964230171073 0.695956454121306],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.80234375 0.802343750000001],...
    [0.840590979782271 0.875583203732504],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.86171875 0.861718750000001],...
    [0.610419906687403 0.645412130637636],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.304687500000001 0.304687500000001],...
    [0.360808709175738 0.395800933125971],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.301562500000005 0.31484375],...
    [0.215396578538103 0.181181959564541],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.304687500000001 0.304687500000001],...
    [0.360808709175738 0.395800933125971],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.613281250000001 0.613281250000002],...
    [0.384914463452566 0.419906687402799],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.614062500000004 0.614062500000004],...
    [0.182737169517885 0.217729393468118],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.613281250000001 0.613281250000002],...
    [0.384914463452566 0.419906687402799],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'textarrow',[0.739843750000001 0.739843750000001],...
    [0.377916018662519 0.346811819595645],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Improved','signal','    energy'},... 
    'LineWidth',2,...
    'HeadWidth',8,...
    'FontWeight','bold');
annotation(gcf,'arrow',[0.546153846153845 0.537427252024287],...
    [0.30892143808256 0.358286512292556],'Color',[1 0 0],'LineWidth',2);

text(n2/200,-0.1,'(a)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/9.5,-0.095,'Raw','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.095,'BP','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.19,-0.095,'BP+SOSVMF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');

text(n2/200,0.18,'(b)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/6.0,0.189,'BP+SOSVMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,0.189,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.19,0.189,'AMF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');

print(gcf,'-depsc','-r300','fig8.pdf');
print(gcf,'-depsc','-r300','fig8.eps');
%% local similarity maps

rect=[10,10,1];niter=20;eps=0;verb=0;
[simid1]=amf_localsimi(dn-d1,d1,rect,niter,eps,verb);
[simid2]=amf_localsimi(dn-d2,d2,rect,niter,eps,verb);
[simid3]=amf_localsimi(dn-d3,d3,rect,niter,eps,verb);
[simid4]=amf_localsimi(dn-d4,d4,rect,niter,eps,verb);
[simid5]=amf_localsimi(dn-d5,d5,rect,niter,eps,verb);

d_sim=[simid1,zeros(n1,ngap),simid2,zeros(n1,ngap),simid3,zeros(n1,ngap),simid4,zeros(n1,ngap),simid5];

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
imagesc(x,t,d_sim);colormap(jet);hold on
text(n2/9,-0.018,'BP','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.25,-0.018,'BP+SOSVMF','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.018,'BP+SOSVMF+FK','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.40,-0.018,'BP+SOSVMF+Curvelet','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.10,-0.018,'AMF','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 16;%c.Label.FontWeight = bold;
caxis([0,0.5]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');

print(gcf,'-depsc','-r300','fig9.pdf');
print(gcf,'-depsc','-r300','fig9.eps');


















