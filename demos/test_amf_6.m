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

clc;clear;close all;
%% addpath
addpath /Users/oboue/Desktop/MATamf/amf_data/
addpath /Users/oboue/Desktop/MATamf/amfsrcs/
addpath /Users/oboue/Desktop/MATamf/seistr/
%% Load the raw DAS seismic data corrupted by strong hogh-frequency noise, high-amplitude 

NOs=[1,20,10,25,11,2];
labels={...                                             %P-arrival sample NO from the SEGY file
    'FORGE\_78-32\_iDASv3-P11\_UTC190423150554.sgy',... %24169
    'FORGE\_78-32\_iDASv3-P11\_UTC190426070723.sgy',... %24811
    'FORGE\_78-32\_iDASv3-P11\_UTC190426062208.sgy',... %26090
    'FORGE\_78-32\_iDASv3-P11\_UTC190426110008.sgy',... %4921
    'FORGE\_78-32\_iDASv3-P11\_UTC190426062553.sgy',... %8934
    'FORGE\_78-32\_iDASv3-P11\_UTC190423182409.sgy'};   %4210
eq=zeros(2000,960);
[n1,n2]=size(eq);
t=[0:n1]*0.0005;

for ii=4
    
if ~ismember(NOs(ii),[14,16,17,27,47,52])
    load(strcat('amf_data/eq-',num2str(NOs(ii)),'.mat'));
end
eq=d1;
din=d1;

%% Denosing using the BP method
%  Parameter tuning for the BP method
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
%% 
%%  Denosing using the BP+SOSVMF method
% Parameter tuning for the BP method

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
ns=15;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
tic
d2=amf_bandpasssosvmf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method
%
w=0.02;                       % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
tic
d3=amf_bandpasssosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
%
%% Denoising using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the key parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.25;              % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
% 
tic
d4=amf_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%
%% AMF
%% Parameter tuning for local orthogonalization operation
% 
rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
rec(1) = 150;
rec(2) = 150;
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
end
%% Plot figures

t=[0:n1]*0.0005;
ngap0=1000;

ngap=50;
x=1:n2*3+2*ngap;
d_1=[din,zeros(n1,ngap0),d1,din-d1,zeros(n1,ngap),d2,din-d2];
d_2=[d3,din-d3,zeros(n1,ngap),d4,din-d4,zeros(n1,ngap),d5,din-d5];

% combined figure
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1);amf_imagesc(d_1,100,2,x,t);caxis([-25,25]);

ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(n2/3.5,-0.03,'Raw','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.64,-0.03,'BP','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.38,-0.03,'BP+SOSVMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',...
    [0.12995951417004 0.754408060453399 0.7758987854251 0.0340050377833753],...
    'Color',[1 0 0],...
    'LineWidth',3);

subplot(2,1,2);amf_imagesc(d_2,100,2,x,t);caxis([-25,25]);
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(n2/2.0,-0.03,'BP+SOSVMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.64,-0.03,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/0.65,-0.03,'AMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(b)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',...
    [0.12995951417004 0.280856423173806 0.7758987854251 0.0340050377833752],...
    'Color',[1 0 0],...
    'LineWidth',3);
print(gcf,'-depsc','-r300','fig13.pdf');
print(gcf,'-depsc','-r300','fig13.eps');
%% Zoomed sections 

inds1=800:1000;

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1);amf_imagesc(d_1(inds1,:),100,2,x,t(inds1));caxis([-25,25]);
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
annotation(gcf,'textarrow',[0.478840115299105 0.496875],...
    [0.763127840546026 0.728615863141523],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'High-amplitude','erratic noise'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14);
annotation(gcf,'textarrow',[0.6980312433909 0.7190312433909],...
    [0.708732717973261 0.758732717973261],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Honrizontal','noise'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14);
subplot(2,1,2);amf_imagesc(d_2(inds1,:),100,2,x,t(inds1));caxis([-25,25]);
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
annotation(gcf,'textbox',...
    [0.159680282796544 0.167798407243139 0.0642003142183817 0.101508916323731],...
    'Color',[1 0 0],...
    'String',{'Random background noise'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textarrow',[0.45219795729137 0.482499619212023],...
    [0.148769549028649 0.169517884914463],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14);
annotation(gcf,'textarrow',[0.725167088079359 0.755468750000012],...
    [0.156602912413093 0.177351248298908],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Improved','signal energy'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14);
text(-200,0.25,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,0.255,'Raw','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.65,0.255,'BP','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.38,0.255,'BP+SOSVMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-200,0.395,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,0.395,'BP+SOSVMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.65,0.395,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.38,0.395,'AMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
print(gcf,'-depsc','-r300','fig14.pdf');
print(gcf,'-depsc','-r300','fig14.eps');
