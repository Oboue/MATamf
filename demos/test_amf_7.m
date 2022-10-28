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
%% Load the DAS data

NOs=[1,3,20,10,25,11,2];
labels={...                                           %P-arrival sample NO from the SEGY file
    'FORGE\_78-32\_iDASv3-P11\_UTC190423150554.sgy',... %24169
    'FORGE\_78-32\_iDASv3-P11\_UTC190423213209.sgy'};   %4210
eq=zeros(2000,960);
[n1,n2]=size(eq);
t=[0:n1]*0.0005;
ngap=50;
x=1:n2*5+4*ngap;

%% Process DAS seismic data with weak signal energy corrupted by strong noise using the AMF method
%%
for ii=1
    
if ~ismember(NOs(ii),[14,16,17,27,47,52])
    load(strcat('amf_data/eq-',num2str(NOs(ii)),'.mat'));
end
eq=d1;
din=d1;
%% Parameter tuning 
%
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
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
ns=10;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
w=0.05;              % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.25;               % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
%
rec = zeros(3, 1);   % 3-D vector denoting smooth radius 
rec(1) = 500;
rec(2) = 500;
rec(3) = 1;
eps1=0;               % regularization parameter, default 0.0
niter2=20;            % number of CG iterations
verb=1;               % verbosity flag (default: 0) 

%
tic
d_amf1=amf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1,rec,eps1,niter2,verb);
toc
%
end

comp1=[din,zeros(n1,ngap),d_amf1,zeros(n1,ngap),din-d_amf1]; 
%% Process the DAS seismic data with strong signal energy corrupted by strong noise using the AMF method
%%
eq=zeros(2000,960);
[n1,n2]=size(eq);

for ii=3
    if ~ismember(ii,[14,16,17,27,47,52])
        strcat('amf_data/eq-',num2str(ii),'.mat')
        load(strcat('amf_data/eq-',num2str(ii),'.mat'));
    end
    d1=d1;
    eq=d1;
din=d1;
%% Parameter tuning 
%
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
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
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
w=0.02;   % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.25;               % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
%
rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
rec(1) = 50;
rec(2) = 50;
rec(3) = 1;
eps1=0;               % regularization parameter, default 0.0
niter2=20;            % number of CG iterations
verb=1;               % verbosity flag (default: 0) 

%
tic
d_amf2=amf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1,rec,eps1,niter2,verb);
toc
%
end
%% Plot figures 

comp2=[din,zeros(n1,ngap),d_amf2,zeros(n1,ngap),din-d_amf2]; 

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1);amf_imagesc(comp1,95,1,x,t);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Channel','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
text(n2/2,-0.03,'Raw data','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(n2/2+ngap+n2,-0.03,'AMF','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(n2/2+ngap*2+n2*2,-0.03,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{1},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,1,2);amf_imagesc(comp2,95,1,x,t);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Channel','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.03,'Raw data','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(n2/2+ngap+n2,-0.03,'AMF','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(n2/2+ngap*2+n2*2,-0.03,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(b)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{2},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

print(gcf,'-depsc','-r300','fig15.pdf');
print(gcf,'-depsc','-r300','fig15.eps');