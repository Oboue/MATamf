% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

%  Copyright (C) Oboue et al., 2022

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
%% Ccomparison between MF, SOMF, SOSVMF, and AMF method
%% Load the DAS data

NOs=[1,3,20,10,25,11,2];
labels={...                                          %P-arrival sample NO from the SEGY file
    'FORGE\_78-32\_iDASv3-P11\_UTC190423150554.sgy',... %24169
    'FORGE\_78-32\_iDASv3-P11\_UTC190423213209.sgy',... 
    'FORGE\_78-32\_iDASv3-P11\_UTC190426070723.sgy',... %24811
    'FORGE\_78-32\_iDASv3-P11\_UTC190426062208.sgy',... %26090
    'FORGE\_78-32\_iDASv3-P11\_UTC190426110008.sgy',... %4921
    'FORGE\_78-32\_iDASv3-P11\_UTC190426062553.sgy',... %8934
    'FORGE\_78-32\_iDASv3-P11\_UTC190423182409.sgy'};   %4210

eq=zeros(2000,960);
[n1,n2]=size(eq);
t=[0:n1]*0.0005;
ngap=50;
x=1:n2*5+4*ngap;

for ii=4

if ~ismember(NOs(ii),[14,16,17,27,47,52])
    load(strcat('amf_data/eq-',num2str(NOs(ii)),'.mat'));
end
eq=d1;
din=d1;
%% Denosing using the MF method 
ns=8;                         % spray radius
%
tic
d1=amf_mf(din,ns*2+1,1,2);
toc
%
%% Denoising using the SOMF method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=50;                   % rect:  smoothing radius (ndim*1)
rect(2)=50;                   % "      "        "
rect(3)=1;                    % "      "        "
verb=1;                       % verbosity flag

ns=8;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=0;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)
%
tic
d2=amf_somf(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%
%% Denosing using the SOSVMF method 
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
verb=1;                       % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=15;                        % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)

%
tic
d3=amf_sosvmf(din,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%
%% Denosing using the AMF method 
%
par.dt=0.0005; % sampling
par.flo=0;     % Low frequency in band, default is 0
par.fhi=200;   % High frequency in band, default is Nyquist
par.nplo=6;    % number of poles for low cutoff
par.nphi=6;    % number of poles for high cutoff
par.phase=0;   % y: minimum phase, n: zero phase
par.verb0=0;   % verbosity flag
%
par.niter=2;                      % number of nonlinear iterations
par.liter=10;                     % liter: number of linear iterations (in divn)
par.order1=3;                     % order: accuracy order
par.eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
par.tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=50;                   % rect:  smoothing radius (ndim*1)
rect(2)=50;                   % "      "        "
par.rect(3)=1;                    % "      "        "
par.verb1=1;                      % verbosity flag

par.adj=0;                        % adjoint flag
par.add=0;                        % adding flag
par.ns=15;                         % spray radius
par.order2=2;                     % PWD order
par.eps=0.01;                     % regularization (default:0.01);
par.ndn=n1*n2;                    % size of dn (n1*n2)
par.nds=n1*n2;                    % size of ds (n1*n2)
par.type_mf=1;                    % 0 (MF) or 1 (SVMF)
par.ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
par.w=0.02;              % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=0.15;               % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration
%
rec = zeros(3, 1);   % 3-D vector denoting smooth radius 
par.rec(1) = 100;
par.rec(2) = 100;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 

%
tic
d_amf=amf(din,par);
d4=d_amf;
toc
%
end
%% Plot figures 
[n1,n2]=size(din);
t=[0:n1]*0.0005;
ngap=50;
x=1:n2*5+4*ngap;

indt1=50:250;indx1=750:900;

d1_z1=d1(indt1,indx1);
d1_z2=d2(indt1,indx1);
d2_z1=d3(indt1,indx1);
d2_z2=d4(indt1,indx1);

comp1=[din,zeros(n1,ngap),d1,din-d1,zeros(n1,ngap),d2,din-d2]; 

comp2=[din,zeros(n1,ngap),d3,din-d3,zeros(n1,ngap),d4,din-d4]; 
% combined figure
figure('units','normalized','Position',[0.0 0.0 0.6, 1],'color','w');
subplot(2,1,1);amf_imagesc(comp1(1:1000,:),98,1,x,t(1:1000));
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Channel','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
text(n2/2,-0.02,'Raw data','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text((n2*2+ngap*3)/2+n2,-0.02,'MF','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text((n2*2+ngap*3)/2+n2*3+ngap*2,-0.02,'SOMF','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.02,'(a)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
text(50,0.46,labels{4},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');
hold on;
plot([indx1(1)+n2+ngap,indx1(1)+n2+ngap],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(end)+n2+ngap,indx1(end)+n2+ngap],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(1)+n2+ngap,indx1(end)+n2+ngap],t([indt1(1),indt1(1)]),'r','linewidth',2);
plot([indx1(1)+n2+ngap,indx1(end)+n2+ngap],t([indt1(end),indt1(end)]),'r','linewidth',2);

plot([indx1(1)+n2*3+ngap*3,indx1(1)+n2*3+ngap*3],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(end)+n2*3+ngap*3,indx1(end)+n2*3+ngap*3],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(1)+n2*3+ngap*3,indx1(end)+n2*3+ngap*3],t([indt1(1),indt1(1)]),'r','linewidth',2);
plot([indx1(1)+n2*3+ngap*3,indx1(end)+n2*3+ngap*3],t([indt1(end),indt1(end)]),'r','linewidth',2);

annotation(gcf,'arrow',[0.4018125 0.3718125],...
    [0.869444790046656 0.830444790046656],'Color',[1 0 0],'LineWidth',2);
annotation(gcf,'arrow',[0.714557291666667 0.684557291666667],...
    [0.871 0.832],'Color',[1 0 0],'LineWidth',2);
%
subplot(2,1,2);amf_imagesc(comp2(1:1000,:),98,1,x,t(1:1000));
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Channel','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
text(n2/2,-0.02,'Raw data','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text((n2*2+ngap*3)/2+n2,-0.02,'SOSVMF','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text((n2*2+ngap*3)/2+n2*3+ngap*2,-0.02,'AMF','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.02,'(b)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
text(50,0.46,labels{4},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

hold on;
plot([indx1(1)+n2+ngap,indx1(1)+n2+ngap],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(end)+n2+ngap,indx1(end)+n2+ngap],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(1)+n2+ngap,indx1(end)+n2+ngap],t([indt1(1),indt1(1)]),'r','linewidth',2);
plot([indx1(1)+n2+ngap,indx1(end)+n2+ngap],t([indt1(end),indt1(end)]),'r','linewidth',2);

plot([indx1(1)+n2*3+ngap*3,indx1(1)+n2*3+ngap*3],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(end)+n2*3+ngap*3,indx1(end)+n2*3+ngap*3],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(1)+n2*3+ngap*3,indx1(end)+n2*3+ngap*3],t([indt1(1),indt1(1)]),'r','linewidth',2);
plot([indx1(1)+n2*3+ngap*3,indx1(end)+n2*3+ngap*3],t([indt1(end),indt1(end)]),'r','linewidth',2);

annotation(gcf,'arrow',[0.402463541666667 0.372463541666667],...
    [0.391217729393468 0.352217729393468],'Color',[1 0 0],'LineWidth',2);
annotation(gcf,'arrow',[0.715208333333333 0.685208333333333],...
    [0.390440124416796 0.351440124416796],'Color',[1 0 0],'LineWidth',2);

a1=axes('Parent',gcf,'Position',[0.287,0.68,0.149,0.15]);
amf_imagesc(d1_z1,10,2);axis off;
a1=axes('Parent',gcf,'Position',[0.600,0.68,0.149,0.15]);
amf_imagesc(d1_z2,10,2);axis off;
a1=axes('Parent',gcf,'Position',[0.287,0.2,0.149,0.15]);
amf_imagesc(d2_z1,10,2);axis off;
a1=axes('Parent',gcf,'Position',[0.600,0.2,0.149,0.15]);
amf_imagesc(d2_z2,10,2);axis off;

print(gcf,'-depsc','-r300','fig16.eps');