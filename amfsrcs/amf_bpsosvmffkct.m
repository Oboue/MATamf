function [dout] = amf_bpsosvmffkct(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1)
% amf_bpsosvmffkct: robust noise attenuation method used to successfully remove various type of strong noise in the advanced median filter (AMF) method.
%
% By Oboue et al., 2022
%
% INPUT
% din:      input data
% dt:       sampling
% flo:      Low frequency in band, default is 0
% fhi:      High frequency in band, default is Nyquist
% nplo:     number of poles for low cutoff
% nphi:     number of poles for high cutoff
% phase:    y: minimum phase, n: zero phase
% verb0:    verbosity flag
%
% niter:    number of nonlinear iterations
% liter:    number of linear iterations (in divn)
% order1:   accuracy order
% eps_dv:   eps for divn  (default: 0.01)
% eps_cg:   eps for CG    (default: 1)
% tol_cg:   tolerence for CG (default: 0.000001)
% rect:     smoothing radius (ndim*1)
% verb1:    verbosity flag
%
% adj:      adjoint flag
% add:      adding flag
% ns:       spray radius
% order2:   PWD order
% eps:      regularization (default:0.01);
% ndn:      size of dn (n1*n2)
% nds:      size of ds (n1*n2)
% type_mf:  0 (MF) or 1 (SVMF)
% ifsmooth: 1 (if smooth) or 0 (only MF)
%
% w:        half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
% c1:       Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
% c2:       Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
% c3:       Thresholding parameter (alpha)
% niter1:   Number of iteration
%%
[d1]=amf_bp(din,dt,flo,fhi,nplo,nphi,phase,verb0);

% amf_bp: Bandpass filtering
%%
[d2] = amf_sosvmf(d1,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);

% amf_sosvmf: structure oriented space-varying median filter (SOSVMF) method 
%%
[d3]=d2-amf_fk_dip(d2,w);

% amf_fk_dip: FK dip filter
%%
[dout]=amf_ct(din,d3,n1,n2,c1,c2,c3,niter1);

% amf_ct: curvelet method
return