function [R,t_out]=SPnP(Xw,uv)
% cann't be used for planar case !!! for 3d case, can works without
% initializaion
n=size(uv,2);
Xc=normc([uv;ones(1,n)]);
% [R0,t0]=initializor(Xw,uv);
% Xw=R0*Xw;
% t0=[88;88;88];
% [R0,t0]=RPnP(Xw,uv);
invl=ones(size(Xc,2),1);
% Xw_t=Xw+repmat(t0,1,n);
% inva=1./sqrt(sum(Xw_t.*Xw_t)).';
it=1;R_old=0;t_old=0;R=888;t=888;
while (norm(R-R_old,'fro')>1e-5 || norm(t-t_old)>1e-5) && it<1000
    sumA2=sum(invl.^2);
    invL=repmat(invl.',3,1);
    sumA2_sumAXc=sum(invL.*Xc,2)/sumA2;
    invAXw=invL.*Xw;
    sumA2_sumA2Xw=sum(invL.*invAXw,2)/sumA2;
    R_old=R;
    [U,S,V]=svd((Xc-invL.*sumA2_sumAXc)*(invAXw-invL.*sumA2_sumA2Xw).');
    R=U*[1 0 0;0 1 0;0 0 det(U*V.')]*V.';
%     if abs(S(3))<1e-4
%         R(:,[1 2])=R;
%     end
%     if det(R)<0, 
%         R=-R; 
%     end
%     R(:,3)=cross(R(:,1),R(:,2));
%     R
    t_old=t;
    t=R.'*sumA2_sumAXc-sumA2_sumA2Xw;
    Xw_t=Xw+repmat(t,1,n);
    invl=1./sqrt(sum(Xw_t.*Xw_t)).';

    it=it+1;
end
% R=R*R0;
t_out=R*t;
% it