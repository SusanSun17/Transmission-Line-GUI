function [ Vr, phiVs, Ir, Is, phiIs ] = Tlinefunc( S, phir, R, XL, XC, Vs, l )
%Tlinefunc( S, phir, R, XL, Vs, l )
%    [ f, ax, a21, a22, a23, a24, a25, a26 ] = Tlinefunc( S, phir, R, XL, XC, Vs, l )
%    Scomplex = S*exp(1i*-phir) - complex load power, MVA
%    R - resistance, ohms per km
%    XL - inductive reactance, ohms per km
%    Vs - supply voltage, kV
%    l - line length, km
%
%    f - figure handle
%    ax - axes handle
%    a11 - arrow Vs
%    a12 - arrow Ir*XL
%    a13 - arrow Ir*R
%    a14 - arrow Vr

const1lr = (1-0.5*l.^2.*XL/XC);
const1 = (const1lr.*l*R);
const2 = (XL.*l+0.25*l.^3.*(R^2-XL^2)/XC);

aa = const1.^2 + const2.^2;

bb = (2*const1.*const1lr + const2*R.*l.^2/XC).*S*cos(phir)...
    + (const1*R*l.^2/XC - 2*const1lr.*const2).*S*sin(phir) - Vs^2;

cc = ((const1lr).^2 + (0.5*l.^2.*R/XC).^2)*S^2;

Irsq = (0.5./aa).*(-bb - sqrt(bb.^2 - 4*aa.*cc) );

% if real(Irsq) < 0
%     disperror1();
% elseif abs(imag(Irsq)) >0
%     disperror1();
% else
%     if exist('errmsg','var')
%         set('errmsg','position',[10 10 0 0])
%         delete(errmsg);
%     end
    
Ir = sqrt( Irsq );
Vr = S./Ir;

const3 = (const1lr + 0.5i*l.^2*R/XC);
Ircomp = Ir.*exp(1i*phir);

phiVs = angle( Vr + 0.5*(R+1i*XL).*l.*((1+const3).*Ircomp +1i*l.*Vr/XC) );

Is = const3.*Ircomp +1i*l.*Vr/XC;
phiIs = angle(Is);
Is = abs(Is);


% sinphiIs = sin(phiIs-phiVs);
% cosphiIs = cos(phiIs-phiVs);
% sinphir = sin(phir-phiVs);
% cosphir = cos(phir-phiVs);
% Isl = 0.5*Is.*l;
% Irl = 0.5*Ir.*l;
% a22x = Vs+XL.*Isl.*sinphiIs;
% a22y = -XL.*Isl.*cosphiIs;
% a23x = a22x -R.*Isl.*cosphiIs;
% a23y = a22y -R.*Isl.*sinphiIs;
% a24x = a23x +XL.*Irl.*sinphir;
% a24y = a23y -XL.*Irl.*cosphir;
% a25x = a24x -R.*Irl.*cosphir;
% a25y = a24y -R.*Irl.*sinphir;
% 
% % a21 = arrow([0 0],[Vs 0],1,'edgecolor','r','facecolor','r');
% % a22 = arrow([Vs 0],[a22x a22y],1,'edgecolor','r','facecolor','r');
% % a23 = arrow([a22x a22y],[a23x a23y],1,'edgecolor','r','facecolor','r');
% % a24 = arrow([a23x a23y],[a24x a24y],1,'edgecolor','r','facecolor','r');
% % a25 = arrow([a24x a24y],[a25x a25y],1,'edgecolor','r','facecolor','r');
% % a26 = arrow([0 0],[Vr.*cos(-phiVs) Vr.*sin(-phiVs)],1,'edgecolor','r','facecolor','r','linewidth',1.5);
% a21 = plot([0 Vs],[0 0],'r');
% a22 = plot([Vs a22x],[0 a22y],'r');
% a23 = plot([a22x a23x],[a22y a23y],'r');
% a24 = plot([a23x a24x],[a23y a24y],'r');
% a25 = plot([a24x a25x],[a24y a25y],'r');
% a26 = plot([0 Vr.*cos(-phiVs)],[0 Vr.*sin(-phiVs)],'r','linewidth',1.5);

end

