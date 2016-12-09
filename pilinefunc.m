function [ Vr, phiVs, Ir, Is, phiIs ] = pilinefunc( S, phir, Scomplex, R, XL, XC, Vs, l )
%pilinefunc( S, phir, Scomplex, R, XL, XC, Vs, l )
%    [ Vr, phiVs, Ir, Is, phiIs, a31, a32, a33, a34 ] = pilinefunc( S, phir, Scomplex, R, XL, XC, Vs, l )
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
% const1 = (const1lr.*l*R);

aa = const1lr.^2;

bb = 2*const1lr.*(R*cos(phir) - XL*sin(phir)).*l*S...
    +(R*S/XC)*(XL*cos(phir) + R*sin(phir)).*l.^3 - Vs^2;

cc = ((R*S)^2 + (XL*S)^2).*l.^2;

Vrsq = (0.5./aa).*(-bb + sqrt(bb.^2 - 4*aa.*cc) );

% if real(Vrsq) < 0
%     disperror1();
% elseif abs(imag(Vrsq)) >0
%     disperror1();
% else
%     if exist('errmsg','var')
%         set('errmsg','position',[10 10 0 0])
%         delete(errmsg);
%     end
Vr = sqrt( Vrsq );
Ir = S./Vr;

phiVs = angle( Vr + l.*(R+1i*XL).*(Scomplex'./Vr + 0.5i*l.*Vr/XC) );

Ircomp = Ir.*exp(1i*phir);
Is = Ircomp + 0.5i*l.*(Vr + Vs.*exp(1i*phiVs))/XC;
phiIs = angle(Is);
Is = abs(Is);

% end

% Ic2 = Vr.*l./(2*XC);
% 
% Isr = Ir.*exp(1i*phir) + 1i*Ic2;
% phisr = angle(Isr);
% Isr = abs(Isr);
% 
% 
% sinphisr = sin(phisr-phiVs);
% cosphisr = cos(phisr-phiVs);
% 
% a32x = Vs+Isr.*XL.*l.*sinphisr;
% a32y = -Isr.*XL.*l.*cosphisr;
% a33x = a32x - Isr.*l.*R.*cosphisr;
% a33y = a32y - Isr.*l.*R.*sinphisr;
% 
% % a31 = arrow([0 0],[Vs 0],1,'edgecolor','b','facecolor','b');
% % a32 = arrow([Vs 0],[a32x a32y],1,'edgecolor','b','facecolor','b');
% % a33 = arrow([a32x a32y],[a33x a33y],1,'edgecolor','b','facecolor','b');
% % a34 = arrow([0 0],[Vr.*cos(-phiVs) Vr.*sin(-phiVs)],1,'edgecolor','b','facecolor','b','linewidth',1.5);
% a31 = plot([0 Vs],[0 0],'b');
% a32 = plot([Vs a32x],[0 a32y],'b');
% a33 = plot([a32x a33x],[a32y a33y],'b');
% a34 = plot([0 Vr.*cos(-phiVs)],[0 Vr.*sin(-phiVs)],'b','linewidth',1.5);


end

