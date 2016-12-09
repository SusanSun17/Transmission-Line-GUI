function [ Vr, phiVs, Ir ] = shortlinefunc( S, phir, Scomplex, R, XL, Vs, l )
%shortlinefunc( S, phir, Scomplex, R, XL, Vs, l )
%    [ f, ax, a11, a12, a13, a14 ] = shortlinefunc( S, phir, Scomplex, R, XL, Vs, l )
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

Vrsq = -S*(cos(phir)*R*l - XL*sin(phir)*l)+0.5*Vs^2 + sqrt( (S*(cos(phir)*R*l - XL*sin(phir)*l)-0.5*Vs^2).^2 - S^2*(R.^2+XL.^2).*l.*l);

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

phiVs = angle(Vr + (Scomplex'./Vr).*(R+1i*XL).*l);

Ir = abs(Scomplex'./Vr);


% sinphir = sin(phir-phiVs);
% cosphir = cos(phir-phiVs);
% a12x = Vs+Ir.*XL.*l.*sinphir;
% a12y = -Ir.*XL.*l.*cosphir;
% a13x = a12x -Ir.*R.*l.*cosphir;
% a13y = a12y -Ir.*R.*l.*sinphir;
% 
% % a11 = arrow([0 0],[Vs 0],1);
% % a12 = arrow([Vs 0],[a12x a12y],1);
% % a13 = arrow([a12x a12y],[a13x a13y],1);
% % a14 = arrow([0 0],[Vr.*cos(-phiVs) Vr.*sin(-phiVs)],1,'linewidth',2);
% a11 = plot([0 Vs],[0 0],'k');
% a12 = plot([Vs a12x],[0 a12y],'k');
% a13 = plot([a12x a13x],[a12y a13y],'k');
% a14 = plot([0 Vr.*cos(-phiVs)],[0 Vr.*sin(-phiVs)],'k','linewidth',2);
% 
% axis equal
% grid on


end

