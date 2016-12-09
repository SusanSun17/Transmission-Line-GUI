
function updatel(f, S, phir, R, XL, XC, VsDropdown, l, arrowscell, labelscell)
%updatel
%  updatel(S, phir, Scomplex, R, XL, XC, Vs, l)
%  Updates phasors when l is changed by slider callback
%
    Scomplex = S.*exp(1i*-phir);
    scrsz = get(groot,'screensize');
    
    if VsDropdown==1
        Vs = 132;
        maxhead = 0.02;
    elseif VsDropdown==2
        Vs = 275;
        maxhead = 0.01;
    elseif VsDropdown==3
        Vs = 400;        
        maxhead = 0.007;
    elseif VsDropdown==4
        Vs = 800;
        maxhead = 0.003;
    end

    %% short line
    [ Vrs(1), phiVss(1), Irs(1)] = shortlinefunc( S, phir, Scomplex, R, XL, Vs, l );
    Vr = Vrs(1);
    phiVs = phiVss(1);
    Ir = Irs(1);
    Iss(1) = Irs(1);
    phiIss(1) = phir(1);
    
    sinphir = sin(phir-phiVs);
    cosphir = cos(phir-phiVs);
    a12x = Vs+Ir.*XL.*l.*sinphir;
    a12y = -Ir.*XL.*l.*cosphir;
    a13x = a12x -Ir.*R.*l.*cosphir;
    a13y = a12y -Ir.*R.*l.*sinphir;

%     set(a11,'xdata',[0 Vs],'ydata',[0 0])
%     set(a12,'xdata',[Vs a12x],'ydata',[0 a12y]);
%     set(a13,'xdata',[a12x a13x],'ydata',[a12y a13y]);
%     set(a14,'xdata',[0 Vr.*cos(-phiVs)],'ydata',[0 Vr.*sin(-phiVs)]);

%     set(a11,'UData',Vs,'maxheadsize',maxhead);
%     set(a12,'XData',Vs,'UData',(a12x-Vs),'VData',a12y);
%     set(a13,'XData',a12x,'YData',a12y,'UData',(a13x-a12x),'VData',(a13y-a12y));
%     set(a14,'UData',Vr.*cos(-phiVs),'VData',Vr.*sin(-phiVs),'maxheadsize',maxhead);
    set(arrowscell{1},'UData',Vs,'maxheadsize',maxhead);
    set(arrowscell{2},'XData',Vs,'UData',(a12x-Vs),'VData',a12y);
    set(arrowscell{3},'XData',a12x,'YData',a12y,'UData',(a13x-a12x),'VData',(a13y-a12y));
    set(arrowscell{4},'UData',Vr.*cos(-phiVs),'VData',Vr.*sin(-phiVs),'maxheadsize',maxhead);
    
%% medium line T model
    [ Vrs(2), phiVss(2), Irs(2), Iss(2), phiIss(2) ] = Tlinefunc( S, phir, R, XL, XC, Vs, l );
    Vr = Vrs(2);
    phiVs = phiVss(2);
    Ir = Irs(2);
    Is = Iss(2);
    phiIs = phiIss(2);
    
    sinphiIs = sin(phiIs-phiVs);
    cosphiIs = cos(phiIs-phiVs);
    sinphir = sin(phir-phiVs);
    cosphir = cos(phir-phiVs);
    Isl = 0.5*Is.*l;
    Irl = 0.5*Ir.*l;
    a22x = Vs+XL.*Isl.*sinphiIs;
    a22y = -XL.*Isl.*cosphiIs;
    a23x = a22x -R.*Isl.*cosphiIs;
    a23y = a22y -R.*Isl.*sinphiIs;
    a24x = a23x +XL.*Irl.*sinphir;
    a24y = a23y -XL.*Irl.*cosphir;
    a25x = a24x -R.*Irl.*cosphir;
    a25y = a24y -R.*Irl.*sinphir;

%     set(a21,'xdata',[0 Vs],'ydata',[0 0])
%     set(a22,'xdata',[Vs a22x],'ydata',[0 a22y]);
%     set(a23,'xdata',[a22x a23x],'ydata',[a22y a23y]);
%     set(a24,'xdata',[a23x a24x],'ydata',[a23y a24y]);
%     set(a25,'xdata',[a24x a25x],'ydata',[a24y a25y]);
%     set(a26,'xdata',[0 Vr.*cos(-phiVs)],'ydata',[0 Vr.*sin(-phiVs)]);
   
%     set(a21,'UData',Vs,'maxheadsize',maxhead);
%     set(a22,'XData',Vs,'UData',(a22x-Vs),'VData',a22y);
%     set(a23,'XData',a22x,'YData',a22y,'UData',(a23x-a22x),'VData',(a23y-a22y));
%     set(a24,'XData',a23x,'YData',a23y,'UData',(a24x-a23x),'VData',(a24y-a23y));
%     set(a25,'XData',a24x,'YData',a24y,'UData',(a25x-a24x),'VData',(a25y-a24y));
%     set(a26,'UData',Vr.*cos(-phiVs),'VData',Vr.*sin(-phiVs),'maxheadsize',maxhead);

    set(arrowscell{5},'UData',Vs,'maxheadsize',maxhead);
    set(arrowscell{6},'XData',Vs,'UData',(a22x-Vs),'VData',a22y);
    set(arrowscell{7},'XData',a22x,'YData',a22y,'UData',(a23x-a22x),'VData',(a23y-a22y));
    set(arrowscell{8},'XData',a23x,'YData',a23y,'UData',(a24x-a23x),'VData',(a24y-a23y));
    set(arrowscell{9},'XData',a24x,'YData',a24y,'UData',(a25x-a24x),'VData',(a25y-a24y));
    set(arrowscell{10},'UData',Vr.*cos(-phiVs),'VData',Vr.*sin(-phiVs),'maxheadsize',maxhead);
    
%% medium line pi model
[ Vrs(3), phiVss(3), Irs(3), Iss(3), phiIss(3) ] = pilinefunc( S, phir, Scomplex, R, XL, XC, Vs, l );
    Vr = Vrs(3);
    phiVs = phiVss(3);
    Ir = Irs(3);

Ic2 = Vr.*l./(2*XC);

Isr = Ir.*exp(1i*phir) + 1i*Ic2;
phisr = angle(Isr);
Isr = abs(Isr);

sinphisr = sin(phisr-phiVs);
cosphisr = cos(phisr-phiVs);

a32x = Vs+Isr.*XL.*l.*sinphisr;
a32y = -Isr.*XL.*l.*cosphisr;
a33x = a32x - Isr.*l.*R.*cosphisr;
a33y = a32y - Isr.*l.*R.*sinphisr;

%     set(a31,'xdata',[0 Vs],'ydata',[0 0])
%     set(a32,'xdata',[Vs a32x],'ydata',[0 a32y]);
%     set(a33,'xdata',[a32x a33x],'ydata',[a32y a33y]);
%     set(a34,'xdata',[0 Vr.*cos(-phiVs)],'ydata',[0 Vr.*sin(-phiVs)]);
%     
%     set(a31,'UData',Vs,'maxheadsize',maxhead);
%     set(a32,'XData',Vs,'UData',(a32x-Vs),'VData',a32y);
%     set(a33,'XData',a32x,'YData',a32y,'UData',(a33x-a32x),'VData',(a33y-a32y));
%     set(a34,'UData',Vr.*cos(-phiVs),'VData',Vr.*sin(-phiVs),'maxheadsize',maxhead);
    
    set(arrowscell{11},'UData',Vs,'maxheadsize',maxhead);
    set(arrowscell{12},'XData',Vs,'UData',(a32x-Vs),'VData',a32y);
    set(arrowscell{13},'XData',a32x,'YData',a32y,'UData',(a33x-a32x),'VData',(a33y-a32y));
    set(arrowscell{14},'UData',Vr.*cos(-phiVs),'VData',Vr.*sin(-phiVs),'maxheadsize',maxhead);
    
   
%     Slabelval.String = [num2str(round(S*10)*0.1),' MVA'];
%     Sphilabelval.String = [num2str(0.1*round(-phir.*1800./pi)),' deg'];
%     llabelval.String = [num2str(round(10*l)*0.1),' km'];
%     Vslabel.String = ['Supply:  Vs = ',num2str(Vs),' V'];
%     Vrlabel0.String = ['End (short line):  Vr = ',num2str(0.1*round(10*Vrs(1))),' V'];
%     Vrlabel1.String = ['End (T model):  Vr = ',num2str(0.1*round(10*Vrs(2))),' V'];
%     Vrlabel2.String = ['End (pi model):  Vr = ',num2str(0.1*round(10*Vrs(3))),' V'];
% 
%     drawnow;

    if VsDropdown==1
        xlim([90 160])
        ylim([-45 10])
    elseif VsDropdown==2
        xlim([240 310])
        ylim([-45 10])
    elseif VsDropdown==3
        xlim([370 440])
        ylim([-45 10])
    elseif VsDropdown==4
        xlim([780 850])
        ylim([-45 10])
    end

    set(labelscell{1},'string',['Load Power Magnitude: ',num2str(round(S*10)*0.1),' MVA']);
    set(labelscell{2},'string',['Load Power Angle: ',num2str(0.1*round(-phir.*1800./pi)),' deg']);
    set(labelscell{3},'string',['Line Length: ',num2str(round(10*l)*0.1),' km']);
    set(labelscell{4},'string',['Supply:  Vs = ',num2str(Vs),' kV']);
    set(labelscell{5},'string',['(short line)  Is = ',num2str(0.1*round(10000*Iss(1))),' A, ',num2str(0.1*round(1800*(phiIss(1)-phiVss(1))/pi)),' deg']);
    set(labelscell{6},'string',['(T model)  Is = ',num2str(0.1*round(10000*Iss(2))),' A, ',num2str(0.1*round(1800*(phiIss(2)-phiVss(2))/pi)),' deg']);
    set(labelscell{7},'string',['(pi model)  Is = ',num2str(0.1*round(10000*Iss(3))),' A, ',num2str(0.1*round(1800*(phiIss(3)-phiVss(3))/pi)),' deg']);
    set(labelscell{8},'string',['End (short line):  Vr = ',num2str(0.1*round(10*Vrs(1))),' kV, ',num2str(0.1*round(-1800*phiVss(1)/pi)),' deg  ;   Ir = ',...
                        num2str(0.1*round(10000*Irs(1))),' A, ',num2str(0.1*round(1800*(phir-phiVss(1))/pi)),' deg']);
    set(labelscell{9},'string',['End (T model):  Vr = ',num2str(0.1*round(10*Vrs(2))),' kV, ',num2str(0.1*round(-1800*phiVss(2)/pi)),' deg  ;   Ir = ',...
                        num2str(0.1*round(10000*Irs(2))),' A, ',num2str(0.1*round(1800*(phir-phiVss(2))/pi)),' deg']);
    set(labelscell{10},'string', ['End (pi model):  Vr = ',num2str(0.1*round(10*Vrs(3))),' kV, ',num2str(0.1*round(-1800*phiVss(3)/pi)),' deg  ;   Ir = ',...
                        num2str(0.1*round(10000*Irs(3))),' A, ',num2str(0.1*round(1800*(phir-phiVss(3))/pi)),' deg']);
    set(labelscell{11},'string',['Line R: ',num2str(round(R*1000)*0.001),' ohm/km']);
    set(labelscell{12},'string',['Line XL: ',num2str(round(XL*1000)*0.001),' ohm/km']);
    set(labelscell{13},'string',['Line XC: ',num2str(round(XC*0.01)*0.1),' kilohm.km']);
    
    %     bgcolor = f.Color;
%     
%     Slabel = uicontrol('Parent',f,'Style','text','Position',[0.55*scrsz(3) 0.6*scrsz(4) 280 20],...
%         'string',['Load Power Magnitude: ',num2str(round(S*10)*0.1),' MVA'],'fontsize',12,'BackgroundColor',bgcolor);
%     Sphilabel = uicontrol('Parent',f,'Style','text','Position',[0.55*scrsz(3) 0.475*scrsz(4) 280 20],...
%         'string',['Load Power Angle: ',num2str(0.1*round(-phir.*1800./pi)),' deg'],'fontsize',12,'BackgroundColor',bgcolor);
%     llabel = uicontrol('Parent',f,'Style','text','Position',[0.59*scrsz(3) 0.35*scrsz(4) 170 20],...
%         'string',['Line Length: ',num2str(round(10*l)*0.1),' km'],'fontsize',12,'BackgroundColor',bgcolor);
% 
%     Vslabel = uicontrol('Parent',f,'Style','text','Position',[0.07*scrsz(3) 0.675*scrsz(4) 170 20],'string',...
%         ['Supply:  Vs = ',num2str(Vs),' kV'],'fontsize',12,'BackgroundColor',[1 1 1]);
%     Islabel0 = uicontrol('Parent',f,'Style','text','Position',[0.2*scrsz(3) 0.7*scrsz(4) 270 20],'string',...
%         ['(short line)  Is = ',num2str(0.1*round(10000*Iss(1))),' A, ',num2str(0.1*round(1800*(phiIss(1)-phiVss(1))/pi)),' deg'],...
%         'fontsize',12,'BackgroundColor',[1 1 1]);
%     Islabel1 = uicontrol('Parent',f,'Style','text','Position',[0.2*scrsz(3) 0.675*scrsz(4) 270 20],'string',...
%         ['(T model)  Is = ',num2str(0.1*round(10000*Iss(2))),' A, ',num2str(0.1*round(1800*(phiIss(2)-phiVss(2))/pi)),' deg'],...
%         'fontsize',12,'foregroundcolor','r','BackgroundColor',[1 1 1]);
%     Islabel2 = uicontrol('Parent',f,'Style','text','Position',[0.2*scrsz(3) 0.65*scrsz(4) 270 20],'string',...
%         ['(pi model)  Is = ',num2str(0.1*round(10000*Iss(3))),' A, ',num2str(0.1*round(1800*(phiIss(3)-phiVss(3))/pi)),' deg'],...
%         'fontsize',12,'foregroundcolor','b','BackgroundColor',[1 1 1]);
%     Vrlabel0 = uicontrol('Parent',f,'Style','text','Position',[0.03*scrsz(3) 0.1*scrsz(4) 600 20],'string',...
%         ['End (short line):  Vr = ',num2str(0.1*round(10*Vrs(1))),' kV, ',num2str(0.1*round(-1800*phiVss(1)/pi)),' deg  ;   Ir = ',...
%         num2str(0.1*round(10000*Irs(1))),' A, ',num2str(0.1*round(1800*(phir-phiVss(1))/pi)),' deg'],...
%         'fontsize',12,'BackgroundColor',bgcolor);
%     Vrlabel1 = uicontrol('Parent',f,'Style','text','Position',[0.03*scrsz(3) 0.065*scrsz(4) 600 20],'string',...
%         ['End (T model):  Vr = ',num2str(0.1*round(10*Vrs(2))),' kV, ',num2str(0.1*round(-1800*phiVss(2)/pi)),' deg  ;   Ir = ',...
%         num2str(0.1*round(10000*Irs(2))),' A, ',num2str(0.1*round(1800*(phir-phiVss(2))/pi)),' deg'],...
%         'fontsize',12,'foregroundcolor','r','BackgroundColor',bgcolor);
%     Vrlabel2 = uicontrol('Parent',f,'Style','text','Position',[0.03*scrsz(3) 0.03*scrsz(4) 600 20],'string',...
%         ['End (pi model):  Vr = ',num2str(0.1*round(10*Vrs(3))),' kV, ',num2str(0.1*round(-1800*phiVss(3)/pi)),' deg  ;   Ir = ',...
%         num2str(0.1*round(10000*Irs(3))),' A, ',num2str(0.1*round(1800*(phir-phiVss(3))/pi)),' deg'],...
%         'fontsize',12,'foregroundcolor','b','BackgroundColor',bgcolor);
%     
%     Rlabel = uicontrol('Parent',f,'Style','text','Position',[0.575*scrsz(3) 0.225*scrsz(4) 0.15*scrsz(3) 20],...
%         'string',['Line R: ',num2str(round(R*1000)*0.001),' ohm/km'],'fontsize',12,'BackgroundColor',bgcolor);
%     XLlabel = uicontrol('Parent',f,'Style','text','Position',[0.575*scrsz(3) 0.15*scrsz(4) 0.15*scrsz(3) 20],...
%         'string',['Line XL: ',num2str(round(XL*1000)*0.001),' ohm/km'],'fontsize',12,'BackgroundColor',bgcolor);
%     XClabel = uicontrol('Parent',f,'Style','text','Position',[0.575*scrsz(3) 0.075*scrsz(4) 0.15*scrsz(3) 20],...
%         'string',['Line XC: ',num2str(round(XC*0.01)*0.1),' kilohm.km'],'fontsize',12,'BackgroundColor',bgcolor);

    
end

