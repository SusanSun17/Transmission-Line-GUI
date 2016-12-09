% Thing to do calculations for transmission lines
% Find: Vr
%       Ir
%       phiIr
%       phiVs
%       Is
%       phiIs

close all
clear all
clc

Scomplex = 50+30i;       % MVA, complex load power, Scomplex = Ir* Vr
S = abs(Scomplex);       % magnitude of complex power
phir = -angle(Scomplex); % phase of Ir (-ve phase of power)

Vs = 132;                % kV, supply voltage magnitude

omega = 2*pi*50;         % rad/s, frequency
R = 0.068;               % ohm/km, resistance
C = 24e-9;               % F/km, capacitance
XL = 0.404;              % ohm/km, inductive reactance
L = XL./omega;
XC = 1./(omega*C);       % ohm.km, capacitive reactance
% l = [20:5:200];                  % km, line length
l = 70;

Vrs = zeros(length(l),3);
phiVss = zeros(length(l),3);
Irs = zeros(length(l),3);
Iss = zeros(length(l),3);
phiIss = zeros(length(l),3);

scrsz = get(groot,'screensize');
f = figure('position',[0.1*scrsz(3) 0.1*scrsz(4) 0.8*scrsz(3) 0.8*scrsz(4)]);
ax = axes('Parent',f,'position',[0.05 0.25 0.55 0.675]);
axis equal
hold on

%% short line
[ Vrs(1), phiVss(1), Irs(1) ] = shortlinefunc( S, phir, Scomplex, R, XL, Vs, l );
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

% a11 = plot([0 Vs],[0 0],'k');
% a12 = plot([Vs a12x],[0 a12y],'k');
% a13 = plot([a12x a13x],[a12y a13y],'k');
% a14 = plot([0 Vr.*cos(-phiVs)],[0 Vr.*sin(-phiVs)],'k','linewidth',2);

a11 = quiver(0, 0, Vs, 0,'k','autoscale','off','maxheadsize',0.02);
a12 = quiver(Vs, 0, (a12x-Vs), a12y,'k','autoscale','off');
a13 = quiver(a12x, a12y, (a13x-a12x), (a13y-a12y),'k','autoscale','off');
a14 = quiver(0,0,Vr.*cos(-phiVs),Vr.*sin(-phiVs),'k','linewidth',2,'autoscale','off','maxheadsize',0.02);
axis equal
grid on
 
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

% a21 = plot([0 Vs],[0 0],'r');
% a22 = plot([Vs a22x],[0 a22y],'r');
% a23 = plot([a22x a23x],[a22y a23y],'r');
% a24 = plot([a23x a24x],[a23y a24y],'r');
% a25 = plot([a24x a25x],[a24y a25y],'r');
% a26 = plot([0 Vr.*cos(-phiVs)],[0 Vr.*sin(-phiVs)],'r','linewidth',1.5);

a21 = quiver(0, 0, Vs, 0,'r','autoscale','off','maxheadsize',0.02);
a22 = quiver(Vs, 0, (a22x-Vs), a22y,'r','autoscale','off');
a23 = quiver(a22x, a22y, (a23x-a22x), (a23y-a22y),'r','autoscale','off');
a24 = quiver(a23x, a23y, (a24x-a23x), (a24y-a23y),'r','autoscale','off');
a25 = quiver(a24x, a24y, (a25x-a24x), (a25y-a24y),'r','autoscale','off');
a26 = quiver(0, 0, Vr.*cos(-phiVs), Vr.*sin(-phiVs),'r','linewidth',1.5,'autoscale','off','maxheadsize',0.02);

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

% a31 = plot([0 Vs],[0 0],'b');
% a32 = plot([Vs a32x],[0 a32y],'b');
% a33 = plot([a32x a33x],[a32y a33y],'b');
% a34 = plot([0 Vr.*cos(-phiVs)],[0 Vr.*sin(-phiVs)],'b','linewidth',1.5);

a31 = quiver(0, 0, Vs, 0,'b','autoscale','off','maxheadsize',0.02);
a32 = quiver(Vs, 0, (a32x-Vs), a32y,'b','autoscale','off');
a33 = quiver(a32x, a32y, (a33x-a32x), (a33y-a32y),'b','autoscale','off');
a34 = quiver(0,0,Vr.*cos(-phiVs),Vr.*sin(-phiVs),'b','linewidth',1.5,'autoscale','off','maxheadsize',0.02);


% xVmin = min([Vrs.*cos(phiVss) Vs]);
% xVmin = floor(xVmin*0.2)*5;
% xVmax = max([Vrs.*cos(phiVss) Vs]);
% xVmax = ceil(xVmax*0.2)*5;
% yVmin = min([2 -Vrs.*sin(phiVss)]);
% yVmin = floor(yVmin*0.5)*2;
% yVmax = max([2 -Vrs.*sin(phiVss)]);
% yVmax = ceil(yVmax*0.5)*2;

axis equal
xlim([90 160])
ylim([-45 10])

title('Voltage Change across a Medium-Length Transmission Line','fontsize',12)
xlabel('Real(Voltage) / kV')
ylabel('Imaginary(Voltage) / kV')
% legend([a11 a21 a31],{'Short Line','T model','Pi Model'})

bgcolor = f.Color;

% supply voltage dropdown
VsDropdown = uicontrol('Parent',f,'Style','popup','String',{'132 kV','275 kV','400 kV','800 kV'},...
    'fontsize',12,'Position',[0.6*scrsz(3) 0.725*scrsz(4) 135 20]);
VsDroplabel = uicontrol('Parent',f,'Style','text','Position',[0.6*scrsz(3) 0.75*scrsz(4) 135 20],...
    'string','Supply Voltage','fontsize',12,'BackgroundColor',bgcolor);
% % load magnitude slider
Sslider = uicontrol('Parent',f,'Style','Slider','Position',[0.55*scrsz(3) 0.625*scrsz(4) 0.2*scrsz(3) 20],'value',S,'min',0,'max',500);
Slabelmin = uicontrol('Parent',f,'Style','text','Position',[0.55*scrsz(3) 0.65*scrsz(4) 20 20],'string','0','fontsize',12,'BackgroundColor',bgcolor);
Slabelmax = uicontrol('Parent',f,'Style','text','Position',[0.725*scrsz(3) 0.65*scrsz(4) 40 20],'string','500','fontsize',12,'BackgroundColor',bgcolor);
Slabel = uicontrol('Parent',f,'Style','text','Position',[0.55*scrsz(3) 0.6*scrsz(4) 280 20],...
    'string',['Load Power Magnitude: ',num2str(round(S*10)*0.1),' MVA'],'fontsize',12,'BackgroundColor',bgcolor);
% Slabelval = uicontrol('Parent',f,'Style','text','Position',[0.6*scrsz(3) 0.625*scrsz(4) 150 20],...
%     'string',[num2str(round(S*10)*0.1),' MVA'],'fontsize',12,'BackgroundColor',bgcolor);
% % load angle slider
Sphislider = uicontrol('Parent',f,'Style','Slider','Position',[0.55*scrsz(3) 0.5*scrsz(4) 0.2*scrsz(3) 20],'value',-phir.*180./pi,'min',-90,'max',90);
Sphilabelmin = uicontrol('Parent',f,'Style','text','Position',[0.5*scrsz(3) 0.525*scrsz(4) 150 20],'string','-90 (capacitive)','fontsize',12,'BackgroundColor',bgcolor);
Sphilabelmax = uicontrol('Parent',f,'Style','text','Position',[0.685*scrsz(3) 0.525*scrsz(4) 150 20],'string','+90 (inductive)','fontsize',12,'BackgroundColor',bgcolor);
Sphilabel = uicontrol('Parent',f,'Style','text','Position',[0.55*scrsz(3) 0.475*scrsz(4) 280 20],...
    'string',['Load Power Angle: ',num2str(0.1*round(-phir.*1800./pi)),' deg'],'fontsize',12,'BackgroundColor',bgcolor);
% Sphilabelval = uicontrol('Parent',f,'Style','text','Position',[0.625*scrsz(3) 0.475*scrsz(4) 80 20],...
%     'string',[num2str(0.1*round(-phir.*1800./pi)),' deg'],'fontsize',12,'BackgroundColor',bgcolor);

% % line length slider
lslider = uicontrol('Parent',f,'Style','Slider','Position',[0.55*scrsz(3) 0.375*scrsz(4) 0.2*scrsz(3) 20],'value',l,'min',0,'max',200);
llabelmin = uicontrol('Parent',f,'Style','text','Position',[0.55*scrsz(3) 0.4*scrsz(4) 20 20],'string','0','fontsize',12,'BackgroundColor',bgcolor);
llabelmax = uicontrol('Parent',f,'Style','text','Position',[0.725*scrsz(3) 0.4*scrsz(4) 40 20],'string','200','fontsize',12,'BackgroundColor',bgcolor);
llabel = uicontrol('Parent',f,'Style','text','Position',[0.59*scrsz(3) 0.35*scrsz(4) 170 20],...
    'string',['Line Length: ',num2str(round(10*l)*0.1),' km'],'fontsize',12,'BackgroundColor',bgcolor);
% llabelval = uicontrol('Parent',f,'Style','text','Position',[0.625*scrsz(3) 0.325*scrsz(4) 80 20],...
%     'string',[num2str(round(10*l)*0.1),' km'],'fontsize',12,'BackgroundColor',bgcolor);

% voltage labels
    Vslabel = uicontrol('Parent',f,'Style','text','Position',[0.07*scrsz(3) 0.675*scrsz(4) 170 20],'string',...
        ['Supply:  Vs = ',num2str(Vs),' kV'],'fontsize',12,'BackgroundColor',[1 1 1]);
    Islabel0 = uicontrol('Parent',f,'Style','text','Position',[0.2*scrsz(3) 0.7*scrsz(4) 270 20],'string',...
        ['(short line)  Is = ',num2str(0.1*round(10000*Iss(1))),' A, ',num2str(0.1*round(1800*(phiIss(1)-phiVss(1))/pi)),' deg'],...
        'fontsize',12,'BackgroundColor',[1 1 1]);
    Islabel1 = uicontrol('Parent',f,'Style','text','Position',[0.2*scrsz(3) 0.675*scrsz(4) 270 20],'string',...
        ['(T model)  Is = ',num2str(0.1*round(10000*Iss(2))),' A, ',num2str(0.1*round(1800*(phiIss(2)-phiVss(2))/pi)),' deg'],...
        'fontsize',12,'foregroundcolor','r','BackgroundColor',[1 1 1]);
    Islabel2 = uicontrol('Parent',f,'Style','text','Position',[0.2*scrsz(3) 0.65*scrsz(4) 270 20],'string',...
        ['(pi model)  Is = ',num2str(0.1*round(10000*Iss(3))),' A, ',num2str(0.1*round(1800*(phiIss(3)-phiVss(3))/pi)),' deg'],...
        'fontsize',12,'foregroundcolor','b','BackgroundColor',[1 1 1]);
    Vrlabel0 = uicontrol('Parent',f,'Style','text','Position',[0.03*scrsz(3) 0.1*scrsz(4) 600 20],'string',...
        ['End (short line):  Vr = ',num2str(0.1*round(10*Vrs(1))),' kV, ',num2str(0.1*round(-1800*phiVss(1)/pi)),' deg  ;   Ir = ',...
        num2str(0.1*round(10000*Irs(1))),' A, ',num2str(0.1*round(1800*(phir-phiVss(1))/pi)),' deg'],...
        'fontsize',12,'BackgroundColor',bgcolor);
    Vrlabel1 = uicontrol('Parent',f,'Style','text','Position',[0.03*scrsz(3) 0.065*scrsz(4) 600 20],'string',...
        ['End (T model):  Vr = ',num2str(0.1*round(10*Vrs(2))),' kV, ',num2str(0.1*round(-1800*phiVss(2)/pi)),' deg  ;   Ir = ',...
        num2str(0.1*round(10000*Irs(2))),' A, ',num2str(0.1*round(1800*(phir-phiVss(2))/pi)),' deg'],...
        'fontsize',12,'foregroundcolor','r','BackgroundColor',bgcolor);
    Vrlabel2 = uicontrol('Parent',f,'Style','text','Position',[0.03*scrsz(3) 0.03*scrsz(4) 600 20],'string',...
        ['End (pi model):  Vr = ',num2str(0.1*round(10*Vrs(3))),' kV, ',num2str(0.1*round(-1800*phiVss(3)/pi)),' deg  ;   Ir = ',...
        num2str(0.1*round(10000*Irs(3))),' A, ',num2str(0.1*round(1800*(phir-phiVss(3))/pi)),' deg'],...
        'fontsize',12,'foregroundcolor','b','BackgroundColor',bgcolor);

% % line resistance slider
Rslider = uicontrol('Parent',f,'Style','Slider','Position',[0.575*scrsz(3) 0.25*scrsz(4) 0.15*scrsz(3) 20],'value',R,'min',0,'max',1);
Rlabelmin = uicontrol('Parent',f,'Style','text','Position',[0.56*scrsz(3) 0.25*scrsz(4) 20 20],'string','0','fontsize',12,'BackgroundColor',bgcolor);
Rlabelmax = uicontrol('Parent',f,'Style','text','Position',[0.73*scrsz(3) 0.25*scrsz(4) 30 20],'string','1.0','fontsize',12,'BackgroundColor',bgcolor);
Rlabel = uicontrol('Parent',f,'Style','text','Position',[0.575*scrsz(3) 0.225*scrsz(4) 0.15*scrsz(3) 20],...
    'string',['Line R: ',num2str(round(R*1000)*0.001),' ohm/km'],'fontsize',12,'BackgroundColor',bgcolor);
% Rlabelval = uicontrol('Parent',f,'Style','text','Position',[0.6*scrsz(3) 0.225*scrsz(4) 150 20],...
%     'string',[num2str(round(R*1000)*0.001),' ohm/km'],'fontsize',12,'BackgroundColor',bgcolor);
% % line reactance inductive slider
XLslider = uicontrol('Parent',f,'Style','Slider','Position',[0.575*scrsz(3) 0.175*scrsz(4) 0.15*scrsz(3) 20],'value',XL,'min',0,'max',1);
XLlabelmin = uicontrol('Parent',f,'Style','text','Position',[0.56*scrsz(3) 0.175*scrsz(4) 20 20],'string','0','fontsize',12,'BackgroundColor',bgcolor);
XLlabelmax = uicontrol('Parent',f,'Style','text','Position',[0.73*scrsz(3) 0.175*scrsz(4) 30 20],'string','1.0','fontsize',12,'BackgroundColor',bgcolor);
XLlabel = uicontrol('Parent',f,'Style','text','Position',[0.575*scrsz(3) 0.15*scrsz(4) 0.15*scrsz(3) 20],...
    'string',['Line XL: ',num2str(round(XL*1000)*0.001),' ohm/km'],'fontsize',12,'BackgroundColor',bgcolor);
% % line reactance inductive slider
XCslider = uicontrol('Parent',f,'Style','Slider','Position',[0.575*scrsz(3) 0.1*scrsz(4) 0.15*scrsz(3) 20],'value',XC,'min',0,'max',2e5);
XClabelmin = uicontrol('Parent',f,'Style','text','Position',[0.56*scrsz(3) 0.1*scrsz(4) 20 20],'string','0','fontsize',12,'BackgroundColor',bgcolor);
XClabelmax = uicontrol('Parent',f,'Style','text','Position',[0.73*scrsz(3) 0.1*scrsz(4) 30 20],'string','200','fontsize',12,'BackgroundColor',bgcolor);
XClabel = uicontrol('Parent',f,'Style','text','Position',[0.575*scrsz(3) 0.075*scrsz(4) 0.15*scrsz(3) 20],...
    'string',['Line XC: ',num2str(round(XC*0.01)*0.1),' kilohm.km'],'fontsize',12,'BackgroundColor',bgcolor);

arrowscell = {a11;a12;a13;a14;a21;a22;a23;a24;a25;a26;a31;a32;a33;a34};
labelscell = {Slabel;Sphilabel;llabel;Vslabel;Islabel0;Islabel1;Islabel2;...
               Vrlabel0;Vrlabel1;Vrlabel2;Rlabel;XLlabel;XClabel};
    
% Sslider.Callback = @(esS,edS) updatel(f, esS.Value, -Sphislider.Value.*pi./180,...
%     Rslider.Value, XLslider.Value, XCslider.Value, VsDropdown.Value, lslider.Value,...
%     a11, a12, a13, a14, a21, a22, a23, a24, a25, a26, a31, a32, a33, a34 ) ;
Sslider.Callback = @(esS,edS) updatel(f, esS.Value, -Sphislider.Value.*pi./180,...
    Rslider.Value, XLslider.Value, XCslider.Value, VsDropdown.Value, lslider.Value,...
    arrowscell, labelscell ) ;
Sphislider.Callback = @(esp,edp) updatel(f, Sslider.Value, -esp.Value.*pi./180,...
    Rslider.Value, XLslider.Value, XCslider.Value, VsDropdown.Value, lslider.Value,...
    arrowscell, labelscell ) ;
lslider.Callback = @(esl,edl) updatel(f, Sslider.Value, -Sphislider.Value.*pi./180,...
    Rslider.Value, XLslider.Value, XCslider.Value, VsDropdown.Value, esl.Value,...
    arrowscell, labelscell ) ;
VsDropdown.Callback = @(esV,edV) updatel(f, Sslider.Value, -Sphislider.Value.*pi./180,...
    Rslider.Value, XLslider.Value, XCslider.Value, esV.Value, lslider.Value,...
    arrowscell, labelscell ) ;
Rslider.Callback = @(esR,edR) updatel(f, Sslider.Value, -Sphislider.Value.*pi./180,...
    esR.Value, XLslider.Value, XCslider.Value, VsDropdown.Value, lslider.Value,...
    arrowscell, labelscell ) ;
XLslider.Callback = @(esXL,edXL) updatel(f, Sslider.Value, -Sphislider.Value.*pi./180,...
    Rslider.Value, esXL.Value, XCslider.Value, VsDropdown.Value, lslider.Value,...
    arrowscell, labelscell ) ;
XCslider.Callback = @(esXC,edXC) updatel(f, Sslider.Value, -Sphislider.Value.*pi./180,...
    Rslider.Value, XLslider.Value, esXC.Value, VsDropdown.Value, lslider.Value,...
    arrowscell, labelscell ) ;

%%







