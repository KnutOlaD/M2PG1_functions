% This code is provided by Mario Veloso who got it from Ira Leifer. Not my code.
%
% function                         buoyvel.m               
%
% Calculates the Buoyancy Velocity for a clean (dirty=0), dirty, or partially
% surfactant covered bubble of radius r (r<1cm) at temperature T (C)
%    Rise velocity for a single radius now works.
%
%    If no arguments in, produces a plot of rise velocities for all Models at 20C,
%       and output of the handles to the MainPlot;Legend;MinorTicks.
%
%     RisePar = Leifer=>  Clift et al., 1978  Fig 7.3  loads from bubrise.dat
%                       T dependence for Re>50 (r>320)
%                                    from polyfit of Re to CD  in Eqn 5-31
%                       T dependence for Re>1 (r>60) & Re<50 from Eqn 5-28
%                       T dependence for Re<1 (r>60)         from Eqn 3-14 (Stokes)
%
%     RisePar = 'Woolf93'            => Only Dirty
%     RisePar = 'Woolf & Thorpe91'   => 
%     RisePar = 'Memery & Merlivat'  => Only Clean
%     RisePar = 'Leifer'
%     RisePar = 'Mendelson'          => Only Clean; only for r>rc.
%     RisePar = 'Leifer&Patro'                 Continental Shelf Research, Accepted
%     RisePar = 'Fan&Tsuchiya'       => Clean and dirty
%     RisePar = 'Off'   yields Rise Velocities = 0 for all radii
%     rad must be a smooth distribution
%
%     Sal is in ppt (default is 35)
%
%    Varying surfactant Transition set at onset of oscillations for
%    Surf=0.5, else at RP*(1+2*(Surf-0.5)) where RP is the onset radius
%
%                Ira Leifer                  4/22/07  MatLab 5.1 
%
%  function Wbuoy = buoyvel(rad, T [,Dirty, RisePar ,Sal]);
   function Wbuoy = buoyvel(rad, T,  Dirty, RisePar ,Sal )

if nargin==0, rad=1e-4*(1:3000);PlotMode=1;  Dirty=0.5; Sal=35;
              radius = (1:3000);T=20; RisePar='Leifer';
else          PlotMode=0;
end

if nargin==2, Dirty=1; RisePar = 'Leifer'; end
if nargin==3
   if ischar(Dirty), RisePar='Leifer'; Dirty=1; end
else
   if ischar(RisePar)~=1, error('Not a Valid Model Type'); end
end
if nargin==4, Sal = 35; end;

[i j]=size(rad); Wbuoy = zeros(i,j);
if strcmp(caps(RisePar),'OFF'); return; end

Visc = viscosity(T,Sal);
DWater = physcon('Density', 'Water', T, Sal);
   
Wbuoy = zeros(size(rad)); WbuoyDirty = Wbuoy; WbuoyStoke = Wbuoy;
                          WbuoyClean = Wbuoy;
if strcmp(caps(RisePar),'LEIFER') || strcmp(caps(RisePar),'LEIFER&PATRO'),
      DataRad   = [25:25:175 200:50:950 1000:100:5900 6000:200:10000]*1e-4;
      DataDirty = [ 0.18   0.65    1.2    1.7    2.2    2.8   3.5     4.0 ...  %  25--200
                    5.40   6.40    7.5    8.5    9.2                      ...  % 250--450
                   10.05  10.90   11.40  12.10  12.8   13.30  13.90  14.30  14.75  15.10 ...  % 500-- 950
                   15.45  15.95   16.30  16.50  16.70  16.80  16.90  17.00  17.10  17.20 ...  %1000--1900
                   17.35  17.45   17.50  17.60  17.65  17.70  17.85  17.95  18.05  18.10 ...  %2000--2900
                   18.25  18.35   18.45  18.55  18.70  18.80  18.95  19.10  19.25  19.35 ...  %3000--3900
                   19.50  19.65   19.75  19.85  20.05  20.20  20.30  20.45  20.60  20.75 ...  %4000--4900
                   20.85  21.00   21.15  21.30  21.45  21.60  21.70  21.85  22.00  22.20 ...  %5000--5900
                   22.35  22.60   22.90  23.20  23.50  23.85  24.15  24.40  24.80  25.30 ...  %6000--7800
                   25.65  26.05   26.45  26.80  27.20  27.60  28.00  28.25  28.55  28.90 ...  %8000--9800
                   29.25];                                                                    %10,000
      DataClean = [ 0.18   0.65    1.2    1.7    2.2    2.8    3.5    4.0 ...  %  25--200
                    5.40   7.10    8.9   10.7   12.3                      ...  % 250--450
                   14.70  20.40   25.00  31.10  33.80  33.20  32.80  31.80  31.00  30.40 ...  % 500-- 950
                   29.85  28.65   27.90  27.10  26.50  26.05  25.60  25.20  24.95  24.60 ...  %1000--1900
                   24.35  24.10   23.95  23.90  23.85  23.80  23.70  23.60  23.50  23.45 ...  %2000--2900
                   23.45  23.45   23.40  23.40  23.40  23.40  23.40  23.45  23.55  23.60 ...  %3000--3900
                   23.65  23.75   23.80  23.85  23.95  24.05  24.20  24.35  24.45  24.60 ...  %4000--4900
                   24.70  24.80   24.85  24.95  25.00  25.20  25.30  25.40  25.50  25.60 ...  %5000--5900
                   25.80  26.00   26.20  26.40  26.70  26.85  27.10  27.25  27.45  27.70 ...  %6000--7800
                   27.95  28.20   28.35  28.60  28.85  29.20  29.45  29.60  29.80  30.00 ...  %8000--8800
                   30.20];                                                                    %10,000
end
                          
%%%%%%%%%%%%%%%%%%% Stokes Rise Velocity %%%%%%%%%%%%%%%%%%%%%%%%% 
h = find(rad>0); Dbub=28/22400; % g/mole / cm3/mole - for air at STP!!!
WbuoyStoke(h) = (2/9)*(980/Visc)* (DWater-Dbub) * rad(h).*rad(h);

if PlotMode==1, RisePar='Stokes'; end;
if strcmp(caps(RisePar),'STOKES'),
   WbuoyClean = WbuoyStoke; WbuoyDirty = WbuoyStoke;  
   if PlotMode==1,
      clf; plot(radius,WbuoyStoke,'r'); axis([0 3000 0 40]); hold on;
      Hndl = gca;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Memery & Merlivat %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Woolf & Thorpe91 %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if PlotMode==1, RisePar='Woolf & Thorpe91'; end;
if strcmp(caps(RisePar),'WOOLF & THORPE91') || strcmp(caps(RisePar),'MEMERY & MERLIVAT'),

   Wbuoy2 = Wbuoy; Wbuoy3 = Wbuoy; y = Wbuoy;
   Chi(h) =  (980*rad(h).^3 / Visc.^2);

   %%%%%%%% Clean
   Wbuoy3(h) = 4.5 * WbuoyStoke(h) ./ (18*(1-2*ones(size(h))./(1+sqrt(1+.091*Chi(h)))));
   Wbuoy80  =  2*(980/Visc)*(80e-4.^2) /9;
   Var = (980/Visc.^2)*(150e-4).^3;
   Wbuoy150 = (980/Visc)*(150e-4).^2 / (18*(1-2/(1+sqrt(1+.091*Var))));
   slope    = (Wbuoy150 - Wbuoy80) / .007; 
   Wbuoy2(h)   =  Wbuoy80 + slope * (rad(h) - .008);       % Inter Radii
   i = find(rad<=80e-4); j = find(rad >80e-4); 
   if length(i)>5 && length(j)>5,
      WbuoyClean = gausjoin(rad, WbuoyStoke , Wbuoy2,  80e-4, 25e-4);
   else
      WbuoyClean = [ WbuoyStoke(i)  Wbuoy2(j) ];
   end
   
   i = find(rad<=150e-4); j = find(rad >150e-4);
   if length(i)>5 && length(j)>5,
      WbuoyClean = gausjoin(rad, WbuoyClean , Wbuoy3, 150e-4,50e-4);
   else
      WbuoyClean = [ WbuoyClean(i) Wbuoy3(j) ];
   end
   
   %%%%%%%% Surfactant Covered
   y(h) = 10.82*ones(size(h)) ./ Chi(h);
   Wbuoy2(h) = WbuoyStoke(h) .* ( sqrt(y(h).^2 + 2*y(h)) - y(h));
   
   i = find(rad<=80e-4); j = find(rad >80e-4); 
   if length(i)>5 && length(j)>5,
      WbuoyDirty = gausjoin(rad, WbuoyStoke, Wbuoy2, 80e-4, 25e-4);
   else
      WbuoyDirty = [ WbuoyStoke(i)  Wbuoy2(j) ];
   end
   if strcmp(caps(RisePar),'MEMERY & MERLIVAT'), WbuoyDirty=WbuoyClean; end;
   
   if PlotMode==1,
      plot(radius,WbuoyClean,'--y');
%     plot(radius,WbuoyDirty,'r'  );  % Woolfe 1991 Dirty superceded by 1993
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Woolf93 %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if PlotMode==1, RisePar='Woolf93'; end
if  strcmp(caps(RisePar),'WOOLF93'),
   WbuoyClean = 0.172*(980^0.76)*(Visc^(-0.56))*(rad(h).^1.28);
   WbuoyDirty = WbuoyClean;
   %%%%%%%%%%%%%%      Vt max = 25 cm/s (Woolf, 1993) %%%%%%%%%%%
   i = find(WbuoyDirty<=25); j = find(WbuoyDirty>25);
   if length(i)>5 && length(j)>5
      r = find(WbuoyDirty<25, 1, 'last' );       
      WbuoyDirty = gausjoin(rad, WbuoyDirty, 25*ones(size(rad)), rad(r), 100e-4);
   else
      i=find(WbuoyDirty>25); WbuoyDirty(i) = 25*ones(size(i));
   end;
   if PlotMode==1, plot(radius,WbuoyDirty,':r'); end; 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Leifer  %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Atkinson %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if PlotMode==1, RisePar = 'Leifer'; end;
if      strcmp(caps(RisePar),'LEIFER'),
   Tdep = zeros(size(Wbuoy));
   i = find(rad<=max(DataRad));
   WbuoyClean(i) = interp1(DataRad,DataClean,rad(i),'spline'); 
   WbuoyDirty(i) = interp1(DataRad,DataDirty,rad(i),'spline');
   i = find(rad>max(DataRad));
   WbuoyClean(i) = DataClean(length(DataRad)) * ones(size(i));  % VB Constant for r>Data
   WbuoyDirty(i) = DataDirty(length(DataRad)) * ones(size(i));
   
   i = find(rad<=320e-4); j = find(rad >320e-4);                % ReNum of 50
   if length(i)>5 && length(j)>5,
      Tdep = gausjoin(rad, -0.64*ones(size(Tdep)), -0.54*ones(size(Tdep)), 320e-4, 150e-4);
   else
      Tdep = [-0.64*ones(size(i))  -0.54*ones(size(j))];
   end
   i = find(rad<=60e-4); j = find(rad >60e-4);                  % ReNum of 1
   if length(i)>5 && length(j)>5,
      Tdep = gausjoin(rad, -0.5*ones(size(Tdep)) , Tdep, 60e-4, 20e-4);
   else
      Tdep = [-0.5*ones(size(i))  Tdep(j)];
   end
    
     i = ones(size(Tdep));                                       % Introduce T-Dependence
   Var = (Visc*i).^Tdep  ./ (0.00998*i).^Tdep;
   WbuoyClean = WbuoyClean.* Var;               
   WbuoyDirty = WbuoyDirty.* Var;
   
   i = find(rad<=80e-4); j = find(rad >80e-4); 
   if length(i)>5 && length(j)>5,
       WbuoyClean = gausjoin(rad, WbuoyStoke , WbuoyClean, 80e-4, 25e-4);
       WbuoyDirty = gausjoin(rad, WbuoyStoke , WbuoyDirty, 80e-4, 25e-4);
   else
       WbuoyClean = [WbuoyStoke(i) WbuoyClean(j)];
       WbuoyDirty = [WbuoyStoke(i) WbuoyDirty(j)];
   end
   if PlotMode
       LHndl=plot(radius,WbuoyDirty,'b' ); set(LHndl,'LineWidth',2); 
       LHndl=plot(radius,WbuoyClean,':r'); set(LHndl,'LineWidth',2);
   end
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Leifer & Patro  %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if PlotMode, RisePar = 'Leifer&Patro'; end
if      strcmp(caps(RisePar),'LEIFER&PATRO')
   Visc = viscosity(T,Sal);
  %Tdep = zeros(size(rad)); Temp = Tdep; Wbuoy2 = Wbuoy; 

   i = find(rad<=max(DataRad));
   WbuoyDirty(i) = interp1(DataRad,DataDirty,rad(i),'spline');

   WbuoyDirty = gausjoin(rad, WbuoyStoke, WbuoyDirty,  50e-4, 20e-4);

   WbuoyDirty2 = polyval([17.8901 11.4222],rad);
   WbuoyDirty = gausjoin(rad, WbuoyDirty, WbuoyDirty2,8000e-4,500e-4);

   TDep = (10.73*rad-1.055); i = find(TDep< -1); TDep(i)=-1*ones(size(i));
 % plot(rad*1e4,TDep,[60 270 580 680 775],[-1 -.69 -.52 -.34 -.17],'o');axis([0 1e3 -1.5 1]);
 % Leifer and Patro 2001
   WbuoyClean  = 326.7 * 1.0000.*(rad.^2.000).*Visc.^TDep;   %    1 -  60 µm
   WbuoyClean2 = 326.7 * 0.7500.*(rad.^2.000).*Visc.^TDep;   %   60 - 275 µm
   WbuoyClean3 = 326.7 * 50.660.*(rad.^3.171).*Visc.^TDep;   %  275 - 475 µm
   WbuoyClean4 = 326.7 * 75148 .*(rad.^5.570).*Visc.^TDep;   %  475 - 665 µm
   WbuoyClean5 = 326.7 * 1277.0.*(rad.^4.068).*Visc.^TDep;   %  665 - 785 µm
   WbuoyClean6 = 326.7 * 5042.1.*(rad.^4.606).*Visc.^TDep;   %  785 - 890 µm
%       Patro et al., 2001
%  WbuoyClean  = 326.7 * 1.0000.*(rad.^2.000).*Visc.^(-1.00);   %    1 -  60 µm
%  WbuoyClean2 = 326.7 * 0.1397.*(rad.^1.372).*Visc.^(-0.64);   %   60 - 500 µm
%  WbuoyClean3 = 326.7 * 11.713.*(rad.^2.851).*Visc.^(-0.64);   %  500 - 660 µm
%  WbuoyClean4 = 326.7 * 0.1563.*(rad.^1.263).*Visc.^(-0.64);   %  660 - 700 µm
%  WbuoyClean5 = 326.7 * 0.0211.*(rad.^0.511).*Visc.^(-0.64);   %  700 - 850 µm

%  WbuoyClean3 = 326.7 * 22.677.*(rad.^2.88 ).*Visc.^(-0.52);   %  500 - 660 µm
%  WbuoyClean4 = 326.7 * 0.6333.*(rad.^1.263).*Visc.^(-0.34);   %  660 - 700 µm
%  WbuoyClean5 = 326.7 * 1.1413.*(rad.^0.320).*Visc.^(-0.17);   %  700 - 850 µm

   WbuoyClean = gausjoin(rad, WbuoyClean, WbuoyClean2,  60e-4, 50e-4);
   WbuoyClean = gausjoin(rad, WbuoyClean, WbuoyClean3, 275e-4, 50e-4);
   WbuoyClean = gausjoin(rad, WbuoyClean, WbuoyClean4, 475e-4, 50e-4); 
   WbuoyClean = gausjoin(rad, WbuoyClean, WbuoyClean5, 665e-4, 50e-4);
   WbuoyClean = gausjoin(rad, WbuoyClean, WbuoyClean6, 785e-4, 50e-4);
  %figure(3);plot(rad,WbuoyClean);pause;
  
%  i = find(WbuoyClean>35); WbuoyClean(i)=35*ones(size(i));

   RP = (T<0)*962 + (T>=0 & T<=35)*(962 -10.2*T) + (T>35)*590;
 % RP = RP*(1+2*(Dirty-0.5));
   
   i = find(rad<RP*1e-4);  j = find(rad>=RP*1e-4);
   if ~isempty(j)
      k =      -4.792e-4* (rad-0.0584) .^ -0.815;
      h = 22.16 + 0.733 * (rad-0.0584) .^ -0.849;
          
      Wbuoy2 = ones(size(i))*h(j(1)) .* exp(k(j(1)).*T);
      
      Wbuoy2 = [Wbuoy2 h(j).*exp(k(j).*T)];
      
      if T>30 && find(rad>.05 & rad < .09)
         disp('Parameterization does not work properly for bubbles 0.05<r<0.09cm, T>35C');
         i = find(WbuoyClean>max(Wbuoy2(j))); WbuoyClean(i) = ones(size(i))*max(Wbuoy2(j));
         i = 1:find(Wbuoy2==max(Wbuoy2(j))); Wbuoy2(i)     = ones(size(i))*max(Wbuoy2(j));
         figure(1);   plot(rad,WbuoyClean,rad,Wbuoy2);pause
      end;
      WbuoyClean = gausjoin(rad, WbuoyClean, Wbuoy2, RP*1e-4, 100e-4);
     %figure(3);   plot(rad,WbuoyClean,'r',rad,Wbuoy2,'k'); pause
   end

   i = find(rad>0.4);
   if length(i)>1,
	   Temp = polyval([11.0478 19.1510], rad);
       WbuoyClean = gausjoin(rad, WbuoyClean, Temp, rad(min(i)), 200e-4);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if      strcmp(caps(RisePar),'MENDELSON'),
	STens = 75.644 -0.1702*T +0.00214*T.*T;           % Surface Tension  (dyn/cm2)
    WbuoyClean = (1.07*STens/DWater./rad + 1.01*980.*rad).^0.5;
end;


%%%%%%%%%%%%%%%%%%%%% Fan & Tsuchiya, 1990 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%function Wbuoy = buoyvel(rad, T,  Dirty, RisePar ,Sal )
%auth: Knut Ola Dølven, https://github.com/KnutOlaD
%Model is described in Leifer & Patro 2002. Original source is 
%Fan, L. and Tsuchiya K., (1990): "Bubble Wake Dynamics in Liquids and LiquidSolid 
%Suspensions", Butterworth-Heinemann, doi: https://doi.org/10.1016/C2009-0-24002-5

if strcmp(caps(RisePar),'FAN&TSUCHIYA')
    
    %Fan & Tsuchiya has tuning input [c,d] which determines if the
    %model is for freshwater/saltwater and clean/dirty bubbles. The c parameter is set
    %to 1.4 for seawater and 1.2 for pure freshwater. The d parameter is set to 1.6 for completely
    %clean bubbles and 0.8 for dirty bubbles .Valid intermediate
    %values for these parameters as well, but see Fan & Tsuchiyua, 1990 and Leifer & Patro, 2002
    %for details.
    %For c, I'm using a linear gradient between S=0 (freshwater, i.e. c=1.2)
    %and S=35 (seawater, i.e. 1.4)
    %For dirty/clean, it is my impression that this is calculated for pure
    %dirty (WbuoyDirty) and pure clean (WbouyClean) bubbles in all instances and then a transition
    %function is applied afterwards
    
    %Calculate input parameters
    %Visc = viscosity(T,Sal); % Viscosity of the seawater [Pa s]
    density = density(T, Sal); % kg/m^3 density of seawater [kg/m^3]
    % Would be better to use gsw, but using this function to avoid dependencies.
    g_acc = 9.81; % m/s^2 Gravitational acceleration [m/s^2]
    surf_tens = 0.0728; % N/m Surface tension of water is 72 dynes/cm, but some suggest 0.0728 N/m -
    % guessing this is the surface tension of seawater then. [N/m]
    %Determine the fresh/sea -water parameter c (linear between 1.2 and 1.4)
    c = 1.2+(Sal/35)*0.2;
    % Calculate the Morton number
    Morton = (g_acc * Visc ^ 4) / (density * surf_tens ^ 3);
    
    
    %Calculate for dirty bubbles
    d = 0.8
    WbuoyDirty = (((density * g_acc * rad .^ 2) / ...
        (3.68 * (Morton ^ -0.038) * Visc)) .^ -d + ...
        ((c * surf_tens) / (density * rad) + ...
        g_acc * rad) .^ (-d / 2)) .^ (-1 / d);
    
    %Calculate for clean bubbles
    d = 1.6
    WbuoyClean = (((density * g_acc * rad .^ 2) / ...
        (3.68 * (Morton ^ -0.038) * Visc)) .^ -d + ...
        ((c * surf_tens) / (density * rad) + ...
        g_acc * rad) .^ (-d / 2)) .^ (-1 / d);

end
%%%%%%%%%%%%%% Assign Wbuoy to Clean, Dirty or Vary @ 400 µm %%%%%%%%%%%%%%%%%%%
RT = 1390e-4;   % Transition radius from clean to dirty

Wbuoy = WbuoyDirty;
if  Dirty==0, Wbuoy = WbuoyClean; end
if (Dirty >0 && Dirty<1) || PlotMode==1
   i = find(rad<=80e-4); j = find(rad>80e-4);
   if length(i)>5 && length(j)>5
      Wbuoy = gausjoin(rad, WbuoyDirty, WbuoyClean, RP*1e-4, 200e-4); % Transition from Tsuchiya et al., 1997
   else
      Wbuoy = [ WbuoyDirty(i)  WbuoyClean(j) ];
   end;
   if PlotMode==1, LHndl=plot(1e4*rad,Wbuoy,'--g'); set(LHndl,'LineWidth',2); end
end

if PlotMode==1
   xlabel('Radius (µm)', 'FontSize',14,'FontWeight','Bold','FontName','Palatino');
   ylabel('Rise Velocity (cm s-1)', 'FontSize',14,'FontWeight','Bold','FontName','Palatino');
   title('Rise Velocities','FontSize',18,'FontWeight','Bold','FontName','Palatino');
   set(Hndl,'XTick',0:250:3000,'YTick',0:5:40,'XGrid','On','YGrid','On');
   hold off;
   Hndl2 = legend('Stokes - Dirty' , 'Memery & Merlivat', ...
                'Woolf 93 - Dirty', ...
                'Leifer - Clean' , 'Leifer - Dirty'  , 'Leifer - Vary');
   set(Hndl2,'Position',[0.605  0.723  0.30  0.204]);
   Wbuoy = [Hndl;Hndl2];
end
