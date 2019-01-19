close all
clear variables
load('/Users/Sjoerd/Desktop/Work/Data/WOA/WOA_gsw_N2LL_plus.mat')

%% CLEAR VARIABLES TO CLEAR SPACE
clear gamma_original ShallowMask Rh0 Rh0_surf
clear b_surf beta_surf alpha_surf Cb_surf Tb_surf
clear gamma_surf Rh0 Chla N2_surf
clear LH_OA SH_OA LW_OA SW_OA E_OA P_OA R_OA LH_Core SH_Core LW_Core SW_Core E_Core P_Core R_Core
clear albedo_diffusive albedo_direct albedo_seferian albedo_seferian_noCHLA albedo_whc 
clear cloud_cover Chla_m
clear rho rho_surf alpha beta Cb Tb N2 b
clear Axy Axz Ayz dx dy dz mlp 
clear cp0 CTmax CTmin deg2m dir dR DR fig_bgc fig_vis gamma lat_A lon_A ml_pres_max ml_pres_min Mycmap
clear N2min Redge rho0 rho_fresh Rl Rmax Rmin Rrange SAmax SAmin savelocation_figures savingtime Slope_max
clear time_stamp vol xl yl zl tl xt_edges yt_edges zt_edges

%% A FUNCTION TO CALCULATE NEUTRAL GRADIENTS AND SLOPES.
% 
% Here we input a 4 dimensional matrix of SA, CT and P and use it to
% calculate isopycnal slopes, at the interfaces between two adjacent casts
% in either the x or y direction. In the vertical we will calculate these
% variables on both the tracer grid-depths and at the interfaces between
% vertical grid points.
%
% INPUT: 
% - SA, CT and P, all are [lon x lat x depth x time] at tracer-grid 
%   locations.
% - xt, yt, and zt are 1D strings of their coordinates.
% - SA_surf, CT_surf and P_surf, all are [lon x lat x time] at zt=0.
%
% OUTPUT
% - Sy_Axz     , Sx_Ayz     , Sx_Axy    , Sy_Axy:
%   These are the Slope of the isopycnal, in either the y or the x 
%   direction, at the interfaces.
%
% - dSAdy_N_Axz, dCTdy_N_Axz, dPdy_N_Axz, DP_Axz
%   The SA, CT and P gradients at the Axz interface (y-direction)
%
% - dSAdx_N_Ayz, dCTdx_N_Ayz, dPdx_N_Axz, DP_Ayz
%   The SA, CT and P gradients at the Ayz interface (x-direction)
%
% - dSAdx_N_Axy, dCTdx_N_Axy, dPdx_N_Axy
%   The SA, CT and P gradients at the Axy interface (z-direction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PREAMBLE
[xl,yl,zl,tl] = size(SA);

%% ADD SURFACE VALUES TO THE SA, CT AND P PROFILES THAT WE USE.
SA1     = NaN(xl,yl,zl+1,tl); SA1(:,:,2:end,:)  =   SA;   SA1(:,:,1,:) = SA_surf     ; SA   =   SA1; clear SA1 SA_surf
CT1     = NaN(xl,yl,zl+1,tl); CT1(:,:,2:end,:)  =   CT;   CT1(:,:,1,:) = CT_surf     ; CT   =   CT1; clear CT1 CT_surf
P1      = NaN(xl,yl,zl+1)   ;  P1(:,:,2:end)    =    p;     P1(:,:,1 ) = zeros(xl,yl); p    =    P1; clear P1  p_surf
%topo1   = NaN(xl,yl,zl+1,tl); topo1(:,:,2:end,:)= topo; topo1(:,:,1,:) = topo_surf   ; topo = topo1; clear topo1 topo_surf
zt      = [0,zt];

%% PRE-ALLOCATION FOR OUTPUT VARIABLES
% Pre-allocate the Sx and Sy values on Tgrid:
Sx     = NaN(xl,yl,zl,tl);
Sy     = NaN(xl,yl,zl,tl);

% Pre-allocate the Axz (in y-direction) interfaces:
     Sy_Axz     = NaN(xl,yl,zl,tl);
dSAdy_N_Axz     = NaN(xl,yl,zl,tl);
dCTdy_N_Axz     = NaN(xl,yl,zl,tl);
 dPdy_N_Axz     = NaN(xl,yl,zl,tl);

% Pre-allocate the Ayz (in x-direction) interfaces:
     Sx_Ayz     = NaN(xl,yl,zl,tl);
dSAdx_N_Ayz     = NaN(xl,yl,zl,tl);
dCTdx_N_Ayz     = NaN(xl,yl,zl,tl);
 dPdx_N_Ayz     = NaN(xl,yl,zl,tl);

% Pre-allocate the Axy (in z-direction) interfaces, on MID-depth:
Sx_Axy          = NaN(xl,yl,zl,tl); % Pre-allocate the Ayz (in x-direction) interfaces, on MID-depth:
Sy_Axy          = NaN(xl,yl,zl,tl); % Pre-allocate the Axz (in y-direction) interfaces, on MID-depth:
dSAdz_N_Axy     = NaN(xl,yl,zl,tl);
dCTdz_N_Axy     = NaN(xl,yl,zl,tl);
 dPdz_N_Axy     = NaN(xl,yl,zl,tl);

% Pre-allocate the (x,y,z) directions, for SA, CT and P, AT the T-grid
dSAdx_N_tg     = NaN(xl,yl,zl,tl);
dSAdy_N_tg     = NaN(xl,yl,zl,tl);
dSAdz_N_tg     = NaN(xl,yl,zl,tl);
dCTdx_N_tg     = NaN(xl,yl,zl,tl);
dCTdy_N_tg     = NaN(xl,yl,zl,tl);
dCTdz_N_tg     = NaN(xl,yl,zl,tl);
dPdx_N_tg      = NaN(xl,yl,zl,tl);
dPdy_N_tg      = NaN(xl,yl,zl,tl);
dPdz_N_tg      = NaN(xl,yl,zl,tl);
     
%% LOOP TO CALCULATE SLOPES AND GRADIENTS:

for nt = 1:tl
    tic
    disp(['Calculating slopes for the month = ', num2str(nt)])
    
    %% ADD MID-DEPTH VALUES BY LINEAR INTERPOLATION OF SA, CT AND P
    SA_nt   = squeeze(SA(:,:,:,nt));
    CT_nt   = squeeze(CT(:,:,:,nt));
    topo_nt = squeeze(topo(:,:,:,nt));
    % Calculate values at the Axy interface.
    SA_mid = 0.5*(SA_nt(:,:,1:end-1) + SA_nt(:,:,2:end));
    CT_mid = 0.5*(CT_nt(:,:,1:end-1) + CT_nt(:,:,2:end));
     p_mid = 0.5*( p(:,:,1:end-1)    +  p(:,:,2:end));
    zt_mid = 0.5*(zt(1:end-1)        + zt(2:end))    ;

    % Preallocation of a new variable to include valies + mid-values.
    zlpm      	= 2*(zl+1)-1;      % Vertical dimension plus mid values.
    SApm      	= NaN(xl,yl,zlpm); % SA plus mid (pm) values.
    CTpm      	= NaN(xl,yl,zlpm); % CT plus mid (pm) values.
     Ppm      	= NaN(xl,yl,zlpm); %  P plus mid (pm) values.
    ztpm        = NaN(1,zlpm);

    % Combine values and mid-values.
    SApm(:,:,1:2:zlpm) = SA_nt;  % Enter original values at odd locations.
    SApm(:,:,2:2:zlpm) = SA_mid; % Enter mid values at even locations.
    CTpm(:,:,1:2:zlpm) = CT_nt;  % Enter original values at odd locations.
    CTpm(:,:,2:2:zlpm) = CT_mid; % Enter mid values at even locations.
     Ppm(:,:,1:2:zlpm) = p;      % Enter original values at odd locations.
     Ppm(:,:,2:2:zlpm) = p_mid;  % Enter mid values at even locations.
    ztpm(1:2:zlpm)     = zt;     % Enter mid values at even locations.
    ztpm(2:2:zlpm)     = zt_mid; % Enter mid values at even locations.
    z                  = -ztpm; % Height on the vertical interfaces.
    clear SA_mid CT_mid p_mid zt_mid CT_nt

    %% THE X-DIRECTION: ON A AYZ SURFACE, FOR DIFFERENT NY:
    % Pre-allocate the Ayz (in x-direction) interfaces, on MID-depth:
         Sx_Ayz_mid = NaN(xl,yl,zl+1);
    dSAdx_N_Ayz_mid = NaN(xl,yl,zl+1);
    dCTdx_N_Ayz_mid = NaN(xl,yl,zl+1);
     dPdx_N_Ayz_mid = NaN(xl,yl,zl+1);
 
    for ny = 1:yl
        disp(['xdir, ny = ',num2str(ny)])

        % Select x-z section, for the Bottles, and make it z by x:
        SA_i    = squeeze(SApm(:,ny,:))'; % was 1:xl-1
        CT_i    = squeeze(CTpm(:,ny,:))';
         P_i    = squeeze( Ppm(:,ny,:))';

        % Select x-z section, for the Cast [EAST of (next to) the bottle], and make it z by x:
        SA_j	= squeeze(SApm([2:xl,1],ny,:))'; % was 2:xl
        CT_j	= squeeze(CTpm([2:xl,1],ny,:))';
         P_j	= squeeze( Ppm([2:xl,1],ny,:))';

        % If all values are NaN, skip this position.
        if sum(isnan(SA_i + CT_i + P_i + SA_j + CT_j + P_j),1)==zlpm; continue; end

        % Latitude:
        lat_i = yt(ny)*ones(1,xl);
        lat_j = lat_i;
        lon_i = xt;
        lon_j = [xt(2:xl), 360 + xt(1)]; % NOTE THIS WORKS ONLY FOR GLOBAL OCEAN.

        % WE NEED TO ADD THE SLOPE CALCULATION
        [~,~,~,Sx_2d,dSAdx_N,dCTdx_N,dPdx_N] = Step2_pressure_iteration(z, lat_i, lat_j, lon_i, lon_j, SA_i, SA_j, CT_i, CT_j, P_i, P_j); % Before this was called "tobenamed_v2"
        clear SA_i CT_i P_i SA_j CT_j P_j lat_i lat_j lon_i lon_j

        % Now allocate the Tracer-Grid points to the right variable:
        % ny,nt    = Current latitude and time.
        % 3:2:zlpm = Selects Values at original Tracer grid location.
        %            However, the first one is AT the surface, which is not
        %            where the first grid point is. So skip the first.
             Sx_Ayz(:,ny,:,nt) = Sx_2d(3:2:zlpm,:)';
        dSAdx_N_Ayz(:,ny,:,nt) = dSAdx_N(3:2:zlpm,:)';
        dCTdx_N_Ayz(:,ny,:,nt) = dCTdx_N(3:2:zlpm,:)';
         dPdx_N_Ayz(:,ny,:,nt) = dPdx_N(3:2:zlpm,:)';

        % Now allocate the Tracer-Grid MID-points to the right variable:
        % 4:2:zlpm = Selects Values at Interfaces. However, the first one 
        %            is supposed to be below the first grid point. While 
        %            here the first one is given between the surface and 
        %            the first grid point. Hence we start at "4".
        
        
        % This is what I did befofe I used multiple months nt = 1:>1
        % Sx_Ayz_mid(:,ny,1:zl-1,nt)      = Sx_2d(4:2:zlpm,:)';
        % dSAdx_N_Ayz_mid(:,ny,1:zl-1,nt) = dSAdx_N(4:2:zlpm,:)';
        % dCTdx_N_Ayz_mid(:,ny,1:zl-1,nt) = dCTdx_N(4:2:zlpm,:)';
        %  dPdx_N_Ayz_mid(:,ny,1:zl-1,nt) = dPdx_N(4:2:zlpm,:)';
         
        Sx_Ayz_mid(:,ny,1:zl-1)      = Sx_2d(4:2:zlpm,:)';
        dSAdx_N_Ayz_mid(:,ny,1:zl-1) = dSAdx_N(4:2:zlpm,:)';
        dCTdx_N_Ayz_mid(:,ny,1:zl-1) = dCTdx_N(4:2:zlpm,:)';
         dPdx_N_Ayz_mid(:,ny,1:zl-1) = dPdx_N(4:2:zlpm,:)';
        clear SAm_x CTm_x Pm_x Sx_2d dSAdx_N dCTdx_N dPdx_N
        
        % Now calculate Sx on Tgrid location, using Sx_Ayz:
        Sx_Tgrid        = squeeze(0.5 .* (Sx_Ayz(:,ny,:,nt) + Sx_Ayz([end,1:end-1],ny,:,nt)));
        index           = (isnan(squeeze(SA_nt(:,ny,2:end))) - isnan(Sx_Tgrid))<0;
        [x2d,z2d]       = ndgrid(1:xl,1:zl);
        values          = index==0; % Where there are values.
        xval            = x2d(values);
        zval            = z2d(values);
        F               = scatteredInterpolant(xval,zval,Sx_Tgrid(values),'nearest','nearest');
        Sx_Tgrid(index) = F(x2d(index),z2d(index)); clear F
        Sx(:,ny,:,nt)   = Sx_Tgrid;            
        clear index x2d z2d values xval zval Sx_Tgrid
    end

    %% THE Y-DIRECTION: ON A AXZ SURFACE, FOR DIFFERENT NX:
    % Pre-allocate the Axz (in y-direction) interfaces, on MID-depth:
         Sy_Axz_mid = NaN(xl,yl,zl+1);
    dSAdy_N_Axz_mid = NaN(xl,yl,zl+1);
    dCTdy_N_Axz_mid = NaN(xl,yl,zl+1);
     dPdy_N_Axz_mid = NaN(xl,yl,zl+1);
     
    for nx = 1:xl
        disp(['ydir, nx = ',num2str(nx)])

        % Select y-z section, for the Bottles, and make it z by Y:
        SA_i    = squeeze(SApm(nx,1:yl-1,:))';
        CT_i    = squeeze(CTpm(nx,1:yl-1,:))';
         P_i    = squeeze( Ppm(nx,1:yl-1,:))';

        % Select y-z section, for the Cast [NORTH of (next to) the bottle], and make it z by y:
        SA_j	= squeeze(SApm(nx,2:yl,:))'; 
        CT_j	= squeeze(CTpm(nx,2:yl,:))';
         P_j	= squeeze( Ppm(nx,2:yl,:))';

        % If all values are NaN, skip this position.
        if sum(isnan(SA_i + CT_i + P_i + SA_j + CT_j + P_j),1)==zlpm; continue; end 

        % Make Lat for the y-direction:
        lat_i   = yt(1:yl-1);
        lat_j   = yt(2:yl);
        lon_i   = xt(nx)*ones(1,yl-1);
        lon_j   = lon_i;

        % WE NEED TO ADD THE SLOPE CALCULATION
        [~,~,~,Sy_2d,dSAdy_N,dCTdy_N,dPdy_N] = Step2_pressure_iteration(z, lat_i, lat_j, lon_i, lon_j, SA_i, SA_j, CT_i, CT_j, P_i, P_j);
        clear SA_i CT_i P_i SA_j CT_j P_j lat_i lat_j lon_i lon_j

        % Now allocate the Tracer-Grid points to the right variable:
        % nx,nt    = Current longitude and time.
        % 1:yl-1   = First service, is the one to the north of the first
        %            T-grid point. Meaning, the last does not have a value.
        % For the rest, see "X direction description".
             Sy_Axz(nx,1:yl-1,:,nt)  = Sy_2d(3:2:zlpm,:)';
        dSAdy_N_Axz(nx,1:yl-1,:,nt)  = dSAdy_N(3:2:zlpm,:)';
        dCTdy_N_Axz(nx,1:yl-1,:,nt)  = dCTdy_N(3:2:zlpm,:)';
         dPdy_N_Axz(nx,1:yl-1,:,nt)  = dPdy_N(3:2:zlpm,:)';

        % This is what I did befofe I used multiple months nt = 1:>1
        % Now allocate the Tracer-Grid points to the right variable:
        %     Sy_Axz_mid(nx,1:yl-1,1:zl-1,:,nt)  = Sy_2d(4:2:zlpm,:)';
        % dSAdy_N_Axz_mid(nx,1:yl-1,1:zl-1,:,nt)  = dSAdy_N(4:2:zlpm,:)';
        % dCTdy_N_Axz_mid(nx,1:yl-1,1:zl-1,:,nt)  = dCTdy_N(4:2:zlpm,:)';
        % dPdy_N_Axz_mid(nx,1:yl-1,1:zl-1,:,nt)  = dPdy_N(4:2:zlpm,:)'; 
         
             Sy_Axz_mid(nx,1:yl-1,1:zl-1,:)  = Sy_2d(4:2:zlpm,:)';
        dSAdy_N_Axz_mid(nx,1:yl-1,1:zl-1,:)  = dSAdy_N(4:2:zlpm,:)';
        dCTdy_N_Axz_mid(nx,1:yl-1,1:zl-1,:)  = dCTdy_N(4:2:zlpm,:)';
         dPdy_N_Axz_mid(nx,1:yl-1,1:zl-1,:)  = dPdy_N(4:2:zlpm,:)'; 
        clear SAm_y CTm_y Pm_y Sy_2d dSAdy_N dCTdy_N dPdy_N
        
        % Now calculate Sx on Tgrid location, using Sx_Ayz:
        Sy_Tgrid            = NaN(yl,zl);
        Sy_Tgrid(2:end,:)   = squeeze(0.5 .* (Sy_Axz(nx,1:end-1,:,nt) + Sy_Axz(nx,2:end,:,nt)));
        index               = (isnan(squeeze(SA_nt(nx,:,2:end))) - isnan(Sy_Tgrid))<0;
        [y2d,z2d]           = ndgrid(1:yl,1:zl);
        values              = index==0; % Where there are values.
        yval                = y2d(values);
        zval                = z2d(values);
        F                   = scatteredInterpolant(yval,zval,Sy_Tgrid(values),'nearest','nearest');
        Sy_Tgrid(index)     = F(y2d(index),z2d(index)); clear F
        Sy(nx,:,:,nt)       = Sy_Tgrid;            
        clear index y2d z2d values yval zval Sy_Tgrid
    end
    
    %% THE Z-DIRECTION: ON A AXY SURFACE, USING THE RESULTS FROM THE X AND Y DIRECTION:
    
    % 1) - Calculate vertical component, at Ayz (x-direction), for Axy interface.
    dSAdz_N_Ayz  = Sx_Ayz_mid .* dSAdx_N_Ayz_mid;
    dCTdz_N_Ayz  = Sx_Ayz_mid .* dCTdx_N_Ayz_mid;
     dPdz_N_Ayz  = Sx_Ayz_mid .*  dPdx_N_Ayz_mid;
    clear dSAdx_N_Ayz_mid dCTdx_N_Ayz_mid dPdx_N_Ayz_mid
    
    % 2) - Calculate slopes/gradients at Axy interface of T-grid, from the x-direction:
            Sx_Axy_3d = 0.5 .* ( Sx_Ayz_mid +  Sx_Ayz_mid([end,1:end-1],:,:));
    dSAdz_N_Ayz_Tgrid = 0.5 .* (dSAdz_N_Ayz + dSAdz_N_Ayz([end,1:end-1],:,:));
    dCTdz_N_Ayz_Tgrid = 0.5 .* (dCTdz_N_Ayz + dCTdz_N_Ayz([end,1:end-1],:,:));
     dPdz_N_Ayz_Tgrid = 0.5 .* ( dPdz_N_Ayz +  dPdz_N_Ayz([end,1:end-1],:,:));
    clear dSAdz_N_Ayz dCTdz_N_Ayz dPdz_N_Ayz Sx_Ayz_mid
    
    % 3) - Calculate vertical component, at Axz (y-direction), for Axy interface.
    dSAdz_N_Axz  = Sy_Axz_mid .* dSAdy_N_Axz_mid;
    dCTdz_N_Axz  = Sy_Axz_mid .* dCTdy_N_Axz_mid;
     dPdz_N_Axz  = Sy_Axz_mid .*  dPdy_N_Axz_mid;
     clear dSAdy_N_Axz_mid dCTdy_N_Axz_mid dPdy_N_Axz_mid
     
    % 4) - Calculate slopes/gradients at Axy interface of T-grid, from the y-direction:
    % Pre-allocate to the right size:
            Sy_Axy_3d = NaN(xl,yl,zl+1);
    dSAdz_N_Axz_Tgrid = NaN(xl,yl,zl+1);
    dCTdz_N_Axz_Tgrid = NaN(xl,yl,zl+1);
     dPdz_N_Axz_Tgrid = NaN(xl,yl,zl+1);

    % Calculate average:
            Sy_Axy_3d(:,2:end,:)  = 0.5 .* ( Sy_Axz_mid(:,1:end-1,:) +  Sy_Axz_mid(:,2:end,:));
    dSAdz_N_Axz_Tgrid(:,2:end,:)  = 0.5 .* (dSAdz_N_Axz(:,1:end-1,:) + dSAdz_N_Axz(:,2:end,:));
    dCTdz_N_Axz_Tgrid(:,2:end,:)  = 0.5 .* (dCTdz_N_Axz(:,1:end-1,:) + dCTdz_N_Axz(:,2:end,:));
     dPdz_N_Axz_Tgrid(:,2:end,:)  = 0.5 .* ( dPdz_N_Axz(:,1:end-1,:) +  dPdz_N_Axz(:,2:end,:));

    % APPROXIMATE first value:
            Sy_Axy_3d(:,1,:)  =  Sy_Axz_mid(:,1,:);
    dSAdz_N_Axz_Tgrid(:,1,:)  = dSAdz_N_Axz(:,1,:);
    dCTdz_N_Axz_Tgrid(:,1,:)  = dCTdz_N_Axz(:,1,:);
     dPdz_N_Axz_Tgrid(:,1,:)  =  dPdz_N_Axz(:,1,:);

    % APPROXIMATE last value:
            Sy_Axy_3d(:,end,:)  =  Sy_Axz_mid(:,end,:);
    dSAdz_N_Axz_Tgrid(:,end,:)  = dSAdz_N_Axz(:,end,:);
    dCTdz_N_Axz_Tgrid(:,end,:)  = dCTdz_N_Axz(:,end,:);
     dPdz_N_Axz_Tgrid(:,end,:)  =  dPdz_N_Axz(:,end,:);
    clear Sy_Axz_mid dSAdz_N_Axz dCTdz_N_Axz dPdz_N_Axz

    % 4) - Calculate vertical component:
    Sx_Axy(:,:,:,nt)       = Sx_Axy_3d(:,:,2:end);
    Sy_Axy(:,:,:,nt)       = Sy_Axy_3d(:,:,2:end);
    dSAdz_N_Axy(:,:,:,nt)  = dSAdz_N_Ayz_Tgrid(:,:,2:end) + dSAdz_N_Axz_Tgrid(:,:,2:end);
    dCTdz_N_Axy(:,:,:,nt)  = dCTdz_N_Ayz_Tgrid(:,:,2:end) + dCTdz_N_Axz_Tgrid(:,:,2:end);
     dPdz_N_Axy(:,:,:,nt)  =  dPdz_N_Ayz_Tgrid(:,:,2:end) +  dPdz_N_Axz_Tgrid(:,:,2:end);
    clear dSAdz_N_Ayz_Tgrid dCTdz_N_Ayz_Tgrid dPdz_N_Ayz_Tgrid
    clear dSAdz_N_Axz_Tgrid dCTdz_N_Axz_Tgrid dPdz_N_Axz_Tgrid 
    clear SApm Ppm CTpm nx ny SA_nt Sx_Axy_3d Sy_Axy_3d CT_nt

    %% NTP TRACER-GRADIENTS, AT T-GRID, FROM GRADIENTS AT INTERFACES
    % FOR SA
    dSAdy_N_3d                           = NaN(xl,yl,zl);
    dSAdz_N_3d                           = NaN(xl,yl,zl);
    dSAdx_N_3d                           = 0.5 .* squeeze(dSAdx_N_Ayz([end,1:end-1],      :,:,nt) + dSAdx_N_Ayz(:,    :,:,nt)); % At the Interface, Deeper/Below of Tgrid.
    dSAdy_N_3d(:,2:end,:)                = 0.5 .* squeeze(dSAdy_N_Axz(            :,1:end-1,:,nt) + dSAdy_N_Axz(:,2:end,:,nt)); % At the Interface, Deeper/Below of Tgrid.
    dSAdz_N_3d(:,2:end,:)                = 0.5 .* squeeze(dSAdz_N_Axy(            :,1:end-1,:,nt) + dSAdz_N_Axy(:,2:end,:,nt)); % At the Interface, Deeper/Below of Tgrid.
    [dSAdx_N_3d, dSAdy_N_3d, dSAdz_N_3d] = interpolate_values_interface_to_tgrid(dSAdx_N_3d, dSAdy_N_3d, dSAdz_N_3d, topo_nt);

    % FOR CT
    dCTdy_N_3d                           = NaN(xl,yl,zl);
    dCTdz_N_3d                           = NaN(xl,yl,zl);
    dCTdx_N_3d                           = 0.5 .* squeeze(dCTdx_N_Ayz([end,1:end-1],      :,:,nt) + dCTdx_N_Ayz(:,    :,:,nt)); % At the Interface, Deeper/Below of Tgrid.
    dCTdy_N_3d(:,2:end,:)                = 0.5 .* squeeze(dCTdy_N_Axz(            :,1:end-1,:,nt) + dCTdy_N_Axz(:,2:end,:,nt)); % At the Interface, Deeper/Below of Tgrid.
    dCTdz_N_3d(:,2:end,:)                = 0.5 .* squeeze(dCTdz_N_Axy(            :,1:end-1,:,nt) + dCTdz_N_Axy(:,2:end,:,nt)); % At the Interface, Deeper/Below of Tgrid.
    [dCTdx_N_3d, dCTdy_N_3d, dCTdz_N_3d] = interpolate_values_interface_to_tgrid(dCTdx_N_3d, dCTdy_N_3d, dCTdz_N_3d, topo_nt);

    % FOR P
    dPdy_N_3d                            = NaN(xl,yl,zl);
    dPdz_N_3d                            = NaN(xl,yl,zl);
    dPdx_N_3d                            = 0.5 .* squeeze(dPdx_N_Ayz([end,1:end-1],      :,:,nt) + dPdx_N_Ayz(:,    :,:,nt)); % At the Interface, Deeper/Below of Tgrid.
    dPdy_N_3d(:,2:end,:)                 = 0.5 .* squeeze(dPdy_N_Axz(            :,1:end-1,:,nt) + dPdy_N_Axz(:,2:end,:,nt)); % At the Interface, Deeper/Below of Tgrid.
    dPdz_N_3d(:,2:end,:)                 = 0.5 .* squeeze(dPdz_N_Axy(            :,1:end-1,:,nt) + dPdz_N_Axy(:,2:end,:,nt)); % At the Interface, Deeper/Below of Tgrid.
    [dPdx_N_3d, dPdy_N_3d, dPdz_N_3d]    = interpolate_values_interface_to_tgrid(dPdx_N_3d, dPdy_N_3d, dPdz_N_3d, topo_nt); clear topo_nt
    
    % ALLOCATE TO TOTAL
	dSAdx_N_tg(:,:,:,nt) = dSAdx_N_3d; clear dSAdx_N_3d
    dSAdy_N_tg(:,:,:,nt) = dSAdy_N_3d; clear dSAdy_N_3d
    dSAdz_N_tg(:,:,:,nt) = dSAdz_N_3d; clear dSAdz_N_3d
    dCTdx_N_tg(:,:,:,nt) = dCTdx_N_3d; clear dCTdx_N_3d
    dCTdy_N_tg(:,:,:,nt) = dCTdy_N_3d; clear dCTdy_N_3d
    dCTdz_N_tg(:,:,:,nt) = dCTdz_N_3d; clear dCTdz_N_3d
     dPdx_N_tg(:,:,:,nt) = dPdx_N_3d; clear dPdx_N_3d
     dPdy_N_tg(:,:,:,nt) = dPdy_N_3d; clear dPdy_N_3d
     dPdz_N_tg(:,:,:,nt) = dPdz_N_3d; clear dPdz_N_3d
    
    
    %% Saving.
    disp(['Saving for month = ', num2str(nt)])
    disp('Saving Slopes on Interface')
    
    Slopes_interface.Sy_Axz           = Sy_Axz;
    Slopes_interface.Sx_Ayz           = Sx_Ayz;
    Slopes_interface.Sx_Axy           = Sx_Axy;
    Slopes_interface.Sy_Axy           = Sy_Axy;
    
    Slopes_interface.dSAdx_N_Ayz      = dSAdx_N_Ayz;
    Slopes_interface.dCTdx_N_Ayz      = dCTdx_N_Ayz;
    Slopes_interface.dPdx_N_Ayz       =  dPdx_N_Ayz;
    
    Slopes_interface.dSAdy_N_Axz      = dSAdy_N_Axz;
    Slopes_interface.dCTdy_N_Axz      = dCTdy_N_Axz;
    Slopes_interface.dPdy_N_Axz       =  dPdy_N_Axz;
    
    Slopes_interface.dSAdz_N_Axy      = dSAdz_N_Axy;
    Slopes_interface.dCTdz_N_Axy      = dCTdz_N_Axy;
    Slopes_interface.dPdz_N_Axy       =  dPdz_N_Axy;
   
    save('/Users/Sjoerd/Desktop/Work/Data/WOA/WOA_Slopes_interface.mat','Slopes_interface','-v7.3')
    
    disp('Saving Slopes on T-grid')
        
    Slopes_Tgrid.Sx       	= Sx;
    Slopes_Tgrid.Sy       	= Sy;

	Slopes_Tgrid.dSAdx_N	= dSAdx_N_tg;
    Slopes_Tgrid.dCTdx_N   	= dCTdx_N_tg;
    Slopes_Tgrid.dPdx_N    	=  dPdx_N_tg;
    
    Slopes_Tgrid.dSAdy_N  	= dSAdy_N_tg;
    Slopes_Tgrid.dCTdy_N  	= dCTdy_N_tg;
    Slopes_Tgrid.dPdy_N   	=  dPdy_N_tg;
    
    Slopes_Tgrid.dSAdz_N  	= dSAdz_N_tg;
    Slopes_Tgrid.dCTdz_N   	= dCTdz_N_tg;
    Slopes_Tgrid.dPdz_N   	=  dPdz_N_tg;
    
    save('/Users/Sjoerd/Desktop/Work/Data/WOA/WOA_Slopes_Tgrid.mat','Slopes_Tgrid','-v7.3')
    
    disp(['Done Saving for month = ', num2str(nt)])
toc
end