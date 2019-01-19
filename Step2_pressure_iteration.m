function [SA_ntp_mid,CT_ntp_mid,P_ntp_mid,dzds,dSAds,dCTds,dPds] = Step2_pressure_iteration(z, lat_i, lat_j, lon_i, lon_j, SA_i, SA_j, CT_i, CT_j, P_i, P_j)

%% SETTINGS PRIOR TO ITERATION LOOP
% Starting Mid-Pressure on the Interface between bottle and cast:
[zl,yl]     = size(SA_i); % Dimensions Bottle.
[zlt,ylt]   = size(SA_j); % Dimensions Cast.
if and(zlt==zl,ylt==yl)==0; disp('Warning: dimensions much agree.'); stop; end; clear zlt ylt

% Calculate Target Pressure.
% (Mid Pressure between two cast, at same depth of the T-grid points).
lat_m       = 0.5 * (lat_i + lat_j);        % Latitude of Target mid-point.
lat_m_2D    = ones(zl,1)*lat_m;             % 2-Dimensional Latitude of Target mid-point [Repeat zl times x latitude] 
z_2D        = z'*ones(1,yl);                % Depth in two dimensions [depth x Repeat yl times].
P_target    = gsw_p_from_z(z_2D,lat_m_2D);  % Target Pressure.
ny_2D       = ones(zl,1)*(1:yl);            % Indeces for bottle and cast.
clear z_2D lat_m_2D

% Pre-Settings for the iteration-loop.
n_max       = 12;                   % Maximal Number of iterations
n           = 1;                    % Start for counting.
DP_crit     = 0.5;                  % The amount of dbars (~meters) close we want to get to Pmid. Note, max = min(dz)/2 ~1.25 in ML, and 25 at depth.
% When decreasing DP too much, the result becomes less neutral. 
% Still need to explain why.

% Create Exponentional Decaying Relaxation Parameter.
a       = log(1e-3);            % End-point of decay.
n0      = 2;                    % Number of "ones" before decaying.
ns      = n0:n_max;             % Number of factors to determine.
dr      = 150;                  % Decay Rate.
Frelax  = [ones(1,n0),exp((ns-n0).*a./dr)]; % Relaxation Number.

% The Bottle-Pressures, which will change iteratively during the iteration-loop:
Pb_new_i   = P_i; % Bottle Pressures on "i" cast, as given by input.
Pb_new_j   = P_j; % Bottle Pressures on "i" cast, as given by input.

%% Bottles and Casts
% Pre-allocation of final variables 
idx           = ones(zl,yl);

% For i:
dSA_ntp_i     = NaN(zl,yl);
dCT_ntp_i     = NaN(zl,yl);
 dP_ntp_i     = NaN(zl,yl);
SA_ntp_mid_i  = NaN(zl,yl);
CT_ntp_mid_i  = NaN(zl,yl);
 P_ntp_mid_i  = NaN(zl,yl);
dz_ntp_i      = NaN(zl,yl);

% For j:
dSA_ntp_j     = NaN(zl,yl);
dCT_ntp_j     = NaN(zl,yl);
 dP_ntp_j     = NaN(zl,yl);
SA_ntp_mid_j  = NaN(zl,yl);
CT_ntp_mid_j  = NaN(zl,yl);
 P_ntp_mid_j  = NaN(zl,yl);
dz_ntp_j      = NaN(zl,yl);

% Initial calculation of NTP variables for i: 
[dSA_ntp_i,dCT_ntp_i,dP_ntp_i,SA_ntp_mid_i,CT_ntp_mid_i,P_ntp_mid_i,dz_ntp_i] = ...
 Step1_ntp_search_support(SA_i,CT_i,P_i,lat_i,SA_j,CT_j,P_j,lat_j,ny_2D,idx,...
 dSA_ntp_i,dCT_ntp_i,dP_ntp_i,SA_ntp_mid_i,CT_ntp_mid_i,P_ntp_mid_i,dz_ntp_i); % Before this was called "tobenamed3"

% Initial calculation of NTP variables for j: 
[dSA_ntp_j,dCT_ntp_j,dP_ntp_j,SA_ntp_mid_j,CT_ntp_mid_j,P_ntp_mid_j,dz_ntp_j] = ...
 Step1_ntp_search_support(SA_j,CT_j,P_j,lat_j,SA_i,CT_i,P_i,lat_i,...
 ny_2D,idx,dSA_ntp_j,dCT_ntp_j,dP_ntp_j,SA_ntp_mid_j,CT_ntp_mid_j,P_ntp_mid_j,dz_ntp_j);

%% Calculate new variables to initate convergence of ntp-variables:
% Difference between Target Pressure and mid-pressure of the NTP:
DP_i           = P_ntp_mid_i - P_target; % Mid-pressure - Target Pressure.
DP_j           = P_ntp_mid_j - P_target; % Mid-pressure - Target Pressure.

% We now select either i2j or j2i, of DP, depending on which one is 
% closest to the Target Pressure. Here DP will be the minimum of both
% versions and "idx_min" will keep track which of the versions (i2j or
% j2i) was used:
DPx             = NaN([size(DP_i),2]);
DPx(:,:,1)      = abs(DP_i);
DPx(:,:,2)      = abs(DP_j);  
[DP,~]          = min(DPx,[],3); clear DPx 
idx             = abs(DP)>DP_crit;    % Index of location with DP>criteria [zl x yl].
redo            = nansum(idx(:));     % Number of profiles that require re-iteration.

% To obtain an better estimate, subtract the difference and
% re-calculate the NTP (making use of the relaxation factor).
% If DP>0, it means that P_avg is deeper, and the new starting 
% pressure should be shallower (lower pressure). The new pressure 
% should be P_new = Pressure_bottle - DP.
Pb_new_i       = Pb_new_i - Frelax(n) .* DP_i; % New bottle Pressure to start ntp calculations from (see below).
Pb_new_j       = Pb_new_j - Frelax(n) .* DP_j; % New bottle Pressure to start ntp calculations from (see below).
    
%% The iterative while loop:    
while and(n<n_max+1,redo>0) % As long as n is smaller than n_max AND there is more than 1 ntp to be redone.
    % The function below uses that idx is a zl x yl string.
    % All output should also be zl x yl
    [SAb_new_i,CTb_new_i] = Interp_SA_CT_to_P(SA_i,CT_i,P_i,Pb_new_i,idx); % Linear interpolation and extrapolation
    [SAb_new_j,CTb_new_j] = Interp_SA_CT_to_P(SA_j,CT_j,P_j,Pb_new_j,idx); % Linear interpolation and extrapolation
    
    % Use new values to re-calculate NTP variables for i:
    [dSA_ntp_i,dCT_ntp_i,dP_ntp_i,SA_ntp_mid_i,CT_ntp_mid_i,P_ntp_mid_i,dz_ntp_i] = ... 
     Step1_ntp_search_support(SAb_new_i,CTb_new_i,Pb_new_i,lat_i,SA_j,CT_j,P_j,lat_j,ny_2D,idx,...
     dSA_ntp_i,dCT_ntp_i,dP_ntp_i,SA_ntp_mid_i,CT_ntp_mid_i,P_ntp_mid_i,dz_ntp_i);

     % Use new values to re-calculate NTP variables for J:
    [dSA_ntp_j,dCT_ntp_j,dP_ntp_j,SA_ntp_mid_j,CT_ntp_mid_j,P_ntp_mid_j,dz_ntp_j] = ... 
     Step1_ntp_search_support(SAb_new_j,CTb_new_j,Pb_new_j,lat_j,SA_i,CT_i,P_i,lat_i,ny_2D,idx,...
     dSA_ntp_j,dCT_ntp_j,dP_ntp_j,SA_ntp_mid_j,CT_ntp_mid_j,P_ntp_mid_j,dz_ntp_j);
    clear SAb_new_i CTb_new_i SAb_new_j CTb_new_j idx DP DP_i DP_j 

    % As above:
    DP_i            = P_ntp_mid_i - P_target; % Mid-pressure - Target Pressure.
    DP_j            = P_ntp_mid_j - P_target; % Mid-pressure - Target Pressure.
    DPx             = NaN([size(DP_i),2]);
    DPx(:,:,1)      = abs(DP_i);
    DPx(:,:,2)      = abs(DP_j);  
    [DP,~]          = min(DPx,[],3); clear DPx  
    idx             = abs(DP)>DP_crit;    % Index of location with DP>criteria [zl x yl].
    redo            = nansum(idx(:));     % Number of profiles that require re-iteration.
    Pb_new_i        = Pb_new_i - Frelax(n) .* DP_i; % New bottle Pressure to start ntp calculations from (see below).
    Pb_new_j        = Pb_new_j - Frelax(n) .* DP_j; % New bottle Pressure to start ntp calculations from (see below).
    n               = n + 1; % Set stage for next iteration.
end % End the while-loop

%% ALLOCATE BEST VALUES
% We must now select the bottles and cast we use, for our 
% estimates for the tracer gradients and slopes, based on which of 
% the two versions (i2j or j2i) provided the most accurate result.
% if exist('min_from_j','var')==0
%    idx_min=0;
% end
% min_from_j    = idx_min==2; % Find indeces for which j2i was most accurate.

% Select
DPx        	= NaN([size(DP_i),2]);         % Allocate Matrix.
DPx(:,:,1)	= abs(P_ntp_mid_i - P_target); % Difference between Target and Mid-pressure.
DPx(:,:,2)  = abs(P_ntp_mid_j - P_target); % Difference between Target and Mid-pressure.
[~,idx_min]	= min(DPx,[],3); clear DPx % Which of two approaches has minimum.
min_from_j  = idx_min==2; % where the "J to i" leads to best values.

% First Allocate the values to be that from "i to j":
dSA_ntp     = dSA_ntp_i;
dCT_ntp     = dCT_ntp_i;
 dP_ntp     =  dP_ntp_i;
SA_ntp_mid  = SA_ntp_mid_i;
CT_ntp_mid  = CT_ntp_mid_i;
 P_ntp_mid  = P_ntp_mid_i;
dz_ntp      = dz_ntp_i;

% Replace where "j to i" was more accurate.
% Note that we define the direction to be diff(x) = xj - xi, so the 
% j2i, which is xi - xj, gets a minus sign. dSA/dCT and dP for j2i:
dSA_ntp(min_from_j)     =-dSA_ntp_j(min_from_j);
dCT_ntp(min_from_j)     =-dCT_ntp_j(min_from_j);
 dP_ntp(min_from_j)     = -dP_ntp_j(min_from_j);
 dz_ntp(min_from_j)     = -dz_ntp_j(min_from_j);
  P_ntp_mid(min_from_j) =   P_ntp_mid_j(min_from_j);
 SA_ntp_mid(min_from_j) =  SA_ntp_mid_j(min_from_j);
 CT_ntp_mid(min_from_j) =  CT_ntp_mid_j(min_from_j);


%% Remove incropping and outcropping
P_too_Small = P_ntp_mid<0; % Values that Outcrop beyond the surface.
P_too_Large = P_ntp_mid > (repmat(max([P_i;P_j]),[zl,1]) + 25); % Values that Incrop into the bottom.
index       = or(P_too_Small,P_too_Large);

% Replacement:
dSA_ntp(index)          = NaN;
dCT_ntp(index)          = NaN;
 dP_ntp(index)          = NaN;
  P_ntp_mid(index)      = NaN;
 SA_ntp_mid(index)      = NaN;
 CT_ntp_mid(index)      = NaN;
 dz_ntp(index)          = NaN;

%% INTERPOLATE NAN AND VALUES THAT DID NOT ITERATE ENOUGH
% Select Profiles that need interpolation.
yni_dp   = sum(idx,1)>0; % Index of cast that require further work, DUE TO ITERATION.
ydx      = 1:yl; % Index of each cast.
cdx      = ydx(yni_dp); % Indeces of cast that have re-iteration requirment..
if isempty(cdx)==0
    for n = 1:length(cdx)
        cdxn    = cdx(n); % Cast to interp.
        zdx     = idx(:,cdxn); % Where this is "1", interpolation is required.
        sel     = ~isnan(P_ntp_mid(:,cdxn)); % Data points that are not NaN.
        
        if sum(sel)>1
            Pq      = P_target(zdx,cdxn); % Get the right cast, get the pressure of ther value.
            Pi      = P_ntp_mid(sel,cdxn);       % Pressures

            % Interp Gradients:
            dSA_ntp(zdx,cdxn)  = interp1(Pi,dSA_ntp(sel,cdxn),Pq,'linear');
            dCT_ntp(zdx,cdxn)  = interp1(Pi,dCT_ntp(sel,cdxn),Pq,'linear');
             dP_ntp(zdx,cdxn)  = interp1(Pi, dP_ntp(sel,cdxn),Pq,'linear');
             dz_ntp(zdx,cdxn)  = interp1(Pi, dz_ntp(sel,cdxn),Pq,'linear');

            % Interp Values:
            SA_ntp_mid(zdx,cdxn)  = interp1(Pi,SA_ntp_mid(sel,cdxn),Pq,'linear');
            CT_ntp_mid(zdx,cdxn)  = interp1(Pi,CT_ntp_mid(sel,cdxn),Pq,'linear');
             P_ntp_mid(zdx,cdxn)  = Pq;
         
        elseif sum(sel)<=1
            % Interp Gradients:
            dSA_ntp(zdx,cdxn)  = NaN;
            dCT_ntp(zdx,cdxn)  = NaN;
             dP_ntp(zdx,cdxn)  = NaN;
             dz_ntp(zdx,cdxn)  = NaN;

            % Interp Values:
            SA_ntp_mid(zdx,cdxn)  = NaN;
            CT_ntp_mid(zdx,cdxn)  = NaN;
             P_ntp_mid(zdx,cdxn)  = NaN;
         
        end
        clear Pi Pq cdxn zdx sel
    end
end
clear cdx

%% TWO-DIMENSIONAL NEAREST NEIGHBOUR INTERPOLATION IN VERTICAL
org_nan     = isnan(SA_i); % Where the original data has NaN
new_nan     = isnan(P_ntp_mid + SA_ntp_mid + CT_ntp_mid); % Where the iterated data has NaN.
dif_nan     = (new_nan - org_nan)>0;  % Where new data has nan, but old data does not.
[z2d,h2d]   = ndgrid(1:zl,1:yl);
values      = new_nan==0; % Where there are values.
zval        = z2d(values);
hval        = h2d(values);

% INTERPOLATE:
% For SA_ntp_mid:
F                   = scatteredInterpolant(zval,hval,SA_ntp_mid(values),'natural','nearest');
SA_ntp_mid(dif_nan) = F(z2d(dif_nan),h2d(dif_nan)); clear F

% For CT_ntp_mid:
F                   = scatteredInterpolant(zval,hval,CT_ntp_mid(values),'natural','nearest');
CT_ntp_mid(dif_nan) = F(z2d(dif_nan),h2d(dif_nan)); clear F

% For P_ntp_mid:
F                   = scatteredInterpolant(zval,hval,P_ntp_mid(values),'natural','nearest');
P_ntp_mid(dif_nan)  = F(z2d(dif_nan),h2d(dif_nan)); clear F

% For dSA_ntp:
F                   = scatteredInterpolant(zval,hval,dSA_ntp(values),'natural','nearest');
dSA_ntp(dif_nan)    = F(z2d(dif_nan),h2d(dif_nan)); clear F

% For dCT_ntp:
F                   = scatteredInterpolant(zval,hval,dCT_ntp(values),'natural','nearest');
dCT_ntp(dif_nan)    = F(z2d(dif_nan),h2d(dif_nan)); clear F

% For dP_ntp:
F                   = scatteredInterpolant(zval,hval,dP_ntp(values),'natural','nearest');
dP_ntp(dif_nan)     = F(z2d(dif_nan),h2d(dif_nan)); clear F

% For dz_ntp:
F                   = scatteredInterpolant(zval,hval,dz_ntp(values),'natural','nearest');
dz_ntp(dif_nan)     = F(z2d(dif_nan),h2d(dif_nan)); clear F

%% Calculate Isopycnal Tracer Gradients (C/m).
% Calculate distance
deg2m        = 6378000*2*pi/360; % meters/degree longitude @ equator.
dlon         = ones(zl,1)*(lon_j - lon_i); % Difference in longitude (degrees).
dlon_per_lat = ones(zl,1)*cos(pi*lat_m/180); % meters/degree as a function of latitude
dx           = dlon .* dlon_per_lat .* deg2m; % Longitudinal distance (m). 
dy           = ones(zl,1)*(lat_j - lat_i)*deg2m; % Difference in latitude (m).
ds           = sqrt(dx.^2 + dy.^2); % Total distance between bottle and cast (m).

% For the structure of ds, I'm asuming that in general either dx or dy is
% zero, depending on which gradient is calculated. To have a generelized
% function, I therefore did it this way, leaving the gradients:
dSAds = dSA_ntp./ds;
dCTds = dCT_ntp./ds;
 dPds =  dP_ntp./ds;
 dzds =  dz_ntp./ds;