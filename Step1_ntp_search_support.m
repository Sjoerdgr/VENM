function [dSA_ntp,dCT_ntp,dP_ntp,SA_ntp_mid,CT_ntp_mid,P_ntp_mid,dz_ntp] = Step1_ntp_search_support(...
          SA_bottle,CT_bottle,P_bottle,lat_bottle,...
          SA_cast,  CT_cast,  P_cast,  lat_cast,...
          ny_2D,idx,dSA_ntp,dCT_ntp,dP_ntp,SA_ntp_mid,CT_ntp_mid,P_ntp_mid,dz_ntp) % Previously called "tobenamed3"

% The Bottles, Casts and "idx" should have input zl x yl.
[zl,~] = size(SA_bottle);

%% Unify NaNs:
% Local "topography", make sure they all have NaNs at the same places.
% For CAST:
cast_sum    = SA_cast + CT_cast + P_cast;
topo_cast   = cast_sum - cast_sum;
SA_cast     = SA_cast + topo_cast;
CT_cast     = CT_cast + topo_cast;
 P_cast     =  P_cast + topo_cast;

 % For BOTTLE:
bottle_sum  = SA_bottle + CT_bottle + P_bottle;
topo_bottle = bottle_sum - bottle_sum;
SA_bottle   = SA_bottle + topo_bottle;
CT_bottle   = CT_bottle + topo_bottle;
 P_bottle   =  P_bottle + topo_bottle;

%% INDEXING
% Take the original cast. If the cast has enough data, it gives a "1".
% For each ny, repeat this with depth, leaving casts of either o's or 1's.
% It should be [zl,yl]:
Enough_Data_in_cast  = repmat(sum((~isnan(cast_sum)),1)>=2,[zl,1]); %  Index where there is more than 2 values in the verticle.
% Enough_Data_in_cast = 1, if there is enough data, and the cast can be used.
% Enough_Data_in_cast = 0, otherwise.

Bottle_is_not_NaN  = (~isnan(SA_bottle) + ~isnan(CT_bottle) + ~isnan(P_bottle))==3;
% Bottle_is_not_NaN = 1, if the bottle is NOT a NaN.
% Bottle_is_not_NaN = 0, otherwise.

index = (idx + Enough_Data_in_cast + Bottle_is_not_NaN)==3; 
% If total == 0, All values are NaN, and the bottle does not reuire replacement.
% If total == 1, Either: There is not enough data in cast, the bottle has already been calculated, or the bottle is a NaN.
% If total == 2, One of the 3, is valid, but the other will nevertheless prohibit it from being correctly calculated.
% if total == 3, We will calculate new values.

%% SELECT BOTTLES AND CORRESPONDING CASTS THAT REQUIRE NTP ESTIMATES.
% Select index:
ny_idx_new_values = ny_2D(index);

% SELECT BOTTLES AND CASTS
% The bottles for SA/CT/P on the "i" cast: 
SAb    = SA_bottle(index)';
CTb    = CT_bottle(index)';
 Pb    =  P_bottle(index)';
 
% The "j Cast" corresponding to the bottles on the "i" cast, for SA/CT/P:
SAc    = SA_cast(:,ny_idx_new_values);
CTc    = CT_cast(:,ny_idx_new_values);
 Pc    =  P_cast(:,ny_idx_new_values);

%% CALCULATING VALUES ON NTP
% Find dSA, dCT, dP and Pmid for "i-bottles" to "j" cast, 
[dSA, dCT, dP, Pmid] = gsw_ntp_bottle_to_cast_gradients_EDITED(SAb,CTb,Pb,SAc,CTc,Pc); % i2j

% Note: "gsw_ntp_bottle_to_cast_gradients" does : dC = Bottle - Cast = I - J = South - North
% We use it the other way around, so we need a MINUS SIGN!!!
% also:
% Cast      = Bottle - dC
% C_ntp_mid = 0.5 x (Bottle + Cast)
%           = 0.5 x (Bottle + Bottle - dC)
%           = Bottle - 0.5 x dC
% This leaves:
dSA_ntp(index)    = -dSA; % ADDED MINUS SIGN, 17:30, 23-10-2017.
dCT_ntp(index)    = -dCT; % ADDED MINUS SIGN, 17:30, 23-10-2017.
 dP_ntp(index)    =  -dP; % ADDED MINUS SIGN, 17:30, 23-10-2017.
SA_ntp_mid(index) = SAb - 0.5 * dSA; % Use  C_ntp_mid = Bottle - 0.5 x dC
 P_ntp_mid(index) = Pmid;
 
% We now obtain the CT at the mid-pressure, not by averaging along the NTP,
% but by using a function that takes into account the nonlineairites of the
% EOS. We do this only for CT, as its nonlineairites are about a factor 10
% stronger than nonlineairites in SA:
LRPD              = gsw_rho(SAb,CTb,Pmid); % The LRPD of the NTP:             
CT_ntp_mid(index) = gsw_CT_from_rho(LRPD,SA_ntp_mid(index)',Pmid); % CT, obtained taking into account the NEOS.

%% Calculate the depth difference between bottle and cast for the NTP:
P_ntp           = Pb - dP;  % Use Cast = Bottle - dC.
latb            = ones(zl,1)*lat_bottle; % Latitude of "i" endpoint.
latc            = ones(zl,1)*lat_cast;   % Latitude of "j" endpoint.
zb              = gsw_z_from_p(Pb,latb(index));    % Depth of Bottle on "i" endpoint, for i2j.
zc              = gsw_z_from_p(P_ntp,latc(index)); % Depth on cast on "j" endpoint, for i2j.
dz_ntp(index)   = zc - zb; % Cast minus Bottle

% DEBUGGING FOR PAUL.
% Send on 11-10-2017, because the ntp code fails on line 
% where we have SA_bottle - SA_tmp, because SA_tmp is resized but SA_bottle
% is not.

% Paul_WOA.SA_bottle      = SAb;
% Paul_WOA.CT_bottle      = CTb;
% Paul_WOA.p_bottle       = Pb;
% Paul_WOA.SA_cast      = SAc;
% Paul_WOA.CT_cast      = CTc;
% Paul_WOA.p_cast       = Pc;
% save('/Users/Sjoerd/Desktop/Paul_WOA.mat','Paul_WOA','-v7.3')