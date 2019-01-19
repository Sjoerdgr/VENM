function [delta_SA_ntp, delta_CT_ntp, delta_p_ntp, p_bar] = gsw_ntp_bottle_to_cast_gradients_EDITED(SA_bottle,CT_bottle,p_bottle,SA,CT,p,drho)

% gsw_ntp_bottle_to_cast        Absolute Salinity, Conservative Temperature
%                                 and pressure on a neighbouring cast using
%                             the neutral tangent plane  (75-term equation)
%==========================================================================
%
% USAGE:
%  [SA_ntp, CT_ntp, p_ntp, p_bar] = ...
%       gsw_ntp_bottle_to_cast(SA_bottle,CT_bottle,p_bottle,SA,CT,p,{drho})
%
% DESCRIPTION:
%  This function finds the Absolute Salinity, Conservative Temperature and 
%  pressure on the casts, which have the same specific volume referred to 
%  the mid-point pressure, p_bar, as the bottle.  
%
%  This programme assumes that the neighbouring cast is stable.  It 
%  searches in one direction only, based on the initial value of the 
%  difference in the locally referenced specific volume.  It searches for
%  the corresponding locally referenced specific volume on the neighbouring 
%  cast starting from the bottle pressure.  
%  
%  Note that this function is intended to work on stabile profiles.  If 
%  there are instabillities then there are mupliple solutions it will
%  generally not find all the crossings, the results may including stable 
%  unstable zero crossings.
%
%--------------------------------------------------------------------------
%                    Still working on this section
%
%  if drho is included it does a seperate function -  see height adjustment 
%  of the surface.
%--------------------------------------------------------------------------
%
% INPUT:  
%  SA_bottle =  Absolute Salinity of the bottle                    [ g/kg ]
%  CT_bottle =  Conservative Temperature of the bottle(ITS-90)    [ deg C ]
%  p_bottle  =  sea pressure of the bottle                         [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  SA  =  Absolute Salinity of the neighbouring cast               [ g/kg ]
%  CT  =  Conservative Temperature of the neighbouring cast       [ deg C ]
%  p   =  sea pressure of the neighbouring cast                    [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  drho
%
%--------------------------------------------------------------------------
%                    Still working on this section
%
% There are two options for inputs:
%  1 s0,ct0 and p0 are scalars (one bottle). Size is (1,1).
%    s,ct and p are vectors of size (nz,1).
%  2 s0,ct0 and p0 are vectors of size (1,nx*ny).
%    s,ct and p are matrices of size (nz,nx*ny).
%---------------------------------------------------------------------------
%
%  SA_bottle, CT_bottle and p_bottle need to be single values, having 
%  dimensions 1x1.
%  SA, CT and p need to have the same dimensions Mx1 or MxN, where M is the
%  number of bottles in the cast and N is the number of profiles.  Note 
%  that they need to be arranged in a column, if they are inputed as a row
%  the programme will assume they are seperate casts. 
%
% OUTPUT:
%  SA_ntp = Absolute Salinity of the neutral tangent plane         [ g/kg ]
%  CT_ntp = Conservative Temperature of the neutral tangent plane [ deg C ]
%  p_ntp  = pressure of the neutral tangent plane                  [ dbar ]
%  p_bar  = the average pressure of the neutral tangent plane      [ dbar ]
%
% AUTHOR: David Jackett
%  Modified by Guillaume Serazin, Stefan Riha and Paul Barker
% Now Sjoerd Groeskamp
%
% VERSION NUMBER: 3.05.7 (20th December, 2016)
%
% REFERENCES: 
%  Jackett, D.R., and T.J. McDougall, 1997: A Neutral Density Variable for 
%   the World's Oceans. J. Phys. Oceanogr., 27, pp. 237–263.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%% PREAMBLE
delta = 3*1e-13; % Originally: 1e-12. The accepted tolerance when calculating the locally 
% referenced specific volume difference between the bottle and the 
% neighbouring cast.  Note that using 1e-12 is approximately equal to 0.1 dbar.

[~, Number_of_profiles] = size(SA); % Size of CAST [number of values in cast x number of casts]
size_dummy = size(SA);

if ~exist('drho','var')
    dspecvol = zeros(size_dummy); % dspecvol is a zero vector
else
    dspecvol = 1./drho;
    dspecvol(drho == 0) = 0;
end

% Pre allocate variables:
Idata         = 1:Number_of_profiles;
size_dummy(1) = 1;
delta_SA_ntp  = NaN(size_dummy);
delta_CT_ntp  = delta_SA_ntp;
delta_p_ntp   = delta_SA_ntp;
p_bar         = delta_SA_ntp;

%% Discard land (Remove NaN) Idata may contain the indices of remaining points
[SA_bottle, CT_bottle, p_bottle, SA, CT, p, dspecvol, Idata] = discard_land(SA_bottle, CT_bottle, p_bottle, SA, CT, p, dspecvol, Idata);

% NEW SIZES AS NaNs are removed:
[profile_length, Number_of_profiles] = size(SA);
Iwet                                 = 1:Number_of_profiles; % Number of bottles.
p0_stacked                           = repmat(p_bottle(:)',[profile_length 1]); % Pressure of bottles.

%% Find index for the bottle above the target bottle is found.
% Then this is turned into 2d.
I_bottle_above                = sum(p0_stacked >= p,1); % Where Pressure bottle >or= Pressure at cast. Sum this to obtain: index of bottle - 1 (shallower location).
bottles_above                 = p_bottle(:)' < p(1,:); % Where Pressure of bottle is smaller than pressure at first bottle of cast. So, if 1, bottle is shallower than the cast reaches.
I_bottle_above(bottles_above) = 1; % If bottle is shallower than cast, replace with "1", to start at first bottle of cast.                                                         % start at the uppermost cast pair
I3d                           = I_bottle_above + profile_length*(0:Number_of_profiles-1); % changing to 2-d index, based on casts supplied.
Iinc                          = NaN(1,Number_of_profiles); % Pre allocation with length of number of bottles=profiles.

% Select SA, CT and p of the bottle-in-cast, shallower than bottle itself:
SA_cast        = SA(I3d);
CT_cast        = CT(I3d);
 p_cast        =  p(I3d);
SA_cast_notNaN = ~isnan(SA_cast); % Index of cases where SA-top-bottle has values.

%% OBTAIN THE VERTICAL INDEX OF THE SHALLOWEST BOTTLE.
% Idata_nn               = ~isnan(SA); % Index of cases where SA has values.
% Idata_intergral        = cumsum(Idata_nn,1); % Number of cases where SA has values. 
% su                     = ones(size(Idata_intergral)); % Pre-allocate variable.
% su(Idata_intergral~=0) = 0; % Where SA has values, this is 0, otherwise (where SA has NaNS) this is 1.
% Ishallowest            = sum(su,1) + 1; % Vertical index of shallowest bottle.
% clear idata_nn Idata_intergral

%% OBTAIN THE VERTICAL INDEX OF THE DEEPEST BOTTLE.
% SA_dummy        = SA(end:-1:1,:,:); % Reverse the SA variable.
% Idata_nn        = ~isnan(SA_dummy); % Where does the reverse have data, and where not.
% Idata_intergral = cumsum(Idata_nn,1); % Number of cases where SA has values. 
% su              = zeros(size(Idata_intergral)); % Pre allocate variables.
% su(Idata_intergral~=0) = 0;
% Ideepest        = sum(su,1) + 1; % Vertical index of shallowest bottle.
% Ideepest        = profile_length - Ideepest + 1; % Vertical index of deepest data point

%% FIND INDEX OF FIRST AND LAST BOTTLE, WITH A VALUE.
idx                  = (1:profile_length)'*ones(1,Number_of_profiles); % Matrix of repeated vertical indeces.
NaN_to_value         = diff(isnan(SA),[],1); % This will give "-1" is there is a cave, and "1" if there is a bottom.
select_shal          = (NaN_to_value==-1); % Find index where there is a cave
Ishallowest = ones(1,Number_of_profiles); % We start of with no caves.
if sum(select_shal(:))>0 % If there are caves,
    Ishallowest(sum(select_shal,1)>0) = idx(select_shal) + 1; % Replace the place where there are caves with the index of the first bottle below the cave.
end
select_deep                                       = (NaN_to_value==1); % Find index where there is a bottom.
select_deep(profile_length,sum(select_deep,1)==0) = profile_length; % If there is "no" bottom, the deepest bottle is the last bottle.
Ideepest                                          = idx(select_deep)'; % Select indeces.

%% OBTAINING TERMS FROM DENSITY POLYNOMIAL:
go_shallow_old        = false(1,Number_of_profiles);
go_deep_old           = false(1,Number_of_profiles);
bisect                = false(1,Number_of_profiles);
done                  = false(1,Number_of_profiles);
Number_of_iterations  = 0;
[vp0,vp1,vp2,vp3,vp4] = gsw_specvol_p_parts(SA_bottle,CT_bottle);
%  vp0 = Absolute Salinity and Conservative Temperature terms for pressure to the 0
%   :
%  vp6 = Absolute Salinity and Conservative Temperature terms for pressure to the 6

%% WHILE LOOP
while any(~done)
    
    Number_of_iterations = Number_of_iterations + 1;
    p_mid                = 0.5*(p_bottle + p_cast); % Mid pressure between bottle above it on Cast, and the Bottle pressure.
    z_mid                = p_mid.*1e-4; % Depth-mid.
    
    % Now it takes the mid-depth and uses the change of density with
    % pressure, to calculate the change in pressure (for the SA and CT
    % given at the bottle) with pressure. This gives the change in density from bottle to mid-pressure.
    v_part_bottle        = vp0 + z_mid.*(vp1 + z_mid.*(vp2 + z_mid.*(vp3 + z_mid.*vp4))); % Specific volume difference

    % Now we do the same, but from the cast to the mid-pressure.
    [vp0_cast,vp1_cast,vp2_cast,vp3_cast,vp4_cast] = gsw_specvol_p_parts(SA_cast,CT_cast);
    v_part_cast                                    = vp0_cast + z_mid.*(vp1_cast + z_mid.*(vp2_cast + z_mid.*(vp3_cast + z_mid.*vp4_cast)));
    clear vp0_cast vp1_cast vp2_cast vp3_cast vp4_cast z_mid
    
    % Differences:
    v_local_diff =  v_part_bottle - v_part_cast + dspecvol; % The difference in density, at the mid-pressure.
    delta_SA_tmp = SA_bottle - SA_cast; % Differences in SA between bottle and Cast.
    delta_CT_tmp = CT_bottle - CT_cast; % Differences in CT between bottle and Cast.
    delta_p_tmp  =  p_bottle -  p_cast; % Differences in  P between bottle and Cast.
    
    %% NEXT ITERATION
    % Next iteration: Deeper or Shallower:
    % If [v_local_diff = 1/rho_bottle - 1/rho_cast > 0], the density of the cast is higher than the bottle: we need to go shallower.
    go_shallow = (v_local_diff >=  delta); % The intercept is shallower
    go_deep    = (v_local_diff <= -delta); % The intercept is deeper 
    SA_cast_notNaN(~isnan(v_local_diff)) = true; % Find where both SA_cast and v_local_diff have values (not NaN).
    
    % If v_local_diff, the locally referenced specific volume difference between
    % the bottle and the neighbouring cast, is NaN at the current pressure 
    % but it is well defined elsewhere in the watercolumn.
    
    % I_bottle_above(Iwet): For each cast, take the INDEX of the bottle-in-the-cast, shallower than the bottle.
    % Ideepest(Iwet): For each cast, take the INDEX of the deepest bottle with a value.
    % Ishallowest(Iwet): For each cast, take the INDEX of the shallowest bottle.
    nskip_up   = I_bottle_above(Iwet) - Ideepest(Iwet); % Difference in INDEX between "shallower bottle" and deepest bottle. if<0: Shallower Bottle = Deeper than data. If >0, deep bottle is deeper, 
    nskip_down = Ishallowest(Iwet) - I_bottle_above(Iwet); % Difference in INDEX between "shallower bottle" and shallowest bottle.
    
    % SELECTING TO GO UP OR DOWN:
    % isnan(v_local_diff): NaN in either bottle or cast.
    % ~SA_cast_notNaN: Where v_local_diff AND SA_cast both NO values.
    % IF nskip_up<0, than deepest bottle of cast, is deeper than the bottle.
    % IF nskip_up>0, than deepest bottle of cast, is shallower than bottle. -NOT GOOD
    % IF nskip_down<0, than shallowest bottle of cast, is shallower than the bottle.
    % IF nskip_down>0, than shallowest bottle of cast, is deeper than bottle. -NOT GOOD
    go_shallow_nan = isnan(v_local_diff) & ~SA_cast_notNaN & (nskip_up>0); % NaN if all options lead to NaN
    go_deep_nan    = isnan(v_local_diff) & ~SA_cast_notNaN & (nskip_down>0); % NaN if all options lead to NaN
    
    % FINAL DECISION TO GO UP OR DOWN THE WATER COLUMN:
    go_shallow     = go_shallow | go_shallow_nan; % If either of these is 1, we go shallower 
    go_deep        = go_deep | go_deep_nan; % If either of these is 1, we go deeper.
    if sum((go_shallow(:) + go_deep(:))==2)>0; disp('Error: We can go both shallower and deeper'); end
    
    % BOTTLES THAT ARE WITHIN LIMITS:
    final                            = abs(v_local_diff) < delta;  % INDEX of bottles that are within tolerence level.
    % Idata = inex of bottles. Iwet = Index of bottles that still need to be done. final = index of bottles that are done.
    % Idata(Iwet(final)): Select bottles that are done, of the bottles that needed to be done, and replaces these for the total bottles.
    delta_SA_ntp(Idata(Iwet(final))) = delta_SA_tmp(final); % Replace SA
    delta_CT_ntp(Idata(Iwet(final))) = delta_CT_tmp(final); % Replace CR
    delta_p_ntp(Idata(Iwet(final)))  = delta_p_tmp(final); % Replace P
    p_bar(Idata(Iwet(final)))        = p_mid(final); % Replace Pmid.
    
    % See if the mid-point pressure is crossed from above or below,
    % starting at the bottle (I think). So is the slope positive or
    % negative?
    cfb            = (go_shallow_old & go_deep & SA_cast_notNaN); % crossed from below
    cfa            = (go_deep_old & go_shallow & SA_cast_notNaN); % crossed from above
    crossed        = cfb | cfa; % If the mid pressure is crossed at all.
    start_bis      = (crossed & ~bisect);   % If crossed and bisect, start bisection here
    bisect         = (bisect | start_bis);  % If bisect, or bisect is started, than bisect here
    search_initial = (go_deep | go_shallow) & ~bisect & ~final; % If need to go deep or shallw, and bisect and not part of final, continue serch.
    
    % CREATE NEW VECTOR OF INDEXCES FOR THE ONES THAT STILL NEED TO BE DONE:
    Iinc(Iwet(go_deep & search_initial))    = 1;  % Select index of bottles (Iinc) and set to 1, for the bottles that need to be done (Iwet), only if they go deep and need to be done.
    Iinc(Iwet(go_shallow & search_initial)) = -1; % Select index of bottles (Iinc) and set to -1, for the bottles that need to be done (Iwet), only if they go deep and need to be done.
    
    % For the ones that need further work (sear_initial) and need to go
    % deeper or shallower, select their indeces and replace them with the
    % difference between bottle and shallowest or deepest in cast, index.
    Iinc(Iwet(go_deep_nan & search_initial))    = nskip_down(Iwet(go_deep_nan & search_initial));
    Iinc(Iwet(go_shallow_nan & search_initial)) = -nskip_up(Iwet(go_shallow_nan & search_initial));
    
    % The new INC has now been adapted so it moves up and down the water column.
    I_bottle_above = I_bottle_above + Iinc; % Get bottle above.
    I3d            = I3d + Iinc; % Shallowest bottle Plus change
    Iwet_tmp       = I_bottle_above(Iwet(search_initial));  % select ones that need to be done.
    I3d_wet        = I3d(Iwet(search_initial)); % here too
    
    % Remove some that are ow beyond the domain given by the cast:
    out                  = (Iwet_tmp < 1) | (Iwet_tmp > profile_length);   % above or below domain 
    out2                 = false(1,length(final)); 
    out2(search_initial) = out; 
    done                 = (isnan(v_local_diff) & SA_cast_notNaN) | final | out2; 
    % Done gives the ones for which differences are NaN, and SA has values,
    % or final or out2? 
    
    I3d_wet         = I3d_wet(~out);
    search_initial  = search_initial(~done);
    bisect          = bisect(~done);
    crossed         = crossed(~done);
    SA_cast_notNaN  = SA_cast_notNaN(~done);
    
    if ~all((search_initial & ~bisect) | (~search_initial & bisect)) % either bisect or keep searching
        error('gsw_ntp_bottle_to_cast: Something is wrong, there may be a NaN in the input data)')
    end
    
    if Number_of_iterations > 1
        SA1 = SA1(~done);
        CT1 = CT1(~done);
         p1 =  p1(~done);
        
        SA2 = SA2(~done);
        CT2 = CT2(~done);
         p2 =  p2(~done);
        
        SA2(crossed) = SA1(crossed);
        CT2(crossed) = CT1(crossed);
         p2(crossed) =  p1(crossed);
    else
        SA2 = SA_cast(~done);
        CT2 = CT_cast(~done);
         p2 =  p_cast(~done);
    end
    
    SA1 = SA_cast(~done);
    CT1 = CT_cast(~done);
     p1 =  p_cast(~done);
    
    % data for next evaluation of the difference in the locally referenced 
    % specific volume between the bottle and the cast
   SA_bottle = SA_bottle(~done);
   CT_bottle = CT_bottle(~done);
    p_bottle = p_bottle(~done);
    dspecvol = dspecvol(~done);
    vp0 = vp0(~done);
    vp1 = vp1(~done);
    vp2 = vp2(~done);
    vp3 = vp3(~done);
    vp4 = vp4(~done);
    
    SA_cast = NaN(1,sum(~done));
    CT_cast = NaN(1,sum(~done));
     p_cast = NaN(1,sum(~done));
    
    SA_cast(search_initial) = SA(I3d_wet);
    CT_cast(search_initial) = CT(I3d_wet);
     p_cast(search_initial) =  p(I3d_wet);
    
    if any(bisect)
        SA_cast(bisect) = 0.5*(SA1(bisect) + SA2(bisect));
        CT_cast(bisect) = 0.5*(CT1(bisect) + CT2(bisect));
         p_cast(bisect) = 0.5*(p1(bisect) + p2(bisect));
    end
    
    go_shallow_old = go_shallow(~done);
    go_deep_old    = go_deep(~done);
    Iwet           = Iwet(~done);
end

end


function [SA_bottle,CT_bottle,p_bottle,SA,CT,p,dspecvol,Idata] = discard_land(SA_bottle,CT_bottle,p_bottle,SA,CT,p,dspecvol,Idata)

% discard_land                             eliminate land from the data set
%==========================================================================
%
% USAGE:
%  [SA_ntp,CT_ntp,p_ntp] = ...
%                 discard_land(SA_bottle,CT_bottle,p_bottle,SA,CT,p,{drho})
%
% DESCRIPTION:
%
%==========================================================================

Iwet = ~isnan(SA_bottle);
Idry = all(isnan(SA));
Iwet = Iwet & ~Idry;

SA_bottle = SA_bottle(Iwet);
CT_bottle = CT_bottle(Iwet);
p_bottle = p_bottle(Iwet);

SA = SA(:,Iwet);
CT = CT(:,Iwet);
p = p(:,Iwet);

Idata = Idata(Iwet);

dspecvol = dspecvol(Iwet);

end
