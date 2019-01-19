function [SA_on_new_pressure,CT_on_new_pressure] = Interp_SA_CT_to_P(SA_bottle,CT_bottle,P_bottle,Pb_new,idx) % Previously called "tobenamed4_v2"
% Pre-allocation
[zl,yl]            = size(SA_bottle);
        
% This loop Interpolates SA, CT of the bottle-cast to their new 
% Pressures (P_iterative):

% Now we will provide a new string of CT and SA, in which certain
% bottles have moved in pressure space (even though they are still
% gven on the zame vertical index):
SA_on_new_pressure = SA_bottle;
CT_on_new_pressure = CT_bottle;
        
for ny = 1:yl

    % For the CAST OF THE BOTTLES, These need replacement:
    idx_ny_all_zl   = idx(:,ny); % string of length zl, with 0 and 1's. Where 1's if it needs replacement.
    if sum(idx_ny_all_zl)>=1

        % SA/CT and P of the "i" cast at which the "i" bottles lives:
     	SA_on_old_pressures = SA_bottle(:,ny); 
      	CT_on_old_pressures = CT_bottle(:,ny); 
       	old_pressures       =  P_bottle(:,ny); 
     	new_pressures       =    Pb_new(idx_ny_all_zl,ny); % New Pressures to interpolate SA and CT onto:

      	% The NOT-NaN values in the "bottle cast":
       	% This is a string of length zl, with 1's where there is NOT a NaN.
      	index     = ~isnan(squeeze(SA_on_old_pressures + CT_on_old_pressures + old_pressures));

      	% Select the NOT-NaN values for  "i bottle cast" for SA/CT and P:
     	SA_on_old_pressures = SA_on_old_pressures(index); 
       	CT_on_old_pressures = CT_on_old_pressures(index); 
      	old_pressures       =       old_pressures(index); 

        % If there are 2 or more values to interpolate between:
       	if sum(index)>1
            % Interpolating SA and CT onto a new P for the "i" bottles:
            SA_on_new_pressure(idx_ny_all_zl,ny) = interp1(old_pressures,SA_on_old_pressures,new_pressures,'linear','extrap');
            CT_on_new_pressure(idx_ny_all_zl,ny) = interp1(old_pressures,CT_on_old_pressures,new_pressures,'linear','extrap');            
        end
	clear old_pressures SA_on_old_pressures CT_on_old_pressures index       
    end
end

% figure
% imagesc(SA_on_new_pressure)
% caxis([33 35.5])
% 
% figure
% imagesc(SA_bottle-SA_on_new_pressure)
% caxis([33 35.5])
% 
% 
% tic
% DP  = P_bottle - Pb_new; % Diffrence between new pressure and bottle pressure.
% A 	= and(idx==1,DP >= 0); % If Pb_new is < P, it means that the new bottle is above the original bottle. Then DP>0.
% B   = and(idx==1,DP < 0); % If Pb_new is > P, it means that the new bottle is below the original bottle. Then DP<0.
% 
% % Below we 1) Calculate difference between new pressure, and the 
% % surrounding pressure. Then we find where is the water column, this 
% % difference is smallest.
% [~,pos]     = min(abs((ones(zl,1) * Pb_new(:)') - max(P_bottle,[],2)*ones(1,yl*zl) )); % Index (vertical) of smallest difference. 
% idx_bottles = (1:zl)'*ones(1,yl); % Index of a bottle.
% idxn        = idx_bottles - reshape(pos,[zl,yl]); % distance in index from bottle.
% bA          = idxn(A); % Select the index-differences, for the bottles above.
% bB          = idxn(B); % Select the index-differences, for the bottles below.
% bA(bA==0)   = 1; % Minimum difference must be 1 (one above).
% bB(bB==0)   = -1; % Minimum difference must be one below.
% 
% idx_A = idx_bottles(A) - bA; % Go up the water column to find bottle for interpolation.
% idx_B = idx_bottles(B) - bB; % Go down the water column to find bottle for interpolation.
% 
% idx_A(idx_A==0)=1; % Can't go into the atmosphere.
% idx_A(idx_A==zl+1)=zl; % Can't go into the bottom.
% idx_B(idx_B==0)=1;
% idx_B(idx_B==zl+1)=zl;
% 
% 
% % % idx is in fact the index of the bottle itself, so for selecting "Bottle above" we need to move these indeces 1 place up:
% % A_top      	= A(1,:); % Bottles at the surface, for which there are no "bottles above".
% % A(1,:)    	= [];
% % A(end+1,:)	= zeros(1,yl);
% % A         	= logical(A);
% % SA_B_above	= [NaN(sum(A_top),1);SA_bottle(A)];
% % CT_B_above	= [NaN(sum(A_top),1);CT_bottle(A)];
% %  P_B_above	= [NaN(sum(A_top),1); P_bottle(A)];
% % clear A A_top
% % 
% % % Bottle below
% % % idx is in fact the index of the bottle itself, so for selecting "Bottle above" we need to move these indeces 1 place down:
% % B_bottom  	= B(end,:);
% % B1         	= zeros(size(idx));
% % B1(2:end,:)	= B(1:end-1,:); 
% % B         	= logical(B1); clear idx_B_below1
% % SA_B_below	= [SA_bottle(B);NaN(sum(B_bottom),1)];
% % CT_B_below	= [CT_bottle(B);NaN(sum(B_bottom),1)];
% %  P_B_below	= [P_bottle(B);NaN(sum(B_bottom),1)];
% % clear B B_top
% 
% % Differences
% SA_on_new_pressure1 = SA_bottle;
% CT_on_new_pressure1 = CT_bottle;
% % y = (dy/dx) * (x-x0) + b
% %
% % We will use Bottle - above
% dSA     = SA_bottle(A) - SA_bottle(idx_A);
% dCT     = CT_bottle(A) - CT_bottle(idx_A);
% dP      =  P_bottle(A) - P_bottle(idx_A);
% x_x0	= Pb_new(A) - P_bottle(A);
% SA_on_new_pressure1(A)  = (dSA./dP) .* x_x0 + SA_bottle(A);
% CT_on_new_pressure1(A)  = (dCT./dP) .* x_x0 + CT_bottle(A);
% 
% clear A dSA dCT dP x_x0
% 
% % For below
% dSA     = SA_bottle(B) - SA_bottle(idx_B);
% dCT     = CT_bottle(B) - CT_bottle(idx_B);
% dP      =  P_bottle(B) - P_bottle(idx_B);
% x_x0	= Pb_new(B) - P_bottle(B);
% SA_on_new_pressure1(B)  = (dSA./dP) .* x_x0 + SA_bottle(B);
% CT_on_new_pressure1(B)  = (dCT./dP) .* x_x0 + CT_bottle(B);
% clear A dSA dCT dP x_x0
% toc
% 
% figure
% imagesc(SA_on_new_pressure1)
% caxis([33 35.5])
% 
% figure
% imagesc(SA_on_new_pressure)
% caxis([33 35.5])