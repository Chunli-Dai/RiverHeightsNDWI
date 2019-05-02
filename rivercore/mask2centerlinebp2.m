function [c]=mask2centerline(data,flagwidth)

%% 
 %Cons: 1\sometimes it is not working. imregionalmin deleted too many points. 2\ need es parameter.
Wn = 10; % nominal width
es = 'SW'; % exit sides
es = 'WE'; % START stream at W, End stream at E.
es ='EW';
plotornot=1;
% 3a. Compute the centerline
[cl, Icl] = centerline_from_mask(data.z,es,Wn,plotornot);
x=data.x(cl(:,1));y=data.y(cl(:,2));
[lat,lon]=polarstereo_inv(x,y,[],[],70,-45);
c.X=lon;c.Y=lat; %.X (longitude), .Y (latitude).

% hold on;plot(lon,lat,'k-')

if nargin>1 && flagwidth==1

%% Get river width; Copied from RivMAP/DEMO.m
Ist=data.z;
resx=mean(data.x(2:end)-data.x(1:end-1));resy=mean(data.y(2:end)-data.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);

if 0 % Index exceeds array bounds.; Error in banklines_from_mask (line 230)
%% 3b. Compute the banklines
close all
banks = banklines_from_mask(Ist, es, plotornot);

%% 3c. Along-channel widths

% Channel width at each pixel from banklines
Wbl = width_from_banklines(cl, banks{1}, banks{2}, Wn);

Wavg = mean(Wbl(isnan(Wbl)==0)); % average of pixel-wise widths - note that it doesn't quite agree with our nominal channel width (Wn), but this is fine since Wn is just used to parameterize buffer boxes
end %if 0 ; width from bankline

% Channel width from mask at specified spacing
% spacing = Wn/2; %spacing of centerline, in pixels (nodes of centerlines).
rescl=resr;
spacing = max([round(100/rescl),1]); %100 m; spacing of centerline, in pixels (nodes of centerlines).
[Wm, SWm] = width_from_mask(Ist, cl, spacing);

% Let's plot the two width methods to see how they compare
% First, we need to compute the streamwise distance along the centerline
S = [0; cumsum(sqrt(diff(cl(:,1)).^2+diff(cl(:,2)).^2))];

% Now plot both widths
figure
% close all
% plot(S,Wbl); hold on
plot(SWm*resr,Wm*resr,'r');
xlabel('streamwise distance (m)'); ylabel('width (m)')
legend('W_m_a_s_k')

%% 3e. Reach-averaged widths

% Finally, we can compute the average width for the entire reach
Wra = sum(sum(Ist))/S(end); %Ist not limited by buffer zone; overestimated.

% This value agrees fairly well with the average of widths computed on a
% pixel-wise basis:
Wra
% Wavg
end %if flagwidth==1

close all
return
end

