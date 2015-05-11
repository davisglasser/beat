% This program obtains motion discrimination thresholds for
% contrast modulated random noise using interleaved QUEST staircases.
% Davis Glasser
% Last Edited: 03/07/2011

%-----Subject Settings------------------
% (Change these for every subject)
stimulus_radius                 = 239;                                  % pixels
mv_length                       = 1;
SF                              = 1/50;                                    % Cycles/deg
TFstep                          = pi/4;                                    % Hz
angle                           = 0;                                    % Degress, 0 = horizontal, 90 = vertical
contrast                        = 20;                                  % Percent
background                      = 128;                                  % Grayscale Units

%-----Leave everything below here alone-
%-----Rig Settings----------------------
scale_factor                    = 1;                                    % Arcmin/pixel
f                               = (SF)*2*pi;
angle                           = angle*pi/180;
a                               = cos(angle)*f;
b                               = sin(angle)*f;
amplitude                       = background*contrast/100;

%-----Spatial Envelope------------------
[x,y]=meshgrid(-stimulus_radius:stimulus_radius,-stimulus_radius:stimulus_radius);
bps = round((stimulus_radius)*2+1);
circle=((stimulus_radius)^2-(x.^2+y.^2));
for i=1:bps; for j =1:bps; if circle(i,j) < 0; circle(i,j) = 0; else circle(i,j) = 1; end; end;
end;
R = (sqrt(x.^2 + y.^2) + eps).*circle;
R = R/max(max(R));
cos2D = (cos(R*pi)+1)/2;
circle = (cos2D.*circle);

% Calculate the grating motion
direction = 2;

motion_step = zeros(1,mv_length);
motion_step(1) = rand*2*pi;
for i=2:mv_length
    motion_step(i) = motion_step(i-1)+TFstep*((direction-1)*(-2)+1);
end

% Make the movie
movie = zeros(bps,bps,1,mv_length);
for i = 1:mv_length;
    moving_grating = sin(a*x+b*y+motion_step(i)).*circle*amplitude+background;
    movie(:,:,1,i) = moving_grating;
end
imwrite(uint8(movie),'LCbig_still.gif','gif','DelayTime',.0312,'Loopcount',inf);