% This program obtains motion discrimination thresholds for
% contrast modulated random noise using interleaved QUEST staircases.
% Davis Glasser
% Last Edited: 03/07/2011

%-----Subject Settings------------------
% (Change these for every subject)
stimulus_radius                 = 47;                                  % pixels
mv_length                       = 17;
px_size                         = 2;                                    % arcmin/pix, plz use even multiple of scale factor
dynamic                         = 1;                                    % 1=new noise each frame.
SF                              = 1/25;                                    % Cycles/deg
TFstep                          = pi/4;                                    % Hz
angle                           = 0;                                    % Degress, 0 = horizontal, 90 = vertical
modulation_depth                = 100;
contrast                        = 99;                                  % Percent
background                      = 128;                                  % Grayscale Units

%-----Leave everything below here alone-
%-----Rig Settings----------------------
scale_factor                    = 1;                                    % Arcmin/pixel
frame_rate                      = 360;                                  % Screen frame rate (hz)
f                               = (SF)*2*pi;
angle                           = angle*pi/180;
a                               = cos(angle)*f;
b                               = sin(angle)*f;
amplitude                       = background*contrast/100;
scale                           = round(px_size/scale_factor);

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
pedestal = imresize(round(rand(bps,bps))*2-1,scale,'nearest');
pedestal = pedestal(1:bps,1:bps);
movie = zeros(1,mv_length);
for i = 1:mv_length;
    moving_grating = sin(a*x+b*y+motion_step(i)).*circle*amplitude+background;
    imwrite(uint8(moving_grating),strcat('firstord_',num2str(i),'.png'),'png');
end
