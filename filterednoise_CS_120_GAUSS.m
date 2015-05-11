% This program tests motion discrimination of luminance-defined
% filtered noise at fixed durations.
% GAUSSIAN TEMPORAL ENVELOPES
% Davis Glasser
% Last Edited: 5/08/2013

clear all;close all;clc;HideCursor;ListenChar(2);                           % Start with a blank slate

try
    %-----Subject Settings------------------
    % (Change these for every subject)
    initials                        = 'test';
    time_sigma                      = [.75 1 1.5 2 2.5 3 3.5 4];
    stimulus_radius                 = 480;                                  % Arcmin
    contrast_rms                    = .15;                                  % RMS right now
    contrast                        = 100;
    
    orientation_mean                = 0;
    orientation_width               = 1;                                   %mean +/- width for sharp, sd for gauss
    orientation_profile             = 1;                                    % 0 = sharp, 1 = gaussian
    
    SF_mean                         = 1;                                    % 0 = broadband;
    SF_width                        = .125;                                    % octaves, moot if broadband;
    SF_profile                      = 0;                                    % 0 = sharp, 2 = gaussian
    
    speed                           = 4;                                    % deg/s
    n_trials                        = 10;                                   % Number of trials per duration, 40 is good
    
    fast                            = 1;                                    % Automatically trigger trials
    feedback                        = 0;
    %-----Experiment Settings---------------
    % (Keep this consistent across subjects)
    experiment_id                   = 'noise_gauss';                        % Used in group filename
    ITI                             = 2;                                    % Intertrial Inverval, Seconds
    background                      = 126;                                  % Grayscale Units
    fixate                          = 0;                                    % Present fixation spot during motion
    % (Bad for fovea, good for periphery)
    H_ecc_stim                      = 0;                                    % Horizontal stimulus ecc (degs, neg is left)
    H_ecc_fix                       = 0;                                    % Horizontal fixation ecc (degs, neg is left)
    V_ecc_stim                      = 0;                                    % Vertical stimulus ecc (degs, neg is up)
    V_ecc_fix                       = 0;                                    % Vertical fixation ecc (degs, neg is up)
    
    linearize                       = 1;                                    % Use calibrated LUT (do this when available)
    spatial_envelope                = 2;                                    % 0 = disk, 1 = Gabor, 2 = raised cosine
    which_envelope                  = 0;
    
    data_path                       = '/Users/tadinlab/Desktop/DGProjector/Data/';                 % Folder for saving data files
    
    %-----Leave everything below here alone-
    %-----Rig Settings----------------------
    scale_factor                    = 2;                                    % Arcmin/pixel
    frame_rate                      = 120;                                  % Screen frame rate (hz)
    
    %-----Housekeeping----------------------
    % Scale things based on viewing distance, and convert other stuff to
    % the units PsychToolbox wants...
    step                            = speed*(60/scale_factor)/frame_rate;
    tme                             = clock;
    n_staircases                    = length(time_sigma);
    results                         = zeros(n_staircases,n_trials);
    stim_dir                        = results;
    count                           = ones(1,n_staircases);
    stimulus_radius                 = stimulus_radius /scale_factor;
    Gaussian_stdev                  = round(stimulus_radius/1.5);
    H_ecc_stim                      = H_ecc_stim*60/scale_factor;
    H_ecc_fix                       = H_ecc_fix*60/scale_factor;
    V_ecc_stim                      = V_ecc_stim*60/scale_factor;
    V_ecc_fix                       = V_ecc_fix*60/scale_factor;
    amplitude                       = background*contrast/100;    
    %-----Randomize Trial Order-------------
    total_trials                    = n_trials*n_staircases;
    perm                            = randperm(total_trials);
    perm                            = mod(perm,n_staircases)+1;
    
    %-----Spatial Envelope------------------
    [x,y]=meshgrid(-stimulus_radius:stimulus_radius,-stimulus_radius:stimulus_radius);
    bps = (stimulus_radius)*2+1;
    circle=((stimulus_radius)^2-(x.^2+y.^2));
    for i=1:bps; for j =1:bps; if circle(i,j) < 0; circle(i,j) = 0; else circle(i,j) = 1; end; end;
    end;
    if spatial_envelope == 1
        circle = (exp(-(((x)/(sqrt(2)*Gaussian_stdev/6)).^2)-((y/(sqrt(2)*Gaussian_stdev/2)).^2)).*circle);
    elseif spatial_envelope == 2
        R = (sqrt(x.^2 + y.^2) + eps).*circle;
        R = R/max(max(R));
        cos2D = (cos(R*pi)+1)/2;
        circle = (cos2D.*circle);
    end
    
    %-----Open Screens----------------------
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    w=Screen('OpenWindow',screenNumber,0,[],8,2);
    screen_rect = Screen('Rect',w);
    if size(screens,2) == 2;
        w2=Screen('OpenWindow',0,0,[],[],2);
        Screen('FillRect',w2, 0); Screen('Flip', w2);
    end
    if linearize
        screen_clut = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
        screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    end
    Screen('FillRect',w, background);
    Screen('Flip', w);
    Screen('FillRect',w, background);
    Screen('TextSize',w,20);Screen('TextFont',w,'Charcoal');
    
    %-----Stimulus Rectangles---------------
    movie_rect= [0,0,bps,bps];
    scr_left_middle = fix(screen_rect(3)/2)-round(bps/2);
    scr_top = fix(screen_rect(4)/2)-round(bps/2);
    screen_rect_middle = movie_rect + [scr_left_middle, scr_top, scr_left_middle, scr_top];
    screen_patch = screen_rect_middle+[H_ecc_stim,V_ecc_stim,H_ecc_stim,V_ecc_stim];
    sr_hor = round(screen_rect(3)/2);
    sr_ver = round(screen_rect(4)/2);
    fix_hor = sr_hor+H_ecc_fix;
    fix_ver = sr_ver+V_ecc_fix;
    
    Screen('DrawText',w,'Motion Discrimination: Filtered Noise (GAUSS)',100,100,0);
    Screen('DrawText',w,'Use the DOWN arrow to start a trial, and LEFT/RIGHT arrows to respond',100,130,0);
    Screen('DrawText',w,[int2str(total_trials),'  trials - press SPACE BAR to start'],100,160,0);
    Screen('Flip',w);
    
    FlushEvents('keyDown');
    validKey = 0;
    while ~validKey
        [secs, keyCode, deltaSecs] = KbWait;
        if keyCode(KbName('space'))
            validKey = 1;
        end
    end
    
    Screen('FillRect', w, background);
    Screen('Flip', w);
    tic;
    %-----Main experimental loop-----------------
    for trial=1:total_trials
        aa = GetSecs;
        % Draw the white fixation cross
        Screen('FillRect',w, background);
        Screen('DrawLine',w,255,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
        Screen('DrawLine',w,255,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
        Screen('Flip',w);
        
        % Calculate the temporal envelope
        [time_gauss,mv_length] = envelope((1000*(time_sigma(perm(trial)))/frame_rate),frame_rate,which_envelope,amplitude);
        
        % Calculate the grating motion
        direction = ceil(2*rand);
        if direction==1
            correct = 'LeftArrow';
            incorrect = 'RightArrow';
        else
            correct = 'RightArrow';
            incorrect = 'LeftArrow';
        end
        % Make noise patch
        variance = 0.16;
        limits = 2;
        
        unfiltered = randn(stimulus_radius*2+ceil(step*mv_length))*sqrt(variance);
        while max(max(abs(unfiltered)))>(limits*sqrt(variance))
            unfiltered = (abs(unfiltered)>(limits*sqrt(variance))).*randn(stimulus_radius*2+ceil(step*mv_length))*sqrt(variance)+(1-(abs(unfiltered)>(limits*sqrt(variance)))).*unfiltered;
        end
        transform = fft2(unfiltered);
        
        % Make filter
        [fx,fy]=meshgrid((-stimulus_radius-ceil(step*mv_length)/2):(stimulus_radius-1+ceil(step*mv_length)/2),(-stimulus_radius-ceil(step*mv_length)/2):(stimulus_radius-1+ceil(step*mv_length)/2));
        angle_map = acosd(abs(fx)./sqrt(fx.^2+fy.^2+eps));
        angle_map(stimulus_radius+1,stimulus_radius+1)=orientation_mean;
        if orientation_profile
            filter_angle = normpdf(angle_map,orientation_mean,orientation_width)/normpdf(orientation_mean,orientation_mean,orientation_width);
        else
            filter_angle = (angle_map<=orientation_width)+(angle_map>=(180-orientation_width));
        end
        sf_map = sqrt(fx.^2+fy.^2)/(stimulus_radius*2)*(60/scale_factor);
        filter_sf = (sf_map>=SF_mean/2^SF_width).*(sf_map<=SF_mean*2^SF_width);
        filtered = ifft2(ifftshift(fftshift(transform).*filter_angle.*filter_sf));
        rms = std(filtered(:));
        %scaled_filtered = real(filtered*contrast_rms/rms+.5);
        scaled_filtered = real(filtered*contrast_rms/rms*2);
        %filtered_output = (filtered*contrast/rms+0.5)*255;
        
        for i=1:mv_length
            moving_grating{i}=max(min((scaled_filtered(1:bps,round((i-1)*step+1):round((i-1)*step+bps))*time_gauss(i).*circle+background),255),0);
            if direction==1
                movie(i) = Screen('MakeTexture',w,moving_grating{i});
            else
                movie(mv_length-(i-1)) = Screen('MakeTexture',w,moving_grating{i});
            end
        end
        
        %Finish the ITI
        WaitSecs(ITI-(GetSecs-aa));
        
        % Draw the black fixation cross
        Screen('DrawLine',w,0,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
        Screen('DrawLine',w,0,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
        Screen('Flip',w);
        
        if ~fast
            FlushEvents('keyDown');
            validKey = 0;
            while ~validKey
                [secs, keyCode, deltaSecs] = KbWait;
                if keyCode(KbName('DownArrow'))&&angle==0
                    validKey = 1;
                elseif keyCode(KbName('LeftArrow'))
                    validKey = 1;
                end
            end
        end
        FlushEvents('keyDown');
        priorityLevel=MaxPriority(w);
        Priority(priorityLevel);
        WaitSecs(0.2);
        
        % Play the movie
        clear StimulusOnsetTime
        for i = 1:mv_length
            Screen('DrawTexture', w, movie(i),movie_rect,screen_patch);
            if fixate
                Screen('DrawLine',w,0,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
                Screen('DrawLine',w,0,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
            end
            [VBLTimestamp StimulusOnsetTime(i) FlipTimestamp Missed Beampos] = Screen('Flip',w);
        end
        Screen('FillRect',w, background);
        Screen('Flip', w);
        
        % Get the response
        validKey = 0;
        while ~validKey
            [secs, keyCode, deltaSecs] = KbWait;
            if keyCode(KbName('ESCAPE'))
                Screen('CloseAll');
                Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
                Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
                ListenChar(1);
                ShowCursor;
                break
            elseif keyCode(KbName(incorrect))
                validKey = 1;
            elseif keyCode(KbName(correct))
                validKey = 1;
                results(perm(trial),count(perm(trial))) = 1;
                if feedback
                    SysBeep;
                end
            else
                Snd('Play',sin(0:1000));
            end
        end
        stim_dir(perm(trial),count(perm(trial)))=direction;
        Priority(0);
        
        for i=1:ceil(mv_length/3)
            Screen('Close',movie(i));
        end
        count(perm(trial)) = count(perm(trial))+1;
    end
    % Shutdown Eyelink:
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    ListenChar(1);
    clc;
    time = toc/60;
    fprintf(1,'\n\n\n');
    fprintf(1,'--------------------------------------------------------\n');
    % Wrap this stuff up.
    for i=1:n_staircases
        fprintf('     -->  %3.0f ms Proportion Correct  =  %5.3f\n',time_sigma(i),mean(results(i,:)));
        fprintf(1,'--------------------------------------------------------\n');
    end
    %fprintf(1,'     -->  Trials with bad timing  =  %3.0f\n',length(find(trialBAD)));
    fprintf(1,'--------------------------------------------------------\n');
    fprintf('     Elapsed time  (minutes)  =  %4.1f\n',time);
    %plot(Q','o- ')
    filename = strcat(data_path,initials,'_',experiment_id,'_',int2str(tme(1)),'_',int2str(tme(2)),'_',int2str(tme(3)),'_',int2str(tme(4)),'_',int2str(tme(5)));
    filename
    clear x y R moving_grating cos2D circle
    save(filename);
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    s = lasterror;
    ddd = psychlasterror;
    msg = ddd.message;
    ListenChar(1);
    ddd = lasterror;
    ddd.message
    ddd.stack(1,1).line
    psychrethrow(lasterror);
    ShowCursor;
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    Priority(0);
    psychrethrow(psychlasterror)
end %try..catch..
