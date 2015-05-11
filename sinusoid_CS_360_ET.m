% This program tests motion discrimination of luminance-defined
% sine-wave gratings at fixed durations, while recording eyemovements.
% Davis Glasser
% Last Edited: 11/03/2011

clear all;close all;clc;HideCursor;ListenChar(2);                           % Start with a blank slate

try
    %-----Subject Settings------------------
    % (Change these for every subject)
    initials                        = 'XX';
    duration                        = [30 40 50 1000];                           % SD of the temporal Gaussian, ms
    stimulus_radius                 = 480;                                  % Arcmin
    n_trials                        = 2;                                    % Number of trials per duration, 40 is good
    end_pad                         = .2;                                   % in seconds
    fast                            = 1;                                    % Automatically trigger trials

    %-----Experiment Settings---------------
    % (Keep this consistent across subjects)
    experiment_id                   = 'pilot';                              % Used in group filename
    ITI                             = 2;                                    % Intertrial Inverval, Seconds
    SF                              = 1;                                    % Cycles/deg
    TF                              = 4;                                    % Hz
    angle                           = 0;                                    % Degress, 0 = horizontal, 90 = vertical
    contrast                        = 100;                                  % Percent
    background                      = 128;                                  % Grayscale Units
    fixate                          = 0;                                    % Present fixation spot during motion
    % (Bad for fovea, good for periphery)
    H_ecc_stim                      = 0;                                    % Horizontal stimulus ecc (degs, neg is left)
    H_ecc_fix                       = 0;                                    % Horizontal fixation ecc (degs, neg is left)
    V_ecc_stim                      = 0;                                    % Vertical stimulus ecc (degs, neg is up)
    V_ecc_fix                       = 0;                                    % Vertical fixation ecc (degs, neg is up)

    linearize                       = 0;                                    % Use calibrated LUT (do this when available)
    spatial_envelope                = 2;                                    % 0 = disk, 1 = Gabor, 2 = raised cosine
    which_envelope                  = 15;

    data_path                       = '/ Experiments/data/';                 % Folder for saving data files

    %-----Leave everything below here alone-
    %-----Rig Settings----------------------
    scale_factor                    = 2;                                    % Arcmin/pixel
    frame_rate                      = 360;                                  % Screen frame rate (hz)

    %-----Housekeeping----------------------
    % Scale things based on viewing distance, and convert other stuff to
    % the units PsychToolbox wants...
    tme                             = clock;
    eye_filename = strcat(initials,'_',int2str(tme(1)),'_',int2str(tme(2)),'_',int2str(tme(3)),'_',int2str(tme(4)),'_',int2str(tme(5)));
    time_sigma                      = polyval([0.2879/2 0.0325],duration);
    return
    n_staircases                    = length(duration);
    results                         = zeros(n_staircases,n_trials);
    count                           = ones(1,n_staircases);
    stimulus_radius                 = stimulus_radius /scale_factor;
    Gaussian_stdev                  = round(stimulus_radius/1.5);
    f                               = (SF*scale_factor/60)*2*pi;
    TFstep                          = (2*pi*TF)/frame_rate;
    H_ecc_stim                      = H_ecc_stim*60/scale_factor;
    H_ecc_fix                       = H_ecc_fix*60/scale_factor;
    V_ecc_stim                      = V_ecc_stim*60/scale_factor;
    V_ecc_fix                       = V_ecc_fix*60/scale_factor;
    angle                           = angle*pi/180;
    a                               = cos(angle)*f;
    b                               = sin(angle)*f;
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
    %     if size(screens,2) == 2;
    %         w2=Screen('OpenWindow',0,0,[],[],2);
    %         Screen('FillRect',w2, 0); Screen('Flip', w2);
    %     end
    if linearize
        fid = fopen('GammaTable_FW900A','r');
        screen_clut = fread(fid,[256 3],'float64');
        fclose(fid);
        screen_clut =   screen_clut -1;
        screen_clut = screen_clut/255;
        screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    end
    Screen('FillRect',w, background);
    Screen('Flip', w);
    Screen('FillRect',w, background);
    Screen('TextSize',w,20);Screen('TextFont',w,'Charcoal');
    
    el=EyelinkInitDefaults(w);
    
    %-----Init Eyetracker----------
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        sca;
        return;
    end
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

    Screen('DrawText',w,'Motion Discrimination: Sine-Wave Gratings',100,100,0);
    if angle==0
        Screen('DrawText',w,'Use the DOWN arrow to start a trial, and LEFT/RIGHT arrows to respond',100,130,0);
    else
        Screen('DrawText',w,'Use the LEFT arrow to start a trial, and UP/DOWN arrows to respond',100,130,0);
    end
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
    % make sure that we get gaze data from the Eyelink
    Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
    
    % open file to record data to
    Eyelink('openfile', eye_filename);
    
    % STEP 4
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
    eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
    % Set background color to 'backgroundcolor' and do flip to show
    % blank screen:
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
        if angle==0
            if direction==1
                correct = 'LeftArrow';
                incorrect = 'RightArrow';
            else
                correct = 'RightArrow';
                incorrect = 'LeftArrow';
            end
        else
            if direction==1
                correct = 'UpArrow';
                incorrect = 'DownArrow';
            else
                correct = 'DownArrow';
                incorrect = 'UpArrow';
            end
        end
        motion_step = zeros(1,mv_length);
        motion_step(1) = rand*2*pi;
        for i=2:mv_length
            motion_step(i) = motion_step(i-1)+TFstep*((direction-1)*(-2)+1);
        end

        % Make the movie
        movie = zeros(1,mv_length);
        for i = 1:mv_length;
            moving_grating{i} = round(((sin(a*x+b*y+ motion_step(i)).*circle*time_gauss(i))+background));
        end
        frame = zeros(bps,bps,3);
        for i = 1:ceil(mv_length/3)
            for j=1:3
                if ((i-1)*3+j)>mv_length
                    %frame(:,:,j) = ones(bps)*background;
                    switch j
                        case 1
                           frame(:,:,3) = ones(bps)*background;
                        case 2
                            frame(:,:,1) = ones(bps)*background;
                        case 3
                            frame(:,:,2) = ones(bps)*background;
                    end
                else
                    switch j
                        case 1
                           frame(:,:,3) = moving_grating{(i-1)*3+j}; 
                        case 2
                            frame(:,:,1) = moving_grating{(i-1)*3+j};
                        case 3
                            frame(:,:,2) = moving_grating{(i-1)*3+j};
                    end
                end
            end
            movie(i) = Screen('MakeTexture',w,frame);
        end
        % do a final check of calibration using driftcorrection
        EyelinkDoDriftCorrection(el);

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
        Eyelink('StartRecording')
        Eyelink('Message',strcat('TRIAL_',num2str(trial),'_',num2str(duration(perm(trial))),'ms_',num2str(direction)));
        WaitSecs(0.2);

        % Play the movie
        clear StimulusOnsetTime
        Eyelink('Message','ONSET');
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
        Eyelink('Message','OFFSET');
        WaitSecs(end_pad);
        % stop eyelink
        Eyelink('StopRecording');
        Priority(0);

        % Get the response
        validKey = 0;
        while ~validKey
            [secs, keyCode, deltaSecs] = KbWait;
            if keyCode(KbName(correct))
                validKey = 1;
                results(perm(trial),count(perm(trial))) = 1;
                if feedback
                   SysBeep; 
                end
            elseif keyCode(KbName(incorrect))
                validKey = 1;
            else
                Snd('Play',sin(0:1000));
            end
        end

        % Make sure stimulus timing was ok
        clear checker
        trialBAD(trial)=0;
        for i=1:mv_length-1
            checker(trial,i) = round(1000*(StimulusOnsetTime(i+1)-StimulusOnsetTime(i)));
            if checker(trial,i)~=round(1000/frame_rate)
                trialBAD(trial)=1;
            end
        end
        count(perm(trial)) = count(perm(trial))+1;
    end
    % Shutdown Eyelink:
    Eyelink('Shutdown');
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
        fprintf('     -->  %3.0f ms Proportion Correct  =  %5.3f\n',duration(i),mean(results(i,:)));
        fprintf(1,'--------------------------------------------------------\n');
    end
    fprintf(1,'     -->  Trials with bad timing  =  %3.0f\n',length(find(trialBAD)));
    fprintf(1,'--------------------------------------------------------\n');
    fprintf('     Elapsed time  (minutes)  =  %4.1f\n',time);
    plot(Q','o- ')
    filename = strcat(data_path,initials,'_',experiment_id,'_',int2str(tme(1)),'_',int2str(tme(2)),'_',int2str(tme(3)),'_',int2str(tme(4)),'_',int2str(tme(5)));
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
