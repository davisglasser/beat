% This program obtains motion discrimination thresholds for
% contrast modulated random noise using interleaved QUEST staircases.
% Davis Glasser
% Last Edited: 10/25/2010

clear all;close all;clc;HideCursor;ListenChar(2);                           % Start with a blank slate

try
    %-----Subject Settings------------------
    % (Change these for every subject)
    initials                        = 'XX';
    starting_point                  = 40;                                  % SD of a temporal gaussian in ms (start high)
    stimulus_radius                 = 480;                                  % Arcmin
    fast                            = 1;                                    % Automatically trigger trials

    %-----Experiment Settings---------------
    % (Keep this consistent across subjects)
    experiment_id                   = 'pilot';                              % Used in group filename
    n_trials                        = 2;                                   % Number of trials per staircase, 25 is good
    n_staircases                    = 2;                                    % Number of staircases, 2 is usually good
    ITI                             = 2;                                    % Intertrial Inverval, Seconds
    
    px_size                         = 2;                                    % arcmin/pix, plz use even multiple of scale factor
    dynamic                         = 1;                                    % 1=new noise each frame.
    SF                              = 1;                                    % Cycles/deg
    TF                              = 4;                                    % Hz
    angle                           = 0;                                    % Degress, 0 = horizontal, 90 = vertical
    modulation_depth                = 10;                              
    contrast                        = 100;                                  % Percent

    background                      = 128;                                  % Grayscale Units
    fixate                          = 0;                                    % Present fixation spot during motion
    % (Bad for fovea, good for periphery)
    feedback                        = 1;                                    
    
    H_ecc_stim                      = 0;                                    % Horizontal stimulus ecc (degs, neg is left)
    H_ecc_fix                       = 0;                                    % Horizontal fixation ecc (degs, neg is left)
    V_ecc_stim                      = 0;                                    % Vertical stimulus ecc (degs, neg is up)
    V_ecc_fix                       = 0;                                    % Vertical fixation ecc (degs, neg is up)

    linearize                       = 0;                                    % Use calibrated LUT (do this when available)
    spatial_envelope                = 2;                                    % 0 = disk, 1 = Gabor, 2 = raised cosine
    which_envelope                  = 15;

    pThreshold                      = 0.82;                                 % Threshold (0.82 works best with Quest)
    tGuessSd                        = .17;                                  % SD of the staircases (bigger = bigger step)
    ceiling                         = 400;                                  % Max duration to allow

    data_path                       = '/ Experiments/data/';                 % Folder for saving data files

    %-----Leave everything below here alone-
    %-----Rig Settings----------------------
    scale_factor                    = 2;                                    % Arcmin/pixel
    frame_rate                      = 60;                                  % Screen frame rate (hz)

    %-----Initialize staircases-------------
    tGuess = starting_point/(1000/frame_rate);
    beta=3.5;delta=0.02;gamma=1/2;
    q = zeros(1,n_staircases);
    for i=1:n_staircases
        Stair(i)=QuestCreate(log10(tGuess),tGuessSd,pThreshold,beta,delta,gamma);
        q(i) = 1;
    end

    %-----Housekeeping----------------------
    % Scale things based on viewing distance, and convert other stuff to
    % the units PsychToolbox wants...
    tme                             = clock;
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
    scale                           = round(px_size/scale_factor);

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
        screen_clut = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
        
%         fid = fopen('GammaTable_FW900A','r');
%         screen_clut = fread(fid,[256 3],'float64');
%         fclose(fid);
%         screen_clut =   screen_clut -1;
%         screen_clut = screen_clut/255;
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

    Screen('DrawText',w,'Duration Threshold: Sine-Wave Gratings',100,100,0);
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
        time_sigma=10^(QuestMean(Stair(perm(trial)))); % get the stimulus duration from Quest
        if time_sigma>ceiling
            time_sigma=ceiling;
        end
        [time_gauss,mv_length] = envelope((1000*(time_sigma)/frame_rate),frame_rate,which_envelope,amplitude);

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
        pedestal = imresize(round(rand(bps,bps))*2-1,scale,'nearest');
        pedestal = pedestal(1:bps,1:bps);
        movie = zeros(1,mv_length);
        for i = 1:mv_length;
            if dynamic
                pedestal = imresize(round(rand(bps,bps))*2-1,scale,'nearest');
                pedestal = pedestal(1:bps,1:bps);
            end
            %grating = (1+sin(a*x+b*y+ motion_step(i)))/2*modulation_depth/100;
            %moving_grating = grating.*pedestal.*circle*time_gauss(i)+background;
            gratinga = (1+sin(a*x+b*y+ motion_step(i)))/2*modulation_depth/100;
            gratingp = gratinga+(1-max(max(gratinga)));
            gratingn = -gratingp;
            moving_grating = (gratingp.*(pedestal==1)+gratingn.*(pedestal==-1)).*circle*time_gauss(i)+background;
            movie(i) = Screen('MakeTexture',w,moving_grating);
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

        WaitSecs(0.2);

        priorityLevel=MaxPriority(w);
        Priority(priorityLevel);
        % Play the movie
        clear StimulusOnsetTime
        for i = 1:mv_length
            Screen('DrawTexture', w, movie(i),movie_rect,screen_patch);
            if fixate
                Screen('DrawLine',w,0,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
                Screen('DrawLine',w,0,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
            end
            Screen('Flip',w);
        end
        Screen('FillRect',w, background);
        Screen('Flip', w);
        Priority(0);

        % Get the response
        FlushEvents('keyDown');
        validKey = 0;
        while ~validKey
            [secs, keyCode, deltaSecs] = KbWait;
            if keyCode(KbName(correct))
                validKey = 1;
                rs = 1;
                if feedback
                    SysBeep;
                end
            elseif keyCode(KbName(incorrect))
                validKey = 1;
                rs = 0;
            else
                Snd('Play',sin(0:1000));
            end
        end

        % Tell Quest what happened
        Stair(perm(trial))=QuestUpdate(Stair(perm(trial)),(log10(time_sigma)),rs);
        Q(perm(trial),q(perm(trial))) = time_sigma*(1000/frame_rate);
        q(perm(trial)) = q(perm(trial))+1;


    end
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
        dur(i)=(10^QuestMean(Stair(i)))*(1000/frame_rate);
        Q(i,q(i)) = dur(i);
        fprintf('     -->  Threshold  %1.0f (ms)  =  %5.3f\n',i,dur(i));
        fprintf(1,'--------------------------------------------------------\n');
    end
    fprintf('     -->  Average        (ms)  =  %5.3f\n',mean(dur));
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
