% This program obtains motion discrimination thresholds for
% contrast modulated random noise using interleaved QUEST staircases.
% Davis Glasser
% Last Edited: 03/07/2011

clear all;close all;clc;HideCursor;ListenChar(2);                           % Start with a blank slate

try
    %-----Subject Settings------------------
    % (Change these for every subject)
    initials                        = 'DT';
    starting_point                  = 25;                                   % SD of a temporal gaussian in ms (start high)
    
    stimulus_duration               = 250;                                  % ms
    ISI                             = 250;
    stimulus_radius                 = 480;                                   % Arcmin
    
    fast                            = 0;                                    % Automatically trigger trials
    
    %-----Experiment Settings---------------
    % (Keep this consistent across subjects)
    experiment_id                   = '2ndOrd_MD';                             % Used in group filename
    n_trials                        = 25;                                   % Number of trials per staircase, 25 is good
    n_staircases                    = 3;                                    % Number of staircases, 2 is usually good
    ITI                             = 2;                                    % Intertrial Inverval, Seconds
    
    px_size                         = 3.5;                                    % arcmin/pix, plz use even multiple of scale factor
    dynamic                         = 1;                                    % 1=new noise each frame.
    SF                              = 1;                                    % Cycles/deg
    TF                              = 4;                                    % Hz
    angle                           = 0;                                    % Degress, 0 = horizontal, 90 = vertical
    contrast                        = 99;                                  % Percent
    supersupressor                  = 0;
    background                      = 128;                                  % Grayscale Units
    fixate                          = 0;                                    % Present fixation spot during motion
    % (Bad for fovea, good for periphery)
    feedback                        = 1;
    
    H_ecc_stim                      = 0;                                    % Horizontal stimulus ecc (degs, neg is left)
    H_ecc_fix                       = 0;                                    % Horizontal fixation ecc (degs, neg is left)
    V_ecc_stim                      = 0;                                    % Vertical stimulus ecc (degs, neg is up)
    V_ecc_fix                       = 0;                                    % Vertical fixation ecc (degs, neg is up)
    
    linearize                       = 1;                                    % Use calibrated LUT (do this when available)
    spatial_envelope                = 2;                                    % 0 = disk, 1 = Gabor, 2 = raised cosine
    which_envelope                  = 15;
    
    pThreshold                      = 0.82;                                 % Threshold (0.82 works best with Quest)
    tGuessSd                        = .2;                                  % SD of the staircases (bigger = bigger step)
    ceiling                         = 100;                                  % Max duration to allow
    
    data_path                       = '/Users/tadinlab/Desktop/DGProjector/Data/';                 % Folder for saving data files
    
    %-----Leave everything below here alone-
    %-----Rig Settings----------------------
    scale_factor                    = 1.75;                                    % Arcmin/pixel
    frame_rate                      = 360;                                  % Screen frame rate (hz)
    
    %-----Initialize staircases-------------
    tGuess = starting_point/100;
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
    time_sigma                      = polyval([0.2879/2 0.0325],stimulus_duration);
    stimulus_radius                 = round(stimulus_radius/scale_factor);
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
    bps = round((stimulus_radius)*2+1);
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
    %         if size(screens,2) == 2;
    %             w2=Screen('OpenWindow',0,0,[],[],2);
    %             Screen('FillRect',w2, 0); Screen('Flip', w2);
    %         end
    if linearize
        screen_clut = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
        
        %         fid = fopen('GammaTable_FW900A','r');
        %         screen_clut = fread(fid,[256 3],'float64');
        %         fclose(fid);
        %         screen_clut =   screen_clut -1;
        %         screen_clut = screen_clut/255;
        Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
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
    
    %Make all the noise
    noise = zeros(bps,bps,360);
    pedestal = imresize(round(rand(bps,bps))*2-1,scale,'nearest');
    pos_range = size(pedestal)-(bps+1);
    for i=1:360
        noise_pos = ceil(pos_range.*rand(1,2));
        noise(:,:,i) = pedestal(noise_pos(1):bps-1+noise_pos(1),noise_pos(2):bps-1+noise_pos(2));
    end
    
    %Premake the noise movie
    [time_gauss,mv_length] = envelope((1000*(time_sigma)/frame_rate),frame_rate,which_envelope,amplitude);
    noise_movie = zeros(1,ceil(mv_length/3));
    frame = zeros(bps,bps,3);
    for i = 1:ceil(mv_length/3)
        for j=1:3
            moving_grating = noise(:,:,i).*circle*time_gauss(i)+background;
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
                        %frame(:,:,3) = moving_grating{(i-1)*3+j};
                        frame(:,:,3) = moving_grating;
                    case 2
                        frame(:,:,1) = moving_grating;
                    case 3
                        frame(:,:,2) = moving_grating;
                end
            end
        end
        noise_movie(i) = Screen('MakeTexture',w,frame);
    end
    
    Screen('DrawText',w,'Modulation Depth Threshold: Second-Order Gratings',100,100,0);
    if angle==0
        Screen('DrawText',w,'Use the DOWN arrow to start a trial, and LEFT (INT=1)/RIGHT (INT-2) arrows',100,130,0);
    else
        Screen('DrawText',w,'Use the LEFT arrow to start a trial, and UP/DOWN arrows to respond',100,130,0);
    end
    Screen('DrawText',w,[int2str(total_trials),'  trials - press SPACE BAR to start'],100,160,0);
    Screen('Flip',w);
    
    FlushEvents('keyDown');
    validKey = 0;
    while ~validKey
        [secs, keyCode, deltaSecs] = KbWait(-1);
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
        
        % Calculate the modulation depth
        log_modulation_depth=QuestMean(Stair(perm(trial))); % get the stimulus duration from Quest
        modulation_depth = 100*10^log_modulation_depth;
        
        if modulation_depth>ceiling
            modulation_depth=ceiling;
        end
        
        % Calculate the grating motion
        
        direction = ceil(2*rand);
        interval = ceil(2*rand);
        if interval==1
            correct = 'LeftArrow';
            incorrect = 'RightArrow';
        else
            correct = 'RightArrow';
            incorrect = 'LeftArrow';
        end
        
        motion_step = zeros(1,mv_length);
        motion_step(1) = rand*2*pi;
        for i=2:mv_length
            motion_step(i) = motion_step(i-1)+TFstep*((direction-1)*(-2)+1);
        end
        movie = zeros(1,ceil(mv_length/3));
        frame = zeros(bps,bps,3);
        for i = 1:ceil(mv_length/3)
            for j=1:3
                gratinga = (1+sin(a*x+b*y+ motion_step(i)))/2*modulation_depth/100;
                gratingp = gratinga+(1-max(max(gratinga)));
                gratingn = -gratingp;
                if dynamic
                    moving_grating = (gratingp.*(noise(:,:,i)==1)+(gratingn.*(noise(:,:,i)==-1))).*circle*time_gauss(i)+background;
                else
                    moving_grating = (gratingp.*(noise(:,:,1)==1)+(gratingn.*(noise(:,:,1)==-1))).*circle*time_gauss(i)+background;
                end
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
                            %frame(:,:,3) = moving_grating{(i-1)*3+j};
                            frame(:,:,3) = moving_grating;
                        case 2
                            frame(:,:,1) = moving_grating;
                        case 3
                            frame(:,:,2) = moving_grating;
                    end
                end
            end
            movie(i) = Screen('MakeTexture',w,frame);
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
        Screen('FillRect',w, background);
        Screen('Flip', w);
        WaitSecs(0.2);
        
        priorityLevel=MaxPriority(w);
        Priority(priorityLevel);
        % Play the movie
        if interval==1
            for i = 1:ceil(mv_length/3)
                Screen('DrawTexture', w, movie(i),movie_rect,screen_patch);
                if fixate
                    Screen('DrawLine',w,0,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
                    Screen('DrawLine',w,0,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
                end
                Screen('Flip',w);
            end
            Screen('FillRect',w, background);
            Screen('Flip', w);
            WaitSecs(ISI/1000);
            for i = 1:ceil(mv_length/3)
                Screen('DrawTexture', w, nosie_movie(i),movie_rect,screen_patch);
                if fixate
                    Screen('DrawLine',w,0,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
                    Screen('DrawLine',w,0,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
                end
                Screen('Flip',w);
            end
            Screen('FillRect',w, background);
            Screen('Flip', w);
        else
            for i = 1:ceil(mv_length/3)
                Screen('DrawTexture', w, nosie_movie(i),movie_rect,screen_patch);
                if fixate
                    Screen('DrawLine',w,0,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
                    Screen('DrawLine',w,0,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
                end
                Screen('Flip',w);
            end
            Screen('FillRect',w, background);
            Screen('Flip', w);
            WaitSecs(ISI/1000);
            for i = 1:ceil(mv_length/3)
                Screen('DrawTexture', w, movie(i),movie_rect,screen_patch);
                if fixate
                    Screen('DrawLine',w,0,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
                    Screen('DrawLine',w,0,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
                end
                Screen('Flip',w);
            end
            Screen('FillRect',w, background);
            Screen('Flip', w);
        end
        
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
        Stair(perm(trial))=QuestUpdate(Stair(perm(trial)),log_modulation_depth,rs);
        Q(perm(trial),q(perm(trial))) = modulation_depth;
        q(perm(trial)) = q(perm(trial))+1;
        
        for i=1:ceil(mv_length/3)
            Screen('Close',movie(i));
        end
        
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
        dur(i)= 100*10^QuestMean(Stair(i));
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
