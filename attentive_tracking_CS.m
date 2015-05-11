% This task requires observers to attentively track an indicated target in
% a circular array. At the end they indicate if a cued target is the
% original target. It follows the basic task from Holcombe et al., 2011,
% Current Biology, reasonably closely.
%
% Doesn't currently monitor fixation, but that would be a good idea in the
% future.
%
% Davis Glasser
% Last Edited: 12/11/2012
clc; clear; HideCursor; ListenChar(2);
try
    %-----Subject Settings------------------
    initials                        = 'test';
    speeds                          = [.5];           % in RPS, 0.5-2 good range, avg ~1.4
    n_trials                        = 10;                                   % Number of trials per duration, 40 is good
    fast                            = 0;                                    % Automatically trigger trials/0 probably best for this task
    feedback                        = 0;
    fixate                          = 1;
    %-----Experiment Settings---------------
    array_radius                    = 4;                                    % Degrees
    element_radius                  = 1.5;
    element_color                   = 0;
    cue_color                       = 255;
    n_elements                      = 4;                                    % even numbers only, please
    cue_duration                    = 400;                                  % ms
    trans_duration                  = 300;
    track_duration                  = 2700;                                 % ms
    
    background                      = 128;
    
    ITI                             = 2;                                    % Intertrial Inverval, Seconds
    
    H_ecc_stim                      = 0;                                    % Horizontal stimulus ecc (degs, neg is left)
    H_ecc_fix                       = 0;                                    % Horizontal fixation ecc (degs, neg is left)
    V_ecc_stim                      = 0;                                    % Vertical stimulus ecc (degs, neg is up)
    V_ecc_fix                       = 0;                                    % Vertical fixation ecc (degs, neg is up)
    
    spatial_envelope                = 2;
    linearize                       = 1;                                    % Use calibrated LUT (do this when available)
    experiment_id                   = 'attn_track';                         % Used in group filename
    data_path                       = '/Users/tadinlab/Desktop/DGProjector/Data/';                 % Folder for saving data files
    
    
    %-----Leave everything below here alone-
    %-----Rig Settings----------------------
    scale_factor                    = 2;                                    % Arcmin/pixel
    frame_rate                      = 60;                                  % Screen frame rate (hz)
    
    %-----Housekeeping----------------------
    % Scale things based on viewing distance, and convert other stuff to
    % the units PsychToolbox wants...
    tme                             = clock;
    n_staircases                    = length(speeds);
    results                         = zeros(n_staircases,n_trials);
    count                           = ones(1,n_staircases);
    array_radius                    = array_radius/scale_factor*60;
    stimulus_radius                 = round(element_radius/scale_factor*60);
    frame_step                      = speeds*360/frame_rate;
    cue_mv_length                   = round(cue_duration/(1000/frame_rate));
    trans_mv_length                 = round(trans_duration/(1000/frame_rate));
    mv_length                       = round(track_duration/(1000/frame_rate));
    H_ecc_stim                      = H_ecc_stim*60/scale_factor;
    H_ecc_fix                       = H_ecc_fix*60/scale_factor;
    V_ecc_stim                      = V_ecc_stim*60/scale_factor;
    V_ecc_fix                       = V_ecc_fix*60/scale_factor;
    
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
    base_circle = circle;
    if spatial_envelope == 1
        circle = (exp(-(((x)/(sqrt(2)*Gaussian_stdev/6)).^2)-((y/(sqrt(2)*Gaussian_stdev/2)).^2)).*circle);
    elseif spatial_envelope == 2
        R = (sqrt(x.^2 + y.^2) + eps).*circle;
        R = R/max(max(R));
        cos2D = (cos(R*pi)+1)/2;
        circle = (cos2D.*circle);
    end
    [x2,y2]=meshgrid((-stimulus_radius/2):(stimulus_radius/2),(-stimulus_radius/2):(stimulus_radius/2));
    circle2=((stimulus_radius)^2-(x2.^2+y2.^2));
    for i=1:bps/2; for j =1:bps/2; if circle2(i,j) < 0; circle2(i,j) = 0; else circle2(i,j) = 1; end; end;
    end;
    
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
        Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    end
    Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Screen('FillRect',w, background);
    Screen('Flip', w);
    Screen('FillRect',w, background);
    Screen('TextSize',w,20);Screen('TextFont',w,'Charcoal');
    
    %-----Stimulus Rectangles---------------
    scr_left_middle = fix(screen_rect(3)/2)-round(bps/2);
    scr_top = fix(screen_rect(4)/2)-round(bps/2);
    sr_hor = round(screen_rect(3)/2);
    sr_ver = round(screen_rect(4)/2);
    fix_hor = sr_hor+H_ecc_fix;
    fix_ver = sr_ver+V_ecc_fix;
    
    %-----Make the Dot Frames
    big_frame = zeros(2*array_radius+bps+10);
    movie_rect= [0,0,size(big_frame,1),size(big_frame,2)];
    screen_patch = CenterRectOnPoint(movie_rect,sr_hor,sr_ver);
    bf_center = round(size(big_frame,1)/2);
    r = ones(bps)*element_color-background;
    dot = r.*circle;
    dot_rect = [0 0 bps bps];
    dot_ang = (0:(n_elements-1))*360/n_elements;
    for i=1:n_elements
        elx = sind(dot_ang(i))*array_radius;
        ely = cosd(dot_ang(i))*array_radius;
        rectdot = CenterRectOnPoint(dot_rect,elx,ely);
        big_frame(round((rectdot(1)+bf_center):(rectdot(3)+bf_center-1)),round((rectdot(2)+bf_center):(rectdot(4)+bf_center-1)))=min(cat(3,big_frame(round((rectdot(1)+bf_center):(rectdot(3)+bf_center-1)),round((rectdot(2)+bf_center):(rectdot(4)+bf_center-1))),dot),[],3 );
    end
    t_rect = big_frame;
    cue = ones(ceil(bps/2))*(cue_color-background).*circle2;
    cue_rect = [0 0 size(cue,2) size(cue,2) ];
    elx = sind(0)*array_radius;
    ely = cosd(0)*array_radius;
    rectcue = CenterRectOnPoint(cue_rect,elx,ely);
    t_rect(round((rectcue(1)+bf_center):(rectcue(3)+bf_center-1)),round((rectcue(2)+bf_center):(rectcue(4)+bf_center-1)))=max(cat(3,t_rect(round((rectcue(1)+bf_center):(rectcue(3)+bf_center-1)),round((rectcue(2)+bf_center):(rectcue(4)+bf_center-1))) ,cue),[],3);
    d_rect = flipud(t_rect);
    sca;
    return;
    
    frame = Screen('MakeTexture',w,big_frame+background);
    t_frame = Screen('MakeTexture',w,t_rect+background);
    d_frame = Screen('MakeTexture',w,d_rect+background);
    
    Screen('DrawText',w,'Attentional Tracking',100,100,0);
    Screen('DrawText',w,'Use the LEFT arrow if the indicated dot was the precued dot.',100,130,0);
    Screen('DrawText',w,'Use the RIGHT arrow if the indicated dot is a different dot.',100,160,0);
    Screen('DrawText',w,[int2str(total_trials),'  trials - press SPACE BAR to start'],100,190,0);
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
    for trial=1:total_trials
        aa = GetSecs;
        % Draw the white fixation cross
        Screen('FillRect',w, background);
        Screen('DrawLine',w,255,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
        Screen('DrawLine',w,255,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
        Screen('Flip',w);
        
        % Calculate the grating motion
        direction = round(rand)*2-1;
        same = round(rand);
        if same
            correct = 'LeftArrow';
            incorrect = 'RightArrow';
        else
            correct = 'RightArrow';
            incorrect = 'LeftArrow';
        end
        
        phase = rand*360;
        
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
                if keyCode(KbName('DownArrow'))
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
        for i = 1:cue_mv_length
            Screen('DrawTexture', w, frame,movie_rect,screen_patch,phase);
            if fixate
                Screen('DrawLine',w,0,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
                Screen('DrawLine',w,0,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
            end
            Screen('DrawDots', w, [sr_hor+ array_radius*sind(90-phase);sr_ver+array_radius*cosd(90-phase)],30,255,[ ],2);
            Screen('Flip', w);
            phase = phase+direction*frame_step(perm(trial));
        end
        for i = 1:trans_mv_length
            Screen('DrawTexture', w, frame,movie_rect,screen_patch,phase);
            if fixate
                Screen('DrawLine',w,0,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
                Screen('DrawLine',w,0,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
            end
            Screen('DrawDots', w, [sr_hor+ array_radius*sind(90-phase);sr_ver+array_radius*cosd(90-phase)],30,[255;255;255;(trans_mv_length-i)/trans_mv_length*255],[ ],2);
            Screen('Flip', w);
            phase = phase+direction*frame_step(perm(trial));
        end
        for i = 1:mv_length
            Screen('DrawTexture', w, frame,movie_rect,screen_patch,phase);
            if fixate
                Screen('DrawLine',w,0,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
                Screen('DrawLine',w,0,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
            end
            Screen('Flip', w);
            phase = phase+direction*frame_step(perm(trial));
        end
        % Get the response
        validKey = 0;
        while ~validKey
            Screen('DrawTexture', w, frame,movie_rect,screen_patch,phase);
            if fixate
                Screen('DrawLine',w,0,fix_hor-2,fix_ver,fix_hor+2,fix_ver,2);
                Screen('DrawLine',w,0,fix_hor,fix_ver-2,fix_hor,fix_ver+2,2);
            end
            if same
                Screen('DrawDots', w, [sr_hor+ array_radius*sind(90-phase);sr_ver+array_radius*cosd(90-phase)],30,255,[ ],2);
            else
                Screen('DrawDots', w, [sr_hor+ array_radius*sind(270-phase);sr_ver+array_radius*cosd(270-phase)],30,255,[ ],2);
            end
            Screen('Flip', w);
            phase = phase+direction*frame_step(perm(trial));
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck;
            if keyIsDown
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
        end
        Priority(0);
        count(perm(trial)) = count(perm(trial))+1;
    end
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    ListenChar(1);
    ShowCursor;
    clc;
    time = toc/60;
    fprintf(1,'\n\n\n');
    fprintf(1,'--------------------------------------------------------\n');
    % Wrap this stuff up.
    for i=1:n_staircases
        fprintf('     -->  %3.1f rps, Proportion Correct  =  %5.3f\n',speeds(i),mean(results(i,:)));
        fprintf(1,'--------------------------------------------------------\n');
    end
    fprintf(1,'--------------------------------------------------------\n');
    fprintf('     Elapsed time  (minutes)  =  %4.1f\n',time);
    %plot(Q','o- ')
    filename = strcat(data_path,initials,'_',experiment_id,'_',int2str(tme(1)),'_',int2str(tme(2)),'_',int2str(tme(3)),'_',int2str(tme(4)),'_',int2str(tme(5)));
    clear x y R moving_grating cos2D circle
    save(filename);
catch ME
    ShowCursor;
    ListenChar(1);
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    Priority(0);
    ME.stack
    ME.message
end