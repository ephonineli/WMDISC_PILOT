% To run: e.g., DACOM1_WMDISC_PILOT('1',1,1) % the first arugment needs to be a string
function WMDISC_PILOT(subName,CurRun,ispractice)

try
    close all;
    clear mex;
    IOPort('Closeall');
    Screen('Closeall');
    AssertOpenGL;
    %test script
    Screen('Preference', 'SkipSyncTests', 1);    
    screenid = max(Screen('Screens'));
    [win,rect] = Screen('OpenWindow',screenid);
    %rect(3) and rect(4) are the screen resolutions
    
    HideCursor;

    Screen('TextFont', win, 'Arial');
    HideCursor;
    white=255;
    black=0;
    gray=(white+black)/2;
    red=[255,0,0];
    green=[0,255,0];
        
    %key buttons    
    KbName('UnifyKeyNames');
    escapeKey = KbName('ESCAPE');
    ltiltKey  = KbName('E');   
    rtiltKey  = KbName('F');   
    leftKey = KbName('LeftArrow');
    rightKey = KbName('RightArrow');
    confirmKey = KbName('UpArrow');
    moveonKey = KbName('space'); 
    xc = rect(3)/2;
    yc = rect(4)/2;
    
    repeatmissedtrials = 0; % flag to repeat missed trials
    radius = 9;
    
    cue_dur=0.2;
    sample_dur = 0.5;
    ISI_dur = 0.5;
    delay_dur = 2.5;
    disc_dur = 0.2;
    TJ_dur = 3.8;
    recall_dur = 4; 
    FB_dur = 0.2;
    ITI_dur_mat=[6,8,10];
    trialnumber=30;% 30 trials, 12 blocks, 360 in total
    if ispractice==1
        trialnumber=8;
    end
    
    fndate = datestr(now, 'yyyymmdd_HHMMSS');
    behavDir = fullfile('data','behavior');
    if ~exist(behavDir,'dir')
        mkdir(behavDir);
    end    
    filename = fullfile(behavDir, sprintf('DACOM1_WMDISC_%s_%d_%s.mat', subName, CurRun, fndate));
    if ispractice==1
        filename= fullfile(behavDir, sprintf('practice_DACOM1_WMDISC_%s_%d_%s.mat', subName, CurRun, fndate));
    end    
     %% visual angle and stimuli sizes
    Distance=50;%cm % EEG room = 53; % 166 = 50
    ScreenWidth=51;%cm % EEG room = 60; % 166 = 51
    ScreenHeight=27.5;%cm % EEG room % 166 = 27.5
    
    %Resolution_IBP=[1200,1920];
    Resolution_IBP=[rect(3), rect(4)]; % rect(3) = width rect(4) = height
    PixelSize_IBP=ScreenWidth/Resolution_IBP(1); %cm per pixel
    degree_meter_IBP=2*Distance*tand(1/2); % cm per degree
    pixel_per_degree=degree_meter_IBP/PixelSize_IBP; % pixels per degree
    
    radius_in_pixel = radius * pixel_per_degree; % convert radius (in degrees) to pixels
    
    % spatial_frequency_in_pixel = spatial_frequency / pixel_per_degree;
    
    size_in_pixel = 2 * radius * pixel_per_degree; 
       
    width = Resolution_IBP(1);
    height = Resolution_IBP(2);
    
    %Eccentricity
    Eccentricity_degree=10;
    Eccentricity_meter_IBP=2*Distance*tand(Eccentricity_degree/2); % convert 'Eccentricity_degree'degree to pixel
    Eccentricity_pixel=Eccentricity_meter_IBP/PixelSize_IBP; % convert cm distance to pixels
    Offset_peri=round(Eccentricity_pixel);
    
    %% -----setting up stimuli--------
    memory_ori = [75 45 15 165 135 105];
    memory_ori_rad = memory_ori*pi/180;
    FixationRect=CenterRect([0,0,12,12],rect);
    FixationRect=OffsetRect(FixationRect,0,0);
    
    wheelrect = CenterRect([0 0 1.5*radius_in_pixel  1.5*radius_in_pixel], rect);
    wheelCenterRect = CenterRectOnPoint(wheelrect, xc, yc);   % centred wheel frame
    xwheel_center = xc;                                     % bar’s x-origin    
    ywheel=yc;
    
    %--------------------s
    % Gabor infor2131mation
    %--------------------

    % Dimension of the region where will draw the Gabor in pixels
    gaborDimPix = round(size_in_pixel);%rect(4) / 2;
    sigma = gaborDimPix / 7;
    % Obvious Parameters
    contrast=0.8;
    aspectRatio = 1.0;
    phase = 0;

    % Spatial Frequency (Cycles Per Pixel)
    % One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
    % gaborDimPix = diameter of gabor (in pixels) // gaborDimPix/pixel_per_degree =
    % stimulus width in degrees 
    % since PTB takes the frequency as cycles per pixel instead of cycles
    % per degree, need to convert the desired cycles per degree to
    % corresponding cycles per pixel

    % using cycles per degree
    % cycles per degree = (# cycles in patch)/(patch diameter) % in
    % Teng2025 it seems ~15.54 per patch with radius of 6 (diameter=12)
    % with r = 9 deg , 15.54/18 ~= (cycles per degree)
    % spatial_frequency = 15.54/(2*radius); % desired cycles per degree
    % freq = spatial_frequency/pixel_per_degree; % cycles per pixel
    
    numCycles = 5;
    freq = numCycles / gaborDimPix;
    backgroundOffset = [0.5 0.5 0.5 0];
    disableNorm = 1;
    preContrastMultiplier = 0.5;
    [gabortex,gabortex_rect] = CreateProceduralGabor(win, gaborDimPix, gaborDimPix, [],...
        backgroundOffset, disableNorm, preContrastMultiplier);
    [gabortex2,gabortex_rect2] = CreateProceduralGabor(win, gaborDimPix, gaborDimPix, [],...
        backgroundOffset, disableNorm, preContrastMultiplier);
    
    center_posrect = CenterRectOnPoint(gabortex_rect, xc, yc);   % JH: this part is added to designate destination for central Gabors    
% 
    % Randomise the phase of the Gabors and make a properties matrix.
    propertiesMat = [phase, freq, sigma, contrast, aspectRatio, 0, 0, 0];    
    %% -----condition matrix-----------
    if CurRun == 1
        mem_cue = [1,2]; %1 = S1, 2 = S2
        %mem_disc_distance=[-72:18:-18,18:18:90];
        mem_disc_distance=[-60 -30 30 60 90]; % possible offsets when congruent??
       
        
        contrast_cond = [0.8]; % contrast levels for the discriminandum
        ori_labels = 1:6; % counterbalance all orientation type per block

        trialcond = []; % will grow to (#reps (18) × #cong (2) × #cue (2) × #contrast (1) × #distance (5)) rows × 4 cols
        cond_temp=[]; % temporary storage per (rep, congruency) block
        for irep = 1:18 % repeat each full set 18 times
            for icong = 1:2 %2 % 1 = “incongruent” block, 2 = “congruent” block
                cond_temp = []; % reset at start of each congruency block
                for icue = 1:length(mem_cue)%2 % loop over each cue (1 and 2)
                    for icontrast = 1:length(contrast_cond)%1 each contrast level
                        for idistance=1:length(mem_disc_distance)%5 % each distance                                   
                            % icong=1 builds incongruent trials and icong=2
                            % for congruent trials. So icong-1 = 0 and icong-2 =1 
                            % builds congruent trials. Later below, in:
                            % run_ori_cong = conditions(:,1), if
                            % run_ori_cong(trial) =1 (congruent), disc_orientation = memory_orientation;
                            % else, disc_orientation = memory_orientation + run_distance(trial);
                            if icong==1 % offset cycles +/-30 +/-60 +90 for congruent trials
                                temp = [icong-1,mem_cue(icue),contrast_cond(icontrast),mem_disc_distance(idistance)];
                            else % offset is forced to 0 for congruent
                                temp = [icong-1,mem_cue(icue),contrast_cond(icontrast),0];
                            end
                            cond_temp = [cond_temp;temp]; % append this one row into the per-congruency block
                        end
                    end
                end

                % now that cond_temp holds all (cue × contrast × distance) rows
                % for one icong, append it into the master list
                trialcond = [trialcond;cond_temp];
            end
        end

        trialcond = Shuffle(trialcond,2); % randomize trial order across *all* reps and conditions
        %1-oricong 2-cue type 3-discriminandum_contrast 4-mem_disc_distance

        ori_labels = (1:6)';                     % column vector [1;2;3;4;5;6]
        nSet = size(trialcond,1) / numel(ori_labels); 
        nRepeats = repmat(ori_labels, nSet, 1);  % each ori appears nSet times
        ITI_dur_repeat = repmat(ITI_dur_mat,size(trialcond,1)/numel(ITI_dur_mat),1);
        ori1 = Shuffle(nRepeats);           
        ori2 = Shuffle(nRepeats);                % independent shuffle
        ITI_durs = Shuffle(ITI_dur_repeat(:));

        % append them as columns 5 & 6 & 7
        trialcond = [trialcond, ori1, ori2,ITI_durs];
        % trialcond = Shuffle(trialcond, 2);        
        save(fullfile('data', ['parameters_WMDISC_' num2str(subName) '.mat']), 'trialcond');
    else
        load(fullfile('data', ['parameters_WMDISC_' num2str(subName) '.mat']));
    end
     
    %% ------show instructions--------
    Screen('TextSize',win,30);
    if ispractice==1
        instructions=['This will be a short practice block. \n\n\n\n Press any button to start.'];   
    else
        instructions=['In this study please always fixate at the center and remember the orientation of the cued item. \n\n During the memory delay, make a left/right judgement.\n\n Press any button to start.'];
    end
    instructions=WrapString(instructions,40);
    Screen(win, 'FillRect', gray);
    %Screen('FillOval', win, white, FixationRect);
    DrawFormattedText(win, instructions, 'center', 'center', black);
    Screen('Flip',win);
    
    KbWait;
    
    %% ------exp start-------------
     
        HideCursor;
        Screen(win, 'FillRect', gray);
        Screen('FillOval', win, white, FixationRect);
        Screen('Flip',win);
        WaitSecs(3);
        
        disc_acc = nan(trialnumber,1);
        disc_resp = nan(trialnumber,1);
        tj_rt = nan(trialnumber,1);
        recall_rt = nan(trialnumber,1);
        recall_rt2 = nan(trialnumber,1); % this rt starts measuring time from the first key press
        wm_error = nan(trialnumber,1);
        target_ori = nan(trialnumber,1);
        probe_ori = nan(trialnumber,1);
        gotRespTJ = nan(trialnumber,1);
        gotRespRecall = nan(trialnumber,1);
        record = [];
        t0=GetSecs;
        
        conditions=trialcond(trialnumber*(CurRun-1)+1:trialnumber*CurRun,:);
        %1-oricong 2-cue type 3-discriminandum_contrast 4-mem_disc_distance
        %5-ori1label 6-ori2label 7-ITI duration
        run_ori_cong=conditions(:,1);
        run_cue_type=conditions(:,2);
        run_disc_contrast=conditions(:,3);
        run_distance=conditions(:,4);
        
        noresptrials = [];
        % initiate empty array in which to save trial numbers in which a
        % response was not recorded
        noresptrialssave = [];
        runningtrialssave = [];        
        finishflag = 0;
        trial = 1;
        % columns: [missed_TJ, missed_recall]
        missed_trials = nan(trialnumber, 2);
    while trial<=trialnumber
        trial_t0 = GetSecs; 
        ori1_label = conditions(trial,5);
        ori2_label = conditions(trial,6);   
        
        orientation1 = memory_ori(ori1_label);
        orientation1 = orientation1 + 3*rand()*(-1)^round(rand());
        orientation2 = memory_ori(ori2_label);
        orientation2 = orientation2 + 3*rand()*(-1)^round(rand());

        ITI_dur = conditions(trial,7);
        
        trial_cue_type = run_cue_type(trial);
        runningtrialssave(end+1) = trial;        

        if trial_cue_type == 1          % S1 is target
            memory_orientation = orientation1;
            cueLabel = '1';
        else                       % S2 is target
            memory_orientation = orientation2;   
            cueLabel = '2';
        end        

        target_ori(trial) = memory_orientation;

        
        if run_ori_cong(trial)==1
            disc_orientation = orientation2; % set disc to be cong/incong wrt ori2
        else
            disc_orientation = orientation2+run_distance(trial);
        end
        % wrap discriminandum orientation into [–90, +90] for computing
        % left/right tilt later
        if disc_orientation<=-90
            disc_orientation=disc_orientation+180;
        elseif disc_orientation>90
            disc_orientation=disc_orientation-180;
        end
        
        run_disc_orientation(trial)=disc_orientation;
       

        contrast = 0.8;
        propertiesMat = [phase, freq, sigma, contrast, aspectRatio, 0, 0, 0];              
% ------- S1 -----------
        % contrast=80;
        % noiselevel=0;
        % [m1]=NoisyGabor(orientation1,size_in_pixel,contrast,noiselevel);
        % sample1=Screen('MakeTexture', win, m1);       
        while (GetSecs - trial_t0) < sample_dur
        % while (GetSecs-t0)<trial*sample_dur+(trial-1)*(ISI_dur+sample_dur+delay_dur+disc_dur+TJ_dur+recall_dur+FB_dur+ITI_dur)
           % Screen('DrawTexture', win, sample1, [], center_posrect);
            Screen('DrawTexture', win, gabortex, gabortex_rect, center_posrect, -orientation1, [], [], [], [],...
                kPsychDontDoRotation, propertiesMat');
            Screen('Flip', win);
        end
        
        
% ------- ISI -----------        
       while (GetSecs - trial_t0) < (sample_dur + ISI_dur)
       % while (GetSecs-t0)<trial*(sample_dur+ISI_dur)+(trial-1)*(sample_dur+delay_dur+disc_dur+TJ_dur+recall_dur+FB_dur+ITI_dur)
            Screen(win, 'FillRect', gray);
            Screen('FillOval', win, white, FixationRect);
            Screen('Flip',win);
        end        
% ------- S2 -----------
        % contrast=80;
        % noiselevel=0;
        % [m2]=NoisyGabor(orientation2,size_in_pixel,contrast,noiselevel);
        % sample2=Screen('MakeTexture', win, m2);

        while (GetSecs - trial_t0) < (sample_dur + ISI_dur + sample_dur)
        % while (GetSecs-t0)<trial*(sample_dur+ISI_dur+sample_dur)+(trial-1)*(delay_dur+disc_dur+TJ_dur+recall_dur+FB_dur+ITI_dur)
            % Screen('DrawTexture', win, sample2, [], center_posrect);
            Screen('DrawTexture', win, gabortex2, gabortex_rect2, center_posrect, -orientation2, [], [], [], [],...
             kPsychDontDoRotation, propertiesMat');
            Screen('Flip', win);
        end
%------- delay 1.1 ------------------
        while (GetSecs - trial_t0) < (sample_dur + ISI_dur + sample_dur + delay_dur)
        % while (GetSecs-t0)<trial*(sample_dur+ISI_dur+sample_dur+delay_dur)+(trial-1)*(disc_dur+TJ_dur+recall_dur+FB_dur+ITI_dur)
            Screen(win, 'FillRect', gray);
            Screen('FillOval', win, white, FixationRect);
            Screen('Flip',win);
        end        
        
% -----Discrimination stimuli---------       
        % % contrast=run_disc_contrast(trial)*100;
        noiselevel=0.2;
        % [d1]=NoisyGabor(disc_orientation,size_in_pixel,contrast,noiselevel);
        % noise=Screen('MakeTexture', win, d1);

        if run_disc_orientation(trial) < 0 
            disc_side(trial) = 1; %left
        else
            disc_side(trial) = 2; %right
        end    
        
        % Precompute a static noise patch for this cue 
        H = gaborDimPix; W = gaborDimPix;
        mask = rand(H,W) < noiselevel;       % Bernoulli(p): which pixels get replaced
        uni  = uint8(255*rand(H,W));           % Uniform(0..255) replacement values
        A    = uint8(mask)*255;                % alpha=255 where we replace, 0 elsewhere
        rgba = cat(3, uni, uni, uni, A);         % RGBA image (noise RGB, mask in A)
        noiseTexCdisc = Screen('MakeTexture', win, rgba);  

        while (GetSecs - trial_t0) < (sample_dur + ISI_dur + sample_dur + delay_dur + disc_dur)
        % while (GetSecs-t0)<trial*(sample_dur+ISI_dur+sample_dur+delay_dur+disc_dur)+(trial-1)*(TJ_dur+recall_dur+FB_dur+ITI_dur)
            Screen('DrawTexture', win, gabortex2, gabortex_rect2, center_posrect, -disc_orientation, [], [], [], [],...
             kPsychDontDoRotation, propertiesMat');

            % overlay the replacement noise ONLY at the masked pixels
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);   % Enable standard alpha blending so Gabor's Gaussian aperture is transparent outside
            Screen('DrawTexture', win, noiseTexCdisc, [], center_posrect);
            Screen('BlendFunction', win, GL_ONE, GL_ZERO);
            
            % Screen('DrawTexture', win, noise,[], center_posrect);
            
            Screen('Flip', win);
        end
        
%---------Tilt Judgement--------------
        count = 1;  
        t1 = GetSecs;      
        gotResponse = 0;
        Screen('FillOval', win, white, FixationRect);
        Screen('Flip', win);   

        while (GetSecs - trial_t0) < (sample_dur + ISI_dur + sample_dur + delay_dur + disc_dur + TJ_dur)
        % while (GetSecs-t0)<trial*(sample_dur+ISI_dur+sample_dur+delay_dur+disc_dur+TJ_dur)+(trial-1)*(recall_dur+FB_dur+ITI_dur)      
            % ----- check for a keypress -----
            % KbCheck returns a boolean keyIsDown plus a vector keyCode indicating which key(s) are pressed.
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
            if keyIsDown == 1
                repskey = find(keyCode,1);
                %    % If the participant happens to hold two keys at once, repskey becomes a vector ([leftKey rightKey], for example). Testing it with if repskey == leftKey will then error or give an unexpected result.
                if repskey ==escapeKey
                    Priority(0);
                    errorMsg = {'ESC pressed, experiment aborted'}
                    ShowCursor;
                    sca;
                else
                    disc_key(count) = repskey;
                    %CCW left arrow keycode80 CW right arrow keycode79
                    if count ==1 && repskey ==ltiltKey && gotResponse == 0
                        disc_resp(trial) = 1;%'left';
                        tj_rt(trial) = GetSecs-t1;
                        gotResponse = 1;
                        if disc_resp(trial)  == disc_side(trial)
                            disc_acc(trial) = 1;
                        else
                            disc_acc(trial) = 0;
                        end
                    elseif count ==1 && repskey ==rtiltKey && gotResponse == 0
                        disc_resp(trial) = 2;%'right';
                        tj_rt(trial) = GetSecs-t1;
                        gotResponse = 1;
                        if disc_resp(trial)  == disc_side(trial)
                            disc_acc(trial) = 1;
                        else
                            disc_acc(trial) = 0;
                        end                    
                    end
                    count = count+1;
                end
                
            end
            if count==1
                Screen('FillOval', win, white, FixationRect);
            else
                Screen('FillOval', win, white, FixationRect);
                if ispractice ==1
                    if disc_acc(trial) == 1
                        Screen('FrameOval', win, green, FixationRect,3,1);
                        
                    elseif disc_acc(trial) == 0
                        Screen('FrameOval', win, red, FixationRect,3,1);
                    end
                else                    
                Screen('FrameOval', win, black, FixationRect,3,1);
                end 
            end
            Screen('Flip', win);
        end
        if gotResponse == 1
            gotRespTJ(trial) = 1;
            if isnan(missed_trials(trial,1))
                missed_trials(trial,1) = 0;  % TJ not missed
            end
        elseif gotResponse==0
                gotRespTJ(trial) = 0;
                tj_rt(trial) = nan;
                disc_resp(trial) = nan;
                disc_acc(trial) = nan;
                if repeatmissedtrials == 1
                    if ~ismember(trial, noresptrials)  % in case subject misses both TJ and Recall within a trial               
                        noresptrials(end+1)      = trial;% noresptrials = [noresptrials trial];
                        noresptrialssave(end+1)  = trial;% noresptrialssave = [noresptrialssave trial];
                    end
                end
                missed_trials(trial,1) = 1;  % TJ missed
        end        

%--------Recall----------
        test_theta = rand()*pi; % pick a new random probe angle for the upcoming recall phase
        probe_ori = test_theta;
        probe_ori(trial) = test_theta/pi*180;
        firstAdjustTime = nan;
        repskey = [];
        gotRecallResponse = 0;
        trial_error = nan;
        Screen('TextSize', win, 40);  
        % centre the cue text horizontally above the wheel
        tb     = Screen('TextBounds', win, cueLabel);
        textX  = xwheel_center - tb(3)/2;            % centred on the wheel
        textY  = ywheel - radius_in_pixel - 50;      % 50 px above the wheel bar        
        
        t2 = GetSecs;

        while (GetSecs - trial_t0) < (sample_dur + ISI_dur + sample_dur + delay_dur + disc_dur + TJ_dur + recall_dur)
        % while (GetSecs-t0)<trial*(sample_dur+ISI_dur+sample_dur+delay_dur+disc_dur+TJ_dur+recall_dur)+(trial-1)*(FB_dur+ITI_dur)      
            Screen('FrameOval', win, black, wheelCenterRect, 2, 1);
            draw_theta = pi/2 - test_theta;
            Screen('DrawLine', win, black, -cos(draw_theta)*radius_in_pixel*0.75+xwheel_center, sin(draw_theta)*radius_in_pixel*0.75+ywheel, ...
                cos(draw_theta)*radius_in_pixel*0.75+xwheel_center, -sin(draw_theta)*radius_in_pixel*0.75+ywheel, 2);
            Screen('DrawText', win, cueLabel, textX, textY, black);            
            Screen('Flip', win);
            %get keyresponse
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
            if keyIsDown == 1
                repskey = find(keyCode,1);
                % repskey=find(keyCode ==1); % If the participant happens to hold two keys at once, repskey becomes a vector ([leftKey rightKey], for example). Testing it with if repskey == leftKey will then error or give an unexpected result.
                if repskey ==escapeKey
                    Priority(0);
                    errorMsg = {'ESC pressed, experiment aborted'}
                    ShowCursor;
                    sca;
                elseif repskey == leftKey || repskey == rightKey
                    if isnan(firstAdjustTime)
                        firstAdjustTime = GetSecs;  % mark first dial movement
                    end            
                    %CCW left arrow keycode80 CW right arrow keycode79
                    if repskey == leftKey
                        test_theta = test_theta - 1.3/180*pi;
                    elseif repskey == rightKey
                        test_theta = test_theta + 1.3/180*pi;
                    end
                    if test_theta > pi
                        test_theta = test_theta - pi;
                    elseif test_theta < 0
                        test_theta = test_theta + pi;
                    end                    
                elseif repskey == confirmKey
                    % RT measured from start of epoch ~ confirm key press
                    recall_rt(trial) = GetSecs - t2;
                    % RT2 measured from start-of-adjustment to confirm key
                    % press
                    if ~isnan(firstAdjustTime)
                        recall_rt2(trial) = GetSecs - firstAdjustTime;
                    else
                        recall_rt2(trial) = nan;  % never adjusted the dial before confirming
                    end
                    gotRecallResponse = 1;
                    break;
                end
            end
            trial_report = test_theta;
        end

        if gotRecallResponse == 1
            gotRespRecall(trial) = 1;
            if isnan(missed_trials(trial,2))
                missed_trials(trial,2) = 0;  % recall not missed
            end
    %         trial_report=trial_report-pi/2;
    %         if trial_report<0
    %             trial_report = trial_report+pi;
    %         end
            trial_error = trial_report - memory_orientation/180*pi;
            if trial_error > pi/2
               trial_error = trial_error-pi;
            elseif trial_error < -pi/2
                trial_error = trial_error+pi;
            end        
      
            wm_error(trial) = trial_error/pi*180;
            wm_report(trial) = trial_report/pi*180;      

        elseif gotRecallResponse==0
                gotRespRecall(trial) = 0;
                recall_rt(trial) = nan;
                wm_report(trial) = nan;
                wm_error(trial) = nan;
                if repeatmissedtrials == 1
                    if ~ismember(trial, noresptrials) % in case subject misses both TJ and Recall within a trial
                        noresptrials(end+1)      = trial;%noresptrials = [noresptrials trial];
                        noresptrialssave(end+1)  = trial;%noresptrialssave = [noresptrialssave trial];
                    end
                end
                missed_trials(trial,2) = 1;  % recall missed
        end                      
        % %calculate bias toward non-matching stimulus
        % if conditions(trial,1)==1
        %     wm_bias(trial) = trial_error/pi*180;
        % else
        %      if conditions(trial,4)<0 %determine if discrimination ori is CW or CCW to the memory ori
        %         wm_bias(trial) = trial_error*-1/pi*180;
        %      else
        %         wm_bias(trial) = trial_error/pi*180;
        %      end
        % end
%--------------Feedback---------------        
        dotColor = [0 0 0];
        if gotRecallResponse == 1 &&  ~isnan(trial_error)
            absErrDeg = abs(trial_error) / pi * 180;      % convert to degrees
           
            if      absErrDeg < 20
                dotColor = [0 255 0];   % green
            elseif  absErrDeg <= 28          % 20 ≤ error ≤ 28
                dotColor = [255 255 0];   % yellow
            else                            % ≥ 28 deg
                dotColor = [255 0 0];   % red
            end
        end

        while (GetSecs - trial_t0) < (sample_dur + ISI_dur + sample_dur + delay_dur + disc_dur + TJ_dur + recall_dur +FB_dur)   
        % while (GetSecs-t0)<trial*(sample_dur+ISI_dur+sample_dur+delay_dur+disc_dur+TJ_dur+recall_dur+FB_dur)+(trial-1)*(ITI_dur)   
            if gotRecallResponse == 1
                if ispractice ==1
                    Screen('DrawTexture', win, gabortex, gabortex_rect, center_posrect, ...
                        -memory_orientation, [], [], [], [], kPsychDontDoRotation, propertiesMat');
                    Screen('FrameOval', win, dotColor, wheelCenterRect, 2, 1);
                    draw_theta = pi/2 - trial_report;  % map 0 rad to 12 o'clock, clockwise positive
                    x0 = -cos(draw_theta) * radius_in_pixel * 0.75 + xwheel_center;
                    y0 =  sin(draw_theta) * radius_in_pixel * 0.75 + ywheel;
                    x1 =  cos(draw_theta) * radius_in_pixel * 0.75 + xwheel_center;
                    y1 = -sin(draw_theta) * radius_in_pixel * 0.75 + ywheel;
                    Screen('DrawLine', win, dotColor, x0, y0, x1, y1, 4);                    
                else
                    Screen('FillOval', win, dotColor, FixationRect);  % color-coded dot
                end 
            elseif gotRecallResponse ==0
                Screen('TextSize',win,30);
                text1='No response was recorded.';
                text1=WrapString(text1,100);
                DrawFormattedText(win, text1, 'center', 'center', black);  
            end
            Screen('Flip', win);
        end
        % --- Press space to proceed to the next trial ---
        if ispractice ==1
    
            Screen('TextSize', win, 30);
            DrawFormattedText(win, 'Press space when ready.', 'center', 'center', black);
            Screen('Flip', win);
            
            while true
                [~,~,keyCode] = KbCheck(-1);
                if keyCode(moveonKey), break; end       % wait for SPACE
                if keyCode(escapeKey), ShowCursor; sca; return; end  % allow ESC to quit
            end        
        end
%--------------ITI---------------     
        if ispractice == 1
            ITI_dur = 3;
        end
        while (GetSecs - trial_t0) < (sample_dur + ISI_dur + sample_dur + delay_dur + disc_dur + TJ_dur + recall_dur +FB_dur + ITI_dur)   
        % while (GetSecs-t0)<trial*(sample_dur+ISI_dur+sample_dur+delay_dur+disc_dur+TJ_dur+recall_dur+FB_dur+ITI_dur)
            Screen(win, 'FillRect', gray);
            Screen('FillOval', win, white, FixationRect);
            Screen('Flip',win);
        end  
%--------------Record Trial Condition---------------        
        %1-oricong 2-cue type 3-discriminandum_contrast 4-mem_disc_distance
    record.ori_cong(trial)         = run_ori_cong(trial); % column 1
    record.run_cue_type(trial)     = run_cue_type(trial); % column 2
    record.disc_contrast(trial)    = run_disc_contrast(trial); % column 3
    record.distance(trial)         = run_distance(trial); % column 4
    record.ori1_label(trial)       = ori1_label; % column 5
    record.ori2_label(trial)       = ori2_label; % column 6        
    record.orientation1(trial)     = orientation1; % column 7
    record.orientation2(trial)     = orientation2; % column 8
            
    record.target_ori(trial)       = memory_orientation; 
    record.disc_orientation(trial) = disc_orientation;       
    record.recall_rt(trial)        = recall_rt(trial); 
    record.recall_rt2(trial)       = recall_rt2(trial); 
    record.tj_rt(trial)            = tj_rt(trial); 
    
    record.disc_resp(trial)        = disc_resp(trial); 
    record.disc_acc(trial)         = disc_acc(trial);  
    record.probe_ori(trial)        = probe_ori(trial); 
    record.wm_error(trial)         = wm_error(trial); 
    record.wm_report(trial)        = wm_report(trial); 
    record.gotRespTJ(trial)        = gotRespTJ(trial); 
    record.gotRespRecall(trial)    = gotRespRecall(trial); 
    % record.wm_bias(trial)          = wm_bias(trial);
    block_duration = GetSecs-t0;


        if repeatmissedtrials ==1 
            if finishflag==0
                if trial<trialnumber % we're still working on the trials
                    trial = trial+1;
                elseif trial==trialnumber && isempty(noresptrials) % we've gone through all of the trials, have no repeats, and can finish up
                    trial = trial+1;
                elseif trial==trialnumber && ~isempty(noresptrials)  % we've gone through all of the trials, but need to repeat some earlier ones
                    finishflag = 1;
                    nrts = unique(noresptrials);
                    trial = nrts(1);
                    nrinds = find(noresptrials==nrts(1));
                    noresptrials(nrinds) = [];
                end
            elseif finishflag == 1
                if trial==trialnumber && isempty(noresptrials)  % we've gone through all of the trials, and repeats, and can finish up
                    trial = trial+1;
                elseif trial<trialnumber && isempty(noresptrials) % we've gone through all of the trials and all of the repeats, and can finish up
                    trial = trialnumber+1;
                elseif trial<trialnumber && ~isempty(noresptrials)
                    nrts = unique(noresptrials);
                    trial = nrts(1);
                    nrinds = find(noresptrials==nrts(1));
                    noresptrials(nrinds) = [];
                    % trial = noresptrials(1);     % first missed trial
                    % noresptrials(1) = [];        % remove it                
                elseif trial==trialnumber && ~isempty(noresptrials)
                    nrts = unique(noresptrials);
                    trial = nrts(1);
                    nrinds = find(noresptrials==nrts(1));
                    noresptrials(nrinds) = [];
                    % trial = noresptrials(1);     % first missed trial
                    % noresptrials(1) = [];        % remove it                               
                end
            end
        else
            trial = trial + 1;
        end
    
    % runningtrialssave = [runningtrialssave runningtrialssave(end)+1]; % record how many trials were run for each block including repeated trials where there were no reponses

    end
        
%--------------end of block---------------    
%     block_data = [conditions(1:trialnumber,:), disc_acc',disc_rt',disc_resp', wm_error', ...
%         wm_bias', target_ori',ntgt_ori', wm_report'/pi*180, ...
%         test_ori',disc_orientation,disc_ori1,disc_ori2];
    % force every field in record to be a column vector
    fn = fieldnames(record);
    for i = 1:numel(fn)
        record.(fn{i}) = record.(fn{i})(:);
    end
    try
        save(filename,'record','runningtrialssave','block_duration','missed_trials');
    catch
        'error'
    end
    mean_abs_error = mean(abs(wm_error), 'omitnan'); % wm_error = per‑trial signed error in degrees exculding no response trials
    msg = sprintf(...
      'End of block.\nMean abs error: %.2f°\n\nPlease call the experimenter.', ...
      mean_abs_error);    

    Screen(win, 'FillRect', gray);                
    Screen('TextSize', win, 30);                       % choose text size
    DrawFormattedText(win, msg, 'center', 'center', black);        
    Screen('Flip', win);                             
    KbWait;                                            % wait for a keypress    
    % Clear screen
    ShowCursor;
    sca;
catch
    sca;
     %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    Priority(0); %Usage: oldPriority=Priority([newPriority])
    psychrethrow(psychlasterror); %same as 'rethrow': reissure error; rethrow(err),usually used in 'try-catch' statement
    ShowCursor;
end
end
