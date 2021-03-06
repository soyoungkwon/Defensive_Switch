% This script will generate modified packman stimuli with electric shock
% Many commands are based on DotDemo(showSprites, waitframes)

clear all; clc
AssertOpenGL;
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1);

screens = Screen('Screens');
screenNumber = max(screens);
[w, rect] = Screen('OpenWindow', screenNumber, 0);
[center(1), center(2)] = RectCenter(rect);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% instruction page
DrawFormattedText(w, 'Attack Mode\n\n You are a WHITE BALL \nMove Up/Down/Left/Right\n to catch the RED BALL', 400, 'center', [255 255 255]);
vbl=Screen('Flip', w);
% Wait for a keystroke on the keyboard:
KbStrokeWait;
DrawFormattedText(w, 'Escape Mode\n\n You are a WHITE BALL  \nMove Up/Down/Left/Right\n to escape the GREEN BALL',  1200,'center', [255 255 255]);
vbl=Screen('Flip', w);
KbStrokeWait;

% key code define
escapeKey = KbName('Escape');
rightKey  = KbName('RightArrow');
leftKey   = KbName('LeftArrow');
upKey     = KbName('UpArrow');
downKey   = KbName('DownArrow');
buttonR   = [];

% open the screen
waitframes = 1;     % Show new dot-images at each waitframes'th monitor refresh.

%====== TIME ==========%
% calculate refresh rate & smooth dots

fps = Screen('FrameRate',w);      % frames per second
ifi = Screen('GetFlipInterval', w);
fps = 1/ifi;
rfr = round(fps);
% moveRate = rfr*0.2;%*2;%0.4;

% mode timing
switchT = repmat([5 6 7.5 6.5 8 9 7 8 7 9 7.5]+20, [1 4]);
switchTime = cumsum(switchT);
switchframe = [round(fps*switchTime)];
minRT = 0.2;
%     modeT = [1 0; 2 5; 1 10; 2 15];
nframes = sum(switchframe);%rfr*50; % number of animation frames in loop

%====== SPACE =========%
% cursor & initial flip & check screen size
HideCursor; % Hide the mouse cursor
Priority(MaxPriority(w));
% vbl=Screen('Flip', w);
[wWidth, wHeight]=Screen('WindowSize', w);
% Clamp point sizes to range supported by graphics hardware:
[minsmooth,maxsmooth] = Screen('DrawDots', w);
s1 = 60;%was 80min(max(s, minsmooth), maxsmooth);% dot size
s2 = 50;
vbl=Screen('Flip', w);

% lines & box positions (in pixels)
% :xPositions (n-1, n-1), yPositions (120,240,360,480,600,720,840,960)
numDemoLines=8; % how many lines in the matrix
halfLines=numDemoLines/2;
lineWidths{1}=6; lineWidths{2}=10;
xPositions=round(linspace(0, wWidth, numDemoLines+2));
xPositions=xPositions(2:end-1);
yPositions=round(linspace(0, wHeight, numDemoLines+2));
yPositions=yPositions(2:end-1);

% one block movement define
% e.g., xGap: 120
xGap = round(mean(diff(xPositions)));
yGap = round(mean(diff(yPositions)));
xMove = xGap; yMove = yGap;

% DOT positions 
moveX = []; moveY = []; resp = []; distArray = []; tArray = [];
moveXYarray = []; %dist4ok = [];
dotMatrix = []; tcaught = []; caughtA = zeros(1,2); allFrames = [];
IJ=[1:(numDemoLines-1)*(numDemoLines-1)]; % all index
DotAll = IJ;
% initial random position of the prey & predator
[xDot,yDot] = ind2sub([numDemoLines-1,numDemoLines-1],DotAll);%DotNoWall);%
xR = xDot(1); yR = yDot(1);
stepGap = (numDemoLines)/2;

% initial position
% xy14 = [xGap*xR-xGap*stepGap yGap*yR-yGap*stepGap];
xy14 = [xDot(1)-stepGap yDot(1)-stepGap];% xR, yR
xy24 = [0 0];% dot positions in Cartesian coordinates (pixels from center)

%---- color code for dots/boxes/lines
myColors{1}=[255 255 255]; % white
myColors{2}=[255 0 0]; % r
myColors{3}=[0 255 0]; % g
myColors{4}=[0 0 255]; % b
myColors{5}=[50 50 50];
myColors{6}=[0 0 0];

polyCenterX{1}=-40; polyCenterY{1}=0;
polyCenterX{2}=0; polyCenterY{2}=-40;
polyCenterX{3}=40; polyCenterY{3}=0;
polyCenterX{4}=0; polyCenterY{4}=40;
numPoints=3;
polyRadius=25;
% angles=(rand(1,numPoints) * 2 * pi);
angles{1}=([0.2 -0.2 -1] * 2) * pi;%)+pi/2;
angles{2}=([0.2 -0.2 -1] * 2+0.5) * pi;%)+pi/2;
angles{3}=([0.2 -0.2 -1] * 2+1)*pi;
angles{4}=([0.2 -0.2 -1] * 2+1.5)*pi;
% angles=([0.2 -0.2 -1] * 2 * pi);
% attack or escape?
AEmode = 1; % attack escape mode (Attack:1, Escape:-1)
m = 1; % ?
p = 0; % press
speedE = 0.25;
speedA = 0.1;
% --------------
% animation loop
% --------------
tic
sw=1; % switch steps
cc=1;
for i = 1:nframes
    telapsed=toc;
    % xGap, yGap multiplied
    xy14M = xy14.*[xGap yGap];
    xy24M = xy24.*[xGap yGap];
    xymatrix1 = transpose(xy14M);% comment out soyoung
    xymatrix2 = transpose(xy24M);
    if (i>0)
        % decide attack or escape mode
        if mod(i, switchframe(sw))==0 %|| sum(abs(xy14-xy24))==0 % switch mode, either time has reached or prey=predator 
            AEmode = -AEmode;
            sw = sw+1;
            xy14 = [xDot(sw)-stepGap yDot(sw)-stepGap];
            xy24 = [xDot(sw) yDot(sw)];
        end

        % lines inside the boarder
        for ll=2:numDemoLines-1
            Screen('DrawLine', w, myColors{5}, xPositions(ll), yPositions(1),  xPositions(ll), yPositions(numDemoLines), lineWidths{2});
            Screen('DrawLine', w, myColors{5}, xPositions(1), yPositions(ll),  xPositions(numDemoLines), yPositions(ll), lineWidths{2});
        end 
        Screen('LineStipple', w, 0);

        %==== boarder ====%
        % zero coordinate (0,0) --> corner
        % e.g., 0 ~ 8
        % vertical lines
        Screen('DrawLine', w, myColors{1}, xPositions(1), yPositions(1),  xPositions(1), yPositions(numDemoLines), lineWidths{2});
        Screen('DrawLine', w, myColors{1}, xPositions(numDemoLines), yPositions(1),  xPositions(numDemoLines), yPositions(numDemoLines), lineWidths{2});
        % horizontal lines
        Screen('DrawLine', w, myColors{1}, xPositions(1), yPositions(1),  xPositions(numDemoLines), yPositions(1), lineWidths{2});
        Screen('DrawLine', w, myColors{1}, xPositions(1), yPositions(numDemoLines),  xPositions(numDemoLines), yPositions(numDemoLines), lineWidths{2});
        
        %==== bar (catch/caught)
        Screen('DrawLine', w, myColors{3}*0.8, 700, yPositions(1)-40,  700+20*m, yPositions(1)-40, lineWidths{2});
        if AEmode ==1,% 1 ==> escape
            if mod(i,30)<28,
                Screen('DrawText', w, 'Escape',  xPositions(round(numDemoLines/2)), 30, [0, 155, 155, 255]);
            end
        elseif AEmode == -1, % -1 ==> attack
            if mod(i,30)<28,
                Screen('DrawText', w, 'Attack',  xPositions(round(numDemoLines/2)), 30, [255, 150, 150, 155]);
            end
        end
        

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%===== Predator moving (automatically) (rarely) ========%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if AEmode == 1, % escape
        moveRate = rfr*speedE;
    elseif AEmode == -1, % attack
        moveRate = rfr*speedA;
    end
    if mod(i,moveRate)==1,% move only occasionally
        % distance current
        dist2 = abs(xy14-xy24);
        dist = dist2(1)^2+dist2(2)^2;
        
        % move all 4 possibility(1 or -1)
        dist4ok = [];
        move4opts = [2*(randperm(2)-1.5) zeros(1,2)];%(randperm(3)-2);
        
        % check all 4 possible move and check its distance & absolute position every 1second!
        rand4 = randperm(4);
        moveLCRX = move4opts(rand4);
        d = 1;
        
        % check location (boarder)
        % this part, make 4 possible movement (NOT added yet), dxdy14C
        % --> choose the movement that
        % makes distance shorter --> choose 1st options among leftovers
        for c=1:4,
            xPredator = moveLCRX(c); % xPredator
            if xPredator == 0, % if X no move, move Y
                moveLCRY = 2*(randperm(2)-1.5); yPredator = moveLCRY(1);
            else, % if X moved, Y should stay
                yPredator = 0;
            end

            % possible changed location (Attack & Escape)
            % 14C are attacker, 24C are escaper
            % C: all possible components
                moveC = [xPredator yPredator];
                xy14C = xy14 + moveC;
                nomoveC = [0 0];
                xy24C = xy24 + nomoveC;
                
                %=== new distance (all 4 possibilities)
                % calculate distance between attacker(14) and escaper (24)
                distt = abs(xy14C-xy24C);
                distnew4 = distt(1)^2+distt(2)^2;
                dist4(c) = distnew4 - dist;
                if dist4(c) < 0, % distance shorter than previous distance
                    dist4ok = [dist4ok; dist4(c)];
                    distSelect = dist4(c);
                    % stamp movement & location (to avoid walls/lines)
                    if d==1, % 1st possiblity from 4 options
                        xPredatorR = xPredator; yPredatorR = yPredator;
%                         dxdy14R = [xPredatorR yPredatorR]; dxdy24R = [0 0];
                    end
                    d = d + 1;
                end
        end
        
        % choose which dxdy!!!
        if AEmode == 1, % escape
            dxdy14R = [xPredatorR yPredatorR]; % choose the shortest!!!
            dxdy24R = [0 0];
        elseif AEmode == -1, % reverse the change
            dxdy24R = [xPredator yPredator];
            dxdy14R = [0 0];
        end

    else, % most of the time, no move!
        dxdy14R = [0 0]; dxdy24R = [0 0];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%===== Prey moving (by pressing button)===================%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [ keyIsDown, time, keyCode ] = KbCheck;
    whichKey = find(keyCode);
    if size(whichKey,1)>1,
    whichKey = whichKey(1);
    end
    % ======= ESCAPE mode ========%
    if AEmode == 1,
        % check any key pressed (to avoid pressed time)
        if keyCode(rightKey) || keyCode(leftKey) || keyCode(downKey) || keyCode(upKey)
            if size(resp,1) == 0 % first pressed
                p = p + 1;
                resp = [resp; AEmode p telapsed whichKey];
            elseif telapsed - resp(end,3) > minRT % avoid stamping continuous press
                p = p + 1;                
                resp = [resp; AEmode p telapsed whichKey];    
            end
            if telapsed == resp(end,3) % update movement only when it's exact moment
                % move by pressing
                if keyCode(rightKey), dxdy24R = [1 0]; 
                elseif keyCode(leftKey), dxdy24R = [-1 0];
                elseif keyCode(upKey), dxdy24R = [0 -1];
                elseif keyCode(downKey), dxdy24R = [0 1];
                end;
            end
%             end
        end
        
    % ======= ATTACK mode ========%
    elseif AEmode == -1,
        % check any key pressed (to avoid pressed time)
        if keyCode(rightKey) || keyCode(leftKey) || keyCode(downKey) || keyCode(upKey)
            if size(resp,1) == 0 % first pressed
                p = p + 1;
                resp = [resp; AEmode p telapsed whichKey];
            elseif telapsed - resp(end,3) > minRT % avoid stamping continuous press
                p = p + 1;
                resp = [resp; AEmode p telapsed whichKey];
            end
            if telapsed == resp(end,3) % update movement only when it's exact moment
                % move by pressing
                if keyCode(rightKey), dxdy14R = [1 0];%[xGap 0;];
                elseif keyCode(leftKey), dxdy14R = [-1 0];%[-xGap 0;];
                elseif keyCode(upKey), dxdy14R = [0 -1];%[0 -yGap;];
                elseif keyCode(downKey), dxdy14R = [0 1];%[0 yGap;];
                end; 
            end
        end
    end
    
    % quit
    if keyCode(escapeKey), break; end
    
    %======= bar =======%
    if sum(abs(xy14-xy24))==0
        tcaught = [tcaught; telapsed];
        if mod(i,rfr) >1 & mod(i,rfr) == 11 % only if mod(i,moveRate)==1  %
            if telapsed -tcaught(end) > minRT
                if AEmode == 1,
                    m = m - 1;
                    fprintf('caught %0.1d %0.2f\n', m, telapsed);
                    Screen('DrawText', w, '- -',  xPositions(round(numDemoLines/2))+200, 50, [0, 0, 255, 255]);
                elseif AEmode == -1,
                    m = m + 1;
                    fprintf('got it %0.1d %0.2f\n', m, telapsed);
                    Screen('DrawText', w, '+ +',  xPositions(round(numDemoLines/2))+200, 50, [0, 0, 255, 255]);
                end
            end
        end
    end
    

    
    %% ===== Draw nice dots =====%
    % zero coordinate (0,0) --> center
    % e.g., -4 ~ +4
    if xPredator == 1 & yPredator == 0, ww=1; elseif xPredator == -1 & yPredator == 0; ww = 2; 
    elseif xPredator == 0 & yPredator == 1; ww=3; elseif xPredator == 0 & yPredator == -1; ww=4; end
    xPos = round([cos(angles{ww})*polyRadius]') + center(1) + polyCenterX{ww};
    yPos = round([sin(angles{ww})*polyRadius]') + center(2) + polyCenterY{ww};
    %         polyPoints = polyPoints1 + center(2);
    if AEmode == 1,%strcmpi(mode, 'escape'),
        Screen('DrawDots', w, xymatrix1, s1, myColors{3}, center, 1);  % xymatrix1 ==> attacker
        Screen('DrawDots', w, xymatrix2, s2, myColors{1}, center, 1);  % xymatrix2 ==> me
        %             Screen('FillRect', w, xymatrix1, s2);%InsetRect(rect, 900, 400));%[255 100 0]
        if mod(i,30) < 15, % mouth of the predator
            Screen('FramePoly', w, myColors{6}, (xymatrix1)'+[xPos yPos],30);
        end
    elseif AEmode == -1, %strcmpi(mode, 'attack')
        Screen('DrawDots', w, xymatrix1, s2, myColors{1}, center, 1);  % xymatrix1 ==> me
        Screen('DrawDots', w, xymatrix2, s1, myColors{2}, center, 1);  % change 1 to 0 or 4 to draw square dots
        if mod(i,30) < 15, % mouth of the predator
            Screen('FramePoly', w, myColors{6}, (xymatrix1)'+[xPos yPos],30);%30);
        end
    end
    Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
        
    % HERE actual location update!!!!
    xy14 = xy14 + dxdy14R; % move dots
    xy24 = xy24 + dxdy24R;
    
    % =====location restrain (not outside of the boarder)
    if xy14(1) > halfLines-1, xy14(1) = halfLines-1; 
    elseif xy14(1) < -halfLines+1, xy14(1) = -halfLines+1; end
    if xy14(2) > halfLines-1, xy14(2) = halfLines-1; 
    elseif xy14(2) < -halfLines+1, xy14(2) = -halfLines+1; end
    if xy24(1) > halfLines-1, xy24(1) = halfLines-1; 
    elseif xy24(1) < -halfLines+1, xy24(1) = -halfLines+1; end
    if xy24(2) > halfLines-1, xy24(2) = halfLines-1; 
    elseif xy24(2) < -halfLines+1, xy24(2) = -halfLines+1; end
    


    % location stamp
       if mod(i,moveRate)==2,
         moveXYarray = [moveXYarray; xPredatorR yPredatorR];
         dotMatrix = [dotMatrix; AEmode telapsed xy14 xy24 sum(abs(xy14-xy24))==0]
         tArray = [tArray; telapsed AEmode];
         if sum(abs(xy14-xy24))==0 %&& (telapsed-caughtA(cc,2) > minRT)% prey=predator
             cc=cc+1;
             caughtA = [caughtA; AEmode telapsed]
             xy14 = [xDot(cc)-stepGap yDot(cc)-stepGap];
             xy24 = [xDot(cc) yDot(cc)];
% %              cc=mod(cc,length(IJ));
         end
       end
       allFrames = [allFrames; AEmode i telapsed xy14 xy24 sum(abs(xy14-xy24))==0];
       vbl=Screen('Flip', w, vbl + (waitframes-0.5)*ifi);
    
    end
end;
Priority(0);
ShowCursor;
sca;
if saveOn
   save([sprintf('testEscape%0.2dAttack%0.2dtrial27', speedE*10, speedA*10)], 'dotMatrix', 'allFrames');
end
