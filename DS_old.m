% This script will generate modified packman stimuli with electric shock
% Many commands are based on DotDemo(showSprites, waitframes)

clear all; clc
AssertOpenGL;

% key code define
escapeKey = KbName('Escape');
rightKey  = KbName('RightArrow');
leftKey   = KbName('LeftArrow');
upKey     = KbName('UpArrow');
downKey   = KbName('DownArrow');
buttonR   = [];

tic
% try
Screen('Preference', 'SkipSyncTests', 1);
% open the screen
waitframes = 1;     % Show new dot-images at each waitframes'th monitor refresh.
screens = Screen('Screens');
screenNumber = 0;%max(screens);
% [w, rect] = Screen('OpenWindow', screenNumber, 0);
[w, rect] = Screen('OpenWindow', screenNumber, 1);

%====== TIME ==========%
% calculate refresh rate & smooth dots
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[center(1), center(2)] = RectCenter(rect);
fps = Screen('FrameRate',w);      % frames per second
ifi = Screen('GetFlipInterval', w);
fps = 1/ifi;
rfr = round(fps);
moveRate = rfr*0.2;%0.4;
nframes = rfr*50; % number of animation frames in loop

% mode timing
switchT = [5 4 5.5 4.5 4 6 7 4 5 6 5]*2;
switchframe = [round(fps*switchT)];
minRT = 0.2;
%     modeT = [1 0; 2 5; 1 10; 2 15];

%====== SPACE =========%
% cursor & initial flip & check screen size
HideCursor; % Hide the mouse cursor
Priority(MaxPriority(w));
vbl=Screen('Flip', w);
[wWidth, wHeight]=Screen('WindowSize', w);
% Clamp point sizes to range supported by graphics hardware:
[minsmooth,maxsmooth] = Screen('DrawDots', w)
s1 = 60;%was 80min(max(s, minsmooth), maxsmooth);% dot size
s2 = 50;

% lines & box positions (in pixels)
% :xPositions (n-1, n-1), yPositions (120,240,360,480,600,720,840,960)
numDemoLines=8; % how many lines in the matrix
lineWidths{1}=6;
lineWidths{2}=10;
xPositions=round(linspace(0, wWidth, numDemoLines+2));
xPositions=xPositions(2:end-1);
yPositions=round(linspace(0, wHeight, numDemoLines+2));
yPositions=yPositions(2:end-1);

% one block movement define
% e.g., xGap: 120
xGap = round(mean(diff(xPositions)));
yGap = round(mean(diff(yPositions)));
xMove = xGap; yMove = yGap;

% WALL positions: IJ(1-n2), could be deleted
i_caught = []; xymatrixall = [];
moveX = []; moveY = []; walls = []; resp = []; distArray = [];
moveXYarray = []; %dist4ok = [];
wI = randperm(numDemoLines-1);
wJ = randperm(numDemoLines-1);
ind=sub2ind([numDemoLines-1,numDemoLines-1],wI,wJ); dotMatrix = [];
IJ=[1:(numDemoLines-1)*(numDemoLines-1)]; % all index

% dot positions (not wall position)
DotNoWall = IJ(~ismember(IJ,ind));
% initial random position of the prey & predator
[xDot,yDot] = ind2sub([numDemoLines-1,numDemoLines-1],DotNoWall);%randperm(numDemoLines-1);%-numDemoLines/2;
xR = xDot(1); yR = yDot(1);
stepGap = (numDemoLines)/2;

% initial position
% xy14 = [xGap*xR-xGap*stepGap yGap*yR-yGap*stepGap];
xy14 = [xR-stepGap yR-stepGap];
xy24 = [0 0];% dot positions in Cartesian coordinates (pixels from center)

%---- color code for dots/boxes/lines
myColors{1}=[255 255 255]; % white
myColors{2}=[255 0 0]; % r
myColors{3}=[0 255 0]; % g
myColors{4}=[0 0 255]; % b
myColors{5}=[50 50 50];
%     white = WhiteIndex(w);
%     colvect=[255 0 0];%white;


% attack or escape?
AEmode = 1; % attack escape mode (Attack:1, Escape:-1)
m = 5; % ?
% --------------
% animation loop
% --------------
sw=1; % switch steps
for i = 1:nframes
    telapsed=toc;
    if (i>1)
        % attack & escape mode
        if mod(i,switchframe(sw))==1
           AEmode = -AEmode;
           sw=sw+1;
           xy14 = [0 0]; xy24 = [3 3];
%            dxdy14 = [3 3];
        end
        Screen('LineStipple',w, 0);
        
        % more lines
        for ll=2:numDemoLines-1
            Screen('DrawLine', w, myColors{5}, xPositions(ll), yPositions(1),  xPositions(ll), yPositions(numDemoLines), lineWidths{2});
            Screen('DrawLine', w, myColors{5}, xPositions(1), yPositions(ll),  xPositions(numDemoLines), yPositions(ll), lineWidths{2});
        end 
        
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
        if AEmode ==1,
            if mod(i,30)<28,
                Screen('DrawText', w, 'Escape',  xPositions(round(numDemoLines/2)), 30, [0, 155, 155, 255]);
            end
        elseif AEmode == -1,
            if mod(i,30)<28,
                Screen('DrawText', w, 'Attack',  xPositions(round(numDemoLines/2)), 30, [255, 150, 150, 155]);
            end
        end
%         % === walls =====%
%         for jj=1:numDemoLines-1
%             Screen('FillRect', w, myColors{1}, [xPositions(wI(jj)), yPositions(wJ(jj)), xPositions(wI(jj)+1), yPositions(wJ(jj)+1)]);  % draw fixation dot (flip erases it)
%         end
        
        %% ===== Draw nice dots =====%
        % zero coordinate (0,0) --> center
        % e.g., -4 ~ +4
        if AEmode == 1,%strcmpi(mode, 'escape'),
            Screen('DrawDots', w, xymatrix1, s1, myColors{3}, center, 1);  % change 1 to 0 or 4 to draw square dots
            Screen('DrawDots', w, xymatrix2, s2, myColors{1}, center, 1);  % change 1 to 0 or 4 to draw square dots
        elseif AEmode == -1, %strcmpi(mode, 'attack')
            Screen('DrawDots', w, xymatrix1, s2, myColors{2}, center, 1);  % change 1 to 0 or 4 to draw square dots
            Screen('DrawDots', w, xymatrix2, s1, myColors{1}, center, 1);  % change 1 to 0 or 4 to draw square dots
            
        end
        Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
        
    end;
    
    %% ===== Predator moving (automatically) (rarely)
    if mod(i,moveRate)==1,
        % distance current
        dist2 = abs(xy14-xy24);
        dist = dist2(1)^2+dist2(2)^2;
        
        % move all 4 possibility(1 or -1)
        dist4ok = [];
        move4opts = [2*(randperm(2)-1.5) zeros(1,2)];%(randperm(3)-2);
        
        % check all 4 possible move and check its distance & absolute
        % position
        % every 1second!
        rand4 = randperm(4);
        moveLCRX = move4opts(rand4);
        d = 1;
        % check location (boarder)
        for c=1:4,
            xPredator = moveLCRX(c); % xPredator
            
            if xPredator == 0, % if X no move, move Y
                moveLCRY = 2*(randperm(2)-1.5); yPredator = moveLCRY(1);
            else, % if X moved, Y should stay
                yPredator = 0;
            end
            
            
            % possible changed location (Attack & Escape)
            % restrain location
            % 14 are attacker, 24 are escaper
            if AEmode == 1,
                dxdy14C = [xPredator yPredator;];% 0 0;]; %moveX = [moveX; xPredator yPredator];
                % check location
%                 xy14C = xy14 + dxdy14C;
                %=== restrain movement (boarder not possible to move)
                % after move (NOT WORKING)
                if (xy14(1)+xPredator) > numDemoLines/2, xPredator = -1;
                elseif (xy14(2)+xPredator) < -numDemoLines/2, xPredator = 1; end
     
%                 if abs(xy14C(1))<numDemoLines/2 & abs(xy14C(2))<numDemoLines/2
                    dxdy14 = [xPredator yPredator];
                    xy14C = xy14 + dxdy14;
                    dxdy24 = [0 0];
%                     xy14 = xy14 + dxdy14;
                    xy24C = xy24 + dxdy24;
%                 end
            elseif AEmode == -1,
                dxdy24C = [xPredator yPredator];
                % check location
%                 xy24C = xy24 + dxdy24C;
                if (xy24(1)+xPredator) > numDemoLines/2, xPredator = -1; 
                elseif (xy24(2)+xPredator) < -numDemoLines/2, xPredator = 1; end

%                 if abs(xy24C(1))<numDemoLines/2 & abs(xy24C(2))<numDemoLines/2
                    dxdy24 = [xPredator yPredator];
                    xy24C = xy24 + dxdy24;
                    dxdy14 = [0 0];
                    xy14C = xy14 + dxdy14;
%                     xy24 = xy24 + dxdy24;
%                 end
            end
            
            %=== new distance (all 4 possibilities)
            % calculate distance between attacker(14) and escaper (24)
            distt = abs(xy14C-xy24C);
            distnew4 = distt(1)^2+distt(2)^2;
            dist4(c) = distnew4 - dist;
            
            % restrain position (not in the wall/box)
            
            
            % restrain distance (only for shorter) & choose last one
            if dist4(c) < 0, % distance shorter than previous distance
                dist4ok = [dist4ok; dist4(c)];
                distSelect = dist4(c);
                % stamp movement & location (to avoid walls/lines)
                if d==1,
                    xPredatorR = xPredator; yPredatorR = yPredator;
                    if AEmode == 1, dxdy14R = [xPredatorR yPredatorR]; dxdy24R = [0 0];
                    elseif AEmode == -1, dxdy24R = [xPredatorR yPredatorR]; dxdy14R = [0 0]; end
                    moveXYarray = [moveXYarray; xPredatorR yPredatorR];
                    dotMatrix = [dotMatrix; xy14 xy24];%
                end
                d = d+1;
            end
            

        end
        
   
    else,
        %         dxdy = [0 0; 0 0]; % no move usually
        dxdy14R = [0 0]; dxdy24R = [0 0];
    end
    
    %% ===== Prey moving (by pressing button)
    [ keyIsDown, time, keyCode ] = KbCheck;
    whichKey = find(keyCode);
    
    % escaping mode
    if AEmode == 1,
        % check any key pressed (to avoid pressed time)
        if keyCode(rightKey) || keyCode(leftKey) || keyCode(downKey) || keyCode(upKey)
            if size(resp,1) == 0 % first pressed
                resp = [resp; AEmode i telapsed whichKey];
            elseif telapsed - resp(end,3) > minRT % avoid immediate pressed
                resp = [resp; AEmode i telapsed whichKey];
                % move by pressing
                if keyCode(rightKey), dxdy24R = [1 0];%[xGap 0];
                elseif keyCode(leftKey), dxdy24R = [-1 0];%[-xGap 0];
                elseif keyCode(upKey), dxdy24R = [0 -1];%[0 -yGap];
                elseif keyCode(downKey), dxdy24R = [0 1];%[0 yGap];
                end;
            end
        end
        
        % Attacking mode
    elseif AEmode == -1,
        % check any key pressed (to avoid pressed time)
        if keyCode(rightKey) || keyCode(leftKey) || keyCode(downKey) || keyCode(upKey)
            if size(resp,1) == 0 % first pressed
                resp = [resp; AEmode i telapsed whichKey];
            elseif telapsed - resp(end,3) > minRT % avoid immediate pressed
                resp = [resp; AEmode i telapsed whichKey];
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
    
    % bar
    if mod(i,rfr)>50,
        if abs(xy14-xy24)==0
            if telapsed - resp(end,3) > minRT
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
    
    % new location
    xy14 = xy14 + dxdy14R; % move dots
    xy24 = xy24 + dxdy24R;
    
    % xGap, yGap multiplied
    xy14M = xy14.*[xGap yGap];
    xy24M = xy24.*[xGap yGap];
    xymatrix1 = transpose(xy14M);% comment out soyoung
    xymatrix2 = transpose(xy24M);
    %         xymatrixall=[xymatrixall; xymatrix(1,:)];
    vbl=Screen('Flip', w, vbl + (waitframes-0.5)*ifi);
    
    
end;
Priority(0);
ShowCursor;
sca;
5
% catch
%     Priority(0);
%     ShowCursor;
%     sca;
%     6
% end