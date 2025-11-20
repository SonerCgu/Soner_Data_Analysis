function fUSI_Live_v1
%% =========================================================================
%   fUSI_Live_v1  —  Stable, polished fUSI Viewer (Documented Version)
%
%   KEY FEATURES:
%   • GUI + image + timecourse panel
%   • Raw / Normalized / % Signal Change modes
%   • Full ROI workflow (preview + permanent curves)
%   • Playback for 4D
%   • Metadata PNG auto export to qc/
%   • RAW + VISUAL NIfTI export
%   • Professional footer
%   • This is PART 1 (file loading + metadata + PSC)
%
%   Author: Soner Caner Cagun
%   Institute: Max Planck Institute for Biological Cybernetics
%   Version: v1 (2025)
% =========================================================================

clear; clc; close all;

%% ========================================================================
% (1) SELECT INPUT FILE
% ========================================================================
[fileName, filePath] = uigetfile('*.mat','Select fUSI MAT file');
if isequal(fileName,0), return; end

S = load(fullfile(filePath, fileName));
vars = fieldnames(S);

fprintf('\nLoaded file: %s\n', fileName);

%% ========================================================================
% (2) READ METADATA
% ========================================================================
if ~isfield(S,'metadata')
    error('This MAT file does not contain metadata structure.');
end

metadata = S.metadata;

dx = metadata.voxelSize(1);
dy = metadata.voxelSize(2);
dz = 1;

%% ========================================================================
% (3) DETECT MAIN DATA ARRAY
% ========================================================================
candidates = {};
for i = 1:numel(vars)
    if isnumeric(S.(vars{i})) && ndims(S.(vars{i})) >= 3
        candidates{end+1} = vars{i};
    end
end
if isempty(candidates)
    error('No numeric fUSI data found in MAT file.');
end

[~, idx] = max(cellfun(@(v) numel(S.(v)), candidates));
dataVar = candidates{idx};

I    = S.(dataVar);
sz   = size(I);
dims = ndims(I);

%% ========================================================================
% (4) SYSTEM TYPE + DIMENSIONS
% ========================================================================
if dims == 3
    systemType = 'Daxasonics (3D Time-Series)';
    Ny = sz(1); Nx = sz(2); T = sz(3); Nz = 1;
else
    systemType = 'Matrix Probe (4D Volumetric)';
    Ny = sz(1); Nx = sz(2); Nz = sz(3); T = sz(4);
end

%% ========================================================================
% (5) DETERMINE TR
% ========================================================================
if isfield(metadata,'time') && numel(metadata.time)>=2
    dt = metadata.time(2) - metadata.time(1);
    if dt > 0.28 && dt < 0.34
        TR = dt;
        fprintf('Detected TR: %.3f s\n', TR);
    else
        TR = 0.320;
        fprintf('Metadata irregular ? TR=0.320 s\n');
    end
else
    answ = inputdlg('Enter TR (ms):','TR Input',1,{'320'});
    TR = str2double(answ{1})/1000;
    fprintf('User TR: %.3f s\n', TR);
end

%% ========================================================================
% (6) PRINT DATASET SUMMARY
% ========================================================================
fprintf('\n================ Dataset Summary ================\n');
fprintf(' File:     %s\n', fileName);
fprintf(' System:   %s\n', systemType);

if dims == 3
    fprintf(' Dim:      %d × %d × 1 × %d\n', Ny, Nx, T);
else
    fprintf(' Dim:      %d × %d × %d × %d\n', Ny, Nx, Nz, T);
end

fprintf(' Voxel:    %.3f × %.3f × %.3f mm\n', dx,dy,dz);
fprintf(' TR:       %.3f s\n', TR);
fprintf(' Duration: %.2f min\n', (T*TR)/60);
fprintf('=================================================\n\n');

%% ========================================================================
% (7) INITIAL PSC CALCULATION
% ========================================================================
Nbaseline = min(T,1000);

if dims == 3
    base = mean(I(:,:,1:Nbaseline),3);
    PSC  = (I - base) ./ base * 100;
    anatomyVol = I;
else
    base = mean(I(:,:,:,1:Nbaseline),4);
    PSC  = bsxfun(@rdivide, bsxfun(@minus, I, base), base) * 100;
    anatomyVol = mean(I,4);
end

anatNorm = normalizeVol(anatomyVol);

%% ========================================================================
%% METADATA
%% ========================================================================
if ~isfield(S,'metadata')
    error('MAT file does not contain metadata.');
end

metadata = S.metadata;

dx = metadata.voxelSize(1);
dy = metadata.voxelSize(2);
dz = 1;

%% ========================================================================
% (9) EXPORT NIFTI (RAW + VISUAL)
% ========================================================================
rawNii = fullfile(filePath,'FUS_raw.nii.gz');
visNii = fullfile(filePath,'FUS_visual.nii.gz');

if ~exist(rawNii,'file') || ~exist(visNii,'file')
    fprintf('Exporting NIfTI files...\n');
    export_fUS_NIfTI_v5_2(I, filePath, dx, dy, dz, TR);
else
    fprintf('NIfTI exists ? skipping.\n');
end

%% ========================================================================
% (10) INTERNAL MODE STATE
% ========================================================================
rawMode = 1;   % 1 = Raw, 2 = Normalized, 3 = PSC

%% ========================================================================
% >>>>>>>>>>>>>>> END OF PART 1 (GUI STARTS IN PART 2) <<<<<<<<<<<<<<<<<
% ========================================================================
%% =========================================================================
%   PART 2 — GUI CREATION  
% =========================================================================

fig = figure('Name',['fUSI Viewer v12.3 — ' systemType], ...
    'Position',[100 40 1600 950],...
    'Color','k',...
    'MenuBar','none',...
    'ToolBar','none',...
    'WindowButtonDownFcn',@mouseClick,...
    'WindowScrollWheelFcn',@mouseWheelScroll,...
    'KeyPressFcn',@keyPressHandler);

drawnow;

%% =========================================================================
%   SIDEBAR PANEL
% =========================================================================
sidebar = uipanel(fig,'Units','normalized',...
    'Position',[0 0 0.20 1],...
    'BackgroundColor',[0.12 0.12 0.12]);

uicontrol(sidebar,'Style','text','String','fUSI Controls',...
    'Units','normalized','Position',[0.05 0.95 0.90 0.035],...
    'ForegroundColor','w','BackgroundColor',[0.12 0.12 0.12],...
    'FontSize',16,'FontWeight','bold');

%% Vertical compact layout engine
Y  = 0.92;      % start high
dY = 0.052;     % compressed spacing to ensure buttons don't get cut

%% ---------------- Display Mode ----------------
uicontrol(sidebar,'Style','text','String','Display Mode',...
    'Units','normalized','Position',[0.05 Y 0.90 0.035],...
    'ForegroundColor','w','BackgroundColor',[0.12 0.12 0.12],...
    'FontSize',12,'FontWeight','bold');
Y = Y - dY;

displayDropdown = uicontrol(sidebar,'Style','popupmenu',...
    'Units','normalized','Position',[0.05 Y 0.90 0.05],...
    'String',{'Raw Intensity','Normalized','% Signal Change'},...
    'Value',1,'FontSize',13,...
    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor','w');
Y = Y - dY;

%% ---------------- ROI Size ----------------
uicontrol(sidebar,'Style','text','String','ROI Size (px)',...
    'Units','normalized','Position',[0.05 Y 0.90 0.035],...
    'ForegroundColor','w','BackgroundColor',[0.12 0.12 0.12],...
    'FontSize',12);
Y = Y - dY;

roiSlider = uicontrol(sidebar,'Style','slider',...
    'Units','normalized','Position',[0.05 Y 0.90 0.03],...
    'Min',2,'Max',150,'Value',10,...
    'SliderStep',[1/150 10/150],...
    'BackgroundColor',[0.3 0.3 0.3]);
Y = Y - dY;

%% ---------------- Brightness ----------------
uicontrol(sidebar,'Style','text','String','Brightness',...
    'Units','normalized','Position',[0.05 Y 0.90 0.035],...
    'ForegroundColor','w','BackgroundColor',[0.12 0.12 0.12],...
    'FontSize',12);
Y = Y - dY;

brightness = uicontrol(sidebar,'Style','slider',...
    'Units','normalized','Position',[0.05 Y 0.90 0.03],...
    'Min',-1,'Max',1,'Value',0,...
    'BackgroundColor',[0.3 0.3 0.3]);
Y = Y - dY;

%% ---------------- Contrast ----------------
uicontrol(sidebar,'Style','text','String','Contrast',...
    'Units','normalized','Position',[0.05 Y 0.90 0.035],...
    'ForegroundColor','w','BackgroundColor',[0.12 0.12 0.12],...
    'FontSize',12);
Y = Y - dY;

contrast = uicontrol(sidebar,'Style','slider',...
    'Units','normalized','Position',[0.05 Y 0.90 0.03],...
    'Min',0.1,'Max',5,'Value',1,...
    'BackgroundColor',[0.3 0.3 0.3]);
Y = Y - dY;

%% ---------------- Gamma ----------------
uicontrol(sidebar,'Style','text','String','Gamma',...
    'Units','normalized','Position',[0.05 Y 0.90 0.035],...
    'ForegroundColor','w','BackgroundColor',[0.12 0.12 0.12],...
    'FontSize',12);
Y = Y - dY;

gammaSlider = uicontrol(sidebar,'Style','slider',...
    'Units','normalized','Position',[0.05 Y 0.90 0.03],...
    'Min',0.1,'Max',5,'Value',1,...
    'BackgroundColor',[0.3 0.3 0.3]);
Y = Y - dY;

%% ---------------- Colormap ----------------
uicontrol(sidebar,'Style','text','String','Colormap',...
    'Units','normalized','Position',[0.05 Y 0.90 0.035],...
    'ForegroundColor','w','BackgroundColor',[0.12 0.12 0.12],...
    'FontSize',12);
Y = Y - dY;

mapDropdown = uicontrol(sidebar,'Style','popupmenu',...
    'Units','normalized','Position',[0.05 Y 0.90 0.05],...
    'String',{'gray','hot'},...
    'FontSize',13,...
    'BackgroundColor',[0.3 0.3 0.3],'ForegroundColor','w');
Y = Y - dY;

%% ---------------- Histogram Equalization ----------------
histEQ = uicontrol(sidebar,'Style','checkbox',...
    'String','Histogram Equalization',...
    'Units','normalized','Position',[0.05 Y 0.90 0.035],...
    'ForegroundColor','w','BackgroundColor',[0.12 0.12 0.12],...
    'FontSize',12);
Y = Y - dY;

%% ---------------- Slice Slider (only for 4D) ----------------
if dims == 4
    uicontrol(sidebar,'Style','text','String','Slice (Z)',...
        'Units','normalized','Position',[0.05 Y 0.90 0.035],...
        'ForegroundColor','w','BackgroundColor',[0.12 0.12 0.12]);
    Y = Y - dY;

    sliceSlider = uicontrol(sidebar,'Style','slider',...
        'Units','normalized','Position',[0.05 Y 0.90 0.03],...
        'Min',1,'Max',Nz,'Value',round(Nz/2),...
        'SliderStep',[1/(Nz-1) 5/(Nz-1)],...
        'BackgroundColor',[0.3 0.3 0.3]);
    Y = Y - dY;
else
    sliceSlider = [];
end

%% ---------------- Baseline Inputs ----------------
uicontrol(sidebar,'Style','text','String','Baseline Start (s)',...
    'Units','normalized','Position',[0.05 Y 0.40 0.035],...
    'ForegroundColor','w','BackgroundColor',[0.12 0.12 0.12]);

uicontrol(sidebar,'Style','text','String','Baseline End (s)',...
    'Units','normalized','Position',[0.55 Y 0.40 0.035],...
    'ForegroundColor','w','BackgroundColor',[0.12 0.12 0.12]);

Y = Y - dY;

baseStart = uicontrol(sidebar,'Style','edit',...
    'Units','normalized','Position',[0.05 Y 0.40 0.04],...
    'String','0',...
    'BackgroundColor',[0.20 0.20 0.20],'ForegroundColor','w');

baseEnd = uicontrol(sidebar,'Style','edit',...
    'Units','normalized','Position',[0.55 Y 0.40 0.04],...
    'String',num2str(min(T*TR,1000*TR)),...
    'BackgroundColor',[0.20 0.20 0.20],'ForegroundColor','w');
Y = Y - dY - 0.02;

%% ---------------- Recalculate PSC ----------------
uicontrol(sidebar,'Style','pushbutton','String','Recalculate PSC',...
    'Units','normalized','Position',[0.05 Y 0.90 0.06],...
    'FontSize',14,'FontWeight','bold',...
    'BackgroundColor',[0.10 0.40 0.10],'ForegroundColor','w',...
    'Callback',@recalcPSCfunc);
Y = Y - 0.085;

%% ---------------- Help + Close Buttons ----------------
uicontrol(sidebar,'Style','pushbutton','String','Help / Info',...
    'Units','normalized','Position',[0.05 Y 0.40 0.06],...
    'FontSize',12,...
    'BackgroundColor',[0.20 0.20 0.60],'ForegroundColor','w',...
    'Callback',@showHelpWindow);

uicontrol(sidebar,'Style','pushbutton','String','Close Viewer',...
    'Units','normalized','Position',[0.55 Y 0.40 0.06],...
    'FontSize',12,...
    'BackgroundColor',[0.40 0 0],'ForegroundColor','w',...
    'Callback',@(src,event) close(fig));

%% =========================================================================
%   PROFESSIONAL FOOTER (bottom right)
% =========================================================================
uicontrol(fig,'Style','text',...
    'Units','normalized','Position',[0.73 0.01 0.26 0.03],...
    'String','Soner Caner Cagun — MPI-BC — fUSI Viewer v1 — 2025',...
    'ForegroundColor',[0.7 0.7 0.7],...
    'BackgroundColor','k',...
    'HorizontalAlignment','right',...
    'FontSize',11);

%% =========================================================================
%   IMAGE DISPLAY AXIS
% =========================================================================
ax1 = axes('Parent',fig,'Units','normalized',...
    'Position',[0.23 0.55 0.72 0.43],...
    'Color','k','XColor','w','YColor','w');
hold(ax1,'on');

if dims==3
    currentZ = 1;
    frame = anatNorm(:,:,1);
else
    currentZ = round(Nz/2);
    frame = anatNorm(:,:,currentZ);
end

frameH = imagesc(rot90(normalizeVol(frame),2));
axis(ax1,'image'); 
axis(ax1,'off');

%% =========================================================================
%   TIMECOURSE AXIS (slightly shifted up)
% =========================================================================
ax2 = axes('Parent',fig,'Units','normalized',...
    'Position',[0.23 0.08 0.72 0.40],...
    'Color','k','XColor','w','YColor','w');
hold(ax2,'on');

xlabel(ax2,'Time (s)','Color','w');
ylabel(ax2,'Raw Intensity','Color','w');

hLive = plot(ax2,(0:T-1)*TR, zeros(1,T),'LineWidth',2);

%% =========================================================================
%   ROI STRUCTURES
% =========================================================================
roiColors  = [1 0 0; 0 1 0; 0.3 0.3 1; 1 1 0; 1 0.5 0; 0 1 1; 1 0 1];
ROI        = struct('x1',{},'x2',{},'y1',{},'y2',{},'z',{},'color',{});
roiHandles = [];
roiPlots   = [];

mousePos = [NaN NaN];

%% --- LIVE ROI PREVIEW ------------------
roiLive = rectangle(ax1,'Position',[1 1 1 1],...
    'EdgeColor',[1 0 0],'LineWidth',2,'Visible','off');

set(fig,'WindowButtonMotionFcn',@(~,~) updateMousePos);

function updateMousePos
    if ~isvalid(ax1), return; end
    C = get(ax1,'CurrentPoint');
    mousePos = round([C(1,1) C(1,2)]);
end
%% =========================================================================
%   PART 3 — UPDATE ENGINE, ROI LOGIC, PLAYBACK, HELP, CLEANUP
% =========================================================================

%% ========================= TIMER SETUP ============================
refreshTimer = timer('ExecutionMode','fixedRate','Period',0.03,...
    'TimerFcn',@updateFrame,'BusyMode','drop');
start(refreshTimer);

if dims==4
    playTimer = timer('ExecutionMode','fixedRate','Period',0.06,...
        'TimerFcn',@playStep,'BusyMode','drop');
end

%% =========================================================================
%   MAIN LIVE UPDATE LOOP (runs every 30 ms)
% =========================================================================
function updateFrame(~,~)

    %% ----- Determine slice -----
    if dims == 4
        currentZ = round(get(sliceSlider,'Value'));
    else
        currentZ = 1;
    end

    %% ----- Get display mode -----
    mode = get(displayDropdown,'Value');  % 1 Raw, 2 Normalized, 3 PSC

    switch mode

        % ============================================================
        %   RAW INTENSITY
        % ============================================================
        case 1
            if dims==4
                F = I(:,:,currentZ,1);
            else
                F = I(:,:,1);
            end

            F = normalizeVol(F);
            F = F * contrast.Value + brightness.Value;
            F = max(0,min(1,F));
            F = F .^ gammaSlider.Value;

            ylabel(ax2,'Raw Intensity','Color','w');

        % ============================================================
        %   NORMALIZED
        % ============================================================
        case 2
            F = anatNorm(:,:,currentZ);
            F = F * contrast.Value + brightness.Value;
            F = max(0,min(1,F));
            F = F .^ gammaSlider.Value;

            ylabel(ax2,'Normalized','Color','w');

        % ============================================================
        %   % SIGNAL CHANGE
        % ============================================================
        case 3
            if dims==4
                F = PSC(:,:,currentZ,1);
            else
                F = PSC(:,:,1);
            end

            F = normalizeVol(F);
            ylabel(ax2,'% Signal Change','Color','w');
    end

    %% Apply histogram equalization
    if histEQ.Value == 1
        F = histeq(F);
    end

    %% Rotate + update image
    F = rot90(F,2);
    set(frameH,'CData',F);

    %% Update colormap
    maps = {'gray','hot'};
    colormap(ax1, maps{mapDropdown.Value});

    %% ===============================================================
    %   ROI LIVE PREVIEW + LIVE TIMECOURSE
    % ===============================================================
    x = mousePos(1);
    y = mousePos(2);

    if isnan(x) || x<1 || y<1 || x>Nx || y>Ny
        roiLive.Visible = 'off';
        drawnow limitrate;
        return;
    end

    rs  = round(roiSlider.Value);
    hlf = floor(rs/2);

    x1 = max(1,x-hlf); 
    x2 = min(Nx,x+hlf);
    y1 = max(1,y-hlf); 
    y2 = min(Ny,y+hlf);

    nextColor = roiColors(mod(numel(ROI), size(roiColors,1))+1, :);

    set(roiLive,'Position',[x1 y1 x2-x1+1 y2-y1+1],...
        'Visible','on','EdgeColor',nextColor);

    %% ----- Extract raw + PSC timecourses -----
    if dims==3
        tc_raw = squeeze(mean(mean(I(y1:y2,x1:x2,:),1),2));
        tc_psc = squeeze(mean(mean(PSC(y1:y2,x1:x2,:),1),2));
    else
        tc_raw = squeeze(mean(mean(I(y1:y2,x1:x2,currentZ,:),1),2));
        tc_psc = squeeze(mean(mean(PSC(y1:y2,x1:x2,currentZ,:),1),2));
    end

    tc_raw = tc_raw(:)';
    tc_psc = tc_psc(:)';

    %% Normalized timecourse
    tc_min = min(tc_raw);
    tc_rng = max(tc_raw) - tc_min;
    if tc_rng == 0, tc_rng = 1; end
    tc_norm = (tc_raw - tc_min) / tc_rng;

    %% ----- Select correct timecourse -----
    switch mode
        case 1, tc = tc_raw;
        case 2, tc = tc_norm;
        case 3, tc = tc_psc;
    end

    set(hLive,'YData',tc,'Color',nextColor);

    drawnow limitrate;
end

%% =========================================================================
%   ROI CLICK HANDLER
% =========================================================================
function mouseClick(~,~)

    type = get(fig,'SelectionType');
    x = mousePos(1); 
    y = mousePos(2);

    if isnan(x) || x<1 || x>Nx || isnan(y) || y<1 || y>Ny
        return;
    end

    rs  = round(roiSlider.Value);
    hlf = floor(rs/2);

    x1 = max(1,x-hlf);  x2 = min(Nx,x+hlf);
    y1 = max(1,y-hlf);  y2 = min(Ny,y+hlf);

    % ---------------------- ADD ROI (left click) ----------------------
    if strcmp(type,'normal')

        col = roiColors(mod(numel(ROI), size(roiColors,1))+1, :);

        ROI(end+1) = struct( ...
            'x1',x1,'x2',x2,'y1',y1,'y2',y2,...
            'z',currentZ,'color',col );

        r = rectangle(ax1,'Position',[x1 y1 x2-x1+1 y2-y1+1],...
                      'EdgeColor',col,'LineWidth',2);
        roiHandles(end+1) = r;

        %% ----- Extract TCs for permanent plot -----
        if dims==3
            tc_raw = squeeze(mean(mean(I(y1:y2,x1:x2,:),1),2));
            tc_psc = squeeze(mean(mean(PSC(y1:y2,x1:x2,:),1),2));
        else
            tc_raw = squeeze(mean(mean(I(y1:y2,x1:x2,currentZ,:),1),2));
            tc_psc = squeeze(mean(mean(PSC(y1:y2,x1:x2,currentZ,:),1),2));
        end

        tc_raw = tc_raw(:)';
        tc_psc = tc_psc(:)';

        %% Normalized version
        mn = min(tc_raw);
        rg = max(tc_raw) - mn;
        if rg==0, rg=1; end
        tc_norm = (tc_raw - mn) / rg;

        %% Select correct (consistent!) mode-based TC
        mode = get(displayDropdown,'Value');
        switch mode
            case 1, tc_final = tc_raw;
            case 2, tc_final = tc_norm;
            case 3, tc_final = tc_psc;
        end

        %% Draw permanent ROI curve
        p = plot(ax2,(0:T-1)*TR, tc_final, 'Color', col, 'LineWidth',1.5);
        roiPlots(end+1) = p;

    % ---------------------- REMOVE ROI (right click) -----------------
    elseif strcmp(type,'alt')

        if isempty(ROI), return; end

        centers = zeros(numel(ROI),2);
        for k=1:numel(ROI)
            centers(k,:) = [(ROI(k).x1+ROI(k).x2)/2 , (ROI(k).y1+ROI(k).y2)/2];
        end

        [~, idxMin] = min(sum((centers - [x y]).^2,2));

        delete(roiHandles(idxMin)); 
        delete(roiPlots(idxMin));

        roiHandles(idxMin) = [];
        roiPlots(idxMin)   = [];
        ROI(idxMin)        = [];
    end
end

%% =========================================================================
%   PLAYBACK (4D only)
% =========================================================================
function playStep(~,~)
    if dims ~= 4, return; end
    v = sliceSlider.Value + 1;
    if v > Nz, v = 1; end
    sliceSlider.Value = v;
end

function startPlayback(~,~)
    if dims == 4 && strcmp(playTimer.Running,'off')
        start(playTimer);
    end
end

function stopPlayback(~,~)
    if dims == 4 && strcmp(playTimer.Running,'on')
        stop(playTimer);
    end
end

%% =========================================================================
%   KEYBOARD SHORTCUTS
% =========================================================================
function keyPressHandler(~,event)
    switch event.Key

        case 'rightarrow'
            if dims==4
                sliceSlider.Value = min(Nz, sliceSlider.Value + 1);
            end

        case 'leftarrow'
            if dims==4
                sliceSlider.Value = max(1, sliceSlider.Value - 1);
            end

        case 'uparrow'
            roiSlider.Value = min(150, roiSlider.Value + 1);

        case 'downarrow'
            roiSlider.Value = max(1, roiSlider.Value - 1);

        case 'space'
            if dims==4
                if strcmp(playTimer.Running,'off')
                    start(playTimer);
                else
                    stop(playTimer);
                end
            end
    end
end

%% =========================================================================
%   MOUSE WHEEL = SLICE SCROLL
% =========================================================================
function mouseWheelScroll(~,event)
    if dims ~= 4, return; end
    v = sliceSlider.Value - event.VerticalScrollCount;
    v = max(1,min(Nz,v));
    sliceSlider.Value = v;
end

%% =========================================================================
%   RECALCULATE PSC
% =========================================================================
function recalcPSCfunc(~,~)

    t1 = str2double(baseStart.String);
    t2 = str2double(baseEnd.String);

    if isnan(t1) || isnan(t2) || t1>=t2
        errordlg('Invalid baseline range.');
        return;
    end

    f1 = max(1, round(t1/TR));
    f2 = min(T, round(t2/TR));

    fprintf('Recomputing PSC: frames %d to %d\n',f1,f2);

    if dims==3
        base = mean(I(:,:,f1:f2),3);
        PSC  = (I - base)./base * 100;
    else
        base = mean(I(:,:,:,f1:f2),4);
        PSC  = bsxfun(@rdivide, bsxfun(@minus,I,base),base) * 100;
    end
end

%% =========================================================================
%   HELP WINDOW
% =========================================================================
function showHelpWindow(~,~)
    helpFig = figure('Name','Help & Info','Color','k',...
        'Position',[350 250 850 620]);

    uicontrol(helpFig,'Style','edit',...
        'Units','normalized','Position',[0.05 0.05 0.90 0.90],...
        'Max',20,'Min',1,...
        'BackgroundColor',[0.15 0.15 0.15],...
        'ForegroundColor','w','FontSize',13,...
        'String',{ ...
        '--------------- HELP ---------------';''; ...
        'DISPLAY MODES:'; ...
        '  Raw = absolute signal'; ...
        '  Normalized = dynamic Z-like normalization'; ...
        '  PSC = baseline-subtracted % change';''; ...
        'ROI:'; ...
        '  Mouse move = live preview'; ...
        '  Left click = add ROI'; ...
        '  Right click = remove ROI';''; ...
        'NAVIGATION:'; ...
        '  ?/? = slice'; ...
        '  ?/? = ROI size'; ...
        '  Mouse wheel = slice scroll'; ...
        '  Space = play/pause'; ...
        '';'------------------------------------' });
end

%% =========================================================================
%   CLEANUP (close window)
% =========================================================================
set(fig,'CloseRequestFcn',@cleanup);

function cleanup(~,~)
    try stop(refreshTimer); delete(refreshTimer); end
    if dims==4
        try stop(playTimer); delete(playTimer); end
    end
    delete(fig);
end
end
%% ========================================================================
%% HELPER — NORMALIZATION FUNCTION
%% ========================================================================
function O = normalizeVol(V)
    V = double(V);
    V = V - min(V(:));
    vmax = max(V(:));
    if vmax > 0
        V = V ./ vmax;
    end
    O = V;
end

function export_fUS_NIfTI_v5_2(I, filePath, dx, dy, dz, TR)

    rawNii = fullfile(filePath,'FUS_raw.nii.gz');
    visNii = fullfile(filePath,'FUS_visual.nii.gz');

    if exist(rawNii,'file') && exist(visNii,'file')
        fprintf('NIfTI exists ? skipping export.\n');
        return;
    end

    dims = ndims(I);
    sz = size(I);

    if dims == 3
        I4 = reshape(I,[sz(1) sz(2) 1 sz(3)]);
    else
        I4 = I;
    end

    % orientation fix
    I4 = rot90(I4,-1);
    I4 = flipud(I4);

    % scale to 0–2000
    I4 = single(I4);
    I4 = I4 - min(I4(:));
    if max(I4(:)) > 0
        I4 = I4 ./ max(I4(:));
    end
    I4 = I4 * 2000;

    % ------- Create a temp header using dummy file -------
    tmp = [tempname '.nii'];
    niftiwrite(zeros(2,2,1,'single'), tmp, 'Compressed', false);
    info = niftiinfo(tmp);
    delete(tmp);

    info.ImageSize = size(I4);
    info.PixelDimensions = [dx dy dz TR];
    info.Datatype = 'single';

    fprintf('Saving RAW NIfTI ? %s\n', rawNii);
    niftiwrite(I4, rawNii(1:end-3), info, 'Compressed', true);

    % ------- Create visual NIfTI (pretty anatomy) -------
    Icon = I4;
    Icon = Icon - min(Icon(:));
    Icon = Icon ./ max(Icon(:));
    Icon = Icon .^ 0.6;
    Icon = Icon * 2000;

    fprintf('Saving VISUAL NIfTI ? %s\n', visNii);
    niftiwrite(Icon, visNii(1:end-3), info, 'Compressed', true);
end
