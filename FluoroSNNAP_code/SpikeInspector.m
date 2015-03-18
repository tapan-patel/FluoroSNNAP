classdef SpikeInspector < handle
    % A class for inspecting spikes from imaging data - adding spikes,
    % removing spikes, etc.
    
    events
        SaveClicked % Thrown when apply button is clicked
    end
    
    properties % Public variables
        
        mask % label matrix
        dF_cell
        F_cell
        Spikes_cell
        N % total number of cells
        s % processed analysis
        frames % total number of frames
        figh = 300; % initial figure height - image is scaled to fit
        % on change of this the window gets resized
        buttons; % handle to uicontrol buttons
        checkbox; % handle to checkbox
        spikeax % fluoresence trace with spikes
        neuron
        spikefig
        spikes
    end
    
    properties(Access=private) % Private variables
        % UI stuff
        guifig % main window
        imax % working image
        
        imag
        
        tl % userinfo bar
        image
        figw % initial window height
        hwar = 2.1 % aspect ratio
        
        % Class stuff
        current % which spike is selected
        
        shapes = {}; % hold all rois
        spike_lines = {}; % hold all spikes for a given cell
        
        % load/save info
        filename
        spike_label;
        thr
    end
    
    %% Public methods
    methods
        
        function this = SpikeInspector(s)
            if(nargin > 0)
                
                % Set public variables
                this.s = s;
                if(isfield(s,'dF_cell'))
                    this.dF_cell = s.dF_cell;
                    this.F_cell = s.F_cell;
                else
                    this.dF_cell = s.F_cell;
                    this.F_cell = s.F_cell;
                    this.s.dF_cell = s.F_cell;
                end
                
                this.N = size(this.s.dF_cell,1);
                this.frames = size(this.s.dF_cell,2);
                this.filename = this.s.filename;
                this.Spikes_cell = this.s.Spikes_cell;
            else
                error('Expected 1 argument, found 0');
            end
            % Constructor. Input is a processed_analysis
            this.figw = this.figh*this.hwar;
            
            % invoke the UI window
            this.createWindow;
            try
                % Try to load the max projection image, if not load first
                % frame, if not, throw an error message
                [folder file] = fileparts(this.filename);
                if(~isfield(this.s,'image'))
                    if(exist([folder '/MAX-' file '.tif']))
                        fullfile = [folder '/MAX-' file '.tif'];
                        I = imread(fullfile);
                    elseif(exist(this.filename,'file'))
                        I = imread(this.filename);
                    else
                        I = this.s.L;
                    end
                    if(size(I,3) ~=3)
                        image = zeros([size(I) 3]);
                        image(:,:,2) = imadjust(I);
                        image = image./max(image(:));
                        this.s.image = image;
                    end
                elseif(size(this.s.image,3)~=3)
                    image = zeros([size(this.s.image) 3]);
                    image(:,:,2) = imadjust(this.s.image);
                    image = image./max(image(:));
                    this.s.image = image;
                end
                this.resetImages;
                this.resizeWindow;
            catch
                error(['Image could not be read from ' this.filename]);
            end
            tmp = load('spikes.mat');
            this.spikes = tmp.spikes;
            msg = sprintf('Loading ROI info. Please wait.\nThis message will self destruct');
            wbh = msgbox(msg,'Please wait');
            delete(findobj(wbh,'string','OK')); drawnow;
            this.current = 1;
            this.neuron = 1;
            % populate shapes with ROIs
            C = regionprops(this.s.L,'BoundingBox','Centroid');
            this.shapes = cell(this.N,1);
            for i=1:numel(C)
                this.shapes{i} = imellipse(this.imax,C(i).BoundingBox);
                setResizable(this.shapes{i}, false);
                fcn = makeConstrainToRectFcn('imellipse',...
                    [C(i).BoundingBox(1) C(i).BoundingBox(1)+C(i).BoundingBox(3)],...
                    [C(i).BoundingBox(2) C(i).BoundingBox(2)+C(i).BoundingBox(4)]);
                setPositionConstraintFcn(this.shapes{i},fcn);
                set(this.shapes{i},'Tag',['imsel_', num2str(i)]);
                % label each neuron
                %                 text(this.imax,C(i).Centroid(1)-5,C(i).Centroid(2)-5,num2str(i),'Color','m');
            end
            this.filename = this.s.filename;
            if(isfield(this.s,'thr') && ~isempty(this.s.thr))
                this.thr = this.s.thr;
            else
                this.thr = .85*ones(this.N,1);
            end
            delete(wbh);
        end
        
        function delete(this)
            % destructor
            delete(this.guifig);
        end
        
        function set.image(this,I)
            % set method for image - read in the first frame in image for
            % roi selection
            
        end
        
        function set.figh(this,height)
            this.figh = height;
            this.figw = this.figh*this.hwar;
            this.resizeWindow;
        end
    end
    %% Private methods
    methods(Access=private)
        function resetImages(this)
            
            % load image
            h = zoom;
            set(h,'Motion','horizontal','Enable','on');
            pan
            this.spikefig = plot(this.s.F_whole,'parent',this.spikeax);
            
            set(this.spike_label,'String','Full frame intensity');
            xlabel(this.spikeax,'Frame #'); ylabel(this.spikeax,'Fluorescence (a.u.)');
            ylim(this.spikeax,'auto');
            this.imag = imshow((this.s.image),'parent',this.imax);
            
        end
        function winpressed(this,h,e,type)
            switch type
                case 'down'
                    SelObj = get(gco,'Parent');
                    Tag = get(SelObj,'Tag');
                    if and(~isempty(SelObj),strfind(Tag,'imsel_'))
                        this.current = str2double(Tag(7:end));
                        this.neuron = this.current;
                        for i=1:numel(this.shapes)
                            if i==this.current
                                setColor(this.shapes{i},'red');
                                updateSpikes(this,i);
                            else
                                setColor(this.shapes{i},'blue');
                            end
                        end
                    end
                    
                    if and(~isempty(SelObj),strfind(Tag,'spike_'))
                        this.current = str2double(Tag(7:end));
                        for i=1:numel(this.spike_lines)
                            if i==this.current
                                setColor(this.spike_lines{i},'w');
                            end
                        end
                        for i=this.current+1:numel(this.spike_lines)
                            set(this.spike_lines{i},'Tag',['spike_', num2str(i-1)]);
                        end
                        this.spike_lines(this.current)=[];
                        spks = this.Spikes_cell{this.neuron};
                        
                        % ask if the spike waveform should be removed from
                        % the library. Need to figure out all transients
                        % that have good similarity to the current spike
                        button = questdlg('Do you want to remove the waveform that auto-detected this spike from the database? This will prevent the automated detection of other spikes that look similar to the current spike being deleted','Remove from database?','Yes','No','No');
                        if(strcmp(button,'Yes'))
                            similarity = zeros(size(this.spikes));
                            for k = 1:length(similarity)
                                R = corrcoef(this.spikes{k},this.dF_cell(this.neuron,spks(this.current):spks(this.current)+length(this.spikes{k})-1));
                                similarity(k) = R(1,2);
                                
                            end
                            
                            id = find(similarity>this.thr(this.neuron));
                            for j = 1:length(id)
                                figure; plot(this.spikes{id(j)});
                            title('This spike waveform in the database will be removed');
                            end
                            button2 = questdlg('Confirm deletion of this spike waveform from the database?','Confirm deletion','Yes','No','No');
                            if(strcmp(button2,'Yes'))
                                this.spikes(id) = [];
                                spikes = this.spikes;
                                save('spikes.mat','spikes');
                                
                            end
                            
                        end
                        spks(this.current) = [];
                        
                        this.Spikes_cell{this.neuron} = spks;
                        %                 updateSpikes(this,this.current_cell);
                    end
                    
            end
        end
        
        function closefig(this,h,e)
            delete(this);
        end;
        
        function updateSpikes(this,i)
            set(this.spike_label,'String',['Trace for cell ' num2str(i)]);
            set(this.buttons(end),'Value',this.thr(i));
            set(this.tl,'String',num2str(this.thr(i)));
            this.spikefig = plot(this.s.dF_cell(i,:),'parent',this.spikeax);
            xlabel(this.spikeax,'Frame #'); ylabel(this.spikeax,'Fluorescence (a.u.)');
            % overlay detected spikes
            this.spike_lines = cell(numel(this.Spikes_cell{i}),1);
            spks = this.Spikes_cell{i};
            spks = sort(spks);
            if(get(this.checkbox,'Value'))
            for j=1:numel(spks)
                axes(this.spikeax)
                yl = ylim(this.spikeax);
                position = [spks(j) yl(1); spks(j) yl(2)];
                
                this.spike_lines{j} = impoly(this.spikeax,position);
                set(this.spike_lines{j},'Tag',['spike_', num2str(j)]);
                setColor(this.spike_lines{j},'black');
                setVerticesDraggable(this.spike_lines{j}, false);
                fcn = makeConstrainToRectFcn('impoly',[min(position(:,1)) max(position(:,1))],[min(position(:,2)) max(position(:,2))]);
                setPositionConstraintFcn(this.spike_lines{j},fcn);
                
                
                %                 hold on
                %                 plot([this.Spikes_cell{i}(j) this.Spikes_cell{i}(j)],yl,'Color','k','LineWidth',2);
                
            end
            end
        end
        
        function slider_cb(this)
            thr = get(this.buttons(end),'Value');
            set(this.tl,'String',['Spike Detection threshold = ' num2str(thr)]);
            msg = sprintf('Updating spikes');
            wbh = msgbox(msg,'Please wait');
            delete(findobj(wbh,'string','OK')); drawnow;
            spks = SpikesLib(this.dF_cell(this.neuron,:),this.spikes,'thr',thr);
            this.Spikes_cell{this.neuron} = spks;
            this.thr(this.neuron) = thr;
            updateSpikes(this,this.neuron);
            close(wbh);
        end
        
        function updateclick(this,h,e)
            % Updates spike detection for all cells using threshold from
            % the slider
            thr = get(this.buttons(end),'Value');
            wb = waitbar(0,['Updating spikes for all neurons using similarity threshold ' num2str(thr)]);
            for i=1:this.N
                spks = SpikesLib(this.dF_cell(i,:),this.spikes,'thr',thr);
                this.Spikes_cell{i} = spks;
                this.thr(i) = thr;
                waitbar(i/this.N);
            end
            close(wb)
            set(this.buttons(2),'Selected','off');
        end
        function addclick(this,h,e)
            % Allow user to click and add spike
            try
                % Current spikes
                spks = this.Spikes_cell{this.neuron};
                % Allow user to select data point
                pt = impoint(this.spikeax);
                vt = pt.getPosition;
                
                % Add a spike to this position
                spks = [spks floor(vt(1))];
                this.Spikes_cell{this.neuron} = spks;
                % Ask user if this spike waveform should be added to the
                % library
                button = questdlg('Do you want to add this spike waveform to the database for improved detection?','Update database?','Yes','No','No');
                switch button
                    case 'Yes'
                        
                        helpdlg('You already selected the onset of the spike. Click on the fluorescence trace and select the endpoint for this transient','Select endpoint for transient');
                        
                        pt = impoint(this.spikeax);
                        vt_end = pt.getPosition;
                        waveform = this.F_cell(this.neuron,floor(vt):floor(vt_end));
                        hw= figure; plot(waveform); title('New waveform added to database. Similar transients will now be detected');
                        button2 = questdlg('Click confirm to add this waveform to the database','Update database?','Confirm','Reject','Confirm');
                        if(strcmp(button2,'Confirm'))
                            this.spikes{end+1} = waveform;
                            spikes = this.spikes;
                            save('spikes.mat','spikes');
                        end
                        close(hw);
                end
                
                updateSpikes(this, this.neuron);
            catch
                uicontrol(this.buttons(1));
            end
            set(this.buttons(1),'Selected','off');
        end
        
        function rasterclick(this,h,e)
            % Plot the rasters
            try
                rfig = figure;
                hold on
                for i=1:this.N
                    spks = this.Spikes_cell{i};
                    if(~isempty(spks))
                        plot(spks/this.s.fps,i,'b.');
                    end
                end
                xlabel('Time (s)'); ylabel('Neuron ID');
                title(['Activity for ' this.s.filename]);
                hold off
            catch
                uicontrol(this.buttons(3));
            end
            set(this.buttons(3),'Selected','off');
        end
        function checkclick(this,h,e)
            updateSpikes(this,this.neuron);
        end
        function scaclick(this,h,e)
            try
                choice = questdlg('Do you want to use pre-processed data or recompute?','Choice','Pre-processed','Recompute','Pre-processed');
                switch choice
                    case 'Recompute'
                        
                        % Show functional connectivity, modularity, rearranged rasters
                        s = this.s;
                        s.Spikes_cell = this.Spikes_cell;
                        L = s.L;
                        N = this.N;
                        % Reperform SCA
                        msg = sprintf('Determining functional connectivity, synchronization and calcium transient kinetics.\nThis may take several minutes.\nThis message will self-destruct');
                        wbh = msgbox(msg,'Please wait');
                        delete(findobj(wbh,'string','OK')); drawnow;
                        try
                            s.phase = [];
                            s.SynchroCluster = [];
                        end
                        C = SCA(s,'DISPLAY',0,'save_flag',0,'wb',0,'analysis_type','phase');
                        this.s.SynchroCluster = C;
                        A = C.A;
                        c = regionprops(L,'Centroid');
                        c = reshape([c.Centroid],2,N);
                        x = c(1,:);
                        y = c(2,:);
                        [DF,rise_time,fall_time,CV] = GetTransients(this.s);
                        this.s.DF = DF;
                        this.s.rise_time = rise_time;
                        this.s.fall_time = fall_time;
                        this.s.CV = CV;
                        
                        [Ci,Q] = modularity_und(A);
                        this.s.modules = Ci;
                        this.s.modularity = Q;
                        
                        % Sychronization by modules
                        SI = zeros(max(Ci),1);
                        for k=1:max(Ci)
                            tmp_s.F_cell = this.s.F_cell(Ci==k,:);
                            tmp_s.Spikes_cell = this.s.Spikes_cell(Ci==k);
                            tmp_s.fps = this.s.fps;
                            c_module = SCA(tmp_s);
                            SI(k) = max(c_module.SI);
                        end
                        this.s.SI_m = SI;
                        Colors = varycolor(max(Ci));
                        
                    case 'Pre-processed'
                        s = this.s;
                        try
                            L = s.L;
                            N = this.N;
                            C = s.SynchroCluster;
                            A = C.A;
                            c = regionprops(L,'Centroid');
                            c = reshape([c.Centroid],2,N);
                            x = c(1,:);
                            y = c(2,:);
                            DF = s.DF; rise_time = s.rise_time; fall_time = s.fall_time; CV = s.CV;
                            Ci = s.modules; Q = s.modularity;
                            SI = s.SI_m;
                            Colors = varycolor(max(Ci));
                        catch
                            errordlg('Could not use pre-processed information. Missing entries. Select "Reprocess" option');
                            return;
                        end
                end
                clu = clustering_coef_bu(A);
                netfig = figure;
                
                hdt = datacursormode(netfig);
                set(hdt,'DisplayStyle','window');
                
                set (netfig, 'Units', 'normalized', 'Position', [0,0,1,1]);
                h1 = subplot(2,2,1);
                set(hdt,'UpdateFcn',{@labeltip,C,c,this.Spikes_cell,DF,...
                    rise_time,fall_time,clu,h1});
                if(size(this.s.image,3)~=3)
                    image = zeros([size(this.s.image) 3]);
                    image(:,:,2) = imadjust(this.s.image);
                else
                    image = this.s.image;
                end
                imagesc(image,'Parent',h1);
                hold on
                
                gplot(A,[x' y']);
                
                set(findobj(gcf,'Type','Line'),'LineStyle',':','LineWidth',.01,'MarkerSize',6,'Marker','o','MarkerFaceColor',[1 1 0],'Color',[1 0 0]);
                
                hold on
                plot(x(1,~logical(sum(A)) & ~logical(sum(A'))),y(1,~logical(sum(A)) & ~logical(sum(A'))),'ro','MarkerFaceColor',[1 0 0])
                set(gca,'YDir','reverse');
                set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
                title({['N = ' num2str(N) ',Dropped = ' num2str(sum(~logical(sum(A)) & ~logical(sum(A'))))]});
                
                
                subplot(2,2,2);
                [~,pos] = max(C.PI,[],2);
                [~,IDX] = sort(pos);
                imagesc(C.C(IDX,IDX),[0 1]); axis square; colorbar;
                xlabel('Neuron ID (rearranged)'); ylabel('Neuron ID (rearranged)');
                msg = sprintf('Global synchronization index: %f\n%d clusters detected',max(C.SI),size(C.PI,2));
                title(msg);
                
                subplot(2,2,3);
                
                set(gca, 'ColorOrder', Colors);
                hold all;
                for i=1:max(Ci)
                    scatter(x(Ci==i),y(Ci==i),'filled','SizeData',10^2);
                end
                set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
                set(gca,'YDir','reverse');
                
                %                 title({sprintf('N = %d. Modules = %d, Modularity = %f',length(x),max(Ci),Q)});
                
                title_string = sprintf('N = %d. Modules = %d. Modularity = %f\n',length(x),max(Ci),Q);
                title_string = sprintf('%sSynchronization by modules: ',title_string);
                for k=1:length(SI)
                    title_string = sprintf('%s %f,',title_string,SI(k));
                end
                title(title_string);
                subplot(2,2,4);
                t = 0:1/s.fps:size(s.F_cell,2)/s.fps-1/s.fps;
                hold all;
                cntr = 1;
                for i=1:max(Ci)
                    modules = find(Ci==i);
                    for j=1:length(modules)
                        try
                            spks = s.Spikes_cell{modules(j)};
                            plot(spks/s.fps,cntr,'.','Color',Colors(i,:));
                            cntr = cntr+1;
                        end
                    end
                end
                title(['Activity rearranged by modules']);
                delete(wbh);
            catch
                uicontrol(this.buttons(4));
            end
            
            set(this.buttons(4),'Selected','off');
            
            
            function out_text = labeltip(h,e,C,c,S,DF,rise_time,fall_time,clu,h1)
                N = size(C.A,1);
                A = C.A;
                pos = get(e,'Position');
                neuron = find(c(1,:)==pos(1) & c(2,:)==pos(2));
                if(isempty(neuron))
                    out_text = '';
                    plot(c(1,:),c(2,:),'ro','MarkerSize',6,'Marker','o','MarkerFaceColor',[1 1 0],'Color',[1 0 0]);
                    plot(c(1,~logical(sum(A)) & ~logical(sum(A'))),c(2,~logical(sum(A)) & ~logical(sum(A'))),'ro','MarkerFaceColor',[1 0 0])
                else
                    
                    % Reset the size and color of all neurons to the way it was
                    % originally
                    
                    %                     cla(h1);
                    %                    imagesc(image,'Parent',h1);
                    %                 hold on
                    %
                    %                 gplot(A,[x' y']);
                    %
                    %                 set(findobj(gcf,'Type','Line'),'LineStyle',':','LineWidth',.01,'MarkerSize',6,'Marker','o','MarkerFaceColor',[1 1 0],'Color',[1 0 0]);
                    %                 set(gca,'XTickLabel',[]);
                    %                 set(gca,'YTickLabel',[]);
                    %                 hold on
                    %                 plot(x(1,~logical(sum(A)) & ~logical(sum(A'))),y(1,~logical(sum(A)) & ~logical(sum(A'))),'ro','MarkerFaceColor',[1 0 0])
                    %                 set(gca,'YDir','reverse');
                    %                 set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
                    %                 title({['N = ' num2str(N) ',Dropped = ' num2str(sum(~logical(sum(A)) & ~logical(sum(A'))))]});
                    
                    
                    out_text = ['Neuron #: ' num2str(neuron)];
                    out_text = sprintf('%s. CI: %.02f',out_text,sum(A(neuron,:))/N);
                    out_text = sprintf('%s,   PI: %.02f,  CC: %.02f',out_text,max(C.PI(neuron,:)),clu(neuron));
                    
                    out_text = sprintf('%s. Events: %d, Amplitude: %.02f, <t_rise> = %.02f (s), <t_fall> = %.02f (s)',...
                        out_text,length(S{neuron}),DF(neuron),mean(rise_time{neuron}),...
                        mean(fall_time{neuron}));
                    
                    connecting = find(A(neuron,:));
                    hold on
                    plot(c(1,:),c(2,:),'ro','MarkerSize',6,'Marker','o','MarkerFaceColor',[1 1 0],'Color',[1 0 0]);
                    plot(c(1,~logical(sum(A)) & ~logical(sum(A'))),c(2,~logical(sum(A)) & ~logical(sum(A'))),'ro','MarkerFaceColor',[1 0 0])
                    plot(c(1,connecting),c(2,connecting),'bo','MarkerFaceColor','b','MarkerSize',6)
                    
                end
            end
        end
        
        
        function saveclick(this,h,e)
            this.s.Spikes_cell = this.Spikes_cell;
            this.s.thr = this.thr;
            notify(this,'SaveClicked');
        end
        
        function summaryclick(this,h,e)
            % Display summary statistics and offer to save them to a .txt
            % file
            choice = questdlg('Do you want to use pre-processed data or recompute?','Choice','Pre-processed','Recompute','Pre-processed');
            switch choice
                case 'Recompute'
                    
                    this.s.Spikes_cell = this.Spikes_cell;
                    msg = sprintf('Gathering summary statistics.\nThis may take several minutes.\n');
                    wbh = msgbox(msg,'Please wait');
                    delete(findobj(wbh,'string','OK')); drawnow;
                    
                    msg = sprintf('Summary for file: %s\n',this.filename);
                    msg = sprintf('%s\nTotal neurons = %d\n',msg,this.N);
                    msg = sprintf('%s\nFrames = %d, frame rate = %d, duration = %.02f (s)\n',msg,this.frames,this.s.fps,this.frames/this.s.fps);
                    %             if(isfield(this.s,'SynchroCluster') & ~isempty(this.s.SynchroCluster))
                    %                 C = this.s.SynchroCluster;
                    %             else
                    
                    C = SCA(this.s,'DISPLAY',0,'save_flag',0,'wb',0,'analysis_type','phase');
                    this.s.SynchroCluster = C;
                    %             end
                    SI = max(C.SI);
                    A = C.A;
                    [Ci,Q] = modularity_und(A);
                    d = sum(A);
                    
                    
                    % Sychronization by modules
                    SI_m = zeros(max(Ci),1);
                    for k=1:max(Ci)
                        tmp_s.F_cell = this.s.F_cell(Ci==k,:);
                        tmp_s.Spikes_cell = this.s.Spikes_cell(Ci==k);
                        tmp_s.fps = this.s.fps;
                        c = SCA(tmp_s);
                        SI_m(k) = max(c.SI);
                    end
                    
                    msg = sprintf('%s\nModularity: %.04f, Number of modules: %d\n',msg,Q,max(Ci(:)));
                    msg = sprintf('%s\nGlobal synchronization index: %.02f\n',msg,SI);
                    msg = sprintf('%s\nSynchronization index by modules: ',msg);
                    
                    for k=1:length(SI_m)
                        msg = sprintf('%s %f,',msg,SI_m(k));
                    end
                    
                    msg = sprintf('%s\n\nMean connectivity: %.04f\n',msg,mean(d)/this.N);
                    % Mean oscillation period (s)
                    ISI = cell(this.N,1);
                    for i=1:this.N
                        ISI{i} = diff(this.Spikes_cell{i});
                    end
                    ISI = [ISI{:}];
                    OP = mean(ISI)/this.s.fps;
                    msg = sprintf('%s\nMean oscillation period: %.02f (s)\n',msg,OP);
                    
                    % Determine deltaF/F0 for all neurons, show mean and 95% CI
                    
                    [DF,rise_time,fall_time,CV] = GetTransients(this.s);
                    ci = bootci(40,@(x) mean(x),DF(~isnan(DF)));
                    msg = sprintf('%s\nMean amplitude: %.04f, 95%% CI: [%.04f, %.04f], CV: %f\n',msg,mean(DF(~isnan(DF))),ci(1),ci(2),mean(CV(~isnan(CV))));
                    
                    Rise = [rise_time{:}]; Fall = [fall_time{:}];
                    ci_r = bootci(40,@(x) mean(x),Rise(~isnan(Rise)));
                    msg = sprintf('%s\nMean rise time: %.02f (s), 95%% CI: [%.02f, %.02f], CV: %f\n',msg,mean(Rise(~isnan(Rise))),ci_r(1),ci_r(2), std(Rise(~isnan(Rise)))./mean(Rise(~isnan(Rise))));
                    
                    ci_tau = bootci(40,@(x) mean(x),Fall(~isnan(Fall)));
                    msg = sprintf('%s\nMean fall time: %.02f (s), 95%% CI: [%.02f, %.02f], CV: %f\n',msg,mean(Fall(~isnan(Fall))),ci_tau(1),ci_tau(2),std(Fall(~isnan(Fall)))./mean(Fall(~isnan(Fall))));
                    close(wbh);
                case 'Pre-processed'
                    C = this.s.SynchroCluster;
                    SI = max(C.SI);
                    A = C.A;
                    Ci = this.s.modules;
                    Q = this.s.modularity;
                    d = sum(A);
                    SI_m = this.s.SI_m;
                    DF = this.s.DF;
                    rise_time = this.s.rise_time;
                    fall_time = this.s.fall_time;
                    CV = this.s.CV;
                    msg = sprintf('Summary for file: %s\n',this.filename);
                    msg = sprintf('%s\nTotal neurons = %d\n',msg,this.N);
                    msg = sprintf('%s\nFrames = %d, frame rate = %d, duration = %.02f (s)\n',msg,this.frames,this.s.fps,this.frames/this.s.fps);
                    msg = sprintf('%s\nModularity: %.04f, Number of modules: %d\n',msg,Q,max(Ci(:)));
                    msg = sprintf('%s\nGlobal synchronization index: %.02f\n',msg,SI);
                    msg = sprintf('%s\nSynchronization index by modules: ',msg);
                    for k=1:length(SI_m)
                        msg = sprintf('%s %f,',msg,SI_m(k));
                    end
                    msg = sprintf('%s\n\nMean connectivity: %.04f\n',msg,mean(d)/this.N);
                    % Mean oscillation period (s)
                    ISI = cell(this.N,1);
                    for i=1:this.N
                        ISI{i} = diff(this.Spikes_cell{i});
                    end
                    ISI = [ISI{:}];
                    OP = mean(ISI)/this.s.fps;
                    msg = sprintf('%s\nMean oscillation period: %.02f (s)\n',msg,OP);
                    ci = bootci(40,@(x) mean(x),DF(~isnan(DF)));
                    msg = sprintf('%s\nMean amplitude: %.04f, 95%% CI: [%.04f, %.04f], CV: %f\n',msg,mean(DF(~isnan(DF))),ci(1),ci(2),mean(CV(~isnan(CV))));
                    Rise = [rise_time{:}]; Fall = [fall_time{:}];
                    ci_r = bootci(40,@(x) mean(x),Rise(~isnan(Rise)));
                    msg = sprintf('%s\nMean rise time: %.02f (s), 95%% CI: [%.02f, %.02f], CV: %f\n',msg,mean(Rise(~isnan(Rise))),ci_r(1),ci_r(2), std(Rise(~isnan(Rise)))./mean(Rise(~isnan(Rise))));
                    
                    ci_tau = bootci(40,@(x) mean(x),Fall(~isnan(Fall)));
                    msg = sprintf('%s\nMean fall time: %.02f (s), 95%% CI: [%.02f, %.02f], CV: %f\n',msg,mean(Fall(~isnan(Fall))),ci_tau(1),ci_tau(2),std(Fall(~isnan(Fall)))./mean(Fall(~isnan(Fall))));
            end
            % Create per neuron report in an table/
            % Number of events, <ISI>, <DF>, rise time, fall time, CI, PI, # assemblies,
            % clustering coef
            dat = zeros(this.N,10);
            clu = clustering_coef_bu(C.A);
            for i=1:this.N
                x = sort(this.s.F_cell(i,:));
                dat(i,1) = mean(x(1:ceil(.1*length(x))));
                dat(i,2) = numel(this.Spikes_cell{i});
                dat(i,3) = mean(diff(this.Spikes_cell{i}))/this.s.fps;
                dat(i,4) = DF(i);
                dat(i,5) = CV(i);
                r = rise_time{i}; tau = fall_time{i};
                dat(i,6) = mean(r(~isnan(r)));
                dat(i,7) = mean(tau(~isnan(tau)));
                dat(i,8) = sum(A(i,:))/this.N;
                dat(i,9) = max(C.PI(i,:));
                dat(i,10) = nnz(C.PI(i,:)>.01);
                dat(i,11) = clu(i);
            end
            
            button = questdlg(msg,'Summary','OK','Cell based report','Export to .txt file','OK');
            
            switch button
                case 'Export to .txt file'
                    % Ask where to save the file
                    [outname, outfolder] = uiputfile('*.txt','Save file name',[this.filename(1:end-4) '_summary.txt']);
                    if(length(outname)>1 && length(outfolder)>1) % Error-checking - if user did not select a name, both outname and outfolder = 0, length = 1
                        fid = fopen(fullfile(outfolder,outname),'w');
                        fprintf(fid,'Summary for file:\t %s\n\n',this.filename);
                        fprintf(fid,'Total cells (ROIs):\t %d\n',this.N);
                        fprintf(fid,'Frames:\t %d\nframe rate:\t %d\nduration (s):\t %.02f\n',...
                            this.frames,this.s.fps,this.frames/this.s.fps);
                        fprintf(fid,'Modularity:\t %.04f\nNumber of modules:\t %d\n',Q,max(Ci(:)));
                        fprintf(fid,'Global synchronization index:\t %.02f\n',SI);
                        fprintf(fid,'Synchronization index by modules:\t ');
                        for k=1:length(SI_m)
                            fprintf(fid,'%f\t',SI_m(k));
                        end
                        fprintf(fid,'\n');
                        fprintf(fid,'Mean global functional connectivity:\t %.04f\n',mean(d)/this.N);
                        fprintf(fid,'Mean oscillation period (s):\t %.02f\n',OP);
                        fprintf(fid,'Mean amplitude (deltaF/F0):\t %.04f\n95%% CI:\t%.04f\t %.04f\n',mean(DF(~isnan(DF))),ci(1),ci(2));
                        fprintf(fid,'Coefficient of variation:\t%f\n',std(CV(~isnan(CV)))./mean(CV(~isnan(CV))))
                        fprintf(fid,'Mean rise time (s):\t %.02f\n95%% CI:\t%.02f\t%.02f\n',mean(Rise(~isnan(Rise))),ci_r(1),ci_r(2));
                        fprintf(fid,'Mean fall time (s):\t %.02f\n95%% CI:\t%.02f\t%.02f\n',mean(Fall(~isnan(Fall))),ci_tau(1),ci_tau(2));
                        %Now output cell based report
                        fprintf(fid,'\n\nNeuron ID\tBaseline fluorescence\tTotal events\t<ISI> (s)\tAmplitude\tCV\t<Rise time> (s)\t<Fall time> (s)\tConnectivity\tParticipation\t# assemblies\tClustering coef\n');
                        for i=1:this.N
                            fprintf(fid,'%d\t%f\t%.02f\t%.04f\t%.02f\t%f\t%.02f\t%.04f\t%.04f\t%d\t%.04f\t%0.04f\n',...
                                i,dat(i,1),dat(i,2),dat(i,3),dat(i,4),dat(i,5),dat(i,6),dat(i,7),dat(i,8),dat(i,9),dat(i,10),dat(i,11));
                        end
                        fclose(fid);
                        msgbox(['Summary output to ' fullfile(outfolder,outname)],'Write complete');
                    else
                        return
                    end
                case 'Cell based report'
                    columnname = {'Baseline fluorescence','Total events','<ISI> (s)','Amplitude','CV','Rise time (s)','Fall time (s)','Connectivity index','Participation index','# assemblies','Clustering coef'};
                    columneditable = [false false false false false false false false false false false];
                    f1 = figure;
                    t = uitable('Units','normalized','Position',...
                        [0.1 0.1 0.9 0.9], 'Data', dat,...
                        'ColumnName', columnname,...
                        'ColumnEditable', columneditable);
            end
        end
        function movieclick(this,h,e)
            
            choice = questdlg('Movie type','Make a selection','Simple (fast)','Elegant (slow)','Simple');
            switch choice
                case 'Elegant (slow)'
                    % make a movie
                    msg = sprintf('Creating an avi movie.\nThis will take approximately %d minutes.',ceil((.7*this.frames)/60));
                    wbh = waitbar(0,msg,'Name','Please wait','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
                    
                    vid = VideoWriter([this.filename(1:end-4) '.avi']);
                    vid.FrameRate = this.s.fps*2;
                    open(vid);
                    %             aviobj = avifile([this.filename(1:end-4) '.avi'],'fps',this.s.fps*2);
                    
                    [~,idx] = min(this.s.F_whole);
                    % idx is the baseline image
                    info = imfinfo(this.filename);
                    Iref = imread(this.filename,'Info',info,'Index',idx);
                    % Find the max frame and adjust image intensity for every frame
                    % such that the max frame is saturating
                    [~,idx] = max(this.s.F_whole);
                    MaxI = imread(this.filename,'Info',info,'Index',idx);
                    MaxIrgb = zeros([size(MaxI) 3],class(MaxI));
                    MaxIrgb(:,:,2) = MaxI;
                    Irefrgb = zeros([size(Iref) 3],class(Iref));
                    Irefrgb(:,:,2) = Iref;
                    lo_high = stretchlim(Irefrgb);
                    
                    % Also determine what the max DF/F0 is
                    DF_max = (double(MaxI)-double(Iref))./double(Iref);
                    high_color = max(DF_max(:));
                    high_color = min([high_color 2.5]);
                    
                    fh = figure('Visible','off');
                    set (fh, 'Units', 'normalized', 'Position', [0,0,.9,.9]);
                    fh1 = subplot(1,2,1);
                    fh2 = subplot(1,2,2);
                    set(fh,'PaperPositionMode','auto');
                    for i=1:numel(info)
                        
                        if(getappdata(wbh,'canceling'))
                            break;
                        end
                        waitbar(i/numel(info),wbh,sprintf('%.2f%% complete.',i/numel(info)*100));
                        I = imread(this.filename,'Index',i,'Info',info);
                        Irgb = zeros([size(I) 3],class(I));
                        Irgb(:,:,2) = I;
                        Irgb = imadjust(Irgb,lo_high);
                        DF = (double(I)-double(Iref))./double(Iref);
                        imshow(Irgb,'Parent',fh1);
                        title(fh1,['t = ' num2str(i/this.s.fps) ' (s)']);
                        image(medfilt2(DF),'Parent',fh2,'CDataMapping','scaled');
                        set(fh2,'CLim',[0 high_color]);
                        axis([fh1 fh2],'image');
                        colorbar(); set(fh2,'XTickLabel',[]); set(fh2,'YTickLabel',[]);
                        title(fh2,'\DeltaF/F_0');
                        print(fh,'temp_render.jpeg','-djpeg');
                        R = imread('temp_render.jpeg');
                        writeVideo(vid,R);
                        %                 aviobj = addframe(aviobj,fh);
                    end
                    %             aviobj = close(aviobj);
                    close(vid);
                    
                    delete(fh);
                    delete(wbh);
                case 'Simple (fast)'
                    msg = sprintf('Creating an avi movie.\nThis will take approximately %d minutes.',ceil((.05*this.frames)/60));
                    wbh = waitbar(0,msg,'Name','Please wait','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
                    
                    vid = VideoWriter([this.filename(1:end-4) '.avi']);
                    vid.FrameRate = this.s.fps*2;
                    open(vid);
                    info = imfinfo(this.filename);
                    I = imread(this.filename);
                    hi_low = stretchlim(im2double(I));
                    for i=1:numel(info);
                        
                        if(getappdata(wbh,'canceling'))
                            break;
                        end
                        waitbar(i/numel(info),wbh,sprintf('%.2f%% complete.',i/numel(info)*100));
                        I = imread(this.filename,'Index',i,'Info',info);
                        I = imadjust(im2double(I),hi_low);
                        writeVideo(vid,I);
                        
                    end
                    close(vid);
                    delete(wbh);
            end
            set(this.buttons(5),'Selected','off');
            msgbox(['Video successfuly created and saved to ' this.filename(1:end-4) '.avi']);
            winopen([this.filename(1:end-4) '.avi']);
        end
        function keyPress(this,h,e)
            
            switch e.Key
                case 'a'
                    set(this.buttons(1),'Selected','on')
                    addclick(this,h,e);
                    
                case 'u'
                    set(this.buttons(2),'Selected','on')
                    updateclick(this,h,e);
                case 'r'
                    set(this.buttons(3),'Selected','on')
                    rasterclick(this,h,e);
                case 's'
                    set(this.buttons(4),'Selected','on')
                    scaclick(this,h,e);
                case 'n'
                    set(this.buttons(5),'Selected','on')
                    networkclick(this,h,e);
            end
        end
        % UI functions
        function createWindow(this, w, h)
            this.guifig=figure('MenuBar','none','Resize','on','Toolbar','none','Name','Spike Inspector', ...
                'NumberTitle','off','Color','white', 'units','pixels','position',[0 0 this.figw this.figh],...
                'CloseRequestFcn',@this.closefig, 'KeyPressFcn',@(h,e)this.keyPress(h,e));
            
            this.imax = axes('parent',this.guifig,'units','normalized','position',[0.15 0.07 0.4 0.87]);
            this.spikeax = axes('parent',this.guifig,'units','normalized','position',[0.59 0.2 0.4 0.7]);
            set(this.guifig,'WindowButtonDownFcn',@(h,e)this.winpressed(h,e,'down'));
            set(this.guifig,'WindowButtonUpFcn',@(h,e)this.winpressed(h,e,'up')) ;
            
            uicontrol('tag','txtimax','style','text','string',sprintf('Working Image\n%s',this.filename),'units','normalized',...
                'position',[0.18 0.95 0.4 0.05], ...
                'BackgroundColor','w');
            this.spike_label = uicontrol('tag','txtroiax','style','text','string','Fluorescence trace','units','normalized',...
                'position',[0.59 0.95 0.4 0.05], ...
                'BackgroundColor','w');
            this.createToolbar(this.guifig);
            %             linkaxes([this.imax this.spikeax]);
            
            this.buttons = gobjects(0);
            this.buttons(end+1) = uicontrol('Parent',this.guifig,'String','Add Spike',...
                'units','normalized',...
                'Position',[0.01 0.65 0.1 0.1], ...
                'Callback',@(h,e)this.addclick(h,e));
            this.buttons(end+1) = uicontrol('Parent',this.guifig,'String','Update spikes for all neurons',...
                'units','normalized',...
                'Position',[0.01 0.52 0.1 0.1],...
                'Callback',@(h,e)this.updateclick(h,e));
            this.buttons(end+1) = uicontrol('Parent',this.guifig,'String','Raster',...
                'units','normalized',...
                'Position',[0.01 0.39 0.1 0.1],...
                'Callback',@(h,e)this.rasterclick(h,e));
            this.buttons(end+1) = uicontrol('Parent',this.guifig,'String','Synchronization & Network',...
                'units','normalized',...
                'Position',[0.01 0.26 0.1 0.1],...
                'Callback',@(h,e)this.scaclick(h,e));
            this.buttons(end+1) = uicontrol('Parent',this.guifig,'String','Make Movie',...
                'units','normalized',...
                'Position',[0.01 0.13 0.1 0.1],...
                'Callback',@(h,e)this.movieclick(h,e));
            this.buttons(end+1) = uicontrol('Parent',this.guifig,'String','Save & Quit',...
                'units','normalized',...
                'Position',[0.01 0.01 0.1 0.1],...
                'Callback',@(h,e)this.saveclick(h,e));
            this.buttons(end+1) = uicontrol('Parent',this.guifig,'String','Summary Statistics',...
                'units','normalized',...
                'Position',[0.01 0.78 0.1 0.1],...
                'Callback',@(h,e)this.summaryclick(h,e));
            this.buttons(end+1) = uicontrol('Style','slider',...
                'Parent',this.guifig,...
                'Max',1,'Min',0,...
                'Value',.85,'Units','normalized',...
                'Position',[.15 .05 .7 .05],...
                'SliderStep',[.005 .002],...
                'String',num2str(5),...
                'Callback',@(src,event)slider_cb(this));
            this.tl = uicontrol('style','text','string',...
                ['Spike Detection threshold = ' num2str(get(this.buttons(end),'Value'))],...
                'units','normalized',...
                'position',[0.15 0.01 0.7 0.03], ...
                'BackgroundColor','g');
            this.checkbox = uicontrol('style','checkbox','Value',1,...
                'string','Overlay spikes','units','normalized','position',[.8 .92 .08 .05],'Callback',@(h,e)this.checkclick(h,e));
        end
        function resizeWindow(this)
            [h,w]=size(this.s.image);
            f = w/h;
            this.figw = this.figh*this.hwar*f;
            
            set(this.guifig,'position',[0 0 this.figw this.figh]);
            movegui(this.guifig,'center');
            set(this.guifig,'visible','on');
            
        end
        function tb=createToolbar(this, fig)
            tb = uitoolbar('parent',fig);
            
            hpt=gobjects(0);
            
            
            %---
            hpt(end+1) = uitoggletool(tb,'CData',localLoadIconCData('tool_zoom_in.png'),...
                'TooltipString','Zoom In',...
                'ClickedCallback',...
                'putdowntext(''zoomin'',gcbo)',...
                'Separator','on');
            hpt(end+1) = uitoggletool(tb,'CData',localLoadIconCData('tool_zoom_out.png'),...
                'TooltipString','Zoom Out',...
                'ClickedCallback',...
                'putdowntext(''zoomout'',gcbo)');
            hpt(end+1) = uitoggletool(tb,'CData',localLoadIconCData('tool_hand.png'),...
                'TooltipString','Pan',...
                'ClickedCallback',...
                'putdowntext(''pan'',gcbo)');
        end
    end
end
function cdata = localLoadIconCData(filename)
% Loads CData from the icon files (PNG, GIF or MAT) in toolbox/matlab/icons.
% filename = info.icon;

% Load cdata from *.gif file
persistent ICONROOT
if isempty(ICONROOT)
    ICONROOT = fullfile(matlabroot,'toolbox','matlab','icons',filesep);
end

if length(filename)>3 && strncmp(filename(end-3:end),'.gif',4)
    [cdata,map] = imread([ICONROOT,filename]);
    % Set all white (1,1,1) colors to be transparent (nan)
    ind = map(:,1)+map(:,2)+map(:,3)==3;
    map(ind) = NaN;
    cdata = ind2rgb(cdata,map);
    
    % Load cdata from *.png file
elseif length(filename)>3 && strncmp(filename(end-3:end),'.png',4)
    [cdata map alpha] = imread([ICONROOT,filename],'Background','none');
    % Converting 16-bit integer colors to MATLAB colorspec
    cdata = double(cdata) / 65535.0;
    % Set all transparent pixels to be transparent (nan)
    cdata(alpha==0) = NaN;
    
    % Load cdata from *.mat file
else
    temp = load([ICONROOT,filename],'cdata');
    cdata = temp.cdata;
end
end
