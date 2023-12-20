function FF_ERN_scoring(wrkdir,savefile,useParallel)

if ~iscell(wrkdir)
    rawfilesloc = dir(fullfile(wrkdir,'*.ept'));
    files = struct2cell(rawfilesloc)';
    files = files(:,1:2);
elseif iscell(wrkdir)
    for ii = 1:size(wrkdir,2)
        if ii == 1
            rawfilesloc = dir(fullfile(wrkdir{ii},'*.ept'));
            files = struct2cell(rawfilesloc)';
            files = files(:,1:2);
        else
            rawfilesloc = dir(fullfile(wrkdir{ii},'*.ept'));
            newfiles = struct2cell(rawfilesloc)';
            newfiles = newfiles(:,1:2);

            files = [files; newfiles];
        end
    end
end

ERPinfo.resp.ern_ROI = {'E6', 'E7', 'E106', 'E129'};
ERPinfo.resp.ern_meanwind = [0 100];
ERPinfo.resp.pe_ROI = {'E62' 'E72' 'E67' 'E77' 'E71' 'E76'};
ERPinfo.resp.pe_meanwind = [200 400];
ERPinfo.resp.events = {'cor' 'err'};


global EPmain EPtictoc

EPtictoc.start=[];
EPtictoc.step=1;
EPtictoc.stop=0;

ep_tictoc('begin');

% Make N by 2 matrix of fieldname + value type
variable_names_types = [...
    ["subjid", "cell"]; ...
    ["task", "cell"]; ...
    ["event", "cell"]; ...
    ["ern", "double"]; ...
    ["pe", "double"]; ...
    ["n_trls", "double"]];
%     ["urtrl", "double"]...



% Make table using fieldnames & value types from above
resp_master = table(...
    'Size',[(size(files,1)*900),size(variable_names_types,1)],...
    'VariableNames', variable_names_types(:,1),...
    'VariableTypes', variable_names_types(:,2));

resp_tbl_ind = 1;
resp_individ = [];


%If user does not indicate whether parallel processing should be used,
% assume that the user does NOT want to process data in parallel
if nargin < 3
    useParallel = 0;
end


%If user indicates that files should be processed in parallel, ensure that
% the parallel processing toolbox is intstalled. If the toolbox is not
% installed, run files serially.
if useParallel == 1
    if exist('parpool.m','file') ~= 2

        dlg = {'Warning: Parallel toolbox not installed';...
            'Data will be run serially rather than in parallel.'};

        for ii = 1:length(dlg)
            fprintf('%s\n',dlg{ii});
        end
        useParallel = 0;
        fprintf('\n\n');
    else
        fprintf('Parallel Toolbox installed\n');
    end
end

fname_out = strcat('FF_ERN_',datestr(now, 'mmddyy'),'.csv');

%Print the date and time that processing began
fprintf(...
    '\n******\nProcessing began on %s at %s \n******\n\n',...
    date, datestr(now, 'HH:MM:SS'));


%start a timer
t1 = tic;

if useParallel == 1

    % Number of files
    nFiles = size(files, 1);

    % Preallocate cell array for temporary filenames
    tempFiles = cell(nFiles, 1);

    % Create a unique temporary directory
    [~,tN] = fileparts(tempname);
    tempDir = fullfile(savefile,tN);
    mkdir(tempDir);

    %use all physical cores (minus 2) for data processing
    nCores = feature('numCores') - 2;

    parpool(nCores);

    %Loop through all subjects in parallel
    parfor ii = 1:size(files, 1)

        fprintf('\n******\nProcessing participant %s\n******\n', files{ii});
        fileloc = fullfile(char(files(ii,2)),char(files(ii,1)));

        % Load ep dataset
        % EPdata = ep_readData('file', fileloc, 'format', 'ep_mat');
        tempVar=load('-mat', fileloc);
        EPdata=tempVar.EPdata;

        resp_individ = resp_scoretrial(EPdata, ERPinfo, EPdata.dataName);

        if ~isempty(resp_individ)

            str = strsplit(EPdata.dataName,'_');
            subjid = str{2};

            % Generate a unique temporary filename
            tempFiles{ii} = fullfile(tempDir, [str{2}, '_', str{3}, '.csv']);

            % Write resp_individ to the temporary file
            writetable(resp_individ, tempFiles{ii});
        end

    end

    %close the pool of workers
    delete(gcp);

    for ii = 1:size(files, 1)
        if ~isempty(tempFiles{ii})
            % Read and append to master file
            data = readtable(tempFiles{ii});
            writetable(data, fullfile(savefile,fname_out), 'WriteMode', 'append');

            % Delete the temporary file
            delete(tempFiles{ii});
        end
    end

    % Remove the temporary directory
    rmdir(tempDir);

else

    for ii = 1:size(files,1)
        %indicate which participant is being processed and the progress
        fprintf('\n******\nProcessing participant %s\n******\n',files{ii});


        fileloc = fullfile(char(files(ii,2)),char(files(ii,1)));

        %load ep dataset
        EPdata = ep_readData('file',fileloc,'format','ep_mat');


        resp_individ = resp_scoretrial(EPdata, ERPinfo, EPdata.dataName);

        if ~isempty(resp_individ)
            tbl_nrow = height(resp_individ);

            resp_master(resp_tbl_ind:(resp_tbl_ind+tbl_nrow-1),:) = resp_individ;
            resp_tbl_ind = resp_tbl_ind + tbl_nrow;

        end



    end


    resp_master = resp_master(1:(resp_tbl_ind-1),:);
    writetable(resp_master,fullfile(savefile,...
        fname_out));
end


%stop timer
t2 = toc(t1);

%incdicate how long it took to process the files
tMinutes = round(t2/60);
disp(['It took ' num2str(tMinutes) ' minutes to process these data.']);


end

function individ = resp_scoretrial(EPdata, ERPinfo, subjid)

str = strsplit(subjid,'_');
subjid = str{2};
task = str{3};

start_ind = 1;

for ii = start_ind:length(EPdata.cellNames)

    event = EPdata.cellNames(ii);

    if strcmp(event,'trial')
        event = EPdata.events{1,ii}.value;
    end

    if (EPdata.analysis.badTrials(ii) == 0)

        ern_mean = zeros(0,length(ERPinfo.resp.ern_ROI));

        for jj = 1:length(ERPinfo.resp.ern_ROI)

            chnind = find(strcmp(EPdata.chanNames,...
                ERPinfo.resp.ern_ROI{jj}));

            ern_mean(jj) = mean(EPdata.data(chnind,...
                knnsearch(EPdata.timeNames,ERPinfo.resp.ern_meanwind(1)):...
                knnsearch(EPdata.timeNames,ERPinfo.resp.ern_meanwind(2)),...
                ii)); %#ok<*FNDSB>

        end

        if length(ern_mean) > 1
            ern_mean = mean(ern_mean);
        end


        pe_mean = zeros(0,length(ERPinfo.resp.pe_ROI));

        for jj = 1:length(ERPinfo.resp.pe_ROI)

            chnind = find(strcmp(EPdata.chanNames,...
                ERPinfo.resp.pe_ROI{jj}));

            pe_mean(jj) = mean(EPdata.data(chnind,...
                knnsearch(EPdata.timeNames,ERPinfo.resp.pe_meanwind(1)):...
                knnsearch(EPdata.timeNames,ERPinfo.resp.pe_meanwind(2)),...
                ii)); %#ok<*FNDSB>

        end

        if length(pe_mean) > 1
            pe_mean = mean(pe_mean);
        end


        if ~isnan(ern_mean) && ~isnan(pe_mean)
            if ~exist('individ','var')

                individ = table;
                individ.subjid = cellstr(subjid);
                individ.task = cellstr(lower(task));
                individ.event = cellstr(event);
                individ.ern = ern_mean;
                individ.pe = pe_mean;

            elseif exist('individ','var')

                row = table;
                row.subjid = cellstr(subjid);
                row.task = cellstr(lower(task));
                row.event = cellstr(event);
                row.ern = ern_mean;
                row.pe = pe_mean;

                individ = vertcat(individ,row); %#ok<AGROW>

            end
        end

    end

end

if ~exist('individ','var')
    individ = [];
end

end
