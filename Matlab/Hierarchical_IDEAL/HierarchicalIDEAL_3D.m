%%%% This function aims to get a 3D reco slice by slice (avoiding some
%%%% problems with the Tsao Jiang original 3D method handling 3D in one
%%%% shot)

% Description: Fat-water separation by hierarchical decomposition.
%              Version of method described, modified to handle multipeak 
%              and arbitrary TE
% Jiang Y, Tsao J. Fast and Robust Separation of Multiple Chemical Species 
% from Arbitrary Echo Times  with Complete Immunity to Phase Wrapping.
% In: Proceedings of the 20th Annual Meeting of ISMRM, Melbourne, Australia
% 2012

% Updates : 30/04/2018 - HALIOT Kylian, PhD Student
%           Process 3D datasets slice by slice

%% Main function
function [outParams] = HierarchicalIDEAL_3D()

% Add path 
[BASEPATH,~] = fileparts(mfilename('fullpath'));
tmp = BASEPATH; addpath(tmp); fprintf('Adding to path: %s\n',tmp); clear tmp;

% Dialog box
str1 = '2D / 3D';
str2 = '3D (Slice by slice)';
str_about = 'About...';
while 1
  button=questdlg({'Hierarchical IDEAL ...','Which reconstruction approach ?'},'Hierarchical IDEAL',str_about,str1,str2,str2);

  if isequal(button,str1),
    outParams = hIDEAL_Sbs('2D / 3D');
    break;
  elseif isequal(button,str2),
    outParams = hIDEAL_Sbs('3D (Slice by slice)');
    break;
  elseif isequal(button,str_about),
    uiwait(msgbox({'Hierarchical IDEAL is an MRI method for fat-water separation',...
                 'by hierarchical decomposition and direct estimation of phase',...
                 'offsets.',...
                 ' ',...
                 'The method is partly described in:',...
                 ' ',...
                 '   Tsao J, Jiang Y. Hierarchical IDEAL: robust water-fat separation',...
                 '   at high field by multiresolution field map estimation. In: ',...
                 '   Proceedings of the 18th Annual Meeting of ISMRM, Toronto, ON, ',...
                 '   Canada, 2008. p 653'...
                 '   Jiang Y, Tsao J. Fast and Robust Separation of Multiple',...
                 '   Chemical Species from Arbitrary Echo Times  with Complete',...
                 '   Immunity to Phase Wrapping. In: Proceedings of the 20th ',...
                 '   Annual Meeting of ISMRM, Melbourne, Australia 2012'}, 'About...','modal'));
  else % Cancel
     break;
  end
end
clear str1 str2 str_about button;
  
if nargout<1, clear outParams; end

end

%% 3D (Slice by slice)
function [outParams] = hIDEAL_Sbs(mode)

    filename=[];
    algoParams=[];
    
    [BASEPATH,~] = fileparts(mfilename('fullpath'));
    tmp = BASEPATH; addpath(tmp); fprintf('Adding to path: %s\n',tmp); clear tmp;

    if isempty(filename),
      p = pwd;
      while 1,
        datapath = fullfile(BASEPATH,'data');
        if exist(datapath,'dir'), cd(datapath); break; end
        datapath = fullfile(BASEPATH,'test','data');
        if exist(datapath,'dir'), cd(datapath); break; end
        break;
      end; clear datapath;
      [tmpfile, tmppath] = uigetfile('*.mat', 'Pick a dataset file to load');
      if isequal(tmpfile,0) | isequal(tmppath,0), cd(p); clear p; return; end
      filename = fullfile(tmppath,tmpfile); clear tmpfile tmppath;
      cd(p); clear p;
    end
    
    tic;
        fprintf('Loading %s',filename);
        load (filename);
        if ~exist('imDataParams','var') && exist('data','var'),  % 2011.09.22 in case variable is named data
          imDataParams = data; clear data;
        end
        fprintf(' (%.2fs)\n',toc);
        fprintf('Matrix: %d',size(imDataParams.images,1));
        tmpsize = size(imDataParams.images); fprintf(' x %d',tmpsize(2:min(end,3))); 
        if numel(tmpsize)>=4,
          if tmpsize(4)>1, fprintf(' x %d coils',tmpsize(4)); else fprintf(' x 1 coil'); end
        end
        if numel(tmpsize)>=5, fprintf(' x %d TE',tmpsize(5)); end
        if numel(tmpsize)>=6, fprintf(' x %d',tmpsize(6:end)); end
        fprintf('\n'); clear tmpsize;
        
        % Set algoParams
        algoParams.MinFractSizeToDivide = 0.01;
        algoParams.MaxNumDiv = 7;
        algoParams = Check_Params(imDataParams, algoParams);
        
        % Run
        if(strcmpi(mode,'2D / 3D'))
            outParams = fw_i2cm0c_3pluspoint_tsaojiang(imDataParams,algoParams);
        end
        
        if(strcmpi(mode,'3D (Slice by slice)'))
            [~, ~, nZ, ~, ~] = size(imDataParams.images);
            im_tmp.TE                    = imDataParams.TE;
            im_tmp.PrecessionIsClockwise = imDataParams.PrecessionIsClockwise;
            im_tmp.FieldStrength         = imDataParams.FieldStrength;
            
            % Slice by slice Water-Fat separation
            for i = 1:nZ
                clc;
                str_msg = fprintf(' Slice : %d/%d\n', i, nZ);
                im_tmp.images            = imDataParams.images(:,:,i,:,:);
                clear out_tmp;
                out_tmp = fw_i2cm0c_3pluspoint_tsaojiang(im_tmp,algoParams);
                
                % Fill with varying parameters
                outParams.phasemap (:,:,i)       = out_tmp.phasemap;
                outParams.r2starmap(:,:,i)       = out_tmp.r2starmap;
                outParams.fiterror (:,:,i)       = out_tmp.fiterror;
                outParams.species(1).amps(:,:,i) = out_tmp.species(1).amps;
                outParams.species(2).amps(:,:,i) = out_tmp.species(2).amps;
            end
            
            % Fill with constant parameters
            outParams.FieldStrength        = out_tmp.FieldStrength;
            outParams.TE                   = out_tmp.TE;
            outParams.species(1).name      = out_tmp.species(1).name;
            outParams.species(1).frequency = out_tmp.species(1).frequency;
            outParams.species(1).relAmps   = out_tmp.species(1).relAmps;
            outParams.species(2).name      = out_tmp.species(2).name;
            outParams.species(2).frequency = out_tmp.species(2).frequency;
            outParams.species(2).relAmps   = out_tmp.species(2).relAmps;
        end
        
        % Display
        if ~isempty(outParams),
          handles.figure = figure; drawnow;
          handles.outParams = outParams;
          if exist('outParamsMP','var'),
            handles.outParamsMP = outParamsMP;
          end
          handles.DataName = filename;
          guidata(handles.figure,handles);
          fw_showresults(handles);
          clear handles;
        end
        if nargout<2, clear outParamsMP; end
        if nargout<1, clear outParams; end

end

%% Check if there are some missing parameters
function [algoParams] = Check_Params(imDataParams, algoParams)

    if ~isfield(algoParams,'species'), algoParams.species = []; end
    if ~isfield(algoParams,'Verbose'), algoParams.Verbose = []; end
    if ~isfield(algoParams,'AlwaysShowGUI'), algoParams.AlwaysShowGUI = []; end
    if ~isfield(algoParams,'Visualize'), algoParams.Visualize = []; end
    if ~isfield(algoParams,'Visualize_FatMapMultipler'), algoParams.Visualize_FatMapMultipler = []; end
    if ~isfield(algoParams,'MinFractSizeToDivide'), algoParams.MinFractSizeToDivide = []; end
    if ~isfield(algoParams,'MaxNumDiv'), algoParams.MaxNumDiv = []; end
    if ~isfield(algoParams,'AssumeSinglePeakAsWater'), algoParams.AssumeSinglePeakAsWater = []; end
    if ~isfield(algoParams,'SnrToAssumeSinglePeak'), algoParams.SnrToAssumeSinglePeak = []; end
    if ~isfield(algoParams,'CorrectAmpForT2star'), algoParams.CorrectAmpForT2star = []; end
    if ~isfield(algoParams,'MaxR2star'), algoParams.MaxR2star = []; end
    if isempty(algoParams.Verbose), algoParams.Verbose = 1; end
    if isempty(algoParams.AlwaysShowGUI), algoParams.AlwaysShowGUI = 0; end
    if isempty(algoParams.Visualize), algoParams.Visualize = 1; end
    if isempty(algoParams.Visualize_FatMapMultipler), algoParams.Visualize_FatMapMultipler = 1.5; end
    if isempty(algoParams.MinFractSizeToDivide), algoParams.MinFractSizeToDivide = 0.05; end
    if isempty(algoParams.MaxNumDiv), algoParams.MaxNumDiv = 6; end
    if isempty(algoParams.AssumeSinglePeakAsWater), algoParams.AssumeSinglePeakAsWater = 1; end
    if isempty(algoParams.SnrToAssumeSinglePeak), algoParams.SnrToAssumeSinglePeak = 2.5; end
    if isempty(algoParams.CorrectAmpForT2star), algoParams.CorrectAmpForT2star = 1; end
    if isempty(algoParams.MaxR2star), algoParams.MaxR2star = 250; end

    % Parameter validation
    %--------------------------------------------------------
    MinNumIDEALechoes = 3;
    InputMissing_FieldStrength         = 0;
    InputMissing_PrecessionIsClockwise = 0;
    InputMissing_species               = 0;
       
    if ~isfield(imDataParams,'TE'),
        if isfield(imDataParams,'TEs'),
            imDataParams.TE = imDataParams.TEs;
            imDataParams = rmfield(imDataParams,'TEs');
        else
            error('TE field is missing in imDataParams');
        end
    end
    if ~isfield(imDataParams,'images'),
        if isfield(imDataParams,'image'),
            imDataParams.images = imDataParams.image;
            imDataParams = rmfield(imDataParams,'image');
        else
            error('images field is missing in imDataParams');
        end
    end
    if numel(imDataParams.TE)<MinNumIDEALechoes, error('Current implementation handles at least %d echoes only.',MinNumIDEALechoes); end
    if numel(imDataParams.TE)<numel(algoParams.species), error('There are fewer TE (%d) than the number of species (%d).',numel(imDataParams.TE),numel(algoParams.species)); end
    if size(imDataParams.images,5)~=numel(imDataParams.TE), error('The data have a number of echoes (%d) different from the number of TE (%d).', size(imDataParams.images,5),numel(imDataParams.TE)); end
    if ndims(imDataParams.images)>5, error('Current implementation handles image data with 5 dimensions only.'); end
    
    %%  % Check manditory parameters
    %--------------------------------------------------------
    if  ~isfield(imDataParams,'FieldStrength') |  isempty(imDataParams.FieldStrength) | ~isnumeric(imDataParams.FieldStrength) | ~isreal(imDataParams.FieldStrength),
      if isfield(imDataParams,'fieldStrength') & ~isempty(imDataParams.fieldStrength) &  isnumeric(imDataParams.fieldStrength) &  isreal(imDataParams.fieldStrength),
        imDataParams.FieldStrength = imDataParams.fieldStrength;
        imDataParams = rmfield(imDataParams,'fieldStrength');
      else
        InputMissing_FieldStrength = 1;
        imDataParams.FieldStrength = [];
      end
    end  
    if ~isfield(imDataParams,'PrecessionIsClockwise') |  isempty(imDataParams.PrecessionIsClockwise) | ~isnumeric(imDataParams.PrecessionIsClockwise) | ~isreal(imDataParams.PrecessionIsClockwise),
      InputMissing_PrecessionIsClockwise =  1;
      imDataParams.PrecessionIsClockwise = [];
    end
    if ~isfield(algoParams,'species') | isempty(algoParams.species),
      InputMissing_species = 1;
      algoParams.species   = [];
    end

%%  % Get missing input
    %--------------------------------------------------------
    if algoParams.AlwaysShowGUI | InputMissing_species | InputMissing_PrecessionIsClockwise | InputMissing_FieldStrength,
      tmpParams = fw_inputparams(algoParams,imDataParams); drawnow;
      if isempty(tmpParams),  % Cancelled
        clear tmpParams;
        outParams = [];
        return;
      end
      if isfield(tmpParams,'species'              ), algoParams.species               = tmpParams.species; end
      if isfield(tmpParams,'PrecessionIsClockwise'), imDataParams.PrecessionIsClockwise = tmpParams.PrecessionIsClockwise; end
      if isfield(tmpParams,'FieldStrength'        ), imDataParams.FieldStrength         = tmpParams.FieldStrength; end
      clear tmpParams;
    end
end