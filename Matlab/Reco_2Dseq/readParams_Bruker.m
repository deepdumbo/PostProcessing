%
% [paramStructure]=backgroundPhaseCorrection(...)
%   [paramStructure]=backgroundPhaseCorrection('dirPath','/Usr/dir/')
%   [paramStructure]=backgroundPhaseCorrection('paramNameList',['PVM_Matrix',..])
%
% Author:   Aurélien TROTIER  (a.trotier@gmail.com)
% Date:     2016-05-26
% Partner:  I worked alone on this program
% Institute: CRMSB (Bordeaux)
%
% Function description:
%   Read parameters from a bruker experiments in acqp, method.
%
%    If optional parameters dirPath is pass, parameters were directly read in this directory.
%    If optional parameters paramNameList is pass, only parameters with the
%    same name were read.
%
% Input:
%
%  Optional:
%    'dirPath': (string) path to the folder of a bruker experiment
%    'paramNameList': (array of string) with the name of each parameters to
%    retrieve from files
%
% Output:
%	paramStructure : structure with the name of each parameters and its value
%
% Algorithm & bibliography:
%
%
% See also :
%   Implementation first done in Goa (Guerbet) and extract to simplify
%   process for research.
%
% To do :
%   read also reco
%

function [paramStructure]=readParams_Bruker(varargin)
%% init
dirPath=[];
paramNameList=[];

% parse optional input parameters
v = 1;
while v < numel(varargin)
    switch varargin{v}
        case 'dirPath'
            assert(v+1<=numel(varargin));
            v = v+1;
            dirPath = varargin{v};
        case 'paramsList'
            assert(v+1<=numel(varargin));
            v = v+1;
            paramNameList=varargin{v};
        otherwise
            fprintf('Unsupported parameter: %s',varargin{v});
    end
    v = v+1;
end

%%
if isempty(dirPath) %choose the directory
    dirPath=uigetdir('*.*','Open rawdata bruker directory');
    cd(dirPath);
    if isequal(dirPath,0)
        error('directory path is not well defined');
    end
end

fileList={'acqp','method'};

for fileIndex = 1 : length(fileList)
    
    fileName=fullfile(dirPath,fileList(fileIndex));
    
    if exist(fileName{1},'file')==2
        % Read parameter file to get a cell array of lines
        lineList = textscan(fileread(fileName{1}), '%s', 'Delimiter', '\n');
        lineList = lineList{1};
        
        paramStructure.dirPath = dirPath;
        
        % For each line
        for lineIndex = 1:length(lineList)
            
            
            
            % Get line
            line = lineList{lineIndex};
            
            % If it is a parameter definition
            if strncmp(line, '##$', 3)
                
                % Parse line to get parameter name and value (or size for arrays)
                index		= find(line == '=', 1);
                paramName	= line(4:(index-1));
                
                paramValue	= line((index+1):end);
                
                % Compare parameter name with those from the list and get its index
                paramIsOK	= strcmp(paramName, paramNameList);
                paramIndex	= find(paramIsOK, 1);
                
                % If parameter name is in the list
                if ~isempty(paramIndex) || isempty(paramNameList)
                    
                    % If parameter is an array (size specification begins with "(")
                    if (paramValue(1) == '(')
                        
                        % If parameter is a string (next line begins with "<")
                        if (lineList{lineIndex+1}(1) == '<') && (paramValue(end) == ')')
                            
                            % Get string
                            lineIndex	= lineIndex + 1;
                            value		= lineList{lineIndex};
                            
                            %                     elseif (lineList{lineIndex+1}(1) == '<')
                            %                         % Else, it is a numeric array
                        else
                            
                            if ((paramValue(end) == ')')==0)
                                paramValue=[paramValue,lineList{lineIndex+1}];
                            end
                            
                            
                            % Get array length
                            arrayLength = str2double(paramValue(2:end-1));
                            
                            if isnan(arrayLength)==0
                                % Get array
                                value = [];
                                while (length(value) < arrayLength)
                                    lineIndex	= lineIndex + 1;
                                    valueList	= textscan(lineList{lineIndex}, '%f');
                                    if isempty(valueList{1})
                                        valueList	= textscan(lineList{lineIndex}, '%s');
                                    end
                                    
                                    value		= [value transpose(valueList{1})];
                                end
                            else %paramValue is a cell
                                
                                idx=1:length(paramValue(1:end));
                                idx=[1,idx(paramValue==','),idx(end)];
                                clear value
                                for i=1:length(idx)-1
                                    value{i}=paramValue(idx(i)+1:idx(i+1)-1);
                                end
                                
                                
                            end
                            
                        end
                        
                        
                        % Else if parameter is a scalar (it begins with a digit or "-")
                    elseif any(paramValue(1) == '-0123456789')
                        
                        % Get scalar
                        value = str2double(paramValue);
                        
                        % Else, it is a string
                    else
                        
                        value = paramValue;
                        
                    end
                    
                    % Update parameter structure
                    paramStructure.(paramName) = value;
                    
                end
            end
        end
    else
        disp([fileName ' does not exist']);
    end
end




end