function [ str_user ] = GetPCname()

current_directory = pwd;
directory_split = strsplit(current_directory,'/');

if(strcmp(directory_split{1,2},'home'))
    
else
    return
end

str_user = sprintf('%s', directory_split{1,3});
rep=['/home/', str_user,'/Tempo/'];

if(exist(rep))
    
else
    return
end

end