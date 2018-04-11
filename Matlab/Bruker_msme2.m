function [already_reconned_image, numb_slices, numb_evolutions , numb_echo, numb_grille ] = Bruker_msme2(img_dir, scan, slice_numb,evolution,echo_numb,grille_numb, recon)
%
% [already_reconned_image, numb_slices, numb_evolutions, numb_echo] = get_precon_Bruker_image(img_dir, scan, slice_numb, evolution, echo_numb, recon);
%
%  This reads in images which have already been reconstructed in Paravision
%
%
%  img_dir [base directory that you input to the function] should contain:
%       1 [directory for scan #1]
%       2 [directory for scan #2]
%       3 [directory for scan #3]
%       4 [directory for scan #4]
%                acqp  [file]
%                method  [file]
%                 fid  [file]
%            TPNAME..  [file(s)]
%               pdata  [directory for reconned data]
%                            1   [directory for recon #1]
%                                       2dseq  [image data binary file]
%                                        reco  [image info file]
%                                      d3proc  [file]
%                                        meta  [file]
%                                       procs  [file]
%                                         roi  [file]
%
%  Assumptions:
%     slices saved sequentially
%     16 or 32 - bit signed integers used
%     2D images



if     (nargin==6), 
            recon=1;
elseif (nargin==5), 
            recon=1;
            grille_numb=1;
elseif (nargin==4),
            recon=1;
            grille_numb=1;
            echo_numb=1;
elseif (nargin==3),
            recon=1;
            grille_numb=1;
            echo_numb=1;
            evolution = 1;
elseif (nargin==2),
            recon=1;
            grille_numb=1;
            echo_numb=1;
            evolution = 1;
            slice_numb = 1;
elseif (nargin==1), 
            recon=1;
            grille_numb=1;
            echo_numb=1;
            evolution = 1;
            slice_numb = 1;
            scan = 1;
            img_dir = 'C:\C_work\matlab\imgdata\biospec';
            fprintf('Assuming directory %s  \n',img_dir);
end;


if length(img_dir)==0,
            img_dir = 'C:\C_work\matlab\imgdata\biospec';
            fprintf('Assuming directory %s  \n',img_dir);
end;
if img_dir(end)=='\', img_dir = img_dir(1:end-1);  end;


evolution=1;


%  Get wordtype, big/little endian, recosize
%     ASSUMING slices saved sequentially -----------------------------

reco_fullname = sprintf('%s\\%d\\pdata\\%d\\reco',img_dir,scan,recon);
fp = fopen(reco_fullname,'r');
if fp<0, 
    fprintf('\nTrouble finding the image directory you specified (for reco):\n');
    fprintf('Argument img_dir = %s\n\n',img_dir); error('See error message above.'); 
end;

cc1=0;cc2=0;cc3=0;endmarker_reached = 0;
while  (~( (cc1)&(cc2)&(cc3) ))   &  (~endmarker_reached),
   linea1=fgetl(fp);
   if linea1~=-1
      texto=cellstr(linea1);
      
      if ~cc1, 
          cc1=strncmp(texto,'##$RECO_wordtype=_',18); % 16 or 32-bit signed integer assumed
          if cc1,
              RECO_WORDTYPE = str2num(linea1(19:20));
          end;
      end;

      
      if ~cc2, 
          cc2=strncmp(texto,'##$RECO_byte_order=',19); 
          if cc2,
              RECO_BYTE_ORDER = linea1(20);
              if RECO_BYTE_ORDER == 'B', RECO_BYTE_ORDER = 'b'; end;
              if RECO_BYTE_ORDER == 'L', RECO_BYTE_ORDER = 'l'; end;
          end;
      end;

      
      if ~cc3, 
          cc3=strncmp(texto,'##$RECO_size',12); 
          if cc3,
             linea1=fgetl(fp);
             space_location = find(linea1==' '); space_location = space_location(1);
             numb_freq_encodes = str2num(linea1(1: (space_location-1)));
             numb_phase_encodes = str2num(linea1((space_location+1):end));
          end;          
      end;

            
      emrtest = strncmp(texto,'##END=',6); 
      if emrtest, endmarker_reached=1; end;
      
     
   end; % if linea1~=-1
end; % while
fclose(fp);

%...............................
%numbber of TE
numb_echo=1;
acqp_fullname = sprintf('%s\\%d\\method',img_dir,scan);
fp = fopen(acqp_fullname,'r');
if fp<0, 
    fprintf('\nTrouble finding the image directory you specified (for method):\n');
    fprintf('Argument img_dir = %s\n\n',img_dir); error('See error message above.'); 
end;

cc1=0;endmarker_reached = 0;
while  (~( (cc1)))   &  (~endmarker_reached),
   linea1=fgetl(fp);
   if linea1~=-1
      texto=cellstr(linea1);
      
      if ~cc1, 
          cc1=strncmp(texto,'##$PVM_NEchoImages=',19);
          if cc1,
              numb_echo = str2double(linea1(20:end));
              %fprintf('ici %f\n', numb_echo);
          end;
      end;
     emrtest = strncmp(texto,'##END=',6); 
      if emrtest, endmarker_reached=1; end;
      
          
   end; % if linea1~=-1
end; % while

fclose(fp);
%......................



%..............................
% all the TE values, and number of slice
% Defaults:
numb_slices = 1;
numb_grille=1;

acqp_fullname = sprintf('%s\\%d\\acqp',img_dir,scan);
fp = fopen(acqp_fullname,'r');
if fp<0, 
    fprintf('\nTrouble finding the image directory you specified (for reco):\n');
    fprintf('Argument img_dir = %s\n\n',img_dir); error('See error message above.'); 
end;

cc1=0;cc2=0;endmarker_reached = 0;
while  (~( (cc1)&(cc2) &(cc3)))   &  (~endmarker_reached),
   linea1=fgetl(fp);
   if linea1~=-1
      texto=cellstr(linea1);
      
       
      if ~cc1, 
          cc1=strncmp(texto,'##$ACQ_echo_time=( ',19); 
          if cc1,
             linea1=fgetl(fp);
             numb_grille=str2num(linea1(1:end));
             %fprintf('ici %f\n',numb_grille);
          end;          
      end;
      
      
      
      if ~cc2, 
          cc2=strncmp(texto,'##$NSLICES=',11); 
          if cc2,
              numb_slices = str2double(linea1(12:end));
          end;
      end;    
      emrtest = strncmp(texto,'##END=',6); 
      if emrtest, endmarker_reached=1; end;
      
     
   end; % if linea1~=-1
end; % while

fclose(fp);
%......................








%  Figure out directions (in absolute x and y)  of read / phase encodes -----------------------------

% Defaults:
normal_order=1;
IM_SIX = numb_freq_encodes;
IM_SIY = numb_phase_encodes;

d3proc_fullname = sprintf('%s\\%d\\pdata\\%d\\d3proc',img_dir,scan,recon);
fp = fopen(d3proc_fullname,'r');
if fp<0, error('Trouble finding the image directory you specified (for d3proc)'); end;
cc1=0;cc2=0;endmarker_reached = 0;
while  (~( (cc1)&(cc2) ))   &  (~endmarker_reached),
   linea1=fgetl(fp);
   if linea1~=-1
      texto=cellstr(linea1);
      
      if ~cc1, 
          cc1=strncmp(texto,'##$IM_SIX=',10); 
          if cc1,
              IM_SIX = str2num(linea1(11:end));
          end;
      end;

      
      if ~cc2, 
          cc2=strncmp(texto,'##$IM_SIY=',10); 
          if cc2,
              IM_SIY = str2num(linea1(11:end));
          end;
      end;

      
      emrtest = strncmp(texto,'##END=',6); 
      if emrtest, endmarker_reached=1; end;
      
     
   end; % if linea1~=-1
end; % while
fclose(fp);

if cc1 & cc2,
    normal_order = (   (IM_SIX == numb_freq_encodes)  &  (IM_SIY == numb_phase_encodes)     );
end;

numb_evolutions=1;

%fprintf('%d bit    %s endian    matrix %d x %d     %d slices   %d evo %d echo  %d echovalues\n',RECO_WORDTYPE,RECO_BYTE_ORDER,numb_freq_encodes,numb_phase_encodes,numb_slices,numb_evolutions,numb_echo,numb_grille);
% We now have the following parameters:
%     RECO_WORDTYPE    16 or 32
           datatype = sprintf('bit%d',RECO_WORDTYPE);
%     RECO_BYTE_ORDER   b or  l
%     numb_freq_encodes
%     numb_phase_encodes
%     numb_slices
%     numb_evolutions
%     numb_echo


if 1==2,
if (evolution > numb_evolutions) ||  (slice_numb > numb_slices)  |  (echo_numb > numb_echo),
    fprintf('Problem:    requested     total available\n');    
    fprintf('evolution     %d              %d\n',evolution,numb_evolutions);
    fprintf('slice         %d              %d\n',slice_numb,numb_slices);
    fprintf('echo          %d              %d\n',echo_numb,numb_echo);
    fprintf('grille echo          %d              %d\n',grille_numb,numb_grille);
    error('Please correct, then try again');
end;
end;


evolution=1;
twodseq_fullname = sprintf('%s\\%d\\pdata\\%d\\2dseq',img_dir,scan,recon); % The only place where "2dseq" appears
fp = fopen(twodseq_fullname,'r',RECO_BYTE_ORDER);
%fseek(fp,0,'eof');taille=ftell(fp);
%fprintf('taille=%f\n',taille);
if fp<0, error('Trouble finding the image directory you specified (for 2dseq)'); end;
clear your_image; 
%fprintf('numb_echo %f, echo_numb %f\n',numb_echo,echo_numb);
% First, skip over images you don't need:
tmp1 = (numb_freq_encodes * numb_phase_encodes) * (RECO_WORDTYPE / 8 ); % Base image size
%               skip evolutions 
tmp2 = (evolution-1)*tmp1*numb_slices*numb_echo + (slice_numb-1)*tmp1*numb_echo  + (echo_numb-1)*tmp1;
%fprintf('tmp2=%f\n',tmp2);
fseek(fp,tmp2,'bof');



% Finally, read in the data:

if normal_order,  [your_image, count] = fread(fp,[numb_freq_encodes, numb_phase_encodes], datatype);
else,             [your_image, count] = fread(fp,[numb_phase_encodes, numb_freq_encodes], datatype);
end;
already_reconned_image = your_image';
fclose(fp);
