function [Co]=getwidthall(odir,c,lateq,loneq)
%Get water profile for all water masks in a folder.
%str=deblank(['ls ',deblank(odir),'/watermask*tif']);
  str=sprintf('find  %s -name ''watermask*.tif''',deblank(odir));
  [status, cmdout]=system(str);
  Co=[];

  if status~=0||isempty(cmdout)
     return;
  end

  fid2 = fopen([deblank(odir),'/gagewidth.txt'], 'w');

  C = strsplit(strtrim(cmdout));
  %wfile=C{end}
  gagewidth=zeros(length(C),1);
  wstd=0;
  parfor i=1:length(C)
    wfile=C{i};
    fprintf(['\n \n Get width for watermask ',num2str(i),':',wfile,'. '])
    [demdir,name,ext] =fileparts([strtrim(wfile)]); %eg.watermaskWV02_20160827bj80.tif
    satname{i}=name(10:13);
    ymdg(i)=str2num(name(15:28));
    ymd=ymdg(i);
    if strcmp(satname{i},'WV01')
            satflag=1;
    else
            satflag=0;
    end
    satflagg(i)=satflag;


    if ~exist(wfile,'file')
      warning([wfile,' not exit']);
      %continue
    else
    %gagewidth(j)=getwidth(wfile,lateq,loneq);
    data=readGeotiff(wfile);
    gagewidthi=0;
    try
    [gagewidthi,widthp]=getwidth(data,wfile,c,lateq,loneq);
    gagewidth(i)=gagewidthi;
    catch e
     fprintf('\nThere was an error! The message was:\n%s',e.message);
     fprintf(['\nwfile is ',wfile])
     %continue
     gagewidthi=0;
    end
    if gagewidthi>0
%     wstd=0;
%   fprintf(fid2,'%d %12.6f %12.6f %d \n',ymd,gagewidthi,wstd,satflag); %1 WV01; 0 otherwise
    end
    end %if

  end %i
  for i=1:length(C)
    wfile=C{i};
    if gagewidth(i) >0
       wstd=0;
       fprintf(fid2,'%d %12.6f %12.6f %d %s \n',ymdg(i),gagewidth(i),wstd,satflagg(i), wfile); %1 WV01; 0 otherwise
    end
  end
  save gagewidth.mat gagewidth
