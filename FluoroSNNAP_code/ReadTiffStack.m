function [FinalImage,time,fps] = ReadTiffStack(FileTif,varargin)
warning off;
% Optional argument to read frames in the range [lo high]
if(nargin>1)
    frange = varargin{1};
end
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);

if(nargin>1)
    frange = varargin{1};
else
    frange = [1 NumberImages];
end
FinalImage=zeros(nImage,mImage,frange(2)-frange(1)+1,'uint16');
time = zeros(frange(2)-frange(1)+1,1);
if(~exist('tifflib','file'))
    
    TifLink = Tiff(FileTif, 'r');
    try
        [~,d] = strtok(InfoImage(frange(1)).DateTime);
        v1 = datevec(d,'HH:MM:SS.FFF');
    end
    cntr = 1;
    for i=frange(1):frange(end)
        TifLink.setDirectory(i);
        FinalImage(:,:,cntr)=TifLink.read();
        try
            if(cntr>1)
                % Read in the time stamp. Assumes that format is YYYYMMDD
                % HH:MM:SS.FFF. Ignore the date part
                [~,d] = strtok(InfoImage(i).DateTime);
                v2 = datevec(d,'HH:MM:SS.FFF');
                time(cntr) = etime(v2,v1);
                
            end
        end
        cntr = cntr+1;
    end
    TifLink.close();
else
    FileID = tifflib('open',FileTif,'r');
    rps = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);
    
    try
        [~,d] = strtok(InfoImage(frange(1)).DateTime);
        v1 = datevec(d,'HH:MM:SS.FFF');
    end
    cntr = 1;
    for i=frange(1):frange(2)
        tifflib('setDirectory',FileID,i);
        % Go through each strip of data.
        rps = min(rps,nImage);
        for r = 1:rps:nImage
            row_inds = r:min(nImage,r+rps-1);
            stripNum = tifflib('computeStrip',FileID,r-1);
            FinalImage(row_inds,:,cntr) = tifflib('readEncodedStrip',FileID,stripNum-1);
        end
        try
            if(i>1)
                % Read in the time stamp. Assumes that format is YYYYMMDD
                % HH:MM:SS.FFF. Ignore the date part
                [~,d] = strtok(InfoImage(i).DateTime);
                v2 = datevec(d,'HH:MM:SS.FFF');
                time(cntr) = etime(v2,v1);
                
            end
        end
        cntr = cntr+1;
    end
    tifflib('close',FileID);
end

% Compute acqusition frame rate
fit = createFit1(frange(1):frange(2),time');
fps = 1./fit.p1;