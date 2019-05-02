function writeRPCCoefficientTag(filename,tagValue)
%writeRPCCoefficientTag Write RPCCoefficient TIFF tag
%
%   writeRPCCoefficientTag(FILENAME,TAGVALUE) sets the value of the
%   RPCCoefficientTag to TAGVALUE in the TIFF file specified by FILENAME.
%   TAGVALUE is a 92 element row vector of doubles. The TIFF file should
%   exist and be writable (have "r+" permissions). If the TIFF file does
%   not exist, a new file is created.

% Copyright 2015 The MathWorks, Inc.

validateattributes(filename,{'char'},{'vector','nonempty'}, ...
    mfilename,'FILENAME',1);

lengthOfTag = 92;
validateattributes(tagValue,{'double'},{'real','vector','numel',lengthOfTag}, ...
    mfilename,'TAGVALUE',2);

mexWriteRPCCoefficientTag(filename,tagValue)
