function extractionstr = getExtractionstr(extraction_type)
%GETEXTRACTIONSTR Summary of this function goes here
%   Detailed explanation goes here
switch extraction_type
    case 'reach duration'
        extractionstr = 'ReachDuration';
    case 'uniform length'
        extractionstr = 'UniformLength';
    case 'percentiles'
        extractionstr = 'Percentiles';
end

end

