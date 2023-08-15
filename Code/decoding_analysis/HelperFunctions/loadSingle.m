function out = loadSingle(myPath)    
% keyboard;
    tic
    strPart = strsplit(myPath, '/');
    fprintf('Loading %s...', strPart{end});
    tmp = load(myPath);
    fldnms = fieldnames(tmp);
    out = tmp.(fldnms{1});
    fprintf('finished in %.2f secs\n', toc);
end