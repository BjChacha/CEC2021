function x = getIGD(fileName)
    x = [];
    fileID = fopen(fileName, 'r');
    tline = fgets(fileID);
    while ischar(tline)
        line = strtrim(tline);
        line = line(2:size(line, 2) - 1);
        varb = regexp(line, ',', 'split');
        x = [x; str2double(varb)];
        tline = fgets(fileID);
    end
    fclose(fileID);
end