function generateSSFinder(SSFinderFullPath, RHSFileFullPath, SRHSFileFullPath)
%
% Latest revision 14.11.2021 
%
% Authors: M.Yu. Khristichenko (INM RAS)
%          Yu.M. Nechepurenko  (INM RAS)
%          E.V. Sklyarova      (MIPT)
%
% Generation of SSFinder file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: 
%        SSFinderFullPath - path of SSFinder;
%        RHSFileFullPath - path of RHS file;
%        SRHSFileFullPath - path of SRHS file;
%        RHS  - a callback function, specifying right hand side of the
%               system.
%
% OUTPUT:   none
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 3 
        isSimplified = true ; 
    else 
        isSimplified = false;
    end
    
    [SSFinderDir,SSFinderName,~] = fileparts(SSFinderFullPath);
    
    if  ~exist(SSFinderDir,'dir')
         error("This directory doesn't exist: " + SSFinderDir)
    end


    if ~isvarname(SSFinderName)
        error("Incorrect SSFinder name is given. Can't create Matlab function with name: " ...
            + SSFinderName)
    end

    if ~exist(RHSFileFullPath,'file') 
        error("This file doesn't exist: " + RHSFileFullPath)
    end

    if isSimplified
        if ~exist(SRHSFileFullPath,'file') 
            error("This file doesn't exist: " + SRHSFileFullPath)
        end
    end


    %We will call RHS func so we add it to the path 
    [RHSDir,RHSFuncName,~] = fileparts(RHSFileFullPath);
    addpath(RHSDir);

    if isSimplified 
        textFilePath = SRHSFileFullPath; 
    else
        textFilePath = RHSFileFullPath;
    end

    [lines, Fname, coefname, coefs] = formatTextFromFile(textFilePath, isSimplified);

    %If coefs were read from SRHS let's compare them with those from RHS file 
    if isSimplified
          [~, ~, ~, RHScoefs] = formatTextFromFile(RHSFileFullPath, false);

          if  ~isempty(setxor(RHScoefs, coefs))
              RHShas  = strjoin(setdiff(RHScoefs, coefs));
              SRHShas = strjoin(setdiff(coefs, RHScoefs));

              error(strcat("Coefs from SRHS and RHS don't match", ...
                          "\nRHS has but SRHS doesn't: ", RHShas,  ...
                          "\nSRHS has but RHS doesn't: ", SRHShas), class(7));
          end 
    end


    try
        %create struct with the found coef names as fields (values are the ones)
        % and try to call RHS function to get the variable names in the right order 
        
        testCoefStruct = cell2struct(num2cell(ones(1, size(coefs, 1))), coefs, 2); 
        [vars, delayvars] =feval(RHSFuncName, testCoefStruct);
        
        if isempty(vars)
            error("Given function " + RHSFuncName + " returns no vars"); 
        end 
        
         if isempty(delayvars)
            error("Given function " + RHSFuncName + " returns no delayed vars"); 
        end 
          
        vars = string(vars); 
    catch ME
        if (strcmp(ME.identifier,'MATLAB:nonExistentField'))
            unmatchedCoef = regexp(ME.message, '\"\w+\"', 'match'); 
            error(strcat("Can't find coef ", unmatchedCoef, " in file ", ...
                textFilePath, ". Check if it is given as a structure element"));
        end
        rethrow(ME)
    end

    vReg = '[a-zA-Z]\w*'; %varname regex
    
    %  '^([a-zA-Z]\w*)=[^\[]'
    [expressedVars, expressedVarsLines] = ...
        getPatternToks(strcat('^(', vReg, ')=[^\[]'), lines);
    expressedVars = unique(expressedVars); 
    
    if ~isempty(expressedVars)
        additionalVars = setdiff(expressedVars, vars); 
    else
        additionalVars = []; 
    end 
    
    % '^Fname\(\d+\)=(.+);*'
    [F, ~] = getPatternToks(['^', Fname, '\(\d+\)=(.+);*'], lines);
    
    if (isempty(F))
        error("Can't find equations to solve in file: " + textFilePath); 
    end
    
    if isSimplified 
          [SRHSDir,SRHSFuncName,~] = fileparts(SRHSFileFullPath);
          addpath(SRHSDir);
          unexprVars = string(feval(SRHSFuncName, testCoefStruct));
          
          if isempty(unexprVars)
              error("Given function " + SRHSFuncName + " returns no unexpressed vars"); 
          end 
          
    else
        unexprVars = vars; 
    end 

    Fnum          = size(F, 1); 
    unexprVarsNum = size(unexprVars, 1);
    
    if isSimplified
        unexprVarInter = intersect(unexprVars,vars);
        
        if size(unexprVarInter,1) ~= unexprVarsNum
            error("Unexpressed vars vector from SRHS func contains vars that are not from RHS func"); 
        end 
    end 

    if unexprVarsNum > Fnum 
        error("Number of unexpressed vars (" + strjoin(unexprVars) + " " + ...
            " - number is: " + unexprVarsNum + ") is less than the number of equations: " + Fnum )
    end

    newTextEqLines = lines(expressedVarsLines);
    newTextEqLines = [newTextEqLines(:); F(:)]; 

    % if we are parsing RHS file there are the variables with delays we have 
    % to substitute with variables without delays   
    if (~isSimplified)
       delayvars=struct2cell(delayvars);
       for i=1:size(delayvars, 1)
           oldnames = string(delayvars{i});
           for j = 1:size(oldnames,1)
               newTextEqLines = subsNewNames(newTextEqLines, oldnames(j), vars(j)); 
           end
       end
    end 


    % Checking if variables are either expressed or from RHS 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    '\<(?<!coefname\.)(?!coefname )([a-zA-Z]\w*)\>'
    foundVars = unique(getPatternToks(strcat('\<(?<!', coefname ,'\.)(?!', coefname, ')(', vReg ,')\>'), ...
       newTextEqLines));
   
    if (isempty(foundVars))
        error("Can't find any vars in equations in file: " + textFilePath); 
    end
    
    unfoundRHSVars = setdiff(vars, foundVars);
    exsessVars     = setdiff(foundVars, [vars; additionalVars]);

    if ~isempty(unfoundRHSVars) || ~isempty(exsessVars)
        unfoundRHSVars =  strjoin(unfoundRHSVars); 
        exsessVars     = strjoin(exsessVars);
        error(strcat("Found vars don't match with vars from RHS", ...
                          "\nRHS has but SRHS doesn't: ", unfoundRHSVars,  ...
                          "\nSRHS has but RHS doesn't: ", exsessVars), class(7)); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %rename all the vars, coefs, and additional vars to avoid any
    %collisions with the user's variables 
    newvars = getNewNames(vars, 'vardidan');
    newAdVars = getNewNames(additionalVars, 'advardidan'); 
    newcoefs = getNewNames(coefs, 'coefdidan');

    %create coef lines for the new func  
    oldCoefNames = fieldnames(newcoefs);
    newTextCoefLines = cellfun(@(a)generateCoefLines(a, newcoefs.(a)), oldCoefNames, ...
        'UniformOutput', false); 
    newTextCoefLines =  vertcat(newTextCoefLines{:}); 

    %delete 'coefname.' prefix 
    newTextEqLines = cellfun(@(a)regexprep(a, strcat('(\<)',coefname, '\.'), '$1'), ...
        newTextEqLines,  'UniformOutput', false);

    newTextEqLines = renameBasedOnStruct(newTextEqLines, newvars); 
    newTextEqLines = renameBasedOnStruct(newTextEqLines, newAdVars); 
    newTextEqLines = renameBasedOnStruct(newTextEqLines, newcoefs);

    %add 'yi=' and '== 0' for the equation lines in new file and wrap them into ' ' 
    %example 'y4=vardidan4*vardidan8*vardidan7-advardidan1 == 0'
    newF = getVarNames('y', Fnum)'; 
    newTextEqLines(end - Fnum + 1: end) = cellfun(@(a, b)strcat(a, '=', b, ' == 0'), ...
       newF,  newTextEqLines(end - Fnum + 1: end),  'UniformOutput', false); 
    newTextEqLines = cellfun(@(a)strcat('''', a, ''''), newTextEqLines, ...
         'UniformOutput', false);

    newUnexprVarsNames = cellfun(@(a)newvars.(a), unexprVars, 'UniformOutput', false);
    newVarNames = struct2cell(newvars); 
    
    inequalities = cellfun(@(a)strcat(a, ' >=0 '), newVarNames, 'UniformOutput', false); 
    % ''res = NSolve[{y1, y2, y3, y4, vardidan1 >=0, vardidan2 >=0, vardidan3 >=0, vardidan4 >=0, vardidan5 >=0, vardidan6 >=0, vardidan7 >=0, vardidan8 >=0, vardidan9 >=0, vardidan10 >=0}, {vardidan2, vardidan7, vardidan5, vardidan1}, Reals]''
    NSolveCall = ['''res = NSolve[{', strjoin(newF, ', '), ', ', strjoin(inequalities, ', ') ...
        '}, {', strjoin(newUnexprVarsNames, ', '), '}, Reals]'''];
   
    %''math2matlab','{vardidan1, vardidan2, vardidan3, vardidan4, vardidan5, vardidan6, vardidan7, vardidan8, vardidan9, vardidan10}/.res''
    ssRes = ['''math2matlab'',', '''{',  strjoin(newVarNames, ', '), '}/.res'''];

    mainNewText = [newTextCoefLines(:); newTextEqLines(:); NSolveCall; ssRes]; 
    mainNewText = cellfun(@(a)wrapMath(a), mainNewText,  'UniformOutput', false);
    mainNewText{end} = ['ss = ', mainNewText{end}];

    newtext{1,1} =['function ss=', SSFinderName, '(coef)'];
    newtext{2,1}= 'math(''$Version'')'; 
    newtext = [newtext(:); mainNewText(:)];
    newtext{end + 1, 1} = 'math(''quit'')'; 

    filename = fullfile(SSFinderDir, strcat(SSFinderName, ".m"));
    fid = fopen(filename, 'w'); 
    
    if fid == -1 
        error("Can't open file: "  + filename); 
    end 
    fprintf(fid, '%s\n', newtext{:}); 
    fclose(fid);

end 



function [wrappedCommand] = wrapMath(command)
   wrappedCommand = strcat('math(', command, ');'); 
end 

function [coefLines] = generateCoefLines(oldname, newname)
    coefLines{1,1} = strcat('''matlab2math'',''', newname,''', coef.', oldname);
    coefLines{2,1} = strcat('''', newname, '=', newname, '[[1,1]]''');
end


function lines = renameBasedOnStruct(lines, namestruct)
    oldnames = fieldnames(namestruct);
    for i = 1:size(oldnames, 1)  
        lines = subsNewNames(lines, oldnames{i}, namestruct.(oldnames{i})); 
    end
end 

function [toks, linescheme] = getPatternToks(pattern, lines) 
    toks = cellfun(@(a)regexp(a, pattern, 'tokens'), lines,  'UniformOutput', false);
    linescheme = ~cellfun('isempty', toks);
    toks = horzcat(toks{:}); 
    toks = vertcat(toks{:});
end

function [names] = getVarNames(prefix, number) 
    names = num2cell(1:number); 
    names = cellfun(@(c)[prefix, num2str(c)],names,'uni',false);
end 

function [newnames] = getNewNames(oldnames, prefix)
    oldnamessize = size(oldnames, 1); 
    if (oldnamessize == 0)
        newnames = struct([]); 
    else 
        newnames = getVarNames(prefix, size(oldnames, 1)); 
        newnames =  cell2struct(newnames, oldnames, 2);
    end
end 


function newlines = subsNewNames(lines, oldname, newname)
    expression = strcat('(\<)',oldname, '(\>)');
    replace = strcat('$1', newname, '$2');    
    newlines = cellfun(@(a)regexprep(a, expression, replace), lines,  'UniformOutput', false);
end

function [lines, Fname, coefname, coefs] = formatTextFromFile(filePath, isSimplified)
    % read RHS file text - line delimiters are ';' and '\n'
    fid = fopen(filePath);
    if fid == -1 
        error("Can't open file: " + filePath); 
    end
    text = textscan(fid,'%s','delimiter',';\n');
    fclose(fid);

    lines = text{1};
    if isempty(lines)
        error('Error. Empty file is given: ' +  filePath);
    end

    %delete all the commented text
    linesWithComment = cellfun(@(s)contains(s, '%'), lines);
    lines(linesWithComment) = ...
        cellfun(@(s)extractBefore(s, '%'), lines(linesWithComment), ...
        'UniformOutput', false);

    %trim leading and tailing whitespaces 
    lines = cellfun(@strtrim, lines, 'UniformOutput', false);

    %delete '...' if the last line contains it 
    lastIdx = size(lines, 1); 
    if contains(lines{lastIdx}, '...')
        lines{lastIdx} = extractBefore(lines{lastIdx}, '...');   
    end

    % combine all the lines devided by '...'
    for i = lastIdx - 1:-1:1
        if contains(lines{i}, '...')
            lines{i} = extractBefore(lines{i}, '...');
            lines{i} = strcat(lines{i}, lines{i+1});
            lines{i + 1} = ''; 
        end
    end

    % delete whitespaces from all the lines except for the lines where syms vars
    % are declared 
    [~, symsLines] = getPatternToks("^syms\s.*", lines);
    lines(~symsLines) = cellfun(@(a)regexprep(a,'\s',''), lines(~symsLines),...
        'UniformOutput', false); 

    % delete all the empty lines 
    lines(cellfun('isempty',lines)) = [];
    
    if (isempty(lines))
        error('Error. Empty file is given: ' +  filePath);
    end 

    vReg = '[a-zA-Z]\w*'; %varname regex 

    if isSimplified
        % '^functionn\[[a-zA-Z]\w*,([a-zA-Z]\w*)\]=[a-zA-Z]\w*\(([a-zA-Z]\w*)\)$'
        funcReg = strcat('^function\[', vReg,',(', vReg, ')\]=', vReg, '\((', vReg ,')\)$'); 
    else
       %'^function\[[a-zA-Z]\w*,[a-zA-Z]\w*,([a-zA-Z]\w*)\]=[a-zA-Z]\w*\(([a-zA-Z]\w*)\)$'
       funcReg = strcat('^function\[', vReg, ',', vReg,',(', vReg, ')\]=',  ...
           vReg, '\((', vReg ,')\)$');
    end
    
    argnames = regexp(lines{1}, funcReg, 'tokens');
    if isempty(argnames) || size(argnames{1}, 2) ~= 2
            error(strcat("In file ",  filePath, ...
                " code doesn't start with the function declaration of the right format"));
    end

    argnames = argnames{1}; 
    Fname = argnames{1}; 
    coefname = argnames{2}; 

    % '\<coef\.([a-zA-Z]\w*)\>'
    coefs = unique(getPatternToks(strcat('\<', coefname, '\.', '(', vReg,')\>'), lines));
    
    if (isempty(coefs))
        error("Can't find coefs in file: " + filePath); 
    end 
end 


