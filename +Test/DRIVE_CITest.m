fprintf('Hello world!');
exit(0);
%result = runtests('tests', 'IncludeSubPackages', true);
%disp(result.table);
%exit(any([result.Failed]));
