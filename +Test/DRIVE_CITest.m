result = runtests('Test.DRIVE_Convergence', 'Test.DRIVE_Mesh');
disp(result.table);
exit(any([result.Failed]));