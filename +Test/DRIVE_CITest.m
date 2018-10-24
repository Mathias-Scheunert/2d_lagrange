result = runtests({'Test.TEST_pick', 'Test.DRIVE_Convergence', 'Test.DRIVE_Mesh'});
disp(result.table);
exit(any([result.Failed]));