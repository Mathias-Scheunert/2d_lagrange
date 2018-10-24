result = runtests({'Test.testLagrange', 'Test.testRaviartThomas'});
disp(result.table);
exit(any([result.Failed]));