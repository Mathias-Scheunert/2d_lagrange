% Script for continuous integration automatic test (used in GIT).
result = runtests({'Test.testLagrange', 'Test.testRaviartThomas'});
disp(result.table);
exit(any([result.Failed]));
