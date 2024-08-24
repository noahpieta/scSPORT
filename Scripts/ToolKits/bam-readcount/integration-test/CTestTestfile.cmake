# CMake generated Testfile for 
# Source directory: /home/zhangyu2/projects/toolkits/bam-readcount/repo/integration-test
# Build directory: /home/zhangyu2/projects/toolkits/bam-readcount/integration-test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(RunBamReadcount "sh" "-ec" "PYTHONPATH=':/home/zhangyu2/projects/toolkits/bam-readcount/repo/build-common/python:/mnt/software/unstowable/root_v5.34.30/lib:/mnt/software/unstowable/PBSuite_15.2.20' /home/zhangyu2/projects/toolkits/bam-readcount/repo/integration-test/bam-readcount_test.py /home/zhangyu2/projects/toolkits/bam-readcount/bin/bam-readcount")
set_tests_properties(RunBamReadcount PROPERTIES  LABELS "integration")
