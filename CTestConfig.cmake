######################################################################
## CTestConfig.txt
## This file is part of the G+Smo library.
##
## Set test timeout (secs) by -DDART_TESTING_TIMEOUT=.. in cmake args
######################################################################

set(CTEST_PROJECT_NAME "Gismo")
set(CTEST_NIGHTLY_START_TIME "00:00:01 CET")
set(CTEST_DROP_METHOD "https")
set(CTEST_CURL_OPTIONS "CURLOPT_SSL_VERIFYPEER_OFF" "CURLOPT_SSL_VERIFYHOST_OFF")
set(CTEST_DROP_SITE "cdash-ci.inria.fr")
set(CTEST_DROP_LOCATION "/submit.php?project=Gismo")
set(CTEST_DROP_SITE_CDASH TRUE)
set(CTEST_OUTPUT_ON_FAILURE TRUE)
set(CTEST_LABELS_FOR_SUBPROJECTS "")
#set(CTEST_USE_LAUNCHERS 1)
