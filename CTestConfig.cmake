## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   ENABLE_TESTING()
##   INCLUDE(CTest)

set(CTEST_PROJECT_NAME "Gismo-stable")
set(CTEST_NIGHTLY_START_TIME "00:00:01 UTC")

set(CTEST_DROP_METHOD "https")
set (# disable certificate (self-signed)
CTEST_CURL_OPTIONS
"CURLOPT_SSL_VERIFYPEER_OFF"
"CURLOPT_SSL_VERIFYHOST_OFF" )
set(CTEST_DROP_SITE "cdash.ricam.oeaw.ac.at")
set(CTEST_DROP_LOCATION "/submit.php?project=Gismo-stable")
set(CTEST_DROP_SITE_CDASH TRUE)
set(CTEST_PROJECT_SUBPROJECTS unittests)
