#!/bin/bash


#bin/brdf_analyze --take_natural_log=true   --whiten_data=true --scale_covariances=true  --whiten_before_log=false  /media/My\ Book/brdf_data/*binary      cov_TTTF.txt
#bin/brdf_analyze --take_natural_log=true   --whiten_data=false --scale_covariances=true  --whiten_before_log=false  /media/My\ Book/brdf_data/*binary      cov_TFTF.txt
#bin/brdf_analyze --take_natural_log=false  --whiten_data=false --scale_covariances=false  --whiten_before_log=false  /media/My\ Book/brdf_data/*binary      cov_FFFF.txt
#bin/brdf_analyze --take_natural_log=false  --whiten_data=true  --scale_covariances=true   --whiten_before_log=false  /media/My\ Book/brdf_data/*binary      cov_FTTF.txt
#bin/brdf_analyze --take_natural_log=false  --whiten_data=true  --scale_covariances=false  --whiten_before_log=false  /media/My\ Book/brdf_data/*binary      cov_FTFF.txt
#bin/brdf_analyze --take_natural_log=true   --whiten_data=false --scale_covariances=false  --whiten_before_log=false  /media/My\ Book/brdf_data/*binary      cov_TFFF.txt
#bin/brdf_analyze --take_natural_log=true   --whiten_data=true  --scale_covariances=false  --whiten_before_log=false  /media/My\ Book/brdf_data/*binary      cov_TTFF.txt
bin/brdf_analyze --take_natural_log=true   --whiten_data=true  --scale_covariances=false  --whiten_before_log=true   /media/My\ Book/brdf_data/*binary      cov_TTFT.txt
bin/brdf_analyze --take_natural_log=true   --whiten_data=true  --scale_covariances=true   --whiten_before_log=true   /media/My\ Book/brdf_data/*binary      cov_TTTT.txt

