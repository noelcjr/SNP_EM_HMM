#!/bin/bash

rm SNP_HMM_MCLUST.out

perl -pi -e 's/XXXX/1kS2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/1kS2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2"

perl -pi -e 's/XXXX/1kS3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/1kS3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3"

perl -pi -e 's/XXXX/1kS4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/1kS4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4"

perl -pi -e 's/XXXX/1kS5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/1kS5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5"

perl -pi -e 's/XXXX/1kS6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/1kS6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6"

####

perl -pi -e 's/XXXX/500S2_500S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S2_500S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S3"

perl -pi -e 's/XXXX/500S2_500S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S2_500S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S4"

perl -pi -e 's/XXXX/500S2_500S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S2_500S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S5"

perl -pi -e 's/XXXX/500S2_500S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S2_500S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S6"

perl -pi -e 's/XXXX/500S3_500S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S3_500S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S2"

perl -pi -e 's/XXXX/500S3_500S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S3_500S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S4"

perl -pi -e 's/XXXX/500S3_500S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S3_500S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S5"

perl -pi -e 's/XXXX/500S3_500S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S3_500S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S6"

perl -pi -e 's/XXXX/500S4_500S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S4_500S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S2"

perl -pi -e 's/XXXX/500S4_500S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S4_500S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S3"

perl -pi -e 's/XXXX/500S4_500S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S4_500S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S5"

perl -pi -e 's/XXXX/500S4_500S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S4_500S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S6"

perl -pi -e 's/XXXX/500S5_500S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S5_500S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S2"

perl -pi -e 's/XXXX/500S5_500S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S5_500S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S3"

perl -pi -e 's/XXXX/500S5_500S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S5_500S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S4"

perl -pi -e 's/XXXX/500S5_500S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S5_500S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S6"

perl -pi -e 's/XXXX/500S6_500S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S6_500S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S2"

perl -pi -e 's/XXXX/500S6_500S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S6_500S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S3"

perl -pi -e 's/XXXX/500S6_500S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S6_500S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S4"

perl -pi -e 's/XXXX/500S6_500S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/500S6_500S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S5"
########

perl -pi -e 's/XXXX/100S2_900S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S2_900S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S3"

perl -pi -e 's/XXXX/100S2_900S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S2_900S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S4"

perl -pi -e 's/XXXX/100S2_900S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S2_900S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S5"

perl -pi -e 's/XXXX/100S2_900S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S2_900S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S6"

perl -pi -e 's/XXXX/100S3_900S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S3_900S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S2"

perl -pi -e 's/XXXX/100S3_900S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S3_900S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S4"

perl -pi -e 's/XXXX/100S3_900S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S3_900S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S5"

perl -pi -e 's/XXXX/100S3_900S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S3_900S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S6"

perl -pi -e 's/XXXX/100S4_900S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S4_900S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S2"

perl -pi -e 's/XXXX/100S4_900S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S4_900S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S3"

perl -pi -e 's/XXXX/100S4_900S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S4_900S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S5"

perl -pi -e 's/XXXX/100S4_900S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S4_900S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S6"

perl -pi -e 's/XXXX/100S5_900S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S5_900S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S2"

perl -pi -e 's/XXXX/100S5_900S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S5_900S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S3"

perl -pi -e 's/XXXX/100S5_900S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S5_900S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S4"

perl -pi -e 's/XXXX/100S5_900S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S5_900S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S6"

perl -pi -e 's/XXXX/100S6_900S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S6_900S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S2"

perl -pi -e 's/XXXX/100S6_900S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S6_900S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S3"

perl -pi -e 's/XXXX/100S6_900S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S6_900S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S4"

perl -pi -e 's/XXXX/100S6_900S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/100S6_900S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S5"

######

perl -pi -e 's/XXXX/900S2_100S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S2_100S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S3"

perl -pi -e 's/XXXX/900S2_100S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S2_100S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S4"

perl -pi -e 's/XXXX/900S2_100S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S2_100S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S5"

perl -pi -e 's/XXXX/900S2_100S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S2_100S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S2S6"

perl -pi -e 's/XXXX/900S3_100S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S3_100S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S2"

perl -pi -e 's/XXXX/900S3_100S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S3_100S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S4"

perl -pi -e 's/XXXX/900S3_100S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S3_100S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S5"

perl -pi -e 's/XXXX/900S3_100S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S3_100S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S3S6"

perl -pi -e 's/XXXX/900S4_100S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S4_100S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S2"

perl -pi -e 's/XXXX/900S4_100S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S4_100S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S3"

perl -pi -e 's/XXXX/900S4_100S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S4_100S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S5"

perl -pi -e 's/XXXX/900S4_100S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S4_100S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S4S6"

perl -pi -e 's/XXXX/900S5_100S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S5_100S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S2"

perl -pi -e 's/XXXX/900S5_100S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S5_100S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S3"

perl -pi -e 's/XXXX/900S5_100S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S5_100S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S4"

perl -pi -e 's/XXXX/900S5_100S6/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S5_100S6/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S5S6"

perl -pi -e 's/XXXX/900S6_100S2/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S6_100S2/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S2"

perl -pi -e 's/XXXX/900S6_100S3/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S6_100S3/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S3"

perl -pi -e 's/XXXX/900S6_100S4/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S6_100S4/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S4"

perl -pi -e 's/XXXX/900S6_100S5/g' SNP_HMM_MCLUST_EM7_4.R
R CMD BATCH SNP_HMM_MCLUST_EM7_4.R
perl -pi -e 's/900S6_100S5/XXXX/' SNP_HMM_MCLUST_EM7_4.R
echo "Ran S6S5"