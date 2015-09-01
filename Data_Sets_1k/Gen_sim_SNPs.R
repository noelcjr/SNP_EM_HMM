library(msm)
##########   Function to simulate SNP data ##########################
baf.simu=function(num,avg,sd){                                      #
  b=c()                                                             #
  for(i in 1:length(avg))                                           #
  {                                                                 #
    print(avg[i])                                                   #
    b=c(b,rtnorm(ceiling(num/length(avg)),avg[i],sd,0,1))           #
  }                                                                 #
  b=sample(b,num)                                                   #
  return(b)                                                         #
}                                                                   #
#####################################################################
S1 <- c()
S2 <- c(0)
S3 <- c(0,0.5,1)
S4 <- c(0,1)
S5 <- c(0,0.4,0.6,1)
S6 <- c(0,0.2,0.4,0.6,0.8,1)
####################################################
mul <- 10;
############   Zero Transitions  ###################
baf0=baf.simu(100*mul,S2,.03); 
write(baf0, file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_1kS2.dat", ncolumns = 1, append = FALSE, sep = " ")  
baf1=baf.simu(100*mul,S3,.03);    
write(baf1, file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_1kS3.dat", ncolumns = 1, append = FALSE, sep = " ")
baf2=baf.simu(100*mul,S4,.03);    
write(baf2, file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_1kS4.dat", ncolumns = 1, append = FALSE, sep = " ")
baf3=baf.simu(100*mul,S5,.03);    
write(baf3, file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_1kS5.dat", ncolumns = 1, append = FALSE, sep = " ")
baf4=baf.simu(100*mul,S6,.03);    
write(baf4, file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_1kS6.dat", ncolumns = 1, append = FALSE, sep = " ")
############   One transition 50/50  #################### 6  ####
bafa=baf.simu(50*mul,S2,.03);
bafb=baf.simu(50*mul,S3,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S2_500S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S2,.03);
bafb=baf.simu(50*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S2_500S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S2,.03);
bafb=baf.simu(50*mul,S5,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S2_500S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S2,.03);
bafb=baf.simu(50*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S2_500S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S3,.03);
bafb=baf.simu(50*mul,S2,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S3_500S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S3,.03);
bafb=baf.simu(50*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S3_500S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S3,.03);
bafb=baf.simu(50*mul,S5,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S3_500S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S3,.03);
bafb=baf.simu(50*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S3_500S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S4,.03);
bafb=baf.simu(50*mul,S2,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S4_500S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S4,.03);
bafb=baf.simu(50*mul,S3,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S4_500S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S4,.03);
bafb=baf.simu(50*mul,S5,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S4_500S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S4,.03);
bafb=baf.simu(50*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S4_500S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S5,.03);
bafb=baf.simu(50*mul,S2,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S5_500S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S5,.03);
bafb=baf.simu(50*mul,S3,.03);
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S5_500S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S5,.03);
bafb=baf.simu(50*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S5_500S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S5,.03);
bafb=baf.simu(50*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S5_500S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S6,.03);
bafb=baf.simu(50*mul,S2,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S6_500S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S6,.03);
bafb=baf.simu(50*mul,S3,.03);
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S6_500S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S6,.03);
bafb=baf.simu(50*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S6_500S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(50*mul,S6,.03);
bafb=baf.simu(50*mul,S5,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_500S6_500S5.dat", ncolumns = 1, append = FALSE, sep = " ")
############   One transition 10/90  #################### 6  ####
bafa=baf.simu(10*mul,S2,.03);
bafb=baf.simu(90*mul,S3,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S2_900S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S2,.03);
bafb=baf.simu(90*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S2_900S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S2,.03);
bafb=baf.simu(90*mul,S5,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S2_900S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S2,.03);
bafb=baf.simu(90*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S2_900S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S3,.03);
bafb=baf.simu(90*mul,S2,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S3_900S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S3,.03);
bafb=baf.simu(90*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S3_900S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S3,.03);
bafb=baf.simu(90*mul,S5,.03);
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S3_900S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S3,.03);
bafb=baf.simu(90*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S3_900S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S4,.03);
bafb=baf.simu(90*mul,S2,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S4_900S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S4,.03);
bafb=baf.simu(90*mul,S3,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S4_900S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S4,.03);
bafb=baf.simu(90*mul,S5,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S4_900S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S4,.03);
bafb=baf.simu(90*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S4_900S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S5,.03);
bafb=baf.simu(90*mul,S2,.03);
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S5_900S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S5,.03);
bafb=baf.simu(90*mul,S3,.03);
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S5_900S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S5,.03);
bafb=baf.simu(90*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S5_900S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S5,.03);
bafb=baf.simu(90*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S5_900S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S6,.03);
bafb=baf.simu(90*mul,S2,.03);
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S6_900S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S6,.03);
bafb=baf.simu(90*mul,S3,.03);
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S6_900S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S6,.03);
bafb=baf.simu(90*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S6_900S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S6,.03);
bafb=baf.simu(90*mul,S5,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S6_900S5.dat", ncolumns = 1, append = FALSE, sep = " ")
############   One transition 90/10  #################### 6  ####
bafa=baf.simu(90*mul,S2,.03);
bafb=baf.simu(10*mul,S3,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S2_100S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S2,.03);
bafb=baf.simu(10*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S2_100S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S2,.03);
bafb=baf.simu(10*mul,S5,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S2_100S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S2,.03);
bafb=baf.simu(10*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S2_100S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S3,.03);
bafb=baf.simu(10*mul,S2,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S3_100S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S3,.03);
bafb=baf.simu(10*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S3_100S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S3,.03);
bafb=baf.simu(10*mul,S5,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S3_100S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S3,.03);
bafb=baf.simu(10*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S3_100S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S4,.03);
bafb=baf.simu(10*mul,S2,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S4_100S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S4,.03);
bafb=baf.simu(10*mul,S3,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S4_100S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S4,.03);
bafb=baf.simu(10*mul,S5,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S4_100S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S4,.03);
bafb=baf.simu(10*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S4_100S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S5,.03);
bafb=baf.simu(10*mul,S2,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S5_100S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S5,.03);
bafb=baf.simu(10*mul,S3,.03);
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S5_100S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S5,.03);
bafb=baf.simu(10*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S5_100S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S5,.03);
bafb=baf.simu(10*mul,S6,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S5_100S6.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S6,.03);
bafb=baf.simu(10*mul,S2,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S6_100S2.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S6,.03);
bafb=baf.simu(10*mul,S3,.03);
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S6_100S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S6,.03);
bafb=baf.simu(10*mul,S4,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S6_100S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(90*mul,S6,.03);
bafb=baf.simu(10*mul,S5,.03);    
write(c(bafa,bafb), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_900S6_100S5.dat", ncolumns = 1, append = FALSE, sep = " ")
############   Two transition 33/34/33  ####################  6 ####
bafa=baf.simu(33*mul,S3,.03);
bafb=baf.simu(34*mul,S4,.03);
bafc=baf.simu(33*mul,S3,.03);   
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_330S3_340S4_330S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(33*mul,S3,.03);
bafb=baf.simu(34*mul,S5,.03);
bafc=baf.simu(33*mul,S3,.03);    
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_330S3_340S5_330S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(33*mul,S4,.03);
bafb=baf.simu(34*mul,S3,.03);
bafc=baf.simu(33*mul,S4,.03);     
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_330S4_340S3_330S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(33*mul,S4,.03);
bafb=baf.simu(34*mul,S5,.03);
bafc=baf.simu(33*mul,S4,.03);    
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_330S4_340S5_330S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(33*mul,S5,.03);
bafb=baf.simu(34*mul,S3,.03);
bafc=baf.simu(33*mul,S5,.03); 
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_330S5_340S3_330S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(33*mul,S5,.03);
bafb=baf.simu(34*mul,S4,.03);
bafc=baf.simu(33*mul,S5,.03);    
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_330S5_340S4_330S5.dat", ncolumns = 1, append = FALSE, sep = " ")
############   Two transition 10/34/56  ####################  6  ###
bafa=baf.simu(10*mul,S3,.03);
bafb=baf.simu(34*mul,S4,.03);
bafc=baf.simu(56*mul,S3,.03);   
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S3_340S4_560S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S3,.03);
bafb=baf.simu(34*mul,S5,.03);
bafc=baf.simu(56*mul,S3,.03);    
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S3_340S5_560S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S4,.03);
bafb=baf.simu(34*mul,S3,.03);
bafc=baf.simu(56*mul,S4,.03);     
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S4_340S3_560S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S4,.03);
bafb=baf.simu(34*mul,S5,.03);
bafc=baf.simu(56*mul,S4,.03);    
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S4_340S5_560S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S5,.03);
bafb=baf.simu(34*mul,S3,.03);
bafc=baf.simu(56*mul,S5,.03); 
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S5_340S3_560S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(10*mul,S5,.03);
bafb=baf.simu(34*mul,S4,.03);
bafc=baf.simu(56*mul,S5,.03);    
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_100S5_340S4_560S5.dat", ncolumns = 1, append = FALSE, sep = " ")
############   Two transition 56/34/10  ####################  6  ###
bafa=baf.simu(56*mul,S3,.03);
bafb=baf.simu(34*mul,S4,.03);
bafc=baf.simu(10*mul,S3,.03);   
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_560S3_340S4_100S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(56*mul,S3,.03);
bafb=baf.simu(34*mul,S5,.03);
bafc=baf.simu(10*mul,S3,.03);    
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_560S3_340S5_100S3.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(56*mul,S4,.03);
bafb=baf.simu(34*mul,S3,.03);
bafc=baf.simu(10*mul,S4,.03);     
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_560S4_340S3_100S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(56*mul,S4,.03);
bafb=baf.simu(34*mul,S5,.03);
bafc=baf.simu(10*mul,S4,.03);    
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_560S4_340S5_100S4.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(56*mul,S5,.03);
bafb=baf.simu(34*mul,S3,.03);
bafc=baf.simu(10*mul,S5,.03); 
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_560S5_340S3_100S5.dat", ncolumns = 1, append = FALSE, sep = " ")
bafa=baf.simu(56*mul,S5,.03);
bafb=baf.simu(34*mul,S4,.03);
bafc=baf.simu(10*mul,S5,.03);
write(c(bafa,bafb,bafc), file = "/home/noelcjr/github/EM_HMM_MCLUST_BAF/Data_Sets_1k/baf_560S5_340S4_100S5.dat", ncolumns = 1, append = FALSE, sep = " ")
