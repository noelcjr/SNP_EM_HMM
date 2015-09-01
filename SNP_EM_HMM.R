# Copyright (c) 2013 Noel Carrascal.
# All rights reserved.

# Redistribution and use in source and binary forms are permitted
# provided that the above copyright notice and this paragraph are
# duplicated in all such forms and that any documentation,
# advertising materials, and other materials related to such
# distribution and use acknowledge that the software was developed
# by Noel Carrascal. The name of Noel Carrascal may not be used to 
# endorse or promote products derived from this software without 
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

# Created by Noel Carrascal, June 6, 2013
# Contact info: noelcarrascal@gmail.com
# Summary: Automated SNP Array Analysis tool.
#
# Known Bugs: Selection of best model is giving wrong states from
#             vitervi's algorithm. Lines 573 to 589. EM, clustering and HMM work.
#             This is key for the working of the program, but not necessary for
#             demostrating machine lerning applicability to this problem since
#             the transitions are selected, but not identified correctly.
library(msm)
library(mclust)
library(RHmm)
#####################################################################
#################    Input File   ########################################################################
dat <- read.table("/home/noelcjr/github/SNP_EM_HMM/Data_Sets_1k/baf_XXXX.dat",header=F);
y <- dat[,1];                
#################    File to output plots ################################################################
png(filename = "/home/noelcjr/github/SNP_EM_HMM/outputs/baf_XXXX.png", width = 800, height = 600, units = "px");
par(mfrow = c(4,1));
#################  Adjustable Parameters ########
n <- 1000;                                      # n = Number of SNP points to be analyzed. Equals the number of lines
incrK <- 0.01;                                  #     in the input file. 
path_to_log_file = "/home/noelcjr/github/SNP_EM_HMM/SNP_EM_HMM.log";
reset_log_file <- TRUE                          # reset_log_file = Results from new tests are appended at the end
iter <- 125;                                    # iter = Some distributions converge faster and an exesive nuber of iterations
                                                #        can give wrong means and standard deviations.
                                                #     of the log file if set to TRUE. This allows to compare new test from previous
threshold <- 0.05;                              #     test results. Setting this to FALSE would delete everything in the log file
meanClusterMethod <- 2;                         #     and continue to do so until reset to TRUE.
x <- 1:n;                                       # The other adjustable parameters have not been tested for different
bins=seq(0,1,by=0.01);                          # values, and should not be changed for now. I future versions these
bins2=seq(0,1,by=0.01);                         # parameters should be tested, and recomended variations described in
bins3=seq(0.01,1,by=0.01);                      # the program's documentation.                                   
#############   Function of a Gaussian  #####################
P <- function(x, mean, variance, pro){                      #
  pro*(exp(-(x-mean)^2/(2*variance)) / sqrt(2*pi*variance));#
}                                                           #
###############  Cluster allocation from hierachical clustering  ##########
GC <- function(t,m,h,cl){
    c1 <- c(0,0,0,0,0,0,0,0);    c2 <- c(0,0,0,0,0,0,0,0);
    if(m[t,1] > 0 && m[t,2] > 0){
       if(h[t] > 0.05){
          c1 <- GC(m[t,1],m,h,cl);   #pass the number of clusters through cl[9]
          c1[9] <- c1[9] + 1;
          c2 <- GC(m[t,2],m,h,c1);
       }else{
       	  c1 <- GC(m[t,1],m,h,cl);
          c2 <- GC(m[t,2],m,h,c1);
       }
       return(c2);
    }else if(m[t,1] < 0 && m[t,2] > 0){
    	if(h[t] > 0.05){
    	   cl[abs(m[t,1])] <- cl[9];
    	   cl[9] <- cl[9] + 1;
    	   c1 <- GC(m[t,2],m,h,cl);
    	   return(c1);
    	}else{
    	   cl[abs(m[t,1])] <- cl[9];
    	   #cl[9] <- cl[9] + 1;
    	   c1 <- GC(m[t,2],m,h,cl);
    	   return(c1);
        }
    }else if(m[t,1] < 0 && m[t,2] < 0){
    	if(h[t] > 0.05){
    	   cl[abs(m[t,1])] <- cl[9];
    	   cl[9] <- cl[9] + 1;
    	   cl[abs(m[t,2])] <- cl[9];
    	   return(cl);
    	}else{
    	   cl[abs(m[t,1])] <- cl[9];
    	   cl[abs(m[t,2])] <- cl[9];
    	   #cl[9] <- cl[9] + 1;
    	   return(cl);
        }
    }else if(m[t,1] > 0 && m[t,2] < 0){
    	if(h[t] > 0.05){  # I haven't seen a merge that meets this criteria. d
    	   cc <- cc + 1;
    	   cl <- GC(m[t,1],m,h,cl);
    	   cl[abs(m[t,2])] <- cc;
    	   return(cl);
    	}else{
    	   cl <- GC(m[t,1],m,h,cl);
    	   cl[abs(m[t,2])] <- cc;
    	   return(cl);
        }
    }
    return(cl);
}
#########  function to renumber monotonically sorted algorithms with GC ##########
numGC <- function(ll){
	srt <- c();
	clstNum <- 1;
	for(i in 1:(length(ll)-1)){
	    if(i == 1){
	       current <- ll[i];
	       srt[i] <- clstNum;
	    }else{
	       if(current == ll[i]){
	       	  srt[i] <- clstNum;
	       }else{
	          current <- ll[i];
	          clstNum <- clstNum + 1;
	          srt[i] <- clstNum;
	       }
	    }
	}
	return(srt)
}
#############   Prior Setup   ####################################################
a <- list();     d <- list();    m <- list();    r <- list();                    #
m[[1]] <- c(0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9);                                    #
r[[1]] <- c(0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04);                            #
a[[1]] <- c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2);                                    #
d[[1]] <- m;     d[[2]] <- r;   d[[3]] <- a;                                     #
a[[1]] <- array(c(P(y,d[[1]][[1]][1],d[[2]][[1]][1],d[[3]][[1]][1]),
                  P(y,d[[1]][[1]][2],d[[2]][[1]][2],d[[3]][[1]][2]),
                  P(y,d[[1]][[1]][3],d[[2]][[1]][3],d[[3]][[1]][3]),
                  P(y,d[[1]][[1]][4],d[[2]][[1]][4],d[[3]][[1]][4]),   
                  P(y,d[[1]][[1]][5],d[[2]][[1]][5],d[[3]][[1]][5]),
                  P(y,d[[1]][[1]][6],d[[2]][[1]][6],d[[3]][[1]][6]),
                  P(y,d[[1]][[1]][7],d[[2]][[1]][7],d[[3]][[1]][7]),     
                  P(y,d[[1]][[1]][8],d[[2]][[1]][8],d[[3]][[1]][8])),dim=c(n,8));#
############   EM Step   #########################################################
ms <- mstep(modelName = "V", data = y, z = a[[1]]);
for(i in 1:iter){
    tempM <- c();  tempS <- c();  tempW <- c();
   	for(j in 1:length(ms$parameters$mean)){
   	    tempM <- c(tempM,ms$parameters$mean[j]);
   	    tempS <- c(tempS,ms$parameters$variance$sigmasq[j]);
   	    tempW <- c(tempW,ms$parameters$pro[j]);
   	}
    d[[1]] <- tempM; 
    d[[2]] <- tempS; 
    d[[3]] <- tempW;
         
    es <- estep(modelName = "V", data = y, parameters = ms$parameters)
    ms <- mstep(modelName = "V", data = y, z = es$z)
}                             
######### Sort the Means, SD and pro in ascending order ################
temp1 <- c();  temp2 <- c();  temp3 <- c();
for(i in 1:(length(d[[1]])-1)){
    for(j in (i+1):length(d[[1]])){
    	if(d[[1]][i] > d[[1]][j]){
    		temp1 <- d[[1]][i];
    		d[[1]][i] <- d[[1]][j];
    		d[[1]][j] <- temp1;
    		temp2 <- d[[2]][i];
    		d[[2]][i] <- d[[2]][j];
    		d[[2]][j] <- temp2;
    		temp3 <- d[[3]][i];
    		d[[3]][i] <- d[[3]][j];
    		d[[3]][j] <- temp3;
    	}
    }	
}
######### Cluster with a threshold and Find Number of Clusters #########
######### This is a poor clustering method that can be improved. #######
if(meanClusterMethod == 1){
   count <- 1;
   cluster <- 1;
   clust <- c();
   clust[count] <- cluster;
   count <- count + 1;
   center <- d[[1]][1];
   for(i in 2:8){
       if((d[[1]][i]-center) < threshold){
	       clust[count] <- cluster;
	       count <- count + 1;
       }else{
           cluster <- cluster + 1;
           clust[count] <- cluster;
           count <- count + 1;
           center <- d[[1]][i];
       }
   }
}else if(meanClusterMethod == 2){
   clust <- c();
   clusters <- c(0,0,0,0,0,0,0,0,1);
   treeNode <- 7;
   finalclusters <- c();
   dst <- dist(d[[1]], method = "euclidean");
   hclst <- hclust(dst, method = "median", members=NULL);
   clust <- GC(treeNode,hclst$merge,hclst$height,clusters);
   clust <- numGC(clust);
}
############# The number of clusters is obtained after grouping the normal distributions. ##################           
numClust <- 1;
for(i in 1:(8-1)){
    if(clust[i] == clust[i+1]){
    }else{numClust <- numClust + 1;}
}
############# Initialize gaussians and calculate possible distributions and probabilities ##################
#############      for all possible cases for the given number of clusters                ##################
#############          This is a first hack, it is possible to improve                    ################## 
gau <-  list();
gau[[1]] <- P(-bins2,d[[1]][1],d[[2]][1],d[[3]][1])+P(+bins2,d[[1]][1],d[[2]][1],d[[3]][1]);
gau[[2]] <- P(+bins2,d[[1]][2],d[[2]][2],d[[3]][2]);
gau[[3]] <- P(+bins2,d[[1]][3],d[[2]][3],d[[3]][3]);
gau[[4]] <- P(+bins2,d[[1]][4],d[[2]][4],d[[3]][4]);
gau[[5]] <- P(+bins2,d[[1]][5],d[[2]][5],d[[3]][5]);
gau[[6]] <- P(+bins2,d[[1]][6],d[[2]][6],d[[3]][6]);
gau[[7]] <- P(+bins2,d[[1]][7],d[[2]][7],d[[3]][7]);
gau[[8]] <- P(-(1-bins2),(1-d[[1]][8]),d[[2]][8],d[[3]][8])+P(+(1-bins2),(1-d[[1]][8]),d[[2]][8],d[[3]][8]);
gau[[9]] <- gau[[1]]+gau[[2]]+gau[[3]]+gau[[4]]+gau[[5]]+gau[[6]]+gau[[7]]+gau[[8]];  
gau[[9]] <- gau[[9]]/sum(gau[[9]]);
dist <- list();
pro <- list();
if(numClust == 1){
   dist[[1]] <- gau[[1]];   dist[[2]] <- gau[[1]];  dist[[3]] <- gau[[1]];
   dist[[4]] <- gau[[1]];   dist[[5]] <- gau[[1]];
   dist[[1]] <- (gau[[1]]+gau[[2]]+gau[[3]]+gau[[4]]+gau[[5]]+gau[[6]]+gau[[7]]+gau[[8]]);
   dist[[1]] <- dist[[1]]/sum(dist[[1]]);   
   dist[[2]] <- dist[[2]]/sum(dist[[2]]);  
   dist[[3]] <- dist[[3]]/sum(dist[[3]]);
   dist[[4]] <- dist[[4]]/sum(dist[[4]]);   
   dist[[5]] <- dist[[5]]/sum(dist[[5]]);
   pro[[1]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
}else if(numClust == 2){
   dist[[1]] <- gau[[1]];   dist[[2]] <- gau[[1]];  dist[[3]] <- gau[[1]];
   dist[[4]] <- gau[[1]];   dist[[5]] <- gau[[1]];
   for(i in 2:length(clust)){
   	   if(clust[i]==1){
   	   	  dist[[1]] <- dist[[1]] + gau[[i]];
   	   	  dist[[3]] <- dist[[3]] + gau[[i]];
   	   }else if(clust[i]==2){
   	      dist[[3]] <- dist[[3]] + gau[[i]];
   	   }
   }
   dist[[1]] <- dist[[1]]/sum(dist[[1]]);   
   dist[[2]] <- dist[[2]]/sum(dist[[2]]);  
   dist[[3]] <- dist[[3]]/sum(dist[[3]]);
   dist[[4]] <- dist[[4]]/sum(dist[[4]]);   
   dist[[5]] <- dist[[5]]/sum(dist[[5]]);
   dist[[6]] <- gau[[8]];
   for(i in (length(clust)-1):1){
   	   if(clust[i]==2){dist[[6]] <- dist[[6]] + gau[[i]];}
   }
   dist[[6]] <- dist[[6]]/sum(dist[[6]]);
   pro[[1]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
   pro[[2]] <- list(dist[[6]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
}else if(numClust == 3){
   dist[[1]] <- gau[[1]];   dist[[2]] <- gau[[1]];  dist[[3]] <- gau[[1]];
   dist[[4]] <- gau[[1]];   dist[[5]] <- gau[[1]];
   for(i in 2:length(clust)){
   	   if(clust[i]==1){
   	   	  dist[[1]] <- dist[[1]] + gau[[i]];
   	   	  dist[[2]] <- dist[[2]] + gau[[i]];
   	   	  dist[[3]] <- dist[[3]] + gau[[i]];
   	   }else if(clust[i]==2){
   	   	  dist[[2]] <- dist[[2]] + gau[[i]];
   	   }else if(clust[i]==3){
   	   	  dist[[2]] <- dist[[2]] + gau[[i]];
   	   	  dist[[3]] <- dist[[3]] + gau[[i]];
   	   }
   }
   dist[[6]] <- gau[[8]];
   for(i in (length(clust)-1):1){
   	   if(clust[i]==3){dist[[6]] <- dist[[6]] + gau[[i]];}
   }
   dist[[1]] <- dist[[1]]/sum(dist[[1]]);   
   dist[[2]] <- dist[[2]]/sum(dist[[2]]);  
   dist[[3]] <- dist[[3]]/sum(dist[[3]]);
   dist[[4]] <- dist[[4]]/sum(dist[[4]]);
   dist[[5]] <- dist[[5]]/sum(dist[[5]]);
   dist[[6]] <- dist[[6]]/sum(dist[[6]]);
   pro[[1]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
   pro[[2]] <- list(dist[[6]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
}else if(numClust == 4){
   dist[[1]] <- gau[[1]];   dist[[2]] <- gau[[1]];  dist[[3]] <- gau[[1]];
   dist[[4]] <- gau[[1]];   dist[[5]] <- gau[[1]];
   for(i in 2:length(clust)){
   	   if(clust[i]==1){
   	   	  dist[[1]] <- dist[[1]] + gau[[i]];
   	   	  dist[[3]] <- dist[[3]] + gau[[i]];
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   }else if(clust[i]==2){
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   }else if(clust[i]==3){
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   }else if(clust[i]==4){
   	      dist[[3]] <- dist[[3]] + gau[[i]];
   	      dist[[4]] <- dist[[4]] + gau[[i]];
   	   }
   }
   dist[[6]] <- gau[[8]];
   for(i in (length(clust)-1):1){
   	   if(clust[i]==numClust){dist[[6]] <- dist[[6]] + gau[[i]];}
   }
   dist[[1]] <- dist[[1]]/sum(dist[[1]]);   
   dist[[2]] <- dist[[2]]/sum(dist[[2]]);  
   dist[[3]] <- dist[[3]]/sum(dist[[3]]);
   dist[[4]] <- dist[[4]]/sum(dist[[4]]);
   dist[[5]] <- dist[[5]]/sum(dist[[5]]);
   dist[[6]] <- dist[[6]]/sum(dist[[6]]);
   pro[[1]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
   pro[[2]] <- list(dist[[6]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
}else if(numClust == 5){
   dist[[1]] <- gau[[1]];   dist[[2]] <- gau[[1]];  dist[[3]] <- gau[[1]];
   dist[[4]] <- gau[[1]];   dist[[5]] <- gau[[1]];
   for(i in 2:length(clust)){
   	   if(clust[i]==1){
   	   	  dist[[1]] <- dist[[1]] + gau[[i]];   
   	   	  dist[[2]] <- dist[[2]] + gau[[i]];
   	   	  dist[[3]] <- dist[[3]] + gau[[i]];
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   }else if(clust[i]==2){
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   }else if(clust[i]==3){
   	   	  dist[[2]] <- dist[[2]] + gau[[i]];
   	   }else if(clust[i]==4){
   	      dist[[4]] <- dist[[4]] + gau[[i]];
   	   }else if(clust[i]==5){
   	      dist[[2]] <- dist[[2]] + gau[[i]];
   	      dist[[3]] <- dist[[3]] + gau[[i]];
   	      dist[[4]] <- dist[[4]] + gau[[i]];
   	   }
   }
   dist[[6]] <- gau[[8]];
   for(i in (length(clust)-1):1){
   	   if(clust[i]==numClust){dist[[6]] <- dist[[6]] + gau[[i]];}
   }
   dist[[1]] <- dist[[1]]/sum(dist[[1]]);   
   dist[[2]] <- dist[[2]]/sum(dist[[2]]);  
   dist[[3]] <- dist[[3]]/sum(dist[[3]]);
   dist[[4]] <- dist[[4]]/sum(dist[[4]]);
   dist[[5]] <- dist[[5]]/sum(dist[[5]]);
   dist[[6]] <- dist[[6]]/sum(dist[[6]]);
   pro[[1]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
   pro[[2]] <- list(dist[[6]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
}else if(numClust == 6){
   dist[[1]] <- gau[[1]];   dist[[2]] <- gau[[1]];  dist[[3]] <- gau[[1]];
   dist[[4]] <- gau[[1]];   dist[[5]] <- gau[[1]];  dist[[6]] <- gau[[1]];
   for(i in 2:length(clust)){
   	   if(clust[i]==1){
   	   	  dist[[1]] <- dist[[1]] + gau[[i]];  
   	   	  dist[[3]] <- dist[[3]] + gau[[i]];
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   	  dist[[6]] <- dist[[6]] + gau[[i]];
   	   }else if(clust[i]==2){
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   }else if(clust[i]==3){
   	   	  dist[[6]] <- dist[[6]] + gau[[i]];
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   }else if(clust[i]==4){
   	      dist[[6]] <- dist[[6]] + gau[[i]];
   	      dist[[5]] <- dist[[5]] + gau[[i]];
   	   }else if(clust[i]==5){
   	      dist[[4]] <- dist[[4]] + gau[[i]];
   	      dist[[5]] <- dist[[5]] + gau[[i]];
   	   }else if(clust[i]==6){ 
   	   	  dist[[3]] <- dist[[3]] + gau[[i]];
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   	  dist[[6]] <- dist[[6]] + gau[[i]];	
   	   }
   }
   dist[[7]] <- gau[[8]];
   for(i in (length(clust)-1):1){
   	   if(clust[i]==numClust){dist[[7]] <- dist[[7]] + gau[[i]];}
   }
   dist[[1]] <- dist[[1]]/sum(dist[[1]]);   
   dist[[2]] <- dist[[2]]/sum(dist[[2]]);  
   dist[[3]] <- dist[[3]]/sum(dist[[3]]);
   dist[[4]] <- dist[[4]]/sum(dist[[4]]);
   dist[[5]] <- dist[[5]]/sum(dist[[5]]);
   dist[[6]] <- dist[[6]]/sum(dist[[6]]);
   dist[[7]] <- dist[[7]]/sum(dist[[7]]);
   pro[[1]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
   pro[[2]] <- list(dist[[7]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
   pro[[3]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[6]],dist[[5]]);
   pro[[4]] <- list(dist[[7]],dist[[2]],dist[[3]],dist[[6]],dist[[5]]);
}else if(numClust == 7){
	print("inside 7");
   dist[[1]] <- gau[[1]];   dist[[2]] <- gau[[1]];  dist[[3]] <- gau[[1]];
   dist[[4]] <- gau[[1]];   dist[[5]] <- gau[[1]];  dist[[6]] <- gau[[1]];
   for(i in 2:length(clust)){
   	   if(clust[i]==1){
   	   	  dist[[1]] <- dist[[1]] + gau[[i]];   
   	   	  dist[[2]] <- dist[[2]] + gau[[i]];
   	   	  dist[[3]] <- dist[[3]] + gau[[i]];
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   	  dist[[6]] <- dist[[6]] + gau[[i]];
   	   }else if(clust[i]==2){
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   }else if(clust[i]==3){
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   	  dist[[6]] <- dist[[6]] + gau[[i]];
   	   }else if(clust[i]==4){
   	      dist[[2]] <- dist[[2]] + gau[[i]];
   	   }else if(clust[i]==5){
   	      dist[[5]] <- dist[[5]] + gau[[i]];
   	      dist[[6]] <- dist[[6]] + gau[[i]];
   	   }else if(clust[i]==6){
   	      dist[[4]] <- dist[[4]] + gau[[i]];
   	      dist[[5]] <- dist[[5]] + gau[[i]];	
   	   }else if(clust[i]==7){
   	      dist[[2]] <- dist[[2]] + gau[[i]];
   	   	  dist[[3]] <- dist[[3]] + gau[[i]];
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   	  dist[[6]] <- dist[[6]] + gau[[i]];
   	   }
   }
   dist[[5]] <- (gau[[1]]+gau[[2]]+gau[[3]]+gau[[4]]+gau[[5]]+gau[[6]]+gau[[7]]+gau[[8]]);
   dist[[7]] <- gau[[8]];
   for(i in (length(clust)-1):1){
   	   if(clust[i]==numClust){dist[[7]] <- dist[[7]] + gau[[i]];}
   }
   dist[[1]] <- dist[[1]]/sum(dist[[1]]);   
   dist[[2]] <- dist[[2]]/sum(dist[[2]]);  
   dist[[3]] <- dist[[3]]/sum(dist[[3]]);
   dist[[4]] <- dist[[4]]/sum(dist[[4]]);
   dist[[5]] <- dist[[5]]/sum(dist[[5]]);
   dist[[6]] <- dist[[6]]/sum(dist[[6]]);
   dist[[7]] <- dist[[7]]/sum(dist[[7]]);
   pro[[1]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
   pro[[2]] <- list(dist[[7]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
   pro[[3]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[6]],dist[[5]]);
}else if(numClust == 8){
   dist[[1]] <- gau[[1]];   dist[[2]] <- gau[[1]];  dist[[3]] <- gau[[1]];
   dist[[4]] <- gau[[1]];   dist[[5]] <- gau[[1]];  dist[[6]] <- gau[[1]];
   dist[[7]] <- gau[[1]];   dist[[8]] <- gau[[1]];  dist[[9]] <- gau[[1]];
   dist[[10]] <- gau[[1]];
   for(i in 2:length(clust)){
   	   if(clust[i]==1){
   	   	  dist[[1]] <- dist[[1]] + gau[[i]];   
   	   	  dist[[2]] <- dist[[2]] + gau[[i]];
   	   	  dist[[3]] <- dist[[3]] + gau[[i]];
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   	  dist[[6]] <- dist[[6]] + gau[[i]];
   	   	  dist[[7]] <- dist[[7]] + gau[[i]];
   	   	  dist[[8]] <- dist[[8]] + gau[[i]];
   	   	  dist[[9]] <- dist[[9]] + gau[[i]];
   	   	  dist[[10]] <- dist[[10]] + gau[[i]];
   	   }else if(clust[i]==2){
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   	  dist[[9]] <- dist[[9]] + gau[[i]];
   	   }else if(clust[i]==3){
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   	  dist[[7]] <- dist[[7]] + gau[[i]];
   	   	  dist[[10]] <- dist[[10]] + gau[[i]];
   	   }else if(clust[i]==4){
   	      dist[[2]] <- dist[[2]] + gau[[i]];
   	      dist[[8]] <- dist[[8]] + gau[[i]];
   	      dist[[9]] <- dist[[9]] + gau[[i]];
   	      dist[[10]] <- dist[[10]] + gau[[i]];
   	   }else if(clust[i]==5){
   	      dist[[6]] <- dist[[6]] + gau[[i]];
   	      dist[[8]] <- dist[[8]] + gau[[i]];
   	      dist[[9]] <- dist[[9]] + gau[[i]];
   	      dist[[10]] <- dist[[10]] + gau[[i]];
   	   }else if(clust[i]==6){
   	      dist[[5]] <- dist[[5]] + gau[[i]];
   	      dist[[7]] <- dist[[7]] + gau[[i]];
   	      dist[[10]] <- dist[[10]] + gau[[i]];	
   	   }else if(clust[i]==7){
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	      dist[[5]] <- dist[[5]] + gau[[i]];
   	      dist[[9]] <- dist[[9]] + gau[[i]];
   	   }else if(clust[i]==8){
   	      dist[[1]] <- dist[[1]] + gau[[i]];   
   	   	  dist[[2]] <- dist[[2]] + gau[[i]];
   	   	  dist[[3]] <- dist[[3]] + gau[[i]];
   	   	  dist[[4]] <- dist[[4]] + gau[[i]];
   	   	  dist[[5]] <- dist[[5]] + gau[[i]];
   	   	  dist[[6]] <- dist[[6]] + gau[[i]];
   	   	  dist[[7]] <- dist[[7]] + gau[[i]];
   	   	  dist[[8]] <- dist[[8]] + gau[[i]];
   	   	  dist[[9]] <- dist[[9]] + gau[[i]];
   	   	  dist[[10]] <- dist[[10]] + gau[[i]];
   	   }
   }
   dist[[5]] <- (gau[[1]]+gau[[2]]+gau[[3]]+gau[[4]]+gau[[5]]+gau[[6]]+gau[[7]]+gau[[8]]);
   dist[[11]] <- gau[[8]];
   for(i in (length(clust)-1):1){
   	   if(clust[i]==numClust){dist[[11]] <- dist[[11]] + gau[[i]];}
   }
   dist[[1]] <- dist[[1]]/sum(dist[[1]]);   
   dist[[2]] <- dist[[2]]/sum(dist[[2]]);  
   dist[[3]] <- dist[[3]]/sum(dist[[3]]);
   dist[[4]] <- dist[[4]]/sum(dist[[4]]);
   dist[[5]] <- dist[[5]]/sum(dist[[5]]);
   dist[[6]] <- dist[[6]]/sum(dist[[6]]);
   dist[[7]] <- dist[[7]]/sum(dist[[7]]);
   dist[[8]] <- dist[[8]]/sum(dist[[8]]);
   dist[[9]] <- dist[[9]]/sum(dist[[9]]);
   dist[[10]] <- dist[[10]]/sum(dist[[10]]);
   pro[[1]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
   pro[[2]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[4]],dist[[9]]);
   pro[[3]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[4]],dist[[10]]);
   pro[[4]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[7]],dist[[5]]);
   pro[[5]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[7]],dist[[9]]);
   pro[[6]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[7]],dist[[10]]);
   pro[[7]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[8]],dist[[5]]);
   pro[[8]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[8]],dist[[9]]);
   pro[[9]] <- list(dist[[1]],dist[[2]],dist[[3]],dist[[8]],dist[[10]]);
   pro[[10]] <- list(dist[[11]],dist[[2]],dist[[3]],dist[[4]],dist[[5]]);
}
###################  Values turned into Strings ######################################################################
bafS1_d <- c()                                        #
for(j in 1:n){                                        #
    ## If full deletion then bafS1_d[j] <- "-0.005";  #
    bafS1_d[j] <- toString(round(y[j],2));            #
}                                                     #
l <- c();         incr <- 0;                          #
for(j in 1:101){                                      #
  l[j] <- toString(incr);                             #
  incr <- incr + incrK;                               #
}                                                     #
############ Initial Probability Generation  ##########################################################################################
numInit <- 1;               inp <- list();                                                                                            #
inp[[1]] <- c(1,0,0,0,0);   inp[[2]] <- c(0,1,0,0,0);  inp[[3]] <- c(0,0,1,0,0);  inp[[4]] <- c(0,0,0,1,0);  inp[[5]] <- c(0,0,0,0,1);#
############   Transition Matrix Generation    ########################################################################################
numTrs <- 1;   trm <- list();   trmN <- list();    indx <- c();                            #
trm[[numTrs]] <- rbind(c(1,0,0,0,0),c(0,1,0,0,0),c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1));  #
trmN[[numTrs]] <- c(0,0);                                                                  #
indx[numTrs] <- numTrs;                                                                    #
numTrs <- numTrs + 1;                                                                      #
for(i in 1:5){                                                                             #
    for(j in 1:5){                                                                         #
    	tempTrs <- rbind(c(1,0,0,0,0),c(0,1,0,0,0),c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1));#
        if(i != j){                                                                        #
           tempTrs[i,j] <- 0.999;                                                          #
           tempTrs[i,i] <- 0.001;                                                          #
           trm[[numTrs]] <- tempTrs;                                                       #
           trmN[[numTrs]] <- c(i,j);
           indx[numTrs] <- numTrs;
           numTrs <- numTrs + 1;
           tempTrs[i,j] <- 0.001;
           tempTrs[i,i] <- 0.999;
           trm[[numTrs]] <- tempTrs;
           trmN[[numTrs]] <- c((i+1),(j+1));
           indx[numTrs] <- numTrs;
           numTrs <- numTrs + 1;    
        }
    }
}
############# Enumerate all combinatios  and pick best viterbi score #########################
hmm <- list();  vit <- list();   track <- list();                                            #
vitLPS <- c();  vitLVS <- c();                                                               #
count <- 1;                                                                                  #
for(i in 1:length(pro)){                                                                     #
    for(j in 1:length(inp)){       # Number Initial Probabilities.                           #
        for(k in 1:length(trm)){    # Number of Transition Matrices.                         #
	    hmm[[count]] <- HMMSet(inp[[j]],trm[[k]],dis="DISCRETE",proba=pro[[i]],labels=l);#
	    vit[[count]] <- viterbi(hmm[[count]],bafS1_d);                                   #
            vitLPS[count] <- vit[[count]]$logProbSeq;                                        #
            vitLVS[count] <- vit[[count]]$logViterbiScore;                                   #
	    track[[count]] <- c(i,j,k,count);                                                #
	    count <- count + 1;                                                              #
        }                                                                                    #
    }                                                                                        #
}                                                                                            #
#####  Best Model bubble sort ################################################################
models <- length(trm)*length(inp)*length(pro);
bestModel <- 1;   offset <- 1;
## vitervi is not well documented in the library's documentation. If my assumptions are right,
## turning NaN in vitLPS to zero and looping for the lowest vitLPS will give the best model.
if(is.nan(vitLPS[bestModel])){vitLPS[bestModel]=0;}
bestModelScore <- vitLPS[bestModel]
for(i in 2:models){
  if(is.nan(vitLPS[i])){vitLPS[i]=0;}  
  if(vitLPS[i] < bestModelScore){
    bestModelScore <- vitLPS[i]; bestModel <- i;
  }
}
## offset is because the first state is actually the empty set,
## no cromosomal information, and no array info.
vit[[bestModel]]$states <- vit[[bestModel]]$states + offset;
################################ Plot ################################################################################
plot(x,y, pch='*', xlim=c(0,n), xlab = "SNPs", ylab = "BAF", main = "Simulated SNP array")                           #
title0 = paste("Histogram of SNP data, clusters = ",toString(numClust)," ",paste(toString(clust)))
hh<-hist(y, breaks=bins, freq=T, col="white", main=title0, xlab="SNPs")
title1 = paste("Normalized histogram and gaussian fit. ",toString(round(d[[1]],3)))
plot(bins3,(hh$density/sum(hh$density)), pch=20, xlab = "SNP Bin", main=title1, ylab = "Probability", xlim=c(0,1), ylim=c(0,max((hh$density/sum(hh$density)))));#
lines(gau[[9]] ~ bins2, col = "red", lwd = 1)                                                                        #
########################   Plot of states   ##########################################################################
title2 = paste("Best Model= (",toString(track[[bestModel]][1]),",",toString(track[[bestModel]][2]),",",toString(track[[bestModel]][3]),") lps=",toString(round(vitLPS[bestModel],3))," S",toString(trmN[track[[bestModel]][3]]),sep="")
plot(x,vit[[bestModel]]$states, type="l", col="red",xlim=c(0,n), ylim=c(0.9,6.1), xlab = "SNP", ylab = "State",main=title2);
texto <- paste(date(),"XXXX, S",toString(trmN[track[[bestModel]][3]]));
write(texto, file = path_to_log_file, ncolumns = 1, append = reset_log_file, sep = " ") 
dev.off()
