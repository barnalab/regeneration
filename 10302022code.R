#==================================================================================================================================================
##DETERMINING ERCC CORRECTION FACTORS
#==================================================================================================================================================
#R version 3.3.2 (2016-10-31) -- "Sincere Pumpkin Patch"
#load library data.table
library(data.table)
#if doesn't load then go to Tools and search for data.table to install
#R version 3.3.2 (2016-10-31) -- "Sincere Pumpkin Patch"
#load your ERCC read counts table (ERCC counts were extracted from the raw read counts data - there are 92 sequences ERCC-00001 through ERCC-000092):
setwd("C:/Users/Olena/Desktop/10252022_code")
ERCC<-read.table("10252022_05192018_set1_161219_set2_170624_expected_counts_raw_data_ERCC_only.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

head(ERCC) 

#This is what my table looks like: 
        
              #Gene   input0h_1 input0h_2 input24h_1 input24h_2 RNP0h_1 RNP0h_2 RNP24h_1 RNP24h_2 light0h_1 light0h_2 light24h_1 light24h_2 heavy0h_1
#1388799 ERCC-00002      2829        13       1381         13     703      11      555       11      2907        42      13671         41      2657
#1388800 ERCC-00003      1362         0        765          4     277       1      260        1       179         1        959          4       175
#1388801 ERCC-00004      1411         0        716          2     329       0      264        0       423         0       2067          5       391
#1388802 ERCC-00009       150         1         72          0      50       2       15        0       201         2        926          4       182
#1388803 ERCC-00012         1         0          0          0       0       0        0        0         0         0          0          0         0
#1388804 ERCC-00013         0         0          1          0       0       0        0        0         0         0          1          0         0

          #heavy0h_2 heavy24h_1 heavy24h_2
#1388799        26      10431         23
#1388800         3        751          3
#1388801         9       1667          0
#1388802         2        774          1
#1388803         0          0          0
#1388804         0          1          0


#I have two time-points 0h and 24h and two library sets for each time-point (labelled _1 and _2). My libraries include input samples (for total RNA) 
#and the pooled free"RNP", "light" and "heavy" polysome fractions at 0h and 24h. This table includes only ERCC read counts extracted from each of these
#libraries. The idea is that you spiked in the SAME amount of ERCC but the read counts here are different. We assume that we added the SAME amount of ERCC so the final reads should be the same.
#Based on this assumption we do a linear transformation of the data sets relative to one set that we pick as normalizer. In my case I normalized all inputs to input0h_1 and all polysome
#libraries to RNP0h_1. The example below is for one libarary each but you would do this for all of them.

nrow(ERCC) 
#[1] 92 
#as expected, I have 92 data points for ERCC in each sample

rownames(ERCC)<-ERCC[,1] #this step is jut to re-name rows based on ERCC name (this is specific to my data table, yours may be different)
ERCC<-ERCC[,-1]
head(ERCC)
#This is what it looks like now:
              #input0h_1 input0h_2 input24h_1 input24h_2 RNP0h_1 RNP0h_2 RNP24h_1 RNP24h_2 light0h_1 light0h_2 light24h_1 light24h_2 heavy0h_1 heavy0h_2 heavy24h_1
#ERCC-00002      2829        13       1381         13     703      11      555       11      2907        42      13671         41      2657        26      10431
#ERCC-00003      1362         0        765          4     277       1      260        1       179         1        959          4       175         3        751
#ERCC-00004      1411         0        716          2     329       0      264        0       423         0       2067          5       391         9       1667
#ERCC-00009       150         1         72          0      50       2       15        0       201         2        926          4       182         2        774
#ERCC-00012         1         0          0          0       0       0        0        0         0         0          0          0         0         0          0
#ERCC-00013         0         0          1          0       0       0        0        0         0         0          1          0         0         0          1
              #heavy24h_2
#ERCC-00002         23
#ERCC-00003          3
#ERCC-00004          0
#ERCC-00009          1
#ERCC-00012          0
#ERCC-00013          0

##Example Normalize one data set by ERCC (my set name is set 1 - #161219); polysome normalized to heavy0, input24h to input0h
#change directory to save where I want it to


#INPUT NORMALIZATION
#input0h_1 vs input24h_1
plot(ERCC$input0h_1, ERCC$input24h_1) #you are plotting the linear relationship between ERCC reads in the two libraries
best_fit <- lm(formula=ERCC$input24h_1 ~ ERCC$input0h_1) #you are fitting a line of best fit for this relationship 
abline(best_fit) #this graph was saved as 10252022_Graph1-ERCC$input0h_1vERCC$input24h_1.jpg
coef(best_fit) #determine the intercept (b) and the slope 
#This is the output for my data:
#(Intercept) ERCC$input0h_1 
#-0.8539127      0.5390024
#this allows you to define the linear relationship (y=mx+b where Intercept = b; and slope m = RCC$input0h_1, x is your "normalizer" library in my case it is 
#input0h_1, and y is the library you are trying to normalize - in my case it is input24h_1
#so y=mx+b is input24h_1=m(input0h_1)+b
#OR:
#input24h_1=0.5390024 (input0h_1)-0.8539127 
#so your normalizer that will make input24h_1 equal to input0h_1 is:
#input24h_1/0.5390024+0.8539127


##Now define the normalizing equation for your polysome libraries (i.e. free/RNP, light and heavy) - I'm normalizng to heavy0h_1

#heavy0h_1 vs RNP0h_1
plot(ERCC$heavy0h_1, ERCC$RNP0h_1) 
best_fit <- lm(formula=ERCC$RNP0h_1 ~ ERCC$heavy0h_1)
abline(best_fit) #this graph was saved as 10252022_Graph2-ERCC$RNP0h_1vERCC$heavy0h_1.jpg
coef(best_fit) 
#(Intercept) ERCC$heavy0h_1 
#10.0502596      0.3748657 
#here y=RNP0h_1, x is heavy0h_1, slope m is 0.3748657 and intercept b is 10.0502596
#so y=mx+b becomes:
#RNP0h_1=0.3748657(heavy0h_1)+10.0502596 
#normalization equation is:
#RNP0h_1/0.3748657-10.0502596 

#Next you go one by one to normalize the rest of the libraries i.e. here it is for normalizing RNP24h_1:
#heavy0h_1 vs RNP24h_1
plot(ERCC$heavy0h_1, ERCC$RNP24h_1) 
best_fit <- lm(formula=ERCC$RNP24h_1 ~ ERCC$heavy0h_1)
abline(best_fit)
coef(best_fit)
#RNP24h_1=0.326368 (heavy0h_1)+7.636469 
#RNP24h_1/0.326368 -7.636469 
#(Intercept) ERCC$heavy0h_1 
#7.636469       0.326368 

#heavy0h_1 vs light0h_1
plot(ERCC$heavy0h_1, ERCC$light0h_1) 
best_fit <- lm(formula=ERCC$light0h_1 ~ ERCC$heavy0h_1)
abline(best_fit)
coef(best_fit)
(Intercept) ERCC$heavy0h_1 
#light0h_1=1.0445083(heavy0h_1)+0.9647055
#light0h_1/1.0445083-0.9647055
#(Intercept) ERCC$heavy0h_1 
#0.9647055      1.0445083 

#heavy0h_1 vs light24h_1
plot(ERCC$heavy0h_1, ERCC$light24h_1) 
best_fit <- lm(formula=ERCC$light24h_1 ~ ERCC$heavy0h_1)
abline(best_fit)
coef(best_fit)
#light24h_1=4.962628(heavy0h_1)+15.289528
#light24h_1/4.962628-15.289528
#(Intercept) ERCC$heavy0h_1 
#15.289528       4.962628 

#heavy0h_1 vs heavy24h_1
plot(ERCC$heavy0h_1, ERCC$heavy24h_1) 
best_fit <- lm(formula=ERCC$heavy24h_1 ~ ERCC$heavy0h_1)
abline(best_fit)
coef(best_fit)
#heavy24h_1=3.856227(heavy0h_1)+16.475864
#heavy24h_1/3.856227-16.475864
#(Intercept) ERCC$heavy0h_1 
#16.475864       3.856227 

##What you would actually do at this point is take these equations and use them to normalize each of the libraries the ERCC reads came from
#i.e. the total library in my case is data$input24h_1 and I would use the same converstion factors so it would look like:
#data$input24h_1 <- ((data$input24h_1/0.5390024)+0.8539127)

#In my case I have two replicates (n=2)
#Each replicate has the following samples: input0h, input24h, RNP0h, RNP24h, light0h, light24h, heavy0h, heavy24h. Replicates were prepared and sequenced
#at different times therefore each replicate is normalized SEPARATELY  to it's own input0h (for input24h) and heavy0h (for RNP, light and heavy) libraries
#Note also if you have a separate "input" library then you can use either ERCC or standard code in limma; if you 
#are using a "total" library reconstructed in silico by adding RNP+light+heavy reads then use the ERCC-normalized reads from RNP+light+heavy to determine the total

#once all of the sets are normalized (see below) you can do analysis in limma as shown below
#for those using this pipeline - please read through the documentation on limma to understand how the calculations work - you will not need all of the equations shown in the example
#i.e. if you just want to look at heavy polysome/total then the equation would be TE=heavy-total

#==================================================================================================================================================
#NORMALIZING DATA
#==================================================================================================================================================
#Note - if you normalize counts for the entire data set (i.e. entire transcriptome set) the correction factors can introduce 'reads' where there were none. Therefore for this analysis
#I'm starting with the subset of raw 'expected counts' file below that has only 473373 lines of data corresponding to all genes that had >0 reads across all 24 libraries 
#(so if any 1 library had a read then it would show up)

library(data.table)
library(edgeR)
setwd("C:/Users/ozhulyn/Desktop/10252022_code")
data<-read.table("10262022_5062018_AS_MATRIX_set_1_161219_set2_170624_expected_counts_matrix_raw_data_REORDERED_FOR_ANALYSIS_rowsum", sep="\t", header=TRUE, stringsAsFactors=FALSE)
nrow(data) #473373 (all columns with '0' reads have been removed)

head(data)
#output looks like this - note c10000_g1 are transcript names in Bryant et al. (2017) transcriptome assembly
                #input0h_1 input0h_2 input24h_1 input24h_2 RNP0h_1 RNP0h_2 RNP24h_1 RNP24h_2 light0h_1 light0h_2 light24h_1 light24h_2 heavy0h_1 heavy0h_2
#c10000_g1           0         0          2          0       0       0        0        0         0         0          0          0         0         1
#c1000000_g1         0         0          0          5       0       0        0        0         0         0          0          0         0         0
#c1000007_g1         1         3          0          4       0       0        0        0         0         0          1          0         0         0
#c1000016_g1         0         0          0          0       0       1        0        0         0         0          0          0         0         0
#c1000016_g2         0         0          0          1       0       0        0        0         0         0          0          0         0         0
#c1000017_g1         0         0          9         11       0       2        0        0         0         0         39         58         0         3
                #heavy24h_1 heavy24h_2
#c10000_g1            0          2
#c1000000_g1          0          0
#c1000007_g1          3          0
#c1000016_g1          0          0
#c1000016_g2          0          0
#c1000017_g1         26         27


#Normalize as follows (normalization is based on the formulas defined in the first section above):

#ERCC normalization sample set #1 161219.txt"
data$input24h_1 <- ((data$input24h_1/0.5390024)+0.8539127)
data$RNP0h_1 <- ((data$RNP0h_1/0.3748657)-10.0502596) 
data$RNP24h_1 <- ((data$RNP24h_1/0.326368)-7.636469) 
data$light0h_1 <- ((data$light0h_1/1.0445083)-0.9647055) 
data$light24h_1 <- ((data$light24h_1/4.962628)-15.289528)
data$heavy24h_1 <- ((data$heavy24h_1/3.856227)-16.475864) 
data$input24h_1M <- (data$input24h_1*(58/88.5)) #correction factor for initial difference in input concentration 


#ERCC normalization sample set #2 170624.txt"
data$input24h_2 <- ((data$input24h_2/1.672407)-0.185723) 
data$RNP0h_2 <- ((data$RNP0h_2/0.5271283)+0.1603392)
data$RNP24h_2 <- ((data$RNP24h_2/0.24806763)-0.01532243) 
data$light0h_2 <- ((data$light0h_2/0.873296370)+0.001509027) 
data$light24h_2 <- ((data$light24h_2/0.94424112)+0.01281649) 
data$heavy24h_2 <- ((data$heavy24h_2/0.8034919)+0.1450031)
data$input24h_2M <- (data$input24h_2*(36.3/30.8)) #correction factor for initial difference in input concentration 

head(data)
#data after normalization
#               input0h_1 input0h_2 input24h_1 input24h_2   RNP0h_1   RNP0h_2  RNP24h_1    RNP24h_2  light0h_1   light0h_2 light24h_1  light24h_2
#c10000_g1           0         0  4.5644713 -0.1857230 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000000_g1         0         0  0.8539127  2.8039799 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000007_g1         1         3  0.8539127  2.2060393 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.088022  0.01281649
#c1000016_g1         0         0  0.8539127 -0.1857230 -10.05026 2.0574106 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000016_g2         0         0  0.8539127  0.4122176 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000017_g1         0         0 17.5514265  6.3916233 -10.05026 3.9544819 -7.636469 -0.01532243 -0.9647055 0.001509027  -7.430789 61.43780506
#             heavy0h_1 heavy0h_2 heavy24h_1 heavy24h_2 input24h_1M input24h_2M
#c10000_g1           0         1 -16.475864  2.6341383   2.9914049  -0.2188878
#c1000000_g1         0         0 -16.475864  0.1450031   0.5596264   3.3046906
#c1000007_g1         0         0 -15.697901  0.1450031   0.5596264   2.5999749
#c1000016_g1         0         0 -16.475864  0.1450031   0.5596264  -0.2188878
#c1000016_g2         0         0 -16.475864  0.1450031   0.5596264   0.4858279
#c1000017_g1         0         3  -9.733522 33.7483288  11.5026298   7.5329846

data<-data[, c(1:2, 17:18, 5:16)] #re-order to keep only input24h_1M and input24h_2M (corrected 24h inputs)
head(data)
#               input0h_1 input0h_2 input24h_1M input24h_2M   RNP0h_1   RNP0h_2  RNP24h_1    RNP24h_2  light0h_1   light0h_2 light24h_1  light24h_2
#c10000_g1           0         0   2.9914049  -0.2188878 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000000_g1         0         0   0.5596264   3.3046906 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000007_g1         1         3   0.5596264   2.5999749 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.088022  0.01281649
#c1000016_g1         0         0   0.5596264  -0.2188878 -10.05026 2.0574106 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000016_g2         0         0   0.5596264   0.4858279 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000017_g1         0         0  11.5026298   7.5329846 -10.05026 3.9544819 -7.636469 -0.01532243 -0.9647055 0.001509027  -7.430789 61.43780506
#               heavy0h_1 heavy0h_2 heavy24h_1 heavy24h_2
#c10000_g1           0         1 -16.475864  2.6341383
#c1000000_g1         0         0 -16.475864  0.1450031
#c1000007_g1         0         0 -15.697901  0.1450031
#c1000016_g1         0         0 -16.475864  0.1450031
#c1000016_g2         0         0 -16.475864  0.1450031
#c1000017_g1         0         3  -9.733522 33.7483288

#create a column that sums the total polysome reads per timepoint per sample (i.e. total0h_1 is RNP0h_1 + light0h_1 + heavy0h_1)
data$total0h_1<-(data$RNP0h_1 + data$light0h_1 + data$heavy0h_1)
data$total0h_2<-(data$RNP0h_2 + data$light0h_2 + data$heavy0h_2)
data$total24h_1<-(data$RNP24h_1 + data$light24h_1 + data$heavy24h_1)
data$total24h_2<-(data$RNP24h_2 + data$light24h_2 + data$heavy24h_2)

head(data)
#               input0h_1 input0h_2 input24h_1M input24h_2M   RNP0h_1   RNP0h_2  RNP24h_1    RNP24h_2  light0h_1   light0h_2 light24h_1  light24h_2
#c10000_g1           0         0   2.9914049  -0.2188878 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000000_g1         0         0   0.5596264   3.3046906 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000007_g1         1         3   0.5596264   2.5999749 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.088022  0.01281649
#c1000016_g1         0         0   0.5596264  -0.2188878 -10.05026 2.0574106 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000016_g2         0         0   0.5596264   0.4858279 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000017_g1         0         0  11.5026298   7.5329846 -10.05026 3.9544819 -7.636469 -0.01532243 -0.9647055 0.001509027  -7.430789 61.43780506
#                 heavy0h_1 heavy0h_2 heavy24h_1 heavy24h_2 total0h_1 total0h_2 total24h_1 total24h_2
#c10000_g1           0         1 -16.475864  2.6341383 -11.01497 1.1618482  -39.40186  2.6316324
#c1000000_g1         0         0 -16.475864  0.1450031 -11.01497 0.1618482  -39.40186  0.1424972
#c1000007_g1         0         0 -15.697901  0.1450031 -11.01497 0.1618482  -38.42239  0.1424972
#c1000016_g1         0         0 -16.475864  0.1450031 -11.01497 2.0589196  -39.40186  0.1424972
#c1000016_g2         0         0 -16.475864  0.1450031 -11.01497 0.1618482  -39.40186  0.1424972
#c1000017_g1         0         3  -9.733522 33.7483288 -11.01497 6.9559909  -24.80078 95.1708114


#re-name columns for consistency
colnames(data)[3] ="input24h_1" #or in data.table setnames(data, 'input24h_1M', 'input24h_1')
colnames(data)[4]="input24h_2" #or in data.table setnames(data, 'input24h_2M', 'input24h_2')
head(data)
#             input0h_1 input0h_2 input24h_1 input24h_2   RNP0h_1   RNP0h_2  RNP24h_1    RNP24h_2  light0h_1   light0h_2 light24h_1  light24h_2
#c10000_g1           0         0  2.9914049 -0.2188878 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000000_g1         0         0  0.5596264  3.3046906 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000007_g1         1         3  0.5596264  2.5999749 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.088022  0.01281649
#c1000016_g1         0         0  0.5596264 -0.2188878 -10.05026 2.0574106 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000016_g2         0         0  0.5596264  0.4858279 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000017_g1         0         0 11.5026298  7.5329846 -10.05026 3.9544819 -7.636469 -0.01532243 -0.9647055 0.001509027  -7.430789 61.43780506
#               heavy0h_1 heavy0h_2 heavy24h_1 heavy24h_2 total0h_1 total0h_2 total24h_1 total24h_2
#c10000_g1           0         1 -16.475864  2.6341383 -11.01497 1.1618482  -39.40186  2.6316324
#c1000000_g1         0         0 -16.475864  0.1450031 -11.01497 0.1618482  -39.40186  0.1424972
#c1000007_g1         0         0 -15.697901  0.1450031 -11.01497 0.1618482  -38.42239  0.1424972
#c1000016_g1         0         0 -16.475864  0.1450031 -11.01497 2.0589196  -39.40186  0.1424972
#c1000016_g2         0         0 -16.475864  0.1450031 -11.01497 0.1618482  -39.40186  0.1424972
#c1000017_g1         0         3  -9.733522 33.7483288 -11.01497 6.9559909  -24.80078 95.1708114


#create a column RNP+Light for calculation of TE (defined as change in heavy/RNP+Light) 
###here have to add other category sums i.e. RNP+Light etc following 
### based on "122017 analysis of three data sets - limma.txt"
data$RL_0h_1<-(data$RNP0h_1+data$light0h_1)
data$RL_24h_1<-(data$RNP24h_1+data$light24h_1)
data$RH_0h_1<-(data$RNP0h_1+data$heavy0h_1) 
data$RH_24h_1<-(data$RNP24h_1+data$heavy24h_1) 
data$LH_0h_1<-(data$light0h_1+data$heavy0h_1) 
data$LH_24h_1<-(data$light24h_1+data$heavy24h_1) 

data$RL_0h_2<-(data$RNP0h_2+data$light0h_2)
data$RL_24h_2<-(data$RNP24h_2+data$light24h_2)
data$RH_0h_2<-(data$RNP0h_2+data$heavy0h_2) 
data$RH_24h_2<-(data$RNP24h_2+data$heavy24h_2) 
data$LH_0h_2<-(data$light0h_2+data$heavy0h_2) 
data$LH_24h_2<-(data$light24h_2+data$heavy24h_2) 


nrow(data) #473373

###re-order so that list matche sorder in Target list
data1<-data[, c(1:20, 21, 27, 22, 28, 23, 29, 24, 30, 25, 31, 26, 32)]
nrow(data1) #473373
head(data1)

#               input0h_1 input0h_2 input24h_1 input24h_2   RNP0h_1   RNP0h_2  RNP24h_1    RNP24h_2  light0h_1   light0h_2 light24h_1  light24h_2
#c10000_g1           0         0  2.9914049 -0.2188878 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000000_g1         0         0  0.5596264  3.3046906 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000007_g1         1         3  0.5596264  2.5999749 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.088022  0.01281649
#c1000016_g1         0         0  0.5596264 -0.2188878 -10.05026 2.0574106 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000016_g2         0         0  0.5596264  0.4858279 -10.05026 0.1603392 -7.636469 -0.01532243 -0.9647055 0.001509027 -15.289528  0.01281649
#c1000017_g1         0         0 11.5026298  7.5329846 -10.05026 3.9544819 -7.636469 -0.01532243 -0.9647055 0.001509027  -7.430789 61.43780506
#               heavy0h_1 heavy0h_2 heavy24h_1 heavy24h_2 total0h_1 total0h_2 total24h_1 total24h_2   RL_0h_1   RL_0h_2  RL_24h_1    RL_24h_2   RH_0h_1
#c10000_g1           0         1 -16.475864  2.6341383 -11.01497 1.1618482  -39.40186  2.6316324 -11.01497 0.1618482 -22.92600 -0.00250594 -10.05026
#c1000000_g1         0         0 -16.475864  0.1450031 -11.01497 0.1618482  -39.40186  0.1424972 -11.01497 0.1618482 -22.92600 -0.00250594 -10.05026
#c1000007_g1         0         0 -15.697901  0.1450031 -11.01497 0.1618482  -38.42239  0.1424972 -11.01497 0.1618482 -22.72449 -0.00250594 -10.05026
#c1000016_g1         0         0 -16.475864  0.1450031 -11.01497 2.0589196  -39.40186  0.1424972 -11.01497 2.0589196 -22.92600 -0.00250594 -10.05026
#c1000016_g2         0         0 -16.475864  0.1450031 -11.01497 0.1618482  -39.40186  0.1424972 -11.01497 0.1618482 -22.92600 -0.00250594 -10.05026
#c1000017_g1         0         3  -9.733522 33.7483288 -11.01497 6.9559909  -24.80078 95.1708114 -11.01497 3.9559909 -15.06726 61.42248263 -10.05026
#           RH_0h_2  RH_24h_1   RH_24h_2    LH_0h_1     LH_0h_2  LH_24h_1   LH_24h_2
#c10000_g1   1.1603392 -24.11233  2.6188159 -0.9647055 1.001509027 -31.76539  2.6469548
#c1000000_g1 0.1603392 -24.11233  0.1296807 -0.9647055 0.001509027 -31.76539  0.1578196
#c1000007_g1 0.1603392 -23.33437  0.1296807 -0.9647055 0.001509027 -30.78592  0.1578196
#c1000016_g1 2.0574106 -24.11233  0.1296807 -0.9647055 0.001509027 -31.76539  0.1578196
#c1000016_g2 0.1603392 -24.11233  0.1296807 -0.9647055 0.001509027 -31.76539  0.1578196
#c1000017_g1 6.9544819 -17.36999 33.7330064 -0.9647055 3.001509027 -17.16431 95.1861338

str(data1) #verify heading order because must match targets file below
data<-data1
head(data)
rm(data1)


#filter the data based on counts per million in input
data1<-subset(data, cpm(data$input0h_1)>1 & cpm(data$input0h_2)>1 & cpm(data$input24h_1)>1 & cpm(data$input24h_2)>1)
nrow(data1) #21140
head(data1)

#filter the data based on each sample set having greater than 0 reads per library and at least 10 reads in total in RNP+light+heavy polysome fraction for each set and time-point.
data2<-subset(data1, RNP0h_1>0 & light0h_1>0 & heavy0h_1>0 & RNP0h_2>0 & light0h_2>0 & heavy0h_2>0 & total0h_1>10 & total0h_2>10 & RNP24h_1>0 & light24h_1>0 & heavy24h_1>0 & RNP24h_2>0 & light24h_2>0 & heavy24h_2>0 & total24h_1>10 & total24h_2>10) 
nrow(data2)#8139

setwd("C:/Users/Olena/Desktop/10252022_code")
write.table(data2, file="05192018_all_libraries_input_cpm_gr1_allgr0_total_gr10_subset_for_TEanalysis.txt", sep="\t")
data<-data2


#==================================================================================================================================================
#NORMALIZING DATA
#==================================================================================================================================================

setwd("C:/Users/Olena/Box Sync/Work_Original/R/R_Yoga/04292018 Final Polysome Seq Data Analysis/05152018 data analysis DE")
targets <- readTargets("09022017_Targets_Original.txt", sep="\t")  
targets

#    Batch      Sample     Type
#1  161219  input0h_1   input0h
#2  170624  input0h_2   input0h
#3  161219  input24h_1 input24h
#4  170624  input24h_2 input24h
#5  161219  RNP0h_1       RNP0h
#6  170624  RNP0h_2       RNP0h
#7  161219  RNP24h_1     RNP24h
#8  170624  RNP24h_2     RNP24h
#9  161219  light0h_1   light0h
#10 170624  light0h_2   light0h
#11 161219  light24h_1 light24h
#12 170624  light24h_2 light24h
#13 161219  heavy0h_1   heavy0h
#14 170624  heavy0h_2   heavy0h
#15 161219  heavy24h_1 heavy24h
#16 170624  heavy24h_2 heavy24h
#17 161219  total0h_1   total0h
#18 170624  total0h_2   total0h
#19 161219  total24h_1 total24h
#20 170624  total24h_2 total24h
#21 161219  RL_0h_1       RL_0h
#22 170624  RL_0h_2       RL_0h
#23 161219  RL_24h_1     RL_24h
#24 170624  RL_24h_2     RL_24h
#25 161219  RH_0h_1       RH_0h
#26 170624  RH_0h_2       RH_0h
#27 161219  RH_24h_1     RH_24h
#28 170624  RH_24h_2     RH_24h
#29 161219  LH_0h_1       LH_0h
#30 170624  LH_0h_2       LH_0h
#31 161219  LH_24h_1     LH_24h
#32 170624  LH_24h_2     LH_24h


x<-data #note this is the data file from section above
group <- factor(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16))
dge<-DGEList(counts=x, group=group)
dge <- calcNormFactors(dge, method="none") 
batch <- factor(targets$Batch)
design<-model.matrix(~0+group+batch) 
factor_names<-factor(targets$Sample, levels=c("input0h", "input24h", "RNP0h", "RNP24h", "light0h", "light24h", "heavy0h", "heavy24h", "total0h", "total24h", "RL_0h", "RL_24h", "RH_0h", "RH_24h", "LH_0h", "LH_24h"))

colnames(design)<-c("input0h", "input24h", "RNP0h", "RNP24h", "light0h", "light24h", "heavy0h", "heavy24h", "total0h", "total24h", "RL_0h", "RL_24h", "RH_0h", "RH_24h", "LH_0h", "LH_24h","batch")

design


contrast.matrix <- makeContrasts(DE_input=input24h-input0h, TEhvl_change_heavy=((heavy24h-heavy0h)-(RL_24h-RL_0h)), TEhvl_change_light=((light24h-light0h)-(RH_24h-RH_0h)), TEhvl_change_RNP=((RNP24h-RNP0h)-(LH_24h-LH_0h)), input0h, input24h, TEhvl_RNP0=RNP0h-LH_0h, TEhvl_RNP24=RNP24h-LH_24h, TEhvl_light0=light0h-RH_0h, TEhvl_light24=light24h-RH_24h, TEhvl_heavy0=heavy0h-RL_0h, TEhvl_heavy24=heavy24h-RL_24h, DE_total=total24h-total0h, TEin_RNP0=RNP0h-input0h, TEin_RNP24=RNP24h-input24h, TEin_light0=light0h-input0h, TEin_light24=light24h-input24h, TEin_heavy0=heavy0h-input0h, TEin_heavy24=heavy24h-input24h, TEin_change_RNP=((RNP24h-RNP0h)-(input24h-input0h)), TEin_change_light=((light24h-light0h)-(input24h-input0h)), TEin_change_heavy=((heavy24h-heavy0h)-(input24h-input0h)), TEtot_RNP0=RNP0h-total0h, TEtot_RNP24=RNP24h-total24h, TEtot_light0=light0h-total0h, TEtot_light24=light24h-total24h, TEtot_heavy0=heavy0h-total0h, TEtot_heavy24=heavy24h-total24h, TEtot_change_RNP=((RNP24h-RNP0h)-(total24h-total0h)), TEtot_change_light=((light24h-light0h)-(total24h-total0h)), TEtot_change_heavy=((heavy24h-heavy0h)-(total24h-total0h)), levels=design)
contrast.matrix


v <- voom(dge, design, plot=TRUE) 
fit <- lmFit(v, design) 
fit2 <-contrasts.fit(fit, contrast.matrix) 
fit3 <- eBayes(fit2)
plotSA(fit3)
results<-decideTests(fit3, method="separate", lfc=1) 

setwd("C:/Users/Olena/Desktop/10252022_code")


hist(fit2$coefficients[,"DE_input"],50)

hist(fit2$coefficients[,"TEhvl_change_heavy"],50)

hist(fit2$coefficients[,"TEhvl_change_light"],50)

hist(fit2$coefficients[,"TEhvl_change_RNP"],50)



head(fit2$coefficients)

summary(results)


hist(fit3$p.value[, "DE_input"], 50)

hist(fit3$p.value[, "TEhvl_change_heavy"], 50)

hist(fit3$p.value[, "TEhvl_change_RNP"], 50)

hist(fit3$p.value[, "TEhvl_change_light"], 50)



setwd("C:/Users/Olena/Box Sync/Work_Original/R/R_Yoga/04292018 Final Polysome Seq Data Analysis/05192018 final DE and TE")
write.table(fit2$coefficients, file="05202018_limma_BOTH_161219_and_170624_TE_8139genes.txt", sep="\t")

DE_input <- topTable(fit3, coef="DE_input", adjust="BH", n=Inf)

TEhvl_change_heavy <- topTable(fit3, coef="TEhvl_change_heavy", adjust="BH", n=Inf)

write.table(DE_input, file="05202018_limmafit3_DE_input.txt", sep="\t", row.names=TRUE, col.names=TRUE)  

write.table(TEhvl_change_heavy,file="05202018_limmafit3_TEhvl_change_heavy.txt", sep="\t", row.names=TRUE, col.names=TRUE)  

DE_input <- topTable(fit3, coef="DE_input", adjust="BH", n=Inf)
#data output of this includes the following: logFC  AveExpr        t       P.Value     adj.P.Val        B

plot(DE_input$AveExpr, DE_input$logFC, main="MA plot", pch=20, axes=TRUE, cex=0.8, xlab="AveExpr", ylab="logFC", col=ifelse(abs(DE_input$logFC)>1 & DE_input$adj.P.Val<0.05, "red", "black"))
abline(h=c(-1,1), lty=2, col="grey")

TEhvl_change_heavy <- topTable(fit3, coef="TEhvl_change_heavy", adjust="BH", n=Inf)
plot(TEhvl_change_heavy$AveExpr, TEhvl_change_heavy$logFC, main="MA plot", pch=20, axes=TRUE, cex=0.8, xlab="AveExpr", ylab="logFC", col=ifelse(abs(TEhvl_change_heavy$logFC)>1 & TEhvl_change_heavy$adj.P.Val<0.05, "red", "black"))
abline(h=c(-1,1), lty=2, col="grey")

#data output of this includes the following: logFC  AveExpr        t       P.Value     adj.P.Val        B



#==================================================================================================================================================
#PLOTTING DATA
#==================================================================================================================================================

#Below is an example of how we subset and plotted the data. For details see the Zhulyn et al. manuscript
setwd("C:/Users/Olena/Box Sync/Work_Original/R/R_Yoga/04292018 Final Polysome Seq Data Analysis/05192018 final DE and TE")
data<-read.table("05202018_SIMPLIFIED_and_ANNOTATED_for_analysis_8139.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char="@")
head(data)
nrow(data)#8139

setwd("C:/Users/Olena/Box Sync/Work_Original/R/R_Yoga/04292018 Final Polysome Seq Data Analysis/05222018 colour groups and gene lists")

#ALL
plot(data$DE_input, data$TEhvl_change_heavy, main="Change in TE (heavy/(light+RNP) vs. change in input mRNA", pch=20, axes=FALSE, cex=0.8, xlim=c(-6,6), ylim=c(-6,6),
     xlab="DE input (log2FC input mRNA)", ylab="change TE heavy vs light+RNP (log2FC TE)")
axis(1, c(-6, -4, -2, 0, 2, 4, 6), pos = 0, cex.axis = 0.8)
axis(2, c(-6, -4, -2, 0, 2, 4, 6), pos = 0, cex.axis = 0.8)
abline(v=c(-1,1), h=c(-1,1), lty=2, col="grey")

#GREEN group = TE change; no DE change
plot(data$DE_input, data$TEhvl_change_heavy, main="Change in TE (heavy/(light+RNP) vs. change in input mRNA", pch=20, axes=FALSE, cex=0.8, xlim=c(-6,6), ylim=c(-6,6),
     xlab="DE input (log2FC input mRNA)", ylab="change TE heavy vs light+RNP (log2FC TE)", col=ifelse(abs(data$TEhvl_change_heavy)>1 & abs(data$DE_input)<1, "green", "black"))
axis(1, c(-6, -4, -2, 0, 2, 4, 6), pos = 0, cex.axis = 0.8)
axis(2, c(-6, -4, -2, 0, 2, 4, 6), pos = 0, cex.axis = 0.8)
abline(v=c(-1,1), h=c(-1,1), lty=2, col="grey")

#ORANGE Group (enhanced DE no change in TE)
plot(data$DE_input, data$TEhvl_change_heavy, main="Change in TE (heavy/(light+RNP) vs. change in input mRNA", pch=20, axes=FALSE, cex=0.8, xlim=c(-6,6), ylim=c(-6,6),
     xlab="DE input (log2FC input mRNA)", ylab="change TE heavy vs light+RNP (log2FC TE)", 
     col=ifelse(abs(data$TEhvl_change_heavy)<1 & data$DE_input>1, "orange", "black"))
axis(1, c(-6, -4, -2, 0, 2, 4, 6), pos = 0, cex.axis = 0.8)
axis(2, c(-6, -4, -2, 0, 2, 4, 6), pos = 0, cex.axis = 0.8)
abline(v=c(-1,1), h=c(-1,1), lty=2, col="grey")


#Data point counts:
GREEN<-subset(data, (abs(data$TEhvl_change_heavy)>1 & abs(data$DE_input)<1))
GREEN_UP<-subset(data, ((data$TEhvl_change_heavy)>1 & abs(data$DE_input)<1))
GREEN_DOWN<-subset(data, ((data$TEhvl_change_heavy)<(-1) & abs(data$DE_input)<1))
nrow(GREEN_UP) #504
nrow(GREEN_DOWN) #521
nrow(GREEN) #1025

ORANGE<-subset(data, (abs(data$TEhvl_change_heavy)<1 & data$DE_input>1))
nrow(ORANGE) #1190

#note changed GREY to just be transcriptionally repressed
GREY<-subset(data, (abs(data$TEhvl_change_heavy)<1 & data$DE_input<(-1)))
nrow(GREY) #156

write.table(GREEN, file="05222018_GREEN.txt", sep="\t")
write.table(GREEN_UP, file="05222018_GREEN_UP.txt", sep="\t")
write.table(GREEN_DOWN, file="05222018_GREEN_DOWN.txt", sep="\t")
write.table(ORANGE, file="05222018_ORANGE.txt", sep="\t")


#no change set
NC3<-subset(data, (abs(data$TEhvl_change_heavy)<1 & abs(data$DE_input)<1) & abs(data$TEhvl_change_RNP)<1 & abs(data$TEhvl_change_light)<1)
nrow(NC3) #5001

write.table(NC3, file="05312018_NC3_no_change_light_heavy_RNP_DE.txt", sep="\t")


