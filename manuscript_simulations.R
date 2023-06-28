#Simulations used for Di Battista et al 2023 JNI manuscript

#script for simulations accompanying the factor models in the tc manuscript

#libraries
library(tidyverse)
library(ggdist)
library(magrittr)
library(modelr)
library(ggplot2)
library(cowplot)
library(rstan)
library(rstanarm)
library(bayesplot)
library(rethinking)
library(tidybayes)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(gtsummary)
library(gt)
library(loo)
library(ggdag)
library(dagitty)
library(MASS)
library(corrplot)

#functions
HDILow<- function(x, HDI=0.9) {
  sortedPts = sort( x)
  ciIdxInc = ceiling( HDI * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc 
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ] }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ] 
  return( HDImin)
}
HDIHigh<- function(x, HDI=0.9) {
  sortedPts = sort( x)
  ciIdxInc = ceiling( HDI * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc 
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ] }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ] 
  return( HDImax)
}
z_score <- function(x) {
  out <- (x - mean(x,na.rm=TRUE))/sd(x,na.rm = TRUE)}
mod_z <- function(x) {
  out <- (x - median(x,na.rm=TRUE))/mad(x,na.rm = TRUE)}

#for speed with mcmc
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#generate a simulated dataset for EFA
#taken from Rick Farouli's website
#https://rfarouni.github.io/assets/projects/BayesianFactorAnalysis/BayesianFactorAnalysis.html
set.seed(10)
D <-1 #number of dimensions
P <- 30  #number of variables (biomarkers)
N <- 60 # number of observations

mu_theta <-rep(0,D) # the mean of eta
mu_epsilon <-rep(0,P) # the mean of epsilon
Phi<-diag(rep(1,D))
Psi <- diag(c(0.20, 0.19, 0.15, 0.20, 0.36, 0.18, 0.18, 1.00, 0.27, 0.27,
              0.65, 0.19, 0.18, 0.15, 0.06, 0.18, 0.22, 0.10, 0.57, 0.24,
              0.31, 0.22, 0.07, 0.10, 0.23, 0.40, 0.19, 0.23, 0.27, 0.37))
l1 <- c(0.85, 0.00, 0.30, 0.00, 0.80, 0.00, -0.50, 0.00, 0.00, 0.00,
        0.45, 0.00, 0.10, 0.00, 0.80, 0.00, -0.50, 0.00, -0.30, 0.00,
        0.85, 0.20, 0.30, 0.00, -0.70, 0.00, 0.50, 0.00, 0.00, 0.00)

Theta <-mvrnorm(N, mu_theta, Phi) # sample factor scores
Epsilon <-mvrnorm(N, mu_epsilon, Psi) # sample error vector
Y<-Theta%*%t(l1)+Epsilon# generate observable data

#are the data scaled?
apply(Y,2,mean)

#scale
Y <- apply(Y,2,z_score)

#visualize
corrplot(cor(Y))

#L and D prior model
dat_mod_ld <- list(m = P, k = D, n = N, n_obs = N*P,y = Y)
str(dat_mod_ld)

#run the model
fact_mod_ld <- stan_model("factor_analysis_cc.stan")
m_ld <- sampling(object = fact_mod_ld, data = dat_mod_ld, chains = 4,
                 iter = 3000 ,init_r = .2) # add some init

#look at posterior
print(m_ld)
post_m_ld <- data.frame(extract.samples(m_ld))
precis(post_m_ld,2)

print(diag(Psi)) 

#are factor scores recovered?
post_z <- post_m_ld[grepl("z",colnames(post_m_ld))]
colnames(post_z)
post_z_means <- apply(post_z,2,mean)

#raw factor scores
Theta

#high correlation!
cor(Theta,post_z_means)

#sample the prior
#update April 6th, 2023 - model priors
dat_ld_prior <- list(m = P, k = D, n = N, n_obs = N*P,y = as.matrix(Y)
)
str(dat_ld_prior)

#compile
fact_ld_prior <-  stan_model("factor_analysis_prior.stan")

#run prior model
m_ld_prior <- sampling(object = fact_ld_prior, data = dat_ld_prior, chains = 4, 
                       iter = 3000 ,init_r = .2) # add some init
prior_m_ld <- data.frame(extract.samples(m_ld_prior))

#Visualize the posterior of a factor model (this is also useful for pub figures)
#trankplots
trankplot(m_ld,pars = "L")
trankplot(m_ld,pars = "z")

#Compare posterior factor loadings
#create a df for factor loadings
post_L <- post_m_ld[grepl("L",colnames(post_m_ld))&
                      !(grepl("Sigma",colnames(post_m_ld)))]
#sigma leftover, just remove.
post_L$sigma_L <- NULL

#create separate dataframes, and then combine
f1 <- post_L[c(1:30)]
colnames(f1) <- colnames(as.data.frame(Y))

#transform to long
f1_long <- gather(f1,key = 'variable', value = 'loading score')
#add factor #
f1_long$factor <- 'Factor 1'


plot_df <- f1_long
plot_df$variable <- factor(plot_df$variable,levels = c('V1','V2','V3','V4','V5',
                                                       'V6','V7','V8','V9','V10',
                                                       'V11','V12','V13','V14',
                                                       'V15','V16','V17','V18',
                                                       'V19','V20','V21','V22',
                                                       'V23','V24','V25','V26',
                                                       'V27','V28','V29','V30'))

L <- as.data.frame(l1)
rownames(L) <- colnames(as.data.frame(Y)) 
L <- as.data.frame(L)
L$variable <- rownames(L)

colnames(L) <- c('Score','Variable')
point_plot <- L
point_plot$Variable <- factor(point_plot$Variable,levels = c('V1','V2','V3','V4','V5',
                                                             'V6','V7','V8','V9','V10',
                                                             'V11','V12','V13','V14',
                                                             'V15','V16','V17','V18',
                                                             'V19','V20','V21','V22',
                                                             'V23','V24','V25','V26',
                                                             'V27','V28','V29','V30'))

#Plot with raw data overlay
theme_set(theme_tidybayes() + panel_border())
plot_loadings <- plot_df %>%
  ggplot(aes(y = fct_rev(variable),
             x = `loading score`,
             fill = stat(x<0))) + 
  geom_vline(xintercept = 0, alpha = 0.8, linetype = 2) + 
  theme(strip.text.x = element_text(face = 'bold', size = 10),
        axis.title.x = element_text(size = 13, face = 'bold'),
        axis.title.y = element_text(size = 13, face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 10),
        legend.position = 'None',
        axis.text.x = element_text(face = 'bold', size = 10)) +
  ylab("Variable") +
  facet_grid(~factor) + xlab("Loading Scores")+
  stat_halfeye(.width = c(0.70,0.90), alpha = 0.7) +
  geom_point(data = point_plot, 
             aes(x = `Score`, y = fct_rev(Variable)),
             shape = 22,fill = 'yellow', color = "black", size = 2)
plot_loadings
ggsave('loadings_plot.jpg',plot_loadings,dpi =600, width = 10)

#plot priors
#Compare posterior factor loadings
#create a df for factor loadings
prior_L <- prior_m_ld[grepl("L",colnames(prior_m_ld))&
                        !(grepl("Sigma",colnames(prior_m_ld)))]
#sigma leftover, just remove.
prior_L$sigma_L <- NULL

#create separate dataframes, and then combine
f1 <- prior_L[c(1:30)]
colnames(f1) <- colnames(as.data.frame(Y))

#transform to long
f1_long <- gather(f1,key = 'variable', value = 'loading score')
#add factor #
f1_long$factor <- 'Factor 1'


prior_df <- f1_long
prior_df$variable <- factor(prior_df$variable,levels = c('V1','V2','V3','V4','V5',
                                                         'V6','V7','V8','V9','V10',
                                                         'V11','V12','V13','V14',
                                                         'V15','V16','V17','V18',
                                                         'V19','V20','V21','V22',
                                                         'V23','V24','V25','V26',
                                                         'V27','V28','V29','V30'))


#Plot 
theme_set(theme_tidybayes() + panel_border())
plot_prior_loadings <- prior_df %>%
  ggplot(aes(y = fct_rev(variable),
             x = `loading score`,
             fill = stat(x<0))) + 
  geom_vline(xintercept = 0, alpha = 0.8, linetype = 2) + 
  theme(strip.text.x = element_text(face = 'bold', size = 10),
        axis.title.x = element_text(size = 13, face = 'bold'),
        axis.title.y = element_text(size = 13, face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 10),
        legend.position = 'None',
        axis.text.x = element_text(face = 'bold', size = 10)) +
  ylab("Variable") +
  facet_grid(~factor) + xlab("Loading Scores")+
  stat_halfeye(.width = c(0.70,0.90), alpha = 0.7) 
plot_prior_loadings

#model together with posterior and raw data
#add sample column to df
plot_df$sample <- 'Posterior' 
prior_df$sample <- "Prior"
combined_ml_df <- rbind.data.frame(plot_df,prior_df)

#plot combined
theme_set(theme_tidybayes() + panel_border())
plot_combined_loadings <- combined_ml_df %>%
  ggplot(aes(y = fct_rev(variable),
             x = `loading score`,
             fill = sample)) + 
  geom_vline(xintercept = 0, alpha = 0.8, linetype = 2) + 
  theme(strip.text.x = element_text(face = 'bold', size = 10),
        axis.title.x = element_text(size = 13, face = 'bold'),
        axis.title.y = element_text(size = 13, face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 10),
        legend.position = 'None',
        axis.text.x = element_text(face = 'bold', size = 10)) +
  ylab("Variable") +
  facet_grid(~factor) + xlab("Loading Scores")+
  stat_halfeye(.width = c(0.70,0.90), alpha = 0.7) +xlim(-1.5,1.5)+
  scale_fill_manual(values = c('red','grey'))+
  geom_point(data = point_plot, 
             aes(x = `Score`, y = fct_rev(Variable)),
             shape = 22,fill = 'yellow', color = "black", size = 2)
plot_combined_loadings
ggsave('combined_loadings_plot.jpg',plot_combined_loadings,dpi =600, width = 10)

#plotting a correlation instead
f1_means <- apply(post_z[c(1:60)],2,mean)
f1_low <- apply(post_z[c(1:60)],2,HDILow)
f1_high <- apply(post_z[c(1:60)],2,HDIHigh)


#create a plot_df2 that contains info for raw and estimated values
fac1 <- cbind.data.frame(Theta,f1_means,f1_low,f1_high)
colnames(fac1) <- c('raw','means','low','high')
fac1$factor <- 'Factor 1'

plot_df2 <- fac1
summary(plot_df2$means)

#plot
factor_plot <-ggplot(plot_df2,aes(x = means,y =raw))+
  geom_point(shape = 21,size = 1.5,stroke =2,colour = 'red')+
  geom_linerange(aes(xmin=low, xmax=high), 
                 colour="red",size = 3,alpha = 0.4,lineend='round') +
  ylab('Raw Factor Scores')+
  xlab('Estimated Factor Scores')+
  theme(title=element_text(face = 'bold'),
        strip.text = element_text(face ='bold',size = 10))+
  geom_smooth(method = 'lm')+ 
  facet_grid(~factor)
factor_plot
ggsave('factor_plot.jpg',factor_plot,dpi = 600)

#data correlation structure (biomarker correlations)
#need to convert the covariances (Sigma) by calculating first the sd of the diagonal...
#error terms for each biomarker (sigma_y)

#first pull the relevant values from the posterior
df_Sigma <- post_m_ld[grepl('Sigma', colnames(post_m_ld))]
df_Sigma <- df_Sigma[!(grepl('L',colnames(df_Sigma)))]
#derive the means
vec_Sigma_means <- apply(df_Sigma,2,mean)

#create a covariance matrix
cov_mat <- matrix(vec_Sigma_means,nrow = 30, ncol = 30)
# Convert to correlation matrix
cor_mat <- cov2cor(cov_mat)

#plot
corrplot(cor_mat)

#how does this compare to the regular data structure?
corrplot(cor(Y))

##Simulating how group differences in latent structures can be uncovered 

#look at Y first
#posterior loading plot of Y suggests Factor 1 is being 'carried' in variance by...
#V1, V5, V15, V17, V21. Therefore, will create two groups based on the values of
#these important variables.

#turn to df first
Y <- as.data.frame(Y)

Y$group <- 1
Y$group[Y$V1>median(Y$V1) & 
          Y$V5>median(Y$V5) &
          Y$V7<median(Y$V7)  ] <- 2

#how many values does that give us?
table(Y$group)

#look at difference in raw factor scores first
Theta
Y$group
latent_score_df <- cbind.data.frame(Y$group,Theta)
by(latent_score_df$Theta,factor(latent_score_df$`Y$group`),mean)

#after pulling z out of the posterior, clean it to feed into the second model
#code to retrieve a subset of random draws of post_z
n_subsets = 100
post_size = 2000
subset_rows <- sample(post_size,n_subsets,replace = F)
z1_subset <- post_z[subset_rows,c(1:60)]

#code to derive a 60 x 1 vector of lv1 draws
z1_samples <- list()
for (i in 1:100) {
  z1_samples[[i]] <- as.vector(t(z1_subset[i,]))
}

#set up stan data list
stan_data1 = list()
for (i in 1:100){
  n = 60 #subjects
  group = Y$group
  z = z1_samples[[i]]
  stan_data1[[i]] <- list(n = n, group = group, z = z)}

#modelling - compile first
linear_from_z <- stan_model('linear_from_z.stan')

#try the model first on the mean of post_z

stan_data1_mean <- list(n = n, group = group, 
                        z = as.vector(z_score(post_z_means) ))
str(stan_data1_mean)
#mean model
test_fit_mean1 <- sampling(linear_from_z,data = stan_data1_mean, chains = 4, iter = 2000)
test_fit_mean1

post_test_fit_mean1 <- as.data.frame(extract.samples(test_fit_mean1))

#try the model on one of the stan_data draws
test_fit_draw1 <- sampling(linear_from_z,data = stan_data1[[1]], chains = 4, iter = 500)
test_fit_draw1

# Run the loop to generate 100 samples
draws1 <- list()
for (i in 1:100) {
  draws1[[i]] <- sampling(linear_from_z, data = stan_data1[[i]],chains = 4, iter = 200)
}

post_draws1<- lapply(draws1,as.data.frame(extract.samples))

post_draws1_df <- do.call(rbind,post_draws1)

# try with smaller # of iterations
draws1_100 <- list()
for (i in 1:100) {
  draws1_100[[i]] <- sampling(linear_from_z, data = stan_data1[[i]],chains = 4, iter = 100)
}

# Check the output of the first iteration
print(summary(draws1_100[[50]]))

trankplot(draws1_100[[10]])
#concatenate draws
post_draws1_100 <- lapply(draws1_100,as.data.frame(extract.samples))

post_draws1_100_df <- do.call(rbind,post_draws1_100)
precis(post_draws1_100_df,3)

post_draws1_100_df$lv1_contrast <- post_draws1_100_df$value.a.2 - 
  post_draws1_100_df$value.a.1

# try with smaller # of iterations
draws1_50 <- list()
for (i in 1:100) {
  draws1_50[[i]] <- sampling(linear_from_z, data = stan_data1[[i]],chains = 4, iter = 50)
}

# Check the output of the first iteration
print(summary(draws1_50[[50]]))

trankplot(draws1_50[[10]])
#concatenate draws
post_draws1_50 <- lapply(draws1_50,as.data.frame(extract.samples))

post_draws1_50_df <- do.call(rbind,post_draws1_50)
precis(post_draws1_50_df,3)

post_draws1_50_df$lv1_contrast <- post_draws1_50_df$value.a.2 - 
  post_draws1_50_df$value.a.1


#single linear secondary model from a 1 LV Factor model
#200 draws
df_ln_2step <- post_draws1_df[c(1,2)]
df_ln_2step$lv1_contrast <- df_ln_2step$value.a.2 - df_ln_2step$value.a.1
colnames(df_ln_2step) <- c('LV1 Group 1','LV1 Group 2','LV1 Contrast')
df_ln_2step_long <- gather(df_ln_2step,key = "Condition",value = "Value")
df_ln_2step_long$Model <- 'Linear 200 Draws'

#100 draws
df_ln_2step_100 <- post_draws1_100_df[c(1,2,5)]
df_ln_2step_100$lv1_contrast <- df_ln_2step_100$value.a.2 - df_ln_2step_100$value.a.1
colnames(df_ln_2step_100) <- c('LV1 Group 1','LV1 Group 2','LV1 Contrast')
df_ln_2step_100_long <- gather(df_ln_2step_100,key = "Condition",value = "Value")
df_ln_2step_100_long$Model <- 'Linear 100 Draws'

#50 draws
df_ln_2step_50 <- post_draws1_50_df[c(1,2,5)]
df_ln_2step_50$lv1_contrast <- df_ln_2step_50$value.a.2 - df_ln_2step_50$value.a.1
colnames(df_ln_2step_50) <- c('LV1 Group 1','LV1 Group 2','LV1 Contrast')
df_ln_2step_50_long <- gather(df_ln_2step_50,key = "Condition",value = "Value")
df_ln_2step_50_long$Model <- 'Linear 50 Draws'

#model of posterior means for single linear model with 1 LV
colnames(post_test_fit_mean1)
df_ln_mean <- post_test_fit_mean1[c(1,2)]
df_ln_mean$lv1_contrast <- df_ln_mean$a.2 - df_ln_mean$a.1
colnames(df_ln_mean) <- c('LV1 Group 1','LV1 Group 2','LV1 Contrast')
df_ln_mean_long <- gather(df_ln_mean,key = "Condition",value = "Value")
df_ln_mean_long$Model <- 'Linear Posterior Means'

#Validating the model works with missing data.
#create a dataset missing the lowest 20% of data for 8 variables.
Y_miss <- Y[-c(31)]
Y_miss[c(1,3,5,7,15,20,25,30)] <-
  apply(Y_miss[c(1,3,5,7,15,20,25,30)],2, 
        function(x) ifelse(x>quantile(x,probs =0.20),x,NA))
#create a subset with just the missing data columns to make it easier to do some...
#downstream math
Y_miss2 <- Y_miss[c(1,3,5,7,15,20,25,30)]

#create bounds
Y_means <- colMeans(Y_miss2,na.rm = TRUE)
Y_sds <- apply(Y_miss2,2,sd,na.rm=TRUE)

#now create a lowerbound vector on the zscale
zeroes <- sapply(1:ncol(Y_miss2),function(i) (0 - Y_means[i])/Y_sds[i] )
cut_values <- sapply(1:ncol(Y_miss2),function(i) (min(Y_miss2[,i],na.rm = TRUE) -
                                                    Y_means[i])/Y_sds[i])

n = nrow(Y_miss)
m = ncol(Y_miss) 
n_miss <- length(Y_miss[is.na(Y_miss)])

#isolate missing from observed
Y_miss <-data.frame(apply(Y_miss,2,z_score))
na_values <- is.na(Y_miss)
y_obs <- Y_miss[!na_values]
n_obs <- length(y_obs)
row_miss <- rep(1:n, times = m)[na_values]
col_miss<- rep(1:m, each = n)[na_values]
miss_check <- cbind(row_miss,col_miss) 

#identifies index where observed values are
row_obs <- rep(1:n, times = m)[!na_values]
col_obs <- rep(1:m, each = n)[!na_values]

bounds <- data.frame(Lo = zeroes,Up = cut_values)
bounds$col_miss <- c(1,3,5,7,15,20,25,30)
col_miss1 <- as.data.frame(col_miss)
bounds1 <- merge(bounds,col_miss1, by = 'col_miss',all =TRUE)
lo <- as.vector(bounds1$Lo)
up <- as.vector(bounds1$Up)

dat_Y_miss = list(N = nrow(Y_miss),
                  P = ncol(Y_miss),
                  row_obs = row_obs,
                  k = 1,
                  m = 30,
                  n = 60,
                  col_obs = col_obs,
                  row_miss = row_miss,
                  col_miss = col_miss,
                  n_miss = n_miss,
                  n_obs = n_obs,
                  y_obs = y_obs,
                  y = as.matrix(Y_miss),
                  Lo = lo,
                  Up = up)

str(dat_Y_miss)

#run the imputed model

mod_impute <- stan_model('factor_analysis_impute_low.stan')

m_miss <- sampling(object = mod_impute, data = dat_Y_miss, chains = 4, cores=4,
                   iter = 2000 ,init_r = .2) # add some init
m_miss

#Run lines 100 - 464 to derive plots, correlations and subsequent linear modelling.

###even more missing data (up to 50% in some columns)
Y_miss50 <- Y[-c(31)]
Y_miss50[c(1,3,5,7,15,20,25,30)] <-
  apply(Y_miss50[c(1,3,5,7,15,20,25,30)],2, 
        function(x) ifelse(x>quantile(x,probs =0.50),x,NA))
#create a subset with just the missing data columns to make it easier to do some...
#downstream math
Y_miss2_50 <- Y_miss50[c(1,3,5,7,15,20,25,30)]

#create bounds
Y_means50 <- colMeans(Y_miss2_50,na.rm = TRUE)
Y_sds50 <- apply(Y_miss2_50,2,sd,na.rm=TRUE)

#now create a lowerbound vector on the zscale
zeroes50 <- sapply(1:ncol(Y_miss2_50),function(i) (0 - Y_means50[i])/Y_sds50[i] )
cut_values50 <- sapply(1:ncol(Y_miss2_50),function(i) (min(Y_miss2_50[,i],na.rm = TRUE) -
                                                         Y_means50[i])/Y_sds50[i])

n = nrow(Y_miss50)
m = ncol(Y_miss50) 
n_miss50 <- length(Y_miss50[is.na(Y_miss50)])

#isolate missing from observed
Y_miss50 <-data.frame(apply(Y_miss50,2,z_score))
na_values50 <- is.na(Y_miss50)
y_obs50 <- Y_miss50[!na_values50]
n_obs50 <- length(y_obs50)
row_miss50 <- rep(1:n, times = m)[na_values50]
col_miss50 <- rep(1:m, each = n)[na_values50]
miss_check50 <- cbind(row_miss50,col_miss50) 

#identifies index where observed values are
row_obs50 <- rep(1:n, times = m)[!na_values50]
col_obs50 <- rep(1:m, each = n)[!na_values50]

bounds50 <- data.frame(Lo = zeroes50,Up = cut_values50)
bounds50$col_miss50 <- c(1,3,5,7,15,20,25,30)
col_miss1_50 <- as.data.frame(col_miss50)
bounds1_50 <- merge(bounds50,col_miss1_50, by = 'col_miss50',all =TRUE)
lo50 <- as.vector(bounds1_50$Lo)
up50 <- as.vector(bounds1_50$Up)

dat_Y_miss50 = list(N = nrow(Y_miss50),
                    P = ncol(Y_miss50),
                    row_obs = row_obs50,
                    k = 1,
                    m = 30,
                    n = 60,
                    col_obs = col_obs50,
                    row_miss = row_miss50,
                    col_miss = col_miss50,
                    n_miss = n_miss50,
                    n_obs = n_obs50,
                    y_obs = y_obs50,
                    y = as.matrix(Y_miss50),
                    Lo = lo50,
                    Up = up50)

str(dat_Y_miss50)

#run the imputed model
m_miss50 <- sampling(object = mod_impute, data = dat_Y_miss50, chains = 4, cores=4,
                     iter = 2000 ,init_r = .2) # add some init
m_miss50

#Run lines 100 - 464 to derive plots, correlations and subsequent linear modelling.


