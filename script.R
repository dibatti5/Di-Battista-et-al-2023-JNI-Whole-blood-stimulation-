#Code to run the analyses in Di Battista et al 2023 JNI

#libraries
library(tidyverse)
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
library(ggdag)

#functions
HDILow<- function(x, HDI=0.95) {
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
HDIHigh<- function(x, HDI=0.95) {
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

#load files
#demo file
demos_df <- read.csv('demos_df.csv',stringsAsFactors = FALSE)
#biomarker file
biomarkers_df <- read.csv('biomarkers_df.csv',stringsAsFactors = FALSE)


#DAGs for scientific models
m1 <- dagitty('dag {
"Conc Hx" [pos="-0.35,-0.252"]
Ex [latent,pos="-0.607,-0.257"]
Imm [outcome,pos="-0.095,-0.254"]
SRC [exposure,pos="-0.610,-0.254"]
Sex [pos="-0.096,-0.257"]
"Conc Hx" -> Imm
"Conc Hx" -> SRC
Ex -> Imm
SRC -> Ex
SRC -> Imm
Sex -> Imm
Sex -> SRC
}
')

drawdag(m1)
impliedConditionalIndependencies(m1) 

#save
d1_plot <- ggdag(m1,text_size = 3)+theme_dag()
d1_plot
ggsave('d1_plot.jpg',d1_plot,dpi = 600)



###LPS modelling
#separate by conditions and add some demos for downstream causal analysis
data_lps_df <- biomarkers_df[biomarkers_df$group=='LPS',]

#find % missing data in a named vector called 'vec'
vec_lps <- apply(data_lps_df[c(5:49)],2,function(x) (length(x[is.na(x)])/length(x))*100)
names_remove <- names(which(vec_lps>50))

#remove less than 50% cytokines
data_lps_df <- data_lps_df[!(colnames(data_lps_df)%in%names_remove)]

####Prep LPS latent model
lps_df <- merge(demos_df,data_lps_df, by = 'id.blood',all.y=TRUE)

#clean up
colnames(lps_df)

#separate the biomarkers to create a latent variable
bio_lps <- lps_df[c(8:32)]

#which columns are missing
#find % missing data in a named vector called 'vec'
miss_lps <- apply(bio_lps,2,function(x) (length(x[is.na(x)])/length(x))*100)
names_miss_lps <- names(which(miss_lps>0))

bio_lps_miss <- bio_lps[colnames(bio_lps)%in%names_miss_lps]

#create bounds for the imputed parameters
bio_lps_means <- colMeans(bio_lps_miss,na.rm = TRUE)
bio_lps_sds <- apply(bio_lps_miss,2,sd,na.rm=TRUE)
#now create a lowerbound vector on the zscale
zeroes_lps <- sapply(1:ncol(bio_lps_miss),function(i) (0 - bio_lps_means[i])/bio_lps_sds[i] )
cut_values_lps <- sapply(1:ncol(bio_lps_miss),function(i) (min(bio_lps_miss[,i],na.rm = TRUE) -
                                                             bio_lps_means[i])/bio_lps_sds[i])

n_lps = nrow(bio_lps)
mc_lps = ncol(bio_lps)
n_miss_lps <- length(bio_lps[is.na(bio_lps)])

#isolate missing from observed
bio_lps_z <-data.frame(apply(bio_lps,2,z_score))
na_values_lps <- is.na(bio_lps_z)
y_obs_lps <- bio_lps_z[!na_values_lps]
n_obs_lps <- length(y_obs_lps)
row_miss_lps <- rep(1:n_lps, times = mc_lps)[na_values_lps]
col_miss_lps<- rep(1:mc_lps, each = n_lps)[na_values_lps]
miss_check_lps <- cbind(row_miss_lps,col_miss_lps) 

#identifies index where observed values are
row_obs_lps <- rep(1:n_lps, times = mc_lps)[!na_values_lps]
col_obs_lps <- rep(1:mc_lps, each = n_lps)[!na_values_lps]


bounds_lps <- data.frame(Lo = zeroes_lps,Up = cut_values_lps)
bounds_lps$col_miss_lps <- unique(col_miss_lps) 
col_miss_lps1 <- as.data.frame(col_miss_lps)
bounds1_lps <- merge(bounds_lps,col_miss_lps1, by = 'col_miss_lps',all =TRUE)
lo_lps <- as.vector(bounds1_lps$Lo)
up_lps <- as.vector(bounds1_lps$Up)

dat_lps = list(n = nrow(bio_lps),
               m = ncol(bio_lps),
               k = 1,
               row_obs = row_obs_lps,
               col_obs = col_obs_lps,
               row_miss = row_miss_lps,
               col_miss = col_miss_lps,
               n_miss = n_miss_lps,
               n_obs = n_obs_lps,
               y_obs = y_obs_lps,
               y = as.matrix(bio_lps_z),
               Lo = lo_lps,
               Up = up_lps)

str(dat_lps)

#compile
factor_impute_low <- stan_model("factor_analysis_impute_low.stan")

#run
m_lps<- sampling(object = factor_impute_low, data = dat_lps, chains = 4,
                 iter = 3000,init_r=.2 )

#extract posterior
m_post_lps <- data.frame(extract.samples(m_lps))

precis(m_post_lps,3)

#check chains
trankplot(m_lps, pars = c('z'))

#model prior
#create fake matrix 
y_prior_lps <- matrix(0,nrow=n_lps,ncol=mc_lps)

dat_prior_lps <-  list(
  k = 1,
  m = mc_lps,
  n = n_lps,
  n_obs = n_lps*mc_lps,
  y = y_prior_lps)

str(dat_prior_lps)

mod_prior <- stan_model('factor_analysis_prior.stan')

bio_prior_lps<- sampling(object = mod_prior, data = dat_prior_lps, chains = 4, cores=4,
                         iter = 3000 ,init_r = .2) # add some init
bio_prior_lps

##########Set up plots
###Posterior
#Compare post factor loadings
#create a df for factor loadings
post_L_lps <- m_post_lps[grepl("L",colnames(m_post_lps))&
                           !(grepl("Sigma",colnames(m_post_lps)))]
#sigma leftover, just remove.
post_L_lps$sigma_L <- NULL

#create separate dataframes, and then combine
f1_lps<- post_L_lps[c(1:mc_lps)]
colnames(f1_lps) <- colnames(bio_lps)

#transform to long
f1_lps_long <- gather(f1_lps,key = 'biomarker', value = 'Loading score')
#add factor #
f1_lps_long$factor <- 'Factor 1'
f1_lps_long$biomarker<- factor(f1_lps_long$biomarker,
                               levels = c('IFNG','IL1B','IL2','IL4','IL6','IL7','IL10','IL13','IL15','IL17A','IL17C','IL17F',
                                          'IL18','IL27','IL33','TNF','TNFSF10','TNFSF12','FLT3LG','LTA','CSF1','CSF2','CSF3','OSM','TSLP',
                                          'CCL2','CCL3','CCL4','CCL7','CCL8','CCL11','CCL13','CCL19','CXCL8','CXCL9',
                                          'CXCL10','CXCL11','CXCL12','HGF','MMP1','MMP12','OLR1','EGF','TGFA','VEGFA'))

#prior
#extract the prior
prior_lps <- data.frame(extract.samples(bio_prior_lps))
precis(prior_lps,2)

#Compare prior factor loadings
#create a df for factor loadings
prior_L_lps <- prior_lps[grepl("L",colnames(prior_lps))&
                           !(grepl("Sigma",colnames(prior_lps)))]
#sigma leftover, just remove.
prior_L_lps$sigma_L <- NULL

#create separate dataframes, and then combine
f1_lps_prior<- prior_L_lps[c(1:mc_lps)]
colnames(f1_lps_prior) <- colnames(bio_lps)

#transform to long
f1_lps_long_prior <- gather(f1_lps_prior,key = 'biomarker', value = 'Loading score')
#add factor #
f1_lps_long_prior$factor <- 'Factor 1'

f1_lps_long_prior$biomarker <- factor(f1_lps_long_prior$biomarker,
                                      levels = c('IFNG','IL1B','IL2','IL4','IL6','IL7','IL10','IL13','IL15','IL17A','IL17C','IL17F',
                                                 'IL18','IL27','IL33','TNF','TNFSF10','TNFSF12','FLT3LG','LTA','CSF1','CSF2','CSF3','OSM','TSLP',
                                                 'CCL2','CCL3','CCL4','CCL7','CCL8','CCL11','CCL13','CCL19','CXCL8','CXCL9',
                                                 'CXCL10','CXCL11','CXCL12','HGF','MMP1','MMP12','OLR1','EGF','TGFA','VEGFA'))

#combine
f1_lps_long$samples <- 'Posterior'
f1_lps_long_prior$samples <-'Prior'

plot_combined_lps <- rbind.data.frame(f1_lps_long,f1_lps_long_prior)

#use this if interpretation is not useful
plot_combined_lps$group <- 'Cytokines'
plot_combined_lps$group[grepl('CC',plot_combined_lps$biomarker)|
                          grepl('CXC',plot_combined_lps$biomarker)] <- 'Chemokines'
plot_combined_lps$group[plot_combined_lps$biomarker=='MMP12'|
                          plot_combined_lps$biomarker=='MMP1'|
                          plot_combined_lps$biomarker=='OLR1'|
                          plot_combined_lps$biomarker=='OSM'|
                          plot_combined_lps$biomarker=='TGFA'|
                          plot_combined_lps$biomarker=='EGF'|
                          plot_combined_lps$biomarker=='HGF'|
                          plot_combined_lps$biomarker=='VEGFA'] <- 'Other'

#combined plot
theme_set(theme_tidybayes() + panel_border())
plot_lps<- plot_combined_lps %>%
  ggplot(aes(y = fct_rev(`biomarker`),
             x = `Loading score`,
             fill = samples)) + 
  geom_vline(xintercept = 0, alpha = 0.8, linetype = 2) + 
  theme(title = element_text(face = 'bold'),
        axis.title.x = element_text(size = 13, face = 'bold'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = 'bold', size = 10),
        legend.position = 'None',
        axis.text.x = element_text(face = 'bold', size = 10)) +
  ggtitle("LPS Stimulation")+xlim(-2,2)+
  scale_fill_manual(values = c('red','grey'))+
  stat_halfeye(aes(color = samples),.width = c(0.70,0.90), 
               alpha = 0.5,linewidth = 1.5,point_size = 1) +
  scale_colour_manual(values = c('red','grey'))

plot_lps

ggsave('plot_lps.jpg',plot_lps,dpi = 600,width = 5,height = 10)

#########Biomarker covariance structure
#need to convert the covariances (Sigma) by calculating first the sd of the diagonal...
#error terms for each biomarker (sigma_y)
#first pull the relevant values from the posterior
df_Sig_lps <- m_post_lps[grepl('Sigma', colnames(m_post_lps))]

#remove the "L"
df_Sig_lps <- df_Sig_lps[!(grepl('L',colnames(df_Sig_lps)))]

#derive the means
vec_Sig_lps_means<- apply(df_Sig_lps,2,mean)

#create a covariance matrix
cov_mat_lps <- matrix(vec_Sig_lps_means,nrow = 25, ncol = 25)

#rename
colnames(cov_mat_lps) <- colnames(bio_lps)
rownames(cov_mat_lps) <- colnames(bio_lps)

# Convert to correlation matrix
cor_mat_lps <- cov2cor(cov_mat_lps)

#plot
corrplot(cor_mat_lps)

#how does this compare to the regular data structure?
corrplot(cor(bio_lps_z,use = "pairwise.complete.obs"))
#by group
corrplot(cor(bio_lps_z[lps_df$id.group.x=='healthy',],use = "pairwise.complete.obs"))

######DOWNSTREAM MODELLING of the latent variable scores created.

#TOTAL EFFECTS

#first create some necessary variables 
#create binary conchx
lps_df$conchx_bin <- as.integer(ifelse(lps_df$demos.hx_conc_number==0,1,2))

#create interaction term
lps_df$int_term <- 'healthy_male_no'   #1
lps_df$int_term[lps_df$id.group.x=='healthy'&
                   lps_df$demos.sex==0&
                   lps_df$conchx_bin==2] <-'healthy_male_yes' #2
lps_df$int_term[lps_df$id.group.x=='healthy'&
                   lps_df$demos.sex==1&
                   lps_df$conchx_bin==1] <-'healthy_female_no' #3
lps_df$int_term[lps_df$id.group.x=='healthy'&
                   lps_df$demos.sex==1&
                   lps_df$conchx_bin==2] <-'healthy_female_yes'#4
lps_df$int_term[lps_df$id.group.x=='concussion'&
                   lps_df$demos.sex==0&
                   lps_df$conchx_bin==1] <-'src_male_no' #5
lps_df$int_term[lps_df$id.group.x=='concussion'&
                   lps_df$demos.sex==0&
                   lps_df$conchx_bin==2] <-'src_male_yes' #6
lps_df$int_term[lps_df$id.group.x=='concussion'&
                   lps_df$demos.sex==1&
                   lps_df$conchx_bin==1] <-'src_female_no' #7
lps_df$int_term[lps_df$id.group.x=='concussion'&
                   lps_df$demos.sex==1&
                   lps_df$conchx_bin==2] <-'src_female_yes' #8
table(lps_df$int_term)
lps_df$int_term <- factor(lps_df$int_term,levels = c('healthy_male_no',
                                                       'healthy_male_yes','healthy_female_no','healthy_female_yes',
                                                       'src_male_no','src_male_yes','src_female_no','src_female_yes'))


#using mean modelling to 1) compare unpooled and pooled estimates, and 2) derive
#the prior predictive checks.

#isolate z scores
colnames(m_post_lps)
post_z_lps <- m_post_lps[grepl("z",colnames(m_post_lps))]

#derive the data list
dat_mean_total_lps = list(
  z = z_score(as.numeric(apply(post_z_lps,2,mean))),
  conchx = lps_df$conchx_bin,
  group = as.integer(as.factor(lps_df$id.group.x)),
  sex = as.integer(lps_df$demos.sex) +1,
  int_term = as.integer(lps_df$int_term),
  n = ncol(post_z_lps)
)

str(dat_mean_total_lps)

#compile the total effects models
student_total_interact_pooled <- stan_model('student_total_interact_pooled.stan')
student_total_interact_unpooled <- stan_model('student_total_interact_unpooled.stan')

#run both mean models
m_stip <- sampling(student_total_interact_pooled,dat = dat_mean_total_lps,
                   chains = 4, iter = 1000)
m_stiu <- sampling(student_total_interact_unpooled,dat = dat_mean_total_lps,
                   chains = 4, iter = 1000)

#evaluate models out of bag
loo_stip <- loo(m_stip)
loo_stiu <- loo(m_stiu)

#compare
loo_compare(loo_stip,loo_stiu)
rethinking::compare(m_stip,m_stiu,func = PSIS)

#extract posteriors
post_stip <- data.frame(extract.samples(m_stip))
post_stiu <- data.frame(extract.samples(m_stiu))

precis(post_stip,2,prob = 0.9)
precis(post_stiu,2,prob = 0.9)

#run priors to overlay on plots
#create data for prior
dat_mean_total_lps_prior = list(
  z = z_score(as.numeric(apply(post_z_lps,2,mean)))*0,
  conchx = lps_df$conchx_bin,
  group = as.integer(as.factor(lps_df$id.group.x)),
  sex = as.integer(lps_df$demos.sex) +1,
  int_term = as.integer(lps_df$int_term),
  n = ncol(post_z_lps)
)

str(dat_mean_total_lps_prior)

student_total_interact_pooled_prior <- 
  stan_model('student_total_interact_pooled_prior.stan')

#run prior model
m_stip_prior <- sampling(student_total_interact_pooled_prior,
                         dat = dat_mean_total_lps_prior,
                         chains = 4, iter = 1000)

#model the total effects from DAG1 with the interaction term pooled
#full model with 100 draws from the latent model posterior
n_subsets = 200
post_size = 6000
subset_rows <- sample(post_size,n_subsets,replace = F)
z1_subset_lps <- post_z_lps[subset_rows,]

z1_samples_lps <- list()
for (i in 1:200) {
  z1_samples_lps[[i]] <- as.vector(t(z1_subset_lps[i,]))
}

#check
z1_samples_lps[[1]]

#set up stan data list
stan_data_total_lps = list()
for (i in 1:200){
  n = ncol(post_z_lps) #subjects
  group = as.integer(as.factor(lps_df$id.group.x))#groups
  z = z_score(z1_samples_lps[[i]])
  sex = as.integer(lps_df$demos.sex) +1
  conchx = lps_df$conchx_bin
  int_term = as.integer(lps_df$int_term)
  stan_data_total_lps[[i]] <- list(n = n, 
                                   group = group,
                                   sex = sex,int_term = int_term,
                                   conchx = conchx,
                                   z = z)}

#check
stan_data_total_lps[[1]]

#Run the loop to generate 200 samples and run each of these 200 models for 3000 
#iterations
draws_total_lps <- list()
for (i in 1:200) {
  draws_total_lps[[i]] <- sampling(student_total_interact_pooled ,
                                   data = stan_data_total_lps[[i]],
                                   chains = 4, iter = 3000)}

post_draws_total_lps<- lapply(draws_total_lps,as.data.frame(extract.samples))
post_draws_total_lps_df <- do.call(rbind,post_draws_total_lps)
precis(post_draws_total_lps_df,2)

#there were scattered divergent transitions, in ~3/200 models, 1 per model
warnings()

#investigate the model performance
#pairs plot of parameters

#check trankplots from various draws
trankplot(draws_total_lps[[2]])

#r hat values 
rhat_df <- list()
for (i in 1:200){
  df_prep <- as.data.frame(summary(draws_total_lps[[i]])$summary)
  rhat_df[[i]]<- df_prep[c(9,10)]}

#pull all together
rhat_lps_df <- do.call(cbind,rhat_df)

#make random draws
draws <- sample(1:nrow(post_draws_total_lps_df),size = 1000)
colnames(post_draws_total_lps_df)

#take a random 1000 draws from the posterior, and just the parameter columns
post_pairs_lps <- post_draws_total_lps_df[draws,c(2:7,123:130)]
pairs(post_pairs_lps)


###R848 Stim
#separate by conditions and add some demos for downstream causal analysis
data_r848_df <- biomarkers_df[biomarkers_df$group=='R848',]

#find % missing data in a named vector called 'vec'
vec_r848 <- apply(data_r848_df[c(5:49)],2,function(x) (length(x[is.na(x)])/length(x))*100)
names_remove <- names(which(vec_r848>50))

#remove less than 50% cytokines
data_r848_df <- data_r848_df[!(colnames(data_r848_df)%in%names_remove)]

####Prep r848 latent model
r848_df <- merge(demos_df,data_r848_df, by = 'id.blood',all.y=TRUE)

#clean up
colnames(r848_df)

#separate the biomarkers to create a latent variable
bio_r848 <- r848_df[c(8:34)]

#which columns are missing
#find % missing data in a named vector called 'vec'
miss_r848 <- apply(bio_r848,2,function(x) (length(x[is.na(x)])/length(x))*100)
names_miss_r848 <- names(which(miss_r848>0))

bio_r848_miss <- bio_r848[colnames(bio_r848)%in%names_miss_r848]

#create bounds for the imputed parameters
bio_r848_means <- colMeans(bio_r848_miss,na.rm = TRUE)
bio_r848_sds <- apply(bio_r848_miss,2,sd,na.rm=TRUE)
#now create a lowerbound vector on the zscale
zeroes_r848 <- sapply(1:ncol(bio_r848_miss),function(i) (0 - bio_r848_means[i])/bio_r848_sds[i] )
cut_values_r848 <- sapply(1:ncol(bio_r848_miss),function(i) (min(bio_r848_miss[,i],na.rm = TRUE) -
                                                               bio_r848_means[i])/bio_r848_sds[i])

n_r848 = nrow(bio_r848)
mc_r848 = ncol(bio_r848)
n_miss_r848 <- length(bio_r848[is.na(bio_r848)])

#isolate missing from observed
bio_r848_z <-data.frame(apply(bio_r848,2,z_score))
na_values_r848 <- is.na(bio_r848_z)
y_obs_r848 <- bio_r848_z[!na_values_r848]
n_obs_r848 <- length(y_obs_r848)
row_miss_r848 <- rep(1:n_r848, times = mc_r848)[na_values_r848]
col_miss_r848<- rep(1:mc_r848, each = n_r848)[na_values_r848]
miss_check_r848 <- cbind(row_miss_r848,col_miss_r848) 

#identifies index where observed values are
row_obs_r848 <- rep(1:n_r848, times = mc_r848)[!na_values_r848]
col_obs_r848 <- rep(1:mc_r848, each = n_r848)[!na_values_r848]


bounds_r848 <- data.frame(Lo = zeroes_r848,Up = cut_values_r848)
bounds_r848$col_miss_r848 <- unique(col_miss_r848) 
col_miss_r8481 <- as.data.frame(col_miss_r848)
bounds1_r848 <- merge(bounds_r848,col_miss_r8481, by = 'col_miss_r848',all =TRUE)
lo_r848 <- as.vector(bounds1_r848$Lo)
up_r848 <- as.vector(bounds1_r848$Up)

dat_r848 = list(n = nrow(bio_r848),
                m = ncol(bio_r848),
                k = 1,
                row_obs = row_obs_r848,
                col_obs = col_obs_r848,
                row_miss = row_miss_r848,
                col_miss = col_miss_r848,
                n_miss = n_miss_r848,
                n_obs = n_obs_r848,
                y_obs = y_obs_r848,
                y = as.matrix(bio_r848_z),
                Lo = lo_r848,
                Up = up_r848)

str(dat_r848)

#run factor model with imputation
m_r848<- sampling(object = factor_impute_low, data = dat_r848, chains = 4, 
                  iter = 3000,init_r=.2 )

#extract posterior
m_post_r848 <- data.frame(extract.samples(m_r848))
precis(m_post_r848,3)

trankplot(m_r848, pars = c('z'))

#model prior
#create fake matrix 
y_prior_r848 <- matrix(0,nrow=n_r848,ncol=mc_r848)

dat_prior_r848 <-  list(
  k = 1,
  m = mc_r848,
  n = n_r848,
  n_obs = n_r848*mc_r848,
  y = y_prior_r848)

str(dat_prior_r848)

#run prior
bio_prior_r848<- sampling(object = mod_prior, data = dat_prior_r848, chains = 4, 
                          iter = 3000 ,init_r = .2) # add some init
bio_prior_r848

##########Set up plots
###Posterior
#Compare post factor loadings
#create a df for factor loadings
post_L_r848 <- m_post_r848[grepl("L",colnames(m_post_r848))&
                             !(grepl("Sigma",colnames(m_post_r848)))]
#sigma leftover, just remove.
post_L_r848$sigma_L <- NULL

#create separate dataframes, and then combine
f1_r848<- post_L_r848[c(1:mc_r848)]
colnames(f1_r848) <- colnames(bio_r848)

#transform to long
f1_r848_long <- gather(f1_r848,key = 'biomarker', value = 'Loading score')
#add factor #
f1_r848_long$factor <- 'Factor 1'
f1_r848_long$biomarker<- factor(f1_r848_long$biomarker,
                                levels = c('IFNG','IL1B','IL2','IL4','IL6','IL7','IL10','IL13','IL15','IL17A','IL17C','IL17F',
                                           'IL18','IL27','IL33','TNF','TNFSF10','TNFSF12','FLT3LG','LTA','CSF1','CSF2','CSF3','OSM','TSLP',
                                           'CCL2','CCL3','CCL4','CCL7','CCL8','CCL11','CCL13','CCL19','CXCL8','CXCL9',
                                           'CXCL10','CXCL11','CXCL12','HGF','MMP1','MMP12','OLR1','EGF','TGFA','VEGFA'))

#prior
#extract the prior
prior_r848 <- data.frame(extract.samples(bio_prior_r848))
precis(prior_r848,2)

#Compare prior factor loadings
#create a df for factor loadings
prior_L_r848 <- prior_r848[grepl("L",colnames(prior_r848))&
                             !(grepl("Sigma",colnames(prior_r848)))]
#sigma leftover, just remove.
prior_L_r848$sigma_L <- NULL

#create separate dataframes, and then combine
f1_r848_prior<- prior_L_r848[c(1:mc_r848)]
colnames(f1_r848_prior) <- colnames(bio_r848)

#transform to long
f1_r848_long_prior <- gather(f1_r848_prior,key = 'biomarker', value = 'Loading score')
#add factor #
f1_r848_long_prior$factor <- 'Factor 1'

f1_r848_long_prior$biomarker <- factor(f1_r848_long_prior$biomarker,
                                       levels = c('IFNG','IL1B','IL2','IL4','IL6','IL7','IL10','IL13','IL15','IL17A','IL17C','IL17F',
                                                  'IL18','IL27','IL33','TNF','TNFSF10','TNFSF12','FLT3LG','LTA','CSF1','CSF2','CSF3','OSM','TSLP',
                                                  'CCL2','CCL3','CCL4','CCL7','CCL8','CCL11','CCL13','CCL19','CXCL8','CXCL9',
                                                  'CXCL10','CXCL11','CXCL12','HGF','MMP1','MMP12','OLR1','EGF','TGFA','VEGFA'))
#combine
f1_r848_long$samples <- 'Posterior'
f1_r848_long_prior$samples <-'Prior'

plot_combined_r848 <- rbind.data.frame(f1_r848_long,f1_r848_long_prior)

#use this if interpretation is not useful
plot_combined_r848$group <- 'Cytokines'
plot_combined_r848$group[grepl('CC',plot_combined_r848$biomarker)|
                           grepl('CXC',plot_combined_r848$biomarker)] <- 'Chemokines'
plot_combined_r848$group[plot_combined_r848$biomarker=='MMP12'|
                           plot_combined_r848$biomarker=='MMP1'|
                           plot_combined_r848$biomarker=='OLR1'|
                           plot_combined_r848$biomarker=='OSM'|
                           plot_combined_r848$biomarker=='TGFA'|
                           plot_combined_r848$biomarker=='EGF'|
                           plot_combined_r848$biomarker=='HGF'|
                           plot_combined_r848$biomarker=='VEGFA'] <- 'Other'

#combined plot
theme_set(theme_tidybayes() + panel_border())
plot_r848<- plot_combined_r848 %>%
  ggplot(aes(y = fct_rev(`biomarker`),
             x = `Loading score`,
             fill = samples)) + 
  geom_vline(xintercept = 0, alpha = 0.8, linetype = 2) + 
  theme(title = element_text(face = 'bold'),
        axis.title.x = element_text(size = 13, face = 'bold'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = 'bold', size = 10),
        legend.position = 'None',
        axis.text.x = element_text(face = 'bold', size = 10)) +
  ggtitle("R848 Stimulation")+xlim(-2,2)+
  scale_fill_manual(values = c('green','grey'))+
  stat_halfeye(aes(color = samples),.width = c(0.70,0.90), 
               alpha = 0.5,linewidth = 1.5,point_size = 1) +
  scale_colour_manual(values = c('green','grey'))

plot_r848

ggsave('plot_r848.jpg',plot_r848,dpi = 600,width = 5,height = 10)


#########Biomarker covariance structure

#need to convert the covariances (Sigma) by calculating first the sd of the diagonal...
#error terms for each biomarker (sigma_y)
#first pull the relevant values from the posterior
df_Sig_r848 <- m_post_r848[grepl('Sigma', colnames(m_post_r848))]

#remove the "L"
df_Sig_r848 <- df_Sig_r848[!(grepl('L',colnames(df_Sig_r848)))]

#derive the means
vec_Sig_r848_means<- apply(df_Sig_r848,2,mean)

#create a covariance matrix
cov_mat_r848 <- matrix(vec_Sig_r848_means,nrow = 27, ncol = 27)

#rename
colnames(cov_mat_r848) <- colnames(bio_r848)
rownames(cov_mat_r848) <- colnames(bio_r848)

# Convert to correlation matrix
cor_mat_r848 <- cov2cor(cov_mat_r848)

#plot
corrplot(cor_mat_r848)

#how does this compare to the regular data structure?
corrplot(cor(bio_r848_z,use = "pairwise.complete.obs"))
#by group
corrplot(cor(bio_r848_z[r848_df$id.group.x=='healthy',],use = "pairwise.complete.obs"))

######DOWNSTREAM MODELLING of the latent variable scores created.

#TOTAL EFFECTS

#first create some necessary variables 
#create binary conchx
r848_df$conchx_bin <- as.integer(ifelse(r848_df$demos.hx_conc_number==0,1,2))

#create interaction term
r848_df$int_term <- 'healthy_male_no'   #1
r848_df$int_term[r848_df$id.group.x=='healthy'&
                    r848_df$demos.sex==0&
                    r848_df$conchx_bin==2] <-'healthy_male_yes' #2
r848_df$int_term[r848_df$id.group.x=='healthy'&
                    r848_df$demos.sex==1&
                    r848_df$conchx_bin==1] <-'healthy_female_no' #3
r848_df$int_term[r848_df$id.group.x=='healthy'&
                    r848_df$demos.sex==1&
                    r848_df$conchx_bin==2] <-'healthy_female_yes'#4
r848_df$int_term[r848_df$id.group.x=='concussion'&
                    r848_df$demos.sex==0&
                    r848_df$conchx_bin==1] <-'src_male_no' #5
r848_df$int_term[r848_df$id.group.x=='concussion'&
                    r848_df$demos.sex==0&
                    r848_df$conchx_bin==2] <-'src_male_yes' #6
r848_df$int_term[r848_df$id.group.x=='concussion'&
                    r848_df$demos.sex==1&
                    r848_df$conchx_bin==1] <-'src_female_no' #7
r848_df$int_term[r848_df$id.group.x=='concussion'&
                    r848_df$demos.sex==1&
                    r848_df$conchx_bin==2] <-'src_female_yes' #8
table(r848_df$int_term)
r848_df$int_term <- factor(r848_df$int_term,levels = c('healthy_male_no',
                                                         'healthy_male_yes','healthy_female_no','healthy_female_yes',
                                                         'src_male_no','src_male_yes','src_female_no','src_female_yes'))


#using mean modelling to 1) compare unpooled and pooled estimates, and 2) derive
#the prior predictive checks.

#isolate z scores
colnames(m_post_r848)
post_z_r848 <- m_post_r848[grepl("z",colnames(m_post_r848))]

#derive the data list
dat_mean_total_r848 = list(
  z = z_score(as.numeric(apply(post_z_r848,2,mean))),
  conchx = r848_df$conchx_bin,
  group = as.integer(as.factor(r848_df$id.group.x)),
  sex = as.integer(r848_df$demos.sex) +1,
  int_term = as.integer(r848_df$int_term),
  n = ncol(post_z_r848)
)

str(dat_mean_total_r848)


#run both mean models
m_stip2 <- sampling(student_total_interact_pooled,dat = dat_mean_total_r848,
                    chains = 4, iter = 1000)
m_stiu2 <- sampling(student_total_interact_unpooled,dat = dat_mean_total_r848,
                    chains = 4, iter = 1000)

#evaluate models out of bag
loo_stip2 <- loo(m_stip2)
loo_stiu2 <- loo(m_stiu2)

#compare
loo_compare(loo_stip2,loo_stiu2)
rethinking::compare(m_stip2,m_stiu2,func = PSIS)

#extract posteriors
post_stip2 <- data.frame(extract.samples(m_stip2))
post_stiu2 <- data.frame(extract.samples(m_stiu2))

precis(post_stip2,2,prob = 0.9)
precis(post_stiu2,2,prob = 0.9)


#model the total effects from DAG1 with the interaction term pooled
#full model with 100 draws from the latent model posterior
n_subsets = 200
post_size = 6000
subset_rows <- sample(post_size,n_subsets,replace = F)
z1_subset_r848 <- post_z_r848[subset_rows,]

z1_samples_r848 <- list()
for (i in 1:200) {
  z1_samples_r848[[i]] <- as.vector(t(z1_subset_r848[i,]))
}

#check
z1_samples_r848[[1]]

#set up stan data list
stan_data_total_r848 = list()
for (i in 1:200){
  n = ncol(post_z_r848) #subjects
  group = as.integer(as.factor(r848_df$id.group.x))#groups
  z = z_score(z1_samples_r848[[i]])
  sex = as.integer(r848_df$demos.sex) +1
  conchx = r848_df$conchx_bin
  int_term = as.integer(r848_df$int_term)
  stan_data_total_r848[[i]] <- list(n = n, 
                                    group = group,
                                    sex = sex,int_term = int_term,
                                    conchx = conchx,
                                    z = z)}

#check
stan_data_total_r848[[1]]


#Run the loop to generate 200 samples and run each of these 200 models for 3000 
#iterations
draws_total_r848 <- list()
for (i in 1:200) {
  draws_total_r848[[i]] <- sampling(student_total_interact_pooled ,
                                    data = stan_data_total_r848[[i]],
                                    chains = 4, iter = 3000)}

post_draws_total_r848<- lapply(draws_total_r848,as.data.frame(extract.samples))
post_draws_total_r848_df <- do.call(rbind,post_draws_total_r848)
precis(post_draws_total_r848_df,2)

#there were scattered divergent transitions, in ~3/200 models, 1 per model
warnings()

#investigate the model performance

#check trankplots from various draws
trankplot(draws_total_r848[[2]])

#r hat values 
rhat_df <- list()
for (i in 1:200){
  df_prep <- as.data.frame(summary(draws_total_r848[[i]])$summary)
  rhat_df[[i]]<- df_prep[c(9,10)]}

#pull all together
rhat_r848_df <- do.call(cbind,rhat_df)

#Pairs plot of random subs amples

#make random draws
draws <- sample(1:nrow(post_draws_total_r848_df),size = 1000)
colnames(post_draws_total_r848_df)

#take a random 1000 draws from the posterior, and just the parameter columns
post_pairs_r848 <- post_draws_total_r848_df[draws,c(2:7,123:130)]
pairs(post_pairs_r848)


####Plotting and Results reporting

#to plot, first we are going to take a few thousand draws from the amalgamated
#posteriors of both models

#generate samples
#lps
lps_sample <- sample(1:nrow(post_draws_total_lps_df),2000)
post_lps_df <- post_draws_total_lps_df[lps_sample,c(1:9,18,123:130)]

#r848
r848_sample <- sample(1:nrow(post_draws_total_r848_df),2000)
post_r848_df <- post_draws_total_r848_df[r848_sample,c(1:9,18,123:130)]

#add a stimulation variable to both and then bind

post_lps_df$stim <-'lps'
post_r848_df$stim <- 'r848'

post_plot_df <- rbind.data.frame(post_lps_df,post_r848_df)

#name columns
colnames(post_plot_df) <- c('nu','src','healthy','no_chx','chx','male','female',
                            'd_mu','d_sig','model_sig','healthy_male_no','healthy_male_yes',
                            'healthy_female_no','healthy_female_yes',
                            'src_male_no','src_male_yes','src_female_no',
                            'src_female_yes','stim')

#add contrasts
#male w conchx
post_plot_df$src_m_y <- post_plot_df$src_male_yes + post_plot_df$male +
  post_plot_df$src + post_plot_df$chx

post_plot_df$hlth_m_y <- post_plot_df$healthy_male_yes + post_plot_df$male +
  post_plot_df$healthy + post_plot_df$chx

#contrast
post_plot_df$contrast_male_yes <- post_plot_df$src_m_y - post_plot_df$hlth_m_y

#male w no conchx
post_plot_df$src_m_n <- post_plot_df$src_male_no + post_plot_df$male +
  post_plot_df$src + post_plot_df$no_chx

post_plot_df$hlth_m_n <- post_plot_df$healthy_male_no + post_plot_df$male +
  post_plot_df$healthy + post_plot_df$no_chx

#contrast
post_plot_df$contrast_male_no <- post_plot_df$src_m_n - post_plot_df$hlth_m_n

#female w conchx
post_plot_df$src_f_y <- post_plot_df$src_female_yes + post_plot_df$female +
  post_plot_df$src + post_plot_df$chx

post_plot_df$hlth_f_y <- post_plot_df$healthy_female_yes + post_plot_df$female +
  post_plot_df$healthy + post_plot_df$chx

#contrast
post_plot_df$contrast_female_yes <- post_plot_df$src_f_y - post_plot_df$hlth_f_y

#female w no conchx
post_plot_df$src_f_n <- post_plot_df$src_female_no + post_plot_df$female +
  post_plot_df$src + post_plot_df$no_chx

post_plot_df$hlth_f_n <- post_plot_df$healthy_female_no + post_plot_df$female +
  post_plot_df$healthy + post_plot_df$no_chx

#contrast
post_plot_df$contrast_female_no <- post_plot_df$src_f_n - post_plot_df$hlth_f_n

#healthy conchx contrasts
#male
post_plot_df$contrast_chx_hlthy_m <- post_plot_df$hlth_m_y - post_plot_df$hlth_m_n

#female
post_plot_df$contrast_chx_hlthy_f <- post_plot_df$hlth_f_y - post_plot_df$hlth_f_n

#have a look (also for results)
precis(post_plot_df[post_plot_df$stim=='r848',],2,probs = 0.9)

#create prior draws to overlay on graphs

#extract from m_stip_prior
post_stip_prior <- data.frame(extract.samples(m_stip_prior))

#match posterior columns
colnames(post_stip_prior)
post_stip_prior <- post_stip_prior[c(1:9,18:26)]
#name columns
colnames(post_stip_prior) <- c('nu','src','healthy','no_chx','chx','male','female',
                               'd_mu','d_sig','model_sig','healthy_male_no','healthy_male_yes',
                               'healthy_female_no','healthy_female_yes',
                               'src_male_no','src_male_yes','src_female_no',
                               'src_female_yes')

#add contrasts
#male w conchx
post_stip_prior$src_m_y <- post_stip_prior$src_male_yes + post_stip_prior$male +
  post_stip_prior$src + post_stip_prior$chx

post_stip_prior$hlth_m_y <- post_stip_prior$healthy_male_yes + post_stip_prior$male +
  post_stip_prior$healthy + post_stip_prior$chx

#contrast
post_stip_prior$contrast_male_yes <- post_stip_prior$src_m_y - post_stip_prior$hlth_m_y

#male w no conchx
post_stip_prior$src_m_n <- post_stip_prior$src_male_no + post_stip_prior$male +
  post_stip_prior$src + post_stip_prior$no_chx

post_stip_prior$hlth_m_n <- post_stip_prior$healthy_male_no + post_stip_prior$male +
  post_stip_prior$healthy + post_stip_prior$no_chx

#contrast
post_stip_prior$contrast_male_no <- post_stip_prior$src_m_n - post_stip_prior$hlth_m_n

#female w conchx
post_stip_prior$src_f_y <- post_stip_prior$src_female_yes + post_stip_prior$female +
  post_stip_prior$src + post_stip_prior$chx

post_stip_prior$hlth_f_y <- post_stip_prior$healthy_female_yes + post_stip_prior$female +
  post_stip_prior$healthy + post_stip_prior$chx

#contrast
post_stip_prior$contrast_female_yes <- post_stip_prior$src_f_y - post_stip_prior$hlth_f_y

#female w no conchx
post_stip_prior$src_f_n <- post_stip_prior$src_female_no + post_stip_prior$female +
  post_stip_prior$src + post_stip_prior$no_chx

post_stip_prior$hlth_f_n <- post_stip_prior$healthy_female_no + post_stip_prior$female +
  post_stip_prior$healthy + post_stip_prior$no_chx

#contrast
post_stip_prior$contrast_female_no <- post_stip_prior$src_f_n - post_stip_prior$hlth_f_n

#healthy conchx contrasts
#male
post_stip_prior$contrast_chx_hlthy_m <- post_stip_prior$hlth_m_y - post_stip_prior$hlth_m_n

#female
post_stip_prior$contrast_chx_hlthy_f <- post_stip_prior$hlth_f_y - post_stip_prior$hlth_f_n

precis(post_stip_prior,probs = 0.9)

###############PLOTS################
#LPS Plot

#save first
pdf('fig1.pdf',width = 12, height = 8, pointsize = 12)
#set multiplot
par(mfrow=c(2,4))
#males no concussion history
dens(x=post_plot_df$src_m_n[post_plot_df$stim=='lps'],ylim = c(0,1.9),xlim=c(-2,2),
     lwd = 3,col = 2, xlab="",
     font.lab=2, cex.lab = 0.9, main = 'Males: No Conc Hx.',
     frame = FALSE)
dens(post_plot_df$hlth_m_n[post_plot_df$stim=='lps'],lwd = 3,col = 4,add = TRUE)
dens(post_stip_prior$src_m_n,lwd = 1,lty = 6,col = 'gray40', add = TRUE)
dens(post_stip_prior$hlth_m_n,lwd = 1,lty = 6,col = 'gray80', add = TRUE)
#legend(-2,1.5, box.lty = 0, cex = 1,
# legend = c("SRC",'Healthy','SRC Prior','Healthy Prior'),
# fill = c("red",'blue','gray40','gray80'))
legend(-2,1.5, box.lty = 0, cex = 1.1,
       legend = c("SRC",'Healthy'),
       fill = c("red",'blue'))
text(-1.7,0.25,'SRC Prior',cex = 0.7,col = 'gray40')
text(-1.6,0.17,'Healthy Prior',cex = 0.7,col = 'gray80')
text(-2,1.5,'A',cex = 0.9,font = 2)

#males concussion history
dens(x=post_plot_df$src_m_y[post_plot_df$stim=='lps'],ylim = c(0,1.9),xlim=c(-2,2),
     lwd = 3,col = 2,xlab = "", ylab = "",
     font.lab=2, cex.lab = 0.8, main = 'Males: Conc Hx.',
     frame = FALSE)
dens(post_plot_df$hlth_m_y[post_plot_df$stim=='lps'],lwd = 3,col = 4,add = TRUE)
dens(post_stip_prior$src_m_y,lwd = 1,lty = 6,col = 'gray40', add = TRUE)
dens(post_stip_prior$hlth_m_y,lwd = 1,lty = 6,col = 'gray80', add = TRUE)
text(-2,1.5,'B',cex = 0.9,font = 2)

#females no concussion history
dens(x=post_plot_df$src_f_n[post_plot_df$stim=='lps'],ylim = c(0,1.9),xlim=c(-2,2),
     lwd = 3,col = 2,main = 'Females: No Conc Hx.',
     font.lab=2, cex.lab = 0.8,xlab = "", ylab = "",
     frame = FALSE)
dens(post_plot_df$hlth_f_n[post_plot_df$stim=='lps'],lwd = 3,col = 4,add = TRUE)
dens(post_stip_prior$src_f_n,lwd = 1,lty = 6,col = 'gray40', add = TRUE)
dens(post_stip_prior$hlth_f_n,lwd = 1,lty = 6,col = 'gray80', add = TRUE)
text(-2,1.5,'C',cex = 0.9,font = 2)
mtext('LPS Reactivity Score',side = 3, line = -29, outer = TRUE,
      font = 2)

#females concussion history
dens(x=post_plot_df$src_f_y[post_plot_df$stim=='lps'],ylim = c(0,1.9),xlim=c(-2,2),
     lwd = 3,col = 2,main = 'Females: Conc Hx.',
     font.lab=2, cex.lab = 0.8,xlab = "", ylab = "",
     frame = FALSE)
dens(post_plot_df$hlth_f_y[post_plot_df$stim=='lps'],lwd = 3,col = 4,add = TRUE)
dens(post_stip_prior$src_f_y,lwd = 1,lty = 6,col = 'gray40', add = TRUE)
dens(post_stip_prior$hlth_f_y,lwd = 1,lty = 6,col = 'gray80', add = TRUE)
text(-2,1.5,'D',cex = 0.9,font = 2)

###CONTRASTS
#contrast male no conchx
dens( post_plot_df$contrast_male_no[post_plot_df$stim=='lps'] , xlim=c(-3,2.5) , 
      ylim=c(0,1.4),lwd=3 , col=1 ,
      cex.lab = 0.9,font.lab = 2,xlab = "",frame = FALSE)
text(-3,1.4,'E',cex = 0.9, font = 2)
mnc <- density(post_plot_df$contrast_male_no[post_plot_df$stim=='lps'], adjust = 0.5)
polygon(c(mnc$x[mnc$x>0], max(mnc$x), 0),
        c(mnc$y[mnc$x>0], 0, 0), col = 4, border = NA )
polygon(c(mnc$x[mnc$x<0], 0, min(mnc$x)), 
        c(mnc$y[mnc$x<0], 0, 0), col = 2, border = NA )

#contrast male conchx
dens( post_plot_df$contrast_male_yes[post_plot_df$stim=='lps'] , xlim=c(-3,2.5) ,
      ylim=c(0,1.4),lwd=3 , col=1 ,
      cex.lab = 0.9 ,font.lab = 2,xlab = "", ylab = "",frame = FALSE)
text(-3,1.4,'F',cex = 0.9, font = 2)
myc <- density(post_plot_df$contrast_male_yes[post_plot_df$stim=='lps'], adjust = 0.5)
polygon(c(myc$x[myc$x>0], max(myc$x), 0),
        c(myc$y[myc$x>0], 0, 0), col = 2, border = NA )
polygon(c(myc$x[myc$x<0], 0, min(myc$x)), 
        c(myc$y[myc$x<0], 0, 0), col = 4, border = NA )

#contrast female no conchx
dens( post_plot_df$contrast_female_no[post_plot_df$stim=='lps'] , xlim=c(-3,2.5) , 
      ylim=c(0,1.4),lwd=3 , col=1 ,
      cex.lab = 0.9,xlab = '', ylab = "",
      font.lab = 2,frame = FALSE)
text(-3,1.4,'G',cex = 0.9, font = 2)
fnc <- density(post_plot_df$contrast_female_no[post_plot_df$stim=='lps'], adjust = 0.5)
polygon(c(fnc$x[fnc$x>0], max(fnc$x), 0),
        c(fnc$y[fnc$x>0], 0, 0), col = 2, border = NA )
polygon(c(fnc$x[fnc$x<0], 0, min(fnc$x)), 
        c(fnc$y[fnc$x<0], 0, 0), col = 4, border = NA )

#contrast female conchx
dens( post_plot_df$contrast_female_yes[post_plot_df$stim=='lps'] , xlim=c(-3,2.5) , 
      ylim=c(0,1.4),lwd=3 , col=1 ,
      cex.lab = 0.9,xlab = '', ylab = "",
      font.lab = 2,frame = FALSE)
text(-3,1.4,'H',cex = 0.9, font = 2)
mtext('Contrast (SRC - Healthy)',side = 3, line = -59, outer = TRUE,
      font = 2)
fyc <- density(post_plot_df$contrast_female_yes[post_plot_df$stim=='lps'], adjust = 0.5)
polygon(c(fyc$x[fyc$x>0], max(fyc$x), 0),
        c(fyc$y[fyc$x>0], 0, 0), col = 4, border = NA )
polygon(c(fyc$x[fyc$x<0], 0, min(fyc$x)), 
        c(fyc$y[fyc$x<0], 0, 0), col = 2, border = NA )

fig1 <- recordPlot()
dev.off() 

#R848 PLOT
#save first
pdf('fig2.pdf',width = 12, height = 8, pointsize = 12)
#set multiplot
par(mfrow=c(2,4))
#males no concussion history
dens(x=post_plot_df$src_m_n[post_plot_df$stim=='r848'],ylim = c(0,1.9),xlim=c(-2,2),
     lwd = 3,col = 3, xlab="",
     font.lab=2, cex.lab = 0.9, main = 'Males: No Conc Hx.',
     frame = FALSE)
dens(post_plot_df$hlth_m_n[post_plot_df$stim=='r848'],lwd = 3,col = 4,add = TRUE)
dens(post_stip_prior$src_m_n,lwd = 1,lty = 6,col = 'gray40', add = TRUE)
dens(post_stip_prior$hlth_m_n,lwd = 1,lty = 6,col = 'gray80', add = TRUE)
legend(-2.2,1.5, box.lty = 0, cex = 1,
       legend = c("SRC",'Healthy'),
       fill = c("green",'blue'))
text(-1.7,0.25,'SRC Prior',cex = 0.7,col = 'gray40')
text(-1.6,0.17,'Healthy Prior',cex = 0.7,col = 'gray80')
text(-2,1.5,'A',cex = 0.9,font = 2)

#males concussion history
dens(x=post_plot_df$src_m_y[post_plot_df$stim=='r848'],ylim = c(0,1.9),xlim=c(-2,2),
     lwd = 3,col = 3,xlab = "", ylab = "",
     font.lab=2, cex.lab = 0.8, main = 'Males: Conc Hx.',
     frame = FALSE)
dens(post_plot_df$hlth_m_y[post_plot_df$stim=='r848'],lwd = 3,col=4,add = TRUE)
dens(post_stip_prior$src_m_y,lwd = 1,lty = 6,col = 'gray40', add = TRUE)
dens(post_stip_prior$hlth_m_y,lwd = 1,lty = 6,col = 'gray80', add = TRUE)
text(-2,1.5,'B',cex = 0.9,font = 2)

#females no concussion history
dens(x=post_plot_df$src_f_n[post_plot_df$stim=='r848'],ylim = c(0,1.9),xlim=c(-2,2),
     lwd = 3,col = 3,main = 'Females: No Conc Hx.',
     font.lab=2, cex.lab = 0.8,xlab = "", ylab = "",
     frame = FALSE)
dens(post_plot_df$hlth_f_n[post_plot_df$stim=='r848'],lwd = 3,col = 4,add = TRUE)
dens(post_stip_prior$src_f_n,lwd = 1,lty = 6,col = 'gray40', add = TRUE)
dens(post_stip_prior$hlth_f_n,lwd = 1,lty = 6,col = 'gray80', add = TRUE)
text(-2,1.5,'C',cex = 0.9,font = 2)
mtext('R848 Reactivity Score',side = 3, line = -29, outer = TRUE,
      font = 2)

#females concussion history
dens(x=post_plot_df$src_f_y[post_plot_df$stim=='r848'],ylim = c(0,1.9),xlim=c(-2,2),
     lwd = 3,col = 3,main = 'Females: Conc Hx.',
     font.lab=2, cex.lab = 0.8,xlab = "", ylab = "",
     frame = FALSE)
dens(post_plot_df$hlth_f_y[post_plot_df$stim=='r848'],lwd = 3,col = 4,add = TRUE)
dens(post_stip_prior$src_f_y,lwd = 1,lty = 6,col = 'gray40', add = TRUE)
dens(post_stip_prior$hlth_f_y,lwd = 1,lty = 6,col = 'gray80', add = TRUE)
text(-2,1.5,'D',cex = 0.9,font = 2)

###CONTRASTS
#contrast male no conchx
dens( post_plot_df$contrast_male_no[post_plot_df$stim=='r848'] , xlim=c(-3,2.5) , 
      ylim=c(0,1.5),lwd=3 , col=1 ,
      cex.lab = 0.9,font.lab = 2,xlab = "",frame = FALSE)
text(-3,1.5,'E',cex = 0.9, font = 2)
mnc <- density(post_plot_df$contrast_male_no[post_plot_df$stim=='r848'], adjust = 0.5)
polygon(c(mnc$x[mnc$x>0], max(mnc$x), 0),
        c(mnc$y[mnc$x>0], 0, 0), col = 4, border = NA )
polygon(c(mnc$x[mnc$x<0], 0, min(mnc$x)), 
        c(mnc$y[mnc$x<0], 0, 0), col = 3, border = NA )

#contrast male conchx
dens( post_plot_df$contrast_male_yes[post_plot_df$stim=='r848'] , xlim=c(-3,2.5) ,
      ylim=c(0,1.5),lwd=3 , col=1 ,
      cex.lab = 0.8 ,font.lab = 2,xlab = "", ylab = "",frame = FALSE)
text(-3,1.5,'F',cex = 0.9, font = 2)
myc <- density(post_plot_df$contrast_male_yes[post_plot_df$stim=='r848'], adjust = 0.5)
polygon(c(myc$x[myc$x>0], max(myc$x), 0),
        c(myc$y[myc$x>0], 0, 0), col = 3, border = NA )
polygon(c(myc$x[myc$x<0], 0, min(myc$x)), 
        c(myc$y[myc$x<0], 0, 0), col = 4, border = NA )

#contrast female no conchx
dens( post_plot_df$contrast_female_no[post_plot_df$stim=='r848'] , xlim=c(-3,2.5) , 
      ylim=c(0,1.5),lwd=3 , col=1 ,
      cex.lab = 0.8,xlab = '', ylab = "",
      font.lab = 2,frame = FALSE)
text(-3,1.5,'G',cex = 0.9, font = 2)
fnc <- density(post_plot_df$contrast_female_no[post_plot_df$stim=='r848'], adjust = 0.5)
polygon(c(fnc$x[fnc$x>0], max(fnc$x), 0),
        c(fnc$y[fnc$x>0], 0, 0), col = 3, border = NA )
polygon(c(fnc$x[fnc$x<0], 0, min(fnc$x)), 
        c(fnc$y[fnc$x<0], 0, 0), col = 4, border = NA )

#contrast female conchx
dens( post_plot_df$contrast_female_yes[post_plot_df$stim=='r848'] , xlim=c(-3,2.5) , 
      ylim=c(0,1.5),lwd=3 , col=1 ,
      cex.lab = 0.8,xlab = '', ylab = "",
      font.lab = 2,frame = FALSE)
text(-3,1.5,'H',cex = 0.9, font = 2)
mtext('Contrast (SRC - Healthy)',side = 3, line = -59, outer = TRUE,
      font = 2)
fyc <- density(post_plot_df$contrast_female_yes[post_plot_df$stim=='r848'], adjust = 0.5)
polygon(c(fyc$x[fyc$x>0], max(fyc$x), 0),
        c(fyc$y[fyc$x>0], 0, 0), col = 4, border = NA )
polygon(c(fyc$x[fyc$x<0], 0, min(fyc$x)), 
        c(fyc$y[fyc$x<0], 0, 0), col = 3, border = NA )

fig2 <- recordPlot()
dev.off() 


#For results reporting, a precis output of contrasts
#lps
precis(post_plot_df[grepl('contrast',
                          colnames(post_plot_df))][post_plot_df$stim=='lps',],2,probs = 0.9)

#r848
precis(post_plot_df[grepl('contrast',
                          colnames(post_plot_df))][post_plot_df$stim=='r848',],2,probs = 0.9)

#pprobs
apply(post_plot_df[grepl('contrast',colnames(post_plot_df))]
      [post_plot_df$stim=='lps',],2,function(x)mean(x<0))

apply(post_plot_df[grepl('contrast',colnames(post_plot_df))]
      [post_plot_df$stim=='r848',],2,function(x)mean(x<0))

#effective sample size for interaction term
table(lps_df$int_term)

colnames(post_plot_df)

#sex differences for discussion
post_plot_df$contrast_sex <- (post_plot_df$healthy + post_plot_df$male) -
  (post_plot_df$healthy + post_plot_df$female)
precis(post_plot_df[grepl('contrast',
                          colnames(post_plot_df))][post_plot_df$stim=='lps',],2,probs = 0.9)
precis(post_plot_df[grepl('contrast',
                          colnames(post_plot_df))][post_plot_df$stim=='r848',],2,probs = 0.9)
mean(post_plot_df$contrast_sex[post_plot_df$stim=='lps']>0)
mean(post_plot_df$contrast_sex[post_plot_df$stim=='r848']>0)



