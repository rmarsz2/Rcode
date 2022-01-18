library(ggplot2)#library for graphics
library(readxl)#library for excel files
library(RColorBrewer)#library for nice colors in graphics
library(broom)#to deal with ugly data
library(robustbase)# for robust non-linear regression
library(plotly)
library(nlstools)
library(ggpubr)
library(plyr)
library(dplyr)#to manipulate data more easily

##### <- comment symbol
#Create portable file path for Windows or Mac or Linux
FILE_SEP <<- .Platform$file.sep
HOME_DIR <<- normalizePath(path = Sys.getenv('HOMEPATH'))
#open file
filename = normalizePath(file.path(HOME_DIR, "Box",
                                   "Richard",
                                   "Data",
                                   "Fiber Data",
                                   "2020-03-13 all mice from monika.xlsx"))

fiber_data1 <- read_excel(filename) #extract data

x_axis <- 2:(length(fiber_data1)-2) #number of pCa concentrations = number of columns minus "two"
number_fibers <- nrow(fiber_data1) #number of rows

#create large Data Frame; requires a precise format for grpaphing
#t() function transposes, colnames() retrieves column names, as.numeric() converts data to numerical format
data_proc <- data.frame("pCa" = as.numeric(t(colnames(fiber_data1)[x_axis])),
                        "Tension" = as.numeric(t(fiber_data1[,x_axis])),
                        "factor" = rep(t(fiber_data1[,(length(fiber_data1)-1)]), each = length(x_axis)),
                        "factor2" = rep(t(fiber_data1[,length(fiber_data1)]), each = length(x_axis)),
                        "fiber" = rep(t(fiber_data1[,1]), each = length(x_axis)))

data_proc$factor <- as.factor(data_proc$factor) #enforce that factor column is composed of factors
data_proc$factor2 <- as.factor(data_proc$factor2) #enforce that factor column is composed of factors
l_factors <- list("factor" = levels(data_proc$factor), "factor2" = levels(data_proc$factor2))
isFactor <- function(x) which(as.logical((x[1] == data_proc$factor)*(x[2] == data_proc$factor2))) #function that returns which data are of a certin type, or factor
combinations <- list(rep(l_factors$factor, length(l_factors$factor2)), rep(l_factors$factor2, each = length(l_factors$factor)))
t_comb <- lapply(1:6, function(x) c(combinations[[1]][x], combinations[[2]][x]))
names(t_comb) <- lapply(1:6, function(x) paste(unlist(t_comb[x]), collapse = " "))

##generate main function
force_Ca_function <- Tension ~ force_Max/(1+10^(nH*(pCa-pCa50)))

##nonlinear regression
##force_Ca_function <- fiber_data1[6] ~ force_Max*(1/(1+10^(nH*(fiber_data1[1]-pCa50))))
##for each different factor, or fiber type, a separate regression is performed
analysis_nls <- lapply(t_comb, function(x) nls(force_Ca_function, data = data_proc[isFactor(x),],  start = list(nH = 3, pCa50 = 6, force_Max = 42)))
names(analysis_nls) <- names(t_comb)
single_fiber_nls <- lapply(1:number_fibers, function(x) nls(force_Ca_function, data = data_proc[((8*x-7):(8*x)),],  start = list(nH = 3, pCa50 = 6, force_Max = 42)))
names(single_fiber_nls) <- t(fiber_data1[,1])

#get bestFit curve
parameters <- lapply(names(t_comb), function(x) tidy(analysis_nls[[x]]))
names(parameters) <- names(t_comb)
single_fiber_parameters <- lapply(1:number_fibers, function(x) single_fiber_nls[[x]]$m$getPars())
names(single_fiber_parameters) <- names(single_fiber_nls)
bestFit <- function(x, i) {parameters[[i]]$estimate[3]*(1/(1+10^(parameters[[i]]$estimate[1]*(x - parameters[[i]]$estimate[2]))))}
avg_pCa50 <- mean(unlist(lapply(1:number_fibers, function(x) single_fiber_parameters[[x]][2])))

##two-way anova
df_params <-as.data.frame(do.call(rbind, single_fiber_parameters))
df_params$ID <- rownames(df_params)
df_params <- merge(df_params, fiber_data1[,c('ID','genotype','age')], by.x = 'ID')
#perform two-way anova on all variables
variables <- names(df_params)[2:4] #returns variables interested in being analyzed
#apply ANOVA over all variables
list_aov <- lapply(variables, function(x) {
  f <- as.formula(paste(df_params[x],"~genotype*age")) #define function with two independant variables (two-way)
  aov(f, data = df_params)})
names(list_aov) <- variables
table_aovs <- lapply(variables, function(x) {
  as_tibble(data.frame(
    "p" = summary(list_aov[[x]])[[1]][["Pr(>F)"]][1:3]))
})
table_aovs <- t(Reduce(function(x, y) cbind(x, y), table_aovs))

rownames(table_aovs) <- variables
colnames(table_aovs) <- rownames(summary(list_aov$nH)[[1]])[1:3]

st.err <- function(x)
  sd(x)/sqrt(length(x))

df_params$genotype <- factor(df_params$genotype, c("NTG", "HCM"))
df_params$age <- factor(df_params$age, c(7,14,28))

list_df_sums <- lapply(variables, function(x) {
  f <- as.formula(paste(df_params[x],"~genotype*age")) #formula with two variables
  y <- aggregate(f, df_params, mean) #compute summary statistics
  y <- merge(y, aggregate(f, df_params, sd), by = c("age", "genotype"))
  y <- merge(y, aggregate(f, df_params, st.err), by = c("age", "genotype"))
  y <- merge(y, aggregate(f, df_params, length), by = c("age", "genotype"))
  names(y) <- c("age","genotype","mean","sd","se","n") #rename columns
  ddply(y, "age") #return dataframe sorted by age
})
names(list_df_sums) <- variables
df_summary <- do.call(rbind, list_df_sums)

list_df_tuk <- lapply(variables, function(x) {
  y <- as.data.frame(do.call(rbind, TukeyHSD(list_aov[[x]])))
  y <- merge(data.frame("permutation" = rownames(y), "lwr" = y$lwr), y)
  y <- merge(data.frame("variable" = x), y)
  ddply(y, c("variable","permutation"))
})
names(list_df_tuk) <- variables
df_tuks <- do.call(rbind, list_df_tuk)

##Percent Max Tension Conversions
df_maxTension <- data_proc
multiple_of_8 <- function(x) { if ((x %% length(x_axis)) != 0) x + (length(x_axis) - x %% length(x_axis)) else x}
df_maxTension$pred_Tension <- lapply(1:nrow(df_maxTension), function(x) df_maxTension$Tension[x] / single_fiber_parameters[[df_maxTension$fiber[x]]][['force_Max']] * 100)
df_maxTension$Tension <- as.numeric(lapply(1:nrow(df_maxTension), function(x) df_maxTension$Tension[x] / df_maxTension$Tension[multiple_of_8(x)] * 100))

##Non-linear Regression of Percent Max Tension
Pforce_Ca_function <- Tension ~ 100*(1/(1+10^(nH*(pCa-pCa50))))
p_max_analysis_nls <- lapply(t_comb, function(x) nls(Pforce_Ca_function, data = df_maxTension[isFactor(x),],  start = list(nH = 1, pCa50 = 6), trace = TRUE))
names(p_max_analysis_nls) <- names(t_comb)
p_max_single_fiber_nls <- lapply(1:number_fibers, function(x) nls(Pforce_Ca_function, data = df_maxTension[((8*x-7):(8*x)),],  start = list(nH = 3, pCa50 = 6), trace = TRUE))
names(p_max_single_fiber_nls) <- t(fiber_data1[,1])

#get bestFit curve of Max Tension
pMax_parameters <- lapply(1:length(p_max_analysis_nls), function(x) p_max_analysis_nls[[x]]$m$getPars())
names(pMax_parameters) <- names(t_comb)
pMax_single_fiber_parameters <- lapply(1:number_fibers, function(x) p_max_single_fiber_nls[[x]]$m$getPars())
names(pMax_single_fiber_parameters) <- names(p_max_single_fiber_nls)
pMax_bestFit <- function(x, i) {100*(1/(1+10^(pMax_parameters[[i]][1]*(x - pMax_parameters[[i]][2]))))}
pMax_avg_pCa50 <- mean(unlist(lapply(1:number_fibers, function(x) pMax_single_fiber_parameters[[x]][2])))

#creating plot with ggplot2
data_proc <- within(data_proc,  factor_c <- paste(factor, factor2, sep=" "))
data_proc$factor_c <- factor(data_proc$factor_c, levels = c("NTG 7","HCM 7","NTG 14","HCM 14","NTG 28","HCM 28"))
mypalette = brewer.pal(7,"Set1") #colors we want to use
df_plot <- subset(data_proc, data_proc$factor2 == 28)
df_plot$factor_c <- factor(df_plot$factor_c, levels = c("NTG 28","HCM 28"))
names(mypalette)[1:2] <- c("HCM 28", "NTG 28")
y_axis_lim <- c(0,30)
p_norm <- ggplot(df_plot, aes(x = pCa, y = Tension, color = as.factor(factor_c), fill = as.factor(factor_c))) + #setup ggplot
  scale_color_manual(values = c(mypalette[1], mypalette[2])) + labs(title = NULL, x = "pCa", y = "Tension (mN/mm2)") +
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 1.25, fun.args = list(mult=1), width = 0.125) +
  stat_summary(fun.y = 'mean', geom = "point", size = 2.5) +
  theme_classic() + scale_x_reverse(breaks = seq(min(df_plot$pCa), max(df_plot$pCa), by = 0.5)) + coord_cartesian(ylim = y_axis_lim) +
  lapply(levels(df_plot$factor_c), function (i) stat_function(fun = bestFit, args = i, color = mypalette[i], size = 1.25)) +
  lapply(levels(df_plot$factor_c), function (i) geom_segment(aes(x = parameters[[i]]$estimate[2], y = 0, xend = parameters[[i]]$estimate[2], yend = parameters[[i]]$estimate[3]/2), linetype = "dashed", color = mypalette[i], size = 1.25)) + 
  lapply(levels(df_plot$factor_c), function (i) geom_segment(aes(x = parameters[[i]]$estimate[2] - parameters[[i]]$std.error[2], y = parameters[[i]]$estimate[3]/2, xend = parameters[[i]]$estimate[2] + parameters[[i]]$std.error[2], yend = parameters[[i]]$estimate[3]/2), size = 1.25, linetype = "solid", color = "black", arrow = arrow(angle = 90, length = unit(0.05, "inches"),ends = "both", type = "open")))  +
  theme(axis.text = element_text(face = "bold", color = "black", size = 24),
        text = element_text(face = "bold", color = "black", size = 24),
        plot.title = element_text(size = 32, hjust = 0.5), legend.title = element_blank())
  

#create max tension plot with ggplot2
p_max <- ggplot(data = df_maxTension, mapping = aes(x = pCa, y = Tension, fill = factor, color = factor)) + #setup ggplot
  scale_color_manual(values = mypalette) +
  theme_classic() + scale_x_reverse() +
  stat_summary(fun.data = mean_se, geom = "errorbar", size = .5, fun.args = list(mult=1), width = 0.125) +
  stat_summary(fun.y = 'mean', geom = "point", size = 1) +
  lapply(l_factors, function (i) stat_function(fun = pMax_bestFit, args = i, color = mypalette[i])) +
  lapply(l_factors, function (i) geom_segment(aes(x = pMax_parameters[[i]][2], y = 0, xend = pMax_parameters[[i]][2], yend = 50), linetype = "dashed", color = mypalette[i])) + 
  lapply(l_factors, function (i) geom_segment(aes(x = pMax_parameters[[i]][2] - summary(p_max_analysis_nls[[i]])$coefficients[2,2], y = 50, xend = pMax_parameters[[i]][2] + summary(p_max_analysis_nls[[i]])$coefficients[2,2], yend = 50), linetype = "solid", color = "black", arrow = arrow(angle = 90, length = unit(0.05, "inches"),ends = "both", type = "open"))) +
  labs(title = "% Force-Calcium Relation", x = "pCa", y = "% Max Tension (mN/mm2)") +
  theme(axis.text = element_text(face = "bold", color = "black", size = 9), text = element_text(face = "bold", color = "black", size = 9), plot.title = element_text(size = 11, hjust = 0.5), legend.title = element_blank())

##Plots residual values
p_resid <- ggplot(aes(x=.fitted, y=.resid), data = augment(analysis_nls[[1]],data_proc[isFactor(1),])) +
  geom_point() + geom_hline(yintercept = 0) + theme_classic() + xlab("fitted values") + ylab("residuals") +
  geom_point(aes(x=.fitted, y=.resid), data = augment(analysis_nls[[2]],data_proc[isFactor(2),]))

##write results into file
export_summary <- function(file_name = "analysis", directory = NA)
{
  if (is.na(directory)) directory = dirname(filename)
  df_parameterSummary <- t(as.data.frame(lapply(1:number_fibers, function(x) unlist(as.vector(t(single_fiber_parameters[x]))))))
  df_pMax_parameterSummary <- t(as.data.frame(lapply(1:number_fibers, function(x) unlist(as.vector(t(pMax_single_fiber_parameters[x]))))))
  df_parameterSummary <- cbind(df_parameterSummary, df_pMax_parameterSummary)
  names(df_parameterSummary) <- c("nH","pCa50","force_Max","nH_pM","pCa50_pM")
  row.names(df_parameterSummary) <- t(fiber_data1[,1])
  write.csv(df_parameterSummary, file = normalizePath(file.path(directory, "Analysis", file_name)))
}

##function that corrects for an error in nlstools
fit_params <- function(nls_object)
{
  x <<- 0
  return(lapply(1:length(nls_object), function(i) {x <<- i
    nlsBoot(analysis_nls[[i]])}))
  remove(x)
}

##function for determining nonlinear regression of samples
nonlinear_regression <- function(df_fibers, method = "std")
{
  force_Ca_function <- Tension ~ force_Max/(1+10^(nH*(pCa-pCa50)))
  if (method == "std")
  {
    analysis_nls <- lapply(1:length(l_factors), function(x) nls(force_Ca_function, data = df_fibers[isFactor(x),],  start = list(nH = 3, pCa50 = 6, force_Max = 42), trace = TRUE))
  }
  if (method == "robust")
  {
    analysis_nls <- lapply(1:length(l_factors), function(x) nlrob(force_Ca_function, data = df_fibers[isFactor(x),],  start = list(nH = 3, pCa50 = 6, force_Max = 42), trace = TRUE))
  }
  names(analysis_nls) <- l_factors
  return(analysis_nls)
}

##This needs more work one day in the future
df_last <- data.frame("pCa50" = t(as.data.frame(lapply(1:number_fibers, function(x) unlist(as.vector(t(single_fiber_parameters[x]))))))[,2],
                      "genotype" = fiber_data1$genotype,
                      "age" = fiber_data1$age)
df_last$genotype <- factor(df_last$genotype, levels = c("NTG","HCM"))
ggplot(data = df_last, aes(x = as.factor(age), y = pCa50, fill = genotype)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1), width = .4, position = position_dodge(width = 0.4), color = "black", size = 1.2) +
  geom_dotplot(binaxis = 'y', position = position_dodge(width = 0.4), stackdir='center', dotsize = 1.25) + theme_classic() +
  scale_fill_brewer(palette = "Set1", direction = -1) + labs(x = "Age (days)", y = "pCa50") + guides(fill=guide_legend(title="genotype")) +
  stat_compare_means(aes(group = df_last[["Genotype"]]), label = "p.signif",
                     label.y = c(max(as.numeric(df_last[(df_last$age == 7),1]))*1.01,
                                 max(as.numeric(df_last[(df_last$age == 14),1]))*1.01,
                     max(as.numeric(df_last[(df_last$age == 28),1]))*1.01))



p_test <- ggplot(df_plot, aes(x = pCa, y = Tension, color = as.factor(factor_c), fill = as.factor(factor_c))) + #setup ggplot
  scale_color_manual(values = c(mypalette[1], mypalette[2])) + labs(title = NULL, x = "pCa", y = "Tension (mN/mm2)") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", size = 1.25, fun.args = list(mult=1), width = 0.125) +
  #stat_summary(fun.y = 'mean', geom = "point", size = 2.5) +
  theme_classic() + scale_x_reverse(breaks = seq(min(df_plot$pCa), max(df_plot$pCa), by = 0.5)) + coord_cartesian(ylim = y_axis_lim) +
  lapply(levels(df_plot$factor_c), function (i) stat_function(fun = bestFit, args = i, color = mypalette[i], size = 1.25)) +
  lapply(levels(df_plot$factor_c), function (i) geom_segment(aes(x = parameters[[i]]$estimate[2], y = 0, xend = parameters[[i]]$estimate[2], yend = parameters[[i]]$estimate[3]/2), linetype = "dashed", color = mypalette[i], size = 1.25)) + 
  lapply(levels(df_plot$factor_c), function (i) geom_segment(aes(x = parameters[[i]]$estimate[2] - parameters[[i]]$std.error[2], y = parameters[[i]]$estimate[3]/2, xend = parameters[[i]]$estimate[2] + parameters[[i]]$std.error[2], yend = parameters[[i]]$estimate[3]/2), size = 1.25, linetype = "solid", color = "black", arrow = arrow(angle = 90, length = unit(0.05, "inches"),ends = "both", type = "open")))  +
  theme(axis.text = element_text(face = "bold", color = "black", size = 24),
        text = element_text(face = "bold", color = "black", size = 24),
        plot.title = element_text(size = 32, hjust = 0.5), legend.title = element_blank())

names(mypalette) <- levels(data_proc$factor_c)
pmax_test <- ggplot(data = df_maxTension, mapping = aes(x = pCa, y = Tension, fill = factor, color = factor)) + #setup ggplot
  scale_color_manual(values = mypalette) +
  theme_classic() + scale_x_reverse() +
  #stat_summary(fun.data = mean_se, geom = "errorbar", size = .5, fun.args = list(mult=1), width = 0.125) +
  #stat_summary(fun.y = 'mean', geom = "point", size = 1) +
  lapply(levels(data_proc$factor_c)[c(3,4,6)], function (i) stat_function(fun = pMax_bestFit, args = i, color = mypalette[i])) +
  lapply(levels(data_proc$factor_c)[c(3,4,6)], function (i) geom_segment(aes(x = pMax_parameters[[i]][2], y = 0, xend = pMax_parameters[[i]][2], yend = 50), linetype = "dashed", color = mypalette[i])) + 
  lapply(levels(data_proc$factor_c)[c(3,4,6)], function (i) geom_segment(aes(x = pMax_parameters[[i]][2] - summary(p_max_analysis_nls[[i]])$coefficients[2,2], y = 50, xend = pMax_parameters[[i]][2] + summary(p_max_analysis_nls[[i]])$coefficients[2,2], yend = 50), linetype = "solid", color = "black", arrow = arrow(angle = 90, length = unit(0.05, "inches"),ends = "both", type = "open"))) +
  labs(title = "% Force-Calcium Relation", x = "pCa", y = "% Max Tension (mN/mm2)") +
  theme(axis.text = element_text(face = "bold", color = "black", size = 9), text = element_text(face = "bold", color = "black", size = 9), plot.title = element_text(size = 11, hjust = 0.5), legend.title = element_blank())
