library(ggplot2)
library(RColorBrewer)
library(readxl)
library(dplyr)
library(Hmisc)
library(Rmisc)
library(ggpubr)
library(gridExtra)
library(tidyr)
library(onewaytests)
library(lme4)
library(nlme)
library(emmeans)
library(cowplot)
library(lmtest) #Breusch-Pagan Test
library(car) #Levene test for variance
library(multcomp)
library(rstatix)
library(GFD)

#FILE_SEP <<- .Platform$file.sep
#HOME_DIR <<- normalizePath(path = 
#open file
filename = normalizePath(file.path(Sys.getenv('HOMEPATH'), "Box",
                                   "Richard",
                                   "Data",
                                   "Echocardiography",
                                   "7 14 28 do echoes.xlsx"))
#open data from excel
echo_data <- lapply(excel_sheets(filename), function(x) read_excel(filename, sheet = x)) #extract data
#make sure names from echo data are correct
names(echo_data) <- excel_sheets(filename)
#create a genotype column as factor from the data file
genotype <- as.factor(c(rep("TnT", nrow(echo_data$`Day28 TnT`)),
                        rep("NTG", nrow(echo_data$`Day28 NTG`)),
                        rep("NTG", nrow(echo_data$`Day14 NTG`)),
                        rep("NTG", nrow(echo_data$`Day7 NTG`)),
                        rep("TnT", nrow(echo_data$`Day7 TnT`)),
                        rep("TnT", nrow(echo_data$`Day14 TnT`))))
#create an age column as factor from the data file
age <- as.factor(c(rep.int(28, nrow(echo_data$`Day28 TnT`)),
                   rep.int(28, nrow(echo_data$`Day28 NTG`)),
                   rep.int(14, nrow(echo_data$`Day14 NTG`)),
                   rep.int(7, nrow(echo_data$`Day7 NTG`)),
                   rep.int(7, nrow(echo_data$`Day7 TnT`)),
                   rep.int(14, nrow(echo_data$`Day14 TnT`))))
#merge all data in one large dataframe
df_echoes <- do.call(rbind,echo_data)
#insert genotype column into this dataframe and set orders of factors to make pretty plots
df_echoes$genotype <- factor(genotype, levels = c("NTG", "TnT"))
#insert age column into dataframe
df_echoes$age <- age
#calculate the AT/FT ratio for coronary flow
df_echoes$'AT/FT' <- df_echoes$`Dias Cor AT`/df_echoes$`Dias Cor FT`
df_echoes$`LXs/LXd` <- df_echoes$`Coronary LXs`/df_echoes$`Coronary LXd`
df_echoes$`SXs/SXd` <- df_echoes$`Coronary SXs`/df_echoes$`Coronary SXd`
#creates more program friendly names for the column names
names(df_echoes) <- make.names(names(df_echoes), unique = TRUE)

#remove accessory variables from memory
remove(genotype)
remove(age)

#perform two-way anova on all variables
variables <- names(df_echoes)[-which(names(df_echoes) == c("age", "genotype"))] #returns variables interested in being analyzed
variables <- variables[5:length(variables)]
#apply ANOVA over all variables

st.err <- function(x)
  sd(x)/sqrt(length(x))

df_echoes$genotype <- factor(df_echoes$genotype, c("NTG", "TnT"))
df_echoes$age <- factor(df_echoes$age, c(7,14,28))

#Shapiro-Wilks test for normality
l_shap_wilk <- lapply(variables, function(x) shapiro.test(resid(lm(df_echoes[[x]] ~ genotype*age, data = df_echoes))))
names(l_shap_wilk) <- variables

deviations <- unlist(lapply(variables, function(x) l_shap_wilk[[x]]$p.value)) < 0.05
dev_subset <- variables[deviations]

########################
#__________________________HOMOGENEITY OF VARIANCE________________________________
########################
#BROWN-FORSYTHE with 1-Way ANOVA
brown_forsythe = vector(mode = "list", length = length(levels(df_echoes$age)))
for(i in 1:length(brown_forsythe))
{
  brown_forsythe[[i]] <- lapply(variables, function(x) {
    f <- as.formula(paste(x, "~genotype")) #define function with two independant variables (two-way)
    bf.test(f, df_echoes[df_echoes$age == levels(df_echoes$age)[i],c("genotype", x)])})
  names(brown_forsythe[[i]]) <- variables
}
names(brown_forsythe) <- levels(df_echoes$age)
#LEVENE TEST & Brown-Forsythe with Multiple Groups
levene_out = vector(mode = "list", length = 2)
center_method <- list(mean, median)
for(i in 1:2)
{
  levene_out[[i]] <- lapply(variables, function(x) {
    f <- as.formula(paste(x, "~genotype*age"))
    leveneTest(f, data = df_echoes, center = center_method[[i]])})
  names(levene_out[[i]]) <- variables
}
names(levene_out) <- c("Levene","Brown.Forsythe")
p_bf <- unlist(lapply(variables, function(x) levene_out$Brown.Forsythe[[x]]$`Pr(>F)`[1]))
names(p_bf) <- variables
hetero_subset <- variables[p_bf < 0.05]
#########################
#_______________________OUTLIER DETECTION_________________________________
#########################
#hatvalues(test_assumptions(test_var = 'LA', df_echoes))
l_cooks_distances <- lapply(variables, function(x) {
  cooks.distance(lm(paste(x, "~ genotype*age"), data = df_echoes))
})
names(l_cooks_distances) <- variables
l_df_out <- lapply(variables, function(x) {
  subset_df <- df_echoes[complete.cases(df_echoes[[x]]),]
  subset_df[(l_cooks_distances[[x]] < 4 / nrow(subset_df)),c("age","genotype",x)]
})
names(l_df_out) <- variables
l_outlier_subsets <- lapply(variables, function(x) l_cooks_distances[[x]] < 4 / nrow(df_echoes[complete.cases(df_echoes[[x]]),]))
names(l_outlier_subsets) <- variables
#########################
#________________________SUMMARIZE DATA________________________________
#########################
list_df_sums <- lapply(variables, function(x) {
  f <- as.formula(paste(df_echoes[x],"~genotype*age")) #formula with two variables
  y <- aggregate(f, df_echoes, mean, subset = l_outlier_subsets[[x]]) #compute summary statistics
  y <- merge(y, aggregate(f, df_echoes, sd, subset = l_outlier_subsets[[x]]), by = c("age", "genotype"))
  y <- merge(y, aggregate(f, df_echoes, st.err, subset = l_outlier_subsets[[x]]), by = c("age", "genotype"))
  y <- merge(y, aggregate(f, df_echoes, length, subset = l_outlier_subsets[[x]]), by = c("age", "genotype"))
  names(y) <- c("age","genotype","mean","sd","se","n") #rename columns
  ddply(y, "age") #return dataframe sorted by age
})
names(list_df_sums) <- variables

##create summary data frame with all variables included
df_summary <- do.call(rbind, list_df_sums)
########################
#__________________________TWO-WAY ANALYSIS AND TABLE GENERATION________________________________
########################
#Declaring variables to be analyzed a certain way
var_2W_Anova <- variables[-which(variables %in% c(dev_subset, hetero_subset))]
var_gfd <- variables[which(variables %in% c(dev_subset, hetero_subset))]

list_aov <- lapply(var_2W_Anova, function(x) {
  lm_aov <- lm(as.formula(paste(x, "~genotype*age")), data = df_echoes,
               contrasts = list(genotype = contr.sum, age = contr.sum), subset = l_outlier_subsets[[x]])
  rstatix::Anova(lm_aov, type = 3)})
names(list_aov) <- var_2W_Anova
table_aovs <- lapply(var_2W_Anova, function(x) {
  as_tibble(data.frame(c(list_aov[[x]]$`Pr(>F)`[2:4], "2W ANOVA")))
})
table_aovs <- Reduce(function(x, y) cbind(x, y), table_aovs)
#summarise_all(table_aovs)
names(table_aovs) <- var_2W_Anova
row.names(table_aovs) <- c(rownames(list_aov$HEART.RATE)[2:4], "Method")
list_gfd <- lapply(var_gfd, function(x) {
  f <- as.formula(paste(x, "~genotype*age"))
  GFD(f, data = l_df_out[[x]], alpha = 0.05, CI.method = "perm")
})
names(list_gfd) <- var_gfd
table_gfd <- lapply(var_gfd, function(x) as_tibble(data.frame(c(list_gfd[[x]]$WTS[,3], "Wald"))))
table_gfd <- Reduce(function(x, y) cbind(x, y), table_gfd)
#summarise_all(table_aovs)
names(table_gfd) <- var_gfd
row.names(table_gfd) <- c(rownames(list_gfd[[1]]$WTS),"method")
sum_table <- t(cbind(table_aovs, table_gfd))
sum_table <- sum_table[variables,]
remove(table_aovs, table_gfd)

combined_2W_list <- c(list_aov,list_gfd)
names(combined_2W_list) <- c(names(list_aov),names(list_gfd))



#############################
#____________________________REDO WITHOUT OUTLIERS
#############################
#l_shap_wilk_out <- lapply(variables, function(x) shapiro.test(resid(lm(l_df_out[[x]][[x]] ~ genotype*age, data = l_df_out[[x]]))))
#names(l_shap_wilk_out) <- variables

#deviations_out <- unlist(lapply(variables, function(x) l_shap_wilk_out[[x]]$p.value)) < 0.05
#dev_subset_out <- variables[deviations_out]

#levene_out_out = vector(mode = "list", length = 2)
#center_method <- list(mean, median)
#for(i in 1:2)
#{
#  levene_out_out[[i]] <- lapply(variables, function(x) {
#    f <- as.formula(paste(x, "~genotype*age"))
#    leveneTest(f, data = df_echoes, center = center_method[[i]])})
#  names(levene_out_out[[i]]) <- variables
#}
#names(levene_out_out) <- c("Levene","Brown.Forsythe")
#p_bf_out <- unlist(lapply(variables, function(x) levene_out_out$Brown.Forsythe[[x]]$`Pr(>F)`[1]))
#names(p_bf_out) <- variables
#hetero_subset_out <- variables[p_bf_out < 0.05]


#var_2W_Anova_out <- variables[-which(variables %in% c(dev_subset_out, hetero_subset_out))]
#var_gfd_out <- variables[which(variables %in% c(dev_subset_out, hetero_subset_out))]

#list_aov_out <- lapply(var_2W_Anova_out, function(x) {
#  lm_aov <- lm(as.formula(paste(x, "~genotype*age")), data = l_df_out[[x]],
#               contrasts = list(genotype = contr.sum, age = contr.sum))
#  rstatix::Anova(lm_aov, type = 3)})
#names(list_aov_out) <- var_2W_Anova_out
#table_aovs_out <- lapply(var_2W_Anova_out, function(x) {
#  as_tibble(data.frame(c(list_aov_out[[x]]$`Pr(>F)`[2:4], "2W ANOVA")))
#})
#table_aovs_out <- Reduce(function(x, y) cbind(x, y), table_aovs_out)
#names(table_aovs_out) <- var_2W_Anova_out
#row.names(table_aovs_out) <- c(rownames(list_aov_out$HEART.RATE)[2:4], "Method")
#list_gfd_out <- lapply(var_gfd_out, function(x) {
#  f <- as.formula(paste(x, "~genotype*age"))
#  GFD(f, data = l_df_out[[x]][,c("age","genotype",x)], alpha = 0.05, CI.method = "perm")
#})
#names(list_gfd_out) <- var_gfd_out
#table_gfd_out <- lapply(var_gfd_out, function(x) as_tibble(data.frame(c(list_gfd_out[[x]]$WTS[,3], "Wald"))))
#table_gfd_out <- Reduce(function(x, y) cbind(x, y), table_gfd_out)
#names(table_gfd_out) <- var_gfd_out
#row.names(table_gfd_out) <- c(rownames(list_gfd_out[[1]]$WTS),"method")
#sum_table_out <- t(cbind(table_aovs_out, table_gfd_out))
#remove(table_aovs_out, table_gfd_out)
#
#combined_2W_list_out <- c(list_aov_out,list_gfd_out)
#names(combined_2W_list_out) <- c(names(list_aov_out),names(list_gfd_out))



#list_df_sums_out <- lapply(variables, function(x) {
#  f <- as.formula(paste(l_df_out[[x]][x],"~genotype*age")) #formula with two variables
#  y <- aggregate(f, l_df_out[[x]], mean) #compute summary statistics
#  y <- merge(y, aggregate(f, l_df_out[[x]], sd), by = c("age", "genotype"))
#  y <- merge(y, aggregate(f, l_df_out[[x]], st.err), by = c("age", "genotype"))
#  y <- merge(y, aggregate(f, l_df_out[[x]], length), by = c("age", "genotype"))
#  names(y) <- c("age","genotype","mean","sd","se","n") #rename columns
#  ddply(y, "age") #return dataframe sorted by age
#})
#names(list_df_sums_out) <- variables

##create summary data frame with all variables included
#df_summary_out <- do.call(rbind, list_df_sums_out)

##########################
#_________________________POST-HOC COMPARISONS
##########################

#list_df_tuk_out <- lapply(variables, function(x) {
#  y <- as.data.frame(do.call(rbind, TukeyHSD(aov(as.formula(paste(x, "~genotype*age")), data = l_df_out[[x]],contrasts = list(genotype = contr.sum, age = contr.sum)))))
#  y <- merge(data.frame("permutation" = rownames(y), "lwr" = y$lwr), y)
#  y <- merge(data.frame("variable" = x), y)
#  ddply(y, c("variable","permutation"))
#})
#names(list_df_tuk_out) <- variables
#df_tuks_out <- do.call(rbind, list_df_tuk_out)
#df_tuks_out_select <- df_tuks_out[c("variable","permutation","p adj")]
#rownames(df_tuks_out_select) <- NULL
#write.table(df_tuks_out_select, "clipboard-256", sep="\t", row.names=FALSE)
#write.table(df_tuks_out_select, file = "df_tuks.txt", sep = ",", quote = FALSE, row.names = F)

list_df_tuk <- lapply(var_2W_Anova, function(x) {
  y <- as.data.frame(do.call(rbind, TukeyHSD(aov(as.formula(paste(x, "~genotype*age")), data = df_echoes,contrasts = list(genotype = contr.sum, age = contr.sum), subset = l_outlier_subsets[[x]]))))
  y <- merge(data.frame("permutation" = rownames(y), "lwr" = y$lwr), y)
  y <- merge(data.frame("variable" = x), y)
  ddply(y, c("variable","permutation"))
})
names(list_df_tuk) <- var_2W_Anova
df_tuks <- do.call(cbind, lapply(var_2W_Anova, function(x) list_df_tuk[[x]]$`p adj`))
colnames(df_tuks) <- var_2W_Anova
v_permutations <- list_df_tuk$HEART.RATE$permutation
rownames(df_tuks) <- v_permutations
select_perms <- v_permutations[c(4,6,7,11,13,15,17,18,19)]
df_tuks_select <- df_tuks[select_perms,]
df_tuks_select <- rbind(df_tuks_select, "Tukey")
rownames(df_tuks_select) <- c(select_perms, "Method")

write.table(df_tuks_select, "clipboard-256", sep="\t", row.names=FALSE)

library(PMCMRplus)

l_tamahane_a <- vector(mode = "list", length = 3)
for(i in 1:3)
{
  l_tamahane_a[[i]] <- lapply(hetero_subset, function(x) tamhaneT2Test(aov(as.formula(paste(x, "~genotype")), data = df_echoes[complete.cases(df_echoes[[x]]),],
                                                                           contrasts = list(genotype = contr.sum), subset = (df_echoes$age == levels(df_echoes$age)[i]) & l_outlier_subsets[[x]]), p.adjust.method = "none"))
  names(l_tamahane_a[[i]]) <- hetero_subset
}
names(l_tamahane_a) <- c("seven","fourteen","twnetyeight")

l_tamahane_g <- vector(mode = "list", length = 2)
for(i in 1:2)
{
  l_tamahane_g[[i]] <- lapply(hetero_subset, function(x) tamhaneT2Test(aov(as.formula(paste(x, "~age")), data = df_echoes[complete.cases(df_echoes[[x]]),],
                                                                           contrasts = list(age = contr.sum), subset = (df_echoes$genotype == levels(df_echoes$genotype)[i]) & l_outlier_subsets[[x]]), p.adjust.method = "none"))
  names(l_tamahane_g[[i]]) <- hetero_subset
}
names(l_tamahane_g) <- c("NTG","TnT")

tamahane_p <- c(unlist(lapply(1:3, function(x) l_tamahane_a[[x]]$ET$p.value)), unlist((lapply(1:2, function(x) l_tamahane_g[[x]]$ET$p.value))))
tamahane_p <- tamahane_p[!is.na(tamahane_p)]
names(tamahane_p) <- select_perms[c(9,4,6,1,3,2,5,8,7)]
p_adj_means_tam <- data.frame("Tamhane" = c(p.adjust(tamahane_p, method = "bonferroni")))
p_adj_means_tam <- c("Tamhane", p_adj_means_tam[select_perms,])

library(FSA)
dunnTest(as.formula(LA~genotype), data = df_echoes[,c("LA","genotype")], method = "none")
l_dunn_a <- vector(mode = "list", length = 3)
for(i in 1:3)
{
  
  l_dunn_a[[i]] <- lapply(dev_subset, function(x) {
    tbl_subset <-  df_echoes[which((df_echoes$age == levels(df_echoes$age)[i]) & l_outlier_subsets[[x]]), c("genotype",x)]
    f <- as.formula(paste(tbl_subset[x],"~genotype")) #formula with two variables
    dunnTest(f, data = tbl_subset, method = "none")
  })
  names(l_dunn_a[[i]]) <- dev_subset
}
names(l_dunn_a) <- c("seven","fourteen","twnetyeight")

l_dunn_g <- vector(mode = "list", length = 2)
for(i in 1:2)
{
  l_dunn_g[[i]] <- lapply(dev_subset, function(x) {
    tbl_subset <-  df_echoes[which((df_echoes$genotype == levels(df_echoes$genotype)[i]) & l_outlier_subsets[[x]]), c("age",x)]
    f <- as.formula(paste(tbl_subset[x],"~age")) #formula with two variables
    dunnTest(f, data = tbl_subset, method = "none")
  })
  names(l_dunn_g[[i]]) <- dev_subset
}
names(l_dunn_g) <-  c("NTG","TnT")


dunn_a <- lapply(dev_subset, function(y) c(unlist(lapply(1:3, function(x) l_dunn_a[[x]][[y]]$res$P.unadj)), unlist((lapply(1:2, function(x) l_dunn_g[[x]][[y]]$res$P.unadj)))))
names(dunn_a) <- dev_subset
for(i in dev_subset) names(dunn_a[[i]]) <- select_perms[c(9,4,6,2,1,3,7,5,8)]
p_adj_means_dunn <-sapply(dev_subset, function(x) p.adjust(dunn_a[[x]], method = "bonferroni"))
p_adj_means_dunn <- rbind(p_adj_means_dunn, "Dunn")
rownames(p_adj_means_dunn)[nrow(p_adj_means_dunn)] <- "Method"

combine_p <- merge(df_tuks_select, p_adj_means_dunn ,by="row.names", all.x=TRUE)
combine_p <- cbind(combine_p, p_adj_means_tam)
names(combine_p)[length(combine_p)] <- hetero_subset
rownames(combine_p) <- combine_p$Row.names
combine_p <- combine_p[-which(names(combine_p) == "Row.names")]
combine_p <- combine_p[,variables]

write.table(combine_p, "clipboard-256", sep="\t", row.names=TRUE)

qqPlot(lm(paste(variables[41],"~genotype*age"), df_echoes), ylab = "Sample Quantiles", xlab = "Theoretical Quantiles")

reformatted_summary <- vector(mode = "list", length = 6)
for(i in 1:3)
{
  reformatted_summary[[2*i-1]] <- paste(signif(df_summary[which(df_summary$age == c(7,14,28)[i] & df_summary$genotype == "NTG"),]$mean,3), " ± ", signif(df_summary[which(df_summary$age == c(7,14,28)[i] & df_summary$genotype == "NTG"),]$se,2))
  reformatted_summary[[2*i]] <- paste(signif(df_summary[which(df_summary$age == c(7,14,28)[i] & df_summary$genotype == "TnT"),]$mean,3), " ± ", signif(df_summary[which(df_summary$age == c(7,14,28)[i] & df_summary$genotype == "TnT"),]$se,2))
}
reformatted_summary <- do.call(cbind, reformatted_summary)
colnames(reformatted_summary) <- c("NTG 7", "TnT 7", "NTG 14", "TnT 14", "NTG 28", "TnT 28")
rownames(reformatted_summary)<- variables

results_2w <- paste("pg=", signif(as.numeric(sum_table[,"genotype"]),2),",pa=", signif(as.numeric(sum_table[,"age"]),2),",pga=", signif(as.numeric(sum_table[,"genotype:age"]),2))
reformatted_summary <- cbind(reformatted_summary, results_2w)
reformatted_summary <- cbind(reformatted_summary, sum_table[,"Method"])

write.table(reformatted_summary, "clipboard-256", sep="\t", row.names=TRUE)


plot(density(df_echoes$SV))
library(mvoutlier)

#Plot code
myPalette <- brewer.pal(3, "Set1")
graphValue <- "HEART.RATE"
title <- "Heart Rate"
y_axis_title <- "Heart Rate (bpm)"
current_df <- df_echoes[l_outlier_subsets[[graphValue]],]
current_df$age <- as.factor(current_df$age)
p_val <- signif(as.numeric(combine_p[[graphValue]][2:nrow(combine_p)]),2)
stat_data <- data.frame("group1" = levels(df_echoes$age)[c(2,3,3,2,2,3,3,3,1)][which(p_val < 0.1)],
                        "group2" = levels(df_echoes$age)[c(1,2,1,2,1,3,2,1,1)][which(p_val < 0.1)],
                        'permutation' = rownames(combine_p)[2:nrow(combine_p)][which(p_val < 0.1)],
                        "pval" = signif(as.numeric(combine_p[[graphValue]][2:nrow(combine_p)]),2)[which(p_val < 0.1)],
                        "y.position" = seq(max(list_df_sums[[graphValue]]$mean + 2*max(list_df_sums[[graphValue]]$se)),
                                           max(list_df_sums[[graphValue]]$mean + 2*max(list_df_sums[[graphValue]]$se))*(1+length(which(p_val < 0.1))/25),
                                           length = length(which(p_val < 0.1))))
plot_stats <- bracket_trans(stat_data)


HR <- ggplot(data = current_df, mapping = aes(x = age, y = HEART.RATE, fill = genotype, color = as.factor(genotype))) +
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 1.75, fun.args = list(mult=1), width = .25, position = "dodge") +
  stat_summary(aes(group = as.factor(genotype)), fun = mean, geom = "line", position = position_dodge(width = 0.25), size = 1.25) +
  theme_classic() + scale_color_manual(values = c(myPalette[2], myPalette[1])) +
  labs(title = title, x = "Age (days)", y = y_axis_title) + scale_fill_manual(values = c(myPalette[2], myPalette[1])) +
  theme(axis.text = element_text(face = "bold", color = "black", size = 11),
        text = element_text(face = "bold", color = "black", size = 13),
        plot.title = element_text(size = 11, hjust = 0.5),
        legend.position = "none", legend.title = element_blank(),
        axis.line = element_line(colour = 'black', size = 1.25),
        axis.ticks = element_line(colour = "black", size = 1.25))
  #lapply(1:(nrow(plot_stats)/2), function(i) geom_segment(data = plot_stats, aes(x = x_position[i], xend = x_position[np+i], y = y.position[i], yend = y.position[np+i], group = genotype),color = "black")) +
  #annotate(geom = "text", x=plot_stats$mid_x[1:(nrow(plot_stats)/2)],
  #         y=plot_stats$y.position[1:(nrow(plot_stats)/2)],
  #         label= stat_data$pval[1:(nrow(plot_stats)/2)],
  #         vjust = -.3, size = 3.5)

bracket_trans <- function(data, min = NA)
{
  library(stringr)
  n_p <- nrow(data)
  indexing_g <- 0
  indexing_a <- 0
  dodging_width = 0.25*.25
  if (is.even(grep("NTG", str_split(data$permutation, "[[:punct:]]")[[1]])))
  {
    indexing_g <- 2
    indexing_a <- 1
  } else {
    indexing_g <- 1
    indexing_a <- 2
  }
  tmp <- data.frame("genotype" = as.factor(c(do.call(rbind,str_split(data$permutation, "[[:punct:]]"))[,c(indexing_g+2,indexing_g)])),
                    "age" = as.factor(c(do.call(rbind,str_split(data$permutation, "[[:punct:]]"))[,c(indexing_a+2,indexing_a)])),
                    "y.position" = rep(data$y.position, 2), "pval" = rep(data$pval, 2))
  exp_group_label <- levels(tmp$genotype)[which(levels(tmp$genotype) != "NTG")]
  if(nrow(data) == 0)
    return(tmp)
  tmp$genotype <- factor(tmp$genotype, levels = c("NTG", exp_group_label))
  tmp$x_position <- 1 - dodging_width
  tmp$x_position[which(tmp$age == "14")] <- 2 - dodging_width
  tmp$x_position[which(tmp$age == "28")] <- 3 - dodging_width
  tmp$x_position[which(tmp$genotype == exp_group_label)] <- tmp$x_position[which(tmp$genotype == exp_group_label)] + 2*dodging_width
  tmp$segment_length <- abs(tmp$x_position[(n_p+1):(2*n_p)] - tmp$x_position[1:n_p])
  tmp$mid_x <- (tmp$x_position[1:n_p]+tmp$x_position[(n_p+1):(2*n_p)])/2
  if(n_p > 1)
  {
    sorting_index <- sort(tmp$segment_length[1:n_p], decreasing = T, index.return = TRUE)$ix
    y_sort <- sort(tmp$y.position[1:n_p], decreasing = T)
    adjust_yn <- n_p+1 - length(which(tmp$segment_length[1:n_p] < dodging_width*4))
    if (adjust_yn <= n_p) {
      new_y <- c(rev(y_sort[n_p:(n_p-adjust_yn+1)]), rep(y_sort[n_p],n_p-adjust_yn))
    } else {
      new_y <- y_sort
    }
    tmp[sorting_index,]$y.position <- new_y
    tmp[sorting_index+n_p,]$y.position <- new_y
  }
  if(!is.na(min))
    tmp$y.position <- tmp$y.position - (min(tmp$y.position) - min)
  tmp$np <- n_p
  return(tmp)
}

ggarrange(SV, EF, cO, HR, nrow = 1, legend = "bottom", common.legend = TRUE)
#v_adj_sh-(max(stat_data$y.position)-min(list_df_sums[[graphValue]]$mean) - max(list_df_sums[[graphValue]]$se))/10