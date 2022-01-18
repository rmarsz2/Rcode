
st.err <- function(x) sd(x)/sqrt(length(x))

iterate_bf.test <- function(df_data, var_data, layer = "age", layer2 = "genotype")
{
  library(onewaytests)
  #var_data <- colnames(df_data)[-which(colnames(df_data) %in% c(layer,layer2))]
  brown_forsythe = vector(mode = "list", length = length(levels(df_data[[layer]])))
  for(i in 1:length(brown_forsythe))
  {
    brown_forsythe[[i]] <- lapply(var_data, function(x) {
      f <- as.formula(paste(x, "~", layer2)) #define function with two independant variables (two-way)
      bf.test(f, df_data[df_data[[layer]] == levels(df_data[[layer]])[i],c("genotype", x)])})
    names(brown_forsythe[[i]]) <- variables
  }
  names(brown_forsythe) <- levels(df_data[[layer]])
  return(brown_forsythe)
}

iterate_levene.test <- function(df_data, var_data, layer = "age", layer2 = "genotype")
{
  library(car)
  #var_data <- colnames(df_data)[-which(colnames(df_data) %in% c(layer,layer2))]
  levene_out = vector(mode = "list", length = 2)
  center_method <- list(mean, median)
  for(i in 1:2)
  {
    levene_out[[i]] <- lapply(var_data, function(x) {
      f <- as.formula(paste(x, "~", layer, "*", layer2, sep = ""))
      leveneTest(f, data = df_data, center = center_method[[i]])})
    names(levene_out[[i]]) <- variables
  }
  names(levene_out) <- c("Levene","Brown.Forsythe")
  
  p_bf <- unlist(lapply(var_data, function(x) levene_out$Brown.Forsythe[[x]]$`Pr(>F)`[1]))
  names(p_bf) <- var_data
  output <- list(levene_out, p_bf)
  names(output) <- c("rawOutput", "pValues")
  return(output)
}

iterate_anova <- function(df_data, var_subset, layer = "age", layer2 = "genotype", subset = NULL)
{
  list_aov <- lapply(var_subset, function(x) {
    lm_aov <- lm(as.formula(paste(x, "~", layer, "*", layer2, sep = "")), data = df_data,
                 contrasts = list(genotype = contr.sum, age = contr.sum), subset[[x]])
    rstatix::Anova(lm_aov, type = 3)})
  names(list_aov) <- var_subset
  table_aovs <- lapply(var_subset, function(x) {
    as_tibble(data.frame(c(list_aov[[x]]$`Pr(>F)`[2:4], "2W ANOVA")))
  })
  table_aovs <- Reduce(function(x, y) cbind(x, y), table_aovs)
  #summarise_all(table_aovs)
  names(table_aovs) <- var_subset
  row.names(table_aovs) <- c(rownames(list_aov[[1]])[2:4], "Method")
  return(table_aovs)
}

iterate_gfd <- function(df_data, var_subset, layer = "age", layer2 = "genotype", subset = NULL)
{
  library(GFD)
  if(is.null(subset))
  {
    subset <- lapply(var_subset, function(x) rep(TRUE, nrow(df_data)))
    names(subset) <- var_subset
  }
  if(class(df_data) != "list")
  {
    df_data <- lapply(var_subset, function(x) {
      subset_df <- df_data[complete.cases(df_data[[x]]),c(layer,layer2,x)]
      subset_df[subset[[x]],]
    })
    names(df_data) <- var_subset
  }
  list_gfd <- lapply(var_subset, function(x) {
    f <- as.formula(paste(x, "~", layer, "*", layer2, sep = ""))
    GFD(f, data = df_data[[x]], alpha = 0.05, CI.method = "perm")
  })
  names(list_gfd) <- var_subset
  table_gfd <- lapply(var_subset, function(x) as_tibble(data.frame(c(list_gfd[[x]]$WTS[,3], "Wald"))))
  table_gfd <- Reduce(function(x, y) cbind(x, y), table_gfd)
  #summarise_all(table_aovs)
  names(table_gfd) <- var_subset
  row.names(table_gfd) <- c(rownames(list_gfd[[1]]$WTS),"method")
  return(table_gfd)
}

iterate_tukey <- function(df_data, var_subset, layer = "age", layer2 = "genotype", subset = NULL)
{
  list_df_tuk <- lapply(var_subset, function(x) {
    y <- as.data.frame(do.call(rbind, TukeyHSD(aov(as.formula(paste(x, "~", layer, "*", layer2, sep = "")), data = df_data, contrasts = list(genotype = contr.sum, age = contr.sum), subset = subset[[x]]))))
    y <- merge(data.frame("permutation" = rownames(y), "lwr" = y$lwr), y)
    y <- merge(data.frame("variable" = x), y)
    ddply(y, c("variable","permutation"))
  })
  names(list_df_tuk) <- var_subset
  df_tuks <- do.call(cbind, lapply(var_subset, function(x) list_df_tuk[[x]]$`p adj`))
  colnames(df_tuks) <- var_subset
  v_permutations <- list_df_tuk[[1]]$permutation
  split_perm <- strsplit(v_permutations, "[[:punct:]]")
  names(split_perm) <- v_permutations
  exclude_perms <- unlist(lapply(1:length(v_permutations), function(x) ("28d" %in% split_perm[[x]]) && (("7d" %in% split_perm[[x]]) || ("14d" %in% split_perm[[x]]))))
  rownames(df_tuks) <- v_permutations
  df_tuks_select <- df_tuks[v_permutations,]
  df_tuks_select <- df_tuks[!exclude_perms,]
  df_tuks_select <- rbind(df_tuks_select, "Tukey")
  rownames(df_tuks_select) <- c(v_permutations[!exclude_perms], "Method")
  return(df_tuks_select)
}

iterate_t.test <- function(data, var_subset, layer = "age", layer2 = "genotype")
{
  levels1 <- levels(data[[layer]])
  levels2 <- levels(data[[layer2]])
  l_ttest <- vector(mode = "list", length = length(levels(data[[layer]])))
  for(i in 1:3)
  {
    l_tamahane_a[[i]] <- lapply(var_subset, function(x) tamhaneT2Test(aov(as.formula(paste(x, "~", layer2)), data = data[complete.cases(data[[x]]),],
                                                                          contrasts = list(genotype = contr.sum), subset = (data$age == levels(data$age)[i]) & l_outlier_subsets[[x]]), p.adjust.method = "none"))
    names(l_tamahane_a[[i]]) <- var_subset
  }
  names(l_tamahane_a) <- levels1
  
  l_tamahane_g <- vector(mode = "list", length = length(levels(data[[layer2]])))
  for(i in 1:2)
  {
    l_tamahane_g[[i]] <- lapply(var_subset, function(x) tamhaneT2Test(aov(as.formula(paste(x, "~", layer)), data = data[complete.cases(data[[x]]),],
                                                                          contrasts = list(age = contr.sum), subset = (data$genotype == levels(data$genotype)[i]) & l_outlier_subsets[[x]]), p.adjust.method = "none"))
    names(l_tamahane_g[[i]]) <- var_subset
  }
  names(l_tamahane_g) <- levels2
  tamahane_pa <- lapply(var_subset, function(y) p.adjust(unlist(lapply(levels1, function(x) l_tamahane_a[[x]][[y]]$p.value)), method = "bonferroni"))
  names(tamahane_pa) <- var_subset
  df_tam_pa <- do.call(cbind, tamahane_pa)
  tamahane_pg <- lapply(var_subset, function(y) p.adjust(unlist(lapply(levels2, function(x) l_tamahane_g[[x]][[y]]$p.value)), method = "bonferroni"))
  names(tamahane_pg) <- var_subset
  df_tam_pg <- do.call(cbind, tamahane_pg)
  names_pa <- c(sapply(levels1, function(x) paste(x, ":", levels2[1],"-", x, ":", levels2[2], sep = "")))
  names_pg <- c(sapply(levels2, function(x) paste(rep(paste(colnames(l_tamahane_g[[x]][[1]]$p.value), ":", x, sep = ""), each = 2), "-",
                                                  rep(paste(rownames(l_tamahane_g[[x]][[1]]$p.value), ":", x, sep = ""), 2), sep = "")))
  split_perm <- strsplit(names_pg, "[[:punct:]]")
  exclude_perms <- unlist(lapply(1:length(split_perm), function(x) ("28d" %in% split_perm[[x]]) && (("7d" %in% split_perm[[x]]) || ("14d" %in% split_perm[[x]]))))
  df_tam_pg <- df_tam_pg[!exclude_perms,]
  df_tam <- rbind(df_tam_pa,df_tam_pg)
  df_tam <- rbind(df_tam, "Tamhane")
  rownames(df_tam) <- c(names_pa, names_pg[!exclude_perms], "Method")
  df_tam <- as.data.frame(df_tam[!is.na(df_tam)[,1],])
  return(df_tam)
}

iterate_tamhane <- function(data, var_subset, layer = "age", layer2 = "genotype", subset = NULL)
{
  library(PMCMRplus)
  levels1 <- levels(data[[layer]])
  levels2 <- levels(data[[layer2]])
  l_tamahane_a <- vector(mode = "list", length = length(levels(data[[layer]])))
  for(i in 1:3)
  {
    l_tamahane_a[[i]] <- lapply(var_subset, function(x) tamhaneT2Test(aov(as.formula(paste(x, "~", layer2)), data = data[complete.cases(data[[x]]),],
                                                                             contrasts = list(genotype = contr.sum), subset = (data$age == levels(data$age)[i]) & l_outlier_subsets[[x]]), p.adjust.method = "none"))
    names(l_tamahane_a[[i]]) <- var_subset
  }
  names(l_tamahane_a) <- levels1
  
  l_tamahane_g <- vector(mode = "list", length = length(levels(data[[layer2]])))
  for(i in 1:2)
  {
    l_tamahane_g[[i]] <- lapply(var_subset, function(x) tamhaneT2Test(aov(as.formula(paste(x, "~", layer)), data = data[complete.cases(data[[x]]),],
                                                                             contrasts = list(age = contr.sum), subset = (data$genotype == levels(data$genotype)[i]) & l_outlier_subsets[[x]]), p.adjust.method = "none"))
    names(l_tamahane_g[[i]]) <- var_subset
  }
  names(l_tamahane_g) <- levels2
  tamahane_pa <- lapply(var_subset, function(y) unlist(lapply(levels1, function(x) l_tamahane_a[[x]][[y]]$p.value)))
  names(tamahane_pa) <- var_subset
  tamahane_pg <- lapply(var_subset, function(y) unlist(lapply(levels2, function(x) l_tamahane_g[[x]][[y]]$p.value)))
  names(tamahane_pg) <- var_subset
  names_pa <- c(sapply(levels1, function(x) paste(x, ":", levels2[1],"-", x, ":", levels2[2], sep = "")))
  names_pg <- c(sapply(levels2, function(x) paste(rep(paste(colnames(l_tamahane_g[[x]][[1]]$p.value), ":", x, sep = ""), each = 2), "-",
                                           rep(paste(rownames(l_tamahane_g[[x]][[1]]$p.value), ":", x, sep = ""), 2), sep = "")))
  split_perm <- strsplit(c(names_pg), "[[:punct:]]")
  exclude_perms <- unlist(lapply(1:length(split_perm), function(x) ("28d" %in% split_perm[[x]]) && (("7d" %in% split_perm[[x]]) || ("14d" %in% split_perm[[x]]))))
  tamhane_combo <- lapply(var_subset, function(x) p.adjust(c(tamahane_pa[[x]],tamahane_pg[[x]][!exclude_perms]), method = "bonferroni"))
  names(tamhane_combo) <- var_subset
  df_tam <- do.call(cbind, tamhane_combo)
  df_tam <- rbind(df_tam, "Tamhane")
  rownames(df_tam) <- c(names_pa,names_pg[!exclude_perms], "Method")
  df_tam <- as.data.frame(df_tam[!is.na(df_tam)[,1],])
  return(df_tam)
}

iterate_dunn <- function(data, var_subset, layer = "age", layer2 = "genotype", perms)
{
  library(FSA)
  levels1 <- levels(data[[layer]])
  levels2 <- levels(data[[layer2]])
  l_dunn_a <- vector(mode = "list", length = 3)
  for(i in 1:3)
  {
    l_dunn_a[[i]] <- lapply(var_subset, function(x) {
      tbl_subset <-  data[which((data[[layer]] == levels(data[[layer]])[i]) & l_outlier_subsets[[x]]), c(layer2,x)]
      f <- as.formula(paste(tbl_subset[x], "~", layer2)) #formula with two variables
      dunnTest(f, data = tbl_subset, method = "none")
    })
    names(l_dunn_a[[i]]) <- var_subset
  }
  names(l_dunn_a) <- levels1
  
  l_dunn_g <- vector(mode = "list", length = 2)
  for(i in 1:2)
  {
    l_dunn_g[[i]] <- lapply(var_subset, function(x) {
      tbl_subset <-  data[which((data$genotype == levels(data$genotype)[i]) & l_outlier_subsets[[x]]), c(layer,x)]
      f <- as.formula(paste(tbl_subset[x],"~", layer)) #formula with two variables
      dunnTest(f, data = tbl_subset, method = "none")
    })
    names(l_dunn_g[[i]]) <- var_subset
  }
  names(l_dunn_g) <-  levels2
  dunn_a <- lapply(var_subset, function(y) c(unlist(lapply(1:3, function(x) l_dunn_a[[x]][[y]]$res$P.unadj)), unlist((lapply(1:2, function(x) l_dunn_g[[x]][[y]]$res$P.unadj)))))
  names(dunn_a) <- var_subset
  #sapply(c("7d","14d","28d"), function(x) paste(x, unlist(strsplit(rownames(test2), " - "))[1:2], sep = "", collapse = "-"))
  #paste(as.data.frame(sapply(c("NTG","HCM"), function(x) paste(unlist(strsplit(rownames(test2), " - "))[7:18], x, sep = "")))$NTG, as.data.frame(sapply(c("NTG","HCM"), function(x) paste(unlist(strsplit(rownames(test2), " - "))[7:18], x, sep = "")))$HCM, sep = "-")
  row_dunn <- c(sapply(levels1, function(x) paste(x, unlist(strsplit(l_dunn_a[[1]][[1]]$res$Comparison, " - ")), sep = "", collapse = "-")),
                c(unlist(lapply(levels2, function(i) (paste(t(sapply(1:3, function(x) paste(unlist(strsplit(l_dunn_g[[1]][[1]]$res$Comparison, " - "))[(2*x-1):(2*x)], ":", i, sep = "")))[,1],
                                                                   t(sapply(1:3, function(x) paste(unlist(strsplit(l_dunn_g[[1]][[1]]$res$Comparison, " - "))[(2*x-1):(2*x)], ":", i, sep = "")))[,2], sep = "-"))))))
  split_perm <- strsplit(row_dunn, "[[:punct:]]")
  exclude_perms <- unlist(lapply(1:length(split_perm), function(x) ("28d" %in% split_perm[[x]]) && (("7d" %in% split_perm[[x]]) || ("14d" %in% split_perm[[x]]))))
  
  #paste(as.data.frame(sapply(c("NTG","HCM"), function(x) paste(unlist(strsplit(l_dunn_g[[1]][[1]]$res$Comparison, " - ")), x, sep = "")))$NTG,
  #      as.data.frame(sapply(c("NTG","HCM"), function(x) paste(unlist(strsplit(l_dunn_g[[1]][[1]]$res$Comparison, " - ")), x, sep = "")))$HCM, sep = "-")
  #for(i in var_subset) names(dunn_a[[i]]) <- select_perms[c(9,4,6,2,1,3,7,5,8)]
  df_dunn <- do.call(cbind, lapply(var_subset, function(x) p.adjust(dunn_a[[x]][!exclude_perms], method = "bonferroni")))
  df_dunn <- rbind(df_dunn, "Dunn")
  colnames(df_dunn) <- var_subset
  rownames(df_dunn) <- c(row_dunn[!exclude_perms], "Method")
  return(df_dunn)
}


bracket_trans <- function(data, min = NA)
{
  library(stringr)
  n_p <- nrow(data)
  dodging_width = 0.25*.25
  tmp <- data.frame("genotype" = as.factor(c(do.call(rbind,str_split(data$permutation, "[[:punct:]]"))[,c(4,2)])),
                    "age" = as.factor(c(do.call(rbind,str_split(data$permutation, "[[:punct:]]"))[,c(3,1)])),
                    "y.position" = rep(data$y.position, 2), "pval" = rep(data$pval, 2))
  exp_group_label <- levels(tmp$genotype)[which(levels(tmp$genotype) != "NTG")]
  if(nrow(data) == 0)
    return(tmp)
  levels_age <- levels(tmp$age)
  tmp$genotype <- factor(tmp$genotype, levels = c("NTG", exp_group_label))
  tmp$x_position <- 1 - dodging_width
  tmp$x_position[which(tmp$age == levels(tmp$age)[grep("14", levels_age)])] <- 2 - dodging_width
  tmp$x_position[which(tmp$age == levels(tmp$age)[grep("28", levels_age)])] <- 3 - dodging_width
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

plot_data <- function(data, graphValue, p_val)
{
  stat_data <- data.frame("group1" = levels(df_data$age)[c(1,2,3,1,1,2,1,1,2)][which(p_val < 0.1)],
                          "group2" = levels(df_data$age)[c(1,2,3,2,3,3,2,3,3)][which(p_val < 0.1)],
                          'permutation' = rownames(combine_p)[1:(nrow(combine_p)-1)][which(p_val < 0.1)],
                          "pval" = signif(as.numeric(combine_p[[graphValue]][1:(nrow(combine_p)-1)]),2)[which(p_val < 0.1)],
                          "y.position" = seq(max(list_norm[[graphValue]]$mean + 2*max(list_norm[[graphValue]]$se)),
                                             max(list_norm[[graphValue]]$mean + 2*max(list_norm[[graphValue]]$se))*(1+length(which(p_val < 0.1))/25),
                                             length = length(which(p_val < 0.1))))
  p <- ggplot(data = current_df, mapping = aes(x = age, y = NOX4, fill = genotype, color = as.factor(genotype))) +
    stat_summary(fun.data = mean_se, geom = "errorbar", size = 1.75, fun.args = list(mult=1), width = .5, position = position_dodge(width = .8), color = "black") +
    stat_summary(aes(group = as.factor(genotype)), fun = mean, geom = "bar", position = position_dodge(width = 0.8), size = 1.25, width = .8, color = "black") +
    geom_dotplot(binaxis = 'y', position = position_dodge(0.8), color = "black", method = "histodot", stackdir = "center", dotsize = 1.5) +
    theme_classic() + scale_color_manual(values = c(myPalette[2], myPalette[1])) +
    labs(title = title, x = "Age (days)", y = y_axis_title) + scale_fill_manual(values = c(myPalette[2], myPalette[1])) +
    scale_y_continuous(sec.axis = sec_axis( trans=~./axis_corrections[[graphValue]])) +
    theme(axis.text = element_text(face = "bold", color = "black", size = 11),
          text = element_text(face = "bold", color = "black", size = 13),
          plot.title = element_text(size = 11, hjust = 0.5),
          legend.position = "none", legend.title = element_blank(),
          axis.line = element_line(colour = 'black', size = 1.25),
          axis.ticks = element_line(colour = "black", size = 1.25)) +
    geom_vline(xintercept = 2.5, linetype = "dashed")
  if(nrow(plot_stats) > 1)
  {
    p <- p +
      lapply(1:(nrow(plot_stats)/2), function(i) geom_segment(data = plot_stats, aes(x = x_position[i], xend = x_position[np+i],y = y.position[i], yend = y.position[np+i], group = genotype),color = "black")) +
      annotate(geom = "text", x=plot_stats$mid_x[1:(nrow(plot_stats)/2)],
               y=plot_stats$y.position[1:(nrow(plot_stats)/2)],
               label= stat_data$pval[1:(nrow(plot_stats)/2)],
               vjust = -.3, size = 3.5)
  }
  return(p)
}