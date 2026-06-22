library(PhyloSim)
library(ggplot2)
library(dplyr)
library(ape)
library(picante)
library(iCAMP)
library(cowplot)
library(tidyr)

# ================================================
# STEP 1: 模型模拟
# ================================================

results_neutral_global    <- list()
results_neutral_local     <- list()
results_density_global    <- list()
results_density_local     <- list()
results_env_global        <- list()
results_env_local         <- list()
results_combined_global   <- list()
results_combined_local    <- list()

for (i in 1:10) {
  cat("Running seed", i, "of 10...\n")
  
  parmBatch <- list(
    createCompletePar(x=100, y=100, dispersal=0, runs=10000,
                      density=0, environment=0, specRate=5, seed=i, type="base"),
    createCompletePar(x=100, y=100, dispersal=1, runs=10000,
                      density=0, environment=0, specRate=5, seed=i, type="base"),
    createCompletePar(x=100, y=100, dispersal=0, runs=10000,
                      density=1, environment=0, specRate=5, seed=i, type="base"),
    createCompletePar(x=100, y=100, dispersal=1, runs=10000,
                      density=1, environment=0, specRate=5, seed=i, type="base"),
    createCompletePar(x=100, y=100, dispersal=0, runs=10000,
                      density=0, environment=1, specRate=5, seed=i, type="base"),
    createCompletePar(x=100, y=100, dispersal=1, runs=10000,
                      density=0, environment=1, specRate=5, seed=i, type="base"),
    createCompletePar(x=100, y=100, dispersal=0, runs=10000,
                      density=1, environment=1, specRate=5, seed=i, type="base"),
    createCompletePar(x=100, y=100, dispersal=1, runs=10000,
                      density=1, environment=1, specRate=5, seed=i, type="base")
  )
  
  simuBatch <- runSimulationBatch(parmBatch, parallel=8)
  
  results_neutral_global[[i]]  <- simuBatch[[1]]
  results_neutral_local[[i]]   <- simuBatch[[2]]
  results_density_global[[i]]  <- simuBatch[[3]]
  results_density_local[[i]]   <- simuBatch[[4]]
  results_env_global[[i]]      <- simuBatch[[5]]
  results_env_local[[i]]       <- simuBatch[[6]]
  results_combined_global[[i]] <- simuBatch[[7]]
  results_combined_local[[i]]  <- simuBatch[[8]]
}

# ================================================
# STEP 2: 数据提取函数
# ================================================

# 提取群落矩阵（每行=1个样方）
make_comm <- function(simu) {
  spec_mat    <- simu$Output[[1]]$specMat
  n_rows      <- nrow(spec_mat)
  all_species <- paste0("s", unique(as.vector(spec_mat)))
  
  comm <- matrix(0, nrow=n_rows, ncol=length(all_species),
                 dimnames=list(paste0("site", 1:n_rows), all_species))
  
  for (i in 1:n_rows) {
    sp  <- paste0("s", spec_mat[i,])
    tab <- table(sp)
    comm[i, names(tab)] <- as.numeric(tab)
  }
  return(comm)
}

# 提取并修剪系统发育树
make_tree <- function(simu, comm) {
  tree       <- simu$Output[[1]]$phylogeny
  current_sp <- colnames(comm)
  tree_pruned <- keep.tip(tree, current_sp)
  return(tree_pruned)
}

# ================================================
# STEP 3: 计算NRI 和 NTI
# ================================================

calc_nri_nti <- function(simu_list, model_label, dispersal_label) {
  results <- data.frame()
  
  for (i in 1:length(simu_list)) {
    
    # 提取数据
    spec_vec        <- as.vector(simu_list[[i]]$Output[[1]]$specMat)
    species_names   <- paste0("s", spec_vec)
    abundance_table <- table(species_names)
    
    # 整个网格作为1个群落
    comm <- matrix(
      as.numeric(abundance_table),
      nrow     = 1,
      dimnames = list("1", names(abundance_table))
    )
    
    tree        <- simu_list[[i]]$Output[[1]]$phylogeny
    current_sp  <- names(abundance_table)
    tree_pruned <- keep.tip(tree, current_sp)
    dist_matrix <- cophenetic(tree_pruned)
    
    # 对齐物种顺序
    comm <- comm[, tree_pruned$tip.label, drop=FALSE]
    
    # 加假行绕过picante bug
    comm_double <- rbind(comm, dummy=rep(0, ncol(comm)))
    
    # 计算NRI
    nri_result <- ses.mpd(comm_double, dist_matrix,
                          null.model         = "taxa.labels",
                          abundance.weighted = TRUE,
                          runs               = 99)
    
    # 计算NTI
    nti_result <- ses.mntd(comm_double, dist_matrix,
                           null.model         = "taxa.labels",
                           abundance.weighted = TRUE,
                           runs               = 99)
    
    results <- rbind(results, data.frame(
      seed      = i,
      model     = model_label,
      dispersal = dispersal_label,
      NRI       = -nri_result[1, "mpd.obs.z"],
      NTI       = -nti_result[1, "mntd.obs.z"]
    ))
  }
  return(results)
}

# Calculating NRI and NTI
cat("Calculating NRI and NTI...\n")
nri_nti_all <- rbind(
  calc_nri_nti(results_neutral_global,  "Neutral",       "Global"),
  calc_nri_nti(results_neutral_local,   "Neutral",       "Local"),
  calc_nri_nti(results_density_global,  "Density Dep",   "Global"),
  calc_nri_nti(results_density_local,   "Density Dep",   "Local"),
  calc_nri_nti(results_env_global,      "Env Filtering", "Global"),
  calc_nri_nti(results_env_local,       "Env Filtering", "Local"),
  calc_nri_nti(results_combined_global, "Combination",   "Global"),
  calc_nri_nti(results_combined_local,  "Combination",   "Local")
)

# ================================================
# STEP 4: 计算 betaNTI
# ================================================

calc_bNTI <- function(simu_list, model_label, dispersal_label) {
  results <- data.frame()
  
  for (i in 1:length(simu_list)) {
    
    comm        <- make_comm(simu_list[[i]])
    tree_pruned <- make_tree(simu_list[[i]], comm)
    dis         <- cophenetic(tree_pruned)
    comm        <- comm[, rownames(dis)]
    
    bNTI_result <- bNTIn.p(
      comm      = comm,
      dis       = dis,
      rand      = 99,
      weighted  = TRUE,
      nworker   = 2,
      sig.index = "bNTI"
    )
    
    bNTI_values <- bNTI_result$index[lower.tri(bNTI_result$index)]
    
    results <- rbind(results, data.frame(
      seed           = i,
      model          = model_label,
      dispersal      = dispersal_label,
      bNTI_mean      = mean(bNTI_values),
      bNTI_sd        = sd(bNTI_values),
      pct_det_pos    = mean(bNTI_values >  2),
      pct_det_neg    = mean(bNTI_values < -2),
      pct_stochastic = mean(abs(bNTI_values) < 2)
    ))
  }
  return(results)
}

# 对所有模型计算betaNTI
cat("Calculating betaNTI...\n")
bNTI_all <- rbind(
  calc_bNTI(results_neutral_global,  "Neutral",       "Global"),
  calc_bNTI(results_neutral_local,   "Neutral",       "Local"),
  calc_bNTI(results_density_global,  "Density Dep",   "Global"),
  calc_bNTI(results_density_local,   "Density Dep",   "Local"),
  calc_bNTI(results_env_global,      "Env Filtering", "Global"),
  calc_bNTI(results_env_local,       "Env Filtering", "Local"),
  calc_bNTI(results_combined_global, "Combination",   "Global"),
  calc_bNTI(results_combined_local,  "Combination",   "Local")
)

# ================================================
# STEP 5: 汇总
# ================================================

# NRI/NTI汇总
summary_nri_nti <- nri_nti_all %>%
  group_by(model, dispersal) %>%
  summarise(
    NRI_mean = round(mean(NRI), 3),
    NRI_sd   = round(sd(NRI),   3),
    NTI_mean = round(mean(NTI), 3),
    NTI_sd   = round(sd(NTI),   3),
    .groups  = "drop"
  )

# betaNTI汇总
summary_bNTI <- bNTI_all %>%
  group_by(model, dispersal) %>%
  summarise(
    bNTI_mean      = round(mean(bNTI_mean),          3),
    bNTI_sd        = round(mean(bNTI_sd),             3),
    pct_det_pos    = round(mean(pct_det_pos)    * 100, 1),
    pct_det_neg    = round(mean(pct_det_neg)    * 100, 1),
    pct_stochastic = round(mean(pct_stochastic) * 100, 1),
    .groups        = "drop"
  )

print(summary_nri_nti)
print(summary_bNTI)

# ================================================
# STEP 6: 作图
# ================================================

my_colors <- c(
  "Neutral"       = "pink",
  "Density Dep"   = "darkblue",
  "Env Filtering" = "darkgreen",
  "Combination"   = "orange"
)

# 设置顺序
for (df in list(summary_nri_nti, summary_bNTI)) {
  df$model <- factor(df$model,
                     levels = c("Neutral", "Density Dep",
                                "Env Filtering", "Combination"))
  df$dispersal <- factor(df$dispersal, levels = c("Global", "Local"))
}

summary_nri_nti$model <- factor(summary_nri_nti$model,
                                levels = c("Neutral", "Density Dep",
                                           "Env Filtering", "Combination"))
summary_nri_nti$dispersal <- factor(summary_nri_nti$dispersal,
                                    levels = c("Global", "Local"))
summary_bNTI$model <- factor(summary_bNTI$model,
                             levels = c("Neutral", "Density Dep",
                                        "Env Filtering", "Combination"))
summary_bNTI$dispersal <- factor(summary_bNTI$dispersal,
                                 levels = c("Global", "Local"))

# NRI图
p_nri <- ggplot(summary_nri_nti,
                aes(x=model, y=NRI_mean, color=model, shape=dispersal)) +
  geom_point(size=4, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=NRI_mean - NRI_sd,
                    ymax=NRI_mean + NRI_sd),
                width=0.2, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray") +
  scale_color_manual(values=my_colors) +
  scale_shape_manual(values=c("Global"=16, "Local"=17)) +
  labs(x="Model", y="NRI", color="Model", shape="Dispersal") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# NTI图
p_nti <- ggplot(summary_nri_nti,
                aes(x=model, y=NTI_mean, color=model, shape=dispersal)) +
  geom_point(size=4, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=NTI_mean - NTI_sd,
                    ymax=NTI_mean + NTI_sd),
                width=0.2, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray") +
  scale_color_manual(values=my_colors) +
  scale_shape_manual(values=c("Global"=16, "Local"=17)) +
  labs(x="Model", y="NTI", color="Model", shape="Dispersal") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# betaNTI mean图
p_bNTI <- ggplot(summary_bNTI,
                 aes(x=model, y=bNTI_mean, color=model, shape=dispersal)) +
  geom_point(size=4, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=bNTI_mean - bNTI_sd,
                    ymax=bNTI_mean + bNTI_sd),
                width=0.2, position=position_dodge(width=0.5)) +
  geom_hline(yintercept= 2, linetype="dashed", color="red") +
  geom_hline(yintercept=-2, linetype="dashed", color="blue") +
  geom_hline(yintercept= 0, linetype="dashed", color="gray") +
  scale_color_manual(values=my_colors) +
  scale_shape_manual(values=c("Global"=16, "Local"=17)) +
  labs(x="Model", y="Mean betaNTI", color="Model", shape="Dispersal") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# betaNTI 比例堆积图
pct_data <- summary_bNTI %>%
  select(model, dispersal, pct_det_pos, pct_det_neg, pct_stochastic) %>%
  pivot_longer(
    cols      = c(pct_det_pos, pct_det_neg, pct_stochastic),
    names_to  = "process",
    values_to = "percentage"
  ) %>%
  mutate(process = recode(process,
                          pct_det_pos    = "Deterministic (+)",
                          pct_det_neg    = "Deterministic (-)",
                          pct_stochastic = "Stochastic"
  ))

p_pct <- ggplot(pct_data,
                aes(x=model, y=percentage, fill=process)) +
  geom_bar(stat="identity", position="stack") +
  facet_wrap(~dispersal) +
  scale_fill_manual(values=c(
    "Deterministic (+)" = "lightblue",
    "Deterministic (-)" = "yellowgreen",
    "Stochastic"        = "pink"
  )) +
  labs(x="Model", y="Percentage (%)", fill="Process") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# 拼图
legend <- get_legend(p_nri + theme(legend.direction="vertical"))

p_nri_nl  <- p_nri  + theme(legend.position="none")
p_nti_nl  <- p_nti  + theme(legend.position="none")
p_bNTI_nl <- p_bNTI + theme(legend.position="none")
p_pct_nl  <- p_pct  + theme(legend.position="none")
p_pct_legend <- get_legend(p_pct)

plot_grid(
  plot_grid(p_nri_nl, p_nti_nl, p_bNTI_nl, p_pct_nl,
            ncol=2, labels=c("A","B","C","D")),
  plot_grid(legend, p_pct_legend, ncol=1),
  ncol       = 2,
  rel_widths = c(1, 0.2)
)

