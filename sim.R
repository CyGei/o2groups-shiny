# Packages
library(tidyverse)
library(o2groups)
library(furrr)

# Helper functions
get_peak <- function(data, by_group = FALSE) {
  if (by_group) {
    
    result <- data %>%
      group_by(simulation, group, date_onset) %>%
      count() %>%
      group_by(simulation, group) %>%
      filter(n == max(n)) %>% 
      ungroup() %>% 
      group_by(group) %>% 
      summarise(
        peak = mean(date_onset),
        lower_ci = quantile(date_onset, 0.025),
        upper_ci = quantile(date_onset, 0.975)
      )

  } else {
    result <- data %>%
      group_by(simulation, date_onset) %>%
      count() %>%
      group_by(simulation) %>%
      filter(n == max(n)) %>% 
      ungroup() %>% 
      summarise(
        peak = mean(date_onset),
        lower_ci = quantile(date_onset, 0.025),
        upper_ci = quantile(date_onset, 0.975)
      )
  }
  
  return(result)
}


plot_stats <- function() {
  if (requireNamespace("ggplot2", quietly = TRUE)) {

    #count the number of cases per day and group and simulation
    count <-data %>%
      group_by(simulation, group, date_onset) %>%
      count()

    # Incidence
    incidence <-
      ggplot2::ggplot(data = count,
                      ggplot2::aes(x = date_onset, y = n, fill = group)) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::geom_vline(data = get_peak(data, TRUE), 
                          aes(xintercept = peak), linewidth = 2*1.5, show.legend = FALSE)+
      ggplot2::geom_vline(data = get_peak(data, TRUE), 
                          aes(xintercept = peak, col = group), linewidth = 2, show.legend = FALSE)+
      ggplot2::labs(x = "Date Onset", y = "Cases", fill = "Group") +
      ggplot2::ggtitle("Incidence (all simulations stacked)") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")

    
    # Proportion of Susceptible
    prop_susceptibles <-
      ggplot2::ggplot(data = stats,
                      ggplot2::aes(x = time, y = prop_susceptible, 
                                   col = group, 
                                   group = interaction(group, simulation))) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::scale_y_continuous(
        breaks = seq(0, 1, 0.1),
        limits = c(0, NA),
        expand = c(0, 0)
      ) +
      ggplot2::labs(x = "Date Infection", y = "Proportion of susceptible", col = "Group") +
      ggplot2::ggtitle("Proportion of susceptible (All simulations)") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    
    return(list(incidence,
                prop_susceptibles))
  } else {
    stop("Please install ggplot2 to use this function")
  }
}


get_mixing <-
  function(data,
           min_t = min(data$date_infection),
           max_t = max(data$date_infection)
  ) {
    
    df <-
      subset(
        data,
        date_infection >= min_t &
          date_infection <= max_t &
          !is.na(source),
        select = c(source_group, group)
      )
    
    counts <- t(table(df$source_group, df$group))
    freqs <- prop.table(counts, margin = 2)
    
    result <-
      expand.grid(source_group = rownames(freqs), group = colnames(freqs))
    result$n <-  as.vector(t(counts))
    result$freq = as.vector(t(freqs))
    
    return(list(Mcol = freqs,
                result = result))
  }

p_Ri <- function(){
  Ri <-
    map(.x = unique(data$simulation), ~ get_Ri(data[data$simulation == .x,])) %>%
    bind_rows(., .id = "simulation")
  
  Ri_long <- Ri %>%
    pivot_longer(cols = name,
                 names_to = "target") %>%
    select(simulation, id, group, target, value, date_infection)
  
  Ri_sims <- Ri_long %>%
    group_by(group, target, date_infection, simulation) %>%
    summarise(Ri = mean(value))
  
  total_Ri_sims <- Ri %>%
    group_by(simulation, group, date_infection) %>%
    summarise(Ri = mean(Ri))
  
  
  p_facet_Ri <- ggplot(data = Ri_sims,
                       aes(x = date_infection,
                           y = Ri,
                           col = group)) +
    facet_grid(target ~ group, scales = "free_y") +
    geom_point(alpha = 0.5) +
    geom_smooth(col = "black", se = FALSE, method = 'lm', formula = y ~ poly(x, 4))+
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none") +
    gg_col +
    gg_x_scale_infection+
    ggtitle("Case Reproduction Number Matrix")+
    labs(caption = "x = source group, y = target group")
  
  
  
  p_total_Ri <- ggplot(total_Ri_sims,
                       aes(x = date_infection,
                           y = Ri,
                           col = group,
                           group = group)) +
    geom_point() +
    geom_smooth(col = "black", se = FALSE, method = 'lm', formula = y ~ poly(x, 4))+ 
    geom_hline(data = tibble(group = name,
                             r0 = r0 ),
               aes(yintercept = r0, group = group),
               col = "#4C4E52", lty = "solid")+
    geom_hline(aes(yintercept = 1), col = "#4C4E52", lty = "dotted")+
    facet_grid(~ group) +
    theme_bw()+
    theme( strip.text.x = element_blank(),
           legend.position = "none")+
    gg_col +
    gg_x_scale_infection+
    labs(caption = "Case Reproduction Number by Source Group")
  
  p_Ri <- patchwork::wrap_plots(p_facet_Ri,p_total_Ri, ncol = 1, guides = "collect", 
                                heights = c(3,1))
  return(p_Ri)
}


p_mixing <- function(){
  # compute the mixing proportions at every timepoint
  mix_t <- purrr::map(.x = unique(data$simulation),
                      function(sim) {
                        sim_data <-  dplyr::filter(data, simulation == sim)
                        purrr::map(.x = unique(sim_data$date_infection),
                                   function(date) {
                                     mix <- sim_data %>% get_mixing(., min_t = date, max_t = date)
                                     return(mix$result)
                                   }) %>%
                          dplyr::bind_rows(.id = "t") %>%
                          dplyr::mutate(t = as.integer(t))
                      }) %>%
    dplyr::bind_rows(., .id = "simulation") %>%
    dplyr::mutate(across(c(source_group, group),
                         ~ factor(., levels = name)))
  
  
  mean_freqs <- mix_t %>%
    group_by(t, source_group, group) %>%
    summarise(freq = mean(freq))
  
  # inital mixing frequencies
  Truth <- generate_Mcol(
    n_groups = n_groups,
    size = size,
    name = name,
    delta = delta
  ) %>%
    as.data.frame(.) %>%
    rownames_to_column(var = "group") %>%
    pivot_longer(-group, names_to = "source_group", values_to = "value") %>%
    select(source_group, group, value) %>%
    arrange(source_group)
  
  
  p_mixing <- ggplot(mix_t) +
    aes(x = t,
        y = freq ,
        col = source_group) +
    facet_grid(group ~ source_group) +
    geom_point() +
    geom_smooth(col = "black", linewidth = 1)+
    geom_hline(data = Truth,
               aes(yintercept = value),  col = "#4C4E52", lty = "solid") +
    theme_bw() +
    theme(legend.position = "none") +
    gg_col +
    gg_x_scale_infection+
    scale_y_continuous(limits = c(0, 1))+
    ggtitle("Case Reproduction Number Matrix")+
    labs(x = "date_infection", y = "Frequency", caption = "x = source group, y = target group")
  
  return(p_mixing)
}


p_delta <- function(){
#estimate delta
truth <- data.frame(truth = delta, name = name)


peaks <- get_peak(data, by_group = TRUE)
est <- map2(peaks$peak, peaks$group, function(peak, group) {
  o2groups::early_delta(
    data = data,
    min_t = 0,
    max_t = peak,
    size = size,
    name = name
  ) %>% as.data.frame() %>% rownames_to_column(var = "name") %>% mutate(group_peak = group)
}) %>% bind_rows()

results <- left_join(truth, est, by = "name")

#plot with error bars
ggplot(results) +
  geom_point(aes(x = name, y = est, col = name)) +
  geom_errorbar(aes(x = name, y = est, col = name,
                    ymin = lower_ci , ymax = upper_ci), width = 0.2) +
  #geom_errorbar(aes(x = name, y=truth,ymin=truth,ymax=truth), color = "black", width = 0.5)+
  geom_hline(aes(yintercept =truth))+
  facet_grid(name ~ group_peak, scales = "free")+
  labs(x = "Group Peak", y = "Estimated value",
       caption = "columns = group's relative peak date, row = group estimate") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  gg_col +
  ggtitle(paste0("Estimated Assortativity Coefficient with 95% CI"))
}






# Parameters
n_simulations = 100
n_cores = parallel::detectCores() - 2
duration = 100
n_groups = 4
size = c(1e3, 1e4, 1e3, 1e3)
name = LETTERS[1:n_groups]
delta = c(10, 2, 1, 4)
intro_group = LETTERS[1:n_groups]
intro_n = rep(10, n_groups)
r0 = c(3, 4, 2, 3.5)
generation_time = simulacr::make_disc_gamma(mean = 5, sd = 2)$d(1:30)
incubation_period = simulacr::make_disc_gamma(mean = 5, sd = 2)$r(1000)

# Simulation
set.seed(123)

out <-
  simulate_groups_furrr(
    n_simulations = n_simulations,
    duration = duration,
    n_groups = n_groups,
    size = size,
    name = name,
    delta = delta,
    intro_n = intro_n,
    r0 = r0,
    generation_time =generation_time,
    incubation_period = incubation_period)

data = out$data
stats = out$stats
param_info <- list(
  n_groups = n_groups,
  name = name,
  size = size,
  delta = delta,
  r0 = r0,
  intro_n = intro_n)
# Plots
pal <- scales::hue_pal()(n_groups)
gg_col <- list(scale_color_manual(values = pal),
               scale_fill_manual(values = pal)) 

max_time_infection <- stats %>%
  filter(new_cases > 0) %>%
  summarise(max_time = max(time) + 1) %>%
  pull(max_time)

gg_x_scale_infection <- list(scale_x_continuous(limits = c(0, max_time_infection)))


p_tree <- plot_tree(subset(data, simulation == 1), pal = pal)
p_tree
p_hist <- plot_stats()[[1]] + gg_col + gg_x_scale_infection
p_hist
p_susceptibles <- plot_stats()[[2]] + gg_col + gg_x_scale_infection
p_susceptibles
p_Ri <- p_Ri()
p_Ri
p_mixing <- p_mixing()
p_mixing
p_delta <- p_delta()
p_delta
