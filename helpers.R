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


plot_stats <- function(data) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    #count the number of cases per day and group and simulation
    count <- data %>%
      group_by(simulation, group, date_onset) %>%
      count()
    peak_df <-  get_peak(data, TRUE)
    # Incidence
    incidence <-
      ggplot2::ggplot(data = count,
                      ggplot2::aes(x = date_onset, y = n, fill = group)) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::geom_vline(
        data = peak_df,
        aes(xintercept = peak),
        linewidth = 2 * 1.5,
        show.legend = FALSE
      ) +
      ggplot2::geom_vline(
        data = peak_df,
        aes(xintercept = peak, col = group),
        linewidth = 2,
        show.legend = FALSE
      ) +
      ggplot2::labs(x = "Date Onset", y = "Cases", fill = "Group") +
      ggplot2::ggtitle("Incidence (all simulations stacked)") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")
    
    return(incidence)
  } else {
    stop("Please install ggplot2 to use this function")
  }
}


get_mixing <-
  function(data,
           min_t = min(data$date_infection),
           max_t = max(data$date_infection)) {
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

p_Ri <- function(data, param_info) {
  Ri <-
    map(.x = unique(data$simulation), ~ get_Ri(data[data$simulation == .x,])) %>%
    bind_rows(., .id = "simulation")
  
  Ri_long <- Ri %>%
    pivot_longer(cols = param_info$name,
                 names_to = "target") %>%
    select(simulation, id, group, target, value, date_infection)
  
  Ri_sims <- Ri_long %>%
    group_by(group, target, date_infection, simulation) %>%
    summarise(Ri = mean(value))
  
  total_Ri_sims <- Ri %>%
    group_by(simulation, group, date_infection) %>%
    summarise(Ri = mean(Ri))
  
  
  p_facet_Ri <- ggplot(data = Ri_sims,
                       aes(
                         x = as.integer(date_infection),
                         y = Ri,
                         col = group
                       )) +
    facet_grid(target ~ group, scales = "free_y") +
    geom_point(alpha = 0.5) +
    geom_smooth(col = "black", se = FALSE) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    ggtitle("Case Reproduction Number Matrix") +
    labs(caption = "x = source group, y = target group")
  
  
  
  p_total_Ri <- ggplot(total_Ri_sims,
                       aes(
                         x = as.integer(date_infection),
                         y = Ri,
                         col = group,
                         group = group
                       )) +
    geom_point() +
    geom_smooth(col = "black", se = FALSE) +
    geom_hline(
      data = tibble(group = param_info$name,
                    r0 = param_info$r0),
      aes(yintercept = r0, group = group),
      col = "#4C4E52",
      lty = "solid"
    ) +
    geom_hline(aes(yintercept = 1), col = "#4C4E52", lty = "dotted") +
    facet_grid(~ group) +
    theme_bw() +
    theme(strip.text.x = element_blank(),
          legend.position = "none") +
    labs(caption = "Case Reproduction Number by Source Group", x = "Date")
  
  p_Ri <-
    patchwork::wrap_plots(
      p_facet_Ri,
      p_total_Ri,
      ncol = 1,
      guides = "collect",
      heights = c(3, 1)
    )
  return(p_Ri)
}


p_mixing <- function(data, param_info) {
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
                         ~ factor(., levels = param_info$name)))
  
  
  mean_freqs <- mix_t %>%
    group_by(t, source_group, group) %>%
    summarise(freq = mean(freq))
  
  # inital mixing frequencies
  Truth <- generate_Mcol(
    n_groups = param_info$n_groups,
    size = param_info$size,
    name = param_info$name,
    delta = param_info$delta
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
    geom_smooth(col = "black", linewidth = 1) +
    geom_hline(
      data = Truth,
      aes(yintercept = value),
      col = "#4C4E52",
      lty = "solid"
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, 1)) +
    ggtitle("Transmission Frequency Matrix") +
    labs(x = "date_infection", y = "Frequency", caption = "x = source group, y = target group")
  
  return(p_mixing)
}


p_delta <- function(data, param_info) {
  #estimate delta
  truth <-
    data.frame(truth = param_info$delta, name = param_info$name)
  
  
  peaks <- get_peak(data, TRUE)
  est <- map2(peaks$peak, peaks$group, function(peak, group) {
    o2groups::early_delta(
      data = data,
      min_t = 0,
      max_t = peak,
      size = param_info$size,
      name = param_info$name
    ) %>% as.data.frame() %>% rownames_to_column(var = "name") %>%
      mutate(group_peak = group)
  }) %>% bind_rows()
  
  results <- left_join(truth, est, by = "name")
  
  #plot with error bars
  ggplot(results) +
    geom_point(aes(x = name, y = est, col = name)) +
    geom_errorbar(aes(
      x = name,
      y = est,
      col = name,
      ymin = lower_ci ,
      ymax = upper_ci
    ),
    width = 0.2) +
    #geom_errorbar(aes(x = name, y=truth,ymin=truth,ymax=truth), color = "black", width = 0.5)+
    geom_hline(aes(yintercept = truth)) +
    facet_grid(name ~ group_peak, scales = "free") +
    labs(x = "Group Peak", y = "Estimated value",
         caption = "columns = group's relative peak date, row = group estimate") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle(paste0("Estimated Assortativity Coefficient with 95% CI"))
}
