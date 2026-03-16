library(tidyverse)
library(brms)
library(marginaleffects)
library(ggdist)

plot_marginal_effects <- function(
    model,
    data,
    component = c("mu", "hu"),
    effect,
    by,
    title,
    filename,
    similarity_values = NULL,
    type = "response",
    ftitle = 40/4,
    faxis.x = 36/4,
    faxis.y = 36/4,
    faxis.title.x = 36/4,
    linewidth = 20/4
) {
  
  component <- match.arg(component)
  
  newdata <- NULL
  by_vars <- by
  
  if (!is.null(similarity_values)) {
    
    if (by == "origin") {
      newdata <- datagrid(
        model = model,
        similarity_z = similarity_values,
        origin = unique(data$origin)
      )
    }
    
    if (by == "position") {
      newdata <- datagrid(
        model = model,
        similarity_z = similarity_values,
        position = unique(data$position)
      )
    }
    
    by_vars <- c(by, "similarity_z")
  }
  
  p <- avg_comparisons(
    model,
    variables = stats::setNames(list("pairwise"), effect),
    by = by_vars,
    type = type,
    dpar = component,
    newdata = newdata
  ) %>%
    get_draws() %>%
    mutate(
      contrast_label = paste0(contrast, " | ", by, ": ", .data[[by]])
    ) %>%
    ggplot(aes(x = draw, y = contrast_label)) +
    stat_interval(
      show.legend = FALSE,
      linewidth = linewidth,
      .width = c(.50, .89, .95),
      point_interval = "median_hdci"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5,
      color = "red"
    ) +
    scale_color_brewer(palette = "Greys") +  # Optional color scale
    ylab("") +
    xlab(paste0("Marginal effect on ", component, " (", type, " scale)")) +
    ggtitle(title) +
    theme_bw() +
    theme(
      plot.title = element_text(size = ftitle, hjust = 0.5),
      axis.text.x = element_text(size = faxis.x),
      axis.text.y = element_text(size = faxis.y),
      axis.title.x = element_text(size = faxis.title.x)
    )
  
  if (!is.null(similarity_values)) {
    p <- p + facet_wrap(~ similarity_z)
  }
  
  #ggsave(filename, plot = p, width = 16, height = 14, dpi = "retina")
  
  return(p)
}

maltese_hurdle <- read_rds("results/maltese_hurdle.rds")

plot_marginal_effects(
  model = maltese_hurdle,
  data = maltese_hurdle$data,
  component = "hu",
  effect = "origin",
  by = "position",
  similarity_values = c(-1, 0, 1),
  title = "Position contrasts within origin at selected similarity levels",
  filename = "hu_position_origin_similarity.png"
)

mu_pos_by_orig <- 
  plot_marginal_effects(
    model = maltese_hurdle,
    component = "mu",
    effect = "position",
    by = "origin",
    title = "Position contrasts within Origin (mu component)",
    filename = "mu_position_by_origin.png"
    )

mu_orig_by_pos <- 
  plot_marginal_effects(
    model = maltese_hurdle,
    component = "mu",
    effect = "origin",
    by = "position",
    title = "Origin contrasts within Position (mu component)",
    filename = "mu_origin_by_position.png"
    )

hu_pos_by_orig <- 
  plot_marginal_effects(
    model = maltese_hurdle,
    component = "hu",
    effect = "position",
    by = "origin",
    title = "Position contrasts within Origin (hu component)",
    filename = "hu_position_by_origin.png"
    )

hu_orig_by_pos <- 
  plot_marginal_effects(
    model = maltese_hurdle,
    component = "hu",
    effect = "origin",
    by = "position",
    title = "Origin contrasts within Position (hu component)",
    filename = "hu_origin_by_position.png"
    )