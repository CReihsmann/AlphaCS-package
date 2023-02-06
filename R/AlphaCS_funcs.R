library(tidyverse)
library(ggplot2)
library(plotly)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(plyr)
library(pals)

#'---Merges data with pathway info
data_pthwy_merge <- function(data_path, comp_lib_path) {
  data <- read_csv(data_path)
  compounds <- read_excel(comp_lib_path, sheet = 2)
  compounds <- compounds %>%
    select(`Item Name`, Pathway, Target, Information, URL)
  data %>%
    left_join(compounds, by =c('Compound' = 'Item Name'))
}

#'---Converts appropriate columns to numeric
data_to_numeric <- function(df) {

  df %>%
    mutate_at(vars(matches('pg')), as.numeric)
}

#'---Make's necessary controls dfs
ctrls_data <- function(df_numeric) {

  controls_list <- c('G 1.7', 'G 1.7 + Arg 20', 'G 16.7')

  ctrls_tbl <- df_numeric %>%
    filter(Compound %in% controls_list) %>%
    filter(!is.na(`pg/mL`)) %>%
    group_by(Compound) %>%
    mutate_at(vars(contains('pg')), .funs = list(sd = ~sd(.))) %>%
    mutate_at(vars(`pg/mL`:`pg/min`), .funs = list(mean = ~mean(.))) %>%
    mutate_at(vars(`pg/mL`:`pg/min`), .funs = list(median = ~median(.))) %>%
    mutate_at(vars(`pg/mL`:`pg/min`), .funs = list(log = ~log(.))) %>%
    mutate_at(vars(`pg/mL`:`pg/min`), .funs = list(z_score = ~((. - mean(.)/sd(.))))) %>%
    mutate_at(vars(contains('log')), .funs = list(z_score = ~((. - mean(.)/sd(.))))) %>%
    mutate_at(vars(`pg/mL`:`pg/min`), .funs = list(std_err = ~sd(.)/sqrt(n())))

  ctrls_avg_plate <- ctrls_tbl %>%
    group_by(donor, plate_numb, .add = T) %>%
    summarise_at(vars(`pg/mL`:`pg/min`), mean) %>%
    mutate(Pathway = 'Controls') %>%
    arrange(`pg/mL`)

  ctrls_avg_stderr <- ctrls_tbl %>%
    select(Compound, `pg/mL_std_err`:`pg/min_std_err`) %>%
    distinct()

  ctrls_avg <- ctrls_tbl %>%
    group_by(donor, .add = T) %>%
    summarise_at(vars(`pg/mL`:`pg/min`), mean) %>%
    left_join(ctrls_avg_stderr)

  list(ctrls_tbl, ctrls_avg_plate, ctrls_avg_stderr, ctrls_avg)

}

#'---Pulls max value from column
max_pull <- function(df, conc) {
  ungroup(df) %>%
    pull(!!sym(conc)) %>%
    max(na.rm = T)
}

#'---Finds max value for given concentration
set_limit <- function(df, conc) {

  if(isTRUE(df %>% ungroup() %>% pull(`pg/mL`) %>% max(na.rm = T) < 700)) {
    if(conc == 'pg/mL'){
      700
    }
    else if(conc == 'pg/0.2mL'){
      140
    }
    else if(conc == 'pg/hr'){
      70
    }
    else{
      1.167
    }
  }
  else {
    if(conc == 'pg/mL') {
      df %>%
        max_pull(conc) %>%
        round_any(100, f = ceiling)
    }
    else if(conc == 'pg/0.2mL') {
      df %>%
        max_pull(conc) %>%
        round_any(20, f = ceiling)
    }
    else if(conc == 'pg/hr') {
      df %>%
        max_pull(conc) %>%
        round_any(10, f = ceiling)
    }
    else {
      df %>%
        max_pull(conc) %>%
        ceiling()
    }
  }

}
#'---produced final tibble needed for graphs
final_df <- function(df_numeric, ctrls_avg_plate, ctrls_avg) {

  controls_list <- c('G 1.7', 'G 1.7 + Arg 20', 'G 16.7')
  arg_glc_ctrl <- ctrl_mean(ctrls_avg, 'pg/mL', 'G 1.7 + Arg 20')

  limit_pgmL <- set_limit(df_numeric, 'pg/mL')
  limit_pg02mL <- set_limit(df_numeric, 'pg/0.2mL')
  limit_pghr <- set_limit(df_numeric, 'pg/hr')
  limit_pgmin <- set_limit(df_numeric, 'pg/min')

  df_numeric_prep <- df_numeric %>%
    filter(!Compound %in% controls_list) %>%
    arrange(Pathway, `pg/mL`)

  final_df <- ctrls_avg_plate %>%
    bind_rows(df_numeric_prep) %>%
    mutate(read_status = if_else(is.na(read_status), 'average', read_status)) %>%
    mutate(`pg/mL` = if_else(read_status == 'na_above', limit_pgmL, `pg/mL`),
           `pg/0.2mL` = if_else(read_status == 'na_above', limit_pg02mL, `pg/0.2mL`),
           `pg/hr` = if_else(read_status == 'na_above', limit_pghr, `pg/hr`),
           `pg/min` = if_else(read_status == 'na_above', limit_pgmin, `pg/min`)) %>%
    mutate(status = if_else(`pg/mL` >= arg_glc_ctrl, 'hit', 'no_hit'))

  only_hits <- final_df %>%
    filter(status == 'hit')

  top_10_perc <- quantile(only_hits$`pg/mL`, probs=0.9, na.rm = T)

  final_df %>%
    mutate(status = if_else(`pg/mL` >= top_10_perc, 'super_hit', status))


}
#'---compiles all tibbles created into list
df_parse <- function(df) {

  df_numeric <- data_to_numeric(df)

  ctrls_tbls <- ctrls_data(df_numeric)

  ctrls_avg_plate <- ctrls_tbls[[2]]
  ctrls_avg <- ctrls_tbls[[4]]

  complete_df <- final_df(df_numeric, ctrls_avg_plate, ctrls_avg)

  c(ctrls_tbls, list(complete_df))

}



#'---sets axis limits
ylim_func <- function(curr_df, orig_df, conc) {
  if(conc == 'pg/mL') {
    ylim(0,max_pull(orig_df, conc))
  }
  else if(conc == 'pg/0.2mL') {
    ylim(0,max_pull(orig_df, conc))
  }
  else if(conc == 'pg/hr') {
    ylim(0,max_pull(orig_df, conc))
  }
  else {
    ylim(0,max_pull(orig_df, conc))
  }
}

#'---pulls ctrl values for vertical lines
ctrl_mean <- function(df, conc, level){
  df %>%
    filter(Compound==level) %>%
    pull(!! sym(conc))
}

#'---creates bar plots for all compounds
all_compounds_plot <- function(df, ctrls_avg, comp_types, conc) {
  low_glc_ctrl <- ctrl_mean(ctrls_avg, conc, 'G 1.7')
  high_glc_ctrl <- ctrl_mean(ctrls_avg, conc, 'G 16.7')
  arg_glc_ctrl <- ctrl_mean(ctrls_avg, conc, 'G 1.7 + Arg 20')

  transpose_df <- t(df)
  transpose_df <- as.data.frame(transpose_df)
  rev_df <- rev(transpose_df)
  rev_df <- t(rev_df)
  rev_df <- as.data.frame(rev_df)

  rev_df %>%
    mutate_at(vars(plate_numb:`pg/min`), ~as.numeric(.)) %>%
    ungroup() %>%
    ggplot(aes(x=fct_inorder(Compound), y = !! sym(conc), fill = Pathway)) +
    geom_col() +
    scale_fill_manual(values=as.vector(polychrome(25))) +
    theme_light() +
    ylim_func(orig_df = df, conc=conc) +
    coord_flip() +
    labs(x = paste(str_to_title(comp_types),'Compound')) +
    guides(fill = guide_legend(ncol = 1)) +
    facet_wrap(vars(plate_numb), nrow = 1, scales = 'free') +
    scale_x_discrete(labels = function(x) str_wrap(x, width=15)) +
    theme(axis.text = element_text(size = 6),
          legend.position = c(0.95, 0.7),
          legend.text = element_text(size = 6)) +
    geom_hline(aes(yintercept = low_glc_ctrl, linetype = 'G 1.7'),
               size = 0.2, color = 'blue') +
    geom_hline(aes(yintercept = high_glc_ctrl, linetype = 'G 16.7'),
               size = 0.2, color = 'red') +
    geom_hline(aes(yintercept = arg_glc_ctrl, linetype = 'G 1.7 + Arg 20'),
               size = 0.2, color = 'orange1') +
    scale_linetype_manual(name = 'Control Means', values = c(2,2,2),
                          guide = guide_legend(override.aes = list(color = c('blue', 'orange1', 'red'))))
}

#'---creates pdf of bar plots
bar_plot_pdf <- function(final_df, ctrls_avg, comp_types) {

  std_err_key = c('pg/mL', "pg/0.2mL", "pg/hr", "pg/min")

  date = format(Sys.Date(), format = '%Y%m%d')

  group_plots = map(std_err_key, ~all_compounds_plot(final_df, ctrls_avg, comp_types, .x))
  group_plots = lapply(group_plots, ggplotGrob)

  if(comp_types == 'ion'){
    ggsave(paste0(date,'-',comp_types,'-all_compounds.pdf'),
           marrangeGrob(group_plots, nrow = 1, ncol = 1),
           path = '../plots/',
           device = 'pdf',
           width = 18,
           height = 20,
           units = 'in')
  }
  else{
    ggsave(paste0(date,'-',comp_types,'-all_compounds.pdf'),
           marrangeGrob(group_plots, nrow = 1, ncol = 1),
           path = '../plots/',
           device = 'pdf',
           width = 35,
           height = 20,
           units = 'in')
  }
}

#'---creates box and bar plots for all controls
controls_bar <- function(df, tbl, conc, conc_st_err, color) {
  df %>%
    ggplot(aes(x=Compound, y=!! sym(conc))) +
    geom_col(fill = color) +
    ylim_func(orig_df = tbl, conc = conc) +
    geom_errorbar(aes(ymin = !! sym(conc) - !! sym(conc_st_err),
                      ymax = !! sym(conc) + !! sym(conc_st_err)),
                  width = 0.2) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45,
                                     vjust = 0.95,
                                     hjust = 0.95))
}

#'---creates boxplot of select control
controls_box <- function(df, conc, color, lim) {
  df %>%
    ggplot(aes(x=Compound, y=!! sym(conc))) +
    geom_boxplot(outlier.shape = NA, fill = color) +
    geom_jitter(width = 0.2) +
    ylim_func(orig_df = df, conc = conc) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45,
                                     vjust = 0.95,
                                     hjust = 0.95))
}

#'---creates bar chart of select control
box_bar_plots <- function(ctrls_avg, ctrls_tbl, comp_types) {
  date = format(Sys.Date(), format = '%Y%m%d')

  pg_mL_bar <- controls_bar(ctrls_avg, ctrls_tbl, 'pg/mL', "pg/mL_std_err", 'dodgerblue3')
  pg_mL_box <- controls_box(ctrls_tbl, 'pg/mL', 'dodgerblue3')
  pg_mL02_bar <- controls_bar(ctrls_avg, ctrls_tbl,'pg/0.2mL', "pg/0.2mL_std_err", 'springgreen3')
  pg_mL02_box <- controls_box(ctrls_tbl, 'pg/0.2mL', 'springgreen3')
  pg_hr_bar <- controls_bar(ctrls_avg, ctrls_tbl, 'pg/hr', "pg/hr_std_err", 'salmon3')
  pg_hr_box <- controls_box(ctrls_tbl, 'pg/hr', 'salmon3')
  pg_min_bar <- controls_bar(ctrls_avg, ctrls_tbl, 'pg/min', "pg/min_std_err", 'gray50')
  pg_min_box <- controls_box(ctrls_tbl, 'pg/min', 'gray50')

  comb_bar_box <- grid.arrange(
    pg_mL_bar, pg_mL_box, pg_mL02_bar, pg_mL02_box,
    pg_hr_bar, pg_hr_box, pg_min_bar, pg_min_box, nrow = 2
  )
  ggsave(paste0(date, '-', comp_types, '-control_comps.pdf'),
         comb_bar_box,
         path = '../plots/',
         device = 'pdf',
         width = 10)
}

