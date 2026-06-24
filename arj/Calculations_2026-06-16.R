export <- FALSE
multip <- 1.96
# Data load ---------------------------------------------------------------
library(tidyverse)
library(parallel)
theme_set(
    theme_bw() +
        theme(
            text = element_text(family = "serif", size = 15),
            legend.position = "bottom"
        )
)
filepath <- dir("data", pattern = "xlsx") %>% 
    sort(decreasing = TRUE) %>% 
    `[`(1) %>% 
    paste0("data/", .)

long <- filepath %>% 
    readxl::read_excel(sheet = "main_data") %>% 
    filter(taxa != "no_insects") %>% 
    select(-age:-end, -duration, -traps, -family, -genus, -species) %>% 
    arrange(desc(site), taxa) %>% # -taxa
    mutate(tur = factor(tur), 
        zone = factor(zone, 
            levels = c("background", "bufer", "impact", "industrial_barren")),
        km = str_extract(site, "[:digit:]{1,}"),
        km = as.numeric(km), 
        year = as.factor(year), 
        .after = "site")

# 1 = turs 1 and 2 are not united
wide <- long %>% 
    select(-num) %>% 
    arrange(taxa) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0) %>% 
    select(-`_no_data`)

nsp_100 <- wide %>% 
    # select(-family) %>% 
    unite("id", year, tur, zone, km, site, plot, sep = "_") %>% 
    column_to_rownames("id") %>% 
    t %>% 
    as.data.frame() %>% 
    mclapply(
        function(x){
        x <- x[x>0]
        x <- as.numeric(as.character(x))
        if(length(x) == 0){0} else 
            if(length(x) == 1){1} else {
                iNEXT::iNEXT(x, size = 100, se = FALSE, q = 0) %>%
                    pluck("iNextEst", "size_based") %>%
                    filter(m == 100) %>%
                    pull(qD)
            }
        },
        mc.cores = detectCores()-1
    ) %>% 
    flatten_dbl()

# 1 = turs 1 and 2 are united
div <- tibble(wide[,1:6], 
               abu = apply(wide[,7:ncol(wide)], 1, sum),
               nsp = apply(wide[,7:ncol(wide)], 1, function(a){length(a[a>0])}),
               nsp100 = nsp_100,
               shan= vegan::diversity(wide[,7:ncol(wide)], 
                                      MARGIN = 1, index = "shannon", base = exp(1))
) %>% 
    mutate(km2 = km^2, kmLog = log(km), .after = km) %>% 
    mutate(abuLog = log10(abu+1), .after = abu)

W <- readxl::read_excel(filepath, "weather")
rm(nsp_100)
res <- list()
plots <- list()
tables <- list()

# Models preset ---------------------------------------------------------
models_fit <- function(y){
    d <- div
    if("tur"  %in% colnames(y)){d <- filter(d, tur  == y$tur)}
    if("year" %in% colnames(y)){d <- filter(d, year == y$year)}
    if("tur"  %in% colnames(y) & "year" %in% colnames(y)){set.seed(1); d <- mutate_if(d, is.numeric, ~`+`(.x, rnorm(nrow(d), sd = 1e-5)))}
    # if("tur"  %in% colnames(y) & "year" %in% colnames(y)){set.seed(1); d <- mutate_at(d, c("abu", "abuLog", "nsp", "nsp100", "shan"), ~`+`(.x, rnorm(nrow(d), sd = 1e-5)))}
    if(str_detect(y$formula, "Segmented")){
        ###
        segmented::segmented(lm(
            str_replace_all(y$formula, "Segmented", ""), 
            data = d), seg.Z = ~km)
        ###
    } else {
        lm(formula = y$formula, data = d)
    }
}

models_pred <- function(x){
    d <- c(if("segmented" %in% class(x)){x$psi[2]}else{numeric()}, seq(1, 32, by = 0.5)) %>% 
        sort() %>% 
        unique() %>%
        tibble(km = ., 
            km2 = km^2, 
            kmLog = log(km)
        )
    if(any(str_detect(names(x$coefficients), "tur"))){
        d <- expand_grid(d, tur = factor(c("1", "2"))) 
    }
    if(any(str_detect(names(x$coefficients), "year"))){
        d <- expand_grid(d, year = factor(c("2009", "2014"))) 
    }
    d %>% 
        mutate(predicted = predict(x, d)) %>% 
        select(any_of(c("km", "year", "tur", "predicted")))
}

manual_aic <- function (a, df_observed = NA, k = 2, simple = TRUE){
    # derived from "AIC.default
    lls <- logLik(a)
    # lls is logLik observer
    if(is.na(df_observed)){
        df_observed <- attr(lls, "df")
    }
    if(simple){
        -2 * as.numeric(lls) + k * df_observed
    } else {
        data.frame(
            df = df_observed, 
            AIC = -2 * as.numeric(lls) + k * df_observed)
    }
}

model_viz_supp <- function(df0){
    df1 <- df0 %>% 
        unite("id", formula, year, tur, sep = " ~ ") %>% 
        pull(id) %>% 
        `names<-`(df0$pred, .) %>% 
        map(~select(.x, km, predicted)) %>% 
        map_dfr(rbind, .id = "model") %>% 
        separate(model, into = c("response", "model", "year", "tur"), sep = " ~ ")
    
    df2 <- div %>% 
        select_at(c("km", predicted = df1$response[1], "year", "tur")) %>% 
        mutate(model = NA)
    
    if(any(str_detect(df0$formula, "Log"))){
        df2$predicted <- 10^df2$predicted
        df1$predicted <- 10^df1$predicted
    }
    
    df1 %>% 
        ggplot(aes(x = km, y = predicted,  linetype = model)) + #
        geom_point(shape = 21, size = 3, data = df2, color = "black", alpha = 0.5) +
        geom_line() + 
        # scale_y_log10() + 
        facet_grid(rows = vars(year), cols = vars(tur), scales = "free_y")  +
        theme(panel.grid = element_blank())
        # guides(linetype = "none")
}

# models count ------------------------------------------------------------
res$models <- expand_grid(
        resp = c("abu", "abuLog", "nsp", "nsp100", "shan"),
        year = c("2009", "2014"),
        # tur = c("1", "2"),
        formula = c("km + tur", "kmSegmented + tur", "km + km2 + tur", "kmLog + tur")
    ) %>% 
    mutate(formula = paste0(resp, " ~ ", formula))

res$fits <- res$models %>% 
    split(1:nrow(.)) %>% 
    lapply(models_fit)

res$pred <- lapply(res$fits, models_pred)

res$models <- res$models %>% 
    mutate(
        fits = unname(res$fits),
        pred = unname(res$pred),
        r2 = map_dbl(res$fits, ~summary(.x)$adj.r.squared),
        aic = map_dbl(res$fits, ~AIC(.x)),
        aic2 = map_dbl(res$fits, ~manual_aic(.x, df_observed = 3)),
        aic2 = case_when(str_detect(formula, "Segm") ~ aic2, TRUE ~ aic)
    ) %>% 
    arrange(resp, aic2)

# models viz --------------------------------------------------------------
df1 <- res$models %>%
    filter(str_detect(formula, "Segmented"), resp != "abu") %>% 
    select(resp, any_of(c("year", "tur", "yeartur")), pred) %>% 
    unnest(pred) %>% 
    mutate(
        yeartur = paste0(year, "_", tur),
        predicted = case_when(
            resp == "abuLog" ~ 10^predicted - 1, 
            TRUE ~ predicted)
    )
df2 <- div %>% 
    transmute(
        year = as.character(year), 
        tur = as.character(tur), 
        yeartur = paste0(year, "_", tur),
        km, abuLog, nsp, nsp100, shan) %>% 
    pivot_longer(names_to = "resp", values_to = "predicted", -km:-year)  %>% 
    mutate(
        predicted = case_when(
            resp == "abuLog" ~ 10^predicted - 1, 
            TRUE ~ predicted)
    )

ggplot(mapping = aes(
    x = km, y = predicted, 
    linetype = year, 
    shape = year,
    fill = year, color = year)) +
    geom_point(size = 2, data = df2, alpha = 0.5, fill = NA) + # color = "black"
    geom_line(data = df1) +
    facet_grid(
        rows = vars(resp), 
        cols = vars(tur),
        labeller = as_labeller(c(
            "1" = "I", "2" = "II",
            "abuLog" = "Обилие, \nэкз. / 100 лов.-сут.",
            "nsp"     = "Число видов \nнаблюдаемое",
            "nsp100"  = "Чисто видов \nна 100 особей", 
            "shan"    = "Индекс Шеннона\nH'"
        )),
        switch = "y",
        scales = "free") +
    scale_shape_manual(values = c(21, 24)) +
    scale_color_manual(values = c("darkgrey", "black")) +
    scale_x_continuous(limits = c(0.1, 32)) +
    theme(
        panel.grid = element_blank(),
        strip.placement = "outside",            # Выравнивает подписи за пределами оси (опционально)
        strip.background = element_blank()
    ) + 
    labs(y = NULL, x = "Расстояние, км", 
         shape = "Год: ", color = "Год: ", linetype = "Год: ")
if(!export){plots$f2_models.selected}




plots$s2_models.all <- map2(models[-1], labs, ~model_viz_supp(.x)+.y)

# Fig. S2. Models all
labs <- list(
    labs(x = NULL, y = "Log10 of individuals per 100 traps-days", 
         subtitle = "log Abundance"), 
    labs(x = NULL, y = "Number of species",
         subtitle = "Number of species"), 
    labs(x = NULL, y = "Number of species per 100 individuals",
         subtitle = "Number of species rarefied to 100 individuals"),
    labs(x = NULL, y = "Shannon index", subtitle = "Diversity")  
)

plots$s2_models.all <- gridExtra::grid.arrange(
    plots$s2_models.all[[1]],
    plots$s2_models.all[[2]], 
    plots$s2_models.all[[3]], 
    plots$s2_models.all[[4]], 
    ncol = 4
)

# Fig. 2. Models selected
# df1 <- models_separated %>% 
#     map_dfr(rbind) %>% 
#     unite("id", formula, year, tur, sep = " ~ ") %>% 
#     pull(id)  %>% 
#     `names<-`(pull(map_dfr(models_separated, ~select(.x, pred)), 1), .) %>% 
#     map_dfr(rbind, .id = "model") %>% 
#     separate(model, into = c("response", "model", "year"), sep = " ~ ") %>% 
#     filter(model == 'kmSegmented', response != "abu") %>% 
#     select(-km2, -kmLog, -model)


#     slice(1) %>% pull(pred)
# %>% select(pred) 
#     map(~.x %>% filter(str_detect(formula, "Segmented")) %>% select(pred) %>% ) %>% 
#     mutate(time_id = paste0(tur, "t_", year))
    




# + Table S3. Models all
tables$ts3_models.all <- res$models %>% 
    select(formula, year, aic = aic2, r2) %>% 
    separate(formula, into = c("predictor", "model"), sep = " ~ ")


    # filter(str_detect(formula, "abu ~ ", negate = TRUE)) %>%
    # map(~dplyr::select(.x, formula, year, aic = aic2)) %>%
    # map_dfr(~mutate_if(.x, is.numeric, ~round(.x, 2))) %>% 
    # filter(str_detect(formula, "abu ~ ", negate = TRUE)) %>% 
   
    # pivot_wider(names_from = predictor, values_from = aic) %>% 
    # arrange(year, model)
if(!export){tables$ts3_models.all }

# Table 2. Models selected comparison
mod_selected <- res$models %>% 
    filter(res$models$resp != "abu" & str_detect(res$models$formula, "Segmented")) %>% 
    unite("id", formula, year, sep = "_")

mod_selected %>% 
    pull(fits) %>% 
    `names<-`(str_remove_all(mod_selected$id, "Segmented \\+ tur" )) %>% 
    map(~.x %>% 
            summary %>% 
            capture.output() %>% 
            as.data.frame()
        ) %>% 
    writexl::write_xlsx("models.selected.xlsx")

tables$t2_mod.comps <- res$models$fits %>% 
    `[`(-1) %>% 
    map(~.x %>% 
            filter(str_detect(formula, "Segment")) %>% 
            select(year, fits)
    ) %>% 
    map(~.x %>% 
            pluck("fits") %>% 
            `names<-`(.x$year) %>% 
            lapply(summary) %>% 
            map(~list(.x$coefficients, .x$psi)) %>% 
            map(~map(.x, ~rownames_to_column(as.data.frame(.x), "pred"))) %>% 
            map_dfr(~rbind(
                select(.x[[1]], pred, `Est.` = 2, `St.Err` = 3), 
                transmute(.x[[2]], pred, `Est.`, `St.Err`)
            ) %>% 
                filter(Est. != 0, pred %in% c("km", "psi1.km")), 
            .id = "year"
            ) %>% 
            split(.$pred) %>% 
            map(~arrange(.x, Est.))
    ) %>% 
    map_dfr(~.x %>% 
                map_dfr(rbind) %>% 
                mutate(E = paste0(round(Est., 3), "±", round(St.Err, 3)), 
                       low_CI = Est. - St.Err * multip, 
                       up_CI  = Est. + St.Err * multip, 
                       CI = paste0("(", round(low_CI, 3), "; ", round(up_CI, 3), ")"),
                       .keep = "unused") %>% 
                select(-ends_with("_CI")) %>% 
                unite("E", E, CI, sep = " || ") %>% 
                arrange(year) %>% 
                pivot_wider(names_from = year, values_from = E),
            .id = "response"
    )
if(!export){tables$t2_mod.comps}


# weather count -----------------------------------------------------------
W %>% 
    mutate(tur = case_when(
        # 2009 г. с 5 по 10 июня; 2014 г. с 3 по 8 июня
        # 2009 г. с 26 по 31 августа; 2014 г. с 18 по 23 августа)
        year == 2009 & month == 6 & day >=  5 & day <= 10 ~ "tur_I",
        year == 2014 & month == 6 & day >=  3 & day <=  8 ~ "tur_I",
        year == 2009 & month == 8 & day >= 26 & day <= 31 ~ "tur_II",
        year == 2014 & month == 8 & day >= 18 & day <= 23 ~ "tur_II"
    )) %>% 
    filter_out(is.na(tur)) %>% 
    mutate(month = tur, .keep = "unused") %>% 
    rbind(mutate(W, month = paste0("month_", month))) %>% 
    # unite("dtm", year, month, sep = "__") %>% 
    select(-day) %>% 
    pivot_longer(names_to = "var", values_to = "val", -1:-2) %>% 
    group_by(year, month, var) %>% 
    summarise(mean = mean(val), sd = sd(val), ci = sd*multip, .groups = "drop") %>% 
    mutate_if(is.numeric, ~round(.x, 1)) %>% 
    mutate(
        i = paste0(mean, "±", sd, " (", round(mean-ci, 1), "...", round(mean+ci, 1), ")"), 
        .keep = "unused") %>% 
    pivot_wider(names_from = var, values_from = i) # %>% writexl::write_xlsx("weater_tmp.xlsx")
    

# Export ------------------------------------------------------------
# Numbers in text 
tables$t0_text <- div %>% 
    group_by(site) %>% 
    summarise(
        mean_abu = mean(abu), 
        min_abu  = min(abu), 
        max_abu  = max(abu),
        .groups = "drop") %>% 
    arrange(mean_abu)
if(!export){tables$t0_text} else {
    # ...
}



if(export){
# Fig. S2. Models all
ggsave(
    plot = plots$s2_models.all, 
    filename = paste0("export/Fig.S2_models.all_", Sys.Date(), ".pdf"), 
    height = 6, width = 11, dpi = 600)
    
# Fig. 2. Models selected
ggsave(
    plot = plots$f2_models.selected, 
    filename = paste0("export/Fig.2_models.selected_", Sys.Date(), ".pdf"), 
    height = 9, width = 11, dpi = 600)

# + Table S3. Models all
    writexl::write_xlsx(
        tables$ts3_models.all,
        paste0("export/Tab.S2_models.selected_", Sys.Date(), ".xlsx")
    )

# Table 2. Models selected
writexl::write_xlsx(tables$t2_mod.comps,  paste0("export/Tab.2_models.comp_", Sys.Date(), ".xlsx"))
# + Fig. 3. Rarefied number of species
ggsave(
    paste0("export/Fig.3_rarefication_", Sys.Date(), ".pdf"), 
    plots$f3_rarefication, 
    width = 9, height = 5.5, dpi = 600)

# Table 3. Dominant species

# Fig. 4. PCoA ordination
ggsave(paste0("export/Fig.3_ord_", Sys.Date(), ".pdf"), 
       plot = plots$f4_pcoa,
       width = 18, height = 13, units = "cm")

# Permanova in text
res$permanova %>% 
    capture.output() %>% 
    `[`(-c(2, 3, 5, 6, 13, 14)) %>% 
    as.data.frame() %>% 
    writexl::write_xlsx(paste0("export/Tab.4a_permanova_", Sys.Date(), ".xlsx"))

# Table 4. Average distances
if(export){
    distances %>% 
        mutate_if(is.numeric, function(a){round(a, 2)}) %>% 
        writexl::write_xlsx(paste0("export/Tab.4_distances_", Sys.Date(), ".xlsx"))
}

}