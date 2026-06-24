export <- TRUE
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

# models viz: manuscript Fig 2  ------------------------------------------------
lns <- res$models %>%
    filter(str_detect(formula, "Segmented"), resp != "abu") %>% 
    select(resp, any_of(c("year", "tur", "yeartur")), pred) %>% 
    unnest(pred) %>% 
    mutate(
        yeartur = paste0(year, "_", tur),
        predicted = case_when(
            resp == "abuLog" ~ 10^predicted - 1, 
            TRUE ~ predicted)
    )
pts <- div %>% 
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

plots$fig2 <- ggplot(mapping = aes(
    x = km, y = predicted, 
    linetype = year, 
    shape = year,
    fill = year)) +
    geom_point(data = pts, size = 2, alpha = 0.5, fill = NA, color = "black") + # 
    geom_line(data = lns) +
    facet_grid(
        rows = vars(resp), 
        cols = vars(tur),
        labeller = as_labeller(c(
            "1" = "I", "2" = "II",
            "abuLog" = "A", #"Обилие, \nэкз. / 100 лов.-сут.",
            "nsp"    = "B", #"Число видов \nнаблюдаемое",
            "nsp100" = "C", #"Чисто видов \nна 100 особей", 
            "shan"   = "D" #"Индекс Шеннона\nH'"
        )),
        switch = "y",
        scales = "free") +
    scale_shape_manual(values = c(21, 24)) +
    # scale_color_manual(values = c("darkgrey", "black")) +
    scale_x_continuous(limits = c(0.1, 32)) +
    theme(
        panel.grid = element_blank(),
        strip.placement = "outside",            # Выравнивает подписи за пределами оси (опционально)
        strip.background = element_blank(),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    ) + 
    labs(y = NULL, x = "Расстояние, км", 
         shape = "Год: ", linetype = "Год: ") # color = "Год: ", 
if(!export){plots$fig2} else {
    ggsave(
        paste0("export/Fig_2_",  Sys.Date(), ".pdf"), 
        device = cairo_pdf,
        plots$fig2,
        height = 20*1.2,
        width = 15*1.2, 
        units = "cm"
    )
}
rm(pts, lns)

# models viz: manuscript Suppl 1  ----------------------------------------------
lns <- c("abuLog", "nsp", "nsp100", "shan") %>% 
    map_dfr(function(R){
        # R = "abuLog"
        lns <- res$models %>% 
            filter(resp == R) %>% 
            separate(formula, into = c("tmp", "type"), sep = " ~ | \\+ tur", extra = "drop") %>% 
            unnest(pred) %>% 
            select(year, tur, km, predicted, resp, type) %>% 
            mutate(yeartur = paste0(tur, " - ", year), .after = tur)
        if(any(str_detect(lns$resp, "Log"))){
            lns$predicted <- 10^lns$predicted-1
            # pts$predicted <- 10^pts$predicted-1
        }
        return(lns)
    })

pts <- c("abuLog", "nsp", "nsp100", "shan") %>% 
    map_dfr(function(R){
        # R = "abuLog"
        pts <- div %>%
            select(year, tur, km, predicted = all_of(R)) %>% 
            mutate(resp = R, type = NA) %>% 
            mutate(yeartur = paste0(tur, " - ", year), .after = tur)
        
        if(any(str_detect(pts$resp, "Log"))){
            pts$predicted <- 10^pts$predicted-1
        }
        return(pts)
    })

plots$supp1 <- ggplot(mapping = aes(x = km, y = predicted,  linetype = type)) + #
        geom_point(data = pts, shape = 21, size = 3, color = "black", alpha = 0.5) +
        geom_line(data = lns) + 
        facet_grid(
            cols = vars(yeartur),
            rows = vars(resp),
            scales = "free_y",
            switch = "y",
            labeller = as_labeller(c(
                "1 - 2009" = "I - 2009", 
                "2 - 2009" = "II - 2009", 
                "1 - 2014" = "I - 2014", 
                "2 - 2014" = "II - 2014",
                "abuLog" = "A",
                "nsp" = "B",
                "nsp100" = "C",
                "shan" = "D"
            )),
            )  +
        theme(
            strip.placement = "outside", 
            strip.background = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12)
        ) + 
    labs(y = NULL, x = "Расстояние, км")
    # guides(linetype = "none")

if(!export){plots$supp1} else {
    ggsave(
        paste0("export/Supp_1_",  Sys.Date(), ".pdf"), 
        device = cairo_pdf,
        plots$supp1,
        height = 20*1.2,
        width = 15*1.2, 
        units = "cm"
    )
}

# Tables ------------------------------------------------------------------
# Models comparison
tables$all_models <- res$models %>% 
    select(formula, year, aic = aic2, r2) %>% 
    separate(formula, into = c("predictor", "model"), sep = " ~ ")
if(export){
    writexl::write_xlsx(
        tables$all_models,
        paste0("export/models.all_", Sys.Date(), ".xlsx")
    )
}

# Modelsselected 
res$selected_models <- res$models %>% 
    filter(res$models$resp != "abu" & str_detect(res$models$formula, "Segmented")) %>% 
    unite("id", formula, year, sep = "_")

if(export){
    res$selected_models %>% 
        pull(fits) %>% 
        `names<-`(str_remove_all(res$selected_models$id, "Segmented \\+ tur" )) %>% 
        map(~.x %>% 
                summary %>% 
                capture.output() %>% 
                as.data.frame()
        ) %>% 
        writexl::write_xlsx(paste0("export/models.selected_", Sys.Date(), ".xlsx"))
}

# weather count -----------------------------------------------------------
res$weather_count <- W %>% 
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
    select(-day) %>% 
    pivot_longer(names_to = "var", values_to = "val", -1:-2) %>% 
    group_by(year, month, var) %>% 
    summarise(mean = mean(val), sd = sd(val), ci = sd*multip, .groups = "drop") %>% 
    mutate_if(is.numeric, ~round(.x, 1)) %>% 
    mutate(
        i = paste0(mean, "±", sd, " (", round(mean-ci, 1), "...", round(mean+ci, 1), ")"), 
        .keep = "unused") %>% 
    pivot_wider(names_from = var, values_from = i) 
    
if(!export){res$weather_count} else {
    writexl::write_xlsx(res$weather_count, paste0("export/weater_", Sys.Date(), ".xlsx"))
}

