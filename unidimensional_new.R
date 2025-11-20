export = T
multip = 2
turs   = "yes"
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

long <- dir("data", pattern = "xlsx") %>% 
    sort(decreasing = TRUE) %>% 
    `[`(1) %>% 
    paste0("data/", .) %>% 
    # преобразовать "sp *" "sp"
    # перед ординацией убирать всех sp (кроме пронумерованных)
    readxl::read_excel(sheet = 1) %>%  
    select(-start, -end) %>% 
    group_by(year, tur, zone, site, plot, taxa) %>% 
    summarise(num = sum(num, na.rm = TRUE), abu = sum(abu, na.rm = TRUE), .groups = "drop") %>% 
    arrange(desc(site), taxa) %>% 
    mutate(
        tur = factor(tur), 
        site = fct_inorder(site),
        zone = factor(zone, 
                      levels = c("background", "bufer", "impact", "industrial_barren", ordered = TRUE)),
        km = as.numeric(str_extract(site, "[:digit:]{1,}")),
        year = as.factor(year), 
        .after = "site")

if(turs == "yes"){
    wide <- long %>% 
        select(-num)
} else {
    wide <- long %>% 
        select(-num, -tur)
}
wide <- wide %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0) %>% 
    select(-`_no_data`)

nsp_30 <- wide %>% 
    select(-year:-plot) %>% 
    t %>% 
    as.data.frame() %>% 
    sapply(function(x){
        x <- x[x>0]
        x <- as.numeric(as.character(x))
        if(length(x) == 0){0} else 
            if(length(x) == 1){1} else {
                iNEXT::iNEXT(x, size = 30, se = FALSE, q = 0) %>%
                    pluck("iNextEst", "size_based") %>%
                    filter(m == 30) %>%
                    pull(qD)
            }
    })

# 1 = turs 1 and 2 are united
div <- wide %>% 
    select(year:plot) %>% 
    mutate(
        abu = apply(select(wide, -year:-plot), 1, sum),
        nsp = apply(select(wide, -year:-plot), 1, function(a){length(a[a>0])}),
        nsp30 = nsp_30,
        shan= vegan::diversity(select(wide, -year:-plot), 
                               MARGIN = 1, index = "shannon", base = exp(1))
    ) %>% 
    mutate(km2 = km^2, kmLog = log(km), .after = km) %>% 
    mutate(abuLog = log10(abu+1), .after = abu)

# div <- tibble(wide[,1:5], 
#               abu = apply(wide[,6:ncol(wide)], 1, sum),
#               nsp = apply(wide[,6:ncol(wide)], 1, function(a){length(a[a>0])}),
#               nsp30 = nsp_30,
#               shan= vegan::diversity(wide[,6:ncol(wide)], 
#                                      MARGIN = 1, index = "shannon", base = exp(1))
# ) %>% 
#     mutate(km2 = km^2, kmLog = log(km), .after = km) %>% 
#     mutate(abuLog = log10(abu+1), .after = abu)

rm(nsp_30)
res <- list()
plots <- list()
tables <- list()

cat("\nData loaded")
# Models functions --------------------------------------------------------
models_fit <- function(x){
    x %>% 
        split(1:nrow(.)) %>%
        lapply(function(y){
            df <-  filter(div, year == y$year)
            if(str_detect(y$formula, "Segmented")){
                ###
                
                
                segmented::segmented(lm(
                    str_replace_all(y$formula, "Segmented", ""), 
                    data = df), seg.Z = ~km)
                    # control = segmented::seg.control(fix.npsi = F))
                ###
            } else {
                lm(formula = y$formula, data = df)
            }
        })
}

models_pred_noturs <- function(fit){
        c(if("segmented" %in% class(fit)){fit$psi[2]}else{numeric()}, seq(1, 32, by = 0.5)) %>%
            sort() %>%
            unique() %>%
            tibble(
                km = .,
                km2 = .^2,
                kmLog = log(.)
            ) %>%
            mutate(predicted = predict(fit, .)) %>% 
        mutate(tur = "", .before = 1)
    }

models_pred_yestur <- function(fit){
    km_list <- c(if("segmented" %in% class(fit)){fit$psi[2]}else{numeric()}, seq(1, 32, by = 0.5)) %>% 
        sort() %>% 
        unique()
    expand_grid(tur = c(1, 2), km = km_list) %>% 
        mutate(
            tur = factor(tur), 
            km = km,
            km2 = km^2, 
            kmLog = log(km)) %>% 
            mutate(predicted = predict(fit, .))
    }

manual_aic <- function (a, df_observed = NA, k = 2, simple = TRUE){
    # derived from "AIC.default"    
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
    
    # old manual AIC
    # manual.aic.old <- function(a){
    #     if(!(a$family$family %in% c("quasipoisson", "poisson"))){return(NA)}
    #     loglik <- sum(dpois(a$model[,1], fitted.values(a), log = TRUE))
    #     phi <- summary(a)$dispersion
    #     -2 * loglik + 2 * summary(a)$df[3] * phi
    # }
}

model_viz_supp_notur <- function(df0){
    df1 <- df0 %>% 
        unite("id", formula, year, sep = " ~ ") %>% 
        pull(id) %>% 
        `names<-`(df0$pred, .) %>% 
        map(~select(.x, km, predicted)) %>% 
        map_dfr(rbind, .id = "model") %>% 
        separate(model, into = c("response", "model", "year"), sep = " ~ ")
    
    df2 <- div %>% 
        select_at(c("km", predicted = df1$response[1], "year")) %>% 
        mutate(model = NA)
    
    df1 %>% 
        ggplot(aes(x = km, y = predicted,  linetype = model)) + #
        geom_point(shape = 21, size = 3, data = df2, color = "black", alpha = 0.5) +
        geom_line() + 
        facet_wrap(~year, scales = "fixed", ncol = 1)  +
        theme(panel.grid = element_blank())
    # guides(linetype = "none")
}

model_viz_supp_tur <- function(df0){
    df0 <- filter(df0, str_detect(formula, "tur"))
    df1 <- df0 %>%
        unite("id", formula, year, sep = " ~ ") %>% 
        pull(id) %>% 
        str_remove_all(" \\+ tur") %>% 
        `names<-`(df0$pred, .) %>% 
        map(~select(.x, tur, km, predicted)) %>% 
        map_dfr(rbind, .id = "model") %>% 
        separate(model, into = c("response", "model", "year"), sep = " ~ ") %>% 
        mutate(tur = fct_relabel(tur, ~paste0("tur ", .x))) 
    
    df2 <- div %>% 
        select_at(c("km", "tur", predicted = df1$response[1], "year")) %>% 
        mutate(model = NA) %>% 
        mutate(tur = fct_relabel(tur, ~paste0("tur ", .x))) 
    
    df1 %>% 
        ggplot(aes(x = km, y = predicted,  color = model, linetype = model)) + #
        geom_point(shape = 21, size = 3, data = df2, color = "black", alpha = 0.5) +
        geom_line() + 
        facet_grid(
            cols = vars(year), 
            rows = vars(tur),
            scales = "fixed")  +
        theme(panel.grid = element_blank())
    # guides(linetype = "none")
}

cat("\nFunctions are ready")
# models count ------------------------------------------------------------
res$dummy <- expand_grid(
        resp = c("abu", "nsp", "nsp30", "shan"), # "abuLog", 
        formula = c("km", "kmSegmented", "km + km2", "kmLog"), 
        tur = switch (turs, yes = c(" + tur", ""), no = ""), 
        year = c(2009, 2014)) %>% 
    mutate(
        formula = paste0(resp, " ~ ", formula, tur)) 
res$fits <- models_fit(res$dummy)
res$predictions <- res$fits %>% 
    # `[`(6:8) %>% 
    lapply(., function(x){
        if(sum(str_detect(names(coefficients(x)), "tur")) > 0){
            models_pred_yestur(x)
        } else {
            models_pred_noturs(x)
        }
    })

res$models <- res$dummy %>% 
    mutate(
        fits = res$fits,
        pred = res$predictions,
        aic = map_dbl(fits, ~AIC(.x)), 
        aic2 = map_dbl(fits, ~manual_aic(.x, df_observed = 3)), # 3 -> ...?
        r2 = map_dbl(fits, ~summary(.x)$adj.r.squared)
    ) %>% 
    mutate(aic = case_when(str_detect(formula, "Segm") ~ aic2, TRUE ~ aic)) %>% 
    select(-aic2) %>% 
    arrange(resp, year, aic) %>% 
    split(.$resp)

if(turs == "yes"){
    res$models <- map(
        res$models, 
        ~filter(.x, str_detect(formula, "tur")))
} else {
    res$models <- map(
        res$models, 
        ~filter(.x, str_detect(formula, "tur", negate = TRUE)))
}

cat("\nModels are calculated")
# Models viz all --------------------------------------------------------------

# Fig. S2. Models all
my_labs <- list(
    labs(x = NULL, y = "Individuals per 100 traps-days", 
         subtitle = "Abundance"), 
    labs(x = NULL, y = "Number of species",
         subtitle = "Number of species"), 
    labs(x = NULL, y = "Number of species per 30 individuals",
         subtitle = "Number of species rarefied to 30 individuals"),
    labs(x = NULL, y = "Shannon index", subtitle = "Diversity")  
)

if(turs == "yes"){
    plots$s2_models.all <- map2(res$models, my_labs, ~model_viz_supp_tur(.x)+.y)
} else {
    plots$s2_models.all <- map2(res$models, my_labs, ~model_viz_supp_notur(.x)+.y)
}

plots$s2_models.all <- gridExtra::grid.arrange(
    plots$s2_models.all[[1]],
    plots$s2_models.all[[2]], 
    plots$s2_models.all[[3]], 
    plots$s2_models.all[[4]], 
    ncol = 4
)
cat("\nFigure 1 done")
# Models selection --------------------------------------------------------
tables$models.all <- res$models %>% 
    map_dfr(rbind) %>% 
    unite("id", resp, year) %>% 
    split(.$id) %>% 
    map(~arrange(select(.x, formula, aic, r2), aic))

# + Table S3. Models all
# tables$ts3_models.all <- res$models %>% 
#     map(~select(.x, formula, year, aic)) %>%
#     map_dfr(~mutate_if(.x, is.numeric, ~round(.x, 2))) %>% 
#     # filter(str_detect(formula, "abu ~ ", negate = TRUE)) %>% 
#     separate(formula, into = c("predictor", "model"), sep = " ~ ") %>% 
#     pivot_wider(names_from = predictor, values_from = aic) %>% 
#     arrange(year, model)
# if(!export){tables$ts3_models.all }

cat("\nTable 1 done")
# Models viz selected  ----------------------------------------------------
if(turs == "yes"){
    plots$models.selected <- res$models %>% 
        map(~filter(.x, str_detect(formula, "Segm"))) %>% 
        map2(my_labs, ~model_viz_supp_tur(.x)+.y)
} else {
    plots$models.selected <- res$models %>% 
        map(~filter(.x, str_detect(formula, "Segm"))) %>% 
        map2(my_labs, ~model_viz_supp_notur(.x)+.y)
}

plots$models.selected <- gridExtra::grid.arrange(
    plots$models.selected[[1]],
    plots$models.selected[[2]], 
    plots$models.selected[[3]], 
    plots$models.selected[[4]], 
    ncol = 4
)

# Fig. 2. Models selected
# df1 <- res$models %>% 
#     map_dfr(rbind) %>% 
#     unite("id", formula, year, sep = " ~ ") %>% 
#     pull(id)  %>% 
#     `names<-`(pull(map_dfr(res$models, ~select(.x, pred)), 1), .) %>% 
#     map_dfr(rbind, .id = "model") %>% 
#     separate(model, into = c("response", "model", "year"), sep = " ~ ") %>% 
#     filter(str_detect(model, 'kmSegmented')) %>% # , response != "abu"
#     select(-km2, -kmLog, -model)

# df2 <- div %>% 
#     transmute(year = as.character(year), km, abu, nsp, nsp30, shan) %>% 
#     pivot_longer(names_to = "response", values_to = "predicted", -km:-year) 
# plots$f2_models.selected <- df1 %>% 
#     ggplot(aes(x = km, y = predicted, fill = year, color = year, shape = year)) +
#     geom_point(size = 3, data = df2, alpha = 0.5, color = "black") +
#     geom_line() +
#     facet_wrap(~response, scales = "free") +
#     scale_shape_manual(values = c(21, 24)) + 
#     scale_x_continuous(limits = c(0.1, 32)) +
#     theme(panel.grid = element_blank()) +
#     labs(y = NULL, x = "Distanse, km", shape = "Year: ", color = "Year: ", fill = "Year: ")
# if(!export){plots$f2_models.selected}

cat("\nFigure 2 done")
# Models param selected ---------------------------------------------------
# Table 2. Models selected comparison
tables$t2_mod.comps <- res$models %>% 
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
# if(!export){tables$t2_mod.comps}
cat("\nTable 2 done")
# export ------------------------------------------------------------------
cat("\nExport in progress...")
if(export){
    if(turs == "yes"){
        turs <- "yes.turs"
    } else {
        turs <- "no.turs"
    }
    writexl::write_xlsx(
        tables$models.all, 
        paste0("export/models_all_", turs, "_", Sys.Date(), ".xlsx")
    )
    writexl::write_xlsx(
        tables$t2_mod.comps, 
        paste0("export/models_selected_", turs, "_", Sys.Date(), ".xlsx")
    )
    
    ggsave(
        paste0("export/models_all_", turs, "_", Sys.Date(), ".pdf"), 
        plots$s2_models.all, 
        height = 210, 
        width = 297*1.5, 
        units = "mm", 
        dpi = 900
    )
    
    ggsave(
        paste0("export/models_selected_", turs, "_", Sys.Date(), ".pdf"), 
        plots$models.selected, 
        height = 210, 
        width = 297*1.5, 
        units = "mm", 
        dpi = 900
    )
}
cat("\nExport is done\n")