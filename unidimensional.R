export = F
multip = 2
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

wide <- long %>% 
    select(-num, -tur) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0)

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
div <- tibble(wide[,1:5], 
              abu = apply(wide[,6:ncol(wide)], 1, sum),
              nsp = apply(wide[,6:ncol(wide)], 1, function(a){length(a[a>0])}),
              nsp30 = nsp_30,
              shan= vegan::diversity(wide[,6:ncol(wide)], 
                                     MARGIN = 1, index = "shannon", base = exp(1))
) %>% 
    mutate(km2 = km^2, kmLog = log(km), .after = km) %>% 
    mutate(abuLog = log10(abu+1), .after = abu)

rm(nsp_30)
res <- list()
plots <- list()
tables <- list()

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
                    data = df), seg.Z = ~km, 
                    control = segmented::seg.control(fix.npsi = F))
                ###
            } else {
                lm(formula = y$formula, data = df)
            }
        })
}

models_pred <- function(fits){
    lapply(fits, function(x){
        c(if("segmented" %in% class(x)){x$psi[2]}else{numeric()}, seq(1, 32, by = 0.5)) %>% 
            sort() %>% 
            unique() %>%
            tibble(
                km = .,
                km2 = .^2, 
                kmLog = log(.)
            ) %>% 
            mutate(predicted = predict(x, .))
    }
    )}

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

model_viz_supp <- function(df0){
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

# models count ------------------------------------------------------------
models <- expand_grid(
    resp = c("abu", "nsp", "nsp30", "shan"), # "abuLog", 
    formula = c("km", "kmSegmented", "km + km2", "kmLog"), 
    year = c(2009, 2014)) %>% 
    mutate(
        formula = paste0(resp, " ~ ", formula)) %>% 
    # slice(-11) %>% 
    mutate(fits = models_fit(.), 
           pred = models_pred(fits),
           aic = map_dbl(fits, ~AIC(.x)), 
           r2 = map_dbl(fits, ~summary(.x)$adj.r.squared)
    ) %>% 
    mutate(aic2 = map_dbl(fits, ~manual_aic(.x, df_observed = 3)), 
           aic2 = case_when(str_detect(formula, "Segm") ~ aic2, TRUE ~ aic),
           .after = aic) %>% 
    arrange(resp, year, aic2) %>% 
    split(.$resp)


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
# plots$s2_models.all <- map2(models[-1], my_labs, ~model_viz_supp(.x)+.y)

plots$s2_models.all <- gridExtra::grid.arrange(
    plots$s2_models.all[[1]],
    plots$s2_models.all[[2]], 
    plots$s2_models.all[[3]], 
    plots$s2_models.all[[4]], 
    ncol = 4
)


# Models selection --------------------------------------------------------
models %>% 
    map_dfr(rbind) %>% 
    unite("id", resp, year) %>% 
    split(.$id) %>% 
    map(~arrange(select(.x, formula, aic = aic2, r2), aic))

# + Table S3. Models all
tables$ts3_models.all <- models %>% 
    map(~select(.x, formula, year, aic = aic2)) %>%
    map_dfr(~mutate_if(.x, is.numeric, ~round(.x, 2))) %>% 
    # filter(str_detect(formula, "abu ~ ", negate = TRUE)) %>% 
    separate(formula, into = c("predictor", "model"), sep = " ~ ") %>% 
    pivot_wider(names_from = predictor, values_from = aic) %>% 
    arrange(year, model)
if(!export){tables$ts3_models.all }


# Models viz selected  ----------------------------------------------------
# Fig. 2. Models selected
df1 <- models %>% 
    map_dfr(rbind) %>% 
    unite("id", formula, year, sep = " ~ ") %>% 
    pull(id)  %>% 
    `names<-`(pull(map_dfr(models, ~select(.x, pred)), 1), .) %>% 
    map_dfr(rbind, .id = "model") %>% 
    separate(model, into = c("response", "model", "year"), sep = " ~ ") %>% 
    filter(model == 'kmSegmented') %>% # , response != "abu"
    select(-km2, -kmLog, -model)

df2 <- div %>% 
    transmute(year = as.character(year), km, abu, nsp, nsp30, shan) %>% 
    pivot_longer(names_to = "response", values_to = "predicted", -km:-year) 
plots$f2_models.selected <- df1 %>% 
    ggplot(aes(x = km, y = predicted, fill = year, color = year, shape = year)) +
    geom_point(size = 3, data = df2, alpha = 0.5, color = "black") +
    geom_line() +
    facet_wrap(~response, scales = "free") +
    scale_shape_manual(values = c(21, 24)) + 
    scale_x_continuous(limits = c(0.1, 32)) +
    theme(panel.grid = element_blank()) +
    labs(y = NULL, x = "Distanse, km", shape = "Year: ", color = "Year: ", fill = "Year: ")
if(!export){plots$f2_models.selected}


# Models param selected ---------------------------------------------------
# Table 2. Models selected comparison
tables$t2_mod.comps <- models %>% 
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

# export ------------------------------------------------------------------


# in progress...