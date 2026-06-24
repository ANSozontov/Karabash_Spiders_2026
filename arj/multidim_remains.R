# Rarefication ------------------------------------------------------------
library(parallel)
cl <- makeCluster(detectCores()-1)
rar <- long %>% 
    group_by(year, zone, taxa) %>% 
    summarise(abu = sum(abu), .groups = "drop") %>% 
    unite("id", year, zone, sep = "_") %>% 
    split(.$id) %>% 
    map(~sort(.x$abu[.x$abu>0])) %>% 
    parLapply(cl = cl, ., function(a){
        a |> 
            as.character() |>
            as.numeric() |>
            iNEXT::iNEXT(
                size = seq(0, 100, by = 2), 
                ### 999
                nboot = 9, 
                ### 999
                se = TRUE, 
                conf = 0.999) |>
            purrr::pluck("iNextEst", "size_based") |>
            dplyr::select(m, Method, qD, qD.LCL,   qD.UCL) 
    }) %>% 
    map_df(rbind, .id = "id") %>% 
    as_tibble() %>% 
    separate(id, into = c("year", "zone"), sep = "_") %>% 
    mutate(zone = factor(zone, ordered = TRUE, levels = levels(long$zone)))

plots$f3_rarefication <- rar %>% 
    filter(m %in% seq(0, 100, by = 2)) %>% 
    mutate(qD.LCL  = case_when(is.na(qD.LCL) ~ qD, TRUE ~ qD.LCL),
           qD.UCL  = case_when(is.na(qD.UCL) ~ qD, TRUE ~ qD.UCL)) %>% 
    ggplot(aes(x = m, y = qD, color = zone, fill = zone)) + 
    facet_wrap(~year) +
    geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.2, color = NA) +
    geom_line() + 
    scale_color_manual(values = c("darkgreen", "#ccff00", "#ff9900", "#FF3300")) +
    scale_fill_manual(values = c("darkgreen", "#ccff00", "#ff9900", "#FF3300")) +
    labs(x = "Individuals", y = "Species") + 
    theme(panel.grid = element_blank())

if(!export){plots$f3_rarefication}

# Multidimensional  --------------------------------------------------
dis <- wide %>% 
    select(-site) %>% 
    unite("ID", zone, year, km, plot, sep = "_") %>% 
    column_to_rownames("ID") %>% 
    vegan::vegdist(method = "bray", binary = FALSE)
pc <- ape::pcoa(dis)
eig <- pc$values$Eigenvalues
eig <- round(eig/sum(eig)*100, 1)
pc <- pc$vectors %>% 
    as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    as_tibble() %>% 
    separate(ID, into = c("zone", "year", "site", "plot"), 
             sep = "_") %>% 
    mutate(zone = factor(zone, levels = levels(long$zone)))

L <- expand_grid(levels(long$zone), c(2009, 2014)) %>% 
    apply(1, function(a){paste0(a, collapse = "_")})

distances <- as.matrix(dis)
distances[upper.tri(distances)] <- NA
diag(distances) <- NA
distances <- distances %>%
    as.data.frame() %>%
    rownames_to_column("id1") %>%
    as_tibble() %>%
    pivot_longer(names_to = "id2", values_to = "dis", -id1) %>%
    filter(!is.na(dis)) %>%
    separate(id1, c("zone1", "year1"), sep = "_", extra = "drop") %>%
    separate(id2, c("zone2", "year2"), sep = "_", extra = "drop") %>%
    mutate(
        id1 = factor(paste0(zone1, "_", year1), levels = L, ordered = TRUE),
        id2 = factor(paste0(zone2, "_", year2), levels = L, ordered = TRUE),
        .keep = "unused") %>% 
    split(1:nrow(.)) %>% 
    lapply(function(x){
        if(x$id1 > x$id2) {
            tibble(dis = x$dis, id1 = x$id2, id2 = x$id1)
        } else {
            tibble(dis = x$dis, id1 = x$id1, id2 = x$id2)
        }
    }) %>% 
    map_dfr(rbind) %>% 
    pivot_wider(
        names_from = id2, values_from = dis, 
        values_fill = NA, values_fn = mean)

set.seed(1); res$permanova <- div %>% 
    vegan::adonis2(
        dis ~ zone * year, 
        data = ., 
        permutations = 999, 
        by = "terms")

if(!export){res$permanova} 

plots$f4_pcoa <- pc %>% 
    ggplot(
        aes(x = Axis.1, y = Axis.2, linetype = year,
            fill = zone, color = zone, shape = year)) +
    geom_point(color = "black", size = 2.5) +  
    stat_ellipse() + 
    scale_shape_manual(values = c(21, 22))+
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_fill_manual(values = c("darkgreen", "greenyellow", "orange", "red")) +
    scale_color_manual(values = c("darkgreen", "greenyellow", "orange", "red")) +
    labs(subtitle = "PCoA\nTurs are united", 
         x = paste0("Axis 1 (", eig[1], " %)"), 
         y = paste0("Axis 2 (", eig[2], " %)"), 
         # fill = "Год", color = "Год", shape = "Год"
    ) +
    theme(panel.grid = element_blank())
if(!export){
    plots$f4_pcoa
}

# PC corelation -----------------------------------------------------------
sites <- readxl::read_excel("data/Carabidae_2025-09-05.xlsx", 
                            "sites", skip = 2) %>% 
    mutate(
        site = as.numeric(str_extract(site, "[:digit:]+")),
        log_km = log10(км)
    ) %>% 
    select(-zone)


way1 <- sites %>% 
    left_join(
        transmute(pc, year, zone, site = as.numeric(site), 
                  plot = as.numeric(plot), Axis.1, Axis.2),
        ., 
        by = c("site", "plot")) %>% 
    select(-zone) %>% 
    unite("id", site, plot, sep = "_")

cor.axis.1 <- way1 %>% 
    select(-Axis.2) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    split(.$param) %>% 
    # `[`(1:2) %>% 
    lapply(function(a){
        b = cor.test(a$Axis.1, a$val)
        tibble(
            Axis.1 = mean(range(a$Axis.1)), 
            val = mean(range(a$val)), 
            r = b$estimate, 
            p = b$p.value) %>% 
            mutate(
                r = round(r, 2), p = round(p, 3), 
                pp = case_when(p <= 0.001 ~ "***", 
                               p<= 0.01 ~ "**", 
                               p <= 0.05 ~ "*",
                               TRUE ~ ""),
                pp = paste0(r, pp))
    }) %>% 
    map_df(rbind, .id = "param")

plots$cor.axis.1 <- way1 %>% 
    select(-Axis.2) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    ggplot(aes(val, Axis.1, color = year)) + 
    geom_point() + 
    geom_smooth(method = "lm", formula = 'y ~ x') + 
    geom_label(data = cor.axis.1, mapping = aes(
        val, Axis.1, label = pp),
        color = "black"
    ) + 
    facet_wrap(~param, scales = "free") + 
    theme_bw() + 
    theme(legend.position = "bottom") + 
    labs(x = NULL, y = NULL, subtitle = "Axis 1, 33%")

cor.axis.2 <- way1 %>% 
    select(-Axis.1) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    split(.$param) %>% 
    # `[`(1:2) %>% 
    lapply(function(a){
        b = cor.test(a$Axis.2, a$val)
        tibble(
            Axis.1 = mean(range(a$Axis.2)), 
            val = mean(range(a$val)), 
            r = b$estimate, 
            p = b$p.value) %>% 
            mutate(
                r = round(r, 2), p = round(p, 3), 
                pp = case_when(p <= 0.001 ~ "***", 
                               p<= 0.01 ~ "**", 
                               p <= 0.05 ~ "*",
                               TRUE ~ ""),
                pp = paste0(r, pp))
    }) %>% 
    map_df(rbind, .id = "param")

plots$cor.axis.2 <- way1 %>% 
    select(-Axis.1) %>% 
    pivot_longer(names_to = "param", values_to = "val", -1:-3) %>% 
    ggplot(aes(val, Axis.2, color = year)) + 
    geom_point() + 
    geom_smooth(method = "lm", formula = 'y ~ x', 
                se = FALSE, ) + 
    geom_label(data = cor.axis.2, mapping = aes(
        val, Axis.1, label = pp),
        color = "black"
    ) + 
    facet_wrap(~param, scales = "free") + 
    theme_bw() + 
    theme(legend.position = "bottom") + 
    labs(x = NULL, y = NULL, subtitle = "Axis 2, 16%")

if(export){
    ggsave("export/corr.axis.1.png", plots$cor.axis.1, 
           height = 9, width = 11, dpi = 300)
    ggsave("export/corr.axis.2.png", plots$cor.axis.2, 
           height = 9, width = 11, dpi = 300)
    lst(axis.1 = cor.axis.1, axis.2 = cor.axis.2) %>% 
        map(~select(.x, Параметр = param, pp)) %>% 
        map_dfr(rbind, .id = "Axis") %>% 
        pivot_wider(names_from = Axis, values_from = pp) %>% 
        writexl::write_xlsx("export/PC, correlation.xlsx")
    
} else { 
    plots$cor.axis.1
    plots$cor.axis.2
}

