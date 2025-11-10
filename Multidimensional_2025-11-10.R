export <- FALSE
# Data load ---------------------------------------------------------------
library(tidyverse)
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
    readxl::read_excel(sheet = 1) %>%  
    group_by(year, tur, zone, site, plot, taxa) %>% 
    summarise(num = sum(num), abu = sum(abu), .groups = "drop") %>% 
    arrange(desc(site), taxa) %>% 
    mutate(
        tur = factor(tur), 
        site = fct_inorder(site),
        zone = factor(zone, 
            levels = c("background", "bufer", "impact", "industrial_barren", ordered = TRUE)),
        km = as.numeric(str_extract(site, "[:digit:]{1,}")),
        year = as.factor(year), 
        .after = "site")

# 1 = turs 1 and 2 are united
wide <- long %>% 
    select(-num, -tur) %>% 
    pivot_wider(names_from = taxa, 
                values_from = abu, 
                # quantitative vartiables come from abundance (inds. per 100 trap-days)
                values_fn = sum, values_fill = 0)

# res <- list()
# plots <- list()
# tables <- list()

# Multidimensional  --------------------------------------------------
my_pcoa <- function(wide, sp.to.remove = "", binary = FALSE){
    if(nchar(sp.to.remove) > 0){
        wide[sp.to.remove] <- NULL
    }
    
    dis <- wide %>% 
        # select(-site) %>% 
        unite("ID", zone, km, plot, year, site, sep = "__") %>% # year, 
        column_to_rownames("ID") %>% 
        vegan::vegdist(method = "bray", binary = binary)
    pc <- ape::pcoa(dis)
    eig <- pc$values$Eigenvalues
    eig <- round(eig/sum(eig)*100, 1)
    pc <- pc$vectors %>% 
        as.data.frame() %>% 
        rownames_to_column("ID") %>% 
        as_tibble() %>% #slice(37, 38, 39, 40, 41, 42) %>% select(1)
        separate(ID, into = c("zone", "km", "plot", "year", "site"), 
                 sep = "__") %>% 
        mutate(zone = factor(zone, levels = levels(long$zone)))
    lst(dis, pc, eig)
}

dis.to.df <- function(d, cols = "zone"){
    distances <- as.matrix(d)
    distances[upper.tri(distances)] <- NA
    diag(distances) <- NA
    # distances <- 
    distances %>%
        as.data.frame() %>%
        rownames_to_column("id1") %>%
        as_tibble() %>%
        pivot_longer(names_to = "id2", values_to = "dis", -id1) %>%
        filter(!is.na(dis)) %>%
        separate(id1, paste0(cols, "1"), sep = "_", extra = "drop") %>%
        separate(id2, paste0(cols, "2"), sep = "_", extra = "drop") 
}

# dis.to.df(dis, c("zone", "km", "plot", "year"))

ref <- wide %>% 
    my_pcoa %>% 
    pluck("dis") %>% 
    dis.to.df() %>% 
    group_by(zone1, zone2) %>% 
    summarise(mean_dis = mean(dis), .groups = "drop") %>% 
    filter(zone1 != zone2) %>%
    unite("zone", zone1, zone2, sep = "-") %>% 
    pivot_wider(names_from = zone, values_from = mean_dis)

sp <- sort(unique(long$taxa))

dif <- sp %>% 
    `[`(1:5) %>% 
    as.list() %>% 
    `names<-`(sp[1:length(.)]) %>% 
    map_dfr(
        ~wide %>% 
            my_pcoa(sp.to.remove = .x) %>% 
            pluck("dis") %>% 
            dis.to.df() %>% 
            group_by(zone1, zone2) %>% 
            summarise(mean_dis = mean(dis), .groups = "drop") %>% 
            filter(zone1 != zone2) %>%
            unite("zone", zone1, zone2, sep = "-") %>% 
            pivot_wider(names_from = zone, values_from = mean_dis), 
        .id = "taxa"
    )

for(i in 1:ncol(ref)){
    
}

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

set.seed(1); res$permanova <- vegan::adonis2(
    dis ~ zone * year, 
    data = wide, 
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


# Export ------------------------------------------------------------
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