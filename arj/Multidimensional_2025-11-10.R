export <- FALSE
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
sp <- long %>% 
    pull(taxa) %>% 
    unique %>% 
    sort %>% 
    str_subset("Gen | sp", negate = TRUE)

wide <- wide %>% 
    select(year:plot, all_of(sp)) %>% 
    mutate(zone = case_when(zone %in% c("bufer", "background") ~ "control", TRUE ~ "impact")) %>% 
    arrange(zone)

ref <- wide %>% 
    my_pcoa %>% 
    pluck("dis") %>% 
    dis.to.df() %>% 
    group_by(zone1, zone2) %>% 
    summarise(mean_dis = mean(dis), .groups = "drop") %>% 
    filter(zone1 != zone2) %>%
    unite("zone", zone1, zone2, sep = "-") %>% 
    pivot_wider(names_from = zone, values_from = mean_dis)



dif <- sp %>% 
    # `[`(1:5) %>%
    as.list() %>% 
    `names<-`(sp[1:length(.)]) %>% 
    mclapply(function(x){
        wide %>% 
            my_pcoa(sp.to.remove = x) %>% 
            pluck("dis") %>% 
            dis.to.df() %>% 
            group_by(zone1, zone2) %>% 
            summarise(mean_dis = mean(dis), .groups = "drop") %>% 
            filter(zone1 != zone2) %>%
            unite("zone", zone1, zone2, sep = "-") %>% 
            pivot_wider(names_from = zone, values_from = mean_dis)
        }, 
        mc.cores = detectCores()-1
    ) %>% 
    map_dfr(rbind, .id = "taxa")
    
    # map_dfr(
    #     ~wide %>% 
    #         my_pcoa(sp.to.remove = .x) %>% 
    #         pluck("dis") %>% 
    #         dis.to.df() %>% 
    #         group_by(zone1, zone2) %>% 
    #         summarise(mean_dis = mean(dis), .groups = "drop") %>% 
    #         filter(zone1 != zone2) %>%
    #         unite("zone", zone1, zone2, sep = "-") %>% 
    #         pivot_wider(names_from = zone, values_from = mean_dis), 
    #     .id = "taxa"
    # )

for(i in 1:ncol(ref)){
    dif[[i+1]] <- dif[[i+1]] - ref[[i]]
}

compared <- "impact-control"

dif %>% 
    select_at(c("taxa", compared)) %>% 
    arrange_at(compared) %>% 
    slice(1:8, (nrow(.)-7):nrow(.)) %>% 
    mutate(
        taxa = factor(taxa),
        taxa = fct_inorder(taxa)
    ) %>% 
    ggplot(aes(x = taxa, ymin = `impact-control`, ymax = 0)) + 
    geom_errorbar() + 
    geom_hline(yintercept = 0, color = "red") + 
    coord_flip() + 
    labs(
        x = NULL,
        y = paste0("Baseline distance = ", round(ref$`impact-control`, 3))) + 
    theme(axis.text.y = element_text(face = "italic"))
    

