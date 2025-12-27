library(tidyverse)
library(sf)
library(terra)
library(parallel)
s <- st_read("2025l5.shp") %>% 
    st_zm %>% 
    st_cast("LINESTRING") %>% 
    mutate(l = st_length(geometry)) %>% 
    arrange(l)

s_3857 <- st_transform(s, 3857)
s_cstm <- st_transform(s, crs = 32600 + floor((st_bbox(s)[1] + 180) / 6) + 1)

res <- list()
# r <- s_3857
r <- s_cstm
st0 <- Sys.time()
    
for(i in 1:nrow(r)){
    st <- Sys.time()
    res[i] <- st_line_sample(
        r$geometry[[i]], 
        density = NULL, 
        sample = seq(0, 1, by = 30 / as.numeric(r$l[[i]]))) 
    cat(i, " ### ", round(Sys.time()-st, 3), "sec\n")
}
cat("\n  ### ### ", round(Sys.time()-st0, 3), "min\n")

# pts_3857 <- res_3857 %>% 
#     map_dfr(~.x %>% as.matrix %>% as.data.frame) %>% 
#     st_as_sf(coords = c("V1", "V2"), crs = st_crs(s_3857))
res_cstm <- res

pts_cstm <- res_cstm %>% 
    map_dfr(~.x %>% as.matrix %>% as.data.frame) %>% 
    st_as_sf(coords = c("V1", "V2"), crs = st_crs(s_cstm))

st_write(pts_cstm, "points.gpkg", "pts_cstm")
st_write(pts_cstm, "pts_cstm.shp")




res <- res[100:105] %>% 
    map_dfr(~st_as_sfc(.x)) %>% 
    `st_crs<-`(st_crs(r))
    # st_cast("POINT")


ggplot() + 
    geom_sf(data = r) + 
    theme_bw()

