##################################################################
######### 2° Workshop: Experimentação na agricultura 4.0 #########
##################################################################

# Instruções:
# Baixar o repositório '2dn_workshop_a4p0'
# Instalar a versão de desenvolvimento do pliman com o comando abaixo
# Descomente para rodar

#devtools::install_github("TiagoOlivoto/pliman")

############################### Pacotes e diretório #####################
library(pliman)
library(tidyverse)
set_wd_here()
################### Contagem e mensuração de objetos ####################
# https://tiagoolivoto.github.io/pliman/reference/analyze_objects.html
# https://tiagoolivoto.github.io/pliman/articles/indexes.html
leaves <- image_import("flax.jpg", plot = TRUE)
image_index(leaves)

leaves_meas <-
  analyze_objects(leaves,
                  show_lw = TRUE,
                  index = "B")

# plotar as medidas
plot_measures(leaves_meas,
              measure = "width",
              col = "green",
              hjust = -90)
plot_measures(leaves_meas, 
              measure = "length", 
              vjust = 60,
              col = "red")


#################### Correção de medidas ############
seg <- 
  image_segment_iter(leaves,
                     index = c("R/(G/B)", "B-R"))

leaves_meas_cor <-
  analyze_objects(leaves,
                  reference = TRUE,
                  reference_area = 20)
image_view(leaves, object = leaves_meas_cor)


################### Análise em lote ##################
set_wd_here("potato_leaves")
potato <-
  analyze_objects(pattern = "G",
                  watershed = FALSE,
                  reference = TRUE,
                  reference_area = 20,
                  marker = "area",
                  save_image = TRUE,
                  dir_processed = "proc",
                  plot = FALSE,
                  efourier = TRUE,
                  object_index = "DGCI")
###### Descritores Fourier ######
img <- image_import("G05.jpg")
ocont <- 
  img |> 
  object_contour(index = "B",
                 watershed = FALSE,
                 filter = 2)
ef <- efourier(ocont) |> efourier_inv()
plot(img)
plot_contour(ef, col = "red", lwd = 2)

###### extraindo insigths #######
meas <- get_measures(potato)
stats <- 
  meas$results |> 
  group_by(img) |> 
  summarise(across(c(area, solidity), .fns = list(m = mean, se = \(x){sd(x)/sqrt(length(x))})))

ggplot(meas$results) +
  geom_point(aes(area, solidity, color = img)) +
  geom_point(data=stats, size = 4,
             aes(x = area_m, y = solidity_m,
                 color = img)) +
  geom_errorbar(data = stats,
                aes(x = area_m,
                    ymin = solidity_m - solidity_se,
                    ymax = solidity_m + solidity_se,
                    color = img)) +
  geom_errorbarh(data = stats,
                 aes(y = solidity_m,
                     xmin = area_m - area_se,
                     xmax = area_m + area_se,
                     color = img)) +
  theme_minimal(base_size = 16)

library(factoextra)
library(FactoMineR)

dfpca <- 
  meas$results |> 
  select(-c(x, y, id, contains("radius_"))) |>
  select(-c(A4:A10,
            B4:B10,
            C4:C10,
            D4:D10)) |> 
  group_by(img) |> 
  summarise(across(where(is.numeric), mean)) |> 
  column_to_rownames("img")

pcam <- PCA(dfpca, graph = FALSE)
fviz_pca_ind(pcam, repel = TRUE)



################### Severidade de doenças ##################
## Interativa
set_wd_here()
img <- image_import("sev_orange.jpg", plot = TRUE)
sev <- measure_disease_iter(img, viewer = "mapview")
sev$results

## Utilizando índices de cores
set_wd_here("sevsoja")
sev_lote <- 
  measure_disease(pattern = "soy",
                  index_lb = "B",
                  index_dh = "NGRDI",
                  threshold = c("Otsu", -0.03),
                  plot =  FALSE,
                  save_image = TRUE,
                  dir_processed = "proc",
                  show_contour = FALSE,
                  parallel = TRUE,
                  col_lesions = "brown")

sev_lote$severity |> 
  ggplot(aes(x = symptomatic)) +
  geom_histogram(bins = 8)



################### Análise de ortomosaicos #################

set_wd_here("orthomosaics")

### Índices de vegetação
# Mosaicos disponíveis no pipeline do FieldimageR (https://github.com/OpenDroneMap/FIELDimageR)
# Grato ao Filipe Matias (https://github.com/filipematias23) por disponibilizar esses materiais
pot <- mosaic_input("potato.tif")
res_pot <- 
  mosaic_analyze(pot,
                 r = 1, g = 2, b = 3,
                 nrow = 9,
                 ncol = 16,
                 buffer_row = -0.15,
                 buffer_col = -0.05,
                 plot_index = c("BGI", "VARI", "NGRDI", "GLI"),
                 summarize_fun = c("mean", "stdev"),
                 plot = FALSE)

mosaic_view(pot,
            r = 1, g = 2, b = 3,
            shapefile = res_pot$result_plot,
            attribute = "mean.NGRDI")

# shapefile usado para a análise
shp <- res_pot$shapefile[[1]] |> shapefile_input()
shapefile_plot(shp)

## modelo digital de superfície
dsm0 <- mosaic_input("dsm_0.tif") 
dsm70 <- mosaic_input("dsm_70.tif")
dsm0 <- mosaic_resample(dsm0, dsm70)
dsm <- dsm70 - dsm0
names(dsm) <- "ph"
ph <- mosaic_analyze(dsm,
                     segment_plot = TRUE,
                     segment_index = "ph",
                     shapefile = shp,
                     plot = FALSE)
mosaic_view(dsm,
            shapefile = ph$result_plot,
            attribute = "coverage")


# Stand de plantas e CV na linha de semeadura
mosaic <- mosaic_input("standcount.tif")
an <- mosaic_analyze(mosaic,
                     r = 1, g = 2, b = 3,
                     ncol = 7,
                     nrow = 1,
                     segment_individuals = TRUE,
                     segment_index = "GLI",
                     check_shapefile = FALSE,
                     map_individuals = TRUE,
                     attribute = "n")
an$result_individ_map

# overlap the two plots
mosaic_view(mosaic,
            r = 1, g = 2, b = 3,
            shapefile = an$result_plot_summ,
            color_regions = custom_palette(c("red", "green", "yellow"), n = 10),
            attribute = c("cv"),
            alpha = 0.6) +
  mosaic_view(mosaic,
              r = 1, g = 2, b = 3,
              shapefile = an$result_indiv,
              attribute = c("area"))


## Índices multiespectrais
# Mosaico disponibilizado pelo Programa Trigo UFV
# Dr. Maicon Nardino (https://www.linkedin.com/in/maicon-nardino-091253268/)
# Caique Machado (https://www.linkedin.com/in/caique-machado-e-silva-404415154/)
trigo <- mosaic_input("trigo_ex.tif")
mosaic_view(trigo)
trigoan <- 
  mosaic_analyze(trigo,
                 nrow = c(5, 1),
                 ncol = c(21, 3),
                 buffer_row = c( -0.15, 0),
                 buffer_col = c(0, -0.05),
                 plot_index = c("NDVI", "EVI"),
                 attribute = "mean.NDVI")
trigoan$map_plot






################ Monitor de colheita  ####################
library(pliman)
set_wd_here("orthomosaics")
moni <- mosaic_input("monitor.jpg")
mosaic_view(moni, r = 1, g = 2, b = 3)
# one block with 10 plot
# 8 plots additional
monires <- 
  mosaic_analyze(moni,
                 r = 1, g = 2, b = 3,
                 grid = c(TRUE, rep(FALSE, 8)),
                 buffer_col = c(-0.05, rep(0, 8)),
                 buffer_row = c(-0.05, rep(0, 8)),
                 nrow = 5,
                 ncol = 2,
                 plot_index = c("R", "G", "B"),
                 summarize_fun = NULL,
                 plot = FALSE)


library(tidyverse)
# Train data
df_train <- 
  monires$result_plot |> 
  as.data.frame() |> 
  select(-geometry) |> 
  filter(block %in% paste0("B0", 2:9)) |> 
  mutate(class = c("130", "125", "115", "105", "95", "85", "75", "70")) |> 
  unnest(cols = data)

# Training a Random Forest model
library(caret)
control <-
  trainControl(method = 'cv',
               p = 0.7,
               number = 5,
               verboseIter = TRUE)
fit <- train(class ~ R + G + B,
             data = df_train,
             method = 'rf',
             trControl = control,
             ntree = 300)

# Data to predict (plots)
df_predict <- 
  monires$result_plot |> 
  as.data.frame() |> 
  select(-geometry) |> 
  filter(block %in% paste0("B01")) |> 
  unnest(cols = data)


# Predict the yield at pixel level and compute the weigthed average
dfpredicted <- 
  df_predict |> 
  mutate(class = as.character(predict(fit, newdata = df_predict))) |> 
  mutate(class = as.numeric(class)) |> 
  group_by(plot_id, class) |> 
  summarise(n = n(), 
            .groups = "drop") |> 
  group_by(plot_id) |> 
  mutate(rgclass = class * n) |> 
  summarise(rgmedio = sum(rgclass) / sum(n)) |> 
  bind_cols(monires$shapefile[[1]]) |> 
  sf::st_as_sf(crs = terra::crs(moni))

# plot the results
mosaic_view(moni,
            r = 1, g = 2, b = 3,
            shapefile = dfpredicted,
            attribute = "rgmedio")






