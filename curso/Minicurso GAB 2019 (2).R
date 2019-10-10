
# Carga de librerías
install.packages (c("randomForest", "caret"))

library(gstat) #
library(geoR) #
library(mapview) #
library(raster) #
library(leaflet) #
library(RColorBrewer) #
library(PerformanceAnalytics) #
library(ggplot2) #
library(caret) #
library(parallel) #

# Seleccion de directorio de trabajo
setwd(choose.dir(getwd(), "Seleccione Directorio de Trabajo"))

# Carga de base de datos
datos <- read.table("datosSF_depurados.txt", header = T)

# Graficos Exploratorios
chart.Correlation(datos[,1:9], histogram=TRUE, pch=19)
chart.Correlation(datos[,c(1,9:18)], histogram=TRUE, pch=19)

plot(datos$x, datos$VUT)
plot(datos$y, datos$VUT)

hist(datos$VUT)
skewness(datos$VUT)   

#Transformación a objeto espacial
coordinates(datos) <- c("x", "y")
crs(datos) <- CRS("+init=epsg:22174")

# Visualización
spplot(datos["VUT"],col.regions=terrain.colors(20))

cols <-colorRampPalette(rev(brewer.pal(11, "Spectral"))) # Genero una paleta
mapview(datos, zcol="VUT",alpha=0, col.regions =cols) 


################################################################
####################### Ajuste Semivariogramas #################
################################################################

# Ajuste de semivariograma experimetal
semi_exp <- variogram(VUT~1, datos)
plot(semi_exp)

# Ajuste de semivariograma teórico, modelo exponencial
mod_exp <- fit.variogram(semi_exp, vgm(3000000, "Exp", 1500, 10000))
mod_exp
plot(semi_exp , mod_exp)

# Ajuste de semivariograma teórico, modelo esférico
mod_esf = fit.variogram(semi_exp, vgm(3000000, "Sph", 1500, 10000))
mod_esf
plot(semi_exp , mod_esf)

# Visualización conjunta 
vgLine <- rbind(
  cbind(variogramLine(mod_exp, maxdist = max(semi_exp$dist)), id ="Exponencial"),
  cbind(variogramLine(mod_esf, maxdist = max(semi_exp$dist)), id ="Esférico")
  )

grafico <-  ggplot(semi_exp, aes(x = dist, y = gamma, color = id)) +
  geom_line(data = vgLine) + 
  geom_point() +
  labs(title ="Semivariograma experimental y teórico ajustado")  +
  xlab("Distancia") + 
  ylab("Semivarianza")  + 
  scale_color_discrete(name = "Modelo", labels = c("Esférico", "Exponencial", "Experimental"))

grafico


# Suma de cuadrado del error de modelos ajustados
attr(mod_exp, 'SSErr')
attr(mod_esf, 'SSErr')

# Ajuste automatico de modelos
v.fit_todos <- fit.variogram(semi_exp,vgm(c("Exp", "Sph", "Gau")))
plot(semi_exp,v.fit_todos)

attr(v.fit_todos, 'SSErr')

################################################################
####################### Interpolación Kriging  #################
################################################################

# Generación de grilla de predicción
Limites <- read.table("bordesSF.txt",header = T)
gr <- pred_grid(Limites, by=100) ;  plot(gr)
gri <- polygrid(gr,bor=Limites)
plot(gri)
names(gri)[1]<-paste("x")
names(gri)[2]<-paste("y")
gridded(gri) = ~x+y
plot(gri)
crs(gri) <- CRS("+init=epsg:22174")

# Carga grilla predicción ejes de manzana
ejes <- read.table("ejesSF.txt", dec = ",",header = T)
names(ejes)
names(ejes)[names(ejes) =="g_perc_edif"] <- "g_porc_edi"
names(ejes)[names(ejes) =="g_perc_baldio"] <- "g_por_bald"
names(datos)

coordinates(ejes) = ~x+y
crs(ejes) <- CRS("+init=epsg:22174")
plot(ejes)
mapview(ejes)

# Kriging Ordinario
kriging_grilla <- krige(VUT~1, datos, gri, model = v.fit_todos)
spplot(kriging_grilla["var1.pred"], main = "Kriging Ordinario: Predicciones", col.regions=terrain.colors(20))
spplot(kriging_grilla["var1.var"],  main = "Kriging Ordinario: Varianza", col.regions=terrain.colors(20))

kriging_ejes <- krige(VUT~1, datos, ejes, model = v.fit_todos)
spplot(kriging_ejes["var1.pred"], main = "Kriging Ordinario: Predicciones", col.regions=terrain.colors(20))
spplot(kriging_ejes["var1.var"],  main = "Kriging Ordinario: Varianza", col.regions=terrain.colors(20))


# Extracción de predichos y varianza de predicción de grilla
pred <- kriging_grilla$var1.pred
varpred <-kriging_grilla$var1.var
cord_X <- kriging_grilla$x
cord_Y <- kriging_grilla$y
tablapred <- cbind(cord_X,cord_Y,pred, varpred)

# Visualización interactiva
cols <-colorRampPalette(rev(brewer.pal(11, "Spectral"))) # Genero una paleta

muestra <- mapview(datos,zcol = "VUT", ceX ="VUT", layer.name ="Muestra", col.regions = cols, alpha=0) 

pred_grilla <- rasterFromXYZ(tablapred[,1:3])
crs(pred_grilla) <-CRS("+init=epsg:22174")
mapapred_grilla <- mapview(pred_grilla,legend=T, col.regions =  cols, layer.name ="Predichos Grilla",na.color = "transparent") 

var_grilla <- rasterFromXYZ(tablapred[,c(1,2,4)])
crs(var_grilla) <-CRS("+init=epsg:22174")
mapavar_grilla <- mapview(var_grilla,legend=T, col.regions =  cols, layer.name ="Varianza Predicción Grilla",na.color = "transparent") 

crs(kriging_ejes) <-CRS("+init=epsg:22174")
mapapred_ejes <- mapview(kriging_ejes, zcol = "var1.pred", ceX ="var1.pred",col.regions =  cols, layer.name ="Predicho Ejes", alpha=0) 
mapavar_ejes <- mapview(kriging_ejes, zcol = "var1.var", ceX ="var1.var",col.regions =  cols, layer.name ="Varianza Predicción Ejes", alpha=0)

muestra +  mapapred_grilla + mapavar_grilla + mapapred_ejes + mapavar_ejes

##### Opciones visualización
muestra_predichos <- muestra +  mapapred_ejes

muestra_predichos %>% 
  leafem::addLogo ("https://idecor.cba.gov.ar/wp-content/uploads/2017/07/cropped-logo-IDECOR.png",position = "topleft" , width = 180, height = 60) %>%
  leafem::addHomeButton(extent(pred_grilla), "Ciudad de San Francisco" , position = "bottomright") %>%
  addMiniMap(position = "bottomleft" , width = 150, height = 150) 


muestra_predichos %>%   
  leafem::addLogo("https://tenor.com/view/gallardo-super-campe%c3%b3n-gif-11367643.gif",position = "topleft",offset.x = 5, offset.y = 100, width = 200, height = 200) %>%
  # leafem::addHomeButton(ext = extent(r), "Ciudad de San Francisco" , position = "bottomright") %>%
  addMiniMap(position = "bottomleft" , width = 100, height = 100)


################################################################
####################### Validación Cruzada Kriging #############
################################################################
mod_exp
mod_esf

attr(mod_exp, 'SSErr')
attr(mod_esf, 'SSErr')

set.seed(17)
kricv_mod_esf <- krige.cv(VUT~1, datos, mod_esf, nfold=10)
set.seed(17)
kricv_mod_exp <- krige.cv(VUT~1, datos, mod_exp, nfold=10)

bubble(kricv_mod_esf, "residual", main = "Residuos Esferico") 
bubble(kricv_mod_exp, "residual", main = "Residuos Exponencial")

# Error medio de predicción (ME), cercano a cero mejor:
mean(kricv_mod_esf$residual)
mean(kricv_mod_exp$residual)

# Error medio absoluto (MAE)
mean(abs(kricv_mod_esf$residual))
mean(abs(kricv_mod_exp$residual))

# Error cuadratico medio de predicción (MSE), mas pequeño mejor
mean(kricv_mod_esf$residual^2)
mean(kricv_mod_exp$residual^2)

# Mean squared deviation ratio (MSDR), Error cuadratico medio normalizado, cercano a 1 mejor
mean(kricv_mod_esf$zscore^2)
mean(kricv_mod_exp$zscore^2)

# RMSE relativo a la media
sqrt(mean(kricv_mod_esf$residual^2))/mean(kricv_mod_esf$observed)*100
sqrt(mean(kricv_mod_exp$residual^2))/mean(kricv_mod_exp$observed)*100

# Correlación lineal entre valores observados y predichos
cor(kricv_mod_esf$observed, kricv_mod_esf$observed - kricv_mod_esf$residual)
cor(kricv_mod_exp$observed, kricv_mod_exp$observed - kricv_mod_exp$residual)

# Correlación lineal entre valores observados y predichos
par(mfrow=c(1,2))
plot(kricv_mod_esf$observed,kricv_mod_esf$observed - kricv_mod_esf$residual, xlab="Observados", ylab="Predichos")
plot(kricv_mod_exp$observed,kricv_mod_exp$observed - kricv_mod_exp$residual,xlab="Observados", ylab="Predichos")



################################################################
####################### Kriging Regresión ######################
################################################################

# Ajuste de modelo de RLM
mlr <- lm(VUT~g_porc_edi+ d_supcom + d_centro + d_indust + d_front , datos)

# Incorporamos los residuos del MLR a la base de datos
datos$residuos <-mlr$residuals
names(datos)

# Ajuste de semivariograma experimetal y teórico a los reiuduos
semiv_rk <- variogram(residuos~1 , datos)
plot(semiv_rk)

semiv_rk <- variogram(residuos~1 , datos, cutoff=2300)
plot(semiv_rk)

v.fit_vut_rk <- fit.variogram(semiv_rk ,vgm(c("Exp","Sph","Gau")))
plot(semiv_rk ,v.fit_vut_rk)

# Kriging sobre residuos
kgres <- krige(residuos~1, datos, ejes, model = v.fit_vut_rk)
spplot(kgres["var1.pred"], main = "Kriging Residual: Predicciones", col.regions=terrain.colors(20))

# Predicción final
ejes$RK_pred <- predict(mlr, newdata=ejes) + kgres$var1.pred

spplot(ejes["RK_pred"], main = "Predicción RK", col.regions=terrain.colors(20))

mapapred_ejesRK <- mapview(ejes, zcol = "RK_pred", ceX ="RK_pred",col.regions =  cols, layer.name ="Predicho Ejes RK", alpha=0)

muestra +  mapapred_grilla + mapapred_ejes + mapapred_ejesRK


################################################################
####################### Random Forest Kriging ##################
################################################################
datos <- read.table("datosSF_depurados.txt", header = T)

# Ajuste del RF con librería caret
seed <-7

# grilla de valores para mtry a evaluar
mtry <-expand.grid(mtry=seq(2,5,1)) 

# opciones de validación
fitControl <- trainControl(method = "cv",number=10, allowParallel = T)

# opciones para praralelizado 
# library(parallel)
# library(doParallel)
# cluster <- makeCluster(detectCores() - 1) 
# registerDoParallel(cluster)

# ajuste del modelo de RF
set.seed(seed)
train_rf <- train(VUT ~ g_porc_edi + d_supcom + d_centro + d_indust + d_front, 
                  data=datos,
                  method = "rf",
                  importance=T,
                  tuneGrid =mtry,
                  trControl = fitControl)

# Incorporamos los residuos del MLR a la base de datos
datos$residuosRF <-datos$VUT - predict(train_rf, newdata=datos)

# Ajuste de semivariograma experimetal y teórico a los residuos del RF
coordinates(datos) <- c("x", "y")
crs(datos) <- CRS("+init=epsg:22174")

semiv_RFk <- variogram(residuosRF~1 , datos)
plot(semiv_RFk)

semiv_RFk <- variogram(residuosRF~1 , datos, cutoff=2300)
plot(semiv_RFk)

v.fit_vut_RFk <- fit.variogram(semiv_RFk ,vgm(c("Exp","Sph","Gau")))
plot(semiv_RFk ,v.fit_vut_RFk)

# Kriging sobre residuos del RF
kgresRF <- krige(residuosRF~1, datos, ejes, model = v.fit_vut_RFk)
spplot(kgresRF["var1.pred"], main = "Kriging Residual (RF): Predicciones", col.regions=terrain.colors(20))

# Predicción final RF
ejes$RFK_pred <- predict(train_rf, newdata=ejes) + kgresRF$var1.pred

spplot(ejes["RFK_pred"], main = "Predicción RFK", col.regions=terrain.colors(20))

mapapredejes_RFK <- mapview(ejes, zcol = "RFK_pred", ceX ="RFK_pred",col.regions =  cols, layer.name ="Predicho Ejes RFK", alpha=0)

muestra +  mapapred_grilla + mapapred_ejes +  mapapred_ejesRK + mapapredejes_RFK


################################################################
####################### Validación Cruzada         #############
################################################################

validacion <-function (fold, base, var.y) {
  require(caret)
  require(gstat)
  require(sp)

  datos <- read.table(base, head=T)
  names(datos)[names(datos) == var.y] <- 'Y'
 
  if (base=="petrel.txt") {
    names(datos)[names(datos) == 'long'] <- 'x'
    names(datos)[names(datos) == 'lat'] <- 'y'
    }
  
  seed <-7
  
  set.seed(seed) 
  datos$id <- sample(rep(1:10, nrow(datos), length.out = nrow(datos)))
  
  list <- 1:10
  prediccion <- data.frame()
  testset<- data.frame()
  
  training<- subset(datos, id %in% list[-fold]) 
  testing<- subset(datos, id %in% c(fold))
  
  # Kriging Ordinario 
  train_ko = training
  test_ko = testing
  coordinates(train_ko)<-~x+y
  coordinates(test_ko)<-~x+y
  vario <- variogram(Y ~1, train_ko)
  VF_vut_KO <- fit.variogram(vario, vgm(c("Sph", "Exp", "Gau"))) 
  KO <- krige(Y~ 1, train_ko, test_ko, VF_vut_KO)
  
  # Regression Kriging  
  train_ko = training
  test_ko = testing
  
  coordinates(train_ko)<-~x+y
  coordinates(test_ko)<-~x+y
  
  mlr <- lm(Y ~. -x -y -id, training)
  
  pred_mlr = predict(mlr, newdata = test_ko)
  
  inside_rk <- predict(mlr, newdata=train_ko)
  train_ko$error_rk <- training$Y - inside_rk
  
  vario_rk <- variogram(error_rk~1, train_ko, cutoff=2300)
  model_rk_ko <- fit.variogram(vario_rk, vgm(c("Sph", "Exp", "Gau"))) 

  test_k <- krige(error_rk~ 1 , train_ko, test_ko, model_rk_ko)
  test_rk_ko <- pred_mlr + test_k$var1.pred
  
  # Random Forest 
  #fitControl <- trainControl(method = "cv", number = 10)
  fitControl <- trainControl(method = "none")
  #mtry <-data.frame(mtry=2) 
  set.seed(seed)
  
  rf <- train(Y ~ . -x -y -id, data=training,
              method = "rf",
              #tuneGrid =mtry,
              trControl = fitControl,
              verbose = FALSE)

  test_rf <- predict(rf, newdata=testing)
  
  # Random Forest + Kriging Ordinario 
  inside_rf <- predict(rf, newdata=training)
  
  train_ko = training
  test_ko = testing
  
  coordinates(train_ko)<-~x+y
  coordinates(test_ko)<-~x+y
  
  train_ko$error_rf <- training$Y - inside_rf
  
  vario_rf <- variogram(error_rf~1, train_ko, cutoff=2300)
  
  model_rf_ko <- fit.variogram(vario_rf, vgm(c("Sph","Exp", "Gau"))) 
  
  test_ko <- krige(error_rf~ 1 , train_ko, test_ko, model_rf_ko) 
  
  test_rf_ko <- test_rf + test_ko$var1.pred
  

  # Tabla observados y predichos
  testset <- rbind(testset, as.data.frame(testing[,"Y"]))
  result <- data.frame(data.frame("x"=testing$x,
                                  "y"=testing$y,
                                  "k-fold"=fold,
                                  "Observado"=testset[,1],
                                  "KO"=KO$var1.pred,
                                  "RK"=test_rk_ko,
                                  "RF"=test_rf,
                                  "RF_KO"=test_rf_ko))
  
  return(result)

}

# correr validacion cruzada
#resultados <- do.call(rbind,lapply(1:10, validacion))

# correr validacion cruzada paralelizado
num_cores <- detectCores()-1
cl <- makeCluster(num_cores)
system.time(resultados <-do.call(rbind,parLapply(cl, 1:10, validacion, base="meuse.txt", var.y="om")))


# Comparación de métodos
head(resultados)
tabla <- resultados[,4:8]
resumen <- function (j) {
  ME <-mean(tabla [,j] - tabla[,"Observado"])
  MAE <- mean(abs(tabla [,j] - tabla[,"Observado"]))
  MAPE <- mean(abs(tabla [,j]-tabla[,"Observado"])/tabla[,"Observado"]) *100
  MSE <- mean((tabla [,j]-tabla[,"Observado"])^2)
  RMSE <-sqrt(mean((tabla [,j]-tabla[,"Observado"])^2))
  nRMSE <-sqrt(MSE)/mean(tabla[,"Observado"]) *100
  rLM <- lm(tabla [,j]~ tabla[,"Observado"])
  R2 <- as.matrix(summary(rLM)$adj.r.squared)
  mx <- mean(tabla[,"Observado"])
  my <- mean(tabla [,j])
  s2x <- var(tabla[,"Observado"])
  s2y <- var(tabla [,j])
  sxy <- mean((tabla[,"Observado"]-mx) * (tabla [,j]-my))
  CCC <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  resumen <- data.frame("Modelo"=names(tabla [j]),ME, MAE,MAPE, MSE, RMSE, nRMSE,R2,CCC)
  return(resumen) 
} 

tablafinal <- do.call("rbind",lapply(2:5,resumen))  
tablafinal

# Otras bases
datosAP <- read.table("datosAP.txt", head=T)
names(datosAP)

datos_meuse <- read.table("meuse.txt", head=T)
names(datos_meuse)

datos_petrel <- read.table("petrel.txt", head=T)
names(datos_petrel) # mud es asimetrica!!!







