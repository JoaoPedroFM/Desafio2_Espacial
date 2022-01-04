## Estatistica Espacial I
## 21/12/2021

library(raster)
library(tidyverse)
library(geoR)
library(readxl)

# Lendo a base de dados
base <- read_excel("seguro_auto.xlsx")
dim(base)
names(base)


# smp_size <- 500
# set.seed(123)
# train_ind <- sample(seq_len(nrow(base)), size = smp_size)
# train <- base[train_ind, ]
# # write.csv2(train,"rio_premio.csv")
# 
# base <- train


# Lendo o SHAPE
library(rgdal)

mapaRJ <- readOGR("Limite_Bairro_si.shp",layer="Limite_Bairro_si")

# Mapa com pontos
par(mar=c(4,4,0,0),cex=1.4)
plot(mapaRJ,xlab="longitude",ylab="latitude")
points(base$long,base$lat,col=2,pch=21,bg=5)
axis(1)
axis(2)

# Histograma do premio
par(mfrow=c(1,2),mar=c(4,4,0.5,0.5))
hist(base$premio,main="",xlab="premio")
hist(log(base$premio),main="",xlab="log(premio)")

base$log_premio <- log(base$premio)

# Transformando em geodata (para usar as funcoes do pacote geoR)
dados <- as.geodata(base, coords.col = 3:2, data.col = 4, covar.col=1)

plot(dados)

par(mfrow=c(1,1))
points(dados,cex.min=1, cex.max=5, pt.sizes="quintiles",
       col=terrain.colors(5))

# Variograma empirico
v <- geoR::variog(dados,trend= ~ dados$coords)

par(mfrow=c(1,1),mar=c(4,4,0.5,0.5),cex=1.5)
plot(v,pch=16,col=4,xlab="distancia",ylab="semivariograma",ylim=c(0,1))

# eyefit(v)

# Estimando via MV
# opcoes para trend: cte (media constante), 1st (tendencia linear) e 2st (tendencia quadratica).

# trend.spatial(~coords + sinistro, dados)[1:5,]

# Estimando o efeito pepita:
mv <- likfit(dados, ini=c(0.1,0.1), trend = trend.spatial(~coords, dados), cov.model="gaussian",
             fix.nugget = TRUE,nugget=0.3)
summary(mv)

plot(v,col=4,xlab="distancia",ylab="semivariograma",ylim=c(0,0.6))
lines(mv,lty=1,col="blue")



### Interpolacao ponderada pelo inverso da distancia

library(gstat)
library(sp)

base.idw <- as.data.frame(cbind(dados$coords,dados$data))
coordinates(base.idw) <- c("long", "lat")

summary(dados$coords)
long = seq(-43.8,-43.1,l=50)
lat = seq(-23.1,-22.79,l=50)
nx = length(long)
ny = length(lat)
# criando o grid de pontos para interpolar
grid = as.data.frame(expand.grid(long,lat))

# coordenadas dos pontos do grid
coordinates(grid) = ~Var1+Var2


# interpolacao
interpol <- idw(V3~1, base.idw, newdata = grid, idp = 2)
interpol.matriz = matrix(interpol$var1.pred,nx,ny)

# grafico
#opcao 1
filled.contour(long, lat,interpol.matriz,
               key.title = title(main = "log premio", cex.main = 1),
               plot.axes={axis(1)
                 axis(2)
                 plot(mapaRJ,xlab="longitude",ylab="latitude",add=TRUE)
               })




### Krigagem

long.pred = seq(-43.8,-43.1,l=50)
lat.pred = seq(-23.1,-22.79,l=50)
grid.pred <- expand.grid(long.pred,lat.pred)


# krigagem
krigagem <- krige.conv(dados, loc=grid.pred,
                       krige=krige.control(type.krige="SK",
                                           trend.d = trend.spatial(~coords, dados),
                                           trend.l = trend.spatial(~Var1+Var2, grid.pred),
                                           beta=mv$beta,
                                           nugget=mv$nugget,
                                           cov.pars=mv$cov.pars))

hist(krigagem$predict)
# media
par(mfrow=c(1,1),mar=c(4,4,0.5,0.5))
image(krigagem, loc=grid.pred,
      col=terrain.colors(256))
plot(mapaRJ,xlab="longitude",ylab="latitude",add=TRUE)

# variancia
par(mfrow=c(1,1),mar=c(4,4,0.5,0.5))
image(krigagem, loc=grid.pred, coords=dados$coords,
      values=krigagem$krige.var, col=terrain.colors(256))
plot(mapaRJ,xlab="longitude",ylab="latitude",add=TRUE)
