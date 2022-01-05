## Estatistica Espacial I
## 21/12/2021

library(raster)
library(tidyverse)
library(geoR)


# Lendo a base de dados
base <- read_csv2("Dados_Rio_new.csv")
dim(base)
names(base)


smp_size <- 500
set.seed(123)
train_ind <- sample(seq_len(nrow(base)), size = smp_size)
train <- base[train_ind, ]
# write.csv2(train,"rio_premio.csv")

base <- train


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
# comportamento assimetrico, variabilidade maior, valores de seguro muito altos (outliers)
hist(log(base$premio),main="",xlab="log(premio)") 
# comportamento se assemelhando a distribuicao normal -> se escolhe modelar o log do premio
# diminui a variabiidade, estibilizando um pouco a incerteza

base$log_premio <- log(base$premio)

# Transformando em geodata (para usar as funcoes do pacote geoR)
dados <- as.geodata(base, coords.col = 8:7, data.col = 17, covar.col=c(1:6,9:16))

plot(dados)
# dado aparentando distribuicao normal, mais pro norte e mais pra direita (dependencia espacial)

par(mfrow=c(1,1))
points(dados,cex.min=1, cex.max=5, pt.sizes="quintiles",
       col=terrain.colors(5))
# quantis: mais pro verde e menor a bola -> menor o valor
# quantis: mais pro branco e maior a bola -> maior o valor
# Indicacoes
# segurados na regiao da Barra sao parecidos entre si, com valores altos de seguro
# em cima -> talvez ha valores altos
# zona a noroeste -> segurados com valores proximos

# Variograma empirico
v <- geoR::variog(dados,trend= ~ dados$coords) # utiliza coordenadas na tendencia

par(mfrow=c(1,1),mar=c(4,4,0.5,0.5),cex=1.5)
plot(v,pch=16,col=4,xlab="distancia",ylab="semivariograma",ylim=c(0,1))
# alcance pequeno, altura = efeito pepita, comportamento nao suave

# eyefit(v)

# Estimando via MV
# opcoes para trend: cte (media constante), 1st (tendencia linear) e 2st (tendencia quadratica).

# trend.spatial(~coords + sinistro, dados)[1:5,]

# Estimando o efeito pepita:
mv <- likfit(dados, ini=c(0.1,0.1), trend = trend.spatial(~coords, dados), 
             cov.model="gaussian",fix.nugget = TRUE,nugget=0.3) # nugget= efeito pepita
# testar modelo da familia matern, esferica e exponencial
summary(mv)
# medidas melhores no modelo de estrutura de dependencia espacial, se ajusta melhor
# do que o modelo que ignora a dependencia espacial na variavel do premio
# (modelo linear que supoe independencia)
# Melhor usar o modelo espacial 

plot(v,col=4,xlab="distancia",ylab="semivariograma",ylim=c(0,0.6))
lines(mv,lty=1,col="blue")
# alcance pequeno, logo chega no patamar, sem longas dependencias no espaco
# nao ha tanta dependencia espacial em pequenas distancias

# se houvesse alcance grande -> processo suave, ate longas distancias apresentam 
# dependencias no espaco


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

# Indicacao
# mapa do log do premio, mais escuro -> maior o log do preco -> maior o preco
# regiao que apresentam precos muito baixos, mas na mesma, ha locais com precos muito altos

# de acordo com o alcance pequeno, na regiao escura, em pequenas distancias nao ha 
# dependencias expressivas entre os dados

# Zona Sul precos mais BAIXOS, regioes com vias, proximos de locais favelizados e 
# divisa de municipios - mais ALTOS



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





########################
### Modelo Bayesiano ###
########################
# FREQUENTISTA

# Gerando os dados
dados <- grf(100, cov.pars=c(10, .4),cov.model="spherical",nug=1)
mv <- likfit(dados, ini=c(10,0.4), nug=1,cov.model="spherical")


plot(dados$coords,pch=16)

variograma=variog(dados)
plot(variograma,col="blue")
lines(mv,lty=1,col="blue")


# Definindo a malha de pontos para predicao
grid <- expand.grid(seq(0, 1, l=30), seq(0, 1, l=30))

plot(dados$coords,pch=16)
points(grid,pch=16,col=2)

krigagem <- krige.conv(dados, loc=grid, krige=krige.control(type.krige="SK", trend.d="cte",beta =mv$beta, cov.model="spherical",cov.pars=mv$cov.pars,
                                                            nugget=mv$nugget))


image(krigagem, loc=grid, coords=dados$coords, col=terrain.colors(256),
      x.leg=c(0.8,0.85), y.leg=c(0.6,0.9),cex.leg=0.7,vertical=TRUE)

image(krigagem, loc=grid, coords=dados$coords, col=terrain.colors(256),
      values=krigagem$krige.var, x.leg=c(0.8,0.85), y.leg=c(0.6,0.9),cex.leg=0.7,vertical=TRUE)


# BAYESIANO

# Fazendo inferencia bayesiana
bayes <- krige.bayes(dados, loc=grid, model = model.control(cov.model = "spherical"),
                     prior = prior.control(beta.prior = "flat", sigmasq.prior = "reciprocal",
                                           phi.prior="uniform", phi.discrete=seq(0,1,l=10),
                                           tausq.rel.prior = "uniform",
                                           tausq.rel.discrete = seq(0, 1, l=10)))

# Priori e posteriori para o parametro phi e tau^2/sigma2
par(mfrow=c(1,2))
plot(bayes, phi.dist = TRUE, tausq.rel.dist = TRUE, type="bars",col=c("red","blue"))
par(mfrow=c(1,1))

# Histogramas da Posteriori
par(mfrow=c(2,2))
hist(bayes,density.est = T)
par(mfrow=c(1,1))

# Plotando o variograma empirico e algumas estimativas bayesianas:

# Variograma empirico
plot(variog(dados))

# Variograma empirico com o verdadeiro valor do variograma no grafico
plot(dados)

# Adicionando a media, mediana e moda a posteriori ao variograma
lines(bayes, summ = mean, lwd=1, lty=2,col=1)
lines(bayes, summ = median, lwd=2, lty=2,col=2)
lines(bayes, summary = "mode", post = "parameters",col=3)

# Plotando o variograma empirico e adicionando a mediana e quantis estimados
plot(variog(dados),ylim=c(0,20))
quantis <- function(x){quantile(x, prob = c(0.05,0.5, 0.95))}
lines(bayes, summ = quantis, ty="l", lty=c(2,1,2), col=c(4,2,4))

# Predicao
image(bayes, main="Valores preditos",col=terrain.colors(256),x.leg=c(0.85,0.9), y.leg=c(0.8,1), cex.leg=0.9,vertical=TRUE)
# points(dados, cex.min=1, cex.max=4, pt.sizes="deciles",col=terrain.colors(10),add=T)

# Variancia de Predicao
image(bayes, val="variance", main="Variancia de predicao",col=terrain.colors(256),x.leg=c(0.85,0.9), y.leg=c(0.8,1), cex.leg=0.9,vertical=TRUE)
