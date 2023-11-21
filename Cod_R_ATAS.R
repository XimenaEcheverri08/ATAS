#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Análisis estadístico multivariado aplicado a datos hidrogeológicos usando R: 
# Acuífero Transfronterizo del Amazonas (ATAS) en Leticia - Amazonas
# 
# Autores: Ximena Echeverri & Breiner Bastidas
# lxecheverrig@unal.edu.co, breiner.bastidas@udea.edu.co
# Noviembre 2023
#
# Bases de datos: matriz concentración química aguas subterráneas
# Muestreo de parámetros fisicoquímicos y microbiológicos: 20 pozos 
# Inventario de puntos de agua: datos de consumos, usos, niveles, condiciones sanitarias
# 22 variables HGQ
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Ver variables precargadas, borrar consola = ctl+L
ls()
rm(list=ls())

# Directorio de trabajo (Importante /) 
Guardar <- 'PC';
if (Guardar == 'USB') {ruta = 'C:/Users/Usuario/Documents/REscritorio'};
if (Guardar == 'PC') {ruta = 'C:/Users/XIMENA ECHEVERRI GTZ/Documents/Rday/Información'};
setwd(ruta);
getwd()

# Cargar librerias y funciones
library(openair);            library(ggplot2);             library(dplyr);          
library(reshape);            library(reshape2);            library(scales);   
library(corrplot);           library(cowplot);             library(gridExtra);
library(ggpubr);             library(moments);             library(reshape2);
library(data.table);         library(factoextra);          library(NbClust); 
library(tidyverse);          library(cluster);             library(tidyr)
library(dplyr);              source('f_varianza.R');       source('ggbiplot.R');

# Leer los datos desde el archivo fuente
Datos <- read.table("BD_Inventario_Puntos_Agua_2015.csv", header = TRUE, sep=";", dec=".")
Pozos <- Datos[,c(1)]; Datos <- Datos[,-c(1,2)]                                # Otra forna de escribir: Datosi$Pozo = NULL; Datosi$Codigo= NULL;
Datos_ICA <- Datos
Datos <- data.matrix(Datos)
Datos_ICA <- as.data.frame(Datos_ICA)

# Analisis estadistico
Media <- apply(Datos, 2, mean); 
Mediana <- apply(Datos, 2, median); 
Desv.estandar <- apply(Datos, 2, sd); 
Kurtosis <- apply(Datos, 2, kurtosis); 
Varianza <- apply(Datos, 2, var); 
Estadisticos <- data.frame(Media, Mediana, Desv.estandar, Kurtosis, Varianza);
Datos <- Datos[,-(which(Varianza < 1))]                                         # Especies con varianza igual a cero 

# Matrices análisis estadisticos
HGQ <-Datos[,c(2:11,13:19)];
Names_HGQ <- colnames(HGQ)
colnames(HGQ)[c(3, 8, 9, 10)] <- c("Dureza total", "Hierro total", "Coliformes totales", "Coliformes fecales")

HGQ_t <- t(HGQ)
names_HGQ_t <- rownames(HGQ_t) 

HGQ_Box <- data.frame(Pozos, HGQ[,-c(9,10)]); 
HGQ_Box <- melt(HGQ_Box, id.var="Pozos", na.rm=TRUE); 
colnames(HGQ_Box) <- c("Fecha", "Categoría", "Concentración");

# Determinacion ICA_agua_subterránea
Datos_ICA$OD_IDEAM <- ifelse(Datos_ICA$Sat_Oxigeno < 0, "Error valor negativo", 
                             ifelse(Datos_ICA$Sat_Oxigeno >= 0 & Datos_ICA$Sat_Oxigeno <= 140, ((0.000000031615*Datos_ICA$Sat_Oxigeno^5)-(0.000010304*Datos_ICA$Sat_Oxigeno^4)+(0.0010076*Datos_ICA$Sat_Oxigeno^3)-(0.027883*Datos_ICA$Sat_Oxigeno^2)+(0.84068*Datos_ICA$Sat_Oxigeno)-(0.1612)),
                                    ifelse(Datos_ICA$Sat_Oxigeno > 140, 50)))

Datos_ICA$DBO5_IDEAM <- ifelse(Datos_ICA$DBO5 < 0, "Error valor negativo", 
                               ifelse(Datos_ICA$DBO5 >= 0 & Datos_ICA$DBO5 <= 30, ((0.00018677*Datos_ICA$DBO5^4) - (0.016615*Datos_ICA$DBO5^3) + (0.59636*Datos_ICA$DBO5^2) - (11.152*Datos_ICA$DBO5) + 100.19),
                                      ifelse(Datos_ICA$DBO5 > 30, 2)))

Datos_ICA$Coli.fec_IDEAM <- ifelse(Datos_ICA$Coliformes_fecales < 0, "Error valor negativo", 
                                   ifelse(Datos_ICA$Coliformes_fecales >= 0 & Datos_ICA$Coliformes_fecales <= 100000, (exp((-0.0152*(log(Datos_ICA$Coliformes_fecales))^2)-(0.1063*(log(Datos_ICA$Coliformes_fecales)))+4.5922)),
                                          ifelse(Datos_ICA$Coliformes_fecales > 100000, 2)))

Datos_ICA$pH_IDEAM <- ifelse(Datos_ICA$pH < 0, "Error", 
                             ifelse(Datos_ICA$pH < 2 & Datos_ICA$pH > 12, 0,
                                    ifelse(Datos_ICA$pH >= 2 & Datos_ICA$pH <= 7.5, (-0.1789*(Datos_ICA$pH^5) + 3.7932*(Datos_ICA$pH^4) - 30.517*(Datos_ICA$pH^3) + 119.75*(Datos_ICA$pH^2) - 224.58*Datos_ICA$pH + 159.46),
                                           ifelse(Datos_ICA$pH >= 7.5 & Datos_ICA$pH <= 12, (-1.11429*(Datos_ICA$pH^4) + 44.50952*(Datos_ICA$pH^3) - 656.6*(Datos_ICA$pH^2) + 4215.34762*(Datos_ICA$pH) - 9840.14286)))))

Datos_ICA$T_IDEAM <- 0.0000019619*(Datos_ICA$D_Temperatura^6) - 0.00013964*(Datos_ICA$D_Temperatura^5) + 0.0025908*(Datos_ICA$D_Temperatura^4) + 0.015398*(Datos_ICA$D_Temperatura^3) - 0.67952*(Datos_ICA$D_Temperatura^2) - 0.67204*Datos_ICA$D_Temperatura + 90.392

Datos_ICA$SDT_IDEAM <- ifelse(Datos_ICA$SDT < 0, "Error", 
                              ifelse(Datos_ICA$SDT >= 0 & Datos_ICA$SDT <= 500, (-0.0000000044289*(Datos_ICA$SDT^4)+0.00000465*(Datos_ICA$SDT^3)-0.0019591*(Datos_ICA$SDT^2)+0.18973*(Datos_ICA$SDT) + 80.608),
                                     ifelse(Datos_ICA$pH > 500, 20)))

Datos_ICA$Turbiedad_IDEAM <- ifelse(Datos_ICA$Turbiedad < 0, "Error", 
                                    ifelse(Datos_ICA$Turbiedad >= 0 & Datos_ICA$Turbiedad <= 100, (0.0000018939*(Datos_ICA$Turbiedad^4) - 0.00049942*(Datos_ICA$Turbiedad^3) + 0.049181*(Datos_ICA$Turbiedad^2) - 2.6284*Datos_ICA$Turbiedad + 98.098),
                                           ifelse(Datos_ICA$Turbiedad > 100, 5)))

Datos_ICA$Fosfatos_IDEAM <- ifelse(Datos_ICA$Fosfatos < 0, "Error", 
                                   ifelse(Datos_ICA$Fosfatos >= 0 & Datos_ICA$Fosfatos <= 10, (0.0046732*(Datos_ICA$Fosfatos^6) - 0.16167*(Datos_ICA$Fosfatos^5) + 2.20595*(Datos_ICA$Fosfatos^4) - 15.0504*(Datos_ICA$Fosfatos^3) + 53.8893*(Datos_ICA$Fosfatos^2) - 99.8933*Datos_ICA$Fosfatos + 99.8311),
                                          ifelse(Datos_ICA$Fosfatos > 100, 5)))

Datos_ICA$Nitratos_IDEAM <- ifelse(Datos_ICA$Nitratos < 0, "Error", 
                                   ifelse(Datos_ICA$Nitratos >= 0 & Datos_ICA$Nitratos <= 100, (0.0000000035603*(Datos_ICA$Nitratos^6) - 0.0000012183*(Datos_ICA$Nitratos^5) + 0.00016238*(Datos_ICA$Nitratos^4) - 0.010693*(Datos_ICA$Nitratos^3) + 0.37304*(Datos_ICA$Nitratos^2) - 7.521*Datos_ICA$Nitratos + 100.95),
                                          ifelse(Datos_ICA$Nitratos > 100, 1)))

Datos_ICA$Coliformes_AS <- ifelse(Datos_ICA$Coliformes_fecales < 0, "Error", 
                                  ifelse(Datos_ICA$Coliformes_fecales == 0, 10,
                                         ifelse(Datos_ICA$Coliformes_fecales > 0 & Datos_ICA$Coliformes_fecales <= 500, 0.00004*(Datos_ICA$Coliformes_fecales^2) - 0.04*(Datos_ICA$Coliformes_fecales) + 10,
                                                ifelse(Datos_ICA$Coliformes_fecales > 500 & Datos_ICA$Coliformes_fecales <= 1000, 0,
                                                       ifelse(Datos_ICA$Coliformes_fecales > 1000 & Datos_ICA$Coliformes_fecales <= 5000, -2,
                                                              ifelse(Datos_ICA$Coliformes_fecales > 5000 & Datos_ICA$Coliformes_fecales <=10000, -3,
                                                                     ifelse(Datos_ICA$Coliformes_fecales > 10000, -4, NA)))))))

Datos_ICA$Color_AS <- ifelse(Datos_ICA$Color < 0, "Error",
                             ifelse(Datos_ICA$Color >= 0 & Datos_ICA$Color <= 15, 10,
                                    ifelse(Datos_ICA$Color > 15 & Datos_ICA$Color <= 35, (-0.5*Datos_ICA$Color + 17.5),
                                           ifelse(Datos_ICA$Color > 35, 0, NA))))

Datos_ICA$Dureza_tl_AS <- ifelse(Datos_ICA$Dureza_total < 0, "Error", 
                                 ifelse(Datos_ICA$Dureza_total >= 0 & Datos_ICA$Dureza_total <= 160, 10,
                                        ifelse(Datos_ICA$Dureza_total > 160 & Datos_ICA$Dureza_total <= 200, -0.25*Datos_ICA$Dureza_total+50,
                                               ifelse(Datos_ICA$Dureza_total > 200, 0))))

Datos_ICA$Hierro_tl_AS <- ifelse(Datos_ICA$Hierro_total < 0, "Error", 
                                 ifelse(Datos_ICA$Hierro_total >= 0 & Datos_ICA$Hierro_total <= 0.3, 10,
                                        ifelse(Datos_ICA$Hierro_total > 0.3 & Datos_ICA$Hierro_total <= 1, 17.191*Datos_ICA$Hierro_total^2 - 36.393*Datos_ICA$Hierro_total + 19.266,
                                               ifelse(Datos_ICA$Hierro_total > 1, 0, NA))))

Datos_ICA$Nitratos_AS <- ifelse(Datos_ICA$Nitratos < 0, "Error", 
                                ifelse(Datos_ICA$Nitratos >= 0 & Datos_ICA$Nitratos <= 10, 10,
                                       ifelse(Datos_ICA$Nitratos > 10 & Datos_ICA$Nitratos <= 25, 0.0177*Datos_ICA$Nitratos^2 - 1.2849*Datos_ICA$Nitratos + 21.068,
                                              ifelse(Datos_ICA$Nitratos > 25, 0))))

Datos_ICA$pH_AS <- ifelse(Datos_ICA$pH < 0, "Error", 
                          ifelse(Datos_ICA$pH > 0 & Datos_ICA$pH <= 2.5, 0,
                                 ifelse(Datos_ICA$pH > 2.5 & Datos_ICA$pH <= 6.5, 2.5*Datos_ICA$pH - 6.25,
                                        ifelse(Datos_ICA$pH > 6.5 & Datos_ICA$pH <= 9, 10,
                                               ifelse(Datos_ICA$pH > 9 & Datos_ICA$pH <= 10, -10 * Datos_ICA$pH + 100,
                                                      ifelse(Datos_ICA$pH > 10, 0))))))

Datos_ICA$Turbiedad_AS <- ifelse(Datos_ICA$Turbiedad < 0, "Error", 
                                 ifelse(Datos_ICA$Turbiedad >= 0 & Datos_ICA$Turbiedad <= 5, 10,
                                        ifelse(Datos_ICA$Turbiedad > 5 & Datos_ICA$Turbiedad <= 25, -0.5*Datos_ICA$Turbiedad + 12.50,
                                               ifelse(Datos_ICA$Turbiedad > 25, 0, NA))))

ICA_IDEAM <- 
  0.17*Datos_ICA$OD_IDEAM + 
  0.16*Datos_ICA$Coli.fec_IDEAM + 
  0.11*Datos_ICA$pH_IDEAM +
  0.11*Datos_ICA$DBO5_IDEAM + 
  0.10*Datos_ICA$Nitratos_IDEAM + 
  0.10*Datos_ICA$Fosfatos_IDEAM + 
  0.08*Datos_ICA$Turbiedad_IDEAM + 
  0.07*Datos_ICA$SDT_IDEAM +  
  0.10*Datos_ICA$T_IDEAM
ICA_I1=ICA_IDEAM/10

ICA_IDEAM_corregido = 
  (0.16 + 0.17/8)*Datos_ICA$Coli.fec_IDEAM + 
  (0.11 + 0.17/8)*Datos_ICA$pH_IDEAM + 
  (0.11 + 0.17/8)*Datos_ICA$DBO5_IDEAM  + 
  (0.10 + 0.17/8)*Datos_ICA$Nitratos_IDEAM + 
  (0.10 + 0.17/8)*Datos_ICA$Fosfatos_IDEAM + 
  (0.10 + 0.17/8)*Datos_ICA$T_IDEAM +
  (0.08 + 0.17/8)*Datos_ICA$Turbiedad_IDEAM + 
  (0.07 + 0.17/8)*Datos_ICA$SDT_IDEAM
ICA_IC=ICA_IDEAM_corregido/10

ICA_AS <- 
  0.3*Datos_ICA$Coliformes_AS + 
  0.12*Datos_ICA$Color_AS +
  0.07*Datos_ICA$Dureza_tl_AS +
  0.07*Datos_ICA$Hierro_tl_AS +
  0.25*Datos_ICA$Nitratos_AS+
  0.07*Datos_ICA$pH_AS +
  0.12*Datos_ICA$Turbiedad_AS 

Datos_ICA$ICA_AS <- ICA_AS
Clasificacion_ICA_AS <- ifelse(Datos_ICA$ICA_AS > 10, "Excelente", 
                               ifelse(Datos_ICA$ICA_AS >= 9 & Datos_ICA$ICA_AS < 10 , "Muy buena",
                                      ifelse(Datos_ICA$ICA_AS >= 6 & Datos_ICA$ICA_AS < 9 , "Buena",
                                             ifelse(Datos_ICA$ICA_AS >= 5 & Datos_ICA$ICA_AS < 6 , "Regular",
                                                    ifelse(Datos_ICA$ICA_AS >= 2.5 & Datos_ICA$ICA_AS < 5 , "Mala",
                                                           ifelse(Datos_ICA$ICA_AS >= 0 & Datos_ICA$ICA_AS < 2.5 , "Muy mala",
                                                                  ifelse(Datos_ICA$ICA_AS < 0 & Datos_ICA$ICA_AS > 10 , "Error", NA)))))))


F_ICA_AS <- data.frame(Pozos, ICA_I1, ICA_IC, ICA_AS); F_ICA_AS
F_ICA_AS <- melt(F_ICA_AS, id.var="Pozos", na.rm=TRUE); 
colnames(F_ICA_AS) <- c("Pozos", "Categoria", "Concentracion");

Bar_Plot <- ggplot(F_ICA_AS, aes(Pozos, Concentracion, fill = Categoria)) + 
  geom_bar(stat="identity", position="dodge") +
  labs( x = quote(""), y = quote("ICA") ) +
  scale_fill_manual(labels = c(expression("ICA OTCA 1"), expression("ICA OTCA 2"), expression("ICA MArtinez")), values = c("#327CB7", "#BCDBEA", "#E98E70")) + ##0A3A70 #9FCBE1 #FDEADE
  theme_void() +
  theme(text = element_text(family = "serif", size = 10),
        legend.text = element_text(size = 10), legend.position = "bottom", legend.direction = "horizontal")+
  theme_minimal()+
  theme(text=element_text(family="serif",  size=10), legend.text = element_text(size = 10), legend.position = "bottom", legend.direction = "horizontal", axis.line = element_line(colour = "black", linewidth =0.1))+
  theme(axis.text.x = element_text(colour = 'black', size = 10, angle = 90, hjust = 1)) +
  theme(axis.text.y = element_text(colour = 'black', size = 10)) +
  theme(axis.title.y  = element_text(colour = 'black', size = 10)) #+ ylim (0, 2.5)#+
  #theme(text = element_text(family = "serif", size = 10),
  #legend.text = element_text(size = 10), legend.position = "bottom", legend.direction = "horizontal")
ggsave(Bar_Plot, filename = "Bar_Plot.png", width = 16, height = 8, units = "cm", dpi = 250)

# Análisis de componentes principales 
PCA <- prcomp(HGQ, center = TRUE, scale. = TRUE)                                # Ejecutar PCA
Resumen_PCA <- summary(PCA); Resumen_PCA                                        # Desviacion, y varianza
Eigenvalores_PCA <-PCA$sdev^2; Eigenvalores_PCA                                 # Eigenvalores (Los importantes son los mayores a 1)
Cargas_PCA <- PCA$rotation; Cargas_PCA                                          # Cargas >0.3 y <-0.3

ggbiplot(PCA,alpha = 0.5, labels = PCA$date) + theme_bw()                       # Graficar resultados PCA
Resultado_PCA = summary(PCA) ; Resultado_PCA

# Prueba rotacion VARIMAX
ncomp <- 4 
rawLoadings     <- PCA$rotation[,1:ncomp] %*% diag(PCA$sdev, ncomp, ncomp);rawLoadings
rotatedLoadings <- varimax(rawLoadings)$loadings; rotatedLoadings

ggbiplot(PCA,alpha = 0.5, labels = PCA$date) + theme_bw()                       # Graficar resultados PCA
Resultado_PCA = summary(PCA) ; Resultado_PCA

# Análisis cluster 
data_clusters=scale(HGQ)
row.names(data_clusters) <- Pozos

# Estimar el número de clústers #Elbow, silhouette o gap_stat  method
fviz_nbclust(data_clusters, kmeans, method = "wss")
fviz_nbclust(data_clusters, kmeans, method = "silhouette")
fviz_nbclust(data_clusters, kmeans, method = "gap_stat")

# Optimizacion para definir cuantos cluster debo utilizar (Regla del codo)
n <- 10
tot_error <- matrix(1:n, ncol = 1)
for(i in 1:n) {
  xmeans <- kmeans(data_clusters,centers = i,nstart = 2000)
  tot_error[i] <- xmeans$tot.withinss
}

tot_error<-data.frame(tot_error); tot_error$n_cluster<-matrix(1:n, ncol = 1)
ggplot(tot_error, aes(y=tot_error, x=n_cluster))+geom_point()+theme_bw()+labs(x="Number of clusters", y="Total error")

paleta <- c("violetred", "forestgreen", "dodgerblue4", "orange1", "grey0")       # Ordenar los colores manualmente para que coincida con los del grafico anterior

#Operacion con el valor optimo que se definio graicamente con el paso anterior
Clusters <- kmeans(data_clusters, centers = 5,nstart = 2000)                    # Solo debe cambiar el valor de centers (numero de clusters)
Clusters$cluster <- as.character(Clusters$cluster)

PCA_Plot <- ggbiplot(PCA, alpha = 1, ellipse=TRUE, groups = Clusters$cluster, labels = Pozos) +
  labs(color="Cluster") +
  theme_bw() + 
  theme(legend.position = "top") + 
  scale_color_manual(values=paleta) +
  theme_void() +
  theme(text = element_text(family = "serif", size = 10),
        legend.text = element_text(size = 10), legend.position = "bottom", legend.direction = "horizontal")+
  theme_minimal()+
  theme(text=element_text(family="serif",  size=10), legend.text = element_text(size = 10), legend.position = "bottom", legend.direction = "horizontal", axis.line = element_line(colour = "black", linewidth = 0.1))+
  theme(axis.text.x = element_text(colour = 'black', size = 10, hjust = 1)) +
  theme(axis.text.y = element_text(colour = 'black', size = 10)) +
  theme(axis.title.y  = element_text(colour = 'black', size = 10)) #+ ylim (0, 2.5)#+
#theme(text = element_text(family = "serif", size = 10),
#legend.text = element_text(size = 10), legend.position = "bottom", legend.direction = "horizontal")
ggsave(PCA_Plot, filename = "PCA_Plot.png", dpi = 300)

Cluster_HGQ <- Clusters$cluster
#aggregate(DATOS[,-c(1,2,3,19)], list(DATOS$Cluster), mean)  # Este le muestra la concentracion promedio de cada cluster

######## Cluster No - Jerarquico
d <- dist(data_clusters)
hc <- hclust(d)
fviz_dend(hc, k=5, cex = 0.7, palette = palette, rect = TRUE, rect_fill = TRUE, labels_track_height = 0.8)
#En lo anterior solo debe cambiar el valor de K, dependiendo del numero de agrupaciones que quiera hacer.

# Otra forma de hacer el grafico de dendrograma
m.distancia <- dist(data_clusters, method = "euclidean") # "pearson"
fviz_dist(m.distancia)
hclust_data_cluster <- hclust(m.distancia, method = "ward.D") 
plot(hclust_data_cluster, labels = Pozos)

# Analisis Spearman: Datos no normalizados
HGQ_cor <- cor(HGQ, method = "pearson")
rownames(HGQ_cor)[c(3, 8, 9, 10)] <- c("Dureza total", "Hierro total", "Coliformes totales", "Coliformes fecales")
colnames(HGQ_cor)[c(3, 8, 9, 10)] <- c("Dureza total", "Hierro total", "Coliformes totales", "Coliformes fecales")
grafico_corr <- corrplot(HGQ_cor, method="color", addCoef.col = "black", type = "upper", diag=FALSE, number.cex= 0.5, tl.col = "black", tl.cex = 0.7)

# Diagrama cajas y bigotes
Box_Plot_1 <- ggplot(HGQ_Box, aes(Categoría, Concentración)) + geom_boxplot(aes(fill = Categoría)) + 
  labs( x = quote(""), y = quote("Concentración ("*mu*"g m"^-3*")")) +  
  #scale_fill_manual(values=c("red", "blue", "green", "green","green"))+ 
  theme_bw() +
  theme(text = element_text(family="serif",  size=11)) +
  #theme(axis.text.x = element_text(colour = 'black', size = 10, angle = 90, hjust = 1)) +
  theme_minimal()+
  theme(axis.text.x = element_text(family="serif", colour = 'black', size = 10, angle = 90, hjust = 1)) +
  theme(axis.text.y = element_text(family="serif", colour = 'black', size = 10)) +
  theme(axis.title.y  = element_text(family="serif", colour = 'black', size = 10)) + #+ ylim (0, 2.5)+
  theme(text = element_text(family = "serif", size = 11))
Box_Plot_1
ggsave(Box_Plot_1, filename = "Box_Plot.png", width = 16, height = 8, units = "cm", dpi = 700)
