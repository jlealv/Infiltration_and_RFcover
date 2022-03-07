######################### packages
library(spdep)
library(ape)
library(sp)
library(MVA)
library(normtest)
library(nortest)
library(readxl)
library(readxl)
library(psych) 
library(spatialreg)
#########################use your own directory
cobsupinf_art1_DATANEW <- read_excel("C:/Users/57316/Desktop/cobsupinf_art1_DATA.xlsx", 
                                     sheet = "XX2")
#View(cobsupinf_art1_DATANEW)
# 
datat=data.frame(cobsupinf_art1_DATA)  #data.frame
infi=datat$INFILKOST                      # response variable  
X=as.matrix(datat[,3:15])                 #Data matrix
XYdata=as.matrix(datat[,1:2])             #coordinates
########################################   creation of the two weight matrices
inf.d=as.matrix(dist(XYdata, diag=T, upper=T))   # distances
inf.d.inv <-as.matrix(1/inf.d)                   #inverse distances
inf.d.inv2 <-as.matrix(1/inf.d**2)               #quadratic inverse distances
diag(inf.d.inv) <- 0                             # zero-diagonal (inverse distances) 
diag(inf.d.inv2) <- 0                            # zero-diagonal (quadratic inverse distances)
W=as.matrix(inf.d.inv)                #weight matrix- inverse distance  W         
W2=as.matrix(inf.d.inv2)              #weight matrix- quadratic inverse distance W2  
sumas=apply(W,1,sum)                      #sums for Standarization W
sumas2=apply(W2,1,sum)                    #sums for Standarization W2
We=W/sumas                               # standarization-W
We2=W2/sumas2                            #standarization-W2
#################### visualization of the weight matrices
plot(We[lower.tri(We)],col="darkred",pch=6,cex=0.5,xlab="Elements in the weight matrix",
     ylab="row-standarized weights",main="Distribution of non-zero weights",
     family="Times New Roman")
points(We2[lower.tri(We2)],col="darkblue",pch=25,cex=0.5)
points(We[upper.tri(We)],col="red",pch=2,cex=0.5)
points(We2[upper.tri(We2)],col="blue",pch=24,cex=0.5)
abline(h=0)
legend(5,0.4, legend=c("id-l", "iqd-l","id-u","iqd-u"),
       col=c("darkred", "darkblue","red","blue"), pch=c(6,25,2,24), cex=0.8)
##############################  exploring spatial dependence in response and  COBFRSUP
Moraninf=(Moran.I(infi, We));Moraninf  #
plot(XYdata[,1],XYdata[,2],col=round(infi,0),pch=19)
grid(10,10)
MoranP=(Moran.I(datat$COBFRSUP, We));MoranP
plot(XYdata[,1],XYdata[,2],col=round(datat$COBFRSUP,0),pch=19)
grid(10,10)
############################# descriptive statistics
options(digits=5)
describe(X)
boxplot(datat$COBFRSUP,main="COBSUP")
boxplot(datat$INFILKOST,main="INFILT")
hist(datat$COBFRSUP,main="COBSUP")
hist(datat$INFILKOST,main="INFILT")
################################ exploring Normality in variables
shapi=c()
for(j in 3:15){
  shapi[j-2]=sf.test(datat[,j])$p.value
}
shapi
###################### Spatial models
##########creation of weights matrix as required by the library Spatialreg
contnb=dnearneigh(coordinates(XYdata),0,380000,longlat = F)
dlist <- nbdists(contnb, XYdata)
dlist <- lapply(dlist, function(x) 1/x)            #inverse distance
Wve=nb2listw(contnb,glist=dlist,style = "W")       #W matriz-standarized
####################pure autoregresive model ########################## without Log
map=spautolm(INFILKOST~1,data=datat,listw=Wve) ############## firts model
summary(map)
resi1=map$fit$residuals                    #extraction of residuals
plot(XYdata[,1],XYdata[,2],col=round(abs(resi1),0)+1,pch=19) #residuals display
grid(10,10)
hist(resi1)             # #histogram of residuals
sf.test(resi1)$p.value  #Normality of residuals
####################pure autoregresive model ########################## with Log
mapln=spautolm(log(INFILKOST)~1,data=datat,listw=Wve) ############## 
summary(mapln)
resi1ln=mapln$fit$residuals                      #extraction of residuals
plot(XYdata[,1],XYdata[,2],col=round(abs(resi1ln),0)+1,pch=19)#residuals display
grid(10,10)
hist(resi1ln)                #histogram of residuals
sf.test(resi1ln)$p.value     #Normality of residuals
inf_e=as.data.frame(mapln$fit["fitted.values"]);inf_e #predicted residuals
DFe=data.frame(log(datat$INFILKOST),inf_e) # Log-response and predicted Log-response         
colnames(DFe) <- c("ln_inf","inf_e")  
plot(DFe$ln_inf ,DFe$inf_e,cex=0.5,pch=19)########scatter - Log-response and predicted Log-response  
cor(DFe$ln_inf ,DFe$inf_e)  #Pearson correlation- Log-response and predicted Log-response  
Moranresi1=(Moran.I(resi1ln,We))  # Spatial dependence in residuals

################################################spatial error model  with Log
mserln =errorsarlm(formula=log(INFILKOST)~COBFRSUP,data=datat,listw=Wve)
summary(mserln)
resi2ln=mserln$residuals       #extraction of residuals
plot(XYdata[,1],XYdata[,2],col=round(abs(resi2ln),0)+1,pch=19)#residuals display
grid(10,10)
hist(resi2ln)                        #histogram of residuals
sf.test(resi2ln)$p.value             #Normality of residuals
inf_e2=as.data.frame(mserln$fitted.values);inf_e2  #predicted response new model
DFe2=data.frame(log(datat$INFILKOST),inf_e2)  
colnames(DFe2) <- c("ln_inf","inf_e2") 
plot(DFe2$ln_inf ,DFe2$inf_e2,cex=0.5,pch=19)#scatter - Log-response and predicted Log-response
cor(DFe2$ln_inf ,DFe2$inf_e2)#Pearson correlation- Log-response and predicted Log-response
Moranresi2=(Moran.I(resi2ln,We))# Spatial dependence in residuals
############################################Spatial Lag model 
meer=lagsarlm(formula=log(INFILKOST)~COBFRSUP,data=datat,listw=Wve)
summary(meer)
resi3ln=meer$residuals         #extraction of residuals
plot(XYdata[,1],XYdata[,2],col=round(abs(resi3ln),0)+1,pch=19) #residuals display
grid(10,10)
hist(resi3)               #histogram of residuals
sf.test(resi3ln)$p.value         #Normality of residuals
inf_e3=as.data.frame(meer$fitted.values);inf_e3  #predicted response new model
DFe3=data.frame(log(datat$INFILKOST),inf_e3)
colnames(DFe3) <- c("ln_inf","inf_e3") 
plot(DFe3$ln_inf ,DFe3$inf_e3,cex=0.5,pch=19)#####scatter - Log-response and predicted Log-response
cor(DFe3$ln_inf ,DFe3$inf_e3)  #Pearson correlation- Log-response and predicted Log-response
Moranresi3=(Moran.I(resi3ln,We))# Spatial dependence in residuals
##############################################SARAR model
contnb=dnearneigh(coordinates(XYdata),0,380000,longlat = F)
dlist <- nbdists(contnb, XYdata)
dlist <- lapply(dlist, function(x) 1/x**2)   #Second weight matrix 
Wve=nb2listw(contnb,glist=dlist,style = "W")
SARAR=sacsarlm(formula=log(INFILKOST)~COBFRSUP+datat$PORTOT,data=datat,listw=Wve) #Sarar model
summary(SARAR)
RESSARAR=SARAR$residuals          #extraction of residuals
plot(XYdata[,1],XYdata[,2],col=round(abs(RESSARAR),0)+1,pch=19)#residuals display
grid(10,10)
hist(RESSARAR)               #histogram of residuals
sf.test(RESSARAR)$p.value #Normality of residuals
shapiro.test(RESSARAR)#Normality of residuals
inf_eSARAR=as.data.frame(SARAR$fitted.values);inf_eSARAR   #predicted response final model
DFeSARAR=data.frame(log(datat$INFILKOST),inf_eSARAR) 
colnames(DFeSARAR) <- c("ln_inf","inf_e_SARAR") 
#####scatter - Log-response and predicted Log-response
plot(DFeSARAR$ln_inf ,DFeSARAR$inf_e_SARAR,cex=0.5,pch=19,
     xlab="Predicted values (log(response))",ylab="Observed values(log(response)) ")
grid(10,10)
cor(DFeSARAR$ln_inf ,DFeSARAR$inf_e_SARAR)#Pearson correlation- Log-response and predicted Log-response
text(0.8,6, "r=0.863")
Moranresi3=(Moran.I(RESSARAR,We))# Spatial dependence in residuals
DFINAL=data.frame(datat,inf_eSARAR)  # original data and predicted response in final model
###################### exporting predictions
library(openxlsx)
write.xlsx(DFINAL, file = "C:/Users/57316/Desktop/avila/pred_finales.xlsx", colNames = TRUE, borders = "surrounding")






