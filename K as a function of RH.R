require(deSolve)

require (scatterplot3d)


#Parameters, probabilities and rates **All rates are annual
################################################################################
#B=uninfected triatomine bugs
#V= Trypanosoma cruzi-infected triatomine bugs (vectors)
#Ia= humans in the acute phase of Chagas disease
#Ii=humans in the chronic indeterminate phase of Chagas disease
#Id=humans in the chronic determinate phase of Chagas disease
#IR= animals, some of which are competent T. cruzi hosts
#Rfv = relative fitness of infected triatomines
#K= triatomine carrying capacity
#r = triatomine birth rate
#mu = triatomine death rate  
#muP= increase in death rate due to vector control 
#mu.R=animal death rate
#muN= background human mortality rate
#Beta.iab = rate of triatomine-host contact*probability of transmission from Ia-B
#Beta.iib = rate of triatomine-host contact*probability of transmission from Ii-B
#Beta.idb =  rate of triatomine-host contact*probability of transmission from Id-B
#Beta.bi= rate of triatomine-host contact*probability of transmission from B-N
#Beta.irb = rate of triatomine-animal host contact*probability of transmission from IR-B
#Beta.bir = rate of triatomine-animal host contact*probability of transmission from B-IR
#RH= number of animals in the system
#N=number of humans
#delta= rate of movement from Ia-Ii
#sigma= rate of movement from Ii-Id  
#alpha =Id death rate
#ViabRH= Proportion of IR that can transmit CHagas- accounts for birds and maintains uninfecteds in population


################################# #############################################
##
## MODEL WITH HUMANS & ANIMALS
##
###############################################################################

ChagasResModel=function(t, x, params){
  
  B=x[1]
  V=x[2]
  Ia=x[3]
  Ii=x[4]
  Id=x[5]
  IR=x[6]
  
  with(as.list(params),{
    
    dB =  r*(B+(RfV*V))*((K - (B+V))/K) - (mu*muP)*B - (B*((Beta.iab*Ia)+(Beta.iib*Ii)+
                                                             (Beta.idb*Id)+(Beta.irb*IR)))/(RH+N)
    
    dV =  B*((Beta.iab*Ia) + (Beta.iib*Ii)+ (Beta.idb*Id)+(Beta.irb*IR))/(RH+N) - ((mu*muP)*V)
    
    dIa =  (Beta.bi*V)*(N-(Ia+Ii+Id))/(N+RH) - Ia*(delta+muN)
    
    dIi =  delta*Ia -  Ii*(sigma+muN)
    
    dId = sigma*Ii - Id*(alpha+muN)
    
    dIR = V*Beta.bir*((ViaRH*RH)-IR)/(N+RH) - mu.R*IR 
    
    res=c(dB, dV, dIa, dIi, dId, dIR)
    list(res)
  }) }

#####################################################################################
##    Bifurcation across changes in abundance of animals with K as a function of RH
##
#ResHosts<-20     ## Density of other hosts that bugs feed on
ViabRH<-0.01      ##  Proportion of these hosts that can transmit CHagas- accounts for birds and maintains uninfecteds in population

RfV = 0.7        ## relative fitness of infected bugs
Host = 10

endt = 50.0
int = 12
times = seq(0, endt, by = 1/int)

reslst = endt*int
resfst = int*(endt-5)

sections = 40

Rmax = 200
Rmin = 5

Rsec=seq(Rmin, Rmax, length=sections+1)

BmatR=matrix(NA, ncol = endt*int+1, nrow = sections+1)
VmatR=matrix(NA, ncol = endt*int+1, nrow = sections+1)
IimatR=matrix(NA, ncol = endt*int+1, nrow = sections+1)
IamatR=matrix(NA, ncol = endt*int+1, nrow = sections+1)
IdmatR=matrix(NA, ncol = endt*int+1, nrow = sections+1)
IRmatR=matrix(NA, ncol = endt*int+1, nrow = sections+1)


for(i in 1:sections){
  
  paras  = c(r = 36, mu=1.73, mu.R = 0.5,muN=0.013, Beta.iab = 25.01, Beta.iib = 1.066, 
             Beta.idb = 8.2,Beta.irb = 19.078, Beta.bi=  0.02378, Beta.bir = 0.8169,
             muP=1, ViaRH = ViabRH, RH=Rsec[i], N=Host, K=100,  delta = 6, sigma = 0.03, alpha =0.1, RfV=RfV) 
  
  #K=((Rsec[i]*25)+100)
  
  xstart = c( B=100, V=10, Ia=0, Ii=0, Id=0, IR=0)
  
  ############################################
  
  outB = as.data.frame(lsoda(xstart, times, ChagasResModel, paras))
  for(j in 1:(endt*int)){
    BmatR[i,j]=outB$B[j]
    VmatR[i,j]=outB$V[j]
    IamatR[i,j]=outB$Ia[j]
    IimatR[i,j]=outB$Ii[j]
    IdmatR[i,j]=outB$Id[j]
    IRmatR[i,j]=outB$IR[j]
  } 
  
  ##
  ##
  
  cat(i, "\n")
}

############ PLOT IT

sel=seq(resfst, reslst, by=int)
plot(1E-4, xlim=c(Rmin,Rmax), ylim=c(1, 700), log='y',
     col='white', xlab='Animal abundance', ylab='Number of individuals',bty='l')
for(i in 1:sections){
  #points(rep(Rsec[i], length(sel)), IamatR[i, sel], col='blue')
  points(rep(Rsec[i], length(sel)), IimatR[i, sel], col='purple')
  points(rep(Rsec[i], length(sel)), IdmatR[i, sel], col='green') 
  points(rep(Rsec[i], length(sel)), BmatR[i,sel], col='orange', pch=20) 
  points(rep(Rsec[i], length(sel)), VmatR[i, sel], col='red', pch=20)
}

par(xpd=TRUE)
leg.text<-c( "Indeterminate", "Determinate", "Bugs", "Vectors")
legend("topright",inset=c(-0.08,-0.14),leg.text,lwd=3,col=c( "purple", "green", "orange","red" ),bty='n')  

leg.text<-c( "Indeterminate", "Determinate", "Bugs", "Vectors")
legend("topright",inset=c(-0.05,0.15),leg.text,lwd=3,col=c( "purple", "green", "orange","red" ),bty='n') 

leg.text<-c( "Indeterminate", "Determinate")
legend("topright",inset=c(-0.05,0.15),leg.text,lwd=3,col=c( "purple", "green"),bty='n') 

par(xpd=TRUE)
leg.text<-c( "Indeterminate", "Determinate")
legend("topright",inset=c(-0.05,-0.1),leg.text,lwd=3,col=c( "purple", "green" ),bty='n')    


