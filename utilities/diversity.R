## Shannon's equitability
local.rnaseq.shannon <- function(exp.mat, pseudoNum = 0){
  # calculate shannon index from transcriptome matrix
  apply(exp.mat, 2, function(x){
    x<-x+pseudoNum
    prop<-x/sum(x)
    #prop<-prop[prop>0]
    shidx = -sum(prop*log(prop), na.rm=T)/log(length(prop))
    shidx
  })
}