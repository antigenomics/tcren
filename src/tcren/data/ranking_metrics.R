# Existing tied-score ROC AUC and average precision; shared unchanged by publication summaries.
metric <- function(y,x) {
 n1<-sum(y==1);n0<-sum(y==0)
 if(!n1 || !n0) return(c(roc_auc=NA_real_,average_precision=NA_real_))
 auc<-(sum(rank(x,ties.method='average')[y==1])-n1*(n1+1)/2)/(n1*n0)
 z<-data.table(x=x,y=y)[,.(pos=sum(y),n=.N),by=x][order(-x)]
 tp<-cumsum(z$pos);rec<-tp/n1;prec<-tp/cumsum(z$n)
 c(roc_auc=auc,average_precision=sum(diff(c(0,rec))*prec))
}
