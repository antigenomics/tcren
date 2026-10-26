#!/usr/bin/env Rscript
# Summarize CLI-derived paired energies; resample whole epitope clusters within each MHC class.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1L)
a<-commandArgs(TRUE)
stopifnot(length(a)==3L);out<-a[1]
e<-fread(file.path(out,'energies.tsv'),colClasses=c(energy='character'));e[,energy:=as.numeric(energy)]; dr<-fread(file.path(out,'matrix_drift.tsv'))
stopifnot(!anyDuplicated(e,by=c('pdb.id','mhc_class','cutoff','potential','replicate')))
B<-2000L;set.seed(20260725L);ans<-list();availability<-list()
for(cl in c('MHCI','MHCII')){
 x<-e[mhc_class==cl];ids<-sort(unique(x$pdb.id))
 co<-unique(x[,.(pdb.id,epitope)])[match(ids,pdb.id)]
 groups<-match(co$epitope,unique(co$epitope));g<-max(groups);n<-length(ids)
 G<-diag(g)[groups,,drop=FALSE]; W<-rmultinom(B,g,rep(1/g,g)); sizes<-colSums(G)
 stopifnot(g>1L,all(x[,uniqueN(pdb.id),by=.(cutoff,potential,replicate)]$V1==n))
 for(cut in sort(unique(x$cutoff))){
  blocks<-lapply(c('Fixed 5 Å','Re-derived','Published'),function(pot){
   d<-dcast(x[cutoff==cut & potential==pot],pdb.id~replicate,value.var='energy')[match(ids,pdb.id)]
   E<-as.matrix(d[,-'pdb.id']);stopifnot(ncol(E)==21L,all(is.finite(E)))
   mu<-rowMeans(E[,-1]);sd<-apply(E[,-1],1,stats::sd)
   list(potential=pot,E=E,z=(E[,1]-mu)/sd,valid=sd>1e-9)
  })
  valid<-Reduce('&',lapply(blocks,`[[`,'valid'))
  availability[[length(availability)+1L]]<-data.table(mhc_class=cl,cutoff=cut,pdb.id=ids,z_defined=valid)
  for(bl in blocks){
   E<-bl$E;U<-matrix(0,n,n)
   # AUC averages all twenty decoys; keep each native and its entire decoy set paired.
   for(rep in 2:21){diff<-outer(-E[,1],-E[,rep],'-');U<-U+(diff>0)+.5*(diff==0)}
   U<-U/20;K<-crossprod(G,U%*%G)
   auc<-sum(U)/n^2;ab<-colSums(W*(K%*%W))/as.vector(crossprod(sizes,W))^2
   z<-ifelse(valid,-bl$z,0);zg<-as.vector(crossprod(G,z));ng<-as.vector(crossprod(G,as.numeric(valid)))
   zb<-as.vector(crossprod(zg,W))/as.vector(crossprod(ng,W))
   ci<-function(v)as.numeric(quantile(v[is.finite(v)],c(.025,.975),names=FALSE,type=7))
   aci<-ci(ab);zci<-ci(zb)
   ans[[length(ans)+1L]]<-data.table(mhc_class=cl,cutoff=cut,potential=bl$potential,n_complexes=n,n_epitopes=g,
    n_z=sum(valid),n_z_epitopes=uniqueN(groups[valid]),n_decoys=20L,n_bootstrap=B,
    auc=auc,auc_lo=aci[1],auc_hi=aci[2],negative_mean_Z=mean(-bl$z[valid]),z_lo=zci[1],z_hi=zci[2])
  }
 }
}
s<-rbindlist(ans);fwrite(s,file.path(out,'summary.tsv'),sep='\t');fwrite(rbindlist(availability),file.path(out,'availability.tsv'),sep='\t')
p<-rbind(s[,.(mhc_class,cutoff,potential,metric='Native/decoy ROC AUC',value=auc,lo=auc_lo,hi=auc_hi)],
 s[,.(mhc_class,cutoff,potential,metric='Negative mean standardized energy',value=negative_mean_Z,lo=z_lo,hi=z_hi)],
 dr[,.(mhc_class,cutoff,potential='Matrix stability',metric='Matrix correlation with 5 Å',value=matrix_R,lo=NA_real_,hi=NA_real_)],
 dr[,.(mhc_class,cutoff,potential='Matrix stability',metric='Mean absolute entry difference',value=matrix_MAD,lo=NA_real_,hi=NA_real_)])
fwrite(p,file.path(out,'cutoff_points.tsv'),sep='\t')
mac<-c('% Formatted from the corrected class-separated calibration.')
for(cl in c('MHCI','MHCII')){
 d<-s[mhc_class==cl];tag<-if(cl=='MHCI')'I' else 'II'
 vals<-list(N=unique(d$n_complexes),Epitopes=unique(d$n_epitopes),ZMin=min(d$n_z),ZMax=max(d$n_z),
  FixedAucFive=d[cutoff==5 & potential=='Fixed 5 Å',auc],PublishedAucFive=d[cutoff==5 & potential=='Published',auc])
 for(nm in names(vals))mac<-c(mac,paste0('\\newcommand{\\Cutoff',tag,nm,'}{',if(grepl('Auc',nm))sprintf('%.3f',vals[[nm]])else vals[[nm]],'}'))
}
mac<-c(mac,'\\newcommand{\\CutoffBootstrapN}{2000}','\\newcommand{\\CutoffDecoyN}{20}')
writeLines(mac,file.path(out,'cutoff_numbers.tex'));writeLines(capture.output(sessionInfo()),file.path(out,'R-session.txt'))
stopifnot(nrow(s)==60L,all(s$auc_lo>=0),all(s$auc_hi<=1),all(s$n_z<=s$n_complexes))
cat('60 class/potential/cutoff summaries; 2000 paired epitope-cluster bootstrap draws per class.\n')

# Preserve historical values and directional differences; protocol/population changes are explicit.
old<-rbindlist(lapply(2:3,function(i){d<-fread(a[i]);rbindlist(lapply(c('Fixed 5 Å','Re-derived','Published'),function(pot){suffix<-switch(pot,'Fixed 5 Å'='','Re-derived'='_rederived','Published'='_legacy');data.table(baseline=if(i==2)'original pooled' else 'pooled replay',cutoff=d$cutoff,potential=pot,old_auc=d[[paste0('auc',suffix)]],old_negative_mean_Z=-d[[paste0('mean_Z',suffix)]],old_n_z=d$n)}))}))
v<-merge(s,old,by=c('cutoff','potential'),allow.cartesian=TRUE)
v[,`:=`(delta_auc=auc-old_auc,delta_negative_mean_Z=negative_mean_Z-old_negative_mean_Z,population_matched=FALSE,
 comparison='Class separation, corrected positions, full-peptide shuffles and twenty-decoy AUC; not a matched gain/loss')]
fwrite(v,file.path(out,'version_differences.tsv'),sep='\t')
