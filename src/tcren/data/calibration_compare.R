# Exact matched historical calibration comparisons; supplied energies only.
suppressPackageStartupMessages(library(data.table));setDTthreads(1L)
a<-commandArgs(TRUE);stopifnot(length(a)==3L);out<-a[3]
read_energy<-function(path){d<-fread(path,colClasses=c(energy='character'));d[,energy:=as.numeric(energy)];d}
prev<-read_energy(a[1]);cur<-read_energy(a[2])
shared<-merge(unique(prev[,.(pdb.id,epitope,mhc_class)]),unique(cur[,.(pdb.id,epitope,mhc_class)]),by=c('pdb.id','epitope','mhc_class'))
fwrite(shared,file.path(out,'matched_previous_cohort.tsv'),sep='\t')
ans<-list()
for(cl in c('MHCI','MHCII'))for(cut in sort(unique(cur$cutoff))){
 ids<-sort(shared[mhc_class==cl,pdb.id]);n<-length(ids);blocks<-list()
 for(version in c('previous','structural'))for(pot in c('Fixed 5 Å','Re-derived','Published')){
  x<-if(version=='previous')prev else cur
  d<-dcast(x[mhc_class==cl & cutoff==cut & potential==pot & pdb.id%in%ids],pdb.id~replicate,value.var='energy')[match(ids,pdb.id)]
  e<-as.matrix(d[,-'pdb.id']);u<-0
  for(j in 2:21){delta<-outer(-e[,1],-e[,j],'-');u<-u+sum(delta>0)+.5*sum(delta==0)}
  sd<-apply(e[,-1],1,sd);blocks[[paste(version,pot)]]<-list(auc=u/(20*n^2),z=(e[,1]-rowMeans(e[,-1]))/sd,valid=sd>1e-9)
 }
 valid<-Reduce('&',lapply(blocks,`[[`,'valid'))
 for(pot in c('Fixed 5 Å','Re-derived','Published')){
  p<-blocks[[paste('previous',pot)]];s<-blocks[[paste('structural',pot)]]
  pz<--mean(p$z[valid]);sz<--mean(s$z[valid])
  ans[[length(ans)+1L]]<-data.table(mhc_class=cl,cutoff=cut,potential=pot,n_matched=n,n_common_z=sum(valid),
   previous_auc=p$auc,structural_auc=s$auc,delta_auc=s$auc-p$auc,auc_regression=s$auc<p$auc,
   previous_negative_mean_Z=pz,structural_negative_mean_Z=sz,delta_negative_mean_Z=sz-pz,z_regression=sz<pz)
 }
}
fwrite(rbindlist(ans),file.path(out,'matched_previous_differences.tsv'),sep='\t')
