# Established tied ranking and continuous-response summaries; explicit input contract.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(1L);setNumericRounding(0L)
a<-commandArgs(TRUE);stopifnot(length(a)==3L)
config<-fromJSON(a[1],simplifyVector=FALSE);out<-a[2]
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
source(file.path(a[3],'tsv_exact.R'));source(file.path(a[3],'ranking_metrics.R'))
write<-function(d,name)write_exact_tsv(d,file.path(out,paste0(name,'.tsv')))
metrics<-curves<-points<-counts<-rescues<-flags<-stats<-list()
for(spec in config$panels){
 ds<-spec$name;d<-read_exact_tsv(spec$rows)
 for(col in unlist(spec$required_true))d<-d[get(col)%in%TRUE]
 for(col in unlist(spec$required_false))d<-d[get(col)%in%FALSE]
 stopifnot(!anyDuplicated(d[,.(id,epitope)]))
 if(!'native_structure_exists'%in%names(d))d[,native_structure_exists:=NA]
 if(!'graded'%in%names(d))d[,graded:=NA_real_]
 if(!is.null(spec$response_column))d[,graded:=get(spec$response_column)]
 if(!is.null(spec$matrices)){
  for(clone in names(spec$matrices)){
   mat<-fread(spec$matrices[[clone]])
   aa<-setdiff(names(mat),'position');m<-as.matrix(mat[,..aa]);shifted<-m-apply(m,1,min)
   ii<-which(d$epitope==clone & nchar(d$id)==nrow(m))
   d$graded[ii]<-vapply(d$id[ii],function(p){
    chars<-strsplit(p,'',fixed=TRUE)[[1]]
    stopifnot(length(chars)==nrow(m),all(chars%in%aa))
    log2(1+sum(shifted[cbind(seq_along(chars),match(chars,aa))]))
   },numeric(1))
  }
  write(d[!is.finite(graded),.(id,epitope,y,iptm,reason='peptide length unsupported by source matrix')],paste0(ds,'_response_unavailable'))
 }
 arms<-unlist(spec$arms);display<-unlist(spec$display)
 stopifnot(all(names(arms)%in%names(d)),all(display%in%names(arms)),all(arms%in%c(-1,1)))
for(ep in sort(unique(d$epitope))){
  part<-d[epitope==ep];coverage<-part$native_structure_exists[1]
  context<-data.table(dataset=ds,epitope=ep,native_structure_exists=coverage)
  counts[[length(counts)+1L]]<-cbind(context,data.table(n_models=nrow(part),n_called=sum(is.finite(part$y)),
    n_positive=sum(part$y==1,na.rm=TRUE),n_negative=sum(part$y==0,na.rm=TRUE),
    n_ambiguous=sum(is.na(part$y)),n_signed_response=sum(is.finite(part$graded)),
    n_screen=sum(is.finite(part$screen_score)),n_iptm=sum(is.finite(part$iptm))))
  finite_common<-Reduce(`&`,lapply(display,function(s)is.finite(part[[s]])))
  for(s in names(arms)){
   for(pop in c('pair','display_common')){
    if(pop=='display_common'&&!s%in%display)next
    finite<-if(pop=='display_common')finite_common else is.finite(part[[s]])&is.finite(part$iptm)
    take<-finite&is.finite(part$y);y<-part$y[take];x<-arms[s]*part[[s]][take];ref<-part$iptm[take]
    now<-metric(y,x);base<-metric(y,ref)
    for(m in names(now))metrics[[length(metrics)+1L]]<-cbind(context,data.table(score=s,population=pop,
      endpoint='call',metric=m,n_rows=nrow(part),n_used=sum(take),n_positive=sum(y),n_negative=sum(y==0),
      value=now[m],iptm_value=base[m],delta=now[m]-base[m]))
    if(pop=='display_common'&&sum(y==1)>0&&sum(y==0)>0){
     z<-data.table(x,y)[,.(pos=sum(y),neg=sum(y==0)),by=x][order(-x)];tp<-cumsum(z$pos);fp<-cumsum(z$neg)
     curves[[length(curves)+1L]]<-cbind(context[rep(1,nrow(z)+1L)],data.table(score=s,n_used=length(y),
       n_positive=sum(y),recall=c(0,tp/sum(y)),precision=c(1,tp/(tp+fp)),fpr=c(0,fp/sum(y==0))))
    }
    # Continuous endpoints do not discard ambiguous classification calls.
    take<-finite&is.finite(part$graded)
    if(sum(take)>2L){
     x<-arms[s]*part[[s]][take];ref<-part$iptm[take];response<-part$graded[take]
     if(length(unique(x))>1L&&length(unique(response))>1L){
      for(method in c('pearson','spearman')){
       value<-cor(x,response,method=method);baseline<-cor(ref,response,method=method)
       metrics[[length(metrics)+1L]]<-cbind(context,data.table(score=s,population=pop,
        endpoint=spec$response_name,metric=method,n_rows=nrow(part),
        n_used=sum(take),n_positive=NA_integer_,n_negative=NA_integer_,value=value,
        iptm_value=baseline,delta=value-baseline))
      }
     }
     if(pop=='display_common'){
      points[[length(points)+1L]]<-cbind(context[rep(1,sum(take))],
       data.table(row_key=part$id[take],score=s,x=x,graded=response,y=part$y[take]))
      # Preserve the existing manuscript's Pearson zero-correlation test; no new contrast.
      test<-cor.test(x,response,method='pearson')
      stats[[length(stats)+1L]]<-cbind(context,data.table(score=s,n_used=sum(take),
        pearson_R=unname(test$estimate),pearson_P=test$p.value))
     }
    }
   }
  }
  # Threshold-free concordance: each positive/negative pair is one comparison.
  take<-is.finite(part$y)&is.finite(part$screen_score)&is.finite(part$iptm)
  p<-part[take];pos<-p$y==1;neg<-p$y==0
  ip<-outer(p$iptm[pos],p$iptm[neg],'-');sc<-outer(p$screen_score[pos],p$screen_score[neg],'-')
  rescues[[length(rescues)+1L]]<-cbind(context,data.table(n_models=nrow(p),n_positive=sum(pos),n_negative=sum(neg),
    n_pairs=length(ip),rescued_strict=sum(ip<0&sc>0),lost_strict=sum(ip>0&sc<0),
    iptm_tie_resolved_correct=sum(ip==0&sc>0),iptm_tie_resolved_wrong=sum(ip==0&sc<0),
    harmonized_tie_from_correct=sum(ip>0&sc==0),harmonized_tie_from_wrong=sum(ip<0&sc==0)))
  fcols<-grep('^flag_',names(part),value=TRUE)
  for(f in fcols)flags[[length(flags)+1L]]<-cbind(context,data.table(flag=f,n_models=nrow(part),
    n_flagged=sum(part[[f]]%in%TRUE),n_available=sum(!is.na(part[[f]]))))
 }
 write(d,paste0(ds,'_rows'))
}
m<-rbindlist(metrics,fill=TRUE);m[,normalized_auc:=ifelse(metric=='roc_auc',2*value-1,NA_real_)];write(m,'metrics')
for(name in c('curves','points','counts','rescues','flags','stats'))write(rbindlist(get(name),fill=TRUE),name)
writeLines(capture.output(sessionInfo()),file.path(out,'R-session.txt'))
