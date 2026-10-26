# Preserve distinct IEEE doubles when reading score TSVs for ranking statistics.
# Some R decimal parsers collapse adjacent doubles;
# the platform strtod preserves the original round-trip decimal representation.
Rcpp::cppFunction('Rcpp::NumericVector exact_double(Rcpp::CharacterVector x) {
 Rcpp::NumericVector y(x.size());
 for (R_xlen_t i=0; i<x.size(); ++i) {
  if (x[i]==NA_STRING) { y[i]=NA_REAL; continue; }
  std::string s=Rcpp::as<std::string>(x[i]);
  if (s.empty()) { y[i]=NA_REAL; continue; }
  char *end;
  y[i]=std::strtod(s.c_str(), &end);
  if (end==s.c_str() || *end) Rcpp::stop("invalid numeric input");
 }
 return y;
}',includes='#include <cstdlib>')
stopifnot(length(unique(exact_double(c('0.98','0.9799999999999999','0.9800000000000001'))))==3L)
read_exact_tsv<-function(path) {
 d<-data.table::fread(path)
 cols<-names(d)[vapply(d,function(x)is.double(x)&&!inherits(x,'integer64'),logical(1))]
 if(length(cols)) {
  raw<-data.table::fread(path,select=cols,colClasses='character')
  for(col in cols)data.table::set(d,j=col,value=exact_double(raw[[col]]))
 }
 d
}

# Round-trip-safe audit exports: fwrite's default decimal precision can merge
# adjacent doubles. Formatting to character preserves the original IEEE value.
write_exact_tsv<-function(d,path) {
 out<-data.table::copy(d)
 cols<-names(out)[vapply(out,function(x)is.double(x)&&!inherits(x,'integer64'),logical(1))]
 for(col in cols) {
  x<-out[[col]]
  data.table::set(out,j=col,value=ifelse(is.na(x),'NaN',sprintf('%.17g',x)))
 }
 data.table::fwrite(out,path,sep='\t',na='NaN')
}
