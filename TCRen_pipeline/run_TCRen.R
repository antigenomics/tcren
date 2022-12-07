library.path <- .libPaths()
library("data.table", lib.loc = library.path)
suppressPackageStartupMessages(library("tidyverse", lib.loc = library.path))
library("magrittr", lib.loc = library.path, warn.conflicts = FALSE)
library("optparse", lib.loc = library.path)
library("stringr", lib.loc = library.path)
 
option_list = list(
  make_option(c("-s", "--structures"), type="character", default="input_structures", 
              help="name of directory with input structures [default= %default]",
              metavar="character"),
  make_option(c("-c", "--candidates"), type="character", default="candidate_epitopes.txt", 
              help="name of file with candidate epitopes [default= %default]", 
              metavar="character"),
  make_option(c("-p", "--potential"), type="character", default="TCRen_potential.csv", 
              help="name of file with energy potential [default= %default]", 
              metavar="character"),
  make_option(c("-o", "--out"), type="character", default="output_TCRen", 
              help="name of directory with TCRen output [default= %default]",
              metavar="character"),
  make_option(c("-m", "--memory"), type="character", default="1G", 
              help="memory allocation [default= %default]",
              metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

MEMORY <- opt$memory

run_mir <- function(input_str = file.path(opt$structures, dir(opt$structures)),
                    mir_path = file.path("mir-1.0-SNAPSHOT.jar"),
                    output_dir = file.path(opt$out, "structures_annotation"),
                    arg, #"annotate-structures", "compute-pdb-geom", "compute-pdb-contacts"
                    print_log = T) {
  .pdb_list <- input_str %>% paste(collapse = " ")
  dir.create(file.path(opt$out), showWarnings = FALSE)
  cmd <- str_glue("java -Xmx{MEMORY} -cp {mir_path} com.milaboratory.mir.scripts.Examples {arg} -I {.pdb_list} -O ./{output_dir}/")
  code <- system(cmd, ignore.stdout = !print_log, ignore.stderr = !print_log)
  if(code != 0) {
    stop(str_glue("Failed to execute '{cmd}'"))
  }
}
system("java -version")

run_mir(arg = "annotate-structures")
run_mir(arg = "compute-pdb-contacts")

contacts <- fread(file.path(opt$out, "structures_annotation", "atomdist.txt")) %>%
  filter(dist <= 5, chain.id.from != chain.id.to) %>%
  select(-atom.from, -atom.to, -dist) %>%
  unique
  
contacts.full <- rbind(
  contacts,
  contacts %>% 
    rename(chain.id.from = chain.id.to,
           chain.id.to = chain.id.from,
           residue.index.from = residue.index.to,
           residue.index.to = residue.index.from)
)

general_resmarkup <- fread(file.path(opt$out, "structures_annotation", "general.txt")) %>% 
  merge(fread(file.path(opt$out, "structures_annotation", "resmarkup.txt")),
        skip.blank.lines = T) %>% 
  group_by(pdb.id, chain.id, region.type) %>% 
  mutate(region.start = min(residue.index)) %>% 
  ungroup 
  
markup <- fread(file.path(opt$out, "structures_annotation", "markup.txt"))

peptides_pdb <- markup %>% 
  filter(region.type == "PEPTIDE") %>% 
  select(pdb.id, peptide = region.sequence) %>% 
  mutate(pep.len = nchar(peptide)) %>% 
  ungroup()

contact.map <- contacts.full %>% 
  merge.data.frame(general_resmarkup %>% 
          select(pdb.id, chain.type.from = chain.type, chain.id.from = chain.id,
                 residue.index.from = residue.index, residue.aa.from = residue.aa,
                 region.type.from = region.type, region.start.from = region.start)) %>% 
  merge.data.frame(general_resmarkup %>% 
          select(pdb.id, chain.type.to = chain.type, chain.id.to = chain.id, 
                 residue.index.to = residue.index, residue.aa.to = residue.aa,
                 region.type.to = region.type, region.start.to = region.start)) %>% 
  mutate(pos.from = residue.index.from - region.start.from,
         pos.to = residue.index.to - region.start.to) %>% 
  filter(chain.type.from %in% c("TRA", "TRB"),
         chain.type.to == "PEPTIDE") %>% 
  select(pdb.id, chain.type.from, region.type.from, 
         pos.from, pos.to, residue.aa.from, residue.aa.to)

peptide_mut <- fread(file.path(opt$candidates)) %>%
  mutate(pep.len = nchar(peptide)) %>% 
  merge(peptides_pdb %>% select(pdb.id, pep.len)) %>% 
  select(-pep.len)

pep_len_max <- nchar(peptide_mut$peptide) %>% max()

contacts.mut.pep <- peptide_mut %>% 
  separate(peptide, into = as.character(c(0:(pep_len_max-1))), sep = 1:pep_len_max, remove = F) %>% 
  melt(id = c("peptide", "pdb.id"), variable.name = "pos.to", value.name = "residue.aa.to") %>%
  mutate(pos.to = as.integer(as.character(pos.to))) %>% 
  merge(contact.map %>% 
          select(-residue.aa.to), allow.cartesian=TRUE) 

potential <- fread(file.path(opt$potential)) %>% 
  melt(id = c("residue.aa.from", "residue.aa.to"), variable.name = "potential")

energy.mut.pep <- contacts.mut.pep %>% 
  merge(potential, by = c("residue.aa.from", "residue.aa.to")) %>% 
  group_by(pdb.id, peptide, potential) %>% 
  summarise(score = sum(value), .groups = "keep") %>% 
  arrange(pdb.id, score) %>%
  select(complex.id = pdb.id, peptide, potential, score)


fwrite(energy.mut.pep,
       file.path(opt$out, "candidate_epitopes_TCRen.csv"))

cat(paste0("The ranked list of candidate epitopes can be found in file ", file.path(opt$out, "candidate_epitopes_TCRen.csv"), "\n"))
