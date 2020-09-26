nef <- list.files(path ="colony/", pattern = "Ne", recursive = T)
out <- data.frame(pop = gsub("/.+", "", nef),
                  ne.rand = numeric(length(nef)), 
                  lcl.rand = numeric(length(nef)),
                  ucl.rand = numeric(length(nef)),
                  ne.nonrand = numeric(length(nef)),
                  lcl.nonrand = numeric(length(nef)),
                  ucl.nonrand = numeric(length(nef)))

for(i in 1:length(nef)){
  dat <- readLines(paste0("colony/", nef[i]))
  dat <- gsub(" ", "", dat)
  dat <- gsub("^.+=", "", dat)
  out[i,-1] <- as.numeric(dat[c(5,6,7,12,13,14)])
}

out[out == 2147483647] <- Inf
write.table(out, "colony_ne.txt", quote = F, sep = "\t", row.names = F)
