out <- c()
for(j in 2196:2440){
    f <- paste0('enc_nt',j,'.Robj')
    if(file.exists(f)){
        load(f)
    }
    out <- rbind(out,val)
}
ids <- unique(out[,'id'])
dmin <- min(out[,'dmin'])
save(out,file='combined_candidates.Robj')
