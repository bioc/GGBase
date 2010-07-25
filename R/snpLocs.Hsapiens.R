snpLocs.Hsapiens = function( rsid, chrtok, spack="SNPlocs.Hsapiens.dbSNP.20090506" ) {
 require(spack, character.only=TRUE)
 ldf = getSNPlocs(chrtok)
 if (nrow(ldf) == 0) stop("chrtok must be wrong")
 locs = ldf$loc
 names(locs) = paste("rs", as.character(ldf$RefSNP_id), sep="")
 chk = intersect(rsid, names(locs))
 if (length(chk) != length(rsid)) warning(paste("some SNP in rsid were not found in location db", spack))
 locs[ chk ]
}
 
 
