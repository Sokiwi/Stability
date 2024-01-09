# the following code can be skipped down to
# just before where it says
# get stability values for each feature
# which is where the created object can be loaded
# from a saved file

# For purposes of computing stabilities the WALS genera and
# language-level glottocode representatives are the
# units of analysis; so all languages in Grambank have to be
# assigned a genus
gb_meta <- read.csv(file="languages_gb.csv")
genus <- rep(NA, nrow(gb_meta))
gb_meta <- cbind(gb_meta, genus)
# remove unwanted "families" from Grambank
# The Language_ID column in the values.csv file corresponds to
# the column ID in languages_gb.csv (originally called languages.csv)
# so the Language_ID's of the excluded languages need to be noted
# such that the values table can be purged of those
IDs_excluded <- gb_meta$ID[which(gb_meta$Family_level_ID%in%c("", "arti1236", "book1242", "mixe1287", "pidg1258", "spee1234", "unat1236", "uncl1493", "sign1238"))]
gb_meta <- gb_meta[!gb_meta$Family_level_ID%in%c("", "arti1236", "book1242", "mixe1287", "pidg1258", "spee1234", "unat1236", "uncl1493", "sign1238"),]
# download.file("https://zenodo.org/record/7079637/files/lexibank/asjp-v20.zip?download=1", "asjp.zip")
# needed is the file languages.csv
# construct ASJP metadata from the current database and a file
# created for the purpose of having metadata for future contributions
listss21 <- read.table(file="listss21_formatted.tab", header=TRUE, sep="\t", quote="", na.strings="", comment.char="")
metadata_missing <- read.table(file="metadata_missing.tab", header=TRUE, sep="\t", quote="", na.strings="", comment.char="")
# download languoid.csv from https://glottolog.org/meta/downloads
asjp_meta <- rbind(listss21, metadata_missing)[1:10]
u <- unique(asjp_meta$iso)
asjp_meta <- asjp_meta[match(u, asjp_meta$iso),]
languoid <- read.csv("languoid.csv")
# go through the Grambank metadata and add genera
for (i in 1:nrow(gb_meta)) {
  if (gb_meta$ID[i]=="atam1239") {
    gb_meta$genus[i] <- "NORTHERN_PAMA-NYUNGAN"
  } else if (gb_meta$ID[i]=="baha1260") {
    gb_meta$genus[i] <- "GERMANIC"
  } else if (gb_meta$ID[i]=="bera1262") {
    gb_meta$genus[i] <- "MALAYO-SUMBAWAN"
  } else if (gb_meta$ID[i]=="daai1236") {
    gb_meta$genus[i] <- "KUKI-CHIN"
  } else if (gb_meta$ID[i]=="east2485") {
    gb_meta$genus[i] <- "NORTH_BORNEO"
  } else if (gb_meta$ID[i]=="lard1243") {
    gb_meta$genus[i] <- "TANGKIC"
  } else if (gb_meta$ID[i]=="hurr1240") {
    gb_meta$genus[i] <- "HURRIAN-URARTIAN"
  } else if (gb_meta$ID[i]=="gana1278") {
    gb_meta$genus[i] <- "SOUTHEASTERN_PAMA-NYUNGAN"
  } else if (gb_meta$ID[i]=="kulo1237") {
    gb_meta$genus[i] <- "NORTHWEST_FORMOSAN"
  } else if (gb_meta$ID[i]=="oldk1249") {
    gb_meta$genus[i] <- "KHMER"
  } else if (gb_meta$ID[i]=="wint1259") {
    gb_meta$genus[i] <- "WINTUAN"
  } else if (gb_meta$ID[i]=="adon1237") {
    gb_meta$genus[i] <- "BIMA-LEMBATA"
  } else if (gb_meta$ID[i]=="aghw1237") {
    gb_meta$genus[i] <- "LEZGIC"
  } else if (gb_meta$ID[i]=="amam1246") {
    gb_meta$genus[i] <- "CORE_GOILALAN"
  } else if (gb_meta$ID[i]=="apol1242") {
    gb_meta$genus[i] <- "BOLIVIA-PARANA"
  } else if (gb_meta$ID[i]=="asur1254") {
    gb_meta$genus[i] <- "MUNDA"
  } else if (gb_meta$ID[i]=="axiy1235") {
    gb_meta$genus[i] <- "BURMESE-LOLO"
  } else if (gb_meta$ID[i]=="azha1235") {
    gb_meta$genus[i] <- "BURMESE-LOLO"
  } else if (gb_meta$ID[i]=="babi1235") {
    gb_meta$genus[i] <- "ATHAPASKAN"
  } else if (gb_meta$ID[i]=="bala1316") {
    gb_meta$genus[i] <- "OCEANIC"
  } else if (gb_meta$ID[i]=="bana1288") {
    gb_meta$genus[i] <- "NORTHERN_LUZON"
  } else if (gb_meta$ID[i]=="bang1369") {
    gb_meta$genus[i] <- "MIJI"
  } else if (gb_meta$ID[i]=="belu1239") {
    gb_meta$genus[i] <- "BANTU"
  } else if (gb_meta$ID[i]=="boka1249") {
    gb_meta$genus[i] <- "TANI"
  } else if (gb_meta$ID[i]=="bula1255") {
    gb_meta$genus[i] <- "WESTERN_PAMA-NYUNGAN"
  } else if (gb_meta$ID[i]=="buta1242") {
    gb_meta$genus[i] <- "TAULIL"
  } else if (gb_meta$ID[i]=="cent2060") {
    gb_meta$genus[i] <- "OCEANIC"
  } else if (gb_meta$ID[i]=="cent2314") {
    gb_meta$genus[i] <- "PEARIC"
  } else if (gb_meta$ID[i]=="cent2322") {
    gb_meta$genus[i] <- "UGRIC"
  } else if (gb_meta$ID[i]=="chas1234") {
    gb_meta$genus[i] <- "BURMESE-LOLO"
  } else if (gb_meta$ID[i]=="chim1312") {
    gb_meta$genus[i] <- "BANTU"
  } else if (gb_meta$ID[i]=="cosa1234") {
    gb_meta$genus[i] <- "BURMESE-LOLO"
  } else if (gb_meta$ID[i]=="east2449") {
    gb_meta$genus[i] <- "OCEANIC"
  } else if (gb_meta$ID[i]=="east2773") {
    gb_meta$genus[i] <- "HIMALAYISH"
  } else if (gb_meta$ID[i]=="east2782") {
    gb_meta$genus[i] <- "BANTU"
  } else if (gb_meta$ID[i]=="eude1234") {
    gb_meta$genus[i] <- "OPATA-EUDEVE"
  } else if (gb_meta$ID[i]=="gang1266") {
    gb_meta$genus[i] <- "KUKI-CHIN"
  } else if (gb_meta$ID[i]=="golo1234") {
    gb_meta$genus[i] <- "UBANGI"
  } else if (gb_meta$ID[i]=="gube1234") {
    gb_meta$genus[i] <- "NYUN"
  } else if (gb_meta$ID[i]=="hier1240") {
    gb_meta$genus[i] <- "ANATOLIAN"
  } else if (gb_meta$ID[i]=="icar1234") {
    gb_meta$genus[i] <- "DARGWIC"
  } else if (gb_meta$ID[i]=="japh1234") {
    gb_meta$genus[i] <- "NA-QIANGIC"
  } else if (gb_meta$ID[i]=="jede1239") {
    gb_meta$genus[i] <- "ASLIAN"
  } else if (gb_meta$ID[i]=="kenu1236") {
    gb_meta$genus[i] <- "NUBIAN"
  } else if (gb_meta$ID[i]=="khan1277") {
    gb_meta$genus[i] <- "TANGKHUL-MARING"
  } else if (gb_meta$ID[i]=="kile1243") {
    gb_meta$genus[i] <- "TUNGUSIC"
  } else if (gb_meta$ID[i]=="koro1313") {
    gb_meta$genus[i] <- "SouthBougainville"
  } else if (gb_meta$ID[i]=="lito1235") {
    gb_meta$genus[i] <- "BANTU"
  } else if (gb_meta$ID[i]=="maip1246") {
    gb_meta$genus[i] <- "NORTHERN_MAIPURAN"
  } else if (gb_meta$ID[i]=="mala1534") {
    gb_meta$genus[i] <- "NORTHERN_LUZON"
  } else if (gb_meta$ID[i]=="malg1251") {
    gb_meta$genus[i] <- "BIU-MANDARA"
  } else if (gb_meta$ID[i]=="mans1260") {
    gb_meta$genus[i] <- "HATAM-MANSIM"
  } else if (gb_meta$ID[i]=="manu1258") {
    gb_meta$genus[i] <- "CENTRAL_MALAYO-POLYNESIAN"
  } else if (gb_meta$ID[i]=="mila1245") {
    gb_meta$genus[i] <- "TANI"
  } else if (gb_meta$ID[i]=="mill1237") {
    gb_meta$genus[i] <- "HUARPE"
  } else if (gb_meta$ID[i]=="minh1238") {
    gb_meta$genus[i] <- "MONGOLIC"
  } else if (gb_meta$ID[i]=="nati1244") {
    gb_meta$genus[i] <- "OCEANIC"
  } else if (gb_meta$ID[i]=="nese1235") {
    gb_meta$genus[i] <- "OCEANIC"
  } else if (gb_meta$ID[i]=="ngal1296") {
    gb_meta$genus[i] <- "UBANGI"
  } else if (gb_meta$ID[i]=="njan1240") {
    gb_meta$genus[i] <- "MAMBILOID"
  } else if (gb_meta$ID[i]=="nuuu1241") {
    gb_meta$genus[i] <- "TU"
  } else if (gb_meta$ID[i]=="ocea1241") {
    gb_meta$genus[i] <- "OCEANIC"
  } else if (gb_meta$ID[i]=="oldm1243") {
    gb_meta$genus[i] <- "MALAYO-SUMBAWAN"
  } else if (gb_meta$ID[i]=="poly1242") {
    gb_meta$genus[i] <- "OCEANIC"
  } else if (gb_meta$ID[i]=="rund1243") {
    gb_meta$genus[i] <- "NORTH_BORNEO"
  } else if (gb_meta$ID[i]=="sirz1237") {
    gb_meta$genus[i] <- "DARGWIC"
  } else if (gb_meta$ID[i]=="sisi1250") {
    gb_meta$genus[i] <- "OCEANIC"
  } else if (gb_meta$ID[i]=="situ1238") {
    gb_meta$genus[i] <- "NA-QIANGIC"
  } else if (gb_meta$ID[i]=="sout3236") {
    gb_meta$genus[i] <- "GUMUZ"
  } else if (gb_meta$ID[i]=="surg1248") {
    gb_meta$genus[i] <- "UGRIC"
  } else if (gb_meta$ID[i]=="tang1373") {
    gb_meta$genus[i] <- "CHINESE"
  } else if (gb_meta$ID[i]=="tang1377") {
    gb_meta$genus[i] <- "TANI"
  } else if (gb_meta$ID[i]=="tief1244") {
    gb_meta$genus[i] <- "TIEFOIC"
  } else if (gb_meta$ID[i]=="tsix1234") {
    gb_meta$genus[i] <- "KHOE-KWADI"
  } else if (gb_meta$ID[i]=="tule1245") {
    gb_meta$genus[i] <- "YOKUTS"
  } else if (gb_meta$ID[i]=="umii1235") {
    gb_meta$genus[i] <- "WORRORRAN"
  } else if (gb_meta$ID[i]=="ware1255") {
    gb_meta$genus[i] <- "JAPURA-COLOMBIA"
  } else if (gb_meta$ID[i]=="winj1235") {
    gb_meta$genus[i] <- "WORRORRAN"
  } else if (gb_meta$ID[i]=="xinc1242") {
    gb_meta$genus[i] <- "XINCAN"
  } else if (gb_meta$ID[i]=="xinc1243") {
    gb_meta$genus[i] <- "XINCAN"
  } else if (gb_meta$ID[i]=="xinc1246") {
    gb_meta$genus[i] <- "XINCAN"
  } else if (gb_meta$ID[i]=="xink1235") {
    gb_meta$genus[i] <- "XINCAN"
  } else if (gb_meta$ID[i]=="yulp1239") {
    gb_meta$genus[i] <- "WESTERN_PAMA-NYUNGAN"
  } else if (gb_meta$ID[i]=="zamb1245") {
    gb_meta$genus[i] <- "BANTU"
  } else if (gb_meta$ID[i]=="lewo1244") {
    gb_meta$genus[i] <- "CENTRAL_MALAYO-POLYNESIAN"
  } else if (gb_meta$ID[i]=="sout2896") {
    gb_meta$genus[i] <- "CENTRAL_MALAYO-POLYNESIAN"
  } else {
    w_id <- which(languoid$id==gb_meta$Language_level_ID[i])
    iso <- languoid$iso639P3code[w_id]
    if (iso!="") {
      w_iso <- which(asjp_meta$iso==iso)
      if (length(w_iso) > 0) {
        gen <- asjp_meta$wls_gen[w_iso]
        if (gen != "GENUS") {
          gb_meta$genus[i] <- gen
        } else {
          s <- strsplit(gb_meta$lineage[i], "\\/")[[1]]
          subgroup_level1 <- s[length(s)-1]
          w_subgroup <- grep(subgroup_level1, gb_meta$lineage)
          genera <- gb_meta$genus[w_subgroup]
          genera <- genera[!is.na(genera)]
          if (length(unique(genera))==1) {
            gb_meta$genus[i] <- unique(genera)
          } else {
            s2 <- strsplit(gb_meta$lineage[i], "\\/")[[1]]
            subgroup_level2 <- s2[length(s2)-2]
            w_subgroup <- grep(subgroup_level2, gb_meta$lineage)
            genera <- gb_meta$genus[w_subgroup]
            genera <- genera[!is.na(genera)]
            if (length(unique(genera))==1) {
             gb_meta$genus[i] <- unique(genera)
            }
          }
        }
      } else if (length(w_iso) == 0) {
          s <- strsplit(gb_meta$lineage[i], "\\/")[[1]]
          subgroup_level1 <- s[length(s)-1]
          w_subgroup <- grep(subgroup_level1, gb_meta$lineage)
          genera <- gb_meta$genus[w_subgroup]
          genera <- genera[!is.na(genera)]
          if (length(unique(genera))==1) {
            gb_meta$genus[i] <- unique(genera)
          } else {
            s2 <- strsplit(gb_meta$lineage[i], "\\/")[[1]]
            subgroup_level2 <- s2[length(s2)-2]
            w_subgroup <- grep(subgroup_level2, gb_meta$lineage)
            genera <- gb_meta$genus[w_subgroup]
            genera <- genera[!is.na(genera)]
            if (length(unique(genera))==1) {
              gb_meta$genus[i] <- unique(genera)
            }
          }
      }
    }
  }
}

# prepare data for values
values <- read.csv(file="values.csv", header=TRUE)
values <- values[!values$Language_ID%in%IDs_excluded,]
values <- values[values$Value!="?",]

# select the best attested doculect representing a language-level languoid;
# first add a column to gb_meta indicating the number of attested values
attestation <- rep(NA, nrow(gb_meta))
gb_meta <- cbind(gb_meta, attestation)
for (i in 1:length(attestation)) {
  gb_meta$attestation[i] <- length(which(values$Language_ID==gb_meta$ID[i]))
}
# make a list of the row indices and other info of recurring language-level IDs
groups <- list()
L_l_ID <- unique(gb_meta$Language_level_ID)
count <- 0
for (i in seq_along(L_l_ID)) {
  m <- which(gb_meta$Language_level_ID==L_l_ID[i])
  if (length(m) > 1) {
    count <- count + 1
    groups[[count]] <- list(L_l_ID[i], m, gb_meta$ID[m], gb_meta$Name[m], gb_meta$attestation[m])
  }
}
# get a vector of doculects to exclude, first by taking out the doculect to
# include and then  make a vector of the doculects to exclude
for (i in seq_along(groups)) {
  w_max <- which(groups[[i]][[5]]==max(groups[[i]][[5]]))[1]
  groups[[i]][[3]] <- groups[[i]][[3]][-w_max]
}
exclude <- c()
for (i in seq_along(groups)) {
  exclude <- c(exclude, groups[[i]][[3]])
}
gb_meta <- gb_meta[!gb_meta$ID%in%exclude,]
values <- values[!values$Language_ID%in%exclude,]

# take out unnecessary columns from the values object and add Parameter_ID,
# and names of language, genus, and family
values <- values[,c("ID", "Language_ID", "Parameter_ID", "Value", "Code_ID")]
parameters <- read.csv(file="parameters.csv", header=TRUE)
get_parameter <- function(x) parameters$Name[which(parameters$ID==x)]
parameter_description <- unlist(lapply(values$Parameter_ID, get_parameter))
get_language_name <- function(x) gb_meta$Name[which(gb_meta$ID==x)]
language_name <- unlist(lapply(values$Language_ID, get_language_name))
get_genus_name <- function(x) gb_meta$genus[which(gb_meta$ID==x)]
genus_name <- unlist(lapply(values$Language_ID, get_genus_name))
get_family_name <- function(x) gb_meta$Family_name[which(gb_meta$ID==x)]
family_name <- unlist(lapply(values$Language_ID, get_family_name))
values <- data.frame(values, parameter_description, language_name, genus_name, family_name)

# take out languages whose WALS genus is CREOLES_AND_PIDGINS
w_c_and_p <- which(values$genus_name=="CREOLES_AND_PIDGINS")
if ( length(w_c_and_p) > 0 ) {
  values <- values[-w_c_and_p,]
}

save(values, file="values.RData")
rm(list=ls())
load("values.RData")  # loading the object values
write.table(values, file="Grambank_stab_data.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# function to verify that genus names are unique, i.e. the same name doesn't
# occur in different families
check_genera <- function() {
  gens <- unique(values$genus_name)
  for (i in 1:length(gens)) {
    w_gen <- which(values$genus_name==gens[i])
    fams_gen <- values$family_name[w_gen]
    u_fams_gen <- unique(fams_gen)
    if (length(u_fams_gen) == 1) {
      # cat(gens[i], "is ok\n")
    } else {
      cat("***ATTENTION***", gens[i], "recurs in", unique(fams_gen), "\n")
    }
  }
}
# can be run as check_genera(); done and verified that genus names are unique

# get stability values for each feature

# function for calculating stability given a feature
stab <- function(feat) {
  # function for getting the proportion of pairs of languages sharing 
  # values within a genus
  Rgen <- function(gen, feat) {
    dat <- values[which(values$genus_name==gen & values$Parameter_ID==feat),]
    L <- nrow(dat)
    if ( L < 2 ) {
      return(NA)
    } else {
      same <- 0
      different <- 0
      n <- 0
      for (i in 1:(L-1)) {
        for (j in (i+1):L) {
          n <- n + 1
          if (dat$Value[i]==dat$Value[j]) {
            same <- same + 1
          } else {
            different <- different + 1 
          }
        }
      }
    }
    prop <- same/n
    weight <- sqrt(n)
    weighted_prop <- prop*weight
    return(c(weighted_prop, weight))
  }
  
  # function for getting the proportion of unrelated languages 
  # sharing values
  Uall <- function(feat) {
    dat <- values[which(values$Parameter_ID==feat),]
    L <- nrow(dat)
    if ( L==0 ) {
      return(NA)
    } else {
      same <- 0
      different <- 0
      for (i in 1:(L-1)) {
        for (j in (i+1):L) {
          if (dat$family_name[i]!=dat$family_name[j]) {
            if (dat$Value[i]==dat$Value[j]) {
              same <- same + 1
            } else {
              different <- different + 1 
            }
          }
        }
      }
      prop <- same/(same + different)
      return(prop)
    }
  }
  
  # function for running Rgen on all genera and outputting the weighted average
  run_Rgen <- function(feat) {
    gens <- unique(values$genus_name)
    numerator <- 0
    sum_weights <- 0
    for (i in 1:length(gens)) {
      Rgen_out <- Rgen(gens[i], feat)
      if ( !is.na(Rgen_out[1]) ) {
        numerator <- numerator + Rgen_out[1]
        sum_weights <- sum_weights + Rgen_out[2]
      }
    }
    weighted_average <- numerator/sum_weights
    return(weighted_average)
  }
  
  # getting stability based on a feature, R, and U
  R <- run_Rgen(feat)
  U <- Uall(feat)
  cat("R: ", R, "\n", sep="")
  cat("U: ", U, "\n", sep="")
  S <- (R - U)/(1 - U)
  S <- round(S * 100, 1)
  return(S)
}

# function for calculating stabilities for all features
runall <- function() {
  feats <- unique(values$Parameter_ID)
  descriptions <- values$parameter_description[match(feats, values$Parameter_ID)]
  cat("feature_ID\tfeature_description\tstability\n", file="Stabilities_Grambank_features.txt")
  for (i in 1:length(feats)) {
    feat <- feats[i]
    description <- descriptions[i]
    cat("doing ", i, " out of ", length(feats), ":\n", sep="")
    cat("description: ", description, "\n", sep="")
    stability <- stab(feat)
    cat("stability: ", stability, "\n\n", sep="")
    cat(feat, "\t", description, "\t", stability, "\n", file="Stabilities_Grambank_features.txt", sep="", append=TRUE)
  }
}
# do runall() to calculate stabilities for all features

# analyzing the output
stabs_g <- read.table(file="Stabilities_Grambank_features.txt", header=TRUE, sep="\t", quote="", comment.char="", na.strings="")
stabs_w <- read.table(file="Stabilities_WALS_features.txt", header=TRUE, sep="\t", quote="", comment.char="", na.strings="")
summary(stabs_g$stability)
summary(stabs_w$stability)
library(ggplot2)
stabs_g$dataset <- "Grambank"
stabs_w$dataset <- "WALS"
stabs <- rbind(stabs_w[,c("stability", "dataset")], stabs_g[,c("stability", "dataset")])
ggplot(stabs, aes(stability, fill = dataset)) +
  geom_density(alpha = 0.2) + 
  xlim(-50, 100) # +
  # theme_bw()

# plot the interplay of R, U, and S
# first make a matrix of S for U and R values running from .1 to 1
s <- seq(.1, 1, .1)
L <- length(s)
m <- matrix(NA, nrow=L, ncol=L, dimnames=list(s,s))
for (i in seq_along(s)) {
  for (j in seq_along(s)) {
    m[i,j] <- (s[i]-s[j])/(1-s[j])
  }
}

# now plot R, U, and S
for (i in 1:10) {
  plot(as.numeric(colnames(m)), m[i,], type="l", xlim=c(0.1,1), ylim=c(-8,1), 
       ylab="S", xlab="U", main="Interplay of R, U, and S", lwd=1.5)
  text(x=0.94, y=m[i,9], paste("R=",colnames(m)[i]))
  grid()
  par(new=TRUE)
}

# prepare boxplot of different functional categories
library(readxl)
gs <- read_excel("Grambank_stab_data_analyzed.xlsx", sheet = 1)
sort(table(gs$functional_category),decreasing=TRUE)  # the most frequent functional categories
argument_marking <- gs$stability[gs$functional_category=="argument marking"]
classification <- gs$stability[gs$functional_category=="classification"]
word_order <- gs$stability[gs$functional_category=="word order"]
number <- gs$stability[gs$functional_category=="number"]
tense.aspect.mood <- gs$stability[gs$functional_category=="tense-aspect-mood"]
possession <- gs$stability[gs$functional_category=="possession"]
func.cats <- list(argument_marking, classification, word_order, number, tense.aspect.mood, possession)
names.func.cats <- c("argument\nmarking", "classification", "word\norder", "number", "tense\n-aspect\n-mood", "possession")
Lfc <- length(func.cats)
for (i in 1:(Lfc-1)) {
  for (j in (i+1):Lfc) {
    p <- t.test(func.cats[[i]], func.cats[[j]])$p.value
    cat(p, "\n")
  }
}

# make boxplot of different functional categories
N <- sapply(func.cats, length)
boxplot(func.cats[[1]], func.cats[[2]], func.cats[[3]], func.cats[[4]], func.cats[[5]], func.cats[[6]],
        ylim = c(-40,90), xlab="Functional categories", ylab="Stability", xaxt="n")
xtick <- c(1:6)
text(x=xtick, par("usr")[3], labels = names.func.cats, pos = 1, xpd = TRUE, cex = .9)
text(x=xtick, par("usr")[3], labels = paste("N=", N, sep=""), pos = 3, xpd = TRUE, cex = .8)

# prepare boxplot of different formal categories
library(readxl)
gs <- read_excel("Grambank_stab_data_analyzed.xlsx", sheet = 1)
sort(table(gs$syntagmatic_category),decreasing=TRUE)  # the most frequent functional categories
morphology <- gs$stability[gs$syntagmatic_category=="morphology"]
syntax <- gs$stability[gs$syntagmatic_category=="syntax"]
lexicon <- gs$stability[gs$syntagmatic_category=="lexicon"]
form.cats <- list(morphology, syntax, lexicon)
names.form.cats <- c("morphology", "syntax", "lexicon")
Lfc <- length(form.cats)
for (i in 1:(Lfc-1)) {
  for (j in (i+1):Lfc) {
    p <- t.test(form.cats[[i]], form.cats[[j]])$p.value
    cat(names.func.cats[i], names.func.cats[j], p, "\n")
  }
}

# make boxplot of different formal categories
N <- sapply(form.cats, length)
boxplot(form.cats[[1]], form.cats[[2]], form.cats[[3]],
        ylim = c(-40,90), xlab="Formal categories", ylab="Stability", xaxt="n")
xtick <- c(1:3)
text(x=xtick, par("usr")[3], labels = names.form.cats, pos = 1, xpd = TRUE, cex = .9)
text(x=xtick, par("usr")[3], labels = paste("N=", N, sep=""), pos = 3, xpd = TRUE, cex = .8)




