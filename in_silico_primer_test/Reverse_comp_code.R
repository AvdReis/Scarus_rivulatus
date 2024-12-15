#R code to reverse complement primers

library(mgsub)
library(stringi)

#https://www.bioinformatics.org/sms/iupac.html
stri_reverse(mgsub::mgsub("Primer_seq",
                                      c(" ","A", "T", "C", "G", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "I"), #no need to sub N
                                      c("","T", "A", "G", "C", "Y", "R", "S", "W", "M", "K", "V", "H", "D", "B", "N")))
