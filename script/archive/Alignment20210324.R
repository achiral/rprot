#Alignment.R
#DECIPHER
#http://www.iu.a.u-tokyo.ac.jp/~kadota/20140612_kadota.pdf
#http://www.iu.a.u-tokyo.ac.jp/~kadota/r_seq.html#about_analysis_general_alignment_multiple
#Biostrings
#http://www.iu.a.u-tokyo.ac.jp/~kadota/bioinfo_ngs_sokushu_2014/20140909_3-4_kadota.pdf
#msa
#https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
####################################################################################################################
#if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("DECIPHER")
#BiocManager::install("msa")
####################################################################################################################
#directory choose関数の定義
dir.choose <- function() {
  system("osascript -e 'tell app \"RStudio\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
         intern = FALSE, ignore.stderr = TRUE)
  p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
  return(ifelse(length(p), p, NA))
}
setwd(dir.choose())
getwd()  #directoryの確認
dir()    #directoryの内容確認
####################################################################################################################
in_f <- "CYP1A2m1F.fasta"                  #入力ファイル名を指定してin_fに格納
out_f <- "out.fasta"                 #出力ファイル名を指定してout_fに格納

#必要なパッケージをロード
library(DECIPHER)                      #パッケージの読み込み

#入力ファイルの読み込み
fasta <- readAAStringSet(in_f, format="fasta")  #in_fで指定したファイルの読み込み
fasta                                  #確認してるだけです

#本番(MSA)
out <- AlignSeqs(fasta)                #MSAを実行した結果をoutに格納
out                                    #確認してるだけです

#ファイルに保存
writeXStringSet(out, file=out_f, format="fasta", width=50)#outの中身を指定したファイル名で保存
####################################################################################################################
library(Biostrings) 
in_f <- "ABCB1AOKES016F1-PREMIX.fasta"
out_f <- "out.fasta"
fasta <- readDNAStringSet(in_f, format="fasta")
writeXStringSet(fasta, file=out_f, format="fasta")
#hoge <- translate(fasta)
#names(hoge) <- names(fasta)
#writeXStringSet(hoge, file=out_f, format="fasta")
writeXStringSet(out, file=out_f, format="fasta", width=50)#outの中身を指定したファイル名で保存
####################################################################################################################
#msa
####################################################################################################################
library(msa)
in_f <- "sequence.txt"
mySequenceFile <- system.file("examples", in_f, package="msa")
mySequences <- readAAStringSet(in_f)
mySequences
myFirstAlignment <- msa(mySequences)
myFirstAlignment
print(myFirstAlignment, show="complete")
msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none", showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
msaPrettyPrint(myFirstAlignment, y=c(164, 213), output="asis", showNames="none", showLogo="none", askForOverwrite=FALSE)
myClustalWAlignment <- msa(mySequences, "ClustalW")
myClustalWAlignment
myClustalOmegaAlignment <- msa(mySequences, "ClustalOmega")
myClustalOmegaAlignment
myMuscleAlignment <- msa(mySequences, "Muscle")
myMuscleAlignment
print(myFirstAlignment)
print(myFirstAlignment, show="complete")
print(myFirstAlignment, showConsensus=FALSE, halfNrow=3)
print(myFirstAlignment, showNames=FALSE, show="complete")
myMaskedAlignment <- myFirstAlignment 
colM <- IRanges(start=1, end=100) 
colmask(myMaskedAlignment) <- colM 
myMaskedAlignment
unmasked(myMaskedAlignment)
conMat <- consensusMatrix(myFirstAlignment) 
dim(conMat)
conMat <- consensusMatrix(myMaskedAlignment) 
conMat[, 95:104]
printSplitString <- function(x, width=getOption("width") - 1){
  starts <- seq(from=1, to=nchar(x), by=width)
  for (i in 1:length(starts))
    cat(substr(x, starts[i], starts[i] + width - 1), "\n")
  }
printSplitString(msaConsensusSequence(myFirstAlignment))
printSplitString(msaConsensusSequence(myFirstAlignment, type="upperlower", thresh=c(40, 20)))
printSplitString(msaConsensusSequence(myMaskedAlignment, type="upperlower", thresh=c(40, 20)))

data(BLOSUM62)
msaConservationScore(myFirstAlignment, BLOSUM62)
msaConservationScore(myFirstAlignment, BLOSUM62, gapVsGap=0, type="upperlower", thresh=c(40, 20))
msaConservationScore(myMaskedAlignment, BLOSUM62, gapVsGap=0, type="upperlower", thresh=c(40, 20))

hemoSeq <- readAAStringSet(system.file("examples/HemoglobinAA.fasta", package="msa"))
hemoAln <- msa(hemoSeq)
hemoAln
hemoAln2 <- msaConvert(hemoAln, type="seqinr::alignment")


install.packages("seqinr")
library(seqinr)
d <- dist.alignment(hemoAln2, "identity") 
as.matrix(d)[2:5, "HBA1_Homo_sapiens", drop=FALSE]


install.packages("ape")
library(ape)
hemoTree <- nj(d)
plot(hemoTree, main="Phylogenetic Tree of Hemoglobin Alpha Sequences")
hemoAln3 <- msaConvert(hemoAln, type="bios2mds::align") 
str(hemoAln3)
hemoAln4 <- as(hemoAln, "BStringSet") 
hemoAln4

msaPrettyPrint(myFirstAlignment, output="asis", y=c(164, 213), subset=c(1:6), showNames="none", showLogo="none",
               consensusColor="ColdHot", showLegend=FALSE, askForOverwrite=FALSE)
msaPrettyPrint(myFirstAlignment, output="asis", y=c(164, 213), subset=c(1:6), showNames="none", showLogo="top",
               logoColors="rasmol", shadingMode="similar", showLegend=FALSE, askForOverwrite=FALSE)
msaPrettyPrint(myFirstAlignment, output="asis", y=c(164, 213), showNames="none", shadingMode="similar",
               shadingColors="blues", showLogo="none", showLegend=FALSE, askForOverwrite=FALSE)
msaPrettyPrint(myFirstAlignment, output="asis", y=c(164, 213), showNames="none", shadingMode="functional",
               shadingModeArg="structure", askForOverwrite=FALSE)
msaPrettyPrint(myFirstAlignment, output="asis", y=c(164, 213), 
               subset=c(1:6), showNames="none", showLogo="none",
               consensusColor="ColdHot", showLegend=FALSE, 
               shadingMode="similar", askForOverwrite=FALSE, 
               furtherCode=c("\\defconsensus{.}{lower}{upper}","\\showruler{1}{top}"))
## how much fits on one page depends on the length of names and the number of sequences;
## change to what suits your needs
chunkSize <- 300 

for (start in seq(1, ncol(aln), by=chunkSize)){
  end <- min(start + chunkSize - 1, ncol(aln))
  alnPart <- DNAMultipleAlignment(subseq(unmasked(aln), start, end))
  msaPrettyPrint(x=alnPart, output="pdf", subset=NULL, file=paste0("aln_", start, "-", end, ".pdf"))
}

toBibtex(citation("msa"))


