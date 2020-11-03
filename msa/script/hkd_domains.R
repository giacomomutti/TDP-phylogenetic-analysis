# dependencies
library(msa)

# this script produces figure of HKD alignment
# unluckily it's not possible to show both subset together
# with \setdomain as this always exceed memory so
# we plot the domains separately and then we then use
# inkscape

# read msa

a <- readRDS("msa/alpha/TDP1alpha_PP.rds")
b <- readRDS("msa/beta/TDP1beta_PP.rds")

# further code for latex figure
further <- c(
  "\\namesrm\\namessl",
  "\\shownumbering{leftright}",
  "\\namesequencelogo{Logo}",
  "\\defconsensus{.}{lower}{upper}",
  "\\showruler{1}{top}",
  "\\bargraphstretch{3}",
  "\\featurestylenamescolor{RoyalBlue}",
  "\\showfeaturestylename{bottom}{Charge}",
  "\\showfeaturestylename{bbottom}{Conservation}",
  "\\showfeaturestylename{bbbottom}{Hydrophob.}",
  "\\showfeaturestylename{bbbbottom}{Weight}",
  "\\bottomspace{0.2\\baselineskip}",
  "\\bbottomspace{-\\baselineskip}",
  "\\bbbottomspace{-0.9\\baselineskip}",
  "\\bbbbottomspace{-\\baselineskip}"
)

# unluckily this function uses texi2pdf which always
# saves in working directory so we'll have to move the files

a_HKD1_prettyprint <- msa::msaPrettyPrint(a,
  file = "msa/alpha/HKD1_alpha.tex",
  output = "tex",
  showNames = "left",
  paperWidth = 20,
  paperHeight = 17,
  showLogo = "top",
  consensusThreshold = c(50, 90),
  logoColors = "rasmol",
  showLogoScale = "none",
  shadingMode = "similar",
  showConsensus = "bottom",
  showLegend = T,
  askForOverwrite = FALSE,
  furtherCode = c(
    further,
    "\\setends{38}{261..291}",
    "\\feature{ttop}{38}{261..291}{brace}{HKD motif I}",
    "\\feature{bottom}{38}{261..291}{color:charge[BlueRed]}{}",
    "\\feature{bbottom}{38}{261..291}{bar:conservation[JungleGreen,Gray10]}{}",
    "\\feature{bbbottom}{38}{261..291}{bar:hydrophobicity[RoyalBlue,Gray10]}{}",
    "\\feature{bbbbottom}{38}{261..291}{bar:molweight[RedViolet,Gray10]}{}"
  )
)

a_HKD2_prettyprint <- msa::msaPrettyPrint(a,
  file = "msa/alpha/HKD2_alpha.tex",
  output = "tex",
  verbose = TRUE,
  showNames = "left",
  paperWidth = 20,
  paperHeight = 17,
  showLogo = "top",
  consensusThreshold = c(50, 90),
  logoColors = "rasmol",
  showLogoScale = "right",
  shadingMode = "similar",
  showConsensus = "bottom",
  showLegend = T,
  askForOverwrite = FALSE,
  furtherCode = c(
    further,
    "\\setends{38}{491..525}",
    "\\feature{ttop}{38}{491..525}{brace}{HKD motif II}",
    "\\feature{bottom}{38}{491..525}{color:charge[BlueRed]}{}",
    "\\feature{bbottom}{38}{491..525}{bar:conservation[JungleGreen,Gray10]}{}",
    "\\feature{bbbottom}{38}{491..525}{bar:hydrophobicity[RoyalBlue,Gray10]}{}",
    "\\feature{bbbbottom}{38}{491..525}{bar:molweight[RedViolet,Gray10]}{}"
  )
)


b_HKD1_prettyprint <- msa::msaPrettyPrint(b,
  file = "msa/beta/HKD1_beta.tex",
  output = "tex",
  verbose = TRUE,
  showNames = "left",
  paperWidth = 20,
  paperHeight = 17,
  showLogo = "top",
  consensusThreshold = c(50, 90),
  logoColors = "rasmol",
  showLogoScale = "none",
  shadingMode = "similar",
  showConsensus = "bottom",
  showLegend = T,
  askForOverwrite = FALSE,
  furtherCode = c(
    further,
    "\\setends{14}{471..499}",
    "\\feature{ttop}{14}{471..499}{brace}{HKD motif I}",
    "\\feature{bottom}{14}{471..499}{color:charge[BlueRed]}{}",
    "\\feature{bbottom}{14}{471..499}{bar:conservation[JungleGreen,Gray10]}{}",
    "\\feature{bbbottom}{14}{471..499}{bar:hydrophobicity[RoyalBlue,Gray10]}{}",
    "\\feature{bbbbottom}{14}{471..499}{bar:molweight[RedViolet,Gray10]}{}"
  )
)

b_HKD1_prettyprint <- msa::msaPrettyPrint(b,
  file = "msa/beta/HKD2_beta.tex",
  output = "tex",
  verbose = TRUE,
  showNames = "left",
  paperWidth = 20,
  paperHeight = 17,
  showLogo = "top",
  consensusThreshold = c(50, 90),
  logoColors = "rasmol",
  showLogoScale = "right",
  shadingMode = "similar",
  showConsensus = "bottom",
  showLegend = T,
  askForOverwrite = FALSE,
  furtherCode = c(
    further,
    "\\setends{14}{890..923}",
    "\\feature{ttop}{14}{890..923}{brace}{HKD motif II}",
    "\\feature{bottom}{14}{890..923}{color:charge[BlueRed]}{}",
    "\\feature{bbottom}{14}{890..923}{bar:conservation[JungleGreen,Gray10]}{}",
    "\\feature{bbbottom}{14}{890..923}{bar:hydrophobicity[RoyalBlue,Gray10]}{}",
    "\\feature{bbbbottom}{14}{890..923}{bar:molweight[RedViolet,Gray10]}{}"
  )
)


# then manueal editing of tex file to do definitive figures

# unusual alpha HKD II species:
sp <- c("sapiens", "truncatula", "clementina", "comosus", "pruriens", "conductrix")
vec_sp <- as.vector(sapply(sp, grep, names(a@unmasked)))


a_HKD2_prettyprint_unsual <- msa::msaPrettyPrint(a,
  subset = vec_sp,
  file = "msa/alpha/HKD2_unusual_alpha.tex",
  output = "tex",
  verbose = TRUE,
  showNames = "left",
  paperWidth = 12,
  paperHeight = 8.3,
  showLogo = "top",
  consensusThreshold = c(50, 90),
  logoColors = "rasmol",
  showLogoScale = "right",
  shadingMode = "similar",
  showConsensus = "bottom",
  showLegend = T,
  askForOverwrite = FALSE,
  furtherCode = c(
    further,
    "\\setends{2}{440..551}",
    "\\feature{ttop}{2}{491..525}{brace}{HKD motif II}",
    "\\feature{bottom}{2}{440..551}{color:charge[BlueRed]}{}",
    "\\feature{bbottom}{2}{440..551}{bar:conservation[JungleGreen,Gray10]}{}",
    "\\feature{bbbottom}{2}{440..551}{bar:hydrophobicity[RoyalBlue,Gray10]}{}",
    "\\feature{bbbbottom}{2}{440..551}{bar:molweight[RedViolet,Gray10]}{}"
  )
)

totlist <- list.files(path = c('./msa/beta','./msa/alpha'), pattern = '*pdf', full.names = T)

# alist <- list.files(path = "./", pattern = "HKD(1|2|2_unusual)_alpha.tex", ignore.case = T)
# blist <- list.files(path = "./", pattern = "HKD(1|2)_beta.tex", ignore.case = T)
# 
# file.copy(alist, to = "msa/alpha/")
# file.copy(blist, to = "msa/beta/")
# 
# file.remove(totlist)

lapply(totlist, function(x)
       pdftools::pdf_convert(x,
       format = 'png', dpi = 400,
       filenames = gsub("pdf","png",x)))
