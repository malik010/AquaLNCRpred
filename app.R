## app.R ##
library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(markdown)
library(seqinr)
library(shinyjs)
library(protr)
#library(RWeka)
library(randomForest)
library(caret)
library(tidyverse)
library(kernlab)
library(ftrCOOL)
library(doParallel)
library(xgboost)
#library(dplyr)

ui <- dashboardPage(title="Welcome to AquaLNCRpred", skin ="green",
   dashboardHeader(title = span("AquaLNCRpred ", style = "color:black; font-size: 28px")),
 
		  
  dashboardSidebar(id="", sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Algorithm", tabName = "algorithm", icon = icon("cogs")),
      menuItem("Download", tabName = "download", icon = icon("download")),
      menuItem("Help", tabName = "help", icon = icon("question-circle")),
      menuItem("Contact", tabName = "contact", icon = icon("envelope")))), 

                dashboardBody(
		 tabItems(
                             # First tab content
                                         tabItem(tabName = "home", h1("Welcome to AquaLNCRpred", align="center", style = "color:green"),
						 p(h3(align="center", "A webserver to predict LncRNA sequences in ", em("Eriocheir sinensis"), "(Chinese Mitten Crab) and related aquaculture organisms")),
                                                fluidRow(
                                                         box(width = 11, height = 500,
                                                               withTags({
                                                         		      div(class="input-area", checked=NA,
                                                            		      p(h3("Enter input sequence(s) in FASTA format. (We highly recommend upto 200 sequences)")),

                                                              		      actionLink("addlink", "Sample Input Sequences"),	
                         			                              textarea(id="Sequence", rows="10", cols="100", style="float: none; width:100%"),
			                                                      fileInput('file1', 'Or upload your input file', accept=c('test/FASTA', 'FASTA', '.fasta', '.txt')),
                                                                              actionButton("submitbutton", "Submit", class = "btn btn-primary"),
                                                                              a(class='btn btn-default', href='', "Reset")
                                                                              ) #div
 
                                                               	      }) #withTags
                                                          ),  # First box

                                                box(title = h4(align="center", "Prediction results will be displayed here when the output is ready. You can also download the output as a CSV file as soon as the results are displayed below."), width = 11,
                                                            downloadButton("downloadData", 'Download output as a CSV file'),
                                                            verbatimTextOutput('contents'))
                                                           ) # Second box

                                           
							 ), #First TabItem

	        	# Second tab content
				     tabItem(tabName = "algorithm", includeMarkdown("markdown/overview.md")),
				     tabItem(tabName = "download", includeMarkdown("markdown/download.md")),
			  	     tabItem(tabName = "help", includeMarkdown("markdown/help.md")),
				     tabItem(tabName = "contact", includeMarkdown("markdown/contact.md"))
                	) # tabItems
                      ) #dashboardBody
      ) #dashboardPage



# Read in the model
mm <- readRDS("models/hybrid_xgb2.rds")

server <- shinyServer(function(input, output, session) {

  observe({
    FASTADATA <- ''
    fastaexample <- '>XR_007761781.1,X
CTTCTCCAGCGGCATGGCGAGGCGGGACGGGGCTTGAACCTTGAACTAATAATTACCTACGCGTGTTCACCAAGAAAGAG
GACCTGTTTTCGAGTGACATCTTGAGGTCCTACCAAGTTTAAAACGTGGAAAACCAGTGTTAATCATGCTGCGAACCAGA
AAACTGTCTCGGCTTTCTTTCCTGATGTCCCATGCATGGGCCTCCAACCACTGCTGGCTGCCGCTCTCCTACTGCAAGAT
GCATTTCATCAGAAAGAAGCATAAGTCATCCAAGGTTGTCAAATACATCCAGTACAAAAAAGAGGATATCAAATCTGAAC
AGAGTAAGAAAGAAGTTAACAAAATACCAAATAAGGTCTCTCTGGAGTTCCAAGGAAACAGTAAATCTTACAAAACATTG
AAAACAATACCTCCAGTTTCTCAAGAGAAAAGCTTTGCCAGAGTCCCCATTGAAACTGCAGCAAAGACAGAAATTGAGCA
GGAAGAGTGTGATGATGTCCCAAAGATACAAGTTGCCTCAGTGGGGGAGGTTCTGGAGAAGGTTGCAGAACTTTATAATG
AAAACAGGAGATTCTCTACTTTGCAATTTGCAGAAAGCAAAATGAAAAATCCAATAACAGGAAATCAGTACTCTGGAAAG
GGAAACAAAATTAGAATTATTCAGAGTCAGAAAAATGACCATGGTAAGATTCCTAATAATAATAAGAAAAAAATTGTTAA
ATTTAGTAATGCAAGAGCAATGAGTGGTGATACATCAGTAATGTGTCATGGTAATTCTGGTCCAGTTGATACTGTTTCAG
AAAATCATATGGTTAATCAAATAGGAAATGATCCAGAAATAAAAGAGGAAAATAATTCAGAGTTTAGCAGAACAACCACC
AGTACAGATGATCTTGATAAGCCTGTGTCAGAACAGGAATTCCACAACATTCGCAGAAATCCTCTTGGCATCCAAATGTT
GTCAGAAGGCATCTATAACCAATTGTTTAAGGAAGTAACCAATGAGGGCTGTGAAAAGGAGGATCTGGAAAGAACCAAAG
AACACTTACAGACTCATGGGTTGTGGAACAAAATTCCATCAGTCATACCCAGTGTGGAGTTTGAGTTGCCTGAAATAGAA
GGTAATGATTTAGACCAGCATTTTCGTGTCATAGCTGAGAATCAGTGTTCTCAGTACAAAGAGCTGCTCAACACGCTGGT
GGAAGATGTTCCAAGTCTCCCAAATAAGTGGGAGTTTGCTCCGGGCTGGACACGTTACCATCCTGATGGAAGTTGTTCAT
CCGTGGATTATCCAGAATGCAGTGCGATAGTTTTTGATGTGGAGGTTTGTGTTACGGAAGGCAACCAGCCTACCATGGCC
ACAGCAGTGTCAAACAAATACTGGTATTCATGGTGCTCTGAAGCACTCATAAATCCAGATGTTCACTTTGATGGGAGTGA
GATTAGAATGGAAGAACTGATTCCACTTGAAACTTCAAGCAGCATGAATCCAGTGCCCTCACACTCTGGGCGGGTTGTTG
TGGGCCATAATGTTAGTTATGATCGCCTGCGAGTCCGGGAGCAGTATTTGTTGAAGGAGACACCCCTTCGGTTTGTTGAT
ACAATGTCCCTCCACATTGCTGTTAGTGGACTAGTGTCAGAGCAGCGGGCCTTGCTGATGAAAAACAAAGGTGAAAAAAA
AATTAGGCTGCCATGGATGTCTGTTGGTTGCCAGAACAGCCTGGTTGAAGTGTATAAATTTTATTGCCGTCCTGAAAAGG
GTCTAGAAAAAAGCACAAGAGATGTGTTTGTTGATGGTACCTTGAGTGATGTAAGAGAAGACTTTCAGAATCTGATGAGT
TATTGTGCTACTGATGTGACAGCAACTCAAAAGGTTCTGGCTAAGTTGTTACCATTGTTTTATGAGCGGTTTCCTCATCC
AGTAACATTTTCTGGAATGCTAGAGATGGGCCTCACATTTCTCCCTGTTACAAAGAACTGGGAAAAGTACATTGAAGCTT
CAGAAGTTCAGTATCACCAAGTTGAGAGGCAGCTGAATGAGGAGCTTGTGAAACAAGTCCAGGCAAGTCTTGGCTCTATG
AAGAATAAAGAGTATGAGAATGATCCCTGGCTATGGAGCCTTGACTGGGCCCAACCCAAAGCTAGGGTGAAAAAACTGCC
TGGCTACCCAAATTGGTATAGAAAGCTCTGTGCACGTACGGGTGAAAGAGAAGGCACTCCAGAGCCTGAAAACATGAGTA
CCAGTTTGCAGATTGTTCCAAAAATTCTAAGATTGACTTGGAATGGATTTCCACTACACCATGAAAGGAAATTTGGATGG
GGTTATTTGAAACCTGTGTACCCCTCTTTCAAGGATATACCACAGTCTGAATGGGATTCCTATGCAGTCAATAACACTAG
TGAACCAGTATTTCCTGTGAAGGCTTTGTATGATATTTGTAATGAAAATGTCACAAGACAGTCAACAGCTCTACAAGATA
TTTCTCAGTATGATGAAAACTTGAATATACATGACCTAGAGTTAGGCAGCAAAAATAAAGGCCCTACAAAGAAAACAGGA
GGCACAACTAGTGAAGAGGGGAGCTTGAAGGATATAGGCATCCCAGGGGTGGGGTTCGTTCCACTGCCACACAAGGACGG
TGCAGGATGTCGGGTTGGCAACCCATTGGCCAAGGATTTCCTGGGAAAGATAGAAGATGGCACATTAACCAGTCACTTGG
GTGATGTTGCAAAGCTTGTGTTGGAGACCAGCAAGTCATTGTCCTATTGGAAGAACAACAGAGATCGCATCCTTTCTCAG
ATGGTTGTGTGGAATGATCATTCTGTCTTACCACATGAAGTCACATCACATGATGCCTACACGAGAGAGAACAAGTATGG
GATCATAGTGCCCATGGTTGTCTCCGCTGGGACCATCACCCGGCGGGCAGTGGAGCGGACCTGGATGACAGCTTCCAATG
CCTACGCGGACCGAATAGGCTCAGAGCTGAAGGCCATGGTTGAGGCTCCGCCGGGTTATAAATTTGTTGGTGCTGATGTA
GACTCCCAGGAACTGTGGATTGCATCATTATTGGGAGATGCATACTTCACTGGAGAGCACGGGGCCACTGCTCTTGGCTG
GATGACTCTCCAAGGGAAGAAGAGTGATGGAACAGACATGCACAGCCACACGGCTCGCAGTGCTGGAATTTCCCGAGACC
ATGCCAAGGTCATTAACTATGGGAGAATTTATGGTGCTGGACTGCGATTCATCCAGAGGCTCCTGAAGCAGTATAATCCC
AAACTGTCTGACCAAGAAACTAAAAGGAAAGCAGAACACTTGTTTGCAGTGACCAAGGGGGAGAAGGGCTGGTACCTGAA
TGAGGTTGGGGAGCAGCTTGCCCTAGATATAGGTCATCCTGTGCCTGAAGAACCCCTGCCTCGAAGGAAGATCACAAAGC
TCTTATGGCAAGCCCGTGAGAGCAACTATGGAGCTAGTTTTTCTGAGATAGTGGAAACCCCTCCTGTGTGGATTGGAGGA
TCTGAGTCTCACATGTTCAACTGCCTGGAGGCCATCGCTCGGTGTGAAGAGCCAAAAACTCCAGTGCTGGGTGCTCGTAT
GACTCGAGCCCTTGAGCCGTACTGGGTGGATAACCAGTTCATGACTAGCAGGGTGAACTGGGTGGTACAGAGCTCAGCGG
TGGACTACCTGCACATCATGCTGGTGTGCATGCGGTGGTTGTTCACCAAGTACGGCATCTCCGGTCGCTTCTGCATCAGC
ATTCATGATGAGGTACGATACTTAGTGGCGGAAGGAGACTGTCATAGAGCAGCGCTGGCCCTCCAGATCACCAACCTTCT
TACACGGGCATACTTTGCCAGCAGACTGGGGTTCAAAGACTTGCCGCAGTCTGTCGCCTTCTTCAGTGGTGTGGACATTG
ACCAGGTGCTAAGAAAGGAGCCCCACATGGACTGTGTTACTCCCTCCAACCCCCAGGGGTTGATGAAGGGGTATAGCATC
AAGCCAGGACTGACACTGGATATACATGAAATCATAAAGTCAACTAATGGGAAACTGACTAAAGCTGATGCTGGCTGCTG
TGAAGAGGCAGCTGAAAGTGAGGATTGACAAAGTGCGCTGTATGCAGAAGAATCAAAAGGAATCCAATTTTGCAAGACGT
GTAGTTTTACTTGCGAAGGCCAAAGTCAACGAAATAGTTACGTTTGAACAGCTATCCTAAATGCAACTTTATAAATGCGA
ATGATTTACTGTTTTTATACCTTTATTTTCACCAAAGATATGGCAACTACGACAAATAAAATGAAGCAAAATACTCACCA
ATGAGTGAGAGAAATGGAAGAAAAAACCCGCACACATGTGAATGCATCTTGAGATGCCAAATACCATCTATAAACTGCAT
TGCCAGTGACAGGAGCACTCAGTAGAGCGTGACGATTCCGCAGCACCCAGCGACAGTAGGGCTGGGCAGACACTGGTACT
GTTAGTACCAATACTGCAAAGTGGTAAGCAAAGGGTCAAGAGTGAGATGTGGGACAGCCCAAGATGCCATGCTGTGGACT
GAGTATCACCTTGTCACCAGGATAGCAC
>XR_007763215.1,P
CTTTAAACGAATATGTTGTGTTTCGATACAGTTTCGAAGCTTTTGTTGAGCGTTTCCTCGGGCGGTTAGTAGTCCAGTGT
GCTGTGGGCCGCGATTCGGGTCAGGGGGTTGTGCAGCAGCGGCGACTGGCGAGGGATGGTTTGTGCTGAGAAAAAGCTAA
GGCGAAGTTTTTATCTACATCTCCTCTAAGAGATAACTTGAAGCCCCACCAACGAGTCTCGAGTGTTAGATGAGGTACCC
TAGAAATTTTGGGCTCCACCAAGCGTGGAAGAACTTGAGAGAGAGACTGCCGAGTAAAGAGGCTGCTACACGGAGTACTG
GAGTCTGGCGGTGTGACCTCGATGAAACTTGCCACCGCCTTCACCTCAACGCCAGAAGAGGGGAAGTGCATTGACTGACC
CAACACCAGTATCTGTGATCAGTCTTTAATAGTTTCCCTAGGTACCCATTTACCCACCAGCTCAAAAGGGAGAATGAACA
GCTGGGTGAGCTGCATGTCAACTGCCCAAGCCAGGGATTTATACCCAGACCCACGATTCATAGTGAGGCATGCCAACCAA
TCCACCATGGAGGCACTTTACAGCAGTGGGAATCTGAAAATCCACAAAATCCAGTTTGGTAATGTAACAAATGAATGTTG
TGTTCACTGCTTATAATGTAACAAACGAATGTTGTGTTCACTGCTTAAGTATCAGTAATAAATTGACAGTTGAGAAAAAA
GTGATTACAGTGCATCATCAGATTGTACTGAATAAATGCATGATTTGCAAACATTTAGCATTATGAACTTTATTATAATA
GAGATATTGATGTGTGATGGTGACTATAGAATTTGATATATTTGTAATTCACTTTATGGCTGAGCAGTCTGGGCTCACAC
CTGTTATGCTTGCAGCCTTAGTCAATTTCAATAAGTCGTTTGATTCAGTGCATCACAAGGCACTCCGAAATATCTTGCGA
CTCCGCAGGAAGGACTATATTAGTTTGCTGACTGGCCTACATTCTGGGACAGAGTGCTGTGAAATGTGGGGGAAGCTTGT
CCAGCTTCTTTCCCATGATGCAGGAGCAAGACGGGGCTCAGTCCTTGCCTCGTCGCTTTTCAACACTTGTAGGGACTGGA
TACTTGGCAGAGTTGTGGACCAGATTCAGTGTAGAGCATCTGTTGGCAATACCAGGGTCTCTGACCTTGTTTTTGCTGAT
GATGCAGTAATCCTGGAGTAGTCGCTGAACGTTCTGGTGATGGCTCTTGAGGCACTGCACAAGGAGGCGAAGCCCTTGGG
CTTCAGGAAAAAGCTCTTGAGGCACTGCACAAGGAGGCGAAGCCCTTGGAACTTCAGGTTTCCTGGGCCAAGACCAAGGT
ACAGTTGTTTGGAGGCTCATTAGATGAAACTGTACAGTTTGTTCATGCGAGTGCCAAGGACAGTGAGATCTTAGAAAATA
TCACATACTTTGGTATCATAGTTCATACTGGCGAGTCTCTTCAGAAAGTCTTACTCTTTTGTTTCTTTTTTAGGAGCAGC
GGGCTTTTTTTTTATTGT
>XR_007763216.1,P
ACCTGTTCCAATTTGCATCGCTTTTTAGGTCCTCCTTGGTACATTTTTATACTTAGATGCTGCATTCAGTGATAGTAAAG
ATTGGGTACCGCGATGCAAAGTTCAGGATGGGATGTCTCACGACCAAGGATGGCAGAGGACTCCACACCGGGCTGCTTCA
AAATAGCCGCCACAGCAGAATGATTGGTGCAACCATAGAAGTTATTTGTTCTGTCAATCCCGTATTTCCGCCGCTCTTGC
CGGGGACACCTGGCGCGGCGGGCGCTCCAACAACACTTATCGTAGTGGAGTAATCTCAGCATGACTTTCGCCGTTGTCTG
GTTTGTGTCCAAACGAAGCATACCTACCTTGCAGTCAATGGTCCAGTCCAGACATAGCCAAGGAGAGTGATTGAGGAGTG
CACGGTGAAGGAATATTCTCGTCCTGAAGGAGTTTGTGTGGGAATGCAGAGAGCTGGGCTGATAGGAATATACTGCACAA
GATAGGTGGATGGCAGGCAGAAGGCAGATGCATGTGTGGTGCACTGTGAACCGTGCCAGGTGTTGGAACACTCCAACAGG
TAACACGTACTTGCCTGGATGGAGATGACTCCCACAGTGGAATCACACTTCAGAAGGTCTATCATAGAGAAGGGAAAATT
GTAGTCAACTTTGAAAAGAACATTTCACTTTCAATAACATATAAGAAATATGTCCA'
     if(input$addlink>0) {
      isolate({
        FASTADATA <- fastaexample
        updateTextInput(session, inputId = "Sequence", value = FASTADATA)
      })
    }
  })

  datasetInput <- reactive({

    inFile <- input$file1
    inTextbox <- input$Sequence

    if (is.null(inTextbox)) {
      return("Please insert/upload sequence in FASTA format")
    } else {
      if (is.null(inFile)) {
        input_fasta <- inTextbox
    
       tmpfile <- tempfile(fileext = ".fasta")
write.fasta(
  sequences = input_fasta,
  names = names(input_fasta),
  nbchar = 80,
  file.out = tmpfile
)

input_fasta <- readFASTA(tmpfile)

       # write.fasta(sequences = input_fasta, names = names(input_fasta),
        #            nbchar = 80, , file.out = "input.fasta")
	#input_fasta <- readFASTA("input.fasta") ################ Here replace input.fasta with your own input file


############################# Generate Features when input is from the text box ###################################
mm <- readRDS("models/hybrid_xgb2.rds")
#input_fasta <- readFASTA("input.fasta")
#input_fasta <- sapply(input_fasta, function(x) {
input_fasta <- readFASTA(tmpfile)

input_fasta <- sapply(input_fasta, function(x) { 
 seq <- paste(x, collapse = "")
  chartr("Tt", "Uu", seq)
})


####################################################k riboNucleotide Composition (kNUComposition_RNA)#############################
######################## k =1 nucleotide composition ########################
k1ncRNA <- kNUComposition_RNA(input_fasta, rng = 1, reverse = FALSE, upto = FALSE, normalized = TRUE, ORF = FALSE, reverseORF = TRUE, label = c())
k1ncRNA <- cbind(ID=rownames(k1ncRNA), k1ncRNA)
k1ncRNA <- as.data.frame(k1ncRNA)
###############################################################################################################################
######################## k =2 di-nucleotide composition ########################
k2ncRNA <- kNUComposition_RNA(input_fasta, rng = 2, reverse = FALSE, upto = FALSE, normalized = TRUE, ORF = FALSE, reverseORF = TRUE, label = c())
k2ncRNA <- cbind(ID=rownames(k2ncRNA), k2ncRNA)
k2ncRNA <- as.data.frame(k2ncRNA)
###############################################################################################################################
######################## k =3 tri-nucleotide composition ########################
k3ncRNA <- kNUComposition_RNA(input_fasta, rng = 3, reverse = FALSE, upto = FALSE, normalized = TRUE, ORF = FALSE, reverseORF = TRUE, label = c())
k3ncRNA <- cbind(ID=rownames(k3ncRNA), k3ncRNA)
k3ncRNA <- as.data.frame(k3ncRNA)
###############################################################################################################################

#####################Codon Usage in RNA (CodonUsage_RNA)#####################################################################
CodonU <- CodonUsage_RNA(input_fasta, ORF = FALSE, reverseORF = TRUE, label = c())
CodonU <- cbind(ID=rownames(CodonU), CodonU)
CodonU <- as.data.frame(CodonU)
colnames(CodonU) <- c("ID", "Codon1", "Codon2", "Codon3", "Codon4", "Codon5", "Codon6", "Codon7", "Codon8", "Codon9", "Codon10", "Codon11", "Codon12", "Codon13", "Codon14", "Codon15", "Codon16", "Codon17", "Codon18", "Codon19", "Codon20", "Codon21", "Codon22", "Codon23", "Codon24", "Codon25", "Codon26", "Codon27", "Codon28", "Codon29", "Codon30", "Codon31", "Codon32", "Codon33", "Codon34", "Codon35", "Codon36", "Codon37", "Codon38", "Codon39", "Codon40", "Codon41", "Codon42", "Codon43", "Codon44", "Codon45", "Codon46", "Codon47", "Codon48", "Codon49", "Codon50", "Codon51", "Codon52", "Codon53", "Codon54", "Codon55", "Codon56", "Codon57", "Codon58", "Codon59", "Codon60", "Codon61", "Codon62", "Codon63", "Codon64")
###############################################################################################################################


#####################G_C content in RNA (G_Ccontent_RNA)#####################################################################
GCCR <- G_Ccontent_RNA(input_fasta, ORF = FALSE, reverseORF = TRUE, normalized = TRUE,  label = c())
#GCCR <- cbind(ID=rownames(GCCR), GCCR)
GCCR <- as.data.frame(cbind(k1ncRNA$ID, GCCR))
names(GCCR)[names(GCCR) == 'V1'] <- 'ID'
GCCR <- as.data.frame(GCCR)
colnames(GCCR) <- c("ID", "GC")
GCCR_tmp <- (as.data.frame(GCCR[, 2:2]))
colnames(GCCR_tmp) <- c("GC.1")
###############################################################################################################################

#####################Maximum Open Reading Frame length in RNA (maxORFlength_RNA)#####################################################################
mORFLen <- maxORFlength_RNA(input_fasta, reverse = TRUE, normalized = TRUE, label = c())
#mORFLen <- cbind(ID=rownames(mORFLen), mORFLen)
mORFLen <- as.data.frame(cbind(k1ncRNA$ID, mORFLen))
names(mORFLen)[names(mORFLen) == 'V1'] <- 'ID'
mORFLen <- as.data.frame(mORFLen)
colnames(mORFLen) <- c("ID", "mORFLen")
mORFLen_tmp <- (as.data.frame(mORFLen[, 2:2]))
colnames(mORFLen_tmp) <- c("mORFLen")

###############################################################################################################################

#####################Z_curve_12bit_RNA (Zcurve12bit_RNA)#####################################################################
Zcurve <- Zcurve12bit_RNA(input_fasta, ORF = FALSE, reverseORF = TRUE, label = c())
Zcurve <- as.data.frame(cbind(k1ncRNA$ID, Zcurve))
names(Zcurve)[names(Zcurve) == 'V1'] <- 'ID'
#Zcurve <- cbind(ID=rownames(Zcurve), Zcurve)
Zcurve <- as.data.frame(Zcurve)
###############################################################################################################################

#####################Z_curve_12bit_RNA (Zcurve12bit_RNA)#####################################################################
asdcRNA <- ASDC_RNA(input_fasta, ORF = FALSE, reverseORF = TRUE, label = c())
asdcRNA <- cbind(ID=rownames(asdcRNA), asdcRNA)
asdcRNA <- as.data.frame(asdcRNA)
colnames(asdcRNA) <- c("ID", "AA.1", "AC.1", "AG.1", "AU.1", "CA.1", "CC.1", "CG.1", "CU.1", "GA.1", "GC.2", "GG.1", "GU.1", "UA.1", "UC.1", "UG.1", "UU.1")

###############################################################################################################################

################################Di riboNucleotide Autocorrelation-Autocovariance (AutoCorDiNUC_RNA)#############################
autoCoDiNRNA <- AutoCorDiNUC_RNA(input_fasta, selectedIdx = c("Rise (RNA)", "Roll (RNA)", "Shift (RNA)", "Slide (RNA)", "Tilt (RNA)", "Twist (RNA)"), maxlag = 3, threshold = 1, type = c("Moran", "Geary", "NormalizeMBorto", "AC", "CC", "ACC"),label = c())
autoCoDiNRNA <- as.data.frame(cbind(k1ncRNA$ID, autoCoDiNRNA))
names(autoCoDiNRNA)[names(autoCoDiNRNA) == 'V1'] <- 'ID'
#autoCoDiNRNA <- cbind(ID=rownames(autoCoDiNRNA), autoCoDiNRNA)
autoCoDiNRNA <- as.data.frame(autoCoDiNRNA)
colnames(autoCoDiNRNA) <- c("ID", "autoCo1", "autoCo2", "autoCo3", "autoCo4", "autoCo5", "autoCo6", "autoCo7", "autoCo8", "autoCo9", "autoCo10", "autoCo11", "autoCo12", "autoCo13", "autoCo14", "autoCo15", "autoCo16", "autoCo17", "autoCo18", "autoCo19", "autoCo20", "autoCo21", "autoCo22", "autoCo23", "autoCo24", "autoCo25", "autoCo26", "autoCo27", "autoCo28", "autoCo29", "autoCo30", "autoCo31", "autoCo32", "autoCo33", "autoCo34", "autoCo35", "autoCo36", "autoCo37", "autoCo38", "autoCo39", "autoCo40", "autoCo41", "autoCo42", "autoCo43", "autoCo44", "autoCo45", "autoCo46", "autoCo47", "autoCo48", "autoCo49", "autoCo50", "autoCo51", "autoCo52", "autoCo53", "autoCo54", "autoCo55", "autoCo56", "autoCo57", "autoCo58", "autoCo59", "autoCo60", "autoCo61", "autoCo62", "autoCo63", "autoCo64", "autoCo65", "autoCo66", "autoCo67", "autoCo68", "autoCo69", "autoCo70", "autoCo71", "autoCo72", "autoCo73", "autoCo74", "autoCo75", "autoCo76", "autoCo77", "autoCo78", "autoCo79", "autoCo80", "autoCo81", "autoCo82", "autoCo83", "autoCo84", "autoCo85", "autoCo86", "autoCo87", "autoCo88", "autoCo89", "autoCo90", "autoCo91", "autoCo92", "autoCo93", "autoCo94", "autoCo95", "autoCo96", "autoCo97", "autoCo98", "autoCo99", "autoCo100", "autoCo101", "autoCo102", "autoCo103", "autoCo104", "autoCo105", "autoCo106", "autoCo107", "autoCo108", "autoCo109", "autoCo110", "autoCo111", "autoCo112", "autoCo113", "autoCo114", "autoCo115", "autoCo116", "autoCo117", "autoCo118", "autoCo119", "autoCo120", "autoCo121", "autoCo122", "autoCo123", "autoCo124", "autoCo125", "autoCo126", "autoCo127", "autoCo128", "autoCo129", "autoCo130", "autoCo131", "autoCo132", "autoCo133", "autoCo134", "autoCo135", "autoCo136", "autoCo137", "autoCo138", "autoCo139", "autoCo140", "autoCo141", "autoCo142", "autoCo143", "autoCo144", "autoCo145", "autoCo146", "autoCo147", "autoCo148", "autoCo149", "autoCo150", "autoCo151", "autoCo152", "autoCo153", "autoCo154", "autoCo155", "autoCo156", "autoCo157", "autoCo158", "autoCo159", "autoCo160", "autoCo161", "autoCo162", "autoCo163", "autoCo164", "autoCo165", "autoCo166", "autoCo167", "autoCo168", "autoCo169", "autoCo170", "autoCo171", "autoCo172", "autoCo173", "autoCo174", "autoCo175", "autoCo176", "autoCo177", "autoCo178", "autoCo179", "autoCo180", "autoCo181", "autoCo182", "autoCo183", "autoCo184", "autoCo185", "autoCo186", "autoCo187", "autoCo188", "autoCo189", "autoCo190", "autoCo191", "autoCo192", "autoCo193", "autoCo194", "autoCo195", "autoCo196", "autoCo197", "autoCo198", "autoCo199", "autoCo200", "autoCo201", "autoCo202", "autoCo203", "autoCo204", "autoCo205", "autoCo206", "autoCo207", "autoCo208", "autoCo209", "autoCo210", "autoCo211", "autoCo212", "autoCo213", "autoCo214", "autoCo215", "autoCo216", "autoCo217", "autoCo218", "autoCo219", "autoCo220", "autoCo221", "autoCo222", "autoCo223", "autoCo224", "autoCo225", "autoCo226", "autoCo227", "autoCo228", "autoCo229", "autoCo230", "autoCo231", "autoCo232", "autoCo233", "autoCo234", "autoCo235", "autoCo236", "autoCo237", "autoCo238", "autoCo239", "autoCo240", "autoCo241", "autoCo242", "autoCo243", "autoCo244", "autoCo245", "autoCo246", "autoCo247", "autoCo248", "autoCo249", "autoCo250", "autoCo251", "autoCo252", "autoCo253", "autoCo254", "autoCo255", "autoCo256", "autoCo257", "autoCo258", "autoCo259", "autoCo260", "autoCo261", "autoCo262", "autoCo263", "autoCo264", "autoCo265", "autoCo266", "autoCo267", "autoCo268", "autoCo269", "autoCo270")
###############################################################################################################################

#####################Amphiphilic Pseudo-k riboNucleotide Composition-di(series) (APkNUCdi_RNA)###############################
APkNUCdiRNA <- APkNUCdi_RNA(input_fasta, selectedIdx = c("Rise (RNA)", "Roll (RNA)", "Shift (RNA)", "Slide (RNA)", "Tilt (RNA)", "Twist (RNA)"), lambda = 3, w = 0.05, l = 2, ORF = FALSE, reverseORF = TRUE, threshold = 1, label = c())
APkNUCdiRNA <- cbind(ID=rownames(APkNUCdiRNA), APkNUCdiRNA)
APkNUCdiRNA <- as.data.frame(APkNUCdiRNA)
colnames(APkNUCdiRNA) <- c("ID", "APk1", "APk2", "APk3", "APk4", "APk5", "APk6", "APk7", "APk8", "APk9", "APk10", "APk11", "APk12", "APk13", "APk14", "APk15", "APk16", "APk17", "APk18", "APk19", "APk20", "APk21", "APk22", "APk23", "APk24", "APk25", "APk26", "APk27", "APk28", "APk29", "APk30", "APk31", "APk32", "APk33", "APk34")
#############################################################################################################################

################ Get hybrids
#hybrid <- as.data.frame(cbind(k1ncRNA, k2ncRNA, k3ncRNA, CodonU, GCCR_tmp, mORFLen_tmp, Zcurve, asdcRNA, autoCoDiNRNA, APkNUCdiRNA))

hybrid <- as.data.frame(cbind(k1ncRNA, k2ncRNA[2:17], k3ncRNA[2:65], CodonU[2:65], GCCR_tmp, mORFLen_tmp, Zcurve[2:13], asdcRNA[2:17], autoCoDiNRNA[2:271], APkNUCdiRNA[2:35]))
hybrid[ , -1] <- lapply(hybrid[ , -1], as.numeric)

Layer1BC_out <- predict(mm, hybrid, type="prob")
Layer1BC_out2 <- predict(mm, hybrid)
Layer1BCpred <- cbind(Layer1BC_out, Layer1BC_out2)
Layer1BCpred <- cbind(ID=rownames(hybrid), Layer1BCpred)
rownames(Layer1BCpred) <- NULL
colnames(Layer1BCpred) <- c("ID", "LncRNA", "mRNA", "Prediction")
is.num <- sapply(Layer1BCpred, is.numeric)
Layer1BCpred[is.num] <- lapply(Layer1BCpred[is.num], round, 3)
Layer1BCpred_Subset <- Layer1BCpred[grep("P", Layer1BCpred$Prediction), ] ###### Get Ids of sequences that were predicted as LncRNAs.
validate(need(nrow(Layer1BCpred_Subset)!=0, "There are potentially NO LncRNA sequences in the dataset."))
Layer1BCpred$Prediction <- ifelse(Layer1BCpred$Prediction == "P", "LncRNA",
                           ifelse(Layer1BCpred$Prediction == "X", "mRNA",
                                  Layer1BCpred$Prediction))

print(Layer1BCpred)
#write.csv(Layer1BCpred, file='prediction_output.csv', row.names=FALSE)
  }
      else {

 input_fasta <- readFASTA(inFile$datapath) ################ Here replace input.fasta with your own input file
input_fasta <- sapply(input_fasta, function(x) {
  seq <- paste(x, collapse = "")
  chartr("Tt", "Uu", seq)
})

#input_fasta <- input_fasta[(sapply(input_fasta, protcheck))] ########## Check if there are non-standard amino acids in the sequences

############################# Generate Features when input is uploaded as a file ###################################
####################################################k riboNucleotide Composition (kNUComposition_RNA)#############################
######################## k =1 nucleotide composition ########################
k1ncRNA <- kNUComposition_RNA(input_fasta, rng = 1, reverse = FALSE, upto = FALSE, normalized = TRUE, ORF = FALSE, reverseORF = TRUE, label = c())
k1ncRNA <- cbind(ID=rownames(k1ncRNA), k1ncRNA)
k1ncRNA <- as.data.frame(k1ncRNA)
###############################################################################################################################
######################## k =2 di-nucleotide composition ########################
k2ncRNA <- kNUComposition_RNA(input_fasta, rng = 2, reverse = FALSE, upto = FALSE, normalized = TRUE, ORF = FALSE, reverseORF = TRUE, label = c())
k2ncRNA <- cbind(ID=rownames(k2ncRNA), k2ncRNA)
k2ncRNA <- as.data.frame(k2ncRNA)
###############################################################################################################################
######################## k =3 tri-nucleotide composition ########################
k3ncRNA <- kNUComposition_RNA(input_fasta, rng = 3, reverse = FALSE, upto = FALSE, normalized = TRUE, ORF = FALSE, reverseORF = TRUE, label = c())
k3ncRNA <- cbind(ID=rownames(k3ncRNA), k3ncRNA)
k3ncRNA <- as.data.frame(k3ncRNA)
###############################################################################################################################

#####################Codon Usage in RNA (CodonUsage_RNA)#####################################################################
CodonU <- CodonUsage_RNA(input_fasta, ORF = FALSE, reverseORF = TRUE, label = c())
CodonU <- cbind(ID=rownames(CodonU), CodonU)
CodonU <- as.data.frame(CodonU)
colnames(CodonU) <- c("ID", "Codon1", "Codon2", "Codon3", "Codon4", "Codon5", "Codon6", "Codon7", "Codon8", "Codon9", "Codon10", "Codon11", "Codon12", "Codon13", "Codon14", "Codon15", "Codon16", "Codon17", "Codon18", "Codon19", "Codon20", "Codon21", "Codon22", "Codon23", "Codon24", "Codon25", "Codon26", "Codon27", "Codon28", "Codon29", "Codon30", "Codon31", "Codon32", "Codon33", "Codon34", "Codon35", "Codon36", "Codon37", "Codon38", "Codon39", "Codon40", "Codon41", "Codon42", "Codon43", "Codon44", "Codon45", "Codon46", "Codon47", "Codon48", "Codon49", "Codon50", "Codon51", "Codon52", "Codon53", "Codon54", "Codon55", "Codon56", "Codon57", "Codon58", "Codon59", "Codon60", "Codon61", "Codon62", "Codon63", "Codon64")
###############################################################################################################################


#####################G_C content in RNA (G_Ccontent_RNA)#####################################################################
GCCR <- G_Ccontent_RNA(input_fasta, ORF = FALSE, reverseORF = TRUE, normalized = TRUE,  label = c())
#GCCR <- cbind(ID=rownames(GCCR), GCCR)
GCCR <- as.data.frame(cbind(k1ncRNA$ID, GCCR))
names(GCCR)[names(GCCR) == 'V1'] <- 'ID'
GCCR <- as.data.frame(GCCR)
colnames(GCCR) <- c("ID", "GC")
GCCR_tmp <- (as.data.frame(GCCR[, 2:2]))
colnames(GCCR_tmp) <- c("GC.1")
###############################################################################################################################

#####################Maximum Open Reading Frame length in RNA (maxORFlength_RNA)#####################################################################
mORFLen <- maxORFlength_RNA(input_fasta, reverse = TRUE, normalized = TRUE, label = c())
#mORFLen <- cbind(ID=rownames(mORFLen), mORFLen)
mORFLen <- as.data.frame(cbind(k1ncRNA$ID, mORFLen))
names(mORFLen)[names(mORFLen) == 'V1'] <- 'ID'
mORFLen <- as.data.frame(mORFLen)
colnames(mORFLen) <- c("ID", "mORFLen")
mORFLen_tmp <- (as.data.frame(mORFLen[, 2:2]))
colnames(mORFLen_tmp) <- c("mORFLen")

###############################################################################################################################

#####################Z_curve_12bit_RNA (Zcurve12bit_RNA)#####################################################################
Zcurve <- Zcurve12bit_RNA(input_fasta, ORF = FALSE, reverseORF = TRUE, label = c())
Zcurve <- as.data.frame(cbind(k1ncRNA$ID, Zcurve))
names(Zcurve)[names(Zcurve) == 'V1'] <- 'ID'
#Zcurve <- cbind(ID=rownames(Zcurve), Zcurve)
Zcurve <- as.data.frame(Zcurve)
###############################################################################################################################

#####################Z_curve_12bit_RNA (Zcurve12bit_RNA)#####################################################################
asdcRNA <- ASDC_RNA(input_fasta, ORF = FALSE, reverseORF = TRUE, label = c())
asdcRNA <- cbind(ID=rownames(asdcRNA), asdcRNA)
asdcRNA <- as.data.frame(asdcRNA)
colnames(asdcRNA) <- c("ID", "AA.1", "AC.1", "AG.1", "AU.1", "CA.1", "CC.1", "CG.1", "CU.1", "GA.1", "GC.2", "GG.1", "GU.1", "UA.1", "UC.1", "UG.1", "UU.1")

###############################################################################################################################

################################Di riboNucleotide Autocorrelation-Autocovariance (AutoCorDiNUC_RNA)#############################
autoCoDiNRNA <- AutoCorDiNUC_RNA(input_fasta, selectedIdx = c("Rise (RNA)", "Roll (RNA)", "Shift (RNA)", "Slide (RNA)", "Tilt (RNA)", "Twist (RNA)"), maxlag = 3, threshold = 1, type = c("Moran", "Geary", "NormalizeMBorto", "AC", "CC", "ACC"),label = c())
autoCoDiNRNA <- as.data.frame(cbind(k1ncRNA$ID, autoCoDiNRNA))
names(autoCoDiNRNA)[names(autoCoDiNRNA) == 'V1'] <- 'ID'
#autoCoDiNRNA <- cbind(ID=rownames(autoCoDiNRNA), autoCoDiNRNA)
autoCoDiNRNA <- as.data.frame(autoCoDiNRNA)
colnames(autoCoDiNRNA) <- c("ID", "autoCo1", "autoCo2", "autoCo3", "autoCo4", "autoCo5", "autoCo6", "autoCo7", "autoCo8", "autoCo9", "autoCo10", "autoCo11", "autoCo12", "autoCo13", "autoCo14", "autoCo15", "autoCo16", "autoCo17", "autoCo18", "autoCo19", "autoCo20", "autoCo21", "autoCo22", "autoCo23", "autoCo24", "autoCo25", "autoCo26", "autoCo27", "autoCo28", "autoCo29", "autoCo30", "autoCo31", "autoCo32", "autoCo33", "autoCo34", "autoCo35", "autoCo36", "autoCo37", "autoCo38", "autoCo39", "autoCo40", "autoCo41", "autoCo42", "autoCo43", "autoCo44", "autoCo45", "autoCo46", "autoCo47", "autoCo48", "autoCo49", "autoCo50", "autoCo51", "autoCo52", "autoCo53", "autoCo54", "autoCo55", "autoCo56", "autoCo57", "autoCo58", "autoCo59", "autoCo60", "autoCo61", "autoCo62", "autoCo63", "autoCo64", "autoCo65", "autoCo66", "autoCo67", "autoCo68", "autoCo69", "autoCo70", "autoCo71", "autoCo72", "autoCo73", "autoCo74", "autoCo75", "autoCo76", "autoCo77", "autoCo78", "autoCo79", "autoCo80", "autoCo81", "autoCo82", "autoCo83", "autoCo84", "autoCo85", "autoCo86", "autoCo87", "autoCo88", "autoCo89", "autoCo90", "autoCo91", "autoCo92", "autoCo93", "autoCo94", "autoCo95", "autoCo96", "autoCo97", "autoCo98", "autoCo99", "autoCo100", "autoCo101", "autoCo102", "autoCo103", "autoCo104", "autoCo105", "autoCo106", "autoCo107", "autoCo108", "autoCo109", "autoCo110", "autoCo111", "autoCo112", "autoCo113", "autoCo114", "autoCo115", "autoCo116", "autoCo117", "autoCo118", "autoCo119", "autoCo120", "autoCo121", "autoCo122", "autoCo123", "autoCo124", "autoCo125", "autoCo126", "autoCo127", "autoCo128", "autoCo129", "autoCo130", "autoCo131", "autoCo132", "autoCo133", "autoCo134", "autoCo135", "autoCo136", "autoCo137", "autoCo138", "autoCo139", "autoCo140", "autoCo141", "autoCo142", "autoCo143", "autoCo144", "autoCo145", "autoCo146", "autoCo147", "autoCo148", "autoCo149", "autoCo150", "autoCo151", "autoCo152", "autoCo153", "autoCo154", "autoCo155", "autoCo156", "autoCo157", "autoCo158", "autoCo159", "autoCo160", "autoCo161", "autoCo162", "autoCo163", "autoCo164", "autoCo165", "autoCo166", "autoCo167", "autoCo168", "autoCo169", "autoCo170", "autoCo171", "autoCo172", "autoCo173", "autoCo174", "autoCo175", "autoCo176", "autoCo177", "autoCo178", "autoCo179", "autoCo180", "autoCo181", "autoCo182", "autoCo183", "autoCo184", "autoCo185", "autoCo186", "autoCo187", "autoCo188", "autoCo189", "autoCo190", "autoCo191", "autoCo192", "autoCo193", "autoCo194", "autoCo195", "autoCo196", "autoCo197", "autoCo198", "autoCo199", "autoCo200", "autoCo201", "autoCo202", "autoCo203", "autoCo204", "autoCo205", "autoCo206", "autoCo207", "autoCo208", "autoCo209", "autoCo210", "autoCo211", "autoCo212", "autoCo213", "autoCo214", "autoCo215", "autoCo216", "autoCo217", "autoCo218", "autoCo219", "autoCo220", "autoCo221", "autoCo222", "autoCo223", "autoCo224", "autoCo225", "autoCo226", "autoCo227", "autoCo228", "autoCo229", "autoCo230", "autoCo231", "autoCo232", "autoCo233", "autoCo234", "autoCo235", "autoCo236", "autoCo237", "autoCo238", "autoCo239", "autoCo240", "autoCo241", "autoCo242", "autoCo243", "autoCo244", "autoCo245", "autoCo246", "autoCo247", "autoCo248", "autoCo249", "autoCo250", "autoCo251", "autoCo252", "autoCo253", "autoCo254", "autoCo255", "autoCo256", "autoCo257", "autoCo258", "autoCo259", "autoCo260", "autoCo261", "autoCo262", "autoCo263", "autoCo264", "autoCo265", "autoCo266", "autoCo267", "autoCo268", "autoCo269", "autoCo270")
###############################################################################################################################

#####################Amphiphilic Pseudo-k riboNucleotide Composition-di(series) (APkNUCdi_RNA)###############################
APkNUCdiRNA <- APkNUCdi_RNA(input_fasta, selectedIdx = c("Rise (RNA)", "Roll (RNA)", "Shift (RNA)", "Slide (RNA)", "Tilt (RNA)", "Twist (RNA)"), lambda = 3, w = 0.05, l = 2, ORF = FALSE, reverseORF = TRUE, threshold = 1, label = c())
APkNUCdiRNA <- cbind(ID=rownames(APkNUCdiRNA), APkNUCdiRNA)
APkNUCdiRNA <- as.data.frame(APkNUCdiRNA)
colnames(APkNUCdiRNA) <- c("ID", "APk1", "APk2", "APk3", "APk4", "APk5", "APk6", "APk7", "APk8", "APk9", "APk10", "APk11", "APk12", "APk13", "APk14", "APk15", "APk16", "APk17", "APk18", "APk19", "APk20", "APk21", "APk22", "APk23", "APk24", "APk25", "APk26", "APk27", "APk28", "APk29", "APk30", "APk31", "APk32", "APk33", "APk34")
#############################################################################################################################

################ Get hybrids
#hybrid <- as.data.frame(cbind(k1ncRNA, k2ncRNA, k3ncRNA, CodonU, GCCR_tmp, mORFLen_tmp, Zcurve, asdcRNA, autoCoDiNRNA, APkNUCdiRNA))

hybrid <- as.data.frame(cbind(k1ncRNA, k2ncRNA[2:17], k3ncRNA[2:65], CodonU[2:65], GCCR_tmp, mORFLen_tmp, Zcurve[2:13], asdcRNA[2:17], autoCoDiNRNA[2:271], APkNUCdiRNA[2:35]))
hybrid[ , -1] <- lapply(hybrid[ , -1], as.numeric)

Layer1BC_out <- predict(mm, hybrid, type="prob")
Layer1BC_out2 <- predict(mm, hybrid)
Layer1BCpred <- cbind(Layer1BC_out, Layer1BC_out2)
Layer1BCpred <- cbind(ID=rownames(hybrid), Layer1BCpred)
rownames(Layer1BCpred) <- NULL
colnames(Layer1BCpred) <- c("ID", "LncRNA", "mRNA", "Prediction")
is.num <- sapply(Layer1BCpred, is.numeric)
Layer1BCpred[is.num] <- lapply(Layer1BCpred[is.num], round, 3)
Layer1BCpred_Subset <- Layer1BCpred[grep("P", Layer1BCpred$Prediction), ] ###### Get Ids of sequences that were predicted as LncRNAs.
validate(need(nrow(Layer1BCpred_Subset)!=0, "There are potentially NO LncRNA sequences in the dataset."))
Layer1BCpred$Prediction <- ifelse(Layer1BCpred$Prediction == "P", "LncRNA",
                           ifelse(Layer1BCpred$Prediction == "X", "mRNA",
                                  Layer1BCpred$Prediction))

print(Layer1BCpred)

#write.csv(LayerBCpred, file='prediction_output.csv', row.names=FALSE)

      }
    }
  })



output$contents <- renderPrint({
    if (input$submitbutton>0) { 
      isolate(datasetInput()) 
    } else {
        return(cat(""))
      }
  })


  output$downloadData <- downloadHandler(
  #  filename = function() { paste('prediction_results', '.csv', sep='') },
filename = function(){"prediction_results.csv"}, 
#filename = function(){Sys.Date(), ".csv"},

 content = function(file) {
      write.csv(datasetInput(), file, row.names=FALSE)
   })

  })

shinyApp(ui, server)
