#library(DT)
Sys.setenv(R_ZIPCMD="/usr/bin/zip")
library(shiny)
#library(shinyIncubator)
library(stringr)
library(RMySQL)
library(knitr)
library(xtable)
library(log4r)
library(shinyBS)
library(digest)
library(curl)
#library(dplyr)
library(stringi)
library(rentrez)

###>>> Additional functions
source("/tmp/additionalFunctions.R")
#source("/Users/schoi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/ICR_2020/cavada/additionalFunctions.R")
###

loggerServer <- create.logger()
logfile(loggerServer) <- 'logs/server.log'
level(loggerServer) <- 'INFO'

options(stringsAsFactors=F)
#host<-'localhost'
host<-'127.0.0.1'
version<-'v15_20140605'

disableActionButton <- function(id,session) {
  session$sendCustomMessage(type="jsCode",
                            list(code=paste0("$('#",id,"').prop('disabled',true)")) )
}
###>>> disable and change color of login button
disableActionButton2 <- function(id,session) {
  session$sendCustomMessage(type='jsCode', list(code=paste0("$('#",id,"').css('background-color','#e80a89').css('color','white').prop('disabled',true)")))
}
###

enableActionButton <- function(id,session) {
  session$sendCustomMessage(type="jsCode",
                            list(code=paste0("$('#",id,"').prop('disabled',false)")))
}

sendAlertMessage <- function(mess,session) {
  session$sendCustomMessage(type="jsCode",
                            list(code=paste0("alert('",mess,"')")))
}

# initShinyVars <- function(session) {
#   session$sendCustomMessage(type="jsCode",list(code="Shiny.unbindAll();
# Shiny.shinyapp.$inputValues.exploreTablePage = undefined;
# Shiny.shinyapp.$inputValues.exploreTableVariation = undefined;
# Shiny.bindAll();")
#   )
# }

shinyServer(function(input, output, session) {

  values <- reactiveValues()

  values$batch <- NULL
  values$sessionid <- paste(format(Sys.time(), "%d%b%Y_%H%M%S"),sep='_')
  values$iname <- ''
  values$inputFileType <- 'simple'
  # can be 'simple' (sample_id gene variation),
  #        'clinical' (Worksheet  InvID	DNA No	Gene	Nomenclature	Het-hom	Class	Exon)
  #     or 'VCF'- to be implemented
  values$sessionschanged <- 0
  values$filesRead <- NULL
  values$filesLoaded <- NULL
  values$lastExploreButton <- 0
  values$lastBatchButton <- 0
  values$lastParseButton <- 0
  values$lastLoadButton <- 0
  values$lastSaveButton <- 0
  values$lastReturnButton <- 0
  values$lastMarkButton <- 0
  values$lastAddNoteButton <- 0
  values$lastSaveNotesButton <- 0
  values$lastExcludeDeep <- TRUE
  values$lastExcludeExtra <- TRUE
  values$lastExcludeCommon <- FALSE
  values$lastExcludeIndels <- FALSE
  values$selectedRow <- NULL
  values$contentType <- 'front'
  values$exploreTable <- NULL
  values$selectedVar <- NULL
  values$selectedVarBaseline <- NULL
  values$selectedVarClinical <- NULL
  values$selectedVarClinicalCurr <- NULL
  values$exploreTableiDL <- 10
  values$exploreTableiDS <- 0
  values$showSaveNotes <- FALSE
  values$showNewNote <- FALSE
  values$newNoteCurr <- ''
  values$offsetLimit <- 10
  values$variantReport <- NULL

  ###>>> additional values ####
  values$showNewNoteExp <- TRUE
  values$showSaveNotesExp <- TRUE
  values$newNoteCurr2 <- ''
  #values$selectedVarClinicalCurr2 <- NULL
  values$lastloginButton <- 0
  #values$variationSel <- ""
  output$USERv2 <- output$USERv <- renderText("Anonymous user")
  values$USER <- ""
  values$lastsaveNotesButton <- 0
  values$logged <- FALSE
  values$curator <- "unknown"
  values$lastsignButton <- 0
  values$lastcreateAccountButton <- 0
  values$addFileName <- ""
  values$dirNotes <- ""
  values$notesupdate <- ""
  values$initVar <- ""
  output$updateAccMessage <-renderText(" ")
  output$registerMessage <-renderText(" ")
  output$registerError <- renderText("Please check the following errors: ")
  output$visibleComments <- renderText("The comments are only visible for the registered users")
  ###

  values$sources <- c('ALAMUT','EVS','TGP','TGP_PHASE3','ICR','HGMD','IARC','DMUDB','BIC','LOVD',
                      'EASTON','LINDOR','HOUDAYER','GUIDUGLI','WALKER','RMH','MUTTASTER',
                      'POLYPHEN2','DROST','UMD','CADD', 'SUSPECT','FINDLAY', 'GUIDUGLI2018',
                      'WESSEX','REVEL', 'GAVIN','LOVD3','TP53_SGE','TP53_Kato','TP53_Fut','BOUWMAN',
                      'BOUWMAN_2020','PARSONS_2019')
  # 'GUIDUGLI2018','WESSEX','REVEL','GAVIN','TP53_SGE','LOVD3')

  disableActionButton('parseButton',session)
  disableActionButton('loadButton',session)
  fixempty <- function(x) {return(ifelse(is.null(x)||is.na(x)||x==0,'',x))}
  freqprint <- function(x){return(sprintf('%5.4f',as.numeric(x)))}
  mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
  rs<-dbSendQuery(mydb,"select table_name,field_name,field_label,variant_report,variant_report_order,analysis_output,analysis_output_order from fields where online_report='Y'")
  f<-fetch(rs,-1)
  fn<-split(f$field_name,f$table_name)
  fl<-split(f$field_label,f$table_name)
  # fa is a structure with the fields for the analysis output
  # the fields from the main classification and dogma_batchlog tables have to be entered manually
  fa<-f[which(f$analysis_output=='Y' & !(f$table_name %in% c('main','classification','dogma_batchlog'))),]
  rownames(fa)<-paste(gsub('tgp_phase3_pop','tgp_phase3',gsub('tgp_pop','tgp',fa[,1])),fa[,2],sep='_')
  fv<-f[which(f$variant_report=='Y' & !(f$table_name %in% c('classification','dogma_batchlog'))),]
  rownames(fv)<-paste(gsub('tgp_phase3_pop','tgp_phase3',gsub('tgp_pop','tgp',fv[,1])),fv[,2],sep='_')
  fan<-split(fa$field_name,fa$table_name)
  fvn<-split(fv$field_name,fv$table_name)

  # fields is a structure for the display panels
  # create list of tables with named lists of fields descriptions:
  fields<-setNames(lapply(names(fl),function(x){return(setNames(paste0(fl[[x]],':'),fn[[x]]))}),names(fl))
  # to access the description fieldX from tableY, use: fields[['tableX']][['fieldX']]

  #hack for tgp_pop & tgp_phase3_pop - these secondary tables store the information in multiple rows, one for each population
  fields[["tgp_pop"]]<-gsub(' .EUR.','',fields[["tgp_pop"]][1:4])
  fields[["tgp_pop"]]<-setNames(fields[["tgp_pop"]],gsub('_EUR','',names(fields[["tgp_pop"]])))
  fields[["tgp_phase3_pop"]]<-gsub(' .ACB.','',fields[["tgp_phase3_pop"]][1:4])
  fields[["tgp_phase3_pop"]]<-setNames(fields[["tgp_phase3_pop"]],gsub('_ACB','',names(fields[["tgp_phase3_pop"]])))

  rs<-dbSendQuery(mydb,"select distinct gene,ensembl transcript, substring_index(refseq,'.',1) rtranscript from cappagenes")
  #                select distinct gene,transcript,substring_index(rtranscript,'.',1) rtranscript from main")
  genelist<-fetch(rs,-1)

  ###>>> remove noncancer genes ###
  noncancergenes  <- c('ACTC1', 'ACTA2', 'APOB','ATP7B',"CACNA1S","COL3A1","DSC2","DSG2","DSP","FBN1","GLA","KCNH2","KCNQ1","LDLR","LMNA","MYBPC3","MYH11","MYH7","MYL2","MYL3","OTC", "PCSK9","PKP2","PRKAG2","RYR1","RYR2","SCN5A","SMAD3", "TGFBR1", "TGFBR2","TMEM43", "TNNI3","TNNT2", "TPM1")
  genelist <- genelist[which(!genelist$gene %in% noncancergenes),]
  ###

  gl<-as.vector(apply(genelist,1,function(x){paste(x[1]," (",x[2],")",sep="",collapse="")}))

  rs<-dbSendQuery(mydb,"select * from tgp_populations")
  dataset<-fetch(rs,-1)
  tp<-as.list(apply(dataset,1,function(x){ifelse(x[3]=='',x[2],paste(x[2]," (Super-population=",x[3],")",sep="",collapse=""))}))
  tp<-setNames(tp,dataset[,1])

  dbDisconnect(mydb)

  observe ({
    if (!is.null(input$offsetLimit) && input$offsetLimit>0)
      values$offsetLimit<-input$offsetLimit
  })

  # explore stuff
  currentgenetrans <- ''
  # Initialize variationText every time we change the gene
  observe({
    if (!is.null(input$genetrans) && (input$genetrans!=currentgenetrans || input$genetrans=='')) {
      updateTextInput(session, 'variationText', value='')

      ###>>> if gene and variant change then the curated class will be reseted.
      #      input$clinicalCIGMAexp <- ""
      #      values$showNewNote <- FALSE
      ###

      currentgenetrans <- input$genetrans
    }
  })

  output$varCountText<-renderText({
    paste(varCount(),'variant(s)')
  })

  # Update the variation select list based on the selected gene
  # and the initial input in the variation text box
  varCount <- reactive({
    splitgene<-unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))
    gene<-splitgene[1]
    transcript<-splitgene[2]
    variation<-input$variationText
    shortp <- substr(variation,4,4)>='0' && substr(variation,4,4)<='9'
    v<-NULL
    if (nchar(gene)>0 && nchar(variation)>3 && grepl('^[cp][.]',variation,perl=T)) {
      # don't bother unless there is a gene selected and at least 4 characters in variation
      mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
      if (substr(variation,1,2) == 'c.') {
        rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                   "select CONCAT(hgvs_cdna,' (',hgvs_prot,')') from main_stable m where m.hgvs_cdna NOT REGEXP 'del|ins|dup' and m.gene='",gene,
                                   "' and m.transcript='",transcript,"' and m.hgvs_cdna like '",variation,"%'",
                                   " order by cdna_pos,offset"))
      }else{ # it has to be 'p.'
        if (shortp) {# short p. notation
          rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                     "select CONCAT(hgvs_cdna,' (',hgvs_prot_code1,')') from main_stable m where m.hgvs_cdna NOT REGEXP 'del|ins|dup' m.gene='",gene,
                                     "' and m.transcript='",transcript,"' and m.hgvs_prot_code1 like '",variation,"%'",
                                     " order by cdna_pos,offset"))
        }else{ # long p notation
          rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                     "select CONCAT(hgvs_cdna,' (',hgvs_prot,')') from main_stable m where m.hgvs_cdna NOT REGEXP 'del|ins|dup' and m.gene='",gene,
                                     "' and m.transcript='",transcript,"' and m.hgvs_prot like '",variation,"%'",
                                     " order by cdna_pos,offset"))
        }
      }
      if(!(grepl('^p[.]',variation,perl=T) && !shortp && nchar(variation)<6)) {
        # don't do long p queries unless at least 6 char long
        v<-fetch(rs,-1)
        dbDisconnect(mydb)
        # Change values for input$variationSel if there are any results
        if(!is.null(v) && class(v)=='data.frame' && dim(v)[1] ){
          updateSelectInput(session, "variationSel", choices = as.vector(v[,1]))
          return(dim(v)[1])
        }else{
          updateSelectInput(session, "variationSel", choices = c(''))
          return(0)
        }
      }else{
        v<-fetch(rs,-1)
        dbDisconnect(mydb)
        updateSelectInput(session, "variationSel", choices = c(''))
        return(0)
      }
    }else{
      updateSelectInput(session, "variationSel", choices = c(''))
      return(0)
    }
    dbDisconnect(mydb)
  })

  # Reactive function collecting loaded sets
  sesslist <- reactive({
    values$sessionschanged
    mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
    rs<-dbSendQuery(mydb,"select investigator_name,session_id,count(*) variations from dogma_batch group by investigator_name,session_id")
    sl<-fetch(rs,-1)
    dbDisconnect(mydb)
    return(as.vector(apply(sl,1,function(x){paste(x[1]," (",x[2],", ",x[3]," variants)",sep="",collapse="")})))
  })

  # Dynamic UI with the loaded sets
  output$loadedSets<-renderUI({
    sl<-sesslist()
    info(loggerServer,paste("Loaded sets='",sl,"'"))
    if (values$sessionschanged == 0) {
      cs<-''
    }else{
      cs<-sl[grep(isolate(values$sessionid),sl)]
    }
    info(loggerServer,paste("Selected set='",cs,"'",sep=""))
    tagList(tags$script(src = "select2-master/select2.js"),
            tags$link(href="select2-master/select2.css",rel="stylesheet"),
            tags$script(src = "js/include_select2.js"),
            do.call(selectInput,list(inputId = "batchid", label="Previously loaded set :", choices = c('',sl), selected=cs, selectize=FALSE))
    )
  })

  ###########
  # BUTTONS #
  ###########

  # Test header of clinical lab file
  fileIsClinical <- function(x) {
    if (grepl('worksheet',x[1],ignore.case=T,perl=T) &&
        grepl('invid',x[2],ignore.case=T,perl=T) &&
        grepl('dna.no',x[3],ignore.case=T,perl=T) &&
        grepl('gene',x[4],ignore.case=T,perl=T) &&
        grepl('nomenclature',x[5],ignore.case=T,perl=T) &&
        grepl('het.hom',x[6],ignore.case=T,perl=T) &&
        grepl('class',x[7],ignore.case=T,perl=T) &&
        grepl('exon',x[8],ignore.case=T,perl=T)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }

  # Upload reactive function
  observe({
    # input$uploadFile will be NA initially. After the user selects and uploads a
    # file, it will be a data frame with 'name', 'size', 'type', and 'datapath'
    # columns. The 'datapath' column will contain the local filenames where the
    # data can be found.
    inFile <- input$uploadFile
    error <- F
    if (!is.null(inFile) && !is.na(inFile) && (is.null(values$fileRead) || is.na(values$fileRead[inFile$datapath])))  {
      info(loggerServer,paste("Reading file:'",inFile$datapath, "'",sep=""))
      values$inputFileType <- 'simple'
      firstline<-scan(inFile$datapath, what=character(), sep="\t", quote='', nlines = 1)
      if (length(which(grepl('^gene$',firstline,ignore.case=T,perl=T)))) { # tab-delimited, no quotes
        if (fileIsClinical(firstline)){values$inputFileType <- 'clinical'}
        batch <- read.csv(inFile$datapath, header=TRUE, sep='\t', quote='')   #header=input$header, sep=input$sep, quote=input$quote)
      }else{
        firstline<-scan(inFile$datapath, what=character(), sep=",", quote='"', nlines = 1)
        if (length(which(grepl('^gene$',firstline,ignore.case=T,perl=T)))) { # comma-delimited, quotes='"'
          if (fileIsClinical(firstline)){values$inputFileType <- 'clinical'}
          batch <- read.csv(inFile$datapath, header=TRUE, sep=',', quote='"')   #header=input$header, sep=input$sep, quote=input$quote)
        }else{
          error(loggerServer,"Error - invalid input file")
          sendAlertMessage("Invalid input file! \\n\\nPlease use a tab-delimited file with no quotes\\nor a comma-delimited file with quotes (like the Excel .csv files)\\ncontaining at least these 3 fields specified on the first line:\\n\\nSampleID, Gene, Variant",session)
          error = T
        }
      }
      if (!error) {
        #debug(loggerServer,paste("batch has these columns:",paste(colnames(batch),collapse=","),sep=""))
        #debug(loggerServer,paste("batch has",dim(batch)[1],"rows"))
        values$fileRead[inFile$datapath]<-1
        values$batch<-batch
        values$contentType<-'table'
        values$sessionid <- paste(format(Sys.time(), "%d%b%Y_%H%M%S"),sep='_')
        enableActionButton('parseButton',session)
        disableActionButton('loadButton',session)
      }
    }
  })
  # Parse reactive function
  observe({
    if (is.null(input$parseButton) || is.na(input$parseButton) || input$parseButton == 0)
      return()
    if (input$parseButton != values$lastParseButton) {
      values$lastParseButton <- input$parseButton
      disableActionButton('parseButton',session)
      batch <- isolate(values$batch)
      if(!is.null(batch)) {
        progress <- Progress$new(session, min=1, max=10)
        on.exit(progress$close())
        progress$set(message = 'Parsing in progress',
                     detail = '  looking for sample_id, gene and variation')
        if (length(which(is.na(batch[,1])))) {batch <- batch[-which(is.na(batch[,1])),]}
        colnames(batch)[which(grepl('dna.no|sample',colnames(batch),ignore.case=T,perl=T))]='sample_id'
        colnames(batch)[which(grepl('gene',colnames(batch),ignore.case=T,perl=T))]='gene'
        colnames(batch)[which(grepl('nomenclature|variation|variant|hgvs',colnames(batch),ignore.case=T,perl=T))]='variation'
        #debug(loggerServer,paste("After parsing batch has these columns:",paste(colnames(batch),collapse=","), sep=""))
        progress$set(detail = ' cleaning up the variation field', value=5)
        batch[,'variation']=sub(';.*$','',batch[,'variation'],perl=T)
        batch[,'variation']=sub('[_\\(]p.*$','',batch[,'variation'],perl=T)
        batch[,'variation']=sub(',c.*$','',batch[,'variation'],perl=T)
        values$batch<-batch#[which(grepl('BRCA',batch$gene)),]
        values$contentType<-'table'
        info(loggerServer," Parsing done")
        progress$set(detail = ' done!', value=10)
        enableActionButton('loadButton',session)
      }
    }
  })


  # function generating SQL code for insert in database of a loaded set
  make_insert <- function(x) {
    i <- "INSERT INTO dogma_batch (session_id,investigator_name,sample_id,gene,variation"
    if (isolate(values$inputFileType == 'clinical')) {
      i <- paste0(i,',worksheet,inv_id,zygosity,exon')
    }
    i <- paste0(i,") VALUES ('",
                paste(isolate(values$sessionid),isolate(values$iname),x["sample_id"],x["gene"],x["variation"],sep="','"),
                "'"
    )
    if(isolate(values$inputFileType=='clinical')){
      i <- paste0(i,",'",
                  paste(x["Worksheet"],x["InvID"],x["Het.hom"],x["Exon"],sep="','"),
                  "'"
      )
    }
    i <- paste0(i,")")
    i
  }

  # Load & Analyze reactive function
  observe({
    if (is.null(input$loadButton) || is.na(input$loadButton) || input$loadButton == 0)
      return()
    if (input$loadButton != values$lastLoadButton) {
      if (isolate(input$iname=='')) {
        session$sendCustomMessage(type="jsCode",
                                  list(code="alert('Please fill in the investigator name!');$('input#iname').focus();"));
      }else{
        values$iname<-isolate(input$iname)
        values$lastLoadButton <- input$loadButton
        disableActionButton('loadButton',session)
        batch <- isolate(values$batch)
        if(!is.null(batch) && (is.null(values$fileLoaded) || is.na(values$fileLoaded[values$sessionid]))) {
          info(loggerServer," Inserting rows in the database")
          if (length(which(grepl('^$|no mutation',batch[,'variation'],ignore.case=T,perl=T)))) {
            batch <- batch[-which(grepl('^$|no mutation',batch[,'variation'],ignore.case=T,perl=T)),]
          }
          progress <- Progress$new(session, min=1, max=dim(batch)[1])
          on.exit(progress$close())
          progress$set(message = 'Loading in progress',
                       detail = paste('  Inserting',max=dim(batch)[1],'rows in the database'))
          mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
          c <- 0
          for (i in apply(batch,1,make_insert)) {
            c <- c + 1
            progress$set(value = c)
            rs<-dbSendQuery(mydb,i)
          }
          rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                     "update dogma_batch b join cappagenes c
             on b.gene=c.gene
            set b.transcript=c.ensembl
          where b.transcript=''
            and b.session_id ='",isolate(values$sessionid),"'"))
          info(loggerServer," Analyzing variants in the database")
          procs<-c("CALL get_alleles_from_hgvs_cdna('dogma_batch','variation')",
                   "CALL get_varType_from_hgvs_cdna('dogma_batch','variation')",
                   "CALL get_varLocation_from_hgvs_cdna('dogma_batch','variation')",
                   "CALL get_cdnacoord_from_hgvs_cdna('dogma_batch','variation')",
                   "CALL get_genomecoord_from_cdnacoord('dogma_batch')",
                   "CALL determine_codingEffect('dogma_batch')",
                   "CALL make_hgvs('dogma_batch')")
          for (i in procs) {
            progress$set(message = 'Analyzing in progress',detail = unlist(str_split(i,'[(]'))[1])
            ###>>> start, make the line below unfunctional or delete. Because there is no mysql procedure created that we can call them...
            #  rs<-dbSendQuery(mydb,i)
            ###
          }
          info(loggerServer," Analyzing done")
          dbDisconnect(mydb)
          values$fileLoaded[values$sessionid]<-1
          values$sessionschanged<-values$sessionschanged + 1
        }
      }
    }
  })

  # Select an already loaded set reactive function
  observe({
    x<-input$batchid
    #debug(loggerServer,paste("Select loaded set changed to '",x,"'",sep=""))
    if (!is.null(x) && !is.na(x)) {
      if (x !='') {
        progress <- Progress$new(session, min=1, max=10)
        on.exit(progress$close())
        progress$set(message = 'Loading in progress',
                     detail = ' querying database')
        sessid<-unlist(str_split(unlist(str_split(x,' [(]'))[2],'[,]'))[1]
        iname<-unlist(str_split(x,' [(]'))[1]
        values$iname <- iname
        values$sessionid <- sessid
        #debug(loggerServer,paste(" sessid=",sessid," iname=",iname))
        exploreTable<-NULL
        mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
        rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                   "create temporary table tmp_",sessid," as
                                   select gene,variation,
                                          max(concat_ws('\t',date_format(creat_date,'%Y-%m-%d'),clinical_cigma_class)) last_clinical,
                                          group_concat(concat_ws('\t',investigator_name,clinical_cigma_class,
                                                                 concat(date_format(creat_date,'%d/%m/%Y'),' (',
                                                                        if(round(datediff(curdate(),creat_date)/7)<9,
                                                                           concat(round(datediff(curdate(),creat_date)/7),' weeks ago'),
                                                                           concat(round(datediff(curdate(),creat_date)/30.4166),' months ago')
                                                                          ),')' ),
                                                                  notes
                                                                 )
                                                         order by creat_date desc separator '\n'
                                                        ) notescontent
                             from dogma_batchlog group by gene,variation"));
        rs<-dbSendQuery(mydb,paste(sep="",collapse="","create index i1 on tmp_",sessid," (gene,variation)"))
        rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                   "select b.sample_id,
                IF(b.hgvs_prot!='',CONCAT(b.gene,':',b.variation,' (',IF(b.codingEffect='synonymous','p.=',b.hgvs_prot),')')
                                  ,CONCAT(b.gene,':',b.variation)) variation,
                c.cigma_class working_cigma_class,
                c.classification_justification,
                IF(b.clinical_cigma_class!='', b.clinical_cigma_class, substring_index(d.last_clinical,'\t',-1)) clinical_cigma_class,
                IF(b.notes!='',substring_index(substring_index(b.notes,'\t',3),'\t',-1),
                               date_format(str_to_date(substring_index(d.last_clinical,' ',1),'%Y-%m-%d'),'%d/%m/%Y')) clinical_cigma_date,
                d.notescontent, b.notes notes_new,
                b.varLocation, b.varType, b.codingEffect, b.offset,
                b.worksheet, b.inv_id, b.zygosity, b.exon, case when cv.reason is not null then 1 else 0 end common
           from dogma_batch b left join dogma_classification c
             on b.gene=c.gene and b.variation=c.hgvs_cdna and c.version='",version,"'
                left join tmp_",sessid," d
             on b.gene=d.gene and b.variation=d.variation
                left join commonvar cv
             on b.gene=cv.gene and b.variation=cv.hgvs_cdna
           where session_id ='",sessid,"'"))
        exploreTable<-fetch(rs,-1)
        exploreTable[is.na(exploreTable[,3]),3]<-''
        exploreTable[is.na(exploreTable[,4]),4]<-''
        exploreTable[is.na(exploreTable[,5]),5]<-''
        exploreTable[is.na(exploreTable[,6]),6]<-''
        exploreTable[is.na(exploreTable[,7]),7]<-''
        exploreTable[is.na(exploreTable[,8]),8]<-''
        exploreTable[is.na(exploreTable[,13]),13]<-''
        exploreTable[is.na(exploreTable[,14]),14]<-''
        exploreTable[is.na(exploreTable[,15]),15]<-''
        exploreTable[is.na(exploreTable[,16]),16]<-''
        dbDisconnect(mydb)
        progress$set(detail =' creating table', value  = 5)
        exploreTable<-cbind(notes='', exploreTable)
        exploreTable[union(which(exploreTable$notescontent!=''),which(exploreTable$notes_new!='')),'notes']<-'<img src="images/details_open.png">'
        exploreTable$visited <-''
        #debug(loggerServer,paste(" exploreTable has: ",dim(exploreTable)," dimensions"))
        values$exploreTable<-exploreTable
        values$inputFileType<-ifelse(exploreTable[1,14]=='','simple','clinical')
        values$contentType<-'exploreTable'
        disableActionButton('parseButton',session)
        disableActionButton('loadButton',session)
        progress$set(detail =' done!', value = 10)
      }else{
        values$batch<-NULL
        values$contentType<-'table'
        disableActionButton('parseButton',session)
        disableActionButton('loadButton',session)
      }
    }
  })


  # function generating SQL code for update of the notes field
  make_update_notes <- function(x) {
    paste0("UPDATE dogma_batch SET notes='",x["notes_new"],"' WHERE session_id='",x["session_id"],"' AND variation='",x["variation"],"'" )
  }

  # function generating SQL code for update of the clinical_cigma_class field
  make_update_clinical <- function(x) {
    paste0("UPDATE dogma_batch SET clinical_cigma_class='",x["clinical_cigma_class"],"' WHERE session_id='",x["session_id"],"' AND variation='",x["variation"],"'")
  }

  # saveButton reactive function
  observe({
    info(loggerServer,"in function for saveButton")
    if (is.null(input$saveButton) || is.na(input$saveButton) || input$saveButton == 0)
      return()
    #debug(loggerServer,paste("Save button pressed:",input$saveButton,sep=""))
    exploreTable <- isolate(values$exploreTable)
    if(!is.null(exploreTable)) {
      x<-isolate(input$batchid)
      #debug(loggerServer,paste(" Saving batch",x))
      if (length(which(exploreTable$clinical_cigma_class!=''))+length(which(exploreTable$notes_new!=''))>0) {
        batch<-unique(exploreTable[union(which(exploreTable$clinical_cigma_class!=''),which(exploreTable$notes_new!='')),c('variation','clinical_cigma_class','notes_new')])
        batch$session_id<-unlist(str_split(unlist(str_split(x,' [(]'))[2],'[,]'))[1]
        batch$variation<-sub(' .*$','',sub('^.*:','',batch$variation,perl=T),perl=T)
        rows1<-length(which(batch$notes_new!=''))
        rows2<-length(which(batch$clinical_cigma_class!=''))
        progress <- Progress$new(session, min=1, max=rows1+rows2)
        on.exit(progress$close())
        progress$set(message = 'Saving in progress',
                     detail = paste('  Updating',max=rows1+rows2,'rows in the database'))
        mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
        c <- 0
        for (i in apply(batch[which(batch$notes_new!=''),],1,make_update_notes)) {
          c <- c + 1
          progress$set(detail = paste('  Updating notes for:',batch$variation[c]),value = c)
          rs<-dbSendQuery(mydb,i)
        }
        for (i in apply(batch[which(batch$clinical_cigma_class!=''),],1,make_update_clinical)) {
          c <- c + 1
          progress$set(detail = paste('  Updating curated class for:',batch$variation[c-rows1]),value = c)
          rs<-dbSendQuery(mydb,i)
        }
      }
      info(loggerServer," Saving done")
      dbDisconnect(mydb)
      values$batch<-NULL
      values$sessionschanged<-0
      values$lastExcludeDeep <- TRUE
      values$lastExcludeExtra <- TRUE
      values$lastExcludeCommon <- FALSE
      values$lastExcludeIndels <- FALSE
      values$selectedVar <- NULL
      values$selectedVarBaseline <- NULL
      values$selectedVarClinical <- NULL
      values$batchid<-NULL
      values$exploreTable<-NULL
      values$exploreTableiDL<-10
      values$exploreTableiDS<-0
      values$contentType<-'table'
    }
  })

  # Explore table click resulting in input$exploreTableVariation change
  observe({
    if(!is.null(input$exploreTableVariation) && input$exploreTableVariation!='') {
      info(loggerServer,"in function for click on variation.")
      x <- unlist(str_split(input$exploreTableVariation, "\t"))
      if (length(x)){
        values$selectedVar <- gsub(' .*$','',x[1],perl=T)
        values$selectedVarBaseline <- x[2]
        values$selectedVarClinical <- x[3]
        values$selectedVarClinicalCurr <- x[3]
        values$exploreTableiDS <- as.numeric(x[4])
      }
      info(loggerServer,paste(" Selected variation:",isolate(values$selectedVar),sep=""))
      isolate({
        if(length(which(values$exploreTable$visited!=''))) {
          values$exploreTable[which(values$exploreTable$visited!=''),'visited']<-'Y'
          #debug(loggerServer,paste("these rows have been set to Y:",which(values$exploreTable$visited!='')))
        }
        values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'visited']<-'L'
        #debug(loggerServer,paste("these rows have been set to L:",which(grepl(values$selectedVar,values$exploreTable$variation))))
        values$newNoteCurr<-unlist(str_split(values$exploreTable[which(values$exploreTable$variation==x[1]),'notes_new'][1],"\t"))[4]
        if(is.null(values$newNoteCurr) || is.na(values$newNoteCurr)){values$newNoteCurr=''}
      })
      values$contentType <- 'browseTabs'
    }
  })

  # Explore table change length resulting in input$exploreTableiDSiDL change
  observe({
    if(!is.null(input$exploreTableiDSiDL) && input$exploreTableiDSiDL!='') {
      info(loggerServer,"in function for change info on exploreTable")
      x <- unlist(str_split(input$exploreTableiDSiDL, "\t"))
      if (length(x)){
        values$exploreTableiDS <- as.numeric(x[1])
        values$exploreTableiDL <- as.numeric(x[2])
      }
    }
  })

  # Explore (undo icon) reactive function
  observe({
    #debug(loggerServer,paste("in function for exploreButton button=",input$exploreButton," last=",isolate(values$lastExploreButton)))
    if (is.null(input$exploreButton) || is.na(input$exploreButton)) {
      return()
    }
    if (input$exploreButton == 0) {
      values$lastExploreButton<-0
      return()
    }
    if ((input$exploreButton != isolate(values$lastExploreButton))) {
      values$lastExploreButton <- input$exploreButton
      values$contentType <- 'exploreSingle'
    }
  })

  # Batch (undo icon) reactive function
  observe({
    #debug(loggerServer,paste("in function for batchButton button=",input$batchButton," last=",isolate(values$lastBatchButton)))
    if (is.null(input$batchButton) || is.na(input$batchButton)) {
      return()
    }
    if (input$batchButton == 0) {
      values$lastBatchButton<-0
      return()
    }
    if ((input$batchButton != isolate(values$lastBatchButton))) {
      values$lastBatchButton <- input$batchButton
      values$contentType <- 'table'
    }
  })


  ###>>> Login Observe

  # LOGIN function
  observe({
    if (is.null(input$loginButton) || is.na(input$loginButton)) {
      return()
    }

    if (input$loginButton == 0) {
      values$lastloginButton <- 0
      return()
    }

    # lets login button acting when user is filled
    if(input$loginButton!=0 & input$USER==""){
      values$lastloginButton <- input$loginButton
      return()
    }

    if(input$loginButton != values$lastloginButton & input$USER!="" & input$PASS!=""){
      mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)

      #     output$userPass<- renderText({ isolate(input$PASS)})

      query <- paste(sep="", collapse="","SELECT auth,email,name,lname,affiliation,institution,usercode from userpass where email='", gsub("'","''",input$USER), "' AND  password='", gsub("'","''",digest(input$PASS, algo = "md5")),"';")
      rs<-dbSendQuery(mydb,query)
      values$curator <- dataset<-fetch(rs,-1)
      values$contentType <- 'front'

      if(length(dataset$auth)>0){
        if(dataset$auth==1 & tolower(dataset$email)==tolower(input$USER)){

          values$USER <- input$USER
          values$logged <- TRUE
          output$USERv <- renderText({paste0(date(), "\nLogged in :",values$USER)})
          output$USERv2 <- renderText("LOGGED IN")
          values$lastloginButton <- input$loginButton
          disableActionButton2("loginButton",session)
          enableActionButton("logoutButton", session)
        }

      }else{
        values$lastloginButton <- input$loginButton
      }
      dbDisconnect(mydb)
      return()

    }
  })

  ###

  ###>>> Signin Observer

  observe({
    #debug(loggerServer,paste("in function for signin button=",input$signButton," last=",isolate(values$lastsignButton)))
    if (is.null(input$signButton) || is.na(input$signButton)) {
      return()
    }
    if (input$signButton == 0) {
      values$lastsignButton<-0
      return()
    }
    if ((input$signButton != isolate(values$lastsignButton))) {
      values$lastsignButton <- input$signButton
      values$contentType <- 'signin'
    }
  })

  ###


  # Return (undo icon) reactive function
  observe({
    #debug(loggerServer,paste("in function for returnButton button=",input$returnButton," last=",isolate(values$lastReturnButton)))
    if (is.null(input$returnButton) || is.na(input$returnButton)) {
      return()
    }
    if (input$returnButton == 0) {
      values$lastReturnButton<-0
      return()
    }
    if ((input$returnButton != isolate(values$lastReturnButton))) {
      values$lastReturnButton <- input$returnButton
      if(values$contentType=='browseTabs'){
        if (isolate(values$showSaveNotes)){saveNotes()}
        values$newNoteCurr <- ''
        values$showNewNote <- FALSE
        values$showSaveNotes <- FALSE
        values$contentType <- 'exploreTable'
      }else{
        values$contentType <-'front'
      }
    }
  })

  # Mark variation reactive function
  observe({
    #debug(loggerServer,paste("in function for markButton button=",input$markButton," last=",isolate(values$lastmarkButton)))
    if (is.null(input$markButton) || is.na(input$markButton)) {
      return()
    }
    if (input$markButton == 0) {
      values$lastmarkButton<-0
      return()
    }
    if ((input$markButton != isolate(values$lastmarkButton))) {
      values$lastmarkButton <- input$markButton
      if (input$clinicalCIGMA==''){
        sendAlertMessage("Please select a proposed class!",session)
      }else{
        values$newNoteCurr <- isolate({paste0(ifelse(values$newNoteCurr!='',values$newNoteCurr,''),
                                              'REVIEWED on ',format(Sys.time(), '%d/%m/%Y\n'))})
        values$showNewNote <- TRUE
        values$showSaveNotes <- TRUE
      }
    }
  })

  # Clinical CIGMA select change function
  observe({
    x<-input$clinicalCIGMA
    if (!is.null(x) && !is.na(x)) {
      info(loggerServer,paste("Proposed classification changed to '",x,"'",sep=""))
      if (x != isolate(values$selectedVarClinical)) {
        values$newNoteCurr <- isolate({paste0(ifelse(values$newNoteCurr!='',values$newNoteCurr,''),
                                              'Proposed classification changed to ',x,'\n')})
        values$showNewNote <- TRUE
        values$showSaveNotes <- TRUE
      }
      values$selectedVarClinicalCurr <- x
    }
  })



  # Add Note reactive function
  observe({
    #debug(loggerServer,paste("in function for addNoteButton button=",input$addNoteButton," last=",isolate(values$lastaddNoteButton)))
    if (is.null(input$addNoteButton) || is.na(input$addNoteButton)) {
      return()
    }
    if (input$addNoteButton == 0) {
      values$lastaddNoteButton<-0
      return()
    }
    if ((input$addNoteButton != isolate(values$lastaddNoteButton))) {
      values$showNewNote <- TRUE
      values$showSaveNotes <- TRUE
    }
  })

  saveNotes <- function() {
    isolate({
      values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'notes_new'] <-
        paste(values$iname,
              input$clinicalCIGMA,
              format(Sys.time(), '%d/%m/%Y'),
              gsub('\n','<br>',input$newNote,perl=T),
              sep="\t"
        )
      values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'notes'] <- '<img src="images/details_open.png">'
      values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'clinical_cigma_class'] <- input$clinicalCIGMA
      values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'clinical_cigma_date'] <- format(Sys.time(), '%d/%m/%Y')
      values$newNoteCurr<-gsub('\n','<br>',input$newNote,perl=T)
    })
  }


  # Save Notes reactive function
  observe({
    debug(loggerServer,paste("in function for saveNotesButton button=",input$saveNotesButton," last=",isolate(values$lastSaveNotesButton)))
    if (is.null(input$saveNotesButton) || is.na(input$saveNotesButton)) {
      return()
    }
    if (input$saveNotesButton == 0) {
      values$lastsaveNotesButton<-0
      return()
    }
    if ((input$saveNotesButton != isolate(values$lastsaveNotesButton))) {
      if(values$contentType == 'browseTabs' || values$contentType == 'exploreTable'){
        saveNotes()
        ###>>>
        #saveNotes() only updates exploreTable and does not save the updates
        #saveButton

        exploreTable <- isolate(values$exploreTable)
        x<-isolate(input$batchid)
        exploreTable$gene <- sapply(exploreTable$variation,function(x) strsplit(x, split=":")[[1]][1])
        batch<-unique(exploreTable[union(which(exploreTable$clinical_cigma_class!=''),which(exploreTable$notes_new!='')),c('variation','clinical_cigma_class','notes_new','gene')])
        batch$session_id<-unlist(str_split(unlist(str_split(x,' [(]'))[2],'[,]'))[1]
        batch$variation<-sub(' .*$','',sub('^.*:','',batch$variation,perl=T),perl=T)
        rows1<-length(which(batch$notes_new!=''))
        rows2<-length(which(batch$clinical_cigma_class!=''))
        mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
        c <- 0
        for (i in apply(batch[which(batch$notes_new!=''),],1,make_update_notes)) {
          c <- c + 1
          rs<-dbSendQuery(mydb,i)
        }
        for (i in apply(batch[which(batch$clinical_cigma_class!=''),],1,make_update_clinical)) {
          c <- c + 1
          rs<-dbSendQuery(mydb,i)
        }
        dbDisconnect(mydb)
      }

      if(values$contentType == 'exploreSingle' && values$lastsaveNotesButton!=isolate(input$saveNotesButton)){
        #saveNotesSingle <- function(){
        info<-cigmainfo()
        geneupdate <- info$summary["gene",]
        variationupdate <- info$summary["hgvs_cdna",]
        classupdate <- isolate(input$clinicalCIGMAexp)
        dateupdate <- Sys.Date()
        notesupdate <-  substr(isolate(input$newNote), start=1, stop=200)
        notesupdate <- gsub(";|>>>|[|]",",",notesupdate)
        inFile <- input$filenote
        if(length(inFile$name)>0){
          if(file.exists(inFile$datapath)){
            notesupdate <- paste0(notesupdate, " +",length(inFile$name)," File(s) uploaded")
          }
        }
        previousnotes <- info$notes$notes

        if(!is.null(variationupdate)){
          if(classupdate!="" || notesupdate!="" || !is.null(inFile)){
            notesupdate <- paste0(">>> ",dateupdate,"; ", paste0(values$curator$name," ",values$curator$lname),"; ",values$curator$institution ,"; ", ifelse(classupdate=="","Curated class not changed",classupdate),"; ", ifelse(notesupdate=="","No comments.",paste0(notesupdate,".")))
            values$notesupdate <- notesupdate <- paste0(notesupdate," | ",previousnotes)
            notesupdate <- gsub("'","''",notesupdate)

            mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
            checknote <- paste(sep="",collapse="","select * from dogma_batchlog where gene='",geneupdate, "' AND  variation='",variationupdate,"'")
            rs<-dbSendQuery(mydb,checknote)
            dataset <- dataset2 <-fetch(rs,-1)

            inFile <- input$filenote
            values$addFileName <- inFile$name

            if(!is.null(input$filenote) && file.exists(inFile$datapath) ){
              splitgene <- unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))
              gene <- splitgene[1]
              transcript <- splitgene[2]
              variation<-gsub(' [(].*$','',input$variationSel,perl=T)
              values$dirNotes <- paste0("/tmp/NOTES/",gene,"/",gsub("[^[:alnum:]\\_\\-]","",variation))
              dir.create(values$dirNotes, showWarnings = F, recursive = T, mode = "0777")
              addFile <- paste0(values$dirNotes,"/", gene,"_",transcript,"_", gsub(">","-",gsub("[^[:alnum:]\\>\\_\\-]","",variation)),
                                "_",paste0(values$curator[,c(3,4)], collapse=""),"_",gsub("[[:punct:]]|\\s","", round(Sys.time(), units = "mins")),"_",inFile$name)
              file.copy(inFile$datapath, file.path(addFile), overwrite=T)
		system(paste0("rm ",inFile$datapath))
            }

            if(nrow(dataset)>0){
              i_update <- paste(sep="",collapse="","UPDATE dogma_batchlog SET notes='",notesupdate,"', investigator_name='", values$curator$usercode, "', clinical_cigma_class='",classupdate,"', creat_date='",dateupdate,"' where gene='",geneupdate, "' AND  variation='",variationupdate,"'")
              rs<-dbSendQuery(mydb,i_update)
            }else{
              insertnotes <- paste(sep="",collapse="","INSERT INTO dogma_batchlog (session_id,creat_date,investigator_name,sample_id,gene,transcript,variation,working_cigma_class,clinical_cigma_class,notes) VALUES ('",values$sessionid,"',CURDATE(),'",values$curator$usercode,"',","'','",geneupdate,"','",info$summary["transcript",],"','",variationupdate,"','",NULL,"','",classupdate,"','",notesupdate,"')")
              rs<-dbSendQuery(mydb,insertnotes)
            }
            output$sucsaved <- renderText(" ")
            sendAlertMessage("SUCCESSFULLY saved ;) \\n\\nYour comment and/or file upload will appear in CanVar-UK Classifications/Notes tab below",session)
            dbDisconnect(mydb)
          }else{
            sendAlertMessage(":o \\n Please select \\n classification or/and \\n make comments or/and \\n upload file(s)",session)
          }
        }else{
          sendAlertMessage(":(( \\n Please select an existing gene and/or variant",session)
        }
      }
      values$lastsaveNotesButton <- input$saveNotesButton
      ###

    }
  })



  ###################
  # DYNAMIC CONTENT #
  ###################
  # Output holding the dynamic content for the control panel

  ###>>> Welcome option
  # #scroll to scroll opens scroll
  welcome <- span(title="Successful login", verbatimTextOutput("USERv"), tags$head(tags$style("#USERv{color:fuchsia; font-size:15px; font-weight:bold; font-style:italic;
border-top:transparent; border-bottom:transparent; border-left:transparent; border-right:transparent;margin:0px,10px;padding-bottom:0px;,
overflow-y:#scroll;  background: #f5f5f5;}")))

  welcome2 <- span(title="Successful login", verbatimTextOutput("USERv2"), tags$head(tags$style("#USERv2{color:fuchsia; font-size:12px; font-weight:bold; font-style:italic;
border-top:transparent; border-bottom:transparent; border-left:transparent; border-right:transparent;margin-bottom: 0px;
overflow-y:#scroll; background: #f5f5f5;}")))

  disableActionButton('logoutButton',session)
  ###

  ###>>> observe before create account
  observe({
    if (is.null(input$lname) || is.na(input$lname)) {
      return()
    }else{
      output$olname  <- renderText({ isolate(input$lname) })
    }
  })



  ###>>> Admin panel

  observe({
    if (is.null(input$updateAccountButton) || is.na(input$updateAccountButton)) {
      return()
    }
    if (input$updateAccountButton == 0) {
      values$lastupdateAccountButton<-0
      return()
    }

    if (input$updateAccountButton != isolate(values$lastupdateAccountButton)) {
      mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
      query <- paste(sep="", collapse="","SELECT email from userpass where email IS NOT NULL;")
      rs<-dbSendQuery(mydb,query)
      dataset<-fetch(rs,-1)
      email <- tolower(trimws(input$Semail3, which = c("both")))
      existuser  <- any(tolower(dataset$email)==email)
      essentialCheck <- existuser & any(input$Sactivate %in% c(1,2,3))

      if(input$Spassnew1!="" & input$Spassnew1==input$Spassnew2 & essentialCheck){
        acc_update <- paste(sep="",collapse="","UPDATE userpass SET password='",digest(input$Spassnew1, algo="md5"),"', auth=",input$Sactivate," where email='",email,"';")
        output$updateAccMessage <- renderText(paste0("Successful - Account is ", ifelse(input$Sactivate==1,"activated", "inactivated")," and password changed. DO NOT FORGET to send this password to the account holder."))
        rs<-dbSendQuery(mydb,acc_update)
      }else if(essentialCheck & input$Spassnew1=="" & input$Spassnew2==""){
        acc_update <- paste(sep="",collapse="","UPDATE userpass SET auth=",input$Sactivate," where email='",email,"';")
        output$updateAccMessage <- renderText(paste0("Successful - Account is ", ifelse(input$Sactivate==1,"activated", "inactivated")))
        rs<-dbSendQuery(mydb,acc_update)
      }else{
        output$updateAccMessage <- renderText("There is no change done yet. Check the fields above and try it again ")
      }

      dbDisconnect(mydb)
      values$lastupdateAccountButton <- input$updateAccountButton
    }else{ return() }
  })


  registerError <- ""

  sendmail<-function(fname,lname){
     recipients <- c('CanVIG@icr.ac.uk','Subin.Choi@icr.ac.uk')
     sender <- 'canvartempnotice@gmail.com'
     # Full email message in RFC2822 format
     message <- paste0('From: "CanVar UK" <temp@noreply.com>
To: "CanVIG" <CanVIG@icr.ac.uk>
Subject: New User Sign-up alert

Dear CanVar UK Admin team,

A new user ', fname, ' ',lname, ' has signed up on CanVar-UK website.
Please update the access authority of the user.')


     # Send the email
     send_mail(sender, recipients, message, smtp_server = 'smtps://smtp.gmail.com',
               username = 'canvartempnotice@gmail.com', password  = 'SubinChoi@@')
   }



  observe({
    #mydb<-dbConnect(MySQL(),user='batch',password='cigma',dbname='cigma2',host=host)
    if (is.null(input$createAccountButton) || is.na(input$createAccountButton)) {
      return()
    }

    if (input$createAccountButton == 0) {
      values$lastcreateAccountButton<-0
      return()
    }
    if (input$createAccountButton != isolate(values$lastcreateAccountButton)) {
      cm2 <- c()
      cm <- c()
      if(isolate(input$Semail1)!=isolate(input$Semail2)){ cm2 <- c(cm2, "Emails do not match") }
      if(isolate(input$Spass1)!=isolate(input$Spass2)){ cm2 <- c(cm2, "Passwords do not match") }
      if(nchar(input$Semail1)>90 || nchar(input$Semail2)>90  || nchar(input$fname)>20 || nchar(input$lname)>20 || nchar(input$Sdept)>500 || nchar(input$Sinst)>500 || nchar(input$city)>20){
        cm <- c(cm, "Max. char length exceed.") }
      if(input$fname==""){ cm <- c(cm, "First name") }
      if(input$lname==""){ cm <- c(cm, "Last name") }
      if(input$Sdept==""){ cm <- c(cm, "Department") }
      if(input$Sinst==""){ cm <- c(cm, "Institution") }
      if(input$city==""){ cm <- c(cm, "City") }
      if(input$Semail1=="" || input$Semail2==""){ cm <- c(cm, "Email") }
      if(input$Spass1=="" || input$Spass2==""){ cm <- c(cm, "Password") }
      if(is.null(c(cm,cm2))){
        mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
        email <- tolower(trimws(input$Semail1, which = c("both")))
        query <- paste(sep="", collapse="","SELECT email from userpass where email IS NOT NULL;")
        rs<-dbSendQuery(mydb,query)
        dataset<-fetch(rs,-1)
        existuser  <- any(tolower(dataset$email)==email)
        passvalid  <- c(length(grep("@|#|[$]|%|\\^|&|<|>|[*]", input$Spass1))>0 && any(grepl("[[:upper:]]", strsplit(input$Spass1, "")[[1]])) && any(grepl("[[:lower:]]", strsplit(input$Spass1, "")[[1]])) && any(grepl("[[:digit:]]", strsplit(input$Spass1, "")[[1]])) && nchar(input$Spass1)<=16 && nchar(input$Spass1)>=8)

        if(!existuser && passvalid){
          #mydb<-dbConnect(MySQL(),user='batch',password='cigma',dbname='cigma2',host=host)
          query <- paste(sep="", collapse="","SELECT userid FROM userpass ORDER BY userid DESC LIMIT 1;")
          rs<-dbSendQuery(mydb,query)
          dataset<-fetch(rs,-1)

          usercode <- gsub(";|>|'","~",paste0(substr(input$fname, start=1, stop=1),substr(input$lname, start=1, stop=1),dataset$userid+1))
          query <- paste(sep="", collapse="", "INSERT INTO userpass (email, name, lname, affiliation, institution, city, password, auth, usercode, registerdate) VALUES ('",gsub("'","''",email),"', '",gsub("'","''",input$fname),"', '",gsub("'","''",input$lname),"', '", gsub("'","''",input$Sdept),"', '",gsub(";|>>>|[|]",",",gsub("'","''",input$Sinst)),"', '",gsub("'","''",input$city),"', '",gsub("'","''",digest(input$Spass1, algo="md5")),"', FALSE, '",usercode,"', CURDATE()",");")
          rs<-dbSendQuery(mydb,query)
          mymessage <- paste0("Dear ",input$fname," ",input$lname,";
					     You have successfully registered in canvaruk.org
					     Your account will be activated in 24-48 hours and a notification will be sent to ",
                              input$Semail1)
          sendmail(input$fname, input$lname)
          output$registerMessage <- renderText(mymessage)

          disableActionButton('createAccount',session)
        }else{
          if(existuser){
            output$registerMessage <- renderText(paste0(email ," is already registered. If you forgot your password please contact us (CanVIG@icr.ac.uk) with #PASSWORD at the email subject line."))
          }
          if(!passvalid){
            output$registerMessage <- renderText("Invalid password.")
          }
          if(existuser && !passvalid){
            output$registerMessage <- renderText(paste0(email ," is already registered. If you forgot your password please contact CanVIG@icr.ac.uk with #PASSWORD at the email subject line.
                                                        And invalid password."))
          }
        }
        dbDisconnect(mydb)
      }else{

        registerError <- span(title="Register errors", verbatimTextOutput("registerError"), tags$head(tags$style("#registerError{color:red; font-size:12px; font-weight:bold; font-style:italic;
	border-top:transparent; border-bottom:transparent; border-left:transparent; border-right:transparent;
	overflow-y:#scroll; max-height: 20px; background: #f5f5f5;}")))
        isare <- ifelse(length(cm)>1, "are", "is")
        mymessage <- ifelse(length(cm)==0, "" , paste0(paste0(cm, collapse="; ")  , " ", isare, " missing. "))
        output$registerMessage <- renderText((paste0(mymessage, paste0(paste0(cm2, collapse="; ")))))
      }
      values$lastcreateAccountButton <- input$createAccountButton
    }
    #dbDisconnect(mydb)
  })

  ###

  output$uiControl<-renderUI({
    #debug(loggerServer,"in renderUI for uiControl")

    #newpass <- tags$input(id="PASS",type="password",placeholder="Password")

    c<-NULL
    if (values$contentType == 'front') {
      c<-div(class="navbar", style="margin-bottom:5px;",
             div(class="well",style="margin-bottom:0;padding:4px;",
                 ###>>> Login option
                 fluidRow(column(width=8, offset = 4,
                                column(width=12, div(align="center",img(src="images/CanVar-UK_Logo_Final_text_Transparent_background.png",style="width:70%"),span(id='version',strong('VERSION 1.2')))),
                                fluidRow(column(width=12, div(align="right", welcome))),

                                fluidRow(column(width=6, offset=6,div(style="margin:0px 10px;",textInput(inputId = "USER", label="", value="", placeholder="Email")))),
                                fluidRow(column(width=6,offset=6,div(style="margin:0px 10px;",textInput(inputId = "PASS",label="", value="",  placeholder="Password", type="password")))),
                                div(style="margin:0px 10px;",
                                fluidRow(column(width=2,offset=7, actionButton(inputId= "loginButton", label="LOG IN ", icon=icon("lock", lib = "font-awesome"))),
                                        column(width=2, actionButton(inputId= "signButton", label="REGISTER", icon=icon("sign-in", lib = "font-awesome"))))),

                                bsTooltip("loginButton", "You can login if you already registered and if your account has been activated", trigger="focus")
                                )
                                )
                 )
             )

                                 ###>>>
                                 # tooltip html code
                                 # Requires wget https://cran.r-project.org/src/contrib/Archive/shinyBS/shinyBS_0.20.tar.gz
                                 # and R CMD INSTALL shinyBS_0.20.tar.gz
                                 # and library(shinyBS)
                                 # and a dummy line to make it active  bsTooltip("bins", "The tooltip can be written here...", "right", options = list(container = "body"))

                                 # HTML("<button id=\"mybutton\" type=\"button\" class=\"btn btn-default action-button\" style=\"padding:4px; font-size:80%\">?</button>"),
                                 # HTML("<script>$(document).ready(function() {setTimeout(function() {shinyBS.addTooltip( \"signButton\", \"tooltip\", {\"container\": 'body\", \"placement\": \"right\", \"trigger\": \"hover\", \"title\": \":) The wait times will be bro\"})}, 500)});</script>")

             ####

      ###>>> signin form

    }else if(values$contentType == 'signin'){

      if(values$logged & input$USER=="admin@canvaruk.org" & input$PASS=="8bda$496e%74b3d6@ICR" & input$signButton){

        # INSERT INTO userpass (email, name, lname, affiliation, institution, city, password, auth, usercode, registerdate) VALUES ('admin@canvaruk.org', 'admin', 'canvaruk', 'managment', 'controlpanel', 'London', '8bda$496e%74b3d6@ICR', FALSE, 'AA', CURDATE());
        mydb<-dbConnect(MySQL(),user='batch',password='cigma',dbname='cigma2',host=host)
        email <- tolower(trimws(input$emailbox3, which = c("both")))
        query <- paste(sep="", collapse="","SELECT email,name,lname,institution,auth from userpass where email IS NOT NULL AND email NOT LIKE 'tmp%' AND email NOT LIKE '%demo%' AND email NOT LIKE 'test%' AND email NOT LIKE '%_DOWN%' AND email NOT IN ('%tmp%','%_DOWN%','admin@canvaruk.org','clare.turnbull@genomicsengland.co.uk','cankutcubuk@gmail.com', 'CanVIG@clare','clare.turnbull@gmail.com');")
        rs<-dbSendQuery(mydb,query)
        dataset<-fetch(rs,-1)
        dbDisconnect(mydb)

        output$downloadData <- downloadHandler(
          filename = function() {
            paste("canvarUK-users-", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            write.csv(dataset, file, row.names=F)
          }
        )

        downlink <- downloadLink("downloadData", label = "Download the list of registered users")
        emailbox3 <- div(aligh="center", tags$textarea(id="Semail3", rows=1, cols=10, placeholder="Email to in/activate",
                                                       style="margin: 0px 0px 10px; width: 32%; height: 22px; resize: none"))
        accountAct <- div(aligh="center", tags$textarea(id="Sactivate", rows=1, cols=10, placeholder="Activation code",
                                                        style="margin: 0px 0px 10px; width: 32%; height: 22px; resize: none"))
        passNew1 <- div(aligh="center", tags$textarea(id="Spassnew1", rows=1, cols=10, placeholder="New password",
                                                      style="margin: 0px 0px 10px; width: 32%; height: 22px; resize: none"))
        passNew2 <- div(aligh="center", tags$textarea(id="Spassnew2", rows=1, cols=10, placeholder="Repeat the new password",
                                                      style="margin: 0px 0px 10px; width: 32%; height: 22px; resize: none"))

        updateAccount <- span(title="Update acconut", do.call(actionButton,list(inputId="updateAccountButton",label="",icon=icon("save","fa-2x",lib="font-awesome"))))

        c<-div(class=" navbar navbar-top", style="margin-bottom:5px;",
               div(class="well",style="margin-bottom:0;;padding:100px;",

                   fluidRow(div(align="center", downlink)),
                   fluidRow(column(align="center",width=12,"----------------------------------------------------")),
                   fluidRow(div(align="center", style="font-style:italic;font-size:17.5px;margin-bottom:10px","Update an account")),
                   # fluidRow(column(offset=0, width=10, textOutput("olname"))),
                   fluidRow(div(align="center", emailbox3)),
                   fluidRow(div(align="center", accountAct)),
                   fluidRow(div(align="center", style="font-style:italic;font-size:17.5px;margin-bottom:10px","If you do not want to update password, leave the the password boxes blank")),
                   fluidRow(div(align="center", passNew1)),
                   fluidRow(div(align="center", passNew2)),
                   fluidRow(div(align="center", updateAccount)),
                   fluidRow(div(align="center", textOutput("updateAccMessage"))),
                   br(),br(),br(), br(),br(),br(),
                   fluidRow(div(align="center", style="font-style:italic;font-size:17.5px;margin-bottom:10px","Return to main page"),
                            div(align="center",span(title="Return to main page",actionButton("returnButton","",icon=icon("home","fa-2x",lib="font-awesome")))),
                            div(align="center",width=12,"----------------------------------------------------"))
               )
        )
        # i_update <- paste(sep="",collapse="","UPDATE dogma_batchlog SET notes='",notesupdate,"', investigator_name='", values$curator$usercode, "', clinical_cigma_class='",classupdate,"', creat_date='",dateupdate,"' where gene='",geneupdate, "' AND  variation='",variationupdate,"'")
        # rs<-dbSendQuery(mydb,i_update)
      }else{

        emailbox1 <- div(aligh="center", tags$textarea(id="Semail1", rows=1, cols=10, placeholder="Email",
                                                       style="margin: 0px 0px 10px; width: 32%; height: 22px; resize: none"))
        emailbox2 <- div(aligh="center", tags$textarea(id="Semail2", rows=1, cols=10, placeholder="Confirm Email",
                                                       style="margin: 0px 0px 10px; width: 32%; height: 22px; resize: none"))
        deptbox <- div(aligh="center", tags$textarea(id="Sdept", rows=1, cols=10, placeholder="Department",
                                                     style="margin: 0px 0px 10px; width: 32%; height: 22px; resize: none"))
        instbox <- div(aligh="center", tags$textarea(id="Sinst", rows=1, cols=10, placeholder="Institution",
                                                     style="margin: 0px 0px 10px; width: 32%; height: 22px; resize: none"))

        createAccount <- span(title="Create your account", do.call(actionButton,list(inputId="createAccountButton",label="",icon=icon("save","fa-2x"))))

        c<-div(class=" navbar navbar-top", style="margin-bottom:5px;",
               div(class="well",style="margin-bottom:0;;padding:100px;",
                   fluidRow(column(align="center",width=12,"----------------------------------------------------")),
                   fluidRow(div(align="center", style="font-style:italic;font-size:17.5px;margin-bottom:10px","Create an account")),
                   fluidRow(column(width=6, fluidRow(div(align="right", textInput(inputId = "fname", label="", value="", placeholder="First Name")))),
                            column(width=6, fluidRow(div(align="left", textInput(inputId = "lname", label="", value="", type="text", placeholder="Last Name"))))

                   ),
                   # fluidRow(column(offset=0, width=10, textOutput("olname"))),
                   fluidRow(div(align="center", deptbox)),
                   fluidRow(div(align="center", instbox)),
                   fluidRow(div(align="center", textInput(inputId = "city", label="", value="", placeholder="City"))),
                   fluidRow(div(align="center", emailbox1)),
                   fluidRow(div(align="center", emailbox2)),
                   fluidRow(column(width=6,
                                   fluidRow(div(align="right", textInput(inputId = "Spass1", label="", value="", placeholder="Password"))),
                                   fluidRow(div(align="right", textInput(inputId = "Spass2", label="", value="", placeholder="Confirm Password")))),
                            column(width=6,helpText(em('Password must include:')),
                                   helpText(em(' - Minimum 8 maximum 16 characters')),
                                   helpText(em(' - Upper and lower case letters')),
                                   helpText(em(' - Numeric character(s)')),
                                   helpText(em(' - Special character(s): @#$%^&<>*')))
                   ),
                   fluidRow(div(align="center", createAccount)),
                   registerError,
                   fluidRow(div(align="center", textOutput("registerMessage"))),
                   fluidRow(column(align="center",width=12,"----------------------------------------------------")),
                   fluidRow(column(width=12,"")),
                   fluidRow(column(width=12,"")),
                   fluidRow(div(align="center", style="font-style:italic;font-size:17.5px;margin-bottom:10px","Return to main page"),
                            div(align="center",span(title="Return to main page",actionButton("returnButton","",icon=icon("home","fa-2x")))),
                            div(align="center",width=12,"----------------------------------------------------"))

               )
        )
      }
      ###
    }else if (values$contentType == 'exploreSingle') {

      ###>>> get the baseline_class and newnote function

      output$baseline_class <- renderText({
        info<-cigmainfo()
        info$baseline_class
      })


      #})

      ###>>> reset the values when the gene is updated on exploreSinglepage.
      # observe({
      #   if (!is.null(input$genetrans){
      #     values$lastgenetrans <- input$genetrans
      #     values$exploreTableiDS <- 0
      #   }
      # })
      ###


      output$downloadadditionalData <- downloadHandler(

        filename = function() {
          # !!! You can move them outside and  make them as global variables
          splitgene <- unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))
          gene <- splitgene[1]
          variation<-gsub(' [(].*$','',input$variationSel,perl=T)
          #            values$dirNotes <- paste0("/tmp/NOTES/",gene"/",gsub("[[:punct:]]","",input$variationText))
          #            dirNotes <- paste0("/tmp/NOTES/",gene,"/",gsub("[[:punct:]]","", variation))
          dirNotes <- paste0(gene,"/",gsub("[^[:alnum:]\\_\\-]","",variation))
          #            paste0(gsub("/","_",dirNotes),".zip")
          paste0(dirNotes,".zip")
        },

        content = function(file) {

          splitgene <- unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))
          gene <- splitgene[1]
          variation<-gsub(' [(].*$','',input$variationSel,perl=T)
          #           values$dirNotes <- paste0("/tmp/NOTES/",gene,"/",gsub("[[:punct:]]","",input$variationText))
          #	    dirNotes <- paste0("/tmp/NOTES/",gene,"/",gsub("[[:punct:]]","",variation))
          dirNotes <- paste0("/tmp/NOTES/",gene,"/",gsub("[^[:alnum:]\\_\\-]","",variation))
#        tared<-paste0(dirNotes,".tar")
	ziped <- paste0(dirNotes,".zip")
          zip(zipfile=ziped, files=list.files(dirNotes, full.names = T, pattern="pdf|doc|docx"), flags="-j")
          file.copy(ziped, file)
#	system(paste0("'","tar -cjf /tmp/NOTES/",gene,"/",gsub("[^[:alnum:]\\_\\-]","",variation),".tar -C ",dirNotes,"'"))
	#	tar(file,dirNotes)
        },
        # https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_types/Complete_list_of_MIME_types
        contentType = "application/zip"
        #     contentType = NULL
        #    contentType = "application/msword"
      )


      ###

      saveNotes <- ''
      additionalFile <- ''

      if(values$showSaveNotesExp){

        if(values$logged  & !values$USER %in% paste0("tmp",c(1:20))){

          saveNotes<-span(title=paste("Save changes made to variant",isolate(values$selectedVar)),
                          do.call(actionButton,list(inputId="saveNotesButton",label="",icon=icon("save","fa-2x"))))


          additionalFile <-  div(style="margin-top:10px;", title="Choose Additional File",
                                  do.call(fileInput,list(inputId="filenote",multiple = TRUE, label="Choose Additional File",accept = c(".docx",".doc",".pdf"))))

          output$sucsaved <- renderText(" ")
          sucsaved <- span(title="Successful save", verbatimTextOutput("sucsaved"),tags$head(tags$style("#sucsaved{color:red; font-size:10px; font-weight:bold; font-style:italic;
	border-top:transparent; border-bottom:transparent; border-left:transparent; border-right:transparent; padding:0px; margin-bottom:0px;
	overflow-y:#scroll; max-height: 20px; background: #f5f5f5;}")))

        }else{

          saveNotes<-span(title=paste("Please LOGIN Save changes made to variant",isolate(values$selectedVar)),
                          do.call(actionButton,list(inputId="saveNotesButtonMessage",label="",icon=icon("save","fa-2x"))))

          output$sucsaved <- renderText("Please log in to save")
          sucsaved <- span(title="Successful save", verbatimTextOutput("sucsaved"), tags$head(tags$style("#sucsaved{color:red; font-size:10px; font-weight:bold; font-style:italic;
	border-top:transparent; border-bottom:transparent; border-left:transparent; border-right:transparent; padding:0px; margin-bottom:0px;
	overflow-y:#scroll; max-height: 20px; background: #f5f5f5;}")))


          #disableActionButton('saveNotes',session)

        }
      }


      if(values$showNewNoteExp) {
        newNote <- tags$textarea(id="newNote", rows=6, cols=70, style="margin: 0px 0px; width:100%; height: 80px; resize: none", gsub('<br>','\n',values$newNoteCurr2),
                                 placeholder="Maximum 200 characters", maxlength="200")
      }else{
        newNote <- span(title=paste("Add note for variant",isolate(values$selectedVar)),actionButton("addNoteButton","",icon=icon("plus-square-o","fa-2x")))
      }
      ###

      c<-div(class=" navbar navbar-top", style="margin-bottom:5px;",
             div(class="well",style="margin-bottom:0;padding:20px;",
                 fluidRow(
                   column(width=4, div(align="center",img(src="images/CanVar-UK_Logo_Final_Transparent_background.png", style="width:100%;margin-top: 10px;"),
                                       div(align="center", welcome2))),

                 column(width=2,


                          helpText(em('Examples: ')),
                                 helpText('BRCA1 c.181T>G'),
                                    # helpText('BRCA1 c.557C>A'),
                                 helpText('BRCA2 c.7994A>G'),
                                    # helpText('BRCA1   p.Arg1203Ter'),
                                   # helpText('BRCA1 c.5207T>C'),
                          helpText('MLH1 c.131C>T'),
                                 helpText('MSH2 c.2500G>A'),
                                 helpText('TP53 p.Arg337His'),

                        span(title="Return to main page",actionButton("returnButton","",icon=icon("home","fa-2x"))),
                        span(title="Download PDF with variant report",downloadButton("downloadVariantReport","")
                        )

                          #

                   ),


                   #####

                   # column(width=2,
                   #        helpText(em("Classes:")),
                   #        helpText("1-Non-pathogenic"),
                   #        helpText("2-Likely non-pathogenic"),
                   #        helpText("3-Uncertain/review"),
                   #        helpText("3*-Uncertain-high suspicion"),
                   #        helpText("4-Likely pathogenic"),
                   #        helpText("5-Pathogenic")
                   # ),
                   column(width=3,
                            selectInput(inputId = "genetrans", label="Gene: ", choices = c('',gl), selectize=FALSE),
                            textInput(inputId = "variationText", label="Variant starting with: " ),
                            selectInput(inputId = "variationSel", label="Select matching variant:", choices = c(''), selectize=FALSE),
                            div(column(offset=0, width=10,textOutput("varCountText")))),
                                   # here the offset was 6 to shift the text


                            #		       column(width=2, "Automated classification:", verbatimTextOutput("baseline_class"),
                            #		      	 	selectInput("clinicalCIGMAexp", "Curated classification:",
                            #			   	choices=c("", "1-Non-pathogenic", "2-Likely non-pathogenic", "3-Uncertain/review", "3*-Uncertain-high suspicion" ,"4-Likely pathogenic", "5-Pathogenic"),
                            #		           selectize=FALSE, selected=values$selectedVarClinicalCurr2)
                            #                  ),


                            column(width=3,


                                   #      fluidRow("Automated classification:"),
                                   #      fluidRow(verbatimTextOutput("baseline_class")),
                                   selectInput("clinicalCIGMAexp", "Proposed classification:",
                                               choices=c("", "1-Non-pathogenic", "2-Likely non-pathogenic", "3-Uncertain/review", "3-Uncertain-high suspicion" ,"4-Likely pathogenic", "5-Pathogenic"),
                                               selectize=FALSE, selected=""),

                                   #			       fluidRow(fileInput("filenote", "Choose Additional File", accept = c("text/csv","text/comma-separated-values,text/plain",".csv",".doc","docx"))),
                                   #			       fluidRow(downloadButton("downloadData", label = "Download"))
                                   #			       fluidRow(downloadLink("downloadData", label = "Download"))
                                   div('Note:'),
					newNote,
                                   #			saveNotes
                                   #              		fluidRow(span(title="Save changes made to the set",actionButton("saveButton","",icon=icon("save","fa-2x"))))
                                   additionalFile,
                                   div(saveNotes, sucsaved)
                                   # 				   fluidRow(saveNotes),
                                   #                                   fluidRow(sucsaved)

                                   )
                            )
                          )
             )




    }
    else if (values$contentType == 'table') {
      c<-div(class=" navbar navbar-top", style="margin-bottom:5px;",
             div(class="well",style="margin-bottom:0;;padding:4px;",
                 fluidRow(
                   column(width=4, div(align="center",img(src="images/CanVar-UK_Logo_Final_Transparent_background.png", style="width:100%;margin-top: 20px;"),
                                       div(align="center", welcome2),
                                       span(title="Return to main page",actionButton("returnButton","",icon=icon("home","fa-2x"))),
                                       )),
                   column(width=8,
                          fluidRow(
                            column(width=4, textInput(inputId = "iname", label="Investigator name:",value=values$iname)
                            ),
                            column(width=4, fileInput('uploadFile', 'Upload variants:',
                                                      accept=c('text/csv','text/comma-separated-value,text/plain','application/csv','.csv')),
                                   actionButton('parseButton',"Parse"),
                                   actionButton('loadButton',"Load & Analyze")
                            ),
                            column(width=4, uiOutput('loadedSets')
                            )

                            )
                          )
                   )
                 )
      )


    }else if(values$contentType == 'exploreTable') {
      c<-div(class=" navbar navbar-top", style="margin-bottom:5px;",
             div(class="well",style="margin-bottom:0;padding:4px;",
                 fluidRow(
                 column(width=2,img(src="images/CanVar-UK_Icon_Transparent.png")),
                   column(width=10,
                          fluidRow(
                            column(width=4,
                                   div(style="font-style:italic;font-size:17.5px;margin-bottom:10px","Set:"),
                                   div(style="font-weight:bold;font-size:17.5px;",textOutput("selectedSet")),
                                   div(style="margin-top:20px;",
                                       span(title="Return to main page",actionButton("returnButton","",icon=icon("home","fa-2x"))),
                                       span(title="Save changes made to the set",actionButton("saveButton","",icon=icon("save","fa-2x"))),
                                       span(title="Download the annotated set",downloadButton("downloadAnnotatedBatch",""))
                                   )
                            ),
                            column(width=5,
                                   checkboxInput("excludeDeep","Exclude intronic variants ",values$lastExcludeDeep),
                                   do.call(numericInput, list("offsetLimit","more than:",isolate(values$offsetLimit))),span(style="width:200px;","nt away from an exon border"),
                                   checkboxInput("excludeExtra","Exclude UTRs and intergenic",values$lastExcludeExtra),
                                   checkboxInput("excludeCommon","Exclude common variants",values$lastExcludeCommon),
                                   checkboxInput("excludeIndels","Exclude insertions/deletions",values$lastExcludeIndels)
                            ),
                            column(width=3,
                                   helpText("1-Non-pathogenic"),
                                   helpText("2-Likely non-pathogenic"),
                                   helpText("3-Uncertain/review"),
                                   helpText("3*-Uncertain-hig suspicion"),
                                   helpText("4-Likely pathogenic"),
                                   helpText("5-Pathogenic")
                            )
                          )
                   )
                 )
             )
      )
    }else if (values$contentType == 'browseTabs') {
      saveNotes <- ''
      if(values$showSaveNotes){
        saveNotes<-span(title=paste("Save changes made to variant",isolate(values$selectedVar)),
                        do.call(actionButton,list(inputId="saveNotesButton",label="",icon=icon("save","fa-2x"))))
      }
      newNote <- span(title=paste("Add note for variant",isolate(values$selectedVar)),actionButton("addNoteButton","",icon=icon("plus-square-o","fa-2x")))
      if(values$showNewNote) {
        newNote <- tags$textarea(id="newNote", rows=6, cols=70, style="margin: 0px 0px 10px; width: 95%; height: 130px; resize: none", gsub('<br>','\n',values$newNoteCurr))
      }
      c<-div(class=" navbar navbar-top", style="margin-bottom:5px;",
             div(class="well",style="margin-bottom:0;padding:4px;",
                 fluidRow(
                   column(width=2,img(src="images/CanVar-UK_Icon_Transparent.png")),
                   column(width=2,
                          div(style="font-style:italic;font-size:17.5px;","Variant:"),
                          div(style="margin-top:5px;font-weight:bold;font-size:17.5px;",textOutput("selectedVar")),
                          div(style="margin-top:20px;",
                              span(title=paste("Return to set:",isolate(input$batchid)),actionButton("returnButton","",icon=icon("reply","fa-2x"))),
                              span(title=paste("Mark variant",isolate(values$selectedVar),'as "REVIEWED on',format(Sys.time(), '%d/%m/%Y"')),actionButton("markButton","",icon=icon("check-square-o","fa-2x"))),
                              saveNotes
                          )
                   ),
                   column(width=2,
                          #      tags$style("#selectedVarBaseline {padding:4px;margin-bottom:5px;} #clinicalCIGMA {margin-bottom:5px;}"),
                          #      "Automated classification:", verbatimTextOutput("selectedVarBaseline"),
                          selectInput("clinicalCIGMA","Curated classification:",choices=c("",
                                                                                          "1-Non-pathogenic",
                                                                                          "2-Likely non-pathogenic",
                                                                                          "3-Uncertain",
                                                                                          "3*-Uncertain-hig suspicion",
                                                                                          "4-Likely pathogenic",
                                                                                          "5-Pathogenic"
                          ),
                          selected=values$selectedVarClinicalCurr, selectize=FALSE)
                   ),
                   column(width=3, newNote
                   ),
                   column(width=2,
                          helpText("1-Non-pathogenic"),
                          helpText("2-Likely non-pathogenic"),
                          helpText("3-Uncertain/review"),
                          helpText("3*-Uncertain-hig suspicion"),
                          helpText("4-Likely pathogenic"),
                          helpText("5-Pathogenic")
                   ),
                   column(width=1,
                          span(title="Download PDF with variant report",downloadButton("downloadVariantReport",""))
                   )
                 )
             )
      )
    }
    return(tagList(tags$script(src = "js/resize_download.js"),c))
  })

  # Output holding the dynamic content for the main panel
  output$uiContents<-renderUI({
    #debug(loggerServer,"in renderUI for output$uiContent")
    #debug(loggerServer,paste(" values$contentType='",values$contentType,"'",sep=''))
    if (values$contentType == 'front') {
      div(class="well",style="padding-top:20px;",
          fluidRow(
        #  column(width=3, 

#fluidRow(
 #                    div(align="center", style="margin-top:48px;",img(src="images/CanVar-UK_Icon_Transparent.png",style="height:200px")))),
            column(width=3,
                   fluidRow(
                     div(align="center", style="margin-top:48px;",
                         helpText('Explore a single variant'),
                         actionButton('exploreButton',"Explore", icon=icon("search", lib = "font-awesome")),
                         div(style="margin-top:16px;",helpText('Explore a batch of variants'),
                             actionButton('batchButton'," Batch ",icon=icon("table", lib = "font-awesome"))
                         )
                     )
                   )
            ),
            column(width=6,
                   div(style="margin-top:6px;",
                       helpText(h4('Instructions:')),
                       helpText('- In ',strong('[Explore]'),' mode, select a gene and type the HGVS name of a single variant you wish to query (beginning with "c." or "p.")'),
                       helpText('- In ',strong('[Batch]'),' mode, upload a tab-delimited .txt or comma-separated .csv file with a header line including: ',em('Sample, Gene, Variant')),
                       fluidRow(column(offset=1, width=11,
                                       helpText('- Click on ',strong('[Parse]'),' to process and ',strong('[Load&Analyze]'),' to load it in the database. A list of previously loaded sets can be viewed via the drop-down: any of these can be selected for re-examination.'),
                                       helpText('- Click on one variant name in the table for detailed information from related data sources'))),
                       helpText('- Click on ',img(src="images/details_open.png"),' for clinically curated notes on variant'),
                       helpText('- The ',strong('automated_class'),' is generated by the gene-specific decision-tree, with an associated ',strong('classification_justification'),'. The ',strong('curated_class'),' is assigned following clinicial review +/- literature search.')
                   )
            )
          ),
          fluidRow(div(style="margin-top:6px; margin-left:20px;",
                       helpText(h4('Disclaimer:')),


                       helpText(em(p('The Cancer Variant Database is a repository of annotations for variants in cancer predisposition genes',br(),'Classifications, both automated and curated, are derived from integration of the variant level data according to objective criteria.'),
                                   p('This resource is currently under development and is for research use only.',br(),'We welcome your feedback at: CanVIG@icr.ac.uk',br(),br(), "CanVIG resources: shorturl.at/eqxPW")))
          )
          )
      )
    }
    else if (values$contentType == 'table' && !is.null(values$batch)) {
      #debug(loggerServer,paste(" values$batch columns:",paste(colnames(values$batch),collapse=","),sep=''))
      #debug(loggerServer,paste(" values$batch has ",dim(values$batch)[1]," rows.",sep=''))
      tagList(
        singleton(
          tags$head(
            tags$style(type="text/css","div.dataTables_length label {float: right !important; text-align: right !important;}
                                    div.dataTables_info {float: right ! important;}
                                    div.dataTables_filter label {float: left !important;}
                                    div.dataTables_paginate {float: left !important; margin: 0;}"
            ))),
        div(class="well",style="padding-top:10px;",dataTableOutput('batch'))
      )
    }else if (values$contentType == 'exploreTable') {
      #debug(loggerServer,paste(" exploreTable has:",dim(values$exploreTable)," dimensions"))
      div(class="well",style="padding-top:10px;",
          selDataTableOutput("exploreTable")
      )
    }else if (values$contentType == 'browseTabs' || values$contentType == 'exploreSingle') {
      #debug(loggerServer,paste("selected variation=",values$selectedVar))
      info <- cigmainfo()
      notes_new <- NULL
      if  (values$contentType == 'browseTabs'){
        notes_new <- values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'notes_new']
      }
      tabs <- vector(mode="list")

      if (!is.null(info) && !is.null(info$summary)) {
        tabs[[length(tabs)+1]] <-  tabPanel("Summary", uiOutput('summaryPanel'))
      }

      if (!is.null(info) && (!is.null(info$alamut) ||
          !is.null(info$muttaster) ||
          !is.null(info$polyphen2) ||
          !is.null(info$cadd) ||
          !is.null(info$suspect) ||
          !is.null(info$revel) ||
          !is.null(info$gavin) ||
          !is.null(info$info$tp53_fut))){
        tabs[[length(tabs)+1]] <- tabPanel("In silico Predictions", uiOutput('insilicoPanel'))
        }

      if (!is.null(info) &&
          (!is.null(info$evs) ||
           !is.null(info$tgp) ||
           !is.null(info$tgp_phase3) ||
           # !is.null(info$icr) ||
           # !is.null(info$bic) ||
           # !is.null(info$lovd3) ||
           # !is.null(info$dmudb) ||
           # !is.null(info$umd) ||
           ###>>> gnomad issue ###
           # otherwise if the variant has only gnomad data then  it does not appear
           # e.g. SDHB c.783A>G
           !is.null(info$gnomad)
           ###
          )
      ) {
        tabs[[length(tabs)+1]] <- tabPanel("Control Frequency", uiOutput('controlfrequencyPanel'))
      }
      if (!is.null(info) &&
          # (!is.null(info$evs) ||
          #  !is.null(info$tgp) ||
          #  !is.null(info$tgp_phase3) ||
           !is.null(info$icr) ||
           !is.null(info$bic) ||
           !is.null(info$lovd3) ||
           !is.null(info$dmudb) ||
           !is.null(info$umd) ||
           ###>>> gnomad issue ###
           # otherwise if the variant has only gnomad data then  it does not appear
           # e.g. SDHB c.783A>G
           !is.null(info$phe)
           ###
          )
       {
        tabs[[length(tabs)+1]] <- tabPanel("Case Frequency", uiOutput('casefrequencyPanel'))
      }

      if (!is.null(info) &&
          (!is.null(info$easton) ||
          !is.null(info[['bouwman_2020_epi']]) ||
          !is.null(info[['parsons_2019']]) ||
          !is.null(info$lindor) ||
          !is.null(info$iarc) ||   #NEWBUG because of insight_class. Here the lovd stays as lovd and not lovd3. Becasue lovd data is not shown in this tab and has no negative effect
          (!is.null(info[['lovd']]) && ifelse(any(rownames(info[['lovd']])=="insight_class"),info[['lovd']]['insight_class',]!='', FALSE))
          )
      ) {
        tabs[[length(tabs)+1]] <- tabPanel("Genetic Epidemiology", uiOutput('geneticepiPanel'))
      }
      ###>>> wessex added below
      if (!is.null(info) &&
          (!is.null(info$houdayer) ||
           !is.null(info$walker) ||
           !is.null(info$wessex)
          )
      ) {
        tabs[[length(tabs)+1]]<- tabPanel("Splicing Analysis", uiOutput('splicingPanel'))
        ###
      }
      if (!is.null(info) &&
          (!is.null(info[["guidugli"]])  ||
           !is.null(info$drost) ||
           !is.null(info[['bouwman']]) ||
           ###>>>  ####
           !is.null(info[['bouwman_2020_func']]) ||
           !is.null(info$findlay) ||
           !is.null(info[["guidugli2018"]]) ||
           !is.null(info$tp53_sge) || !is.null(info$tp53_kato) || !is.null(info$tp53_fut))

          ####

      ) {
        tabs[[length(tabs)+1]] <- tabPanel("Functional Analysis", uiOutput('functionalPanel'))
      }


      # notestabs <- vector(mode="list")
      if (!is.null(info) && !is.null(info$notes) && !is.na(info$notes)) { #&& !is.null(info$notes)) || (!is.null(notes_new) &&!is.null(notes_new[1]) && notes_new[1]!='')) {
        tabs[[length(tabs)+1]]<- tabPanel("CanVar-UK Classifications/notes", uiOutput('notesPanel'))
      }

      div(class="well",style="padding-top:10px;",
          # fluidRow(column(width=12, div(style="font-size:15.0px;border-color:white;border-width:3px;border-style:solid;",  do.call(tabsetPanel, c(id="notespanel", notestabs))))),
          # hr(),
          fluidRow(
            ###>>> font size (13.0px)
            # the top:270px was before 170 but to see panel names (like Summary, Frequency) we changed it
            column(width=12, div(style="font-size:15.0px", do.call(tabsetPanel, c(id="tabpanel", tabs)))),
            #	    column(width=8,do.call(tabsetPanel, c(id="tabpanel", tabs))),
            #	    column(width=4,do.call(tabsetPanel, c(id="notespanel", notestabs)))
            ###
          )
      )
    }else{
      NULL
    }
  })

  #################################################
  # 1. Plain dataTable with loaded/parsed content #
  #    contentType = "table"                      #
  #################################################
  output$batch<-shiny::renderDataTable({
    values$batch
  }, options = function () {
    list(sDom = "<'row-fluid'<'span6'p><'span6'i>>t<'row-fluid'<'span6'f><'span6'l>r>",
         bDestroy = TRUE,
         bLengthChange = TRUE,
         iDisplayLength = values$exploreTableiDL,
         aLengthMenu= list(10, 25, 50, 100),
         fnInfoCallback = I("(function( oSettings, iStart, iEnd, iMax, iTotal, sPre ) {Shiny.onInputChange('exploreTableiDSiDL',oSettings._iDisplayStart+'\t'+oSettings._iDisplayLength);return(sPre);})")
    )}
  )

  ####################################
  # 2. dataTable with hidden fields  #
  #    contentType = "exploreTable"  #
  ####################################

  output$exploreTable <- shiny::renderDataTable({
    v<-values$exploreTable
    info(loggerServer,"in renderDataTable for output$exploreTable")
    #debug(loggerServer,paste(" exploreTable has:",dim(v)," dimensions"))
    if (is.null(v)) return(v);
    if(input$excludeDeep && length(which(v$varLocation == 'intron' & abs(v$offset)>input$offsetLimit))) {
      v<-v[-which(v$varLocation == 'intron' & abs(v$offset)>input$offsetLimit),]
    }
    if(input$excludeExtra && length(which(v$varLocation %in% c("3'UTR","5'UTR","intergenic")))) {
      v<-v[-which(v$varLocation %in% c("3'UTR","5'UTR","intergenic")),]
    }
    if(input$excludeCommon && length(which(as.logical(v$common)))) {
      v<-v[-which(as.logical(v$common)),]
    }
    if(input$excludeIndels && length(which(v$varType %in% c('insertion','deletion','complex')))) {
      v<-v[-which(v$varType %in% c('insertion','deletion','complex')),]
    }
    colnames(v)[1]<-'Notes'
    colnames(v)[2]<-'SampleID'
    colnames(v)[3]<-'Gene:Variation'
    colnames(v)[4]<-'Automated classification'
    colnames(v)[5]<-'Justification (automated class)'
    colnames(v)[6]<-'Proposed classification'
    colnames(v)[7]<-'Date'
    v
  },options = function () {
    list(sDom = "<'row-fluid'<'span5'p><'span3'T><'span4'i>>t<'row-fluid'<'span5'f><'span3'T><'span4'l>r>",
         bDestroy = TRUE,
         bLengthChange = TRUE,
         iDisplayLength = values$exploreTableiDL,
         iDisplayStart = values$exploreTableiDS,
         aLengthMenu= list(10, 25, 50, 100),
         bSortClasses = TRUE,
         bAutoWidth = FALSE,
         aoColumnDefs = list(list(bVisible=FALSE, aTargets=sapply(c(7:18),list)),
                             list(sWidth='275px', aTargets=sapply(2,list)),
                             list(sWidth='50px' , aTargets=sapply(0:1,list)),
                             list(sWidth='225px' , aTargets=sapply(4,list)),
                             list(sWidth='160px', aTargets=sapply(c(3,5),list))),
         fnRowCallback  = I("(function (nRow,aData,iDisplayIndex,iDisplayIndexFull){if(aData[18]=='Y'){$(nRow).addClass('visited');};if(aData[18]=='L'){$(nRow).addClass('visitedLast')}})"),
         fnInfoCallback = I("(function( oSettings, iStart, iEnd, iMax, iTotal, sPre ) {Shiny.onInputChange('exploreTableiDSiDL',oSettings._iDisplayStart+'\t'+oSettings._iDisplayLength);return(sPre);})")
         #                                                          var p=0;
         #                                                          if(typeof Shiny.shinyapp.$inputValues.exploreTablePage != 'undefined'){p=Shiny.shinyapp.$inputValues.exploreTablePage};
         #                                                          var t=oSettings.oInstance;
         #                                                          if(Math.ceil(oSettings._iDisplayStart/oSettings._iDisplayLength)!=p){var tId=setTimeout(function(){t.fnPageChange(p,true)},100);}})")#,
         #            oTableTools = I("({'sSwfPath': 'TableTools-2.2.0/swf/copy_csv_xls_pdf.swf',
         #                              'aButtons': [{'sExtends':'text',
         #                                            'sButtonText':'Show all notes',
         #                                            'fnClick':function(button,config){
         #                                                           var oTable=getExploreTable();
         #                                                           $('#exploreTable tbody tr td img').each(function(){
         #                                                             var nTr=$(this).parents('tr')[0];
         #                                                             this.src='images/details_close.png';
         #                                                             oTable.fnOpen(nTr,fnFormatDetails(nTr),'details');
         #                                                           })
         #                                                      }
         #                                           },
         #                                           {'sExtends':'text',
         #                                            'sButtonText':'Hide all notes',
         #                                            'fnClick':function(button,config){
         #                                                           var oTable=getExploreTable();
         #                                                           $('#exploreTable tbody tr td img').each(function() {
         #                                                             var nTr=$(this).parents('tr')[0];
         #                                                             this.src='images/details_open.png';
         #                                                             oTable.fnClose(nTr);
         #                                                           })
         #                                                        }
         #                                            }
         #                                           ]
         #                               })")
    )
  }
  )

  observe({
    if (!is.null(input$excludeDeep) && (input$excludeDeep != isolate(values$lastExcludeDeep))) {
      values$lastExcludeDeep <- input$excludeDeep
      values$exploreTableiDS <- 0
    }
  })

  observe({
    if (!is.null(input$excludeExtra) && (input$excludeExtra != isolate(values$lastExcludeExtra))) {
      values$lastExcludeExtra <- input$excludeExtra
      values$exploreTableiDS <- 0
    }
  })

  observe({
    if (!is.null(input$excludeCommon) && (input$excludeCommon != isolate(values$lastExcludeCommon))) {
      values$lastExcludeCommon <- input$excludeCommon
      values$exploreTableiDS <- 0
    }
  })

  observe({
    if (!is.null(input$excludeIndels) && (input$excludeIndels != isolate(values$lastExcludeIndels))) {
      values$lastExcludeIndels <- input$excludeIndels
      values$exploreTableiDS <- 0
    }
  })


  ############################################
  # 3. Dynamic tabsetPanel for browsing info #
  #    contentType = "browseTabs"            #
  ############################################

  # Reactive function returning the info from CIGMA2 for the requested variation
  cigmainfo <- reactive({
    if (values$contentType=='exploreSingle' && !is.null(input$genetrans)){
      if (!nchar(input$genetrans)) {return()}
      gene<-unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))[1]
      transcript<-unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))[2]
      variation<-gsub(' [(].*$','',input$variationSel,perl=T)
    }else if (values$contentType=='browseTabs'){
      gene<-unlist(str_split(values$selectedVar,':'))[1]
      variation<-unlist(str_split(values$selectedVar,':'))[2]
      transcript <- NULL
    }else{
      return()
    }
    info(loggerServer,paste("in cigmainfo, gene=",gene," variation=",variation,sep=""))
    i<-NULL
    i$gene <- gene
    i$transcript <- transcript
    i$variation <- variation
    dataset<-NULL
    if (nchar(gene) && nchar(variation)) {
      mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
      rs<-dbSendQuery(mydb,paste(sep="",collapse="",
"select m.source,m.gene,m.transcript,c.transcript,r.transcript,hgvs_cdna,hgvs_prot,hgvs_prot_code1,
                m.altname,varLocation,r.name,varType,codingEffect,hg19_chr,hg19_pos,rsID,
                m.cdna_pos,offset,ntwt,ntmut,codon,aawt,aamut,c.firstorlast3,m.flags,pubmed,google_search_new,google_search_old
           from main_stable m
      left join cdna2genomic c on m.transcript=c.transcript and m.gene=c.gene and m.cdna_pos=c.cdna_pos and m.varLocation='exon'
      left join capparegions r on m.transcript=r.transcript and m.gene=r.gene and m.hg19_chr=r.chr and m.hg19_pos between r.exonstart and r.exonend          
where m.gene='",gene,"' and m.hgvs_cdna='",variation,"'"))
      dataset<-fetch(rs,-1)
      if (!is.null(dataset) && dim(dataset)[1]>0) {

        revelDat <- dataset[,c("hg19_chr","hg19_pos","hgvs_prot_code1")]

        i$rsid<-dataset[,'rsID']
        i$google_search_old<-dataset[,dim(dataset)[2]]
        dataset<-dataset[,-dim(dataset)[2]]
        i$google_search_new<-dataset[,dim(dataset)[2]]
        dataset<-dataset[,-dim(dataset)[2]]
        i$pubmed<-dataset[,dim(dataset)[2]]
        dataset<-dataset[,-dim(dataset)[2]]
        i$summary<-as.data.frame(t(dataset),stringsAsFactors=F)
        rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                   "SELECT type,IF(aa_from=aa_to,CONCAT(description,' [',aa_from,']'),CONCAT(description,' [',aa_from,'-',aa_to,']')) description
           from main_stable m
           join swissprot s on m.gene=s.gene and m.codon between s.aa_from and s.aa_to
          where m.gene='",gene,"' and m.hgvs_cdna='",variation,"'"))
        dataset<-fetch(rs,-1)
        if (!is.null(dataset) && dim(dataset)[1]>0) {
          i$swissprot<-as.data.frame(dataset,stringsAsFactors=F)
        }
        varsources <-strsplit(i$summary['source',],' ')[[1]]
        for (source in values$sources) {
          table <- tolower(source)
          if (source=='EVS') table <- 'esp_cappa'
          if (source=='ICR') table <- 'icr_controls'
          if (source=='HGMD') table <- 'hgmd_cappa'
          prefix <- table
          prefix <- sub('_cappa$','',prefix)
          prefix <- paste0(prefix,'_')

          ###>>> Findlay data query ####
          if(source=='FINDLAY'){
            query <- paste(sep="", collapse="","select score,class,probability,meanRNA from Findlay where gene='",gene,
                           "' AND  hgvs_cdna='",variation,"'")
            rs <- dbSendQuery(mydb,query)
            Findlay <- fetch(rs,-1)
            if(is.na(Findlay[1,1])){ Findlay <- NULL }
            i[[tolower(source)]] <- Findlay
          }
          ###
          ###>>> Bouwman2020 data query ####
          if(source=='BOUWMAN_2020'){
            query <- paste(sep="", collapse="","select Cisplatin_probability_deleterious, Cisplatin_prediction, Olaparib_probability_deleterious, Olaparib_prediction, DR_GFP_probability_deleterious, DR_GFP_prediction, Combined_LR, Posterior_Probability, Bouwman_Parsons_combined_LR, Bouwman_Parsons_combined_Posterior_Probability from bouwman2020 where gene='",gene,
                           "' AND  hgvs_cdna='",variation,"'")
            rs <- dbSendQuery(mydb,query)
            bouwman_2020 <- fetch(rs,-1)
            if(nrow(bouwman_2020)==0){bouwman_2020  <- NULL}
            bouwman_2020_func<-bouwman_2020[1:6]
            bouwman_2020_epi <-bouwman_2020[7:10]
            if(which(bouwman_2020_func != "NA") %>% length() ==0){bouwman_2020_func  <- NULL}

            if(which(bouwman_2020_epi != "NA") %>% length() ==0){ bouwman_2020_epi  <- NULL }
            i[['bouwman_2020_func']] <- bouwman_2020_func
            i[['bouwman_2020_epi']] <- bouwman_2020_epi
          }
          ###



          ###>>> Parsons2019 data query ####
          if(source=='PARSONS_2019'){
            query <- paste(sep="", collapse="","select Combined_LR, Prior_Probability_of_Pathogenicity, Posterior_Probability,IARC_Class, Comment from parsons2019 where gene='",gene,
                           "' AND  hgvs_cdna='",variation,"'")
            rs <- dbSendQuery(mydb,query)
            parsons_2019 <- fetch(rs,-1)
            if(is.na(parsons_2019[1,1])){ parsons_2019  <- NULL }
            i[[tolower(source)]] <- parsons_2019
          }
          ###



          ###>>> gnomad query ####
          i$gnomad <- "show gnomad always"
          ###

          ###>>> PHE query ####
          i$phe <- "show phe always"
          ###

          ###>>> Guidugli2018_BRCA2 data query ####

          if(source=='GUIDUGLI2018'){

            query <- paste(sep="", collapse="","select probability,classification from guidugli2018 where gene='",gene,
                           "' AND  hgvs_cdna='",variation,"'")
            rs <- dbSendQuery(mydb,query)
            guidugli2018 <- fetch(rs,-1)
            if(is.na(guidugli2018[1,1])){ guidugli2018 <- NULL }
            i[[tolower(source)]] <- guidugli2018
          }

          ###

          ###>>> wessex data
          if(source=='WESSEX'){
            query <- paste(sep = "",collapse = "", "select specimen, spliceabnormality, sampletype from WESSEX where gene='",gene,"' and hgvs_cdna='",variation,"'")
            rs<-dbSendQuery(mydb, query)
            wessex <- fetch(rs,-1)
            if(nrow(wessex)==0){ wessex <- NULL }
            i[[tolower(source)]] <- wessex
          }
          ###

          ###>>> add info REVEL and GAVIN
          if(source=='REVEL'){
            query <- paste(sep = "",collapse = "", "select REVEL from REVEL where chr='",revelDat$hg19_chr,"' and hg19_pos='",revelDat$hg19_pos,
                           "' and aachange='",gsub("p.|[0-9]+","",revelDat$hgvs_prot_code1),"'")
            rs<-dbSendQuery(mydb, query)
            revel <- fetch(rs,-1)
            if(nrow(revel)==0){ revel <- NULL }
            # if(nrow(revel)>1){ revel =revel$REVEL[which(paste0(revel$ref,revel$alt)==gsub("c.|[0-9]+|>","",i$variation)
            #                                             |paste0(revel$ref,revel$alt) ==gsub("c.|[0-9]+|>","",i$variation) %>% stringi::stri_reverse())]
            #  }
            i[[tolower(source)]] <- as.data.frame(unique(revel))
          }

          if(source=='GAVIN'){
            query <- paste(sep = "",collapse = "", "select gavinClass from GAVIN  where gene='",gene,"' and hgvs_cdna='",variation,"'")
            rs<-dbSendQuery(mydb, query)
            gavin <- fetch(rs,-1)
            if(nrow(gavin)==0){ gavin <- NULL }
            i[[tolower(source)]] <- gavin
          }

          if(source=='LOVD3'){
            query <- paste(sep = "",collapse = "", "select times_reported from LOVD3  where gene='",gene,"' and hgvs_cdna='",variation,"'")
            rs<-dbSendQuery(mydb, query)
            LOVD3 <- fetch(rs,-1)
            if(nrow(LOVD3)==0){ LOVD3 <- NULL }
            i[[tolower(source)]] <- LOVD3
          }


          ###

          ###>>> Giacomelli2018 (TP53SGE), Kato2003 and Fortuno2019

          if(gene=="TP53"){
            if(source=='TP53_SGE'){
              query <- paste(sep = "",collapse = "", "select Combined_Model,ClassSGE from TP53_SGE where hgvs_prot_code1='",gsub("p.","",revelDat$hgvs_prot_code1),"'")
              rs<-dbSendQuery(mydb, query)
              tp53sge <- fetch(rs,-1)

              if(nrow(tp53sge)==0){ tp53sge <- NULL }
              i[[tolower(source)]] <- tp53sge
            }

            if(source=='TP53_Kato'){
              query <- paste(sep = "",collapse = "", "select TransactivationClass from TP53_Kato where hgvs_prot_code1='",gsub("p.","",revelDat$hgvs_prot_code1),"'")
              rs<-dbSendQuery(mydb, query)
              tp53kato <- fetch(rs,-1)
              if(nrow(tp53kato)==0){ tp53kato <- NULL }
              i[[tolower(source)]] <- tp53kato
            }

            if(source=='TP53_Fut'){
              query <- paste(sep = "",collapse = "", "select BayesDel,SuggestedPrediction from TP53_Fut where hgvs_prot_code1='",gsub("p.","",revelDat$hgvs_prot_code1),"'")
              rs<-dbSendQuery(mydb, query)
              tp53fut <- fetch(rs,-1)
              if(nrow(tp53fut)==0){ tp53fut <- NULL }
              i[[tolower(source)]] <- tp53fut
            }

          }

          ###


          if(source %in% varsources) {
            dataset<-NULL
            rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                       "select * from ",table," where gene='",gene,"' and hgvs_cdna='",variation,"'"))
            dataset<-fetch(rs,-1)
            dataset<-as.data.frame(t(dataset),stringsAsFactors=F)
            selectedrows<-grep(prefix,rownames(dataset))
            i[[tolower(source)]]<-as.data.frame(row.names=sub(paste0('^',prefix),'',rownames(dataset)[selectedrows]),
                                                dataset[selectedrows,],stringsAsFactors=F)


            if (source == 'TGP') {
              rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                         "select tgp_pop,tgp_AF,tgp_alleleCount,tgp_alleleTotal,tgp_genotypeCount,tgp_genotypePopSize from tgp_pop where tgp_rsID='",i$tgp["rsID",1],"'"))
              dataset<-fetch(rs,-1)
              i$tgp_pop<-as.data.frame(dataset,stringsAsFactors=F)
              colnames(i$tgp_pop)<-sub('^tgp_','',colnames(dataset))
              i$tgp_pop<-t(i$tgp_pop[order(i$tgp_pop$genotypePopSize,decreasing=T),])
            }
            if (source == 'TGP_PHASE3') {
              rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                         "select tgp_phase3_pop,tgp_phase3_AF,tgp_phase3_alleleCount,tgp_phase3_alleleTotal,tgp_phase3_genotypeCount,tgp_phase3_genotypePopSize from tgp_phase3_pop where tgp_phase3_id='",i$tgp_phase3["id",1],"'"))
              dataset<-fetch(rs,-1)
              i$tgp_phase3_pop<-as.data.frame(dataset,stringsAsFactors=F)
              colnames(i$tgp_phase3_pop)<-sub('^tgp_phase3_','',colnames(dataset))
              i$tgp_phase3_pop<-t(i$tgp_phase3_pop[order(i$tgp_phase3_pop$genotypePopSize,decreasing=T),])
            }
          }
        }
      }
      rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                 "select investigator_name,
             clinical_cigma_class,
             concat(date_format(creat_date,'%d/%m/%Y'),' (',
                    if(round(datediff(curdate(),creat_date)/7)<9,
                      concat(round(datediff(curdate(),creat_date)/7),' weeks ago'),
                      concat(round(datediff(curdate(),creat_date)/30.4166),' months ago')
                    ),')') date,
             notes
        from dogma_batchlog b
       where b.gene='",gene,"' and b.variation='",variation,"'
       order by creat_date desc"));
      dataset<-fetch(rs,-1)
      if (!is.null(dataset) && dim(dataset)[1]>0) {
        i$notes<-as.data.frame(dataset,stringsAsFactors=F)
      }
      rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                 "select cigma_class, classification_justification
        from dogma_classification
       where gene='",gene,"' and hgvs_cdna='",variation,"' and version='",version,"'"));
      dataset<-fetch(rs,-1)
      if (!is.null(dataset) && dim(dataset)[1]>0) {
        i$baseline_class<-dataset[1,1]
        i$justification<-dataset[1,2]
      }
      dbDisconnect(mydb)
    }
    return(i)
  })

  # Output for the Summary panel
  output$summaryPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in summaryPanel")
    scholarlink1<-a(class="btn btn-primary",href=info$google_search_new,target="_blank","Google Scholar")
    scholarlink2<-a(class="btn btn-primary",href=info$google_search_old,target="_blank","Google Scholar (old nomen.)")
    search_new<-gsub('/scholar','/search',gsub('http://scholar.','http://www.',info$google_search_new))
    searchlink<-a(class="btn btn-primary",href=search_new,target="_blank","Google Search")
    pubmedlink <- ''
    if (!is.null(info) && !is.null(info$pubmed) && info$pubmed!=''){
      pubmedlink <- a(class="btn btn-primary",href=info$pubmed,target="_blank","HGMD Pubmed IDs")
    }
    dbsnplink <- ''
    clinvarReview <- NULL
    # if (!is.null(info) && !is.na(info$rsid) && info$rsid!='character(0)'){
    #   dbsnplink <- a(class="btn btn-primary",href=paste('http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',substr(info$rsid,3,nchar(info$rsid)),'#Diversity',sep=''),target="_blank","dbSNP")
    # }

    ###>>> clinvar update and add review status
    ### !!! rs numbers seems to be same for all substitutions in a codon
    # clinvarlink <- a(class="btn btn-primary",href=paste('http://www.ncbi.nlm.nih.gov/clinvar?term=',info$rsid,sep=''),target="_blank","ClinVar")
    mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
    #  rs<-dbSendQuery(mydb,paste(sep="",collapse="", "select * from CLINVAR where gene='",info$gene,"' and rsID='",info$rsid,"' and hgvs_cdna='",info$variation,"'"))
    # rs<-dbSendQuery(mydb,paste(sep="",collapse="", "select * from CLINVAR where gene='",info$gene,"' and hgvs_cdna='",info$variation,"'"))
    rs<-dbSendQuery(mydb,paste(sep="",collapse="", "select VariationID from CLINVAR where gene='",info$gene,"' and hgvs_cdna='",info$variation,"'"))

    dataset<-fetch(rs,-1)
    dbDisconnect(mydb)


    # clinvarDate <- trimws(strsplit(dataset$Stars[1], split="[(]")[[1]][1], which="right")
    
if(nrow(dataset)==1){

clinvarlink <- a(class="btn btn-primary",href=paste('https://www.ncbi.nlm.nih.gov/clinvar/variation/',dataset$VariationID,sep=''),target="_blank","Clinvar")
}else {clinvarlink=""}

   # if(nrow(dataset)==1){
      ###>>> Clinvar clinvarlink clinvarReview cStarEmpty cStarFull
    #  cStarEmpty <- HTML("<span  style=\"color:lightgrey\" <i class=\"fa fa-star\" aria-hidden=\"false\"></i> </span>")
     # cStarFull <- HTML("<span  style=\"color:black\" <i class=\"fa fa-star\" aria-hidden=\"false\"></i> </span>")
      ## cStarFull <- icon("star", lib = "font-awesome")
      #  clinvarsummary=entrez_summary(db="clinvar", id=dataset$VariationID)['clinical_significance']
       # clinvar_last_evaluated <-  clinvarsummary$clinical_significance$last_evaluated %>% substr(1,10)
       # clinvar_review_status <-clinvarsummary$clinical_significance$review_status
       # clinvar_description <-clinvarsummary$clinical_significance$description



      #if(clinvar_review_status=="practice guideline"){
       # clinvarReview <- list(p(em(paste0("ClinVar status (live sync): ")),strong(clinvar_description)," (",cStarFull,cStarFull,cStarFull,cStarFull,")",paste0(" - Last evaluated by ClinVar:", clinvar_last_evaluated)))
     # }else if(clinvar_review_status=="reviewed by expert panel"){
      #  clinvarReview <- list(p(em(paste0("ClinVar status (live sync): ")),strong(clinvar_description)," (",cStarFull,cStarFull,cStarFull, cStarEmpty,")",paste0(" - Last evaluated by ClinVar:", clinvar_last_evaluated)))
     # }else if(clinvar_review_status=="criteria provided, multiple submitters, no conflicts"){
      #  clinvarReview <- list(p(em(paste0("ClinVar status (live sync): ")),strong(clinvar_description), " (",cStarFull,cStarFull,cStarEmpty,cStarEmpty,")",paste0(" - Last evaluated by ClinVar:", clinvar_last_evaluated)))
     # }else if(clinvar_review_status=="criteria provided, conflicting interpretations"
      #         |clinvar_review_status=="criteria provided, single submitter"){
       # clinvarReview <- list(p(em(paste0("ClinVar status (live sync): ")),strong(clinvar_description), " (",cStarFull,cStarEmpty,cStarEmpty,cStarEmpty,")",paste0(" - Last evaluated by ClinVar:", clinvar_last_evaluated)))
     # }else{
      #  clinvarReview <- list(p(em(paste0("ClinVar status (live sync): ")),strong(clinvar_description), " (",cStarEmpty,cStarEmpty,cStarEmpty,cStarEmpty,")",paste0(" - Last evaluated by ClinVar:", clinvar_last_evaluated)))
     # }
#}
##old
    #}else if(nrow(dataset)>1){
      # e.g https://www.ncbi.nlm.nih.gov/clinvar/?term=TP53%3Ac.1009C%3ET
     # clinvarlink <- a(class="btn btn-primary",href=paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=",info$gene,"%3A",info$variation),target="_blank","ClinVar")
      #clinvarReview <- list(p(em(paste0("ClinVar status (live sync): ")), strong("This variant (HGVS cDNA) has multiple entries in ClinVar. Please check it on ClinVar")))
    #}else{
      # clinvarlink <- a(class="btn btn-primary",href=paste('https://ncbi.nlm.nih.gov/clinvar?term=',info$rsid,sep=''),target="_blank","ClinVar")
     # clinvarlink <- ""
      #clinvarReview <- list(p(em("ClinVar status (live sync): "), strong("this variant does not exist in ClinVar")))
    #}

    #   if(length(grep("100",dataset$Stars))>0){
    #     clinvarReview <- list(p(em(paste0("Clinvar status ( Last evaluated:", clinvar_last_evaluated,"): ")),strong(dataset$Class)," (",cStarFull,cStarFull,cStarFull,cStarFull,")"))
    #   }else if(length(grep("75",dataset$Stars))>0){
    #     clinvarReview <- list(p(em(paste0("Clinvar status ( Last evaluated:", clinvar_last_evaluated,"): ")),strong(dataset$Class)," (",cStarFull,cStarFull,cStarFull, cStarEmpty,")"))
    #   }else if(length(grep("50",dataset$Stars))>0){
    #     clinvarReview <- list(p(em(paste0("Clinvar status ( Last evaluated:", clinvar_last_evaluated,"): ")),strong(dataset$Class), " (",cStarFull,cStarFull,cStarEmpty,cStarEmpty,")"))
    #   }else if(length(grep("25",dataset$Stars))>0){
    #     clinvarReview <- list(p(em(paste0("Clinvar status ( Last evaluated:", clinvar_last_evaluated,"): ")),strong(dataset$Class), " (",cStarFull,cStarEmpty,cStarEmpty,cStarEmpty,")"))
    #   }else{
    #     clinvarReview <- list(p(em(paste0("Clinvar status ( Last evaluated:", clinvar_last_evaluated,"): ")),strong(dataset$Class), " (",cStarEmpty,cStarEmpty,cStarEmpty,cStarEmpty,")"))
    #   }
    #
    # }else if(nrow(dataset)>1){
    #   # e.g https://www.ncbi.nlm.nih.gov/clinvar/?term=TP53%3Ac.1009C%3ET
    #   clinvarlink <- a(class="btn btn-primary",href=paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=",info$gene,"%3A",info$variation),target="_blank","Clinvar")
    #   clinvarReview <- list(p(em(paste0("Clinvar status (", clinvarDate,"): ")), strong("This variant (HGVS cDNA) has multiple entries in Clinvar. Please check it on Clinvar")))
    # }else{
    #   # clinvarlink <- a(class="btn btn-primary",href=paste('https://ncbi.nlm.nih.gov/clinvar?term=',info$rsid,sep=''),target="_blank","ClinVar")
    #   clinvarlink <- ""
    #   clinvarReview <- list(p(em(paste0("Clinvar status (", clinvarDate,"): ")), strong("this variant does not exist in Clinvar")))
    # }

    ###

    info.main<-as.list(info$summary[,1])
    names(info.main)<-rownames(info$summary)
   gene=info.main[['gene']]    
info.main[['name']] = gsub(info.main[['gene']],'',info.main[['name']])
    info.main[['name']] = gsub('on','on ',info.main[['name']])
    info.main[['name']] = gsub('of',' of ',info.main[['name']])
    firstorlast3<-fixempty(info.main[["firstorlast3"]])
    if (firstorlast3>0) firstorlast3<-paste0('+',as.character(firstorlast3))
    rmhclass<-''
    if (!is.null(info) && !is.null(info$rmh)) rmhclass<- info$rmh['class',1]
    sources<-gsub('ALAMUT', 'INSILICO', info.main[["source"]])
    sources<-gsub(' POLYPHEN2| MUTTASTER| CADD| SUSPECT','',sources)
    l<-list(tags$table(border="0",cellspacing="5",
                       tags$tr(tags$td(em(fields[["main"]][["gene"]]),             align="right"),tags$td(strong(info.main[["gene"]])),
                               tags$td(em(fields[["main"]][["hgvs_cdna"]]),        align="right"),tags$td(strong(info.main[["hgvs_cdna"]])),
                               tags$td(em(fields[["main"]][["varLocation"]]),      align="right"),tags$td(strong(info.main[["varLocation"]]))),
                       tags$tr(tags$td(em(fields[["main"]][["transcript"]]),       align="right"),tags$td(strong(info.main[["transcript"]])),
                               tags$td(em(fields[["main"]][["hgvs_prot"]]),        align="right"),tags$td(strong(fixempty(info.main[["hgvs_prot"]]))),
                               tags$td(em(fields[["main"]][["varType"]]),          align="right"),tags$td(strong(info.main[["varType"]]))),
                       tags$tr(tags$td(em(fields[["main"]][["rtranscript"]]),      align="right"),tags$td(strong(info.main[["rtranscript"]])),
                               tags$td(em(fields[["main"]][["hgvs_prot_code1"]]),  align="right"),tags$td(strong(fixempty(info.main[["hgvs_prot_code1"]]))),
                               tags$td(em(fields[["main"]][["codingEffect"]]),     align="right"),tags$td(strong(fixempty(info.main[["codingEffect"]])))),
                       tags$tr(tags$td(em(fields[["main"]][["hg19_chr"]]),         align="right"),tags$td(strong(info.main[["hg19_chr"]])),
                               tags$td(em(fields[["main"]][["altname"]]),          align="right"),tags$td(strong(fixempty(info.main[["altname"]]))),
                               tags$td(em(fields[["main"]][["cdna_pos"]]),         align="right"),tags$td(strong(info.main[["cdna_pos"]]))),
                       tags$tr(tags$td(em(fields[["main"]][["hg19_pos"]]),         align="right"),tags$td(strong(info.main[["hg19_pos"]])),
                               tags$td(em(fields[["main"]][["rsID"]]),             align="right"),tags$td(strong(fixempty(info.main[["rsID"]]))),
                               tags$td(em(fields[["main"]][["offset"]]),           align="right"),tags$td(strong(fixempty(info.main[["offset"]])))),
                       tags$tr(tags$td(em(fields[["main"]][["ntwt"]]),             align="right"),tags$td(strong(info.main[["ntwt"]])),
                               tags$td(em("Exon/intron:"),                         align="right"),tags$td(strong(info.main[["name"]])),
                               tags$td(em(fields[["main"]][["codon"]]),            align="right"),tags$td(strong(fixempty(info.main[["codon"]])))),
                       tags$tr(tags$td(em(fields[["main"]][["ntmut"]]),            align="right"),tags$td(strong(info.main[["ntmut"]])),
                               tags$td(),tags$td(),
                               tags$td(em(fields[["main"]][["firstorlast3"]]),     align="right"),tags$td(strong(firstorlast3)))
    ),
    br(),
    p(em('Data sources:'),strong(sources)),
    br(),
    p(em('Literature links:'),scholarlink1,scholarlink2,searchlink,pubmedlink),
    p(em('Database links:'),dbsnplink,clinvarlink),
    # p(em('Flags:'),strong(fixempty(info.main[["flags"]]))),
    # p(em('RMH class:'),strong(rmhclass)),
    br()
    )
    l1<-list()
    ###>>>
    if(!is.null(clinvarReview)){ l1<-c(l1,clinvarReview) }
    #
    if (!is.null(info$hgmd) && info$hgmd['tag',]!='') {l1<-c(l1,list(p(em('HGMD classification (2015):'),strong(info$hgmd['tag',]))))}
    if (!is.null(info$dmudb) && info$dmudb['interpretation',]!='') {l1<-c(l1,list(p(em('DMuDB classification (2015):'),strong(info$dmudb['interpretation',]))))}
    if (!is.null(info$umd) && info$umd['significance',]!='') {l1<-c(l1,list(p(em('UMD classification (2015):'),strong(info$umd['significance',]))))}
    l<-c(l,l1)
    return(as.character(tagList(l)))
  })

  # Output for the In silico predictions panel
  output$insilicoPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in insilicoPanel")
    l1a<-vector(mode="list")
    l2c<-vector(mode="list")
    l2s<-vector(mode="list")
    if (!is.null(info) && !is.null(info$alamut)) {
      fields.alamut<-names(fields[["alamut"]])
      info.alamut<-as.list(info$alamut[,1])
      names(info.alamut)<-rownames(info$alamut)
      fields.missense<-fields.alamut[which(grepl('AGVGD|SIFT|MAPP',fields.alamut))]
      fields.conserv<-fields.alamut[which(grepl('Orthos|Conserved|BLOSUM',fields.alamut))]
      fields.splicing<-fields.alamut[which(grepl('SS|MaxEnt|NNS|SSF',fields.alamut))]
      l1a<-lapply(fields.missense,function(x){
        v<-info.alamut[[x]]
        if(v!='') {
          return(tags$tr(class=ifelse(x %in% c('AGVGDclass','SIFTprediction','MAPPprediction'),'highlight','normal'),tags$td(em(fields[["alamut"]][[x]]),align="right"),tags$td(strong(v))))
        }else{
          return('')
        }
      }
      )
      l2c<-lapply(fields.conserv,function(x){
        v<-info.alamut[[x]]
        if(v!='') {
          return(tags$tr(tags$td(em(fields[["alamut"]][[x]]),align="right"),tags$td(strong(v))))
        }else{
          return('')
        }
      }
      )
      l2s<-lapply(fields.splicing,function(x){
        v<-info.alamut[[x]]
        ###>>> !is.na(v) added to control the NAs that come from SSF. e.g. TP53 c.794T>A
        if(v!='' && !is.na(v)){
          return(tags$tr(tags$td(em(fields[["alamut"]][[x]]),align="right"),tags$td(strong(v))))
        }else{
          return('')
        }
      }
      )
    }
    ###>>> edited for tooltips example
    #    muttasterInfo  <- icon("info", lib = "font-awesome")   class=\"btn btn-default action-button\"
    muttasterButton <- HTML("<button id=\"muttasterBtn\" type=\"button\" class=\"btn btn-default fa fa-info\" style=\"padding:0.8px 5px; font-size:80%; border-radius: 60%; color: #4CAF50; border: 1px dashed #4CAF50\"></button>")
    muttasterButtonTip <- bsTooltip("muttasterBtn", "Tips for the Mutation Taster tool",trigger="hover")
    l1m<-vector(mode="list")
    if (!is.null(info) && !is.null(info$muttaster)) {
      fields.muttaster<-names(fields[["muttaster"]])
      info.muttaster<-as.list(info$muttaster[,1])
      names(info.muttaster)<-rownames(info$muttaster)
      l1m<-lapply(fields.muttaster,function(x){
        return(
          if(x=="prediction"){
            tags$tr(class=ifelse(x == 'prediction','highlight','normal'), tags$td(em(fields[["muttaster"]][[x]]), muttasterButton, align="right"),  muttasterButtonTip, tags$td(strong(info.muttaster[[x]])))
          }else{
            tags$tr(class=ifelse(x == 'prediction','highlight','normal'), tags$td(em(fields[["muttaster"]][[x]]), align="right"), tags$td(strong(info.muttaster[[x]])))
          }
        )
      }
      )
    }
    ###
    l1p<-vector(mode="list")
    if (!is.null(info) && !is.null(info$polyphen2)) {
      fields.polyphen2<-names(fields[["polyphen2"]])
      info.polyphen2<-as.list(info$polyphen2[,1])
      names(info.polyphen2)<-rownames(info$polyphen2)
      l1p<-lapply(fields.polyphen2,function(x){
        return(tags$tr(class=ifelse(x == 'hvar_prediction','highlight','normal'),tags$td(em(fields[["polyphen2"]][[x]]),align="right"),tags$td(strong(info.polyphen2[[x]]))))
      }
      )
    }
    l1c<-vector(mode="list")
    if (!is.null(info) && !is.null(info$cadd)) {
      fields.cadd<-names(fields[["cadd"]])
      info.cadd<-as.list(info$cadd[,1])
      names(info.cadd)<-rownames(info$cadd)
      l1c<-lapply(fields.cadd,function(x){
        return(tags$tr(class=ifelse(x == 'Cscore','highlight','normal'),tags$td(em(fields[["cadd"]][[x]]),align="right"),tags$td(strong(info.cadd[[x]]))))
      }
      )
    }
    l1s<-vector(mode="list")
    if (!is.null(info) && !is.null(info$suspect)) {
      fields.suspect<-names(fields[["suspect"]])
      info.suspect<-as.list(info$suspect[,1])
      names(info.suspect)<-rownames(info$suspect)
      l1s<-lapply(fields.suspect,function(x){
        return(tags$tr(class=ifelse(x == 'score','highlight','normal'),tags$td(em(fields[["suspect"]][[x]]),align="right"),tags$td(strong(info.suspect[[x]]))))
      }
      )
    }
    ###>>> REVEL and GAVIN added
    # below l1r also added to l1

    ###>>> tool tip
    revelButton <- HTML("<button id=\"revelBtn\" type=\"button\" class=\"btn btn-default fa fa-info\" style=\"padding:0.8px 5px; font-size:80%; border-radius: 60%; color: #4CAF50; border: 1px dashed #4CAF50\"></button>")
    revelButtonTip <- bsTooltip("revelBtn", "Higher score is reflecting greater likelihood that the variant is disease-causing. For UK diagnostic labs thresholds are: <0.4 benign (BP1), >0.7 pathogenic (PP3)",trigger="hover", placement="top")

    l1r<-vector(mode="list")
    if (!is.null(info) && !is.null(info$revel)) {
      l1r <-  tags$tr(class='highlight',tags$td(em("Revel [0-1]: "), revelButton ,align="right"), revelButtonTip, tags$td(strong(info$revel[1,])))
    }

    l1g<-vector(mode="list")
    if (!is.null(info) && !is.null(info$gavin)) {
      l1g <-  tags$tr(class='highlight',tags$td(em("Gavin class: "),align="right"),tags$td(strong(info$gavin[,1])))
    }else{
      l1g <-  tags$tr(class='highlight',tags$td(em("Gavin class: "),align="right"),tags$td(strong("No result for this genomic seq. change")))
    }

    lbd<-vector(mode="list")
    if (!is.null(info) && !is.null(info$tp53_fut)) {
      lbd <-  tags$tr(class='highlight',tags$td(em("BayesDel: "),align="right"),tags$td(strong(round(as.numeric(info$tp53_fut$BayesDel),3))))
    }
    ###

    ###>>>
    missenseButton <- HTML("<button id=\"missenseBtn\" type=\"button\" class=\"btn btn-default action-button fa fa-info\" style=\"padding:0.8px 5px; font-size:80%; border-radius: 50%; color: #4CAF50; border: 1px dashed #4CAF50\"></button>")

    l1<-list(h4('Missense/nonsense predictions',missenseButton),
             tags$table(border="0",cellspacing="5",l1a,l1m,l1p,l1c,l1s,l1r,lbd,l1g),

             ###>>>
             bsTooltip("missenseBtn", "here you see some tips here you see some tips here you see some tips here you see some tips here you see some tips here",trigger="click")
    )
    l2<-vector(mode="list")
    if (!is.null(info) && !is.null(info$swissprot)) {
      l2<-list(h4('Swissprot features'),
               tags$table(border="0",cellspacing="5",
                          tags$tr(tags$td(em('Type')),tags$td(em('Description'))),
                          lapply(1:dim(info$swissprot)[1],function(x){
                            return(tags$tr(tags$td(info$swissprot[x,1]),tags$td(info$swissprot[x,2])))
                          }
                          )
               )
      )
    }

    spliceButton  <- HTML("<button id=\"spliceBtn\" type=\"button\" class=\"btn btn-default fa fa-info\" style=\"padding:0.8px 5px; font-size:80%; border-radius: 60%; color: #4CAF50; border: 1px dashed #4CAF50\"></button>")
    spliceButtonTip <- bsTooltip("spliceBtn", "PP3 for splicing is awarded if MaxEntScan >= 15% difference AND SSF >= 5% difference",trigger="hover")

    l2<-c(l2,list(h4('Conservation'),tags$table(border="0",cellspacing="5",l2c),
                  h4('Splicing predictions',spliceButton),spliceButtonTip,tags$table(border="0",cellspacing="5",l2s)
    )
    )
    return(as.character(tagList(fluidRow(column(width=8,l1),column(width=4,l2)))))
  })




  # Output for the Frequency panel
  output$controlfrequencyPanel <- renderUI ({

    ###>>>
    ### to arrange text size, h4 changed with  h6
    ### to have more space we changed the align arguments
    mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)

    info<-cigmainfo()
    info(loggerServer,"in controlfrequencyPanel")
    l<-vector(mode="list")


    if (!is.null(info) && !is.null(info$tgp)) {
      fields.tgp<-names(fields[["tgp"]])
      info.tgp<-as.list(info$tgp[,1])
      names(info.tgp)<-rownames(info$tgp)
      fields.tgp_pop<-names(fields[["tgp_pop"]])
      l<-c(l,list(h5(a('1000 Genome project (TGP / 1KG)',

                                 ###>>>
                                 #                       href=paste0('http://browser.1000genomes.org/Homo_sapiens/Variation/Population?db=core;v=',info$rsid,';vdb=variation'),
                                 href=paste0("http://phase3browser.1000genomes.org/Homo_sapiens/Variation/Population?db=core;r=",info$summary["hg19_chr",],":",info$summary["hg19_pos",],"-",info$summary["hg19_pos",],";v=",info$summary["rsID",],";vdb=variation;vf=9043431"),

                                 ###
                                 target='_blank')),
                            tags$table(border="0",cellspacing="5",
                                       lapply(fields.tgp,function(x){
                                         v<-info.tgp[[x]]
                                         if (grepl('AF',x)) v<-freqprint(v)
                                         return(tags$tr(tags$td(em(fields[["tgp"]][[x]]),align="right"),tags$td(strong(v), align="left")))
                                       }
                                       )
                            ),
                            tags$table(border="1",
                                       lapply(fields.tgp_pop,function(x){
                                         return(tags$tr(tags$td(em(fields[["tgp_pop"]][[x]]),align="right"),
                                                        lapply(info$tgp_pop[x,],function(y,fld=x){
                                                          v<-y
                                                          if (grepl('AF',fld)) v<-freqprint(v)
                                                          if (fld=='pop'){
                                                            return(tags$td(strong(v),align="left",title=tp[[v]]))
                                                          }else{
                                                            return(tags$td(strong(v),align="left"))
                                                          }
                                                        }
                                                        )
                                         )
                                         )
                                       }
                                       )
                            )
      )
      )
    }
    if (!is.null(info) && !is.null(info$tgp_phase3)) {
      fields.tgp_phase3<-names(fields[["tgp_phase3"]])
      info.tgp_phase3<-as.list(info$tgp_phase3[,1])
      names(info.tgp_phase3)<-rownames(info$tgp_phase3)
      fields.tgp_phase3_pop<-names(fields[["tgp_phase3_pop"]])
      ###>>>
      #      l<-c(l,list(h5('1000 Genome project - phase3 (TGP-phase3 / 2.5KG)'),

      l<-c(l,list(h5(a('1000 Genome project - phase3 (TGP-phase3 / 2.5KG)',
                       href=paste0("http://phase3browser.1000genomes.org/Homo_sapiens/Variation/Population?db=core;r=",info$summary["hg19_chr",],":",info$summary["hg19_pos",],"-",info$summary["hg19_pos",],";v=",info$summary["rsID",],";vdb=variation;vf=9043431"),
                       target='_blank')),

                  ####
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.tgp_phase3,function(x){
                               v<-info.tgp_phase3[[x]]
                               if (grepl('AF',x)) v<-freqprint(v)
                               return(tags$tr(tags$td(em(fields[["tgp_phase3"]][[x]]),align="right"),tags$td(strong(v))))
                             }
                             )
                  ),
                  tags$table(border="1",
                             lapply(fields.tgp_phase3_pop,function(x){
                               return(tags$tr(tags$td(em(fields[["tgp_phase3_pop"]][[x]]),align="right"),
                                              lapply(info$tgp_phase3_pop[x,],function(y,fld=x){
                                                v<-y
                                                if (grepl('AF',fld)) v<-freqprint(v)
                                                if (fld=='pop'){
                                                  return(tags$td(strong(v),align="left",title=tp[[v]]))
                                                }else{
                                                  return(tags$td(strong(v),align="left"))
                                                }
                                              }
                                              )
                               )
                               )
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$evs)) {
      fields.evs<-names(fields[["esp"]])
      info.evs<-as.list(info$evs[,1])
      names(info.evs)<-rownames(info$evs)
      info.main<-as.list(info$summary[,1])
      names(info.main)<-rownames(info$summary)
      l<-c(l,list(h5(a('Exome Sequencing Project (ESP) / Exome Variant Server (EVS)',
                       href=paste0('http://evs.gs.washington.edu/EVS/PopStatsServlet?searchBy=chromosome&chromosome=',
                                   info.main[["hg19_chr"]],'&chromoStart=',info.main[["hg19_pos"]],
                                   '&chromoEnd=',info.main[["hg19_pos"]],'&x=0&y=0'),target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.evs,function(x){
                               v<-info.evs[[x]]
                               if (grepl('MAF',x)) v<-freqprint(v)
                               return(tags$tr(tags$td(em(fields[["esp"]][[x]]),align="right"),tags$td(strong(v), align="left")))
                             }
                             )
                  )
      )
      )
    }




    ###>>> gnomAD data
    if(1==1){

      # dataset <- list()
      #
      # for(gp in c("","_count","_no","_homo")){
      #   toselect <- paste0("AFR",gp, ", AMR",gp, ", ASJ",gp, ", EAS",gp, ", FIN",gp, ", NFE",gp, ", OTH",gp, ", SAS",gp, ", Total", gp)


      #Cancer GnomAD v2.1.1

        myquery <- paste(sep = "",collapse = "", "select * from gnomadtest where gene='",info$summary["gene",],"' and hgvs_cdna='",info$summary["hgvs_cdna",],"' and transcript='",info$summary["transcript",],"'")
        rs<-dbSendQuery(mydb,myquery)
        dataset<-as.list(fetch(rs,-1))
        dataset[is.na(dataset) | dataset=="NA" | dataset=="character(0)"|dataset=="numeric(0)"] <- 0

        Allele=dataset[names(dataset) %>% endsWith("no")]  %>% as.data.frame() %>% t()
        AlleleCount=dataset[names(dataset) %>% endsWith("count")]    %>% as.data.frame() %>% t()
        Homozygotes=dataset[names(dataset) %>% endsWith("homo")]  %>% as.data.frame() %>% t()
        gnomad_cancer=as.data.frame(cbind(AlleleCount,Homozygotes,Allele))
        rownames(gnomad_cancer)=rownames(Allele) %>% substr(1,3)
        gnomad_cancer=gnomad_cancer[order(row.names(gnomad_cancer)),]
        fields.gnomad2 <- list("AFR"="African","AMR"="Latino","ASJ"="Ashkenazi Jewish","EAS"="East Asian","FIN"="European (Finnish)","NFE"="European (Non-Finnish)","OTH"="Other","SAS"="South Asian","Tot"="Total")
        rownames(gnomad_cancer)=lapply(rownames(gnomad_cancer),function(x){fields.gnomad2[[x]]})
        colnames(gnomad_cancer)=c('Allele Count',	'# Homozygotes',	'# Alleles sequenced')

        # gnomad_cancer$"# Total" <- c("12487","17720","5185","9977","12562","64603","3614","15308","141456")
        gnomad_cancer$"# WES in dataset" <- c("8128","17296","5040","9197","10824","56885","3070","15308","125748")
        # gnomad_cancer$"# WGS" <- c("4359","424","145","780","1738","7718","544","0","15708")

        #Non-Cancer GnomAD v2.1.1
          myquery <- paste(sep = "",collapse = "", "select * from gnomadnoncancer where gene='",info$summary["gene",],"' and hgvs_cdna='",info$summary["hgvs_cdna",],"' and transcript='",info$summary["transcript",],"'")
          rs<-dbSendQuery(mydb,myquery)
          dataset<-as.list(fetch(rs,-1))
          dataset[is.na(dataset) | dataset=="NA" | dataset=="character(0)"|dataset=="numeric(0)"] <- 0

          Allele=dataset[names(dataset) %>% endsWith("no")]  %>% as.data.frame() %>% t()
          AlleleCount=dataset[names(dataset) %>% endsWith("count")]    %>% as.data.frame() %>% t()
          Homozygotes=dataset[names(dataset) %>% endsWith("homo")]  %>% as.data.frame() %>% t()
          gnomad_NC=as.data.frame(cbind(AlleleCount,Homozygotes,Allele))
          rownames(gnomad_NC)=rownames(Allele) %>% substr(1,3)
          gnomad_NC=gnomad_NC[order(row.names(gnomad_NC)),]
          fields.gnomad2 <- list("AFR"="African","AMR"="Latino","ASJ"="Ashkenazi Jewish","EAS"="East Asian","FIN"="European (Finnish)","NFE"="European (Non-Finnish)","OTH"="Other","SAS"="South Asian","Tot"="Total")
          rownames(gnomad_NC)=lapply(rownames(gnomad_NC),function(x){fields.gnomad2[[x]]})
          colnames(gnomad_NC)=c('Allele Count',	'# Homozygotes',	'# Alleles sequenced')

          # gnomad_NC$"# Total" <- c("12487","17720","5185","9977","12562","64603","3614","15308","141456")
         gnomad_NC$"# WES in dataset" <- c("7451","17130","4786","8846","10816","51377","2810","15263","118479")
          # gnomad_NC$"# WGS" <- c("4359","424","145","780","1738","7718","544","0","15708")
          #
          #

          #Non-Cancer(Female only) GnomAD v2.1.1
          myquery <- paste(sep = "",collapse = "", "select * from gnomadnoncancer_female where gene='",info$summary["gene",],"' and hgvs_cdna='",info$summary["hgvs_cdna",],"' and transcript='",info$summary["transcript",],"'")
          rs<-dbSendQuery(mydb,myquery)
          dataset<-as.list(fetch(rs,-1))
          dataset[is.na(dataset) | dataset=="NA" | dataset=="character(0)"|dataset=="numeric(0)"] <- 0

          Allele=dataset[names(dataset) %>% endsWith("no")]  %>% as.data.frame() %>% t()
          AlleleCount=dataset[names(dataset) %>% endsWith("count")]    %>% as.data.frame() %>% t()
          Homozygotes=dataset[names(dataset) %>% endsWith("homo")]  %>% as.data.frame() %>% t()
          gnomad_NC_female=as.data.frame(cbind(AlleleCount,Homozygotes,Allele))
          rownames(gnomad_NC_female)=rownames(Allele) %>% substr(1,3)
          gnomad_NC_female=gnomad_NC_female[order(row.names(gnomad_NC_female)),]
          fields.gnomad2 <- list("AFR"="African","AMR"="Latino","ASJ"="Ashkenazi Jewish","EAS"="East Asian","FIN"="European (Finnish)","NFE"="European (Non-Finnish)","OTH"="Other","SAS"="South Asian","Tot"="Total")
          rownames(gnomad_NC_female)=lapply(rownames(gnomad_NC_female),function(x){fields.gnomad2[[x]]})
          colnames(gnomad_NC_female)=c('Allele Count',	'# Homozygotes',	'# Alleles sequenced')

          # gnomad_NC_female$"# Total" <- c("12487","17720","5185","9977","12562","64603","3614","15308","141456")
         #gnomad_NC_female$"# WES in dataset" <- c("8128","17296","5040","9197","10824","56885","3070","15308","125748")
          # gnomad_NC_female$"# WGS" <- c("4359","424","145","780","1738","7718","544","0","15708")

        # output$gnomad_frequency<-renderDT(if (input$gnomad_cancerdatatoshow=='gnomAD v2.1.1')
          # {DT::datatable(gnomad_cancer, options = list(bSort=FALSE,bPaginate=FALSE,dom='t',
           #                               columnDefs = list(list(className = 'dt-center', targets = 0:4)))
            #                                ) %>%formatStyle(
             #                              '# WES in dataset',
              ###                             `border-left` = "solid 2px")}

                 #                          else if (input$gnomad_cancerdatatoshow=='gnomAD v2.1.1 (non-cancer)')
                  #                           {DT::datatable(gnomad_NC, options = list(bSort=FALSE,bPaginate=FALSE,dom='t',
                   ##                                     columnDefs = list(list(className = 'dt-center',  targets = 0:4)))
                     #                        ) %>%formatStyle(

                       #                      '# WES in dataset',
                      #                       `border-left` = "solid 2px")


                        #                     }
                         #                  else {DT::datatable(gnomad_NC_female, options = list(bSort=FALSE,bPaginate=FALSE,dom='t',
                          #                                                                      columnDefs = list(list(className = 'dt-center', targets = 0:4)))
                           #                ) %>%formatStyle(
                            #                 '# WES in dataset',
                             #                `border-left` = "solid 2px")}
          # )

 gnomad_list <- list('a'=gnomad_cancer,'b'=gnomad_NC,'c'=gnomad_NC_female)
          observeEvent(input$gnomad_cancerdatatoshow,{
            # observe element change and render table
            output$gnomad_frequency <- renderTable(
           gnomad_list[input$gnomad_cancerdatatoshow]
            ,na = NA, digits = 0,align = 'c', width = "70%",rownames =TRUE,colnames = T, striped=T, hover=T, bordered=T, spacing = "s")
          })
      lgnomad <- list(h5(a('gnomAD Population Frequencies (v2.1 - Jan 2019)',
                           href=ifelse(dataset$link!="0", paste0('http://gnomad.broadinstitute.org/variant/',dataset$link), paste0('http://gnomad.broadinstitute.org/gene/',info$summary["gene",])),
                           # href=ifelse(gnomadlink!="character(0)", paste0('http://gnomad.broadinstitute.org/variant/',gnomadlink), paste0('http://gnomad.broadinstitute.org/gene/',info$summary["gene",])),
                           target='_blank')))



      dbDisconnect(mydb)
      return(tagList(
  fluidRow(column(12,lgnomad,selectInput(inputId = "gnomad_cancerdatatoshow",
                                      label = "Choose a dataset:",
                                      choices = c("gnomAD v2.1.1"="a", "gnomAD v2.1.1 (non-cancer)"="b",
                                                  "gnomAD v2.1.1 (non-cancer female)"="c"),
                                      selected = "b"),
                      style="margin-left:5px;margin-right:5px;font-size:15px;",
               div( tags$head(tags$style("#gnomad_frequency table tr td{background-color: white;word-wrap: break-word;}", media="screen", type="text/css"))),
                      
                      div(tableOutput("gnomad_frequency"),align="center"),l)
               
               )))
     # fluidRow(column(12,lgnomad,selectInput(inputId = "gnomad_cancerdatatoshow",
      #                                label = "Choose a dataset:",
       #                               choices = c("gnomAD v2.1.1", "gnomAD v2.1.1 (non-cancer)",
        #                                          "gnomAD v2.1.1 (non-cancer female)"),
         #                             selected = "gnomAD v2.1.1 (non-cancer)"),
          #            style="margin-left:5px;margin-right:5px;font-size:15px;",
          #            div(DT::dataTableOutput("gnomad_frequency"),align="center"),l)
           #    )))
    #return(as.character(tagList(l)))
    }
    })

  output$casefrequencyPanel <- renderUI ({

    ###>>>
    ### to arrange text size, h4 changed with  h6
    ### to have more space we changed the align arguments
    mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)

    info<-cigmainfo()
    info(loggerServer,"in casefrequencyPanel")
    case_freq_l<-vector(mode="list")


    if(1==1){
      ###>>> Total probands (fullscreen) ; Total probands tested (total) olarak degisti, "pvalue" cikarildi

      fields.phe2 <- list("BIRq3"="Birmingham", "BRIrvj"="Bristol", "CAMrgt"="Cambridge",
                          "GUYrj1"="Guy's", "LEEDrr8"="Leeds", "NEWrtd"="Newcastle", "MANrw3"="Manchester", "NOTrx1"="Nottingham",
                          "SALrnz"="Salisbury", "SHEFrcu"="Sheffield","scrCount"="Total proband count", "all_tested"="Total probands tested",
                          "WhiteCount"="White eth.(use as NFE): total proband count", "white_tested"="White eth.(use as NFE):total probands tested")


      toselect <- paste0(names(fields.phe2), collapse = ", ")
      myquery <- paste(sep = "",collapse = "", "select ",toselect  ," from PHE LEFT JOIN PHE_probands ON PHE.gene = PHE_probands.gene where PHE.gene='",info$summary["gene",],"' and transcript='",info$summary["transcript",],"' and hgvs_cdna='", info$summary["hgvs_cdna",],"'")
      rs<-dbSendQuery(mydb,myquery)
      phest <- as.list(fetch(rs,-1))
      t= as.data.frame(phest) %>% t()

      if(!values$logged){
        fields.phe2 <- list("scrCount"="Total proband count", "all_tested"="Total probands tested"
)        ###>>> Total probands (fullscreen) ; Total probands tested (total) olarak degisti, "pvalue" cikarildi
       t=as.data.frame(phest)[c('scrCount','all_tested')] %>% t()
         }
      rownames(t)= fields.phe2
      output$protable2 =renderTable(na = NA, digits = 0,align = 'c',expr = t, width = "70%",rownames =TRUE,colnames = F, striped=T, hover=T, bordered=T, spacing = "s")
      # myquery <- paste(sep = "",collapse = "", "select * from PHE_probands where gene='",info$summary["gene",],"'")
      # rs<-dbSendQuery(mydb,myquery)
      # probands <- as.data.frame(fetch(rs,-1)) %>% select(-'gene') %>% t()
      #
      # output$protable =renderTable(probands,rownames =TRUE,
      #                              colnames = F)
      #



    ###>>> the line below has 1==0
    if (!is.null(info) && !is.null(info$icr) && 1==0) {
      fields.icr<-names(fields[["icr"]])
      info.icr<-as.list(info$icr[,1])
      names(info.icr)<-rownames(info$icr)
      case_freq_l<-c(case_freq_l,list(h5('1958 Birth Cohort (ICR)'),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.icr,function(x){
                               return(tags$tr(tags$td(em(fields[["icr"]][[x]]),align="right"),tags$td(strong(info.icr[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$bic)) {
      fields.bic<-names(fields[["bic"]])
      info.bic<-as.list(info$bic[,1])
      names(info.bic)<-rownames(info$bic)
      ###>>>
      case_freq_l<-c(case_freq_l,list(h5(a('BIC (2015)',href='https://research.nhgri.nih.gov/bic/',target='_blank')),
                  #
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.bic,function(x){
                               return(tags$tr(tags$td(em(fields[["bic"]][[x]]),align="right"),tags$td(strong(info.bic[[x]]))))
                             }
                             )
                  )
      )
      )
    }

    ###>>> edit for LOVD3
    if (!is.null(info) && !is.null(info$lovd3)) {
      # fields.lovd<-names(fields[["lovd"]])
      # fields.lovd<-fields.lovd[-which(fields.lovd=='insight_class')]
      # info.lovd<-as.list(info$lovd[,1])
      # names(info.lovd)<-rownames(info$lovd)
      gene<-info$summary['gene',1]
      variant<-info$summary['hgvs_cdna',1]
      lovd_link<-'LOVD'
      if (grepl('BRCA',gene)){
        lovd_link <- a(href=paste0('http://chromium.liacs.nl/LOVD2/cancer/variants.php?select_db=',gene,'&action=search_unique&search_Variant%2FDNA=',variant),'LOVD',target="_blank")
      }
      if (gene %in% c('MLH1','MSH2','MSH6','PMS2')){
        lovd_link <- a(href=paste0('http://chromium.liacs.nl/LOVD2/colon_cancer/variants.php?select_db=',gene,'&action=search_unique&search_Variant%2FDNA=',variant),'LOVD',target="_blank")
      }

      lovd_link <- a(href=paste0('https://databases.lovd.nl/shared/variants/in_gene?search_geneid=%3D"',gene,'"&search_VariantOnTranscript/DNA=%3D"',variant,'"'),'LOVD3',target="_blank")

      # l<-c(l,list(h5(lovd_link),
      #             tags$table(border="0",cellspacing="5",
      #                        lapply(fields.lovd,function(x){
      #                          return(tags$tr(tags$td(em(fields[["lovd"]][[x]]),align="right"),tags$td(strong(info.lovd[[x]]))))
      #                        }
      #                        )
      #             )
      # )
      # )

      l1lovd3 <-  list(h5(lovd_link),tags$tr(class='highlight',tags$td(em("Number of times reported in LOVD3 (July 2019): "),align="right"),tags$td(strong(info$lovd3[,1]))))
      case_freq_l<-c(case_freq_l,l1lovd3)
    }

    if (!is.null(info) && !is.null(info$dmudb)) {
      fields.dmudb<-names(fields[["dmudb"]])
      info.dmudb<-as.list(info$dmudb[,1])
      names(info.dmudb)<-rownames(info$dmudb)
      case_freq_l<-c(case_freq_l,list(h5('DMuDB'),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.dmudb,function(x){
                               return(tags$tr(tags$td(em(fields[["dmudb"]][[x]]),align="right"),tags$td(strong(info.dmudb[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$umd)) {
      fields.umd<-names(fields[["umd"]])
      info.umd<-as.list(info$umd[,1])
      names(info.umd)<-rownames(info$umd)
      case_freq_l<-c(case_freq_l,list(h5('UMD'),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.umd,function(x){
                               return(tags$tr(tags$td(em(fields[["umd"]][[x]]),align="right"),tags$td(strong(info.umd[[x]]))))
                             }
                             )
                  )
      )
      )
    }

    dbDisconnect(mydb)
    return(tagList(

      fluidRow(column(12,
                      h5('UK diagnostic labs (Release date: July 2020)'),
                      div( tags$head(tags$style("#protable2 table tr td{background-color: white; }", media="screen", type="text/css")),
                           tableOutput("protable2"),align="center")
                      ,case_freq_l))))
    #return(as.character(tagList(l)))
  }
    })


  # Output for the Genetic Epidemiology panel
  output$geneticepiPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in geneticepiPanel")
    l<-vector(mode="list")
    if (!is.null(info) && !is.null(info$easton)) {
      fields.easton<-names(fields[["easton"]])
      info.easton<-as.list(info$easton[,1])
      names(info.easton)<-rownames(info$easton)
      l<-c(l,list(h4(a('Easton et al, Am J Hum Genet. 2007',href='http://www.ncbi.nlm.nih.gov/pubmed/17924331',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.easton,function(x){
                               return(tags$tr(tags$td(em(fields[["easton"]][[x]]),align="right"),tags$td(strong(info.easton[[x]]))))
                             }
                             )
                  )
      )
      )
    }

    if (!is.null(info) && !is.null(info$lindor)) {
      fields.lindor<-names(fields[["lindor"]])
      info.lindor<-as.list(info$lindor[,1])
      names(info.lindor)<-rownames(info$lindor)
      l<-c(l,list(h4(a('Lindor et al, Hum Mutat. 2011',href='http://www.ncbi.nlm.nih.gov/pubmed/21990134',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.lindor,function(x){
                               v<-info.lindor[[x]]
                               if (grepl('reference',x) && info$lindor['pubmed',1]!='') v<-a(v,href=paste0('http://www.ncbi.nlm.nih.gov/pubmed/',info$lindor['pubmed',1]),target='_blank')
                               return(tags$tr(tags$td(em(fields[["lindor"]][[x]]),align="right"),tags$td(strong(v))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$iarc)) {
      fields.iarc<-names(fields[["iarc"]])
      info.iarc<-as.list(info$iarc[,1])
      names(info.iarc)<-rownames(info$iarc)
      l<-c(l,list(h4(a('IARC: Vallee et al, Hum Mutat. 2012',href='http://www.ncbi.nlm.nih.gov/pubmed/21990165',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.iarc,function(x){
                               return(tags$tr(tags$td(em(fields[["iarc"]][[x]]),align="right"),tags$td(strong(info.iarc[[x]]))))
                             }
                             )
                  )
      )
      )
    }

    if (!is.null(info) && !is.null(info[['lovd']]) && info[['lovd']]['insight_class',]!="") {
      gene<-info$summary['gene',1]
      variant<-info$summary['hgvs_cdna',1]
      l<-c(l,list(h4(a('InSiGHT: Thompson et al, Nat Genet. 2014',href='http://www.ncbi.nlm.nih.gov/pubmed/24362816',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             tags$tr(tags$td(em('InSiGHT class:'),align="right"),tags$td(a(strong(info$lovd['insight_class',]),href=paste0('https://googledrive.com/host/0B8HVsr5izQxJUi1XTzEtWFlRc00/index.html?gene=',gene,'&protein=&variant=',variant),target='_blank')))
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info[["parsons_2019"]])) {
      parsons_2019<- info[["parsons_2019"]]
      parsons_2019 <- parsons_2019[,c(1:5)]
      colnames(parsons_2019) <- c("Combined LR (Odds for Causality):","Prior probability of pathogenicity:","Posterior probability:", "IARC class", "Comment:")


      l<-c(l,list(h4(a('Parsons et al, Human Mutation 2019',href='https://pubmed.ncbi.nlm.nih.gov/31131967/',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(colnames(parsons_2019),function(x){
                               return(tags$tr(tags$td(em(x),align="right"),tags$td(strong(parsons_2019[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info[["bouwman_2020_epi"]])) {
      bouwman_2020_epi<- info[["bouwman_2020_epi"]]
      colnames(bouwman_2020_epi) <- c("Combined LR:","Posterior probability:","Bouwman et al., 2020 and Parsons et al., 2019 combined LR:", "Bouwman et al., 2020 and Parsons et al., 2019 combined posterior probability:")


      l<-c(l,list(h4(a('Bouwman et al, Clin Cancer Res 2020',href='https://pubmed.ncbi.nlm.nih.gov/32546644/', target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(colnames(bouwman_2020_epi),function(x){
                               return(tags$tr(tags$td(em(x),align="right"),tags$td(strong(bouwman_2020_epi[[x]]))))
                             }
                             )
                  )
      )
      )
    }



    return(as.character(tagList(l)))
  })

  # Output for the Splicing Analysis panel
  output$splicingPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in splicingPanel")
    l<-vector(mode="list")
    if (!is.null(info) && !is.null(info$houdayer)) {
      fields.houdayer<-names(fields[["houdayer"]])
      info.houdayer<-as.list(info$houdayer[,1])
      names(info.houdayer)<-rownames(info$houdayer)
      l<-c(l,list(h4(a('Houdayer et al, Hum Mutat. 2012',href='http://www.ncbi.nlm.nih.gov/pubmed/22505045',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.houdayer,function(x){
                               return(tags$tr(tags$td(em(fields[["houdayer"]][[x]]),align="right"),tags$td(strong(info.houdayer[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$walker)) {
      fields.walker<-names(fields[["walker"]])
      l<-c(l,list(h4(a('Walker et al, Hum Mutat. 2013',href='http://www.ncbi.nlm.nih.gov/pubmed/23893897',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.walker,function(x){
                               return(tags$tr(tags$td(em(fields[["walker"]][[x]]),align="right"),
                                              lapply(info$walker[x,],function(y,fld=x){
                                                v<-y
                                                if (grepl('pubmed',fld)) v<-a(v,href=paste0('http://www.ncbi.nlm.nih.gov/pubmed/',v),target='_blank')
                                                return(tags$td(strong(v),align="right"))
                                              }
                                              )
                               )
                               )
                             }
                             )
                  )
      )
      )
    }

    ###>>> Wessex data
    if (!is.null(info) && !is.null(info$wessex)){
      if(values$logged){

        dataset <- as.data.frame(t(do.call("cbind", as.list(info$wessex))),stringsAsFactors=F)
        #fields.wessex <- colnames(dataset) <- as.character(dataset["specimen",])
        #dataset <- cbind(c("Specimen","Splice abnormality","Sample type"), dataset)
        dataset <- cbind("headers"=c("Specimen","Splice abnormality","Sample type"), dataset, stringsAsFactors=F)
        fields.wessex <- colnames(dataset) <- c("headers",as.character(dataset["specimen",-1]))
        dataset <- lapply(dataset, as.list)

        lwsx <- list(h5('UK diagnostic labs (Wessex data)'),
                     tags$table(border="1", cellspacing="5",
                                lapply(fields.wessex,function(x){
                                  v <- dataset[[x]]
                                  return(tags$tr(tags$td(v[[1]],align="center"),
                                                 tags$td(v[[2]],align="center"),
                                                 tags$td(v[[3]],align="center")))
                                }
                                )
                     )
        )
      }else{
        lwsx <- list(h5('UK diagnostic labs'),em("Data available"))
      }
    }else{
      lwsx <- list()
    }
    ###

    ###>>> Findlay appears in 2 tabs "Functional analysis" and "Splicing analysis"

    if (!is.null(info) && !is.null(info$findlay) && !is.na(info$findlay[1,1])) {

      findlayButton  <- HTML("<button id=\"findlayBtn\" type=\"button\" class=\"btn btn-default fa fa-info\" style=\"padding:0.8px 5px; font-size:80%; border-radius: 60%; color: #4CAF50; border: 1px dashed #4CAF50\"></button>")
      findlayButtonTip <- bsTooltip("findlayBtn", "RNA score is a mRNA expression score derived from (log2) normalized SNV frequencies in cDNA to their frequencies in gDNA. RNA score between -2 and -3 indicate around 75% reduction in mRNA expression. RNA scores < -3 equates to a high (>75%) reduction in mRNA expression.", trigger="hover", placement="top")

      findlayButton3  <- HTML("<button id=\"findlayBtn3\" type=\"button\" class=\"btn btn-default fa fa-info\" style=\"padding:0.8px 5px; font-size:80%; border-radius: 60%; color: #4CAF50; border: 1px dashed #4CAF50\"></button>")
      findlayButtonTip3 <- bsTooltip("findlayBtn3", "Functional classifications were made by setting thresholds for Pnf (function probability):  Pnf > 0.99: non-functional, 0.01< Pnf < 0.99: intermediate, Pnf <0.01: functional", trigger="hover",  placement="top")


      Findlay <- info$findlay
      Findlay$score <- round(as.numeric(Findlay$score),3)
      Findlay$probability <- round(as.numeric(Findlay$probability),3)
      Findlay$meanRNA <- round(as.numeric(Findlay$meanRNA),3)
      colnames(Findlay) <- c("Function score:","Functional classifications:","Function probability:", "Mean RNA score:")
      Findlay <- Findlay[,c(1,3,2,4)]

      l<-c(l,list(h4(findlayButton3,a('Saturation genome editing BRCA1 haploid cell assay (Findlay 2018)',href='https://www.nature.com/articles/s41586-018-0461-z',target='_blank'),findlayButton),findlayButtonTip,findlayButtonTip3,
                  tags$table(border="0",cellspacing="5",
                             lapply(colnames(Findlay),function(x){
                               return(tags$tr(tags$td(em(x),align="right"),tags$td(strong(Findlay[[x]]))))
                             }
                             )
                  )
      )
      )

    }
    ###

    return(as.character(tagList(fluidRow(column(width=6,l), column(width=5, offset=1, lwsx)))))
    #return(as.character(tagList(l)))
  })

  # Output for the Functional Analysis panel
  output$functionalPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in functionalPanel")
    l<-vector(mode="list")

    if (!is.null(info) && !is.null(info[["guidugli"]])) {
      fields.guidugli<-names(fields[["guidugli"]])
      info.guidugli<-as.list(info$guidugli[,1])
      names(info.guidugli)<-rownames(info$guidugli)
      l<-c(l,list(h4(a('Guidugli et al, Cancer Res. 2013',href='http://www.ncbi.nlm.nih.gov/pubmed/23108138',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.guidugli,function(x){
                               return(tags$tr(tags$td(em(fields[["guidugli"]][[x]]),align="right"),tags$td(strong(info.guidugli[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$drost)) {
      fields.drost<-names(fields[["drost"]])
      info.drost<-as.list(info$drost[,1])
      names(info.drost)<-rownames(info$drost)
      if (grepl('MSH',info$summary['gene',1])) {
        title<-a('Drost et al, Hum Mutat. 2012',href='http://www.ncbi.nlm.nih.gov/pubmed/22102614',target='_blank')
      }else{
        title<-a('Drost et al, Hum Mutat. 2010',href='http://www.ncbi.nlm.nih.gov/pubmed/20020535',target='_blank')
      }
      l<-c(l,list(h4(title),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.drost,function(x){
                               return(tags$tr(tags$td(em(fields[["drost"]][[x]]),align="right"),
                                              lapply(info.drost[[x]],function(y,fld=x){
                                                v<-y
                                                if (grepl('reference|ficient',fld)) v<-gsub('([^(,]+)\\(PMID:(\\d+)\\)','<A HREF="http://www.ncbi.nlm.nih.gov/pubmed/\\2" target="_blank">\\1</A>',y,perl=T)
                                                return(tags$td(strong(HTML(v))))
                                              }
                                              )
                               )
                               )
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info[['bouwman']])) {
      fields.bouwman<-names(fields[["bouwman"]])
      info.bouwman<-as.list(info[["bouwman"]][,1])
      names(info.bouwman)<-rownames(info[["bouwman"]])
      l<-c(l,list(h4(a('Bouwman et al, Cancer Discov. 2013',href='http://www.ncbi.nlm.nih.gov/pubmed/23867111',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.bouwman,function(x){
                               return(tags$tr(tags$td(em(fields[["bouwman"]][[x]]),align="right"),tags$td(strong(info.bouwman[[x]]))))
                             }
                             )
                  )
      )
      )
    }




    #####>>> adding Findlay data

    #   info <- cigmainfo()
    #   query <- paste(sep="", collapse="","select score,class,probability from Findlay where gene='",info$summary["gene",],
    #         "' AND  hgvs_cdna='",info$summary["hgvs_cdna",],"'")
    #   mydb<-dbConnect(MySQL(),user='batch',password='cigma',dbname='cigma2',host=host)
    #   rs <- dbSendQuery(mydb,query)
    #   Findlay <- fetch(rs,-1)
    #   dbDisconnect(mydb)

    if (!is.null(info) && !is.null(info$findlay) && !is.na(info$findlay[1,1])) {

      findlayButton1  <- HTML("<button id=\"findlayBtn1\" type=\"button\" class=\"btn btn-default fa fa-info\" style=\"padding:0.8px 5px; font-size:80%; border-radius: 60%; color: #4CAF50; border: 1px dashed #4CAF50\"></button>")
      findlayButtonTip1 <- bsTooltip("findlayBtn1", "RNA score is a mRNA expression score derived from (log2) normalized SNV frequencies in cDNA to their frequencies in gDNA. RNA score between -2 and -3 indicate around 75% reduction in mRNA expression. RNA scores < -3 equates to a high (>75%) reduction in mRNA expression. ", trigger="hover",  placement="top")

      findlayButton2  <- HTML("<button id=\"findlayBtn2\" type=\"button\" class=\"btn btn-default fa fa-info\" style=\"padding:0.8px 5px; font-size:80%; border-radius: 60%; color: #4CAF50; border: 1px dashed #4CAF50\"></button>")
      findlayButtonTip2 <- bsTooltip("findlayBtn2", "Functional classifications were made by setting thresholds for Pnf (function probability):  Pnf > 0.99: non-functional, 0.01< Pnf < 0.99: intermediate, Pnf <0.01: functional", trigger="hover",  placement="top")


      Findlay <- info$findlay
      Findlay$score <- round(as.numeric(Findlay$score),3)
      Findlay$probability <- round(as.numeric(Findlay$probability),3)
      Findlay$meanRNA <- round(as.numeric(Findlay$meanRNA),3)
      colnames(Findlay) <- c("Function score:","Functional classifications:","Function probability:", "Mean RNA score:")
      Findlay <- Findlay[,c(1,3,2,4)]

      l<-c(l,list(h4(findlayButton2, a('Saturation genome editing BRCA1 haploid cell assay (Findlay 2018)',href='https://www.nature.com/articles/s41586-018-0461-z',target='_blank'),findlayButton1),findlayButtonTip1,findlayButtonTip2,
                  tags$table(border="0",cellspacing="5",
                             lapply(colnames(Findlay),function(x){
                               return(tags$tr(tags$td(em(x),align="right"),tags$td(strong(Findlay[[x]]))))
                             }
                             )
                  )
      )
      )
    }

    ###

    ###>>> adding guidugli2018 data

    if (!is.null(info) && !is.null(info$guidugli2018) && !is.na(info$guidugli2018[1,1])) {

      Guidugli2018 <- info$guidugli2018
      colnames(Guidugli2018) <- c("Probability deleterious (Guidugli):", "Classification (Guidugli):")

      l<-c(l,list(h4(a('Guidugli et al, Am J Hum Genet. 2018',href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5985401/',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(colnames(Guidugli2018),function(x){
                               return(tags$tr(tags$td(em(x),align="right"),tags$td(strong(Guidugli2018[[x]]))))
                             }
                             )
                  )
      )
      )
    }

    ###

    ###>>> adding TP53_SGE
    if (!is.null(info) && !is.null(info$tp53_sge)) {

      # c.1177G>T
      tp53sgeButton <- HTML("<button id=\"tp53sgeBtn\" type=\"button\" class=\"btn btn-default fa fa-info\" style=\"padding:0.8px 5px; font-size:80%; border-radius: 60%; color: #4CAF50; border: 1px dashed #4CAF50\"></button>")
      tp53sgeButtonTip <- bsTooltip("tp53sgeBtn", "IARC classification based on growth suppression assays in A549 human cells; DNE+LOF is p53WTNutlin3 Z-score >= 0.61 and Etoposide Z-score <= -0.21; noDNE+noLOF is p53WTNutlin3 Z-score < 0.61 and Etoposide Z-score > -0.21",trigger="hover", placement="top")

      l<-c(l,list(h4(a('Giacomelli et al, Nat Genet. 2018',href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6168352/',target='_blank'),tp53sgeButton),tp53sgeButtonTip,
                  tags$tr(class='highlight',tags$td(em("Combined model score [0:+Inf]:  "),align="right"),tags$td(strong(round(as.numeric(info$tp53_sge$Combined_Model),4)))),
                  br(),
                  tags$tr(class='highlight',tags$td(em("IARC classification:  "),align="right"),tags$td(strong(info$tp53_sge$ClassSGE)))
      )
      )

    }

    if (!is.null(info) && !is.null(info$tp53_kato)) {

      tp53katoButton <- HTML("<button id=\"tp53katoBtn\" type=\"button\" class=\"btn btn-default fa fa-info\" style=\"padding:0.8px 5px; font-size:80%; border-radius: 60%; color: #4CAF50; border: 1px dashed #4CAF50\"></button>")
      tp53katoButtonTip <- bsTooltip("tp53katoBtn", "Functional classification based on the overall transcriptional activity (TA) on 8 different promoters as measured in yeast assays by Kato et al(2003). For each mutant, the median of the 8 promoter-specific activities (expressed as percent of the wild-type protein) is calculated and missense mutations are classified as non-functional if the median is <=20%, partially functional if the median is >20% and <=75%, functional if the median is >75% and <=140%, and supertrans if the median is >140%.",trigger="hover", placement="top")

      l<-c(l,list(h4(a('Kato et al, PNAS, 2003',href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC166245/',target='_blank'),tp53katoButton),tp53katoButtonTip,
                  tags$tr(class='highlight',tags$td(em("Transactivation Class: "),align="right"),tags$td(strong(info$tp53_kato$TransactivationClass)))
      )
      )
    }


    if (!is.null(info) && !is.null(info$tp53_fut)) {

      l<-c(l,list(h4(a('Fortuno et al, Hum. Mutat., 2018',href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6043381/',target='_blank')),
                  tags$tr(class='highlight',tags$td(em("BayesDel: "),align="right"),tags$td(strong(round(as.numeric(info$tp53_fut$BayesDel),3)))),
                  br(),
                  tags$tr(class='highlight',tags$td(em("Suggested prediction: "),align="right"),tags$td(strong(info$tp53_fut$SuggestedPrediction)))
      )
      )
    }


    ###
    if (!is.null(info) && !is.null(info$bouwman_2020_func)) {
      bouwman_2020_func<- info[["bouwman_2020_func"]]
      bouwman_2020_func$Cisplatin_probability_deleterious <- signif(as.numeric(bouwman_2020_func$Cisplatin_probability_deleterious),digits = 3)
      bouwman_2020_func$Olaparib_probability_deleterious<- signif(as.numeric(bouwman_2020_func$Olaparib_probability_deleterious),digits = 3)
      bouwman_2020_func$DR_GFP_probability_deleterious <- signif(as.numeric(bouwman_2020_func$DR_GFP_probability_deleterious),digits =3)

      colnames(bouwman_2020_func) <- c("Cisplatin_probability_deleterious:","Cisplatin_prediction:","Olaparib_probability_deleterious:", "Olaparib_prediction:","DR_GFP_probability_deleterious:","DR_GFP_prediction:")
      l<-c(l,list(h4(a('Bouwman et al, Cancer Res. 2020',href='https://pubmed.ncbi.nlm.nih.gov/32546644/', target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(colnames(bouwman_2020_func),function(x){
                               return(tags$tr(tags$td(em(x),align="right"),tags$td(strong(bouwman_2020_func[[x]]))))
                             }
                             )
                  )
      )
      )
    }

    return(as.character(tagList(l)))
  })


  varNotes <- reactive({
    i <- list()
    mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
    splitgene<-unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))
    gene<-splitgene[1]
    transcript<-splitgene[2]
    variation<-gsub(' [(].*$','',input$variationSel,perl=T)
    rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                               "select investigator_name,
                                 clinical_cigma_class,
                                 concat(date_format(creat_date,'%d/%m/%Y'),' (',
                                 if(round(datediff(curdate(),creat_date)/7)<9,
                                 concat(round(datediff(curdate(),creat_date)/7),' weeks ago'),
                                 concat(round(datediff(curdate(),creat_date)/30.4166),' months ago')
                                 ),')') date,
                                 notes
                                 from dogma_batchlog b
                                 where b.gene='",gene,"' and b.variation='",variation,"'
                                 order by creat_date desc"));
    dataset<-fetch(rs,-1)
    if (!is.null(dataset) && dim(dataset)[1]>0) {
      i$notes <- as.data.frame(dataset,stringsAsFactors=F)
    }else{
      i$notes <- NULL
    }

    dirNotesOld <- paste0("/tmp/NOTES/",gene,"/", gsub("[^[:alnum:]\\_\\-]","",variation))
    if(dir.exists(dirNotesOld)){ i$downlink <- downloadLink("downloadadditionalData", label = "Download the additional notes") }else{ i$downlink <- ""}
    dbDisconnect(mydb)
    return(i)
  })

  output$notesPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in notesPanel")
    l<-list();
    notes_new <- NULL
    if (values$contentType=='browseTabs') {
      notes_new <- values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'notes_new']
      if (notes_new[1] != '') {
        nn <- unlist(str_split(notes_new[1],"\t"))
        ntext <- gsub('PMID[:]?(\\s+)?(\\d+)','<A HREF="http://www.ncbi.nlm.nih.gov/pubmed?term=\\2" target="_blank">PMID:\\2</A>',nn[4],perl=T)
        l<-c(l,
             list(div(class='highlight',
                      fluidRow(column(width=4,em('Provisioned class:')),column(width=8,strong(nn[2]))),
                      fluidRow(column(width=4,em('Date: ')),column(width=8,nn[3])),
                      fluidRow(column(width=4,em('Name: ')),column(width=8,nn[1]))),
                  fluidRow(div(class="details",style="margin-bottom:20px;",em("Notes:"),lapply(unlist(str_split(ntext,'<br>')),function(x){div(HTML(x))})))
             )
        )
      }
    }
    if (!is.null(info) && !is.null(info$notes) && !is.na(info$notes)){
      ###>>> serial save
      ### and show comments only to the logged in users
      ### and above line shows the comments if the user logged in values$logged

      info$notes <- t(as.matrix(gsub(" [|] ","<br>",info$notes)))

      #      l<-c(l,apply(info$notes,1,function(x){
      #        ntext <- gsub('PMID:(\\s+)?(\\d+)','<A HREF="http://www.ncbi.nlm.nih.gov/pubmed?term=\\2" target="_blank">PMID:\\2</A>',x[4],perl=T)
      #        list(div(class='highlight',
      #                 fluidRow(column(width=4,em('Curated classification: ')),column(width=8,strong(x[2]))),
      #                 fluidRow(column(width=4,em('Date: ')),column(width=8,x[3])),
      #                 fluidRow(column(width=4,em('Name: ')),column(width=8,x[1]))),
      #                 fluidRow(div(class="details",style="margin-bottom:20px",em('Notes (date; reviewer ID; affiliation; curated classification; note):'),lapply(unlist(str_split(ntext,'<br>')),function(x){div(HTML(x))})))
      #        )
      #      }
      #      )
      #      )

      # remove the orange part in the notes part
      if(1==0){
        l<-c(l,apply(info$notes,1,function(x){
          list(div(class='highlight',
                   fluidRow(column(width=4,em('Proposed class: ')),column(width=8,strong(x[2]))),
                   fluidRow(column(width=4,em('Date: ')),column(width=8,x[3]))
                   # fluidRow(column(width=4,em('Name: ')),column(width=8,x[1]))
          )
          )}))
      }
      if(values$logged){
        mynotes_tmp <- varNotes()
        if(values$notesupdate!=""){
          #  mynotes <- as.matrix(gsub(" [|] ","<br>",values$notesupdate)) #mynotes, info$notes ile degistirildi, x[4] x yapildi ntext satirinda
          mynotes_tmp <- varNotes()
          mynotes <- t(as.matrix(gsub(" [|] ","<br>", mynotes_tmp$notes[,4])))
          nvals <- values$notesupdate
        }else{
          mynotes <- as.matrix(info$notes[,4])
        }

        #    splitgene <- unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))
        #	gene <- splitgene[1]
        #	dirNotesOld <- paste0("/tmp/NOTES/",gene,"/",gsub("[[:punct:]]","",input$variationText))
        #	if(dir.exists(dirNotesOld) || values$dirNotes!=""){ values$downlink <- downloadLink("downloadData", label = "Download the additional notes") }else{values$downlink <- ""}

        l<-c(l,apply(mynotes,1,function(x){
          ntext <- gsub('PMID:(\\s+)?(\\d+)','<A HREF="http://www.ncbi.nlm.nih.gov/pubmed?term=\\2" target="_blank">PMID:\\2</A>',x,perl=T)
          fluidRow(div(class="details",style="margin:20px;",em('Notes (date; reviewer; affiliation; proposed classification; note): '),
                       list(div(HTML(paste0(mynotes_tmp$downlink)))), list(br()),# list(div(HTML(values$notesupdate))),
                       #				lapply(unlist(str_split(ntext,'<br>')),function(x){div(HTML(x))})
                       ###>>> Make Date, reviewer, classification bold
                       lapply(unlist(strsplit(ntext,'<br>')),function(x){  w <-strsplit(x,";")[[1]]; w<-paste0("<strong>",paste0(w[c(1,2)], collapse=";"),";</strong>", w[3], "; ",paste0("<strong>",w[4],"; </strong>"), w[5], collapse=""); div(HTML(w))})))
        }))
      }else{

        visibleComments <- span(title="Visible comments",
                                verbatimTextOutput("visibleComments"))

        l <- c(l, list(fluidRow(visibleComments)))
      }
    }
    #    l<-c(l,list(div(class='highlight',
    #                    fluidRow(column(width=4,em('Automated class:')),column(width=8,strong(info$baseline_class))),
    #                    fluidRow(column(width=4,em('Justification:')),column(width=8,strong(info$justification)))
    #    )
    #    )
    #    )
    return(as.character(tagList(l)))
  })




  ###################
  # Downloaded data #
  ###################
  output$downloadVariantReport = downloadHandler(
    filename =  function() {info<-cigmainfo()
    paste0(info$summary['gene',1],'_',info$summary['hgvs_cdna',1],'_variantReport.pdf') },
    content = function(file) {
      #debug(loggerServer, 'VariantReport: getting cigmainfo')
      info<-cigmainfo()
      #debug(loggerServer, '               getting VariantReport')
      variantReport<-getVariantReport()
      #debug(loggerServer, '               getting VariantNotes')
      variantNotes<-getVariantNotes()
      #debug(loggerServer, '          done getting VariantNotes')
      savewd<-getwd()
      #debug(loggerServer, paste('current directory=',savewd))
      #debug(loggerServer, '          done getting VariantNotes')
      #debug(loggerServer, 'dumping variables')
      dump(c('info','variantReport','variantNotes'),file=paste(savewd,'reports','dump.txt',sep='/'))
      #debug(loggerServer, 'dumping done')
      #      setwd('reports') # knitr can only write in folder reports
      #      tryCatch({
      #debug(loggerServer, '               creating .tex from .Rnw')
      knit('reports/cavada_variant_report.Rnw')
      #debug(loggerServer, '               creating .pdf from .tex')
      system('/usr/bin/pdflatex cavada_variant_report.tex')

      #debug(loggerServer, paste('               renaming cavada_variant_report.pdf to',file))
      file.rename('cavada_variant_report.pdf', file) # move pdf to file for downloading
      #      },finally = {
      #        setwd(savewd)
      #debug(loggerServer, 'done creating .pdf from .Rnw')
      #      })
    },

    contentType = 'application/pdf'
  )

  output$downloadAnnotatedBatch <- downloadHandler(
    filename = function() { paste0(values$iname,'_',values$sessionid,'.curated.csv') },
    content = function(file) {write.csv(getAnnotatedBatch(), file, row.names=F, quote=T, eol="\r\n")}
  )

  getAnnotatedBatch <- reactive({
    if (values$inputFileType=='simple') {
      v<-data.frame(sampleID=values$exploreTable[,'sample_id'],
                    gene=sapply(values$exploreTable[,'variation'],function(x){return(unlist(str_split(x,':'))[1])}),
                    hgvs_cdna=sapply(values$exploreTable[,'variation'],function(x){return(unlist(str_split(unlist(str_split(x,':'))[2],' '))[1])}),
                    hgvs_protein=sapply(values$exploreTable[,'variation'],function(x){return(gsub('[()]','',unlist(str_split(unlist(str_split(x,':'))[2],' '))[2]))}),
                    varType=values$exploreTable[,'varType'],
                    varLocation=values$exploreTable[,'varLocation'],
                    codingEffect=values$exploreTable[,'codingEffect'],
                    automated_class=values$exploreTable[,'working_cigma_class'],
                    classification_justification=values$exploreTable[,'classification_justification'],
                    curated_class=values$exploreTable[,'clinical_cigma_class'])
      s<-lapply(1:dim(v)[1],
                function(y){
                  return(sapply(names(fan),
                                function(x,yy=y){
                                  table<-x
                                  if(table=='esp') table<-'esp_cappa'
                                  if(table=='hgmd') table<-'hgmd_cappa'
                                  prefix<-x
                                  if(table=='tgp_pop') {
                                    prefix<-'tgp'
                                    table<-'tgp_pop_main'
                                  }
                                  if(table=='tgp_phase3_pop') {
                                    prefix<-'tgp_phase3'
                                    table<-'tgp_phase3_pop_main'
                                  }
                                  fields <- paste(paste0(prefix,'_',fan[[x]]),collapse=',')
                                  if(table=='walker') fields <- paste0('GROUP_CONCAT(',fields,') ',fields)
                                  select_statement<-paste0("select ",fields," from ",table," where gene='",v[yy,2],"' and hgvs_cdna='",v[yy,3],"'")
                                  return(select_statement)
                                }
                  )
                  )
                }
      )
      mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
      extra<-lapply(1:length(s),
                    function(y) {
                      return(
                        unlist(sapply(1:length(s[[y]]),function(x){rs<-dbSendQuery(mydb,s[[y]][x]);dataset<-fetch(rs,-1);return(as.vector(dataset))}))
                      )
                    }
      )
      dbDisconnect(mydb)
      for(i in 1:length(extra)) {
        if(length(extra[[i]])>0) {
          for (j in 1:length(extra[[i]])) {
            v[i,names(extra[[i]])[j]]<-extra[[i]][j]
          }
        }
      }
      for (j in 1:dim(v)[2]) {
        if (length(which(is.na(v[,j])))>0) {
          v[which(is.na(v[,j])),j]<-''
        }
      }
      colnames(v)=c('Sample ID',
                    'Gene',
                    'HGVS cDNA',
                    'HGVS protein',
                    'Variant type',
                    'Variant location',
                    'Variant coding effect',
                    'Automated classification',
                    'Class justification',
                    'Proposed classification',
                    fa[colnames(v)[11:length(colnames(v))],3])
      ordered_fields<-fa[order(fa$analysis_output_order),'field_label']
      return(v[,c(colnames(v)[1:10],ordered_fields[which(ordered_fields %in% colnames(v))])])
    }
    if (values$inputFileType=='clinical') {
      v<-data.frame(Worksheet=values$exploreTable[,'worksheet'],
                    InvID=values$exploreTable[,'inv_id'],
                    DNA.No=values$exploreTable[,'sample_id'],
                    Gene=sapply(values$exploreTable[,'variation'],function(x){return(unlist(str_split(x,':'))[1])}),
                    Nomenclature=sapply(values$exploreTable[,'variation'],function(x){return(unlist(str_split(x,':'))[2])}),
                    Het.Hom=values$exploreTable[,'zygosity'],
                    Class=values$exploreTable[,'clinical_cigma_class'],
                    Exon=values$exploreTable[,'exon'])
      return(v)
    }
  })

  getVariantReport <- reactive({
    if (values$contentType=='exploreSingle' && !is.null(input$genetrans)){
      if (!nchar(input$genetrans)) {return()}
      gene<-unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))[1]
      transcript<-unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))[2]
      variation<-gsub(' [(].*$','',input$variationSel,perl=T)
    }else if (values$contentType=='browseTabs'){
      gene<-unlist(str_split(values$selectedVar,':'))[1]
      variation<-unlist(str_split(values$selectedVar,':'))[2]
    }else{
      return()
    }
    #debug(loggerServer,paste("in getVariantReport, gene=",gene," variation=",variation,sep=""))
    ordered_fields<-rownames(fv)[order(fv$variant_report_order)]
    info<-cigmainfo()
    #debug(loggerServer,'   got cigmainfo')
    info$summary['name',1] = gsub(info$summary['gene',1],'',info$summary['name',1])
    info$summary['name',1] = gsub('on','on ',info$summary['name',1])
    info$summary['name',1] = gsub('of',' of ',info$summary['name',1])
    variantReport<-data.frame(
      Description=c(
        'Gene',
        'Ensembl transcript',
        'RefSeq transcript',
        'Chromosome',
        'Genomic position',
        'Reference allele',
        'Variant allele',
        'HGVS cDNA',
        'HGVS protein',
        'HGVS protein (1 letter)',
        'Alternative notation',
        'rsID',
        'Variant location',
        'Exon/intron number',
        'Variant type',
        'Variant coding effect',
        'cDNA position',
        'Codon number',
        'Intronic position',
        'First/last 3 bases of exon'
      ),
      Value=c(
        gene,
        info$summary['transcript',1],
        info$summary['rtranscript',1],
        info$summary['hg19_chr',1],
        info$summary['hg19_pos',1],
        info$summary['ntwt',1],
        info$summary['ntmut',1],
        variation,
        info$summary['hgvs_prot',1],
        info$summary['hgvs_prot_code1',1],
        info$summary['altname',1],
        info$summary['rsID',1],
        info$summary['varLocation',1],
        info$summary['name',1],
        info$summary['varType',1],
        info$summary['codingEffect',1],
        info$summary['cdna_pos',1],
        info$summary['codon',1],
        info$summary['offset',1],
        info$summary['firstorlast3',1]
      ),
      row.names=c(
        'main_gene',
        'main_transcript',
        'main_rtranscript',
        'main_hg19_chr',
        'main_hg19_pos',
        'main_ntwt',
        'main_ntmut',
        'main_hgvs_cdna',
        'main_hgvs_prot',
        'main_hgvs_prot_code1',
        'main_altname',
        'main_rsID',
        'main_varLocation',
        'main_name',
        'main_varType',
        'main_codingEffect',
        'main_cdna_pos',
        'main_codon',
        'main_offset',
        'main_firstorlast3'
      )
    )
    # hack to insert the exon/intron number after varLocation
    ordered_fields<-c(ordered_fields[1:13],'main_name',ordered_fields[14:length(ordered_fields)])
    if (!is.null(info) && !is.null(info$swissprot)) {
      swissprots<-lapply(1:dim(info$swissprot)[1],function(x){
        return(c(paste0('Swissprot feature',x,' - ',info$swissprot[x,1]),info$swissprot[x,2]))
      })
      swissprot_fields<-NULL
      for(i in 1:length(swissprots)) {
        variantReport<-rbind(variantReport,data.frame(Description=swissprots[[i]][1],Value=swissprots[[i]][2],row.names=paste0('swissprot_feature_',i)))
        swissprot_fields<-c(swissprot_fields,paste0('swissprot_feature_',i))
      }
      # hack to insert the swissprot features after firstorlast3
      ordered_fields<-c(ordered_fields[1:17],swissprot_fields,ordered_fields[18:length(ordered_fields)])
    }

    t<-names(fvn)
    t<-t[-which(t=='main')]
    s<-lapply(t,
              function(x){
                table<-x
                if(table=='esp') table<-'esp_cappa'
                if(table=='hgmd') table<-'hgmd_cappa'
                prefix<-x
                if(table=='icr') {
                  table<-'icr_controls'
                  prefix<-'icr_controls'
                }
                if(table=='tgp_pop') {
                  prefix<-'tgp'
                  table<-'tgp_pop_main'
                }
                if(table=='tgp_phase3_pop') {
                  prefix<-'tgp_phase3'
                  table<-'tgp_phase3_pop_main'
                }
                fields <- paste(paste0(prefix,'_',fvn[[x]]),collapse=',')
                if(table=='walker') fields <- paste0('GROUP_CONCAT(',fields,') ',fields)
                select_statement<-paste0("select ",fields," from ",table," where gene='",gene,"' and hgvs_cdna='",variation,"'")
                return(select_statement)
              }
    )
    s<-setNames(s,names(t))
    mydb<-dbConnect(MySQL(max.con=200, fetch.default.rec=5000),user='batch',password='cigma',dbname='cigma2',host=host)
    extra<-lapply(1:length(s),function(x){
      #debug(loggerServer,paste('trying to run SQL statement:',s[[x]],sep='\n'))
      rs<-dbSendQuery(mydb,s[[x]]);
      dataset<-fetch(rs,-1);
      return(as.vector(dataset))
    }
    )
    extra<-setNames(extra,names(t))
    #   dbDisconnect(mydb)
    #debug(loggerServer,'adding extra fields')
    for(i in 1:length(extra)) {
      if(length(extra[[i]])>0) {
        for (j in 1:length(extra[[i]])) {
          field_label <- fv[names(extra[[i]])[j],'field_label']
          if (!is.na(field_label) && !is.na(extra[[i]][1,j])) {
            variantReport<-rbind(variantReport,data.frame(Description=field_label,Value=extra[[i]][1,j],row.names=names(extra[[i]])[j]))
          }
        }
      }
    }
    #debug(loggerServer,'done adding extra fields')
    if (!(gene %in% c('MSH2','MSH6','MLH1','PMS2')) && ('lovd_insight_class' %in% row.names(variantReport))){
      variantReport<-variantReport[-which(row.names(variantReport)=='lovd_insight_class'),]
    }
    variantReport<-variantReport[ordered_fields[which(ordered_fields %in% row.names(variantReport))],]
    variantReport<-variantReport[-which(is.na(variantReport$Value)),]
    variantReport[,1]<-paste0(variantReport[,1],':')
    variantReport1<-data.frame(Values=variantReport[,2],row.names=variantReport[,1])
    #debug(loggerServer,'done fixing the variantReport')

    ###>>> added to include frequenc, ..., etc. data to pdf report
    varReportExtra1 <- c()
    varReportExtra2 <- c()

    fields.phe2 <- list("BIRq3"="Birmingham", "BRIrvj"="Bristol", "CAMrgt"="Cambridge",
                        "GUYrj1"="Guy's", "LEEDrr8"="Leeds", "NEWrtd"="Newcastle", "MANrw3"="Manchester", "NOTrx1"="Nottingham",
                        "SALrnz"="Salisbury", "SHEFrcu"="Sheffield","scrCount"="Total proband count", "totCount"="Total probands tested")


    toselect <- paste0(names(fields.phe2), collapse = ", ")
    myquery <- paste(sep = "",collapse = "", "select ",toselect  ," from PHE where gene='",info$summary["gene",],"' and transcript='",info$summary["transcript",],"' and hgvs_cdna='", info$summary["hgvs_cdna",],"'")
    rs<-dbSendQuery(mydb,myquery)
    phest <- as.list(fetch(rs,-1))
    PHEprobands <- ifelse(info$summary["gene",]=="BRCA1" | info$summary["gene",]=="BRCA2",25773,0)
    ###>>>
    if(values$logged){

      if(all(sapply(phest,is.na))){
        varReportExtra1 <- c(varReportExtra1, "UK diagnostic labs (probands):")
        varReportExtra2 <- c(varReportExtra2, "NA")
      }else{
        varReportExtra1 <- c(varReportExtra1, "UK diagnostic labs (probands):")
        varReportExtra2 <- c(varReportExtra2, paste0(phest$scrCount," / ",PHEprobands))
      }
    }else{
      #varReportExtra1 <- c(varReportExtra1, "UK diagnostic labs frequency:")
      #varReportExtra2 <- c(varReportExtra2, "Please login to access this data")
      varReportExtra1 <- c(varReportExtra1, "UK diagnostic labs (probands):")
      varReportExtra2 <- c(varReportExtra2, paste0(phest$scrCount," / ",PHEprobands))
    }

    toselect <- paste0("NFE",", NFE_count",", NFE_no",", NFE_homo",", Total",", Total_count",", Total_no",", Total_homo")
    myquery <- paste(sep = "",collapse = "", "select ",toselect  ," from gnomadtest where gene='",info$summary["gene",],"' and hgvs_cdna='",info$summary["hgvs_cdna",],"'")
    rs<-dbSendQuery(mydb,myquery)
    dtst <- as.list(fetch(rs,-1))

    if(all(sapply(dtst,length)==0)){
      varReportExtra1 <- c(varReportExtra1, "gnomAD population frequencies:")
      varReportExtra2 <- c(varReportExtra2, "Variant not found in gNomad database")
    }else{
      varReportExtra1 <- c(varReportExtra1, "gnomAD population frequencies:")
      varReportExtra2 <- c(varReportExtra2, paste0(paste0(gsub("_"," ",names(dtst)),": ", signif(as.numeric(dtst),2)), collapse=", "))

      varReportExtra1 <- c(varReportExtra1, "")
      varReportExtra2 <- c(varReportExtra2, paste0("NFE=Non-Finnish European, ", "The GNOMAD dataset contains 125,748 exome", "and 15,708 whole-genome sequences", "from unrelated individuals", "(141,456 individuals, 282,912 alleles)"))
    }

    variantReport1 <- data.frame(Values=c(variantReport1[,1],varReportExtra2),row.names=c(rownames(variantReport1),varReportExtra1))
    ###
    dbDisconnect(mydb)
    return(variantReport1)
  })

  sanitize <- function(str) {
    Sys.setlocale('LC_ALL','C')
    result <- str
    result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
    result <- gsub("$", "\\$", result, fixed = TRUE)
    result <- gsub(">", "$>$", result, fixed = TRUE)
    result <- gsub("<", "$<$", result, fixed = TRUE)
    result <- gsub("|", "$|$", result, fixed = TRUE)
    result <- gsub("{", "\\{", result, fixed = TRUE)
    result <- gsub("}", "\\}", result, fixed = TRUE)
    result <- gsub("%", "\\%", result, fixed = TRUE)
    result <- gsub("_", "\\_", result, fixed = TRUE)
    result <- gsub("#", "\\#", result, fixed = TRUE)
    result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
    result <- gsub("~", "\\~{}", result, fixed = TRUE)
    result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$", result, fixed = TRUE)
    return(result)
  }

  getVariantNotes<-reactive({
    info<-cigmainfo()
    info(loggerServer,'in getVariantNotes')
    variantNotes<-data.frame()
    if (!is.null(info) && !is.null(info$notes) && !is.na(info$notes)) {
      mynotes<-info$notes[,c('date','notes')]
      #debug(loggerServer,'sanitizing mynotes')
      mynotes$notes<-sanitize(mynotes$notes)
      #debug(loggerServer,'adding PMID URLs')
      mynotes$notes<-gsub('PMID[:]?(\\s+)?(\\d+)','\\\\href{http://www.ncbi.nlm.nih.gov/pubmed?term=\\2}{PMID:\\2}',gsub('<br>','\\baselineskip',mynotes$notes))
      colnames(mynotes)<-c('Date','Notes')
      variantNotes<-rbind(variantNotes,mynotes)
      if(!values$logged){
        variantNotes[,2]  <-  "The comments are only visible for the registered users. Please login !!!"
      }
    }
    #debug(loggerServer,'done with notes')
    return(variantNotes)
  })
  #################
  # Other outputs #
  #################

  output$exploreTableVariation <- renderText({
    input$exploreTableVariation
  })

  # Output variation name (for the browseTab screen)
  output$selectedVar <- renderText({
    values$selectedVar
  })

  output$selectedVarBaseline <- renderText({
    ifelse(values$selectedVarBaseline!='',values$selectedVarBaseline,' ')
  })

  output$selectedSet <- renderText({
    input$batchid
  })

})
