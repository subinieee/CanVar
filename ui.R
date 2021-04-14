library(shiny)
#library(shinyIncubator)
library(RMySQL)
options(stringsAsFactors=F)
library(markdown)
shinyUI(
 
  fluidPage(
#theme="css/united.bootstrap.css",
            title="CanVar-UK",
            #progressInit(),
            tagList(
              singleton(
                tags$head(
                 includeHTML(("www/google-analytics.html")),
                  tags$script(src = "js/jsCodeMessage.js"), 
                  #tags$script(src = "js/disable_back.js"), 
			tags$style("@import url(https://use.fontawesome.com/releases/v5.7.2/css/all.css);"),                
                  tags$style(type="text/css", '#version {color:rgb(2, 182, 230); font-size: 12px; margin-left:10px;}'),
                  tags$style(type="text/css", '.form-group, .progress {margin-bottom:0px;}'),
                  tags$style(type="text/css", '.btn {margin-top:0px;}'),
                  tags$style(type="text/css", '#visibleComments {margin-top: 20px;
                            margin-bottom: 30px;
                            padding-top: 0px;
                            padding-bottom: 0px;"}'),
                  tags$style(type="text/css", 'table.dataTable tr.odd{ background-color: white !important}'),
                  tags$style(type="text/css", 'table.dataTable tr.even{ background-color: white !important}'),
                  
                  tags$style(type="text/css", '.row {margin-left:auto; margin-right:auto; !important;}'),
                  
                  tags$style(type="text/css", '.help-block {color: #A0A0A0; margin-bottom:5px;}'),
                  tags$style(type="text/css", 'select,input[type="text"],input[type="password"] {height:33px;margin-bottom:5px;}'),
                  tags$style(type="text/css", 'label[for=offsetLimit] {margin-left:20px; display:inline-block;}'),
                  tags$style(type="text/css", '#offsetLimit {width:50px; margin-bottom:0px; margin-left:5px; margin-right:5px;}'),
                  tags$style(type="text/css", '.thin {min-height:20px !important;}'),
		  # for orange color set below line rgba(232, 101, 55, 0.6)      
                  tags$style(type="text/css",".highlight {background-color: rgb(0, 191, 255) !important}"),
                  tags$style(type="text/css",'.row-fluid [class*="span"] {min-height: 20px !important}'),
		            tags$style(type="text/css", "#visibleComments{color:red; font-size:12px; font-weight:bold; font-style:italic; 
	border-top:transparent; border-bottom:transparent; border-left:transparent; border-right:transparent;
	overflow-y:#scroll; max-height: 20px; background: #f5f5f5; margin-left:20px;}"),
	
                  tags$style(type="text/css",'.details {background-color: white !important;}'),
		              HTML('<style type="text/css"> td > strong {margin-left:5px;}</style>')
                )
              )
            ),

		  uiOutput('uiControl'),

		  uiOutput('uiContents')
  )
)
