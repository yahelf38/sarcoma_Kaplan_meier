
library(dplyr)
library(tidyverse)
library(tibble)
library(readxl)
library(survival)
library(ggplot2)
library(survminer)
library(ggpubr)
library(shiny)
library(data.table)


#load data
clinicalData <- read_excel("data_clinical_patient.xlsx")
data_mRNA_median_all_sample_Zscores<-read_excel("data_RNA_Seq_v2_mRNA_median_Zscores.xlsx") 

Expression<-data_mRNA_median_all_sample_Zscores[,-2]
Expression<-as.data.table(Expression)
Expression<-Expression[ ,lapply(.SD, mean), by = Hugo_Symbol] #remove duplicate genes and replace Expression values as means 
v<-as.vector(Expression$Hugo_Symbol)
Expression<-as.data.frame(Expression)
rownames(Expression)<-v
Expression<-Expression[,-1]
Expression<-t(Expression)
Expression<-as.data.frame(Expression)
patient_names<-rownames(Expression)
patient_names<- substr(patient_names,1,12)
rownames(Expression)<-patient_names
Expression<-as.data.frame(Expression)
##### choose method for KM P-value calcuatlion 
# "survdiff", log-rank;
# "1": log-rank, LR; –> Regular log-rank test, sensitive to detect late differences.
# "n": Gehan-Breslow (generalized Wilcoxon), GB; –> detect early differences.
# "sqrtN": Tarone-Ware, TW; –> detect early differences.
# "S1": Peto-Peto's modified survival estimate, PP; –> more robust than Tharone-Whare or Gehan-Breslow, detect early differences
# "S2": modified Peto-Peto (by Andersen), mPP
# "FH_p=1_q=1": Fleming-Harrington(p=1, q=1), FH
#("log-rank", "log-rank-LR","Gehan-Breslow (generalized Wilcoxon)","Tarone-Ware",
#"Peto-Peto's"," modified Peto-Peto (by Andersen)","Fleming-Harrington")

method_km<- c("survdiff","1","n","sqrtN","S1","S2","FH_p=1_q=1")
downloadObjUI <- function(id) {
  ns <- NS(id)
  
  downloadButton(ns("data_download"), label = "Download Data", class = "btn-primary")
}
downloadObj <- function(input, output, session, data) {
  
  output$data_download <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(data(), file) # add parentheses to data arg if reactive
    }
  )
}
ui <- fluidPage(
  
  
  titlePanel("Making gene table"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput("Main_Expression","Choose Driver Gene Expression Level (20%)", choices = c("HIGH","LOW")),
      selectInput("Secondery_Expression","Choose secondery Gene Expression Level (20%)", choices = c("HIGH","LOW")),
      selectInput("driver_gene","choose driver gene",choices = names(Expression),selected = "MET"),
      textInput("list_of_genes","choose secondery genes" , value = "MET"),
      selectInput("method_km","Choose p-value calcualation method for kaplan meiers", choices = method_km),
      p("choose method for KM P-value calcuatlion",br(),
        "survdiff, log-rank;  ",br(),
        "1: log-rank, LR; –> Regular log-rank test",br(),
        "n: .Gehan-Breslow (generalized Wilcoxon), GB; –> detect early differences",br(),
        "sqrtN: Tarone-Ware, TW; –> detect early differences.",br(),
        "S1: Peto-Peto's modified survival estimat,",br(),
        "PP; –> more robust than Tharone-Whare or Gehan-Breslow, detect early differences",br(),
        "S2: modified Peto-Peto (by Andersen), mPP",br(),
        "FH_p=1_q=1: Fleming-Harrington(p=1, q=1), FH"),br(),
      
    ),
    mainPanel(
      tableOutput("table"),
      textOutput("driver"),
      downloadObjUI(id = "download1")
    ))
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  Making_the_table<-reactive({
    b<-input$list_of_genes
    b<-strsplit(b, " ")[[1]]
    b<-as.character(b)
    b <- b[b %in% colnames(Expression)]
    b<-as.character(b)
    
  
    
    driver<-input$driver_gene
    preData<-dplyr::select(Expression,c(driver,b))
    
    ##merge data with survival analysis
    data<-merge(preData,clinicalData,by.x = 0, by.y ='#Patient Identifier')
    names(data)<-gsub(" ","",names(data))
    names(data)<-gsub("(","",names(data),fixed = TRUE)
    names(data)<-gsub(")","",names(data))
    names(data)<-gsub("'","",names(data))
    data<-subset(data,OverallSurvivalStatus == "DECEASED" | OverallSurvivalStatus == "LIVING")
    data$OverallSurvivalMonths<-as.numeric(data$OverallSurvivalMonths,na.rm=T)
    data$Status[data$OverallSurvivalStatus=="DECEASED"]<-1
    data$Status[data$OverallSurvivalStatus=="LIVING"]<- 0
    data$Status<-as.numeric(data$Status)
    
    
    ## change level of Expression of the driver gene based on user input also based on 20/80 ratio
    if (input$Main_Expression == "HIGH"){
      ReplacebaleDB <- data.frame(row.names = data$Row.names , DriverGene = data[,driver],
                                  High_Low = ifelse(data[,driver]>quantile(data[,driver],probs = 0.8),"HIGH","LOW"),
                                  OverallSurvivalMonths=data$OverallSurvivalMonths,
                                  OverallSurvivalStatus=data$Status)
    }
    
    if (input$Main_Expression == "LOW"){
      ReplacebaleDB <- data.frame(row.names = data$Row.names , DriverGene = data[,driver],
                                  High_Low = ifelse(Expression[,driver]>quantile(Expression[,driver],probs = 0.2),"HIGH","LOW"),
                                  OverallSurvivalMonths=data$OverallSurvivalMonths,
                                  OverallSurvivalStatus=data$Status)
    }
    
    Vector_for_High_High_Expression<-vector()
    Vector_for_High_Low_Expression<-vector()
    Vector_for_Low_High_Expression<-vector()
    Vectorfor_Low_Low_Expression<-vector()
    Vector_for_KM<-vector()
    
    if(input$Secondery_Expression=="HIGH"){
      for (i in b) {
        ReplacebaleDB$variablegene<-data[,i]
        ReplacebaleDB$variablegene_Expression_Level <- ifelse(ReplacebaleDB$variablegene>quantile(ReplacebaleDB$variablegene,probs = 0.8),"HIGH","LOW")
        ReplacebaleDB$CEL <- ifelse(ReplacebaleDB$High_Low == "HIGH" & ReplacebaleDB$variablegene_Expression_Level == "HIGH","HIGH_HIGH",
                                    ifelse(ReplacebaleDB$High_Low == "HIGH" & ReplacebaleDB$variablegene_Expression_Level == "LOW","HIGH_LOW",
                                           ifelse(ReplacebaleDB$High_Low == "LOW" & ReplacebaleDB$variablegene_Expression_Level == "HIGH","LOW_HIGH",
                                                  ifelse(ReplacebaleDB$High_Low == "LOW" & ReplacebaleDB$variablegene_Expression_Level == "LOW","LOW_LOW"," "))))
        Var_High_High<-length(ReplacebaleDB$CEL[ReplacebaleDB$CEL=="HIGH_HIGH"])
        Var_High_Low<-length(ReplacebaleDB$CEL[ReplacebaleDB$CEL=="HIGH_LOW"])
        Var_Low_High<-length(ReplacebaleDB$CEL[ReplacebaleDB$CEL=="LOW_HIGH"])
        Var_Low_Low<-length(ReplacebaleDB$CEL[ReplacebaleDB$CEL=="LOW_LOW"])
        ## run KM models 
        data_For_Km<-data.frame(OverallSurvivalMonths=ReplacebaleDB$OverallSurvivalMonths,OverallSurvivalStatus = ReplacebaleDB$OverallSurvivalStatus, CEL= ReplacebaleDB$CEL)
        km <- survfit(Surv(time = OverallSurvivalMonths,event = OverallSurvivalStatus) ~ CEL,data=data_For_Km)
        text<-as.character(surv_pvalue(km,method = input$method_km,data = data_For_Km)$pval.txt)
        
        
        Vector_for_High_High_Expression <- append(Vector_for_High_High_Expression,Var_High_High)
        Vector_for_High_Low_Expression <- append(Vector_for_High_Low_Expression,Var_High_Low)
        Vector_for_Low_High_Expression <- append(Vector_for_Low_High_Expression,Var_Low_High)
        Vectorfor_Low_Low_Expression <- append(Vectorfor_Low_Low_Expression,Var_Low_Low)
        Vector_for_KM <- append(Vector_for_KM,text)
      }
    }
    
    if(input$Secondery_Expression=="LOW"){
      for (i in b) {
        ReplacebaleDB$variablegene<-data[,i]
        ReplacebaleDB$variablegene_Expression_Level <- ifelse(ReplacebaleDB$variablegene>quantile(ReplacebaleDB$variablegene,probs = 0.2),"HIGH","LOW")
        ReplacebaleDB$CEL <- ifelse(ReplacebaleDB$High_Low == "HIGH" & ReplacebaleDB$variablegene_Expression_Level == "HIGH","HIGH_HIGH",
                                    ifelse(ReplacebaleDB$High_Low == "HIGH" & ReplacebaleDB$variablegene_Expression_Level == "LOW","HIGH_LOW",
                                           ifelse(ReplacebaleDB$High_Low == "LOW" & ReplacebaleDB$variablegene_Expression_Level == "HIGH","LOW_HIGH",
                                                  ifelse(ReplacebaleDB$High_Low == "LOW" & ReplacebaleDB$variablegene_Expression_Level == "LOW","LOW_LOW"," "))))
        Var_High_High<-length(ReplacebaleDB$CEL[ReplacebaleDB$CEL=="HIGH_HIGH"])
        Var_High_Low<-length(ReplacebaleDB$CEL[ReplacebaleDB$CEL=="HIGH_LOW"])
        Var_Low_High<-length(ReplacebaleDB$CEL[ReplacebaleDB$CEL=="LOW_HIGH"])
        Var_Low_Low<-length(ReplacebaleDB$CEL[ReplacebaleDB$CEL=="LOW_LOW"])
        ## run KM models 
        
        data_For_Km<-data.frame(OverallSurvivalMonths=ReplacebaleDB$OverallSurvivalMonths,OverallSurvivalStatus = ReplacebaleDB$OverallSurvivalStatus, CEL= ReplacebaleDB$CEL)
        km <- survfit(Surv(time = OverallSurvivalMonths,event = OverallSurvivalStatus) ~ CEL,data=data_For_Km)
        text<-as.character(surv_pvalue(km,method = input$method_km,data = data_For_Km)$pval.txt)
        km <- survfit(Surv(time = OverallSurvivalMonths,event = OverallSurvivalStatus) ~ CEL,data=ReplacebaleDB)
        text<-as.character(surv_pvalue(km,data = ReplacebaleDB)$pval.txt)
        
        Vector_for_High_High_Expression <- append(Vector_for_High_High_Expression,Var_High_High)
        Vector_for_High_Low_Expression <- append(Vector_for_High_Low_Expression,Var_High_Low)
        Vector_for_Low_High_Expression <- append(Vector_for_Low_High_Expression,Var_Low_High)
        Vectorfor_Low_Low_Expression <- append(Vectorfor_Low_Low_Expression,Var_Low_Low)
        Vector_for_KM <- append(Vector_for_KM,text)
      }
    }
    
    #make Expression table 
    finalDataBase <- data.frame("Genes"= b,
                                "High_High"=Vector_for_High_High_Expression,
                                "High_Low"=Vector_for_High_Low_Expression,
                                "Low_High"=Vector_for_Low_High_Expression,
                                "Low_Low"=Vectorfor_Low_Low_Expression,
                                "Co_ExpressionP_value"=Vector_for_KM)
    return(finalDataBase)
  })
  
  
  
  whats_driver<-reactive({
    x<-input$driver_gene
    return(x)
  })
  
  ##plot functions 
  output$table <- renderTable({
    Making_the_table()
  })
  ## 
  
  callModule(downloadObj, id = "download1", data = Making_the_table)
  
}

# Run the application 
shinyApp(ui = ui, server = server)

