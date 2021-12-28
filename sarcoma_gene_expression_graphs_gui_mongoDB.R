## load libraries
library(dplyr)
library(tidyverse)
library(tibble)
library(readxl)
library(survival)
library(ggplot2)
library(survminer)
library(ggpubr)
library(shiny)
library(survMisc)
library(janitor)
library(data.table)
library(mongolite)



connection_string <- "mongodb://sarcoma-manager:MongoDb-1dfae4e9912c99eea55ac@132.66.207.18:80/sarcomadb?authSource=admin"
db_clinicalData <-  mongo(collection="ClinicalData", db="sarcomadb", url=connection_string)
db_Gene_Expression <-  mongo(collection="RNASeq", db="sarcomadb", url=connection_string)
#filtering clinical data 
OS_Months=paste(shQuote("OS_MONTHS", type="cmd"), collapse=", ")
OS_Status=paste(shQuote("OS_STATUS", type="cmd"), collapse=", ")
OS_Months_Query =paste0('{"name":{"$in":[',OS_Months,']}}')
OS_Status_Query =paste0('{"name":{"$in":[',OS_Status,']}}')
clinicalData_OS_Months <- db_clinicalData$find(fields = '{
}',query =OS_Months_Query)
clinicalData_OS_Status <- db_clinicalData$find(fields = '{
}',query =OS_Status_Query)
clinical<- merge(clinicalData_OS_Months,clinicalData_OS_Status[,c("patient","value")],by.x = "patient", by.y = "patient")
colnames(clinical)[colnames(clinical) == 'value.x'] <- 'OverallSurvivalMonths'
colnames(clinical)[colnames(clinical) == 'value.y'] <- 'OverallSurvivalStatus'



# functions 
Gene_Query<-function(Gene){
  db_Gene_Expression <-  mongo(collection="RNASeq", db="sarcomadb", url=connection_string)
  Gene=paste(shQuote(Gene, type="cmd"), collapse=", ")
  OS_Gene =paste0('{"name":{"$in":[',Gene,']}}')
  Gene_expression <- db_Gene_Expression$find(fields = '{
}',query =OS_Gene)
}

data_Expression_level_cut_off<-function(Expresion_table1,Expresion_table2,clinical_table,input_main_expressoin,input_secondery_Expression,Expression_level){
   if (input_main_expressoin == "HIGH") {
     Expresion_table1$Main_Expression_Level <- ifelse(Expresion_table1[,"value"]>quantile(Expresion_table1[,"value"],probs = Expression_level),"HIGH","LOW")
    }
    if (input_main_expressoin == "LOW") {
      Expresion_table1$Main_Expression_Level <- ifelse(Expresion_table1[",value"]>quantile(Expresion_table1[,"value"],probs = 1 - Expression_level),"HIGH","LOW")
    }
    if (input_secondery_Expression == "HIGH") {
      Expresion_table2$Secondary_Expression_Level <- ifelse(Expresion_table2[,"value"]>quantile(Expresion_table2[,"value"],probs =  Expression_level),"HIGH","LOW")
    }
    if (input_secondery_Expression == "LOW") {
      Expresion_table2$Secondary_Expression_Level <- ifelse(Expresion_table2[,"value"]>quantile(Expresion_table2[,"value"],probs = 1 - Expression_level),"HIGH","LOW")
    }
  data<-merge(Expresion_table1,Expresion_table2[,c("patient","Secondary_Expression_Level")],by.x = "patient", by.y ="patient")
  data$CEL <- ifelse(data$Main_Expression_Level == "HIGH" & data$Secondary_Expression_Level == "HIGH","HIGH_HIGH",
                                          ifelse(data$Main_Expression_Level == "HIGH" & data$Secondary_Expression_Level == "LOW","HIGH_LOW",
                                                 ifelse(data$Main_Expression_Level == "LOW" & data$Secondary_Expression_Level == "HIGH","LOW_HIGH",
                                                        ifelse(data$Main_Expression_Level == "LOW" & data$Secondary_Expression_Level == "LOW","LOW_LOW"," "))))

  data<-merge(data,clinical_table,by.x = "patient", by.y ="patient")
  names(data)<-gsub(" ","",names(data))
  names(data)<-gsub("(","",names(data),fixed = TRUE)
  names(data)<-gsub(")","",names(data))
  names(data)<-gsub("'","",names(data))
  data<-subset(data,OverallSurvivalStatus == "DECEASED" | OverallSurvivalStatus == "LIVING")
  data$OverallSurvivalMonths<-as.numeric(data$OverallSurvivalMonths,na.rm=T)
  data$Status[data$OverallSurvivalStatus=="DECEASED"]<-1
  data$Status[data$OverallSurvivalStatus=="LIVING"]<- 0
  data$Status<-as.numeric(data$Status)
  return(data)
}

method_km<- c("survdiff","1","n","sqrtN","S1","S2","FH_p=1_q=1")
high_Expression_80 <-0.8
low_Expression_20<-0.2
high_Expression_70 <-0.7
low_Expression_30<-0.3

ui <- fluidPage(
  
  # Application title
  titlePanel("KM by Expression"),
  
  
  sidebarLayout(
    sidebarPanel(
      selectInput("Main_gene","choose driver gene",choices = db_Gene_Expression$distinct(key = "name")),
      selectInput("Main_Expression","Choose Driver Gene Expression Level (20%)", choices = c("HIGH","LOW")),
      selectInput("Secondery_gene","choose secondery gene",choices = db_Gene_Expression$distinct(key = "name")),
      selectInput("Secondery_Expression","Choose secondery Gene Expression Level (20%)", choices = c("HIGH","LOW")),
      selectInput("method_km","Choose p-value calcualation method for kaplan meiers", choices = method_km),
      sliderInput(inputId = "x_axis",label = "months axis KM",min = 5,max = 400,value = 20)
    ),
    
    
    mainPanel(
      tabsetPanel(type = "tab",
                  tabPanel("80%-20%",
                           textOutput("text_main"),
                           textOutput("text_second"),
                           tableOutput("Driver_gene"),
                           tableOutput("secondary_gene"),
                           tableOutput("Combined_genes"),
                           plotOutput("CombinedKM")),
                  tabPanel("70%-30%",
                           textOutput("text_main_30_70"),
                           textOutput("text_second_30_70"),
                           tableOutput("Driver_gene_30_70"),
                           tableOutput("secondary_gene_30_70"),
                           tableOutput("Combined_genes_30_70"),
                           plotOutput("combinedKm30_70")))  
    )
  )
)

server <- function(input, output,session) {
  

  # prepare data for 80-20 cross between genes
  choosen_gene<-reactive({
    data<-data_Expression_level_cut_off(Expresion_table1 = Gene_Query(input$Main_gene),
                                        Expresion_table2 = Gene_Query(input$Secondery_gene),
                                        clinical_table = clinical,
                                        Expression_level = high_Expression_80,
                                        input_main_expressoin = input$Main_Expression,
                                        input_secondery_Expression = input$Secondery_Expression)
    return(data)
  })
  
  # prepare data for 70-30 cross between genes
  choosen_gene_30_70<-reactive({
    data<-data_Expression_level_cut_off(Expresion_table1 = Gene_Query(input$Main_gene),
                                        Expresion_table2 = Gene_Query(input$Secondery_gene),
                                        clinical_table = clinical,
                                        Expression_level = high_Expression_70,
                                        input_main_expressoin = input$Main_Expression,
                                        input_secondery_Expression = input$Secondery_Expression)
    return(data)
  })
  obsB <- observe({
    (choosen_gene())
  })
  
  km_fit1 <- reactive({
    #print driver gene p value for kaplan meier analysis betwwen himself 80-20 retio 
    km <- survfit(Surv(OverallSurvivalMonths, Status) ~ Main_Expression_Level , data=choosen_gene())
    text<-surv_pvalue(km,method = input$method_km,data = choosen_gene())$pval.txt
    text<- paste("Driver Gene cross Kaplan meier P_value:", text)
    return(text)
    
    
  })
  #print driver gene p value for kaplan meier analysis betwwen himself 70-30 retio 
  km_fit1_30_70 <- reactive({
    
    km <- survfit(Surv(OverallSurvivalMonths, Status) ~ Main_Expression_Level , data=choosen_gene_30_70())
    text<-surv_pvalue(km,method = input$method_km,data = choosen_gene_30_70())$pval.txt
    text<- paste("Driver Gene cross Kaplan meier P_value:", text)
    return(text)
    
  })
  #print secondary gene p value for kaplan meier analysis betwwen himself 80-20 retio 
  second_gene <- reactive({
    km <- survfit(Surv(OverallSurvivalMonths, Status) ~ Secondary_Expression_Level , data=choosen_gene())
    text<-surv_pvalue(km,  method = input$method_km ,data = choosen_gene())$pval.txt
    text<- paste("Second Gene cross Kaplan meier P_value:", text)
    return(text)
  })
  
  #print secondary gene p value for kaplan meier analysis betwwen himself 70-30 retio 
  second_gene_30_70 <- reactive({
    km <- survfit(Surv(OverallSurvivalMonths, Status) ~ Secondary_Expression_Level , data=choosen_gene_30_70())
    text<-surv_pvalue(km,method = input$method_km, data = choosen_gene_30_70())$pval.txt
    text<- paste("Second Gene cross Kaplan meier P_value:", text)
    return(text)
  })
  
  first_Expression_table <- reactive({
    data <- data.frame(
      High_Driver_Gene = length(choosen_gene()$Main_Expression_Level[choosen_gene()$Main_Expression_Level=="HIGH"]),
      Low_Driver_Gene = length(choosen_gene()$Main_Expression_Level[choosen_gene()$Main_Expression_Level=="LOW"])
    )
    return(as.data.frame(data))
  })
  # count how many patients in Driver Gene cross
  first_Expression_table30_70 <- reactive({
    data <- data.frame(
      High_Driver_Gene = length(choosen_gene_30_70()$Main_Expression_Level[choosen_gene_30_70()$Main_Expression_Level=="HIGH"]),
      Low_Driver_Gene = length(choosen_gene_30_70()$Main_Expression_Level[choosen_gene_30_70()$Main_Expression_Level=="LOW"])
    )
    return(as.data.frame(data))
  })
  second_Expression_table <- reactive({
    data<-data.frame( 
      High_Secondary_Gene = length(choosen_gene()$Secondary_Expression_Level[choosen_gene()$Secondary_Expression_Level=="HIGH"]),
      Low_Secondary_Gene = length(choosen_gene()$Secondary_Expression_Level[choosen_gene()$Secondary_Expression_Level=="LOW"])
    )
    return(as.data.frame(data))
  })
  # count how many patients in Secondery Gene cross
  second_Expression_table30_70 <- reactive({
    data<-data.frame( 
      High_Secondary_Gene = length(choosen_gene_30_70()$Secondary_Expression_Level[choosen_gene_30_70()$Secondary_Expression_Level=="HIGH"]),
      Low_Secondary_Gene = length(choosen_gene_30_70()$Secondary_Expression_Level[choosen_gene_30_70()$Secondary_Expression_Level=="LOW"])
    )
    return(as.data.frame(data))
  })
  
  CombinedGenes <- reactive({
    df<-data.frame( HIGH_HIGH = length(choosen_gene()$CEL[choosen_gene()$CEL=="HIGH_HIGH"]),
                    HIGH_LOW = length(choosen_gene()$CEL[choosen_gene()$CEL=="HIGH_LOW"]),
                    LOW_HIGH = length(choosen_gene()$CEL[choosen_gene()$CEL=="LOW_HIGH"]),
                    LOW_LOW = length(choosen_gene()$CEL[choosen_gene()$CEL=="LOW_LOW"]))
    return(df)
  })
  
  CombinedGenes30_70 <- reactive({
    df<-data.frame( HIGH_HIGH = length(choosen_gene_30_70()$CEL[choosen_gene_30_70()$CEL=="HIGH_HIGH"]),
                    HIGH_LOW = length(choosen_gene_30_70()$CEL[choosen_gene_30_70()$CEL=="HIGH_LOW"]),
                    LOW_HIGH = length(choosen_gene_30_70()$CEL[choosen_gene_30_70()$CEL=="LOW_HIGH"]),
                    LOW_LOW = length(choosen_gene_30_70()$CEL[choosen_gene_30_70()$CEL=="LOW_LOW"]))
    return(df)
  })
  
  Combined_km <- reactive({
    km <- survfit(Surv(OverallSurvivalMonths, Status) ~ CEL , data=choosen_gene())
    GGkm<-ggsurvplot(km,data = choosen_gene(),pval = TRUE, pval.method = TRUE,
                     log.rank.weights = input$method_km, pval.method.coord = c(5, 0.1),
                     pval.method.size = 3,
                     break.x.by=20,
                     xlim = c(0,input$x_axis),
                     axes.offset = FALSE
                     
    )
    return(GGkm)
  })
  
  Combined_km30_70 <- reactive({
    km <- survfit(Surv(OverallSurvivalMonths, Status) ~ CEL , data=choosen_gene_30_70())
    GGkm<-ggsurvplot(km,data = choosen_gene_30_70(),pval = TRUE, pval.method = TRUE,
                     log.rank.weights = input$method_km, pval.method.coord = c(5, 0.1),
                     pval.method.size = 3,
                     break.x.by=20,
                     xlim = c(0,input$x_axis),
                     axes.offset = FALSE
                     
                     
    )
    return(GGkm)
  })
  
  ## output funcitons 
  output$text_main <- renderText({
    km_fit1()
  })
  output$text_second <- renderText({
    second_gene()
  })
  output$Driver_gene <- renderTable({
    table<-first_Expression_table()
  })
  output$Combined_genes <- renderTable({
    table <- CombinedGenes()
  })
  output$secondary_gene <- renderTable({
    second_Expression_table()
  })
  output$CombinedKM <- renderPlot({
    Combined_km()
  })
  output$text_main_30_70 <- renderText(
    km_fit1_30_70()
  )
  output$text_second_30_70 <- renderText(
    second_gene_30_70()
  )
  output$Driver_gene_30_70<- renderTable({
    first_Expression_table30_70()
  })
  output$secondary_gene_30_70 <- renderTable({
    second_Expression_table30_70()
  })
  output$Combined_genes_30_70 <- renderTable({
    CombinedGenes30_70()
  })
  output$combinedKm30_70 <- renderPlot({
    Combined_km30_70()
  })
}
 
# Run the application 
shinyApp(ui = ui, server = server)
