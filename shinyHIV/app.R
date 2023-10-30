library(shiny)

source("HIV.R")
source("RTI.R")

###Add dashed red line on CD4+ to mark AIDS threshold
###summary includes selected options, if 1 year selected "demonstrates withing host dynamics 365 post infection"
###make option to view cell counts side by side-- only if scale allows i.e. cd4+ and infected cd4+ on same plot


ui <- fluidPage(

  
    navbarPage("Modeling the Within-Host Dynamics of HIV",

               tabPanel("How HIV Works",
                        titlePanel("How HIV Works: Between & Within Host(s) Dynamics"),
                        
                        mainPanel(
                            "Acquired Immunodeficiency Syndrome, also known as AIDS, was recognized
                             as a new disease in 1981 after a series of young, previously healthy
                             men died from unusual, and in many cases, rare, infections
                             and malignancies [20], [7], [4].
                             On April 23, 1984 
                             Human Immunodeficiency Virus Type 1 was identified as the virus that causes AIDS [20] https://www.hiv.gov/hiv-basics/overview/history/hiv-and-aids-timeline/#year-1984.
                             At least 85.6 million people have been infected by HIV over the last 4
                             decades and 40.4 million have died [23].",
                            br(),
                            br(),
                            #"As a sexually transmitted disease",
                            "HIV is spread, primarily, through sex. Additionally, HIV can
                             be contracted through exposure to an infected persons bodily fluids. 
                             Transmission most commonly occurs through contact with an open wound,
                             perinatally (during pregnancy), from the
                             use, or accidental use of, an infected persons needle(s) [20].",
                            br(),
                            br(),
                            "HIV results in a chronic, progressive disease that has
                             no cure. HIV is hallmarked by its high viral reproduction rate and
                             destruction of the immune system. The virus targets CD4+ T-cells, a type
                             of lymphocyte (lymphocytes are white blood cells responsible for riding
                             the body of invading antigens and infected cells). CD4+ T-cells are
                             responsible for detecting attacks on the bodies own cells and sending
                             signals to other 'killer' lymphocytes triggering an immune response.",
                            br(),
                            br(),
                            "Overtime HIV reduces the CD4+ T-cell population to roughly 20% of normal
                             levels [19]. Such a low T-cell count leaves the body susceptible to
                             opportunistic infections and cancers. The near total destruction of the
                             immune system, perpetrated by HIV, makes treating such infections exceptionally
                             difficult and survival rates plumit.",
                            br(),
                            br(),
                            h3("Stages of HIV"),
                            #br(),
                            "HIV has 3 distinct stages of infection: Acute HIV
                             Infection, Chronic HIV Infection, and Aquired Immune Deficiency Syndrome. The
                             state of infection is determined based upon a patient’s CD4+ T-cell count,
                             viral load, and presenting symptoms.",
                            br(),
                            br(),
                            "Acute HIV Infection: begins 2-4 weeks after the initial
                             infection, lasting for ~8-10 weeks [12], [8]. Acute infection is hallmarked by
                             flu-like symptoms that are often accompanied by rash and swollen lymph nodes.
                             In this time the viral load is very high, as is the likelihood of transmission
                             [12]. During acute infection the patients CD4+ T-cell count will drop
                             dramatically[8]. It is at this time that the patient will benefit the most
                             from beginning antiretroviral treatment.",
                            br(),
                            br(),
                            "Chronic HIV Infection (also called Asymptomatic HIV or Clinical Latency):
                             lasts for approxiamtely 10 years post
                             Acute Infection, although, for some, Clinical Latency can be as short as
                             5 years [12]. Patients taking effective treatment regimens can prolong this
                             stage for up to 30 years. Clinical Latency patients will have very minimal
                             symptoms. The viral load reaches 'steady state', meaning the number of viruses
                             in the body stays relatively stable-- typically at lower numbers--
                             CD4+ T-cell counts 
                             near normal levels [12], [8]. Towards the end of Chronic Infection, CD4+ T-cell counts will
                             decline rapidly as the viral load begins to rise. People taking effective CARV’s
                             (combination antiretrovirals) can arrive at and maintain an undetectable viral
                             load during this stage, making the risk of transmitting HIV through sex
                             nearly 0 [12].",
                            br(),
                            br(),
                            "Acquired Immunodeficiency Syndrome (AIDS): Eventually, an infected person’s CD4+
                             T-cell will become so low,
                             < 200 cells/mm3 (a normal count is roughly 1,000 cells/mm3), that they will be
                             classified as having AIDS [8]. Acquired Immunodeficiency Syndrome is the final
                             stage of HIV. The patient’s immune system is severely damaged and can no longer
                             fight off infection. Without treatment, a patient can be expected to live for 3
                             years after an AIDS diagnosis [12].",
                            br(),
                            br(),
                            fluidRow(
                                column(7,
                                       h3(" HIV Viral Lifecycle"),
                                       "The life cycle of an
                                       HIV viral particle (either a virus or virion) is as follows: [14]",
                                       br(),
                                       tags$ol(type = "1",
                                               tags$li(" Unbound virus or virion circulates in the blood stream;"),
                                               tags$li(" Virus attaches itself to an uninfected cell;"),
                                               tags$li(" Virus empties contents (genetic code) into this cell;"),
                                               tags$li(" HIV RNA is used by the reverse transcriptase enzyme to build HIV DNA;"),
                                               tags$li(" HIV DNA is then inserted into the cell’s chromosome by the HIV
                                                        integrase enzyme, establishing HIV infection of the cell;"),
                                               tags$li(" The infected cell reproduces, activating HIV DNA which makes
                                                        new raw material for new viral particles (virions);"),
                                               tags$li(" Material packets for the new virus come together;"),
                                               tags$li(" The immature virus breaks free of infected cell and enters the blood stream;"),
                                               tags$li(" New virus matures and is now capable of infecting healthy cells."))),
                                column(5,
                                       img(src = "lifeCycle.png", height = 450, width = 600)))
                            )),
  ## explain antiretroviral
### make drop down panel explaining each different kind
               tabPanel("Treatment Options",
                        mainPanel(
                            h3("Background"),
                            "In 1987, the first drug was approved to treat HIV
                             [7]. AZT (full name azidothymidine), originally used
                             to treat cancer, is a reverse transcriptase inhibitor 
                             (RTI). It’s success was short lived as drug resistant
                             variants of HIV took over. Since then, combination
                             antiretroviral therapies have become the standard
                             of HIV care. There are 5 classes of antiretroviral
                             therapies (ART’s) used to treat HIV, table III outlines
                             the purpose of each drug class. Combination
                             antiretrovirals (cARV’s) typically include 3 or more
                             drugs from more than one class, in an attempt to
                             overcome HIV’s adaptive and drug resistant nature
                             [13].",
                            br(),
                            br(),
                            "3 enzymes vital to the productive infection of a
                             cell with HIV, were quickly identified as the reverse
                             transcriptase, integrase, and protease enzymes [7].
                             Reverse transcriptase takes viral RNA and turns it
                             into DNA, so that viral DNA can be joined with
                             the cell’s DNA, this is step 4 in section II-C [14].
                             The integrase enzyme completes the second half of
                             this process, by, as the name suggests, integrating
                             the viral DNA into the cell’s (step 5, II-C) [14]. The
                             protease enzyme has a different job. As the last step
                             in the process, protease finishes the maturation of an
                             HIV virion that has burst from the host cell [14]. The
                             identification of these three enzymes have provided
                             opportunities for targeted drug therapies.",
                            br(),
                            br(),
                            
                            
                        )),
               
               tabPanel("Model",
                        fluidRow(
                            column(3,
                                   h3("Model Modifications"),

                                   helpText("words to brief explain"),
                                   #checkboxGroupInput("drugTherapy",
                                   radioButtons("drugTherapy",
                                                      h4("Anti-retroviral Therapy"),
                                                      choices = c("No drug therapy" = 1,
                                                                  "Reverse Transcriptase Inhibitor (RTI)" = 2,
                                                                  "Protease Inhibitor (PI)" = 3,
                                                                  "Integrase Inhibitor" = 4),
                                                                  selected = 1),
                                   

                                   helpText("words to brief explain"),
                                   #need different time for ART?
                                   sliderInput("time",
                                               h4("Time (years)"),
                                               min = 0, max = 15, value = c(0,1)),

                                   #advanced: downward arrow: select plots
                                   #default is all 5, can select individual
                                   #inquire about how to do 1,23,4,5 and pass info onto HIV.R as "CD4+" etc.
                                   #(not have to write out each plot name in selected)
                                   checkboxGroupInput("isolatedPlots",
                                                      h4("Select Plot"),
                                                      choices = c("CD4+ T-Cells",
                                                                  "Infected CD4+ T-Cells",
                                                                  "Viral Load",
                                                                  "CD8+ T-Cells",
                                                                  "Activated CD8+ T-Cells"),
                                                      selected = c("CD4+ T-Cells",
                                                                   "Infected CD4+ T-Cells",
                                                                   "Viral Load",
                                                                   "CD8+ T-Cells",
                                                                   "Activated CD8+ T-Cells"))),
                            column(6,
                                   plotOutput("shinyGraphs")
                                   #plotOutput("shinyRTI")
                                   ),
                            column(3,
                                   wellPanel(
                                       
                                       "Define",
                                       "T-cell",
                                       "Viral load",
                                       "Activated t-cells",
                                       "etc."
                                   ))
                        )
                        ),


               tabPanel("Resources",
                        mainPanel(
                        )),
               tabPanel("About the Model",
                        mainPanel(
                            "summary, mission, etc.",
                            br(),
                            "Full paper",
                            br(),
                            "repo"))
               

               ))

server <- function(input, output, session) {
   # out_plot <- reactive({
                
    #    shinyPlot(choicesVar = input$isolatedPlots)}#, choicesDT = input$drugTherapy)
#       shinyRTI(choicesVar = input$isolatedPlots)
                                        # })
        # output$shinyRTI <- renderPlot({
     #       out_plot()
                                        # })
    out_plot_rti <- reactive({
       shinyRTI(choicesVar = input$isolatedPlots)
    })

    out_plot_hiv <- reactive ({
        shinyPlot(choicesVar = input$isolatedPlots)
    })

    output$shinyGraphs <- renderPlot({
        #out_plot_hiv()
        
        if(input$drugTherapy == 2)
        {
            out_plot_rti()
        }else{
            out_plot_hiv()
        }
    })
    
        
    
#    out_plot <- reactive({
        
        #shinyPlot(choicesVar = input$isolatedPlots)
                  
#         if(input$drugTherapy == "Reverse Transcriptase Inhibitor (RTI)")
#         {
#            print("me too")
#            shinyRTI(choicesVar = input$isolatedPlots)
#         }else{
#             print("im also running")
             #shinyRTI(choicesVar = NA)
#            shinyPlot(choicesVar = input$isolatedPlots)
#         }
#     })
    
#    output$shinyGraphs <- renderPlot({
#        out_plot()
#        })
    
}

shinyApp(ui = ui, server = server)
