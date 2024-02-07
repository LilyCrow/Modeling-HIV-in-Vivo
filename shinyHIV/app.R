library(shiny)

source("HIV.R")
source("RTI.R")

###Add dashed red line on CD4+ to mark AIDS threshold
###summary includes selected options, if 1 year selected "demonstrates withing host dynamics 365 post infection"
###make option to view cell counts side by side-- only if scale allows i.e. cd4+ and infected cd4+ on same plot
###scroll down to learn more about/ view intracellular delay, plots of drug concentration, two compartment (blood and lymphoid tissue) CD8+ focused models
#####Do i need a vocab section? Where should I put it, new panel? or elsewhere within an existing panel

#####when done, ask Sruthi/ nichols etc. to give to their students, do work sheet, provide feedback

ui <- fluidPage(

  
    navbarPage("Modeling the Within-Host Dynamics of HIV",

               tabPanel("How HIV Works",
                        titlePanel("How HIV Works: Between & Within Host(s) Dynamics"),
                        
                        mainPanel(
                            h3("Background"),
                            "Acquired Immunodeficiency Syndrome, also known as AIDS, was recognized
                             as a new disease in 1981 after a series of young, previously healthy
                             men died from unusual, and in many cases, rare, infection and malignancies [20], [7], [4].
                             On April 23, 1984 Human Immunodeficiency Virus Type 1 was identified as the virus that causes
                             AIDS [20] https://www.hiv.gov/hiv-basics/overview/history/hiv-and-aids-timeline/#year-1984.
                             Since then, at least 85.6 million people have tested positive for HIV
                             and 40.4 million have died [23].",
                            br(),
                            br(),
                            #"As a sexually transmitted disease",
                            "HIV is spread, primarily, through sex. Additionally, HIV can be contracted through
                             contact with an infect persons bodily fluids like blood, etc. HIV", em("cannot"), "be
                             contracted through contact with saliva, etc. 
                             In addition to sex, transmission most commonly occurs through contact with an open wound,
                             perinatally (during pregnancy), or from sharing needles [20].",
                            br(),
                            br(),
                            "It is important to acknowledge common misconceptions about HIV, and other STI's/STD's in
                             general. With the rise of the HIV epidemic, homophobic rhetoric and assumptions led to
                             dangerous behaviors that impacted both the LGBT+ community and the population at large.
                             While HIV was initially prominent among the LGBT+ community, due to a multitude of factors,
                             the virus was equally dangerous for others. Anti-LGBTQ thinking, promoted by interest
                             groups, schools, and then president, Ronald Reagan, led to heightened cases of HIV/ AIDS
                             as inaccurate and not scientific information was spread. No matter your sexual orientation,
                             you can still be affected by HIV, just as you can be by any other disease or infection.",
                            br(),
                            br(),
                            h3("How HIV works"),
                            "HIV results in a chronic, progressive disease that has
                             no cure. HIV is hallmarked by its high viral reproduction rate and
                             destruction of the immune system. The virus targets CD4+ T-cells, a type
                             of lymphocyte. CD4+ T-cells are responsible for detecting attacks on the
                             bodies own cells and sending signals to other 'killer' lymphocytes
                             triggering an immune response.",
                            br(),
                            #move (lymphocytes are white blood cells responsible for riding
                            #the body of invading antigens and infected cells) to vocabulary 
                            br(),
                            "Overtime HIV reduces the CD4+ T-cell population to roughly 20% of normal
                             levels (normal: ~1,000 cells/mm3 final stage of hiv ~200 cells/mm3) [19].
                             Such a low T-cell count leaves the body susceptible to opportunistic infections
                             and cancers. The near total destruction of the immune system, perpetrated by HIV,
                             makes treating such infections exceptionally difficult and survival rates plummit.",
                            br(),
                            br(),
                            h3("Stages of HIV"),
                            #br(),
                            "HIV has 3 distinct stages of infection: Acute HIV Infection, Chronic HIV
                             Infection, and Aquired Immune Deficiency Syndrome. The state of infection
                             is determined based upon a patient’s CD4+ T-cell count, viral load, and
                             presenting symptoms.",
                            br(),
                                        #make stages more organized:?
                                        #name
                                        #duration and time start (timeline)
                                        #cd4+ tcell count
                                        #viral load
                                        #presenting symptoms
                            br(),
                            "Acute HIV Infection: begins 2-4 weeks after the initial infection, lasting
                             for ~8-10 weeks [12], [8]. Acute infection is hallmarked by
                             flu-like symptoms that are often accompanied by rash, persistent cough,
                             and swollen lymph nodes. In this time the viral load is very high, as is
                             the likelihood of transmission [12]. During acute infection the patients
                             CD4+ T-cell count will drop dramatically[8]. It is at this time that the
                             patient will benefit most from beginning antiretroviral treatment.",
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
                                       h3("Understanding the HIV Viral Lifecycle"),
                                       "When a virus enters the body, the first thing it has to do is begin reproducing.
                                        How a virus replicates looks different for every virus, but for HIV, it looks like
                                        the following: [14]",
                                       br(),
                                       tags$ol(type = "1",
                                               tags$li("An unbound, or free floating, HIV viral particle, sometimes referred
                                                        to as a virion, circulates in the blood stream;"),
                                               tags$li("The unbound virus attaches itself to an uninfected CD4+ T-cell;"),
                                                        #, the cell is now infected
                                               tags$li("The virus empties its contents— genetic code in the form of RNA—
                                                        into the uninfected cell;"), #, the cell is now infected
                                               tags$li("An enzyme called reverse transcriptase takes the HIV RNA and
                                                        turns it into DNA, preparing the viral genetic code to be built
                                                        into the cells DNA;"), #\emph{cells}?
                                               tags$li("The integrase enzyme inserts the HIV DNA into the cell’s chromosome,
                                                        establishing HIV infection in the cell;"),
                                               tags$li("The infected cell reproduces, activating HIV DNA which makes
                                                        material for new viral particles (virions), instead of new
                                                        T-cells;"),
                                               tags$li("Material packets for the new virus come together;"),
                                               tags$li("The immature viral particles break free of the infected cell
                                                        and enter the blood stream;"),
                                               tags$li("The new viral particles mature and are now capable of infecting
                                                        healthy cells."))),
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
