##########################################
# Isotopic enrichment calculator 
# Author: Juris Meija, NRC Canada
# Date: March-April 2018

require(ecipex)
require(CHNOSZ)
require(shiny)
require(shinyjs)
require(rhandsontable)
require(mvtnorm)

#Isotope abundance dataframe (TICE-2013)
x.data <- read.csv("data/tice-2013-dataframe.csv")
#Isotope abundance covariance dataframe (TICE-2013) NOT YET IMPLEMENTED
x.cov.data <- read.csv("data/tice-2013-dataframe-covariances.csv")

# Helper functions to extract (1) nuclide masses, (2) isotopic abundances and (3) abundance covariances
mass <- function(el)   {x.data[which(x.data$element==el),]$mass}
mass.u <- function(el) {x.data[which(x.data$element==el),]$'mass.6u'}
x.mean <- function(el) {x.data[which(x.data$element==el),]$abundance}
x.cov <- function(el) {
  cov = matrix(c(2*sin((pi/6)*x.cov.data[which(x.cov.data$element==el),]$abundance.covariance)), nrow=length(x.mean(el)), ncol=length(x.mean(el)))
  (cov + t(cov))/2 
}

# consolidate isotopologues that are within the specified mass window
compact = function(z, tol=0.1){
  mask = rep(0, length(z$mass))
  k = 0
  for (i in 1:length(z$mass)) { k = k+1; mask[abs(z$mass[i]-z$mass) <= 2*tol] = k }
  d = data.frame(cbind('mask'=mask, 'm'=z$mass, 'x'=z$abundance))
  # abundance-weighted mean for mass
  m.aggregate = by(d, d$mask, function(x) weighted.mean(x$m, x$x))
  # sum isotopic abundances
  x.aggregate = by(d, d$mask, function(x) sum(x$x))
  list(m=m.aggregate, x=x.aggregate)
}

# isotope pattern fit without considering the instrumental mass bias
  ipfit0 = function(p, element, formula, charge, mtolerance, user.data, diso){
   diso$abundance[diso$element==element] <- c(1-p[1], p[1])
   z=ecipex(formula, isoinfo = diso, sortby = "mass")[[1]]
   #### convert masses to mass-to-charge ratios
   z$mass <- (z$mass - charge*0.0005486)/abs(charge)
   #### compact the mass spectrum 
   df=compact(z, tol=mtolerance)
   
   #### join the experimental and theoretical spectra
   m.j = outer(df$m, user.data$mass, function(x, y) abs(x-y))
   matches = which(m.j < mtolerance, arr.ind=TRUE)
   res = data.frame(model.m = df$m[matches[,1]], 
                    model.x = df$x[matches[,1]], 
                    data.m = user.data$mass[matches[,2]], 
                    data.x = user.data$abundance[matches[,2]])
   #### least squares fit 
   ifelse(length(res$data.x)<2, 1e35, c(sum(lm( res$data.x ~ 0 + res$model.x )$residuals^2)))
 }

## TEST DATA of microcystin [Asp3]MC-RR (C48H73N13O12) with x(15N) = 0.90 (z=+2)
   init.df = data.frame(mass=c(517.2693, 517.7680, 518.2666, 518.7652, 519.2642),
                       abundance=c(11698799,52250179,150346678,278752635,293253478))
  
### server file
server <- function(input, output, session) {
  
  output$A <- renderText(paste0("Isotopic abundance of ", c('carbon', 'nitrogen')[c(input$element=='C', input$element=='N')],"-",c(13, 15)[c(input$element=='C', input$element=='N')]))

  dfiso <- x.data

  inputformula <- reactive({
    ch = as.double(substr(input$charge,1,3))
    z = makeup(input$formula)
    z[which(names(z)=="H")] <- z[which(names(z)=="H")] + ch
    list(f=paste0(names(z),z, collapse=""), 
         ch=ch, fo=substr(input$charge, 3, nchar(input$charge)))
  })
  
  # standard molar mass
  s <- reactive({ q=ecipex(input$formula, isoinfo = x.data, limit = 1e-6)[[1]]; q[,1]%*%q[,2] })
  output$text0 <- renderText({ paste('Standard molar mass of M (natural isotopic composition):', formatC(s(), digits=2, format="f"), " g/mol") })

  # hot data input table (enrichment calculator)
  output$hot <- renderRHandsontable({
    if (is.null(input$hot)) { DF = init.df } else { DF = hot_to_r(input$hot) }
    rhandsontable(DF, readOnly = FALSE, stretchH = "all", selectCallback = TRUE) %>%
      hot_context_menu(allowColEdit = FALSE ) %>%
      hot_validate_numeric(cols = 1:2, min = 0) %>%
      hot_cols('float', format = '0.00')
  })

  user.data  <- reactive({
    if (is.null(input$hot)) { z = init.df } else { z = hot_to_r(input$hot) }
    list(mass=z[,1], abundance=z[,2])
  })
  
  v <- reactiveValues(enr = NULL, sd = NULL, all = NULL, M = NULL)
  
  observeEvent(input$button, {
            
              if(T){
               v$enr <- optimise(ipfit0, interval = c(0, 1), 
                        tol = 1e-10,
                        element = input$element, 
                        formula = inputformula()$f,
                        charge = inputformula()$ch, 
                        mtolerance = input$tolerance,
                        diso = x.data,
                        user.data = user.data() )$minimum
               v$sd <- NULL 
               }
             
            diso = x.data
            diso$abundance[diso$element==input$element] <- c(1-v$enr, v$enr)
            zz = ecipex(inputformula()$f, isoinfo = diso, sortby = "mass", limit=1e-6)[[1]]
            v$M <- zz[,1] %*% zz[,2]
            
            # update input$abundance slider
            updateSliderInput(session, "abundance", "Explore data by-hand", value = v$enr, min = 0, max = 1, step = 0.01)
            
            # error : no matching peaks
            #diso = x.data
            #diso$abundance[diso$element==input$element] = c(1-v$enr, v$enr)
            z=ecipex(inputformula()$f, isoinfo = diso, sortby = "mass")[[1]]
            z$mass <- (z$mass - inputformula()$ch*0.0005486)/abs(inputformula()$ch)
            df=compact(z, tol=input$tolerance)
            m.j = outer(df$m, user.data()$mass, function(x, y) abs(x-y))
            matches = which(m.j < input$tolerance, arr.ind=TRUE)
            res = data.frame(model.m = df$m[matches[,1]], 
                             model.x = df$x[matches[,1]], 
                             data.m = user.data()$mass[matches[,2]], 
                             data.x = user.data()$abundance[matches[,2]])
            
            if(length(res$data.x)<2){
              output$text2  <- renderText({ paste("<font color=\"#FF0000\">", "NO MATCHING PEAKS", "</font>") })
              output$text3  <- renderText({ paste("<font color=\"#FF0000\">", "Check the input parameters or try increasing the m/z tolerance", "</font>") })
            }
            
            if(length(res$data.x)>=2){
            # refresh plot upon the input$button
            output$pdf <- renderPlot({ 
              #diso = x.data
              #diso$abundance[diso$element==input$element] <- c(1-v$enr, v$enr)
              z=ecipex(inputformula()$f, isoinfo = diso, sortby = "mass", limit=1e-5)[[1]]
              #### convert masses to mass-to-charge ratios
              z$mass <- (z$mass - inputformula()$ch*0.0005486)/abs(inputformula()$ch)
              #### compact the mass spectrum 
              df=compact(z, tol=input$tolerance)
              #### rescale the two spectra (model and experimental)
              q=which(user.data()$abundance==max(user.data()$abundance))
              mass.max = user.data()$mass[q][1]
              slope = user.data()$abundance[q][1] / df$x[ abs(df$m - mass.max)==min(abs(df$m - mass.max)) ] [1]
              
              xy = user.data()
              max.x = max(xy$mass, df$m)
              min.x = max(0, min(xy$mass, df$m))
              max.y = max(xy$abundance, df$x*slope)
              par(mar=c(4.5, 4.1, 1.5, 3.1))
              plot(xy$mass, xy$abundance, type="h", main="", xlim=c(min.x, max.x), ylim=c(0, max.y), 
                   xlab="Mass-to-charge ratio", ylab="Abundance", cex.axis=1.3, cex.lab=1.3)
              rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "gray95")
              segments(df$m, 0, df$m, df$x*slope, lwd=0.75, lty=2, col='red')
              segments(xy$mass, 0, xy$mass, xy$abundance, lwd=3.5)
              points(df$m, df$x*slope, 
                     pch=21, col="red", bg="red", cex=1.3)
              box()
              })
            
            output$text2  <- renderText({ paste0("Best fit for the isotopic abundance of ",
                                                 c(13,15)[c(input$element=='C', input$element=='N')],
                                                 input$element,": ", formatC(v$enr, digits=3, format="f")," mol/mol") })
            #output$text2.u  <- renderText({ paste0(v$all) })
            output$text3  <- renderText({ paste0("Best fit for the molar mass of M: ",
                                                        formatC(v$M, digits=2, format="f")," g/mol") })
            }
                })
  
  observeEvent(input$abundance, {
    
    diso = x.data
    diso$abundance[diso$element==input$element] <- c(1-input$abundance, input$abundance)
    zz = ecipex(inputformula()$f, isoinfo = diso, sortby = "mass", limit=1e-6)[[1]]
    v$M <- zz[,1] %*% zz[,2]
    
    output$pdf <- renderPlot({ 
      diso = x.data
      diso$abundance[diso$element==input$element] <- c(1-input$abundance, input$abundance)
      z=ecipex(inputformula()$f, isoinfo = diso, sortby = "mass", limit=1e-6)[[1]]
      #### convert masses to mass-to-charge ratios
      z$mass <- (z$mass - inputformula()$ch*0.0005486)/abs(inputformula()$ch)
      #### compact the mass spectrum 
      df=compact(z, tol=input$tolerance)
      #### rescale the two spectra (model and experimental)
      q=which(user.data()$abundance==max(user.data()$abundance))
      mass.max = user.data()$mass[q][1]
      slope = user.data()$abundance[q][1] / df$x[ abs(df$m - mass.max)==min(abs(df$m - mass.max)) ] [1]
      
      xy = user.data()
      max.x = max(xy$mass, df$m)
      min.x = max(0, min(xy$mass, df$m))
      max.y = max(xy$abundance, df$x*slope)
      par(mar=c(4.5, 4.1, 1.5, 3.1))
      plot(xy$mass, xy$abundance, type="h", main="", xlim=c(min.x, max.x), ylim=c(0, max.y), 
           xlab="Mass-to-charge ratio", ylab="Abundance", cex=1.8, cex.axis=1.3, cex.lab=1.3)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "gray95")
      segments(df$m, 0, df$m, df$x*slope, lwd=0.75, lty=2, col='red')
      segments(xy$mass, 0, xy$mass, xy$abundance, lwd=3.5)
      points(df$m, df$x*slope, 
             pch=21, col="red", bg="red", cex=1.8)
      box()
    })
    
  })
  
  output$text1   <- renderText({ paste0("Neutral substance, M: ", input$formula) })
  output$text1a  <- renderText({ paste0("Measured ion, ", inputformula()$fo,": ", inputformula()$f, ' (charge = ',ifelse(inputformula()$ch>0,paste0("+",inputformula()$ch),inputformula()$ch),")") })

  shinyjs::onclick("toggleextra", shinyjs::toggle(id = "filterextra", anim = TRUE))
  shinyjs::onclick("togglenotes", shinyjs::toggle(id = "filternotes", anim = TRUE))
  shinyjs::onclick("toggleextra.mol", shinyjs::toggle(id = "filterextra", anim = TRUE))
  shinyjs::onclick("togglenotes.mol", shinyjs::toggle(id = "filternotes", anim = TRUE))
  
}

### user interface
ui <- fluidPage(

  titlePanel( title="Isotopic enrichment calculator" ),
  shinyjs::useShinyjs(),
  sidebarLayout(
    
    sidebarPanel(
      textInput("formula", label = "Elemental composition of the (neutral) substance", value = "C48H73N13O12", width = '85%'),
      selectInput("element", label = "Select the isotopic element", choices = c('carbon'='C','nitrogen'='N'), selected = 'N', width = '85%'),
      selectInput("charge", label = "Select the charge-state", choices=c("-2 [M-2H]2-","-1 [M-H]-","+1 [M+H]+","+2 [M+2H]2+"), selected="+2 [M+2H]2+", width = '85%'),
      br(),
      wellPanel(tags$style(type="text/css", '.well {width: 85%}'),
        sliderInput("abundance", label = "Explore data by-hand", min=0, max=1, value=0.9, round=TRUE, step=0.01, width='100%'),
        helpText(textOutput("A")), align="center"),
        
      h5("Additional parameters", a(id = "toggleextra", "show/hide")),
      shinyjs::hidden(div(id = "filterextra",
                          fluidRow(
                            column(12, numericInput("tolerance", label = "Enter m/z tolerance", value = 0.1, max = 0.5, min=0.0001, step = 0.001, width = '85%'))
                            #column(12, radioButtons("unc", "Should the uncertainty about the natural isotopic composition be considered?", choices=c("Yes","No"), selected="No", inline = TRUE, width = '85%'))
      ))),
      br(),
      h5(tags$b("Enter (paste) the observed isotopic pattern")),
      rHandsontableOutput("hot"),
      helpText("right-click to add or delete rows"),
      br(),
      conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",   
                         actionButton("button", label = "Perform least squares fitting", icon = icon('bar-chart-o'))),
      conditionalPanel(condition = "$('html').hasClass('shiny-busy')",   
                         actionButton("button", label = "busy...", icon = icon('hourglass-half')))      
    ),
    
    mainPanel(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"),
      fluidRow(column(helpText("This calculator evaluates the experimental mass spectra for a given substance and determines the isotopic composition of a selected element that best matches the experimental data."),width=11)),
      br(),textOutput("text1"),textOutput("text1a"),textOutput("text0"),
      br(),
      fluidRow(column(plotOutput(("pdf")),width=11)),
      conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",   
                       strong(htmlOutput("text2")),
                       strong(htmlOutput("text2.u")),
                       strong(htmlOutput("text3")),
                       htmlOutput("text4"),
                       strong(htmlOutput("textU")),
                       htmlOutput("textUU")
                       ),
      br(),
      h5("Notes and explanations", a(id = "togglenotes", "show/hide")),
      shinyjs::hidden(div(id = "filternotes",
                          fluidRow(column(
                            p("Black lines: experimental data; red lines: theoretical model."),
                            p("The calculator performs simple fitting of the isotope patterns and, at this stage, does not consider the uncerrtainty of the natural isotopic composition of the elements nor any instrumental isotope ratio fractionation effects during the mass spectrometry measurements."),
                            p("Created using R and Shiny. Juris Meija (2019) NRC Canada"),width=11)
                          ))),
      br(),
      p("NRC Isotopic Enrichment Calculator (2019) v.1.8")
      
    )
  )
)

shinyApp(ui = ui, server = server)
