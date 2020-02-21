##########################################
# Isotopic enrichment calculator 
# Author: Juris Meija, NRC Canada
# Date: 2019-2020, version 2.11 (Feb 2020)
##########################################

# Brief description
# This calculator obtains the best fit of the isotopic composition of nitrogen/carbon/oxygen
# from an observed mass spectrum of the molecule with known identity (molecular formula)

require(ecipex)
require(CHNOSZ)
require(shiny)
require(shinyjs)
require(rhandsontable)
#require(mvtnorm)
#setwd('F:/SHINY-METROLOGY/isotopic-enrichment-calculator')
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
  if(element=='H') {
    diso$abundance[diso$element=='Hx'] <- c(1-p[1], p[1])
    diso$abundance[diso$element=='Hy'] <- c(1-p[2], p[2])
    }
  if(element=='C') {
    diso$abundance[diso$element=='Cx'] <- c(1-p[1], p[1])
    diso$abundance[diso$element=='Cy'] <- c(1-p[2], p[2])
  }
  if(element=='N') {
    diso$abundance[diso$element=='Nx'] <- c(1-p[1], p[1])
    diso$abundance[diso$element=='Ny'] <- c(1-p[2], p[2])
  }
  if(element=='O') {
    diso$abundance[diso$element=='Ox'] <- c(1-p[1], 0, p[1])
    diso$abundance[diso$element=='Oy'] <- c(1-p[2], 0, p[2])
  }
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

## TEST DATA of microcystin [Asp3]MC-RR (C48H73N13O12) with x(15N) = 0.98 (z=+2)
init.df = data.frame(mass=c(518.264, 518.763, 519.261, 519.763, 520.264),
                     abundance=c(2534313,14690051,58837471,28533149,8711512))
#init.df = data.frame(mass=c(828.4899,829.4936,830.4936,831.4966,832.4989,833.5016,834.5037,835.5064),
#                     abundance=c(7.06,4.7,100,51.65,22.42,6.82,1.79,0.48))

### server file
server <- function(input, output, session) {
  
  # Text for input$abundance slider
  output$A0 <- renderText(paste0("Isotopic abundance of ", c('hydrogen-2','carbon-13', 'nitrogen-15', 'oxygen-18')[c(input$element=='H', input$element=='C', input$element=='N', input$element=='O')]," (pool 1)"))
  output$A1 <- renderText(paste0("Isotopic abundance of ", c('hydrogen-2','carbon-13', 'nitrogen-15', 'oxygen-18')[c(input$element=='H', input$element=='C', input$element=='N', input$element=='O')]," (pool 2)"))
  
  dfiso <- x.data
  
  # Parse the chemical formula from input parameters and data
  inputformula <- reactive({
    #input=list(formula='C48H73N13O12',charge=2,adduct='H2',element='N',atoms1=1,atoms2=12)
    z = makeup(input$formula)
    
    # add placeholders for C, H, N, O and their enriched counterparts Cx Hx Nx Ox
    if(all(names(z)!='C')) z = c(z,C=0)
    if(all(names(z)!='H')) z = c(z,H=0)
    if(all(names(z)!='N')) z = c(z,N=0)
    if(all(names(z)!='O')) z = c(z,O=0)
    
    if(input$element=='H') z = c(z, 'Hx'=input$atoms1, 'Hy'=input$atoms2)
    if(input$element=='C') z = c(z, 'Cx'=input$atoms1, 'Cy'=input$atoms2)
    if(input$element=='N') z = c(z, 'Nx'=input$atoms1, 'Ny'=input$atoms2)
    if(input$element=='O') z = c(z, 'Ox'=input$atoms1, 'Oy'=input$atoms2)
    
    # readjust the total number of atoms
    if(input$element=='H') z[names(z)=='H'] <- z[names(z)=='H'] - z[names(z)=='Hx'] - z[names(z)=='Hy']
    if(input$element=='C') z[names(z)=='C'] <- z[names(z)=='C'] - z[names(z)=='Cx'] - z[names(z)=='Cy']
    if(input$element=='N') z[names(z)=='N'] <- z[names(z)=='N'] - z[names(z)=='Nx'] - z[names(z)=='Ny']
    if(input$element=='O') z[names(z)=='O'] <- z[names(z)=='O'] - z[names(z)=='Ox'] - z[names(z)=='Oy']
    
    # combine formula with the adduct while considering the charge (add positive, subtract negative)
    # validate input$adduct
      an.error.occured <- FALSE
      tryCatch( { mm = makeup(input$adduct) }, error = function(e) {an.error.occured <<- TRUE} )
      if(an.error.occured) mm = makeup('H0')
    
    q=merge(z,mm*ifelse(input$charge<0,-1,+1),by.x="row.names", by.y="row.names", all.x=TRUE, all.y=TRUE)
    q[is.na(q)]<-0
    
    # adjust the total number of atoms for the presence of secified enriched atoms
    #q[q$'Row.names'==input$element,'x'] <- q[q$'Row.names'==input$element,'x'] - q[q$'Row.names'==paste0(input$element,'x'),'x']
    # combine
    q[,'total']<-q[,'x']+q[,'y']
    q <- q[q[,'total']!=0,]
    
    q.a = cbind(q[q[,'y']!=0,c('Row.names')],abs(q[q[,'y']!=0,c('y')]))
    
    foo.t.text=gsub('Oy','[<sup>**</sup>O]',gsub('Hy','[<sup>**</sup>H]',gsub('Ny','[<sup>**</sup>N]',gsub('Cy','[<sup>**</sup>C]',gsub('Ox','[<sup>*</sup>O]',gsub('Hx','[<sup>*</sup>H]',gsub('Nx','[<sup>*</sup>N]',gsub('Cx','[<sup>*</sup>C]', 
                          paste0(q[,'Row.names'],'<sub>',q[,'total'],'</sub>',collapse='') ))))))))
    foo.m.text=gsub('Oy','[<sup>**</sup>O]',gsub('Hy','[<sup>**</sup>H]',gsub('Ny','[<sup>**</sup>N]',gsub('Cy','[<sup>**</sup>C]',gsub('Ox','[<sup>*</sup>O]',gsub('Hx','[<sup>*</sup>H]',gsub('Nx','[<sup>*</sup>N]',gsub('Cx','[<sup>*</sup>C]', 
                          paste0(q[,'Row.names'],'<sub>',q[,'x'],'</sub>',collapse='') ))))))))
    
    list(f=paste0(q[,'Row.names'], q[,'total'], collapse=""), 
         fo=paste0('[M',c('-','+')[c(input$charge<0,input$charge>0)],input$adduct,']',input$charge,c('-','+')[c(input$charge<0,input$charge>0)]),
         foo.t=foo.t.text,
         foo.m=foo.m.text,
         foo.a=paste0(q.a[,1],'<sub>',q.a[,2],'</sub>',collapse=''),
         ch=as.double(input$charge),
         h.max=sum(z[which(names(z) %in% c('H','Hx','Hy'))]),
         c.max=sum(z[which(names(z) %in% c('C','Cx','Cy'))]),
         n.max=sum(z[which(names(z) %in% c('N','Nx','Ny'))]),
         o.max=sum(z[which(names(z) %in% c('O','Ox','Oy'))]),
         enr.atoms=input$atoms1+input$atoms2
         )
  })
  
  # update input$atoms1
  observe({
    updateSliderInput(session, "atoms1", max = unname(c(inputformula()$h.max, inputformula()$c.max, inputformula()$n.max, inputformula()$o.max)[c('H','C','N','O')==input$element]) - input$atoms2)
  })
  # update input$atoms2
  observe({
    updateSliderInput(session, "atoms2", max = unname(c(inputformula()$h.max, inputformula()$c.max, inputformula()$n.max, inputformula()$o.max)[c('H','C','N','O')==input$element]) - input$atoms1)
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
      hot_cols('float', format = '0.000')
  })
  
  user.data  <- reactive({
    if (is.null(input$hot)) { z = init.df } else { z = hot_to_r(input$hot) }
    list(mass=z[,1], abundance=z[,2])
  })
  
  v <- reactiveValues(enr = NULL, sd = NULL, all = NULL, M = NULL, atoms = NULL, atoms.max = NULL)
  
  observeEvent(input$element, {
    
    v$atoms.max <- unname(c(inputformula()$h.max, inputformula()$c.max, inputformula()$n.max, inputformula()$o.max)[c('H','C','N','O')==input$element])
    
    # update input$atoms
    updateSliderInput(session, "atoms1", max = v$atoms.max - input$atoms2)
    updateSliderInput(session, "atoms2", max = v$atoms.max - input$atoms1)
    
  })
  
  observeEvent(input$button, {
    
      v$atoms.max <- unname(c(inputformula()$h.max, inputformula()$c.max, inputformula()$n.max, inputformula()$o.max)[c('H','C','N','O')==input$element])
      # update input$atoms
      updateSliderInput(session, "atoms1", max = v$atoms.max)
      updateSliderInput(session, "atoms2", max = v$atoms.max)
      
      v$enr <- optim(c(0.95, 0.95), ipfit0, method = "L-BFGS-B", 
                        lower= c(0, 0), upper = c(1, 1), 
                        control = list(factr = 1e-10),
                        element = input$element, 
                        formula = inputformula()$f,
                        charge = inputformula()$ch, 
                        mtolerance = input$tolerance,
                        diso = x.data,
                        user.data = user.data() )$par
      v$sd <- NULL 
  
      diso = x.data
      if(input$element=='H') {
        diso$abundance[diso$element=='Hx'] <- c(1-v$enr[1], v$enr[1])
        diso$abundance[diso$element=='Hy'] <- c(1-v$enr[2], v$enr[2])
      }
      if(input$element=='C') {
        diso$abundance[diso$element=='Cx'] <- c(1-v$enr[1], v$enr[1])
        diso$abundance[diso$element=='Cy'] <- c(1-v$enr[2], v$enr[2])
      }
      if(input$element=='N') {
        diso$abundance[diso$element=='Nx'] <- c(1-v$enr[1], v$enr[1])
        diso$abundance[diso$element=='Ny'] <- c(1-v$enr[2], v$enr[2])
      }
      if(input$element=='O') {
        diso$abundance[diso$element=='Ox'] <- c(1-v$enr[1], 0, v$enr[1])
        diso$abundance[diso$element=='Oy'] <- c(1-v$enr[2], 0, v$enr[2])
      }

    zz = ecipex(inputformula()$f, isoinfo = diso, sortby = "mass", limit=1e-6)[[1]]
    v$M <- zz[,1] %*% zz[,2]
    
    # update input$abundance sliders
    updateSliderInput(session, "abundance1", value = v$enr[1])
    updateSliderInput(session, "abundance2", value = v$enr[2])
    
    # error : no matching peaks
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
      output$text2A  <- renderText({ paste("<font color=\"#FF0000\">", "NO MATCHING PEAKS", "</font>") })
      output$text2B  <- renderText({ paste("<font color=\"#FF0000\">", "Check the input parameters or try increasing the m/z tolerance", "</font>") })
    }
    
    if(length(res$data.x)>=2){
      # refresh plot upon the input$button
      output$pdf <- renderPlot({ 
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
        segments(df$m, 0, df$m, df$x*slope, lwd=4)
        points(x=xy$mass, y=xy$abundance, pch=21, cex=2, col='black', bg='white')
        box()
      })
      
      if(input$atoms2>0){
        output$text2A  <- renderText({ paste0("Best fit for the isotopic abundance of ",
                                              "<sup>", c(2,13,15,18)[c(input$element=='H',input$element=='C', input$element=='N', input$element=='O')],
                                              "</sup>", input$element," in ", 
                                              paste0('[<sup>*</sup>',input$element,']')," and ",
                                              paste0('[<sup>**</sup>',input$element,']')," : ",
                                              formatC(v$enr[1], digits=3, format="f")," mol/mol (pool 1) and ", formatC(v$enr[1], digits=3, format="f"), " mol/mol (pool 2)") })
      }
      if(input$atoms2==0){
        output$text2A  <- renderText({ paste0("Best fit for the isotopic abundance of ",
                                              "<sup>", c(2,13,15,18)[c(input$element=='H',input$element=='C', input$element=='N', input$element=='O')],
                                              "</sup>", input$element," in ", paste0('[<sup>*</sup>',input$element,']'),": ", formatC(v$enr[1], digits=3, format="f")," mol/mol") })    
      }
        
      output$text2B  <- renderText({ paste0("Best fit for the molar mass of M: ",
                                           formatC(v$M, digits=2, format="f")," g/mol") })
    }
  })
  
  observeEvent(input$abundance1, {
    
    diso = x.data
    if(input$element=='H') {
      diso$abundance[diso$element=='Hx'] <- c(1-input$abundance1, input$abundance1)
      diso$abundance[diso$element=='Hy'] <- c(1-input$abundance2, input$abundance2)
      }
    if(input$element=='C') {
      diso$abundance[diso$element=='Cx'] <- c(1-input$abundance1, input$abundance1)
      diso$abundance[diso$element=='Cy'] <- c(1-input$abundance2, input$abundance2)
    }
    if(input$element=='N') {
      diso$abundance[diso$element=='Nx'] <- c(1-input$abundance1, input$abundance1)
      diso$abundance[diso$element=='Ny'] <- c(1-input$abundance2, input$abundance2)
    }
    if(input$element=='O') {
      diso$abundance[diso$element=='Ox'] <- c(1-input$abundance1, 0, input$abundance1)
      diso$abundance[diso$element=='Oy'] <- c(1-input$abundance2, 0, input$abundance2) 
    }
    
    zz = ecipex(inputformula()$f, isoinfo = diso, sortby = "mass", limit=1e-6)[[1]]
    v$M <- zz[,1] %*% zz[,2]
    
    output$pdf <- renderPlot({ 
      diso = x.data
      if(input$element=='H') {
        diso$abundance[diso$element=='Hx'] <- c(1-input$abundance1, input$abundance1)
        diso$abundance[diso$element=='Hy'] <- c(1-input$abundance2, input$abundance2)
      }
      if(input$element=='C') {
        diso$abundance[diso$element=='Cx'] <- c(1-input$abundance1, input$abundance1)
        diso$abundance[diso$element=='Cy'] <- c(1-input$abundance2, input$abundance2)
      }
      if(input$element=='N') {
        diso$abundance[diso$element=='Nx'] <- c(1-input$abundance1, input$abundance1)
        diso$abundance[diso$element=='Ny'] <- c(1-input$abundance2, input$abundance2)
      }
      if(input$element=='O') {
        diso$abundance[diso$element=='Ox'] <- c(1-input$abundance1, 0, input$abundance1)
        diso$abundance[diso$element=='Oy'] <- c(1-input$abundance2, 0, input$abundance2)
      }
      z=ecipex(inputformula()$f, isoinfo = diso, sortby = "mass", limit=1e-6)[[1]]
      #### convert masses to mass-to-charge ratios
      z$mass <- (z$mass - as.double(inputformula()$ch)*0.0005486)/abs(as.double(inputformula()$ch))
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
      segments(df$m, 0, df$m, df$x*slope, lwd=3.5)
      points(x=xy$mass, y=xy$abundance, pch=19, cex=2)
      box()
    })
    
  })

  observeEvent(input$abundance2, {
    
    diso = x.data
    if(input$element=='H') {
      diso$abundance[diso$element=='Hx'] <- c(1-input$abundance1, input$abundance1)
      diso$abundance[diso$element=='Hy'] <- c(1-input$abundance2, input$abundance2)
    }
    if(input$element=='C') {
      diso$abundance[diso$element=='Cx'] <- c(1-input$abundance1, input$abundance1)
      diso$abundance[diso$element=='Cy'] <- c(1-input$abundance2, input$abundance2)
    }
    if(input$element=='N') {
      diso$abundance[diso$element=='Nx'] <- c(1-input$abundance1, input$abundance1)
      diso$abundance[diso$element=='Ny'] <- c(1-input$abundance2, input$abundance2)
    }
    if(input$element=='O') {
      diso$abundance[diso$element=='Ox'] <- c(1-input$abundance1, 0, input$abundance1)
      diso$abundance[diso$element=='Oy'] <- c(1-input$abundance2, 0, input$abundance2) 
    }
    
    zz = ecipex(inputformula()$f, isoinfo = diso, sortby = "mass", limit=1e-6)[[1]]
    v$M <- zz[,1] %*% zz[,2]
    
    output$pdf <- renderPlot({ 
      diso = x.data
      if(input$element=='H') {
        diso$abundance[diso$element=='Hx'] <- c(1-input$abundance1, input$abundance1)
        diso$abundance[diso$element=='Hy'] <- c(1-input$abundance2, input$abundance2)
      }
      if(input$element=='C') {
        diso$abundance[diso$element=='Cx'] <- c(1-input$abundance1, input$abundance1)
        diso$abundance[diso$element=='Cy'] <- c(1-input$abundance2, input$abundance2)
      }
      if(input$element=='N') {
        diso$abundance[diso$element=='Nx'] <- c(1-input$abundance1, input$abundance1)
        diso$abundance[diso$element=='Ny'] <- c(1-input$abundance2, input$abundance2)
      }
      if(input$element=='O') {
        diso$abundance[diso$element=='Ox'] <- c(1-input$abundance1, 0, input$abundance1)
        diso$abundance[diso$element=='Oy'] <- c(1-input$abundance2, 0, input$abundance2)
      }
      z=ecipex(inputformula()$f, isoinfo = diso, sortby = "mass", limit=1e-6)[[1]]
      #### convert masses to mass-to-charge ratios
      z$mass <- (z$mass - as.double(inputformula()$ch)*0.0005486)/abs(as.double(inputformula()$ch))
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
      segments(df$m, 0, df$m, df$x*slope, lwd=3.5)
      points(x=xy$mass, y=xy$abundance, pch=19, cex=2)
      box()
    })
    
  })
  
  output$text1A  <- renderText({ paste0("Molecular ion, M : ", inputformula()$foo.m) })
  output$text1B  <- renderText({ paste0("Measured ion, [M",ifelse(inputformula()$ch>0,'+','-'), inputformula()$foo.a,"]","<sup>",abs(inputformula()$ch),ifelse(inputformula()$ch>0,'+','-'),'</sup> : ',inputformula()$foo.t) })
  
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
      selectInput("element", label = "Isotopic element", choices = c('hydrogen-2'='H','carbon-13'='C','nitrogen-15'='N', 'oxygen-18'='O'), selected = 'N', width = '85%'),
      fluidRow(
        column(5, sliderInput("atoms1", label = "Number of isotopic atoms (pool 1)", min=1, max=13, value=13, round=TRUE, step = 1, ticks = TRUE)),
        column(5, sliderInput("atoms2", label = "Number of isotopic atoms (pool 2)", min=0, max=13, value=0, round=TRUE, step = 1, ticks = TRUE)),
        width = '85%'
      ),
      fluidRow(
        column(5, textInput("adduct", label = "Specify adduct", value="H2") ),
        column(5, selectInput("charge", label = "Charge-state", choices=c('-2'=-2,'-1'=-1,'+1'=1,'+2'=2), selected=2) ),
        width = '85%' ),
      h5("Additional settings", a(id = "toggleextra", "show/hide")),
      shinyjs::hidden(div(id = "filterextra",
                          fluidRow(
                            column(12, sliderInput("abundance1", label = "Explore data by-hand", min=0, max=1, value=0.98, round=TRUE, step=0.01, width='85%') ),
                            column(11, helpText(textOutput("A0"))),
                            column(12, sliderInput("abundance2", label = "", min=0, max=1, value=0.98, round=TRUE, step=0.01, width='85%') ),
                            column(11, helpText(textOutput("A1"))),
                            column(12, numericInput("tolerance", label = "Enter m/z tolerance", value = 0.1, max = 0.5, min=0.0001, step = 0.001, width = '85%'))
                            ))),
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
      fluidRow(column(helpText("This calculator evaluates the experimental mass spectra for a given substance and determines the isotopic composition of a selected element that best matches the experimental data."), width=11)),
      br(),
      strong(htmlOutput("text1A")),
      strong(htmlOutput("text1B")),
      br(),
      fluidRow(column(plotOutput(("pdf")),width=11)),
      conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",   
                       strong(htmlOutput("text2A")),
                       strong(htmlOutput("text2B")),
                       htmlOutput("text4"),
                       strong(htmlOutput("text1b")),
                       strong(htmlOutput("textU")),
                       htmlOutput("textUU")
      ),
      br(),
      h5("Notes and explanations", a(id = "togglenotes", "show/hide")),
      shinyjs::hidden(div(id = "filternotes",
                          fluidRow(column(
                            p("Circles: experimental data; Lines: theoretical model."),
                            p("The calculator performs simple fitting of the isotope patterns and, at this stage, does not consider the uncerrtainty of the natural isotopic composition of the elements nor any instrumental isotope ratio fractionation effects during the mass spectrometry measurements."),
                            p("Created using R and Shiny. Juris Meija (2019-2020) NRC Canada"),width=11)
                          ))),
      br(),
      p("NRC Isotopic Enrichment Calculator (2019-2020) v.2.11")
      
    )
  )
)

shinyApp(ui = ui, server = server)