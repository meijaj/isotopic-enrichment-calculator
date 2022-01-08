##########################################
# Isotopic enrichment calculator 
# Author: Juris Meija, NRC Canada
# Date: 2019-2021, version December 2021
##########################################

# Brief description
# This calculator obtains the best fit of the isotopic composition of nitrogen/carbon/oxygen
# from an observed mass spectrum of the molecule with known identity (molecular formula)

require(shiny)
require(shinyjs)       # Javascript for shiny
require(rhandsontable) # Hot tables in shiny

require(ecipex)    # FT-based calculation of molecular formulas
require(CHNOSZ)    # parsing chemical formulas
require(textclean) # mgsub function

require(optimx)    # multi-dimensional non-linear optimization
require(writexl)   # write xlsx files
require(tibble)    # write xlsx files

#Isotope abundance dataframe (TICE-2013)
x.data <- read.csv("data/tice-2013-dataframe.csv")
names(x.data)[2]<-'nucleons'
x.data <- x.data[,c(1,3,4,2)]

# Helper functions to extract (1) nuclide masses, (2) isotopic abundances and (3) abundance covariances
mass <- function(el)   {x.data[which(x.data$element==el),]$mass}
mass.u <- function(el) {x.data[which(x.data$element==el),]$'mass.6u'}
x.mean <- function(el) {x.data[which(x.data$element==el),]$abundance}

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
    diso$abundance[diso$element=='Hz'] <- c(1-p[3], p[3])
    diso$abundance[diso$element=='Hw'] <- c(1-p[4], p[4])
    }
  if(element=='C') {
    diso$abundance[diso$element=='Cx'] <- c(1-p[1], p[1])
    diso$abundance[diso$element=='Cy'] <- c(1-p[2], p[2])
    diso$abundance[diso$element=='Cz'] <- c(1-p[3], p[3])
    diso$abundance[diso$element=='Cw'] <- c(1-p[4], p[4])
  }
  if(element=='N') {
    diso$abundance[diso$element=='Nx'] <- c(1-p[1], p[1])
    diso$abundance[diso$element=='Ny'] <- c(1-p[2], p[2])
    diso$abundance[diso$element=='Nz'] <- c(1-p[3], p[3])
    diso$abundance[diso$element=='Nw'] <- c(1-p[4], p[4])
  }
  if(element=='O') {
    diso$abundance[diso$element=='Ox'] <- c(1-p[1], 0, p[1])
    diso$abundance[diso$element=='Oy'] <- c(1-p[2], 0, p[2])
    diso$abundance[diso$element=='Oz'] <- c(1-p[3], 0, p[3])
    diso$abundance[diso$element=='Ow'] <- c(1-p[4], 0, p[4])
  }
  z=ecipex(formula, isoinfo = diso, sortby = "mass", limit = 1e-13)[[1]]
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
                   data.x = user.data$intensity[matches[,2]])
  # Least squares fit 
  ifelse(length(res$data.x)<2, 1e35, c(sum(lm( res$data.x ~ 0 + res$model.x )$residuals^2)))
  
}

## TEST DATA of microcystin [Asp3]MC-RR (C48H73N13O12) with x(15N) = 0.98 (z=+2)
init.df = data.frame(mass=c(518.264, 518.763, 519.261, 519.763, 520.264),
                     intensity=c(2534313,14690051,58837471,28533149,8711512))

### server file
server <- function(input, output, session) {
  
  # Text for input$abundance sliders
  ia = c('hydrogen-2','carbon-13', 'nitrogen-15', 'oxygen-18')
  output$A1 <- renderText(paste0("Isotopic abundance of ", ia[c(input$element=='H', input$element=='C', input$element=='N', input$element=='O')], " (pool 1)"))
  output$A2 <- renderText(paste0("Isotopic abundance of ", ia[c(input$element=='H', input$element=='C', input$element=='N', input$element=='O')], " (pool 2)"))
  output$A3 <- renderText(paste0("Isotopic abundance of ", ia[c(input$element=='H', input$element=='C', input$element=='N', input$element=='O')], " (pool 3)"))
  output$A4 <- renderText(paste0("Isotopic abundance of ", ia[c(input$element=='H', input$element=='C', input$element=='N', input$element=='O')], " (pool 4)"))
  
  dfiso <- x.data
  
  # Parse the chemical formula from input parameters and data
  inputformula <- reactive({
  z = makeup(input$formula)
    
    # create placeholders for C, H, N, O (if needed)
    if(all(names(z)!='C')) z = c(z, C=0)
    if(all(names(z)!='H')) z = c(z, H=0)
    if(all(names(z)!='N')) z = c(z, N=0)
    if(all(names(z)!='O')) z = c(z, O=0)
    
    # create placeholders for the labeled elements Ex, Ey, Ez, Ew (E = H, C, N, or O)
    if(input$element=='H') z = c(z, 'Hx'=input$atoms1, 'Hy'=input$atoms2, 'Hz'=input$atoms3, 'Hw'=input$atoms4)
    if(input$element=='C') z = c(z, 'Cx'=input$atoms1, 'Cy'=input$atoms2, 'Cz'=input$atoms3, 'Cw'=input$atoms4)
    if(input$element=='N') z = c(z, 'Nx'=input$atoms1, 'Ny'=input$atoms2, 'Nz'=input$atoms3, 'Nw'=input$atoms4)
    if(input$element=='O') z = c(z, 'Ox'=input$atoms1, 'Oy'=input$atoms2, 'Oz'=input$atoms3, 'Ow'=input$atoms4)
    
    # readjust the total number of atoms
    if(input$element=='H') z[names(z)=='H'] <- z[names(z)=='H'] - z[names(z)=='Hx'] - z[names(z)=='Hy'] - z[names(z)=='Hz'] - z[names(z)=='Hw']
    if(input$element=='C') z[names(z)=='C'] <- z[names(z)=='C'] - z[names(z)=='Cx'] - z[names(z)=='Cy'] - z[names(z)=='Cz'] - z[names(z)=='Cw']
    if(input$element=='N') z[names(z)=='N'] <- z[names(z)=='N'] - z[names(z)=='Nx'] - z[names(z)=='Ny'] - z[names(z)=='Nz'] - z[names(z)=='Nw']
    if(input$element=='O') z[names(z)=='O'] <- z[names(z)=='O'] - z[names(z)=='Ox'] - z[names(z)=='Oy'] - z[names(z)=='Oz'] - z[names(z)=='Ow']
    
    # combine formula with the adduct while considering the charge (add positive, subtract negative)
    
    # validate input$adduct
      an.error.occured <- FALSE
      tryCatch( { mm = makeup(input$adduct) }, error = function(e) {an.error.occured <<- TRUE} )
      if(an.error.occured) mm = makeup('H0')
    
    q=merge(z,mm*ifelse(input$charge<0, -1, +1), by.x="row.names", by.y="row.names", all.x=TRUE, all.y=TRUE)
    q[is.na(q)]<-0
    
    # combine everything
    
    q[,'total'] <- q[,'x'] + q[,'y'] # x and y here refer to the two joined datasets 
    q <- q[q[,'total']!=0,]
    
    # adduct
    q.a = cbind(q[q[,'y']!=0, c('Row.names')], abs(q[q[,'y']!=0,c('y')]))
    
    # parse molecular formulas
    pp1 = c('Cx','Cy','Cz','Cw','Hx','Hy','Hz','Hw','Nx','Ny','Nz','Nw','Ox','Oy','Oz','Ow')
    pp2 = c('[<sup>*</sup>C]','[<sup>**</sup>C]','[<sup>***</sup>C]','[<sup>****</sup>C]','[<sup>*</sup>H]','[<sup>**</sup>H]','[<sup>***</sup>H]','[<sup>****</sup>H]',
            '[<sup>*</sup>N]','[<sup>**</sup>N]','[<sup>***</sup>N]','[<sup>****</sup>N]','[<sup>*</sup>O]','[<sup>**</sup>O]','[<sup>***</sup>O]','[<sup>****</sup>O]')
    
    # measured adduct
    foo.t.text=mgsub(paste0(q[,'Row.names'],'<sub>',q[,'total'],'</sub>',collapse=''), pp1, pp2)
      
    # molecular ion
    foo.m.text=mgsub(paste0(q[,'Row.names'],'<sub>',q[,'x'],'</sub>',collapse=''), pp1, pp2)
    
    x.natural = x.data[x.data[,'element']==input$element,][c(2,2,2,3)[input$element==c('H','C','N','O')],'abundance']
    
    result = list(f=paste0(q[,'Row.names'], q[,'total'], collapse=""), 
         fo=paste0('[M', c('-', '+')[c(input$charge<0, input$charge>0)], input$adduct,']', input$charge, c('-', '+')[c(input$charge<0, input$charge>0)]),
         foo.t=foo.t.text,
         foo.m=foo.m.text,
         foo.a=paste0(q.a[,1], '<sub>', q.a[,2], '</sub>', collapse=''),
         ch=as.double(input$charge),
         h.max=sum(z[which(names(z) %in% c('H','Hx','Hy','Hz','Hw'))]),
         c.max=sum(z[which(names(z) %in% c('C','Cx','Cy','Cz','Cw'))]),
         n.max=sum(z[which(names(z) %in% c('N','Nx','Ny','Nz','Nw'))]),
         o.max=sum(z[which(names(z) %in% c('O','Ox','Oy','Oz','Ow'))]),
         enr.atoms=input$atoms1 + input$atoms2 + input$atoms3 + input$atoms4,
         x.natural=x.natural
         )
    
    return(result)
    
  })
  
  # update input$atoms1 slider
  observe({
    updateSliderInput(session, "atoms1", max = unname(c(inputformula()$h.max, inputformula()$c.max, inputformula()$n.max, inputformula()$o.max)[c('H','C','N','O')==input$element]) - input$atoms2 - input$atoms3 - input$atoms4)
  })
  # update input$atoms2 slider
  observe({
    updateSliderInput(session, "atoms2", max = unname(c(inputformula()$h.max, inputformula()$c.max, inputformula()$n.max, inputformula()$o.max)[c('H','C','N','O')==input$element]) - input$atoms1 - input$atoms3 - input$atoms4)
  })
  # update input$atoms3 slider
  observe({
    updateSliderInput(session, "atoms3", max = unname(c(inputformula()$h.max, inputformula()$c.max, inputformula()$n.max, inputformula()$o.max)[c('H','C','N','O')==input$element]) - input$atoms1 - input$atoms2 - input$atoms4)
  })
  # update input$atoms4 slider
  observe({
    updateSliderInput(session, "atoms4", max = unname(c(inputformula()$h.max, inputformula()$c.max, inputformula()$n.max, inputformula()$o.max)[c('H','C','N','O')==input$element]) - input$atoms1 - input$atoms2 - input$atoms3)
  })
  
  # hot data input table (enrichment calculator)
  output$hot <- renderRHandsontable({
    file1 <- input$upload
    if(is.null(input$hot) & is.null(file1)) { DF = init.df } else { DF = hot_to_r(input$hot) }
    if(!is.null(file1)) {DF_read = data.frame(readxl::read_xlsx(file1$datapath))
        read_start = which(DF_read[,1]=='DATA')
        read_stop = which(DF_read[,1]=='OUTPUT')
        read_formula = gsub('Neutral formula: ', '', DF_read[grep('Neutral formula', DF_read[,1]),1])
        read_element = gsub('Isotopic element: ', '', DF_read[grep('Isotopic element', DF_read[,1]),1])
        read_adduct  = gsub('Adduct: ', '', DF_read[grep('Adduct', DF_read[,1]),1])
        read_charge  = gsub('Charge state: ', '', DF_read[grep('Charge state', DF_read[,1]),1])
        read_tolerance = gsub('Mass tolerance: ', '', DF_read[grep('Mass tolerance', DF_read[,1]),1])
        read_pools   = strsplit(gsub('Isotopic pool sizes: ', '', DF_read[grep('Abundance', DF_read[,1]),1]),',')[[1]]
        read_abund1  = DF_read[grepl('Abundance', DF_read[,1])&grepl('pool1', DF_read[,1]),2]
        read_abund2  = DF_read[grepl('Abundance', DF_read[,1])&grepl('pool2', DF_read[,1]),2]
        read_abund3  = DF_read[grepl('Abundance', DF_read[,1])&grepl('pool3', DF_read[,1]),2]
        read_abund4  = DF_read[grepl('Abundance', DF_read[,1])&grepl('pool4', DF_read[,1]),2]
        
        if(is.numeric(read_start)) {
          # read data frame
          DF_r = data.frame(unname(DF_read[(2+read_start):(read_stop-2),1:2]))
          names(DF_r) = c('mass','intensity')
          DF_r$mass <- as.double(DF_r$mass)
          DF_r$intensity <- as.double(DF_r$intensity)
          DF = DF_r
          
          # update other info from the input file
          updateTextInput(session, 'formula', value = read_formula)
          updateTextInput(session, 'adduct', value = read_adduct)
          updateTextInput(session, 'tolerance', value = read_tolerance)
          updateSelectInput(session, 'element', selected = read_element)
          updateSelectInput(session, 'charge', selected = read_charge)
          if(!is.na(as.numeric(read_pools[1]))) updateSliderInput(session, 'atoms1', label = paste0("Pool 1 [*", input$element,"]"), value = read_pools[1])
          if(!is.na(as.numeric(read_pools[2]))) updateSliderInput(session, 'atoms2', label = paste0("Pool 2 [**", input$element,"]"), value = read_pools[2])
          if(!is.na(as.numeric(read_pools[3]))) updateSliderInput(session, 'atoms3', label = paste0("Pool 3 [***", input$element,"]"), value = read_pools[3])
          if(!is.na(as.numeric(read_pools[4]))) updateSliderInput(session, 'atoms4', label = paste0("Pool 4 [****", input$element,"]"), value = read_pools[4])
          if(!is.na(as.numeric(read_abund1)))   updateSliderInput(session, "abundance1", value = read_abund1)
          if(!is.na(as.numeric(read_abund2)))   updateSliderInput(session, "abundance2", value = read_abund2)
          if(!is.na(as.numeric(read_abund3)))   updateSliderInput(session, "abundance3", value = read_abund3)
          if(!is.na(as.numeric(read_abund4)))   updateSliderInput(session, "abundance4", value = read_abund4)
        }
    }
    rhandsontable(DF, readOnly = FALSE, stretchH = "all", selectCallback = TRUE) %>%
      hot_context_menu(allowColEdit = FALSE ) %>%
      hot_validate_numeric(cols = 1:2, min = 0) %>%
      hot_cols('float', format = '0.000')
  })
  
  user.data  <- reactive({
    if (is.null(input$hot)) { z = init.df } else { z = hot_to_r(input$hot) }
    list(mass=z[,1], intensity=z[,2])
  })
  
  v <- reactiveValues(enr = NULL, rss = NULL, all = NULL, M = NULL, atoms = NULL, atoms.max = NULL, abundances = NULL, R2 = NULL, Nmatches = 0, fitted = NULL)
  
  observeEvent(input$element, {
    
    v$atoms.max <- unname(c(inputformula()$h.max, inputformula()$c.max, inputformula()$n.max, inputformula()$o.max)[c('H','C','N','O')==input$element])
    
    # update input$atoms
    updateSliderInput(session, "atoms1", label = paste0("Pool 1 [*", input$element,"]"), max = v$atoms.max - input$atoms2 - input$atoms3 - input$atoms4)
    updateSliderInput(session, "atoms2", label = paste0("Pool 2 [**", input$element,"]"), max = v$atoms.max - input$atoms1 - input$atoms3 - input$atoms4)
    updateSliderInput(session, "atoms3", label = paste0("Pool 3 [***", input$element,"]"), max = v$atoms.max - input$atoms1 - input$atoms2 - input$atoms4)
    updateSliderInput(session, "atoms4", label = paste0("Pool 4 [****", input$element,"]"), max = v$atoms.max - input$atoms1 - input$atoms2 - input$atoms3)
    
  })
  
  observeEvent(input$button, {
    
      v$atoms.max <- unname(c(inputformula()$h.max, inputformula()$c.max, inputformula()$n.max, inputformula()$o.max)[c('H','C','N','O')==input$element])
      
      # update input$atoms
      updateSliderInput(session, "atoms1", label = input$element, max = v$atoms.max)
      updateSliderInput(session, "atoms2", max = v$atoms.max)
      updateSliderInput(session, "atoms3", max = v$atoms.max)
      updateSliderInput(session, "atoms4", max = v$atoms.max)
      
      qw = optimx(c(max(0.000001, min(input$abundance1,0.99999)), 
                    max(0.000001, min(input$abundance2,0.99999)), 
                    max(0.000001, min(input$abundance3,0.99999)), 
                    max(0.000001, min(input$abundance4,0.99999))), 
                        ipfit0, method = "L-BFGS-B", 
                        lower = c(0, 0, 0, 0) + 1e-6, upper = c(1, 1, 1, 1) - 1e-5, 
                        itnmax = 1000, 
                        element = input$element, 
                        formula = inputformula()$f,
                        charge = inputformula()$ch, 
                        mtolerance = input$tolerance,
                        diso = x.data,
                        user.data = user.data() )
      
      v$enr <- c(qw$p1, qw$p2, qw$p3, qw$p4)
      v$rss <- log(qw$value) # log of residual sum squares
      
      diso = x.data
      if(input$element=='H') {
        diso$abundance[diso$element=='Hx'] <- c(1-v$enr[1], v$enr[1])
        diso$abundance[diso$element=='Hy'] <- c(1-v$enr[2], v$enr[2])
        diso$abundance[diso$element=='Hz'] <- c(1-v$enr[3], v$enr[3])
        diso$abundance[diso$element=='Hw'] <- c(1-v$enr[4], v$enr[4])
      }
      if(input$element=='C') {
        diso$abundance[diso$element=='Cx'] <- c(1-v$enr[1], v$enr[1])
        diso$abundance[diso$element=='Cy'] <- c(1-v$enr[2], v$enr[2])
        diso$abundance[diso$element=='Cz'] <- c(1-v$enr[3], v$enr[3])
        diso$abundance[diso$element=='Cw'] <- c(1-v$enr[4], v$enr[4])
      }
      if(input$element=='N') {
        diso$abundance[diso$element=='Nx'] <- c(1-v$enr[1], v$enr[1])
        diso$abundance[diso$element=='Ny'] <- c(1-v$enr[2], v$enr[2])
        diso$abundance[diso$element=='Nz'] <- c(1-v$enr[3], v$enr[3])
        diso$abundance[diso$element=='Nw'] <- c(1-v$enr[4], v$enr[4])
      }
      if(input$element=='O') {
        diso$abundance[diso$element=='Ox'] <- c(1-v$enr[1], 1e-8, v$enr[1])
        diso$abundance[diso$element=='Oy'] <- c(1-v$enr[2], 1e-8, v$enr[2])
        diso$abundance[diso$element=='Oz'] <- c(1-v$enr[3], 1e-8, v$enr[3])
        diso$abundance[diso$element=='Ow'] <- c(1-v$enr[4], 1e-8, v$enr[4])
      }

    zz = ecipex(inputformula()$f, isoinfo = diso, sortby = "mass", limit=1e-5, id=TRUE)[[1]]
    v$M <- zz[,1] %*% zz[,2]
    
    # update input$abundance sliders
    updateSliderInput(session, "abundance1", value = v$enr[1])
    updateSliderInput(session, "abundance2", value = v$enr[2])
    updateSliderInput(session, "abundance3", value = v$enr[3])
    updateSliderInput(session, "abundance4", value = v$enr[4])
    
    # error : no matching peaks
    z=ecipex(inputformula()$f, isoinfo = diso, sortby = "mass")[[1]]
    z$mass <- (z$mass - inputformula()$ch*0.0005486)/abs(inputformula()$ch)
    df=compact(z, tol=input$tolerance)
    m.j = outer(df$m, user.data()$mass, function(x, y) abs(x-y))
    matches = which(m.j < input$tolerance, arr.ind=TRUE)
    res = data.frame(model.m = df$m[matches[,1]], 
                     model.x = df$x[matches[,1]], 
                     data.m = user.data()$mass[matches[,2]], 
                     data.x = user.data()$intensity[matches[,2]])
    v$Nmatches = length(res$data.x)
    if(v$Nmatches < 2){
      output$text2A  <- renderText({ paste("<font color=\"#FF0000\">", "NO MATCHING PEAKS. Check the input parameters or try increasing the m/z tolerance.", "</font>") })
    }
    
    if(v$Nmatches >= 2){
    
      # pools
      pools = 1
      if(input$atoms2>0 & input$atoms3==0 & input$atoms4==0) pools = 2
      if(input$atoms2>0 & input$atoms3>0 & input$atoms4==0) pools = 3
      if(input$atoms2>0 & input$atoms3>0 & input$atoms4>0) pools = 4
      
      # Proportions of multiple enriched atoms in the molecule ######
      reg = c('^2H','^13C','^15N','^18O')[c(input$element=='H',input$element=='C',input$element=='N',input$element=='O')]
      # Number of enriched atoms and the corresponding abundances
      zz$tot <- apply(as.matrix(zz[,grepl(reg, colnames(zz))], ncol=length(which(grepl(reg, colnames(zz)))  )), 1, sum)
      
      ###############################################################
      v$abundances <- sapply(0:max(zz$tot), function(i) sum(zz$abundance[zz$tot == i]))
      results.enr.table <- data.frame(
        'Number of isotopic atoms' = 0:max(zz$tot),
        'Abundance' = formatC(v$abundances, digits=input$digits, format='f')
        )
      
      results.table = data.frame('element'=c(input$element,c(paste0('[<sup>*</sup>',input$element,']'),
                                           paste0('[<sup>**</sup>',input$element,']'),
                                           paste0('[<sup>***</sup>',input$element,']'),
                                           paste0('[<sup>****</sup>',input$element,']'))[1:pools]),
                                 'isotopic abundance'=formatC(c(inputformula()$x.natural, v$enr[1:pools]), digits=input$digits, format="f")
                                 )
      names(results.table) <- c('Isotopic element',
                                paste0('Abundance of <sup>',c(2,13,15,18)[c(input$element=='H',input$element=='C', input$element=='N', input$element=='O')],'</sup>',input$element)
                                )
      names(results.enr.table) <- c(paste0('Number of <sup>',c(2,13,15,18)[c(input$element=='H',input$element=='C', input$element=='N', input$element=='O')],'</sup>',input$element,' atoms'), 'Abundance')
      
      output$res.table <- renderTable({ results.table }, sanitize.text.function = function(x) x)
      output$enr.table <- renderTable({ results.enr.table }, sanitize.text.function = function(x) x)
      
      output$text2C  <- renderText({ paste0("Molar mass of M: ", formatC(v$M, digits=2, format="f")," g/mol") })
    }
  })
  
  toListen <- reactive({ list(input$abundance1, input$abundance2, input$abundance3, input$abundance4) })
  
  observeEvent(toListen(), {
    
    diso = x.data
    if(input$element=='H') {
      diso$abundance[diso$element=='Hx'] <- c(1-input$abundance1, input$abundance1)
      diso$abundance[diso$element=='Hy'] <- c(1-input$abundance2, input$abundance2)
      diso$abundance[diso$element=='Hz'] <- c(1-input$abundance3, input$abundance3)
      diso$abundance[diso$element=='Hw'] <- c(1-input$abundance4, input$abundance4)
      }
    if(input$element=='C') {
      diso$abundance[diso$element=='Cx'] <- c(1-input$abundance1, input$abundance1)
      diso$abundance[diso$element=='Cy'] <- c(1-input$abundance2, input$abundance2)
      diso$abundance[diso$element=='Cz'] <- c(1-input$abundance3, input$abundance3)
      diso$abundance[diso$element=='Cw'] <- c(1-input$abundance4, input$abundance4)
    }
    if(input$element=='N') {
      diso$abundance[diso$element=='Nx'] <- c(1-input$abundance1, input$abundance1)
      diso$abundance[diso$element=='Ny'] <- c(1-input$abundance2, input$abundance2)
      diso$abundance[diso$element=='Nz'] <- c(1-input$abundance3, input$abundance3)
      diso$abundance[diso$element=='Nw'] <- c(1-input$abundance4, input$abundance4)
    }
    if(input$element=='O') {
      diso$abundance[diso$element=='Ox'] <- c(1-input$abundance1, 1e-8, input$abundance1)
      diso$abundance[diso$element=='Oy'] <- c(1-input$abundance2, 1e-8, input$abundance2)
      diso$abundance[diso$element=='Oz'] <- c(1-input$abundance3, 1e-8, input$abundance3)
      diso$abundance[diso$element=='Ow'] <- c(1-input$abundance4, 1e-8, input$abundance4) 
    }
    
    zz = ecipex(inputformula()$f, isoinfo = diso, sortby = "mass", limit=1e-6)[[1]]
    v$M <- zz[,1] %*% zz[,2]
    v$fitted <- zz
    
    output$pdf <- renderPlot({ 
      
      z = v$fitted
      #### convert masses to mass-to-charge ratios
      z$mass <- (z$mass - as.double(inputformula()$ch)*0.0005486)/abs(as.double(inputformula()$ch))
      #### compact the mass spectrum 
      df=compact(z, tol=input$tolerance)
      #### re-scale the two spectra (model and experimental)
      q=which(user.data()$intensity==max(user.data()$intensity))
      mass.max = user.data()$mass[q][1]
      slope = user.data()$intensity[q][1] / df$x[ abs(df$m - mass.max)==min(abs(df$m - mass.max)) ][1]
      
      xy = user.data()
      max.x = max(xy$mass, df$m)
      min.x = max(0, min(xy$mass, df$m))
      max.y = max(xy$intensity, df$x*slope)
      par(mar=c(4.5, 4.1, 2, 3.1))
      plot(xy$mass, xy$intensity, type="h", main="", xlim=c(min.x, max.x), ylim=c(0, max.y), 
           xlab="Mass-to-charge ratio", ylab="Intensity", cex=1.8, cex.axis=1.3, cex.lab=1.3)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "gray95")
      points(x=xy$mass, y=xy$intensity, pch=21, cex=3, col='tomato', bg=NULL, lwd=2)
      # segments(xy$mass-2*input$tolerance, xy$intensity, xy$mass+2*input$tolerance, xy$intensity, col='tomato', lwd=2)
      segments(df$m, 0, df$m, df$x*slope, lwd=4)
      
      v$fit = list(mass = df$m, intensity = df$x*slope)
      
      #### join the experimental and theoretical spectra
      m.j = outer(df$m, xy$mass, function(x, y) abs(x-y))
      matches = which(m.j < input$tolerance, arr.ind=TRUE)
      v$rss <- log(sum(slope*df$x[matches[,1]] - xy$intensity[matches[,2]])^2)
      if(v$Nmatches>=2) v$R2  <- summary(lm(I(df$x[matches[,1]]) ~ 0 + I(xy$intensity[matches[,2]])))$r.squared
      mtext(side=3,line=0.2,adj=1,text=paste0('ln(RSS) = ',formatC(v$rss, digits=3, format='e'),', R2 = ',formatC(v$R2, digits=3,format='f')), cex=0.9, col='tomato',font=2)
      box()
    })
    
  })

  download.out <- reactive({
    
    z=v$fitted
    #### convert masses to mass-to-charge ratios
    z$mass <- (z$mass - as.double(inputformula()$ch)*0.0005486)/abs(as.double(inputformula()$ch))
    #### compact the mass spectrum 
    df=compact(z, tol=input$tolerance)
    #### rescale the two spectra (model and experimental)
    q=which(user.data()$intensity==max(user.data()$intensity))
    mass.max = user.data()$mass[q][1]
    slope = user.data()$intensity[q][1] / df$x[ abs(df$m - mass.max)==min(abs(df$m - mass.max)) ][1]
    
    NN = length(user.data()$mass)
    Nen = length(v$abundances)
    Nfit = length(df$m)
    
    iso = paste0(c(2,13,15,18)[c(input$element=='H',input$element=='C', input$element=='N', input$element=='O')],input$element)
    res = matrix('', nrow = 29+NN+Nen+Nfit, ncol=2)
    res[1,1] = 'DESCRIPTION'
    res[2,1] = 'id:'
    res[3,1] = 'notes:'
    res[4,1] = paste('Neutral formula:', input$formula)
    res[5,1] = paste('Isotopic element:', input$element)
    res[6,1] = paste('Adduct:', input$adduct)
    res[7,1] = paste('Charge state:', input$charge)
    res[8,1] = paste('Isotopic pool sizes:', paste(c(input$atoms1,input$atoms2,input$atoms3,input$atoms4), collapse=','))
    res[9,1] = paste('Mass tolerance:', input$tolerance)
    res[11,1] = 'DATA'
    res[12,1] = 'm/z'
    res[12,2] = 'Intensity'
    for(i in 1:NN) res[12+i, 1] = user.data()$mass[i]
    for(i in 1:NN) res[12+i, 2] = user.data()$intensity[i]
    res[14+NN,1] = 'OUTPUT'
    res[16+NN,1] = 'Goodness-of-fit, ln(RSS):'
    res[16+NN,2] = formatC(v$rss, digits=3, format='e')
    res[17+NN,1] = 'Goodness-of-fit, R^2:'
    res[17+NN,2] = formatC(v$R2, digits=3, format='e')
    
    res[18+NN,1] = gsub('<sub>','',gsub('</sub>','',gsub('</sup>','',gsub('<sup>','',paste0("Measured ion, [M",ifelse(inputformula()$ch>0,'+','-'), inputformula()$foo.a,"]","<sup>",abs(inputformula()$ch),ifelse(inputformula()$ch>0,'+','-'),'</sup> : ',inputformula()$foo.t)))))
    res[19+NN,1] = paste0('Abundance of ',iso,' in ',input$element,' (natural):')
    res[19+NN,2] = inputformula()$x.natural
    res[20+NN,1] = paste0('Abundance of ',iso,' in [*',input$element,'] (pool1):')
    res[20+NN,2] = formatC(v$enr[1],digits=4,format='f')
    res[21+NN,1] = paste0('Abundance of ',iso,' in [**',input$element,'] (pool2):')
    res[21+NN,2] = formatC(v$enr[2],digits=4,format='f')
    res[22+NN,1] = paste0('Abundance of ',iso,' in [***',input$element,'] (pool3):')
    res[22+NN,2] = formatC(v$enr[3],digits=4,format='f')
    res[23+NN,1] = paste0('Abundance of ',iso,' in [****',input$element,'] (pool4):')
    res[23+NN,2] = formatC(v$enr[4],digits=4,format='f')
    res[25+NN,1] = paste0('Number of ',iso,' atoms')
    res[25+NN,2] = paste0('Abundance')
    res[seq(26+NN,length.out = Nen),] = cbind(seq(0,length.out=length(v$abundances)), v$abundances)
    
    res[27+NN+Nen,1] = paste0('FITTED VALUES')
    res[28+NN+Nen,1] = paste0('mass')
    res[28+NN+Nen,2] = paste0('intensity')
    res[seq(29+NN+Nen,length.out = Nfit),] = cbind(df$m, df$x*slope)
    
    return( as_tibble(res, .name_repair = 'minimal') )
  })
  
  output$text1A  <- renderText({ paste0("Molecular ion, M : ", inputformula()$foo.m) })
  output$text1B  <- renderText({ paste0("Measured ion, [M",ifelse(inputformula()$ch>0,'+','-'), inputformula()$foo.a,"]","<sup>",abs(inputformula()$ch),ifelse(inputformula()$ch>0,'+','-'),'</sup> : ',inputformula()$foo.t) })
  
  # Download the configuration file
  
  output$download <- downloadHandler(
    filename = function() { "out.xlsx" },
    content = function(file) { write_xlsx(download.out(), col_names = FALSE, path = file) }
  )
  
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
        column(10, h5(tags$b("Number of isotopic atoms"))),
        width = '85%'
      ),
      fluidRow(
      column(5, sliderInput("atoms1", label = "Pool 1", min=1, max=13, value=13, round=TRUE, step = 1, ticks = TRUE)),
      column(5, sliderInput("atoms2", label = "Pool 2", min=0, max=13, value=0, round=TRUE, step = 1, ticks = TRUE)),
      width='85%'
            ),
      fluidRow(
        column(5, sliderInput("atoms3", label = "Pool 3", min=0, max=13, value=0, round=TRUE, step = 1, ticks = TRUE)),
        column(5, sliderInput("atoms4", label = "Pool 4", min=0, max=13, value=0, round=TRUE, step = 1, ticks = TRUE)),
        width='85%'
      ),
      fluidRow(
        column(5, textInput("adduct", label = "Specify adduct", value="H2") ),
        column(5, selectInput("charge", label = "Charge-state", choices=c('-2'=-2,'-1'=-1,'+1'=1,'+2'=2), selected=2) ),
        width = '85%' ),
      h5("Additional settings", a(id = "toggleextra", "show/hide")),
      shinyjs::hidden(div(id = "filterextra",
                          fluidRow(column(12, numericInput("digits", label = "Number of decimal digits", value = 4, min = 2, max = 8, step = 1, width='85%')) ),
                          fluidRow(
                            column(12, sliderInput("abundance1", label = "Explore data by-hand", min=0, max=1, value=0.98, round=TRUE, step=0.01, width='85%') ),
                            column(11, helpText(textOutput("A1"))),
                            column(12, sliderInput("abundance2", label = "", min=0, max=1, value=0.98, round=TRUE, step=0.01, width='85%') ),
                            column(11, helpText(textOutput("A2"))),
                            column(12, sliderInput("abundance3", label = "", min=0, max=1, value=0.98, round=TRUE, step=0.01, width='85%') ),
                            column(11, helpText(textOutput("A3"))),
                            column(12, sliderInput("abundance4", label = "", min=0, max=1, value=0.98, round=TRUE, step=0.01, width='85%') ),
                            column(11, helpText(textOutput("A4"))),
                            column(12, numericInput("tolerance", label = "Enter m/z tolerance", value = 0.1, max = 0.5, min=0.0001, step = 0.001, width = '85%'))
                            ))),
      h5(tags$b("Enter (paste) the observed isotopic pattern")),
      rHandsontableOutput("hot"),
      helpText("right-click to add or delete rows"),
      br(),
      #fluidRow(
      #  column(11, fileInput("upload", "Upload configuration file", multiple = FALSE, accept = ".xlsx"),
      #         ), width = '85%' ),
      br(),
      conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",   
                       actionButton("button", label = "Perform least squares fitting", icon = icon('bar-chart-o'))),
      conditionalPanel(condition = "$('html').hasClass('shiny-busy')",   
                       actionButton("button", label = "busy...", icon = icon('hourglass-half'))),
      br()
    ),
    
    mainPanel(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"),
      fluidRow(column(helpText("This calculator evaluates the experimental mass spectra for a given substance and determines the isotopic composition of a selected element that best matches the experimental data."), width=11)),
      br(),
      fluidRow(column(6, tableOutput("res.table")),
          column(6, fluidRow(
          strong(htmlOutput("text1A")),
          strong(htmlOutput("text1B")),
          strong(htmlOutput("text2A")),
          strong(htmlOutput("text2B")),
          strong(htmlOutput("text2C"))
        ))  
      ),
      fluidRow(column(plotOutput(("pdf")),width=11)),
      conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",   
                       br(), strong(htmlOutput("text1b"))
      ),
      fluidRow(column(6, tableOutput("enr.table"))),
      br(),
      br(),
      h5("Notes and explanations", a(id = "togglenotes", "show/hide")),
      shinyjs::hidden(div(id = "filternotes",
                          fluidRow(column(
                            p("Red circles: experimental data; Vertical black lines: theoretical model; ln(RSS): the natural logarithm of the residual sum squares, a measure of how well the model fits the data."),
                            p("The calculator performs simple fitting of the isotope patterns and, at this stage, does not consider the uncerrtainty of the natural isotopic composition of the elements nor any instrumental isotope ratio fractionation effects during the mass spectrometry measurements."),
                            p("Created using R and Shiny. Juris Meija (2019-2021) NRC Canada"), width=11)
                          ))),br(),
      downloadButton("download", "Download configuration file (xlsx)"),br(),br(),
      fileInput("upload", label=NULL, placeholder="Upload configuration file (xlsx)", multiple = FALSE, accept = ".xlsx"),
      br(),br(),
      p("NRC Isotopic Enrichment Calculator (2019-2021) December 2021")
      
    )
  )
)

shinyApp(ui = ui, server = server)