#
# Run the Shiny app: Single Stock Calculator 
#
SingleStockCalculator <- function() {
  appDir <- system.file("SingleStockCalculator", package="fishsizespectrum")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}


SingleStockAssessment <- function(W=10000, param=baseparameters(W), F=0.3) {
  # Calculate reference spectrum:
  N0 = spectrum(param)
  
  # Set fishing mortality:
  p$F = F
  
  # Calc fished spectrum:
  N = spectrum(param)
  
  # Calculate reference points:
  ref = calcRefpoints(param)
  
  # Print them out:
  cat("Fisheries impact assessment on a stock with asymptotic weight", p$W, "g\n",
      " Fishing mortality:",p$F,"yr^-1\n",
      " \nState:\n",
      " Spawning stock biomass (SSB) relative to max SSB:", N$SSBperR[1]/N0$SSBperR[1], "\n",
      " SSB/Bmsy:", N$SSBperR[1]/ref$Bmsy,"\n",
      " SSB/Blim", N$SSBperR[1]/ref$Blim,"\n",
      " Recruitment relative to max recruitment:", N$R[1], "\n",
      " \nReference points:\n",
      " Fmsy: ", ref$Fmsy, "yr^-1\n",
      " Fmax: ", ref$Fmax, "yr^-1\n",
      " Flim: ", ref$Flim, "yr^-1\n",
      " Fcrash:", ref$Fcrash, "yr^-1\n"
  )
  
  # Add fisheries induced evolution results:
  
  # Add the QG parameters to the param structure
  paramQG = baseparamQG(p=param)
  
  S = calcSelectionResponse(paramQG, W, F)
  cat("\n",
      "Selection responses:\n",
      " Weight at maturation:", 100*S$dwmdt, "%/yr^-1\n",
      " Growth rate:",100*S$dAdt, "%/yr^-1\n",
      " Investment in reproduction:", 100*S$dkrdt, "%/yr^-1\n"
  )
}