quick_plot <- function(ODEinput){
  
  
  par(mfrow = c(2, 2))
  plot(ODEinput$CS ~ ODEinput$times, type = 'l', ylim = c(0, 55), 
       lwd = 3, col = 'blue', main = "Cattle (no Prophylaxis)", 
       xlab = "Time", ylab = "Number")
  lines(ODEinput$CEs ~ ODEinput$times, lwd = 3, col = 'orange') # Exposed
  lines(ODEinput$CEr ~ ODEinput$times, lwd = 3, col = 'orange', lty = 2) # Exposed
  lines(ODEinput$CIs ~ ODEinput$times, lwd = 3, col = 'red') # Infected
  lines(ODEinput$CIr ~ ODEinput$times, lwd = 3, col = 'red', lty = 2) # Infected
  lines(ODEinput$CTs ~ ODEinput$times, lwd = 3, col = 'green') # Treated
  lines(ODEinput$CTr ~ ODEinput$times, lwd = 3, col = 'green', lty = 2) # Treated
  lines(( ODEinput$CEs + ODEinput$CEr + ODEinput$CIs + ODEinput$CIr + ODEinput$CTs + ODEinput$CTr + ODEinput$CS) ~ODEinput$times,lty = 2)
  lines(( ODEinput$CEs + ODEinput$CEr + ODEinput$CIs + ODEinput$CIr + ODEinput$CTs + ODEinput$CTr + ODEinput$CS +
            ODEinput$PEs + ODEinput$PEr + ODEinput$PIs + ODEinput$PIr + ODEinput$PTs + ODEinput$PTr + ODEinput$PS + ODEinput$PF+ ODEinput$PPr + ODEinput$PPs) ~ODEinput$times,lty = 2, col = "pink")
  #legend("topright", col = c("blue", "orange", "darkorange", "red","darkred", "green", "darkgreen", "grey"),
  #       y = c("CS", "CEs", "CEr", "CIs", "CIr","CTs", "CTr"), pch = 15)
  
  plot( ODEinput$PF ~ ODEinput$times, type = 'l', ylim = c(0, 55),col = 'purple',
        lwd = 3,main = "Cattle (with Prophylaxis)",xlab = "Time", ylab = "Number")
  lines(ODEinput$PS ~ ODEinput$times, lwd = 3, col = 'blue') # Exposed
  lines(ODEinput$PEs ~ ODEinput$times, lwd = 3, col = 'orange') # Exposed
  lines(ODEinput$PEr ~ ODEinput$times, lwd = 3, col = 'orange', lty = 2) # Exposed
  lines(ODEinput$PIs ~ ODEinput$times, lwd = 3, col = 'red') # Infected
  lines(ODEinput$PIr ~ ODEinput$times, lwd = 3, col = 'red', lty = 2) # Infected
  lines(ODEinput$PTs ~ ODEinput$times, lwd = 3, col = 'green') # Treated
  lines(ODEinput$PTr ~ ODEinput$times, lwd = 3, col = 'green', lty = 2) # Treated
  lines(ODEinput$PPs ~ ODEinput$times, lwd = 3, col = 'grey') # Treated
  lines(ODEinput$PPr ~ ODEinput$times, lwd = 3, col = 'grey', lty = 2) # Treated
  lines((ODEinput$PEs + ODEinput$PEr + ODEinput$PIs + ODEinput$PIr + ODEinput$PTs + ODEinput$PTr  + ODEinput$PS + ODEinput$PF + ODEinput$PPr + ODEinput$PPs) ~ODEinput$times, lty = 2)
  lines(( ODEinput$CEs + ODEinput$CEr + ODEinput$CIs + ODEinput$CIr + ODEinput$CTs + ODEinput$CTr + ODEinput$CS +
            ODEinput$PEs + ODEinput$PEr + ODEinput$PIs + ODEinput$PIr + ODEinput$PTs + ODEinput$PTr + ODEinput$PS+ ODEinput$PF+ ODEinput$PPr + ODEinput$PPs) ~ODEinput$times,lty = 2, col = "pink")
  #legend("topright", col = c("purple","blue", "orange", "darkorange", "red","darkred", "green", "darkgreen", "grey", "black", "grey"),
  #       y = c("PF","PS", "PEs", "PEr", "PIs", "PIr","PTs", "PTr", "PPs", "PPr"), pch = 15)
  
  plot(ODEinput$WS ~ ODEinput$times,type = 'l', ylim = c(0, max(ODEinput$WS+10)),col = 'blue',lwd = 3,
       main = "Wildlife", xlab = "Time", ylab = "Number")
  lines(ODEinput$WEs ~ ODEinput$times, lwd = 3, col = 'orange') # Exposed
  lines(ODEinput$WEr ~ ODEinput$times, lwd = 3, col = 'orange', lty = 2) # Exposed
  lines(ODEinput$WIs ~ ODEinput$times, lwd = 3, col = 'red') # Infected
  lines(ODEinput$WIr ~ ODEinput$times, lwd = 3, col = 'red', lty = 2) # Infected
  lines((ODEinput$WEs + ODEinput$WEr + ODEinput$WIs + ODEinput$WIr + ODEinput$WS) ~ ODEinput$times,lty = 2)
  #legend("topright", col = c("blue", "orange", "darkorange", "red","darkred"),
  #       y = c("WS", "WEs", "WEr", "WIs", "WIr"), pch = 15)
  
  plot(ODEinput$VSt ~ ODEinput$times, type = 'l', ylim = c(0, max(ODEinput$VSt)+100), col = 'blue',
       lwd = 3, main = "Vector", xlab = "Time", ylab = "Number")
  lines(ODEinput$VSf ~ ODEinput$times, lwd = 3, col = 'lightblue') # Exposed
  lines(ODEinput$VEs ~ ODEinput$times, lwd = 3, col = 'orange') # Exposed
  lines(ODEinput$VEr ~ ODEinput$times, lwd = 3, col = 'orange', lty = 2) # Exposed
  lines(ODEinput$VIs ~ ODEinput$times, lwd = 3, col = 'red') # Infected
  lines(ODEinput$VIr ~ ODEinput$times, lwd = 3, col = 'red', lty = 2) # Infected
  lines((ODEinput$VEs + ODEinput$VEr + ODEinput$VIs + ODEinput$VIr + ODEinput$VSt + ODEinput$VSf) ~
          ODEinput$times, lty = 2)
  #legend("topright", col = c("blue", "lightblue","orange", "darkorange", "red","darkred"),
  #  €   y = c("VSt","VSf", "VEs", "VEr", "VIs", "VIr"), pch = 15)
  
}