#!/bin/bash
### timer ---------------------------------------------------------------------------------
# Function to start the timer 
startimer () {
  start=$(date +%s) # Take the system time at the start 
}
# Usage: 
# > startimer

# Function to end the timer AFTER startimer function started within a terminal session
endtimer () {
  end=$(date +%s) # Take the system time at the end 
  
  seconds=$(echo "$end - $start" | bc) # bc is a "basic calculator"
  # echo $seconds' sec' # Shows difference in seconds 
  
  echo 'Formatted:' # Prints a title for the formated time 
  awk -v t=$seconds 'BEGIN{t=int(t*1000); printf "%d:%02d:%02d\n", t/3600000, t/60000%60, t/1000%60}'
  # This last command will show the time in 00:00:00
}