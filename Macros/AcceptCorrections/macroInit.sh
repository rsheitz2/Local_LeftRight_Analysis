#!/bin/bash

if [ $# -ne 1 ]; then
     echo "" 
     echo "Script gets macros ready for git commits"
     echo ""
     echo "Initalization is:"
     echo "    toWrite =false"
     echo ""
     echo "Enter any character to initialized macros (i.e. macroInit.sh 1)"
     echo ""
     
else

    echo "sed -i.bak 's/toWrite =true/toWrite =false/' *.C"
    echo "mv *.bak Tests"
    echo "sed -i.bak 's/toWrite =true/toWrite =false/' FalseAsym/*.C"
    echo "mv FalseAsym/*.bak Tests"
    
    sed -i.bak 's/toWrite =true/toWrite =false/' *.C
    mv *.bak Tests
    sed -i.bak 's/toWrite =true/toWrite =false/' FalseAsym/*.C
    mv FalseAsym/*.bak Tests
fi
