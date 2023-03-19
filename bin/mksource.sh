#!/bin/ksh
ls *.F90 2>/dev/null |\
awk  '{
    CC[NR]=$1
}
END { 
    if (NR>0) {
       printf("F90SOURCE = \\\n")
       for (i = 1; i < NR; i++) printf("\t%s\\\n", CC[i])
       printf("\t%s\n\n", CC[NR])
    }
}'

ls *.F 2>/dev/null |\
awk  '{
    CC[NR]=$1
}
END { 
    if (NR>0) {
       printf("F77SOURCE = \\\n")
       for (i = 1; i < NR; i++) printf("\t%s\\\n", CC[i])
       printf("\t%s\n\n", CC[NR])
    }
}'

ls *.H  2>/dev/null |\
awk  '{
    CC[NR]=$1
}
END { 
    if (NR > 0) {
       printf("INCSOURCE = \\\n")
       for (i = 1; i < NR; i++) printf("\t%s\\\n", CC[i])
       printf("\t%s\n\n", CC[NR])
    }
}'

