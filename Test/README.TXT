OLD SCHEME

-m module (clean, trinity/spade, map, align, subset)
-s [modules=trinity,map,align] name of txt file containg list of PE samples to process
-i [modules=trinity,map,align] name of txt file containg list of SE samples to process
-e [modules=ALL] name of reference fasta file (aa)
-r [module=clean] full path to folder containing the raw sequenes
-p [modules=ALL] number of threads to use for multithreading steps
-x [module=subset] name of txt file containing list of subset exons  
-t [module=subset] name of txt file containing list of subset samples
-n [module=align,subset] suffix for alignment and subset runs

NEW SCHEME

-m module (clean, trinity/spade, map, align, subset)
-p [modules=trinity,map,align] name of txt file containg list of PE samples to process 
-s [modules=trinity,map,align] name of txt file containg list of SE samples to process
-a [modules=ALL] name of reference fasta file (aa)
-r [module=clean] full path to folder containing the raw sequenes
-t [modules=ALL] number of threads to use for multithreading steps
-e [module=subset] name of txt file containing list of subset exons  
-f [module=subset] name of txt file containing list of subset samples
-x [module=align,subset] suffix for alignment and subset runs  

-m -m
-s -p
-i -s
-e -a
-r -r
-p -t
-x -e
-t -f
-n -x
