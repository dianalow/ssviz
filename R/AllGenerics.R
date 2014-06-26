#############
# Accessors #
#############


#########
# ssviz #
#########

setGeneric("readBam",
           function(file_name,tags=character(0)){
            standardGeneric("readBam")
})

setMethod("readBam",signature=c(file_name="character"),
          definition=function(file_name,tags=character(0)){
            #reads bam files together with any additional tags and outputs the dataframe formatted version
            #returns a data.frame version as described in the Rsamtools manual
            bam<-scanBam(file_name,param=ScanBamParam(what=scanBamWhat(),tag=tags))
            tagnames<-names(bam[[1]]$tags)
            bam <- unname(bam) # names not useful in unlisted result
            bam_field<-names(bam[[1]])
            lst <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
            bam.df<-do.call("DataFrame", lst)
            if(is.null(tagnames)) {
              names(bam.df)<-bam_field
            } else {names(bam.df)<-c(bam_field[-length(bam_field)],tagnames)}
            return(bam.df)
          })

setGeneric("getCountMatrix",
           function(bam_file,pseudo=FALSE){
             standardGeneric("getCountMatrix")
           })

setMethod("getCountMatrix",signature=c(bam_file="DataFrame"),
          definition=function(bam_file,pseudo=FALSE){
            # returns the bam data.frame with an additional column counts
            # only relevant if the fasta file used for mapping input was previously collapsed via fastx_toolkit
            # to return a fasta read name in the format of : readnumber-totalcounts
            # pseudo option adds "1" as count
            if(pseudo){
              counts<-1
            } else {
              bsplit<-data.frame(matrix(unlist(strsplit(bam_file$qname,"-")),ncol=2,byrow=TRUE))
              colnames(bsplit)<-c("readname","counts")
              counts<-as.numeric(as.character(bsplit$counts))
            }
            #bam1<-cbind(bam_file,DataFrame(counts))
            bam1<-bam_file
            bam1$counts<-counts
            return(bam1)
          })

setGeneric("ntfreq",
           function(bam_file,ntlength,toRNA=TRUE,count_type="total"){
             standardGeneric("ntfreq")
           })

setMethod("ntfreq",signature=c(bam_file="DataFrame",ntlength="numeric"),
          definition=function(bam_file,ntlength,toRNA=TRUE,count_type="total"){
            # calculates nucleotide frequency of reads in bam file
            
            fun<-function(data,counts){
              bases<-c("A","C","G","T")
              res<-as.matrix(aggregate(counts,by=list(data),sum))
              mm<-match(bases,res[,1],nomatch=0)
              wh<-which(mm==0)
              
              if(length(wh)>0){
                res<-matrix(rbind(res,matrix(cbind(bases[wh],rep(0,length(wh))),ncol=2)),ncol=2)
              }
              
              res[,2]<-as.numeric(res[,2])/sum(counts)
              return(res)
            }
            
            if(!("counts" %in% colnames(bam_file))) {
              message("No counts column provided in bam file. Running unique counts.\nOtherwise provide a count matrix or run getCountMatrix.")
              count_type<-"unique" #assume no counts provided
              bam_file$counts<-1
            }
            
            #get counts and reverse complement the opposite strand
            b1<-unique(bam_file[c("qname","strand","seq","counts")]) # to remove multimapped
            #b1<-getCountMatrix(b0)
            
            if(count_type=="unique") {
              b1$counts<-1 #easy-out for same workflow
            }
            
            b1.ss<-cbind(b1[b1$strand=='+',]$counts,sapply(as.list(b1[b1$strand=='+',]$seq),toString))
            b1.as<-cbind(b1[b1$strand=='-',]$counts,sapply(as.list(reverseComplement(b1[b1$strand=='-',]$seq)),toString))
            if(nrow(b1.ss)==0) x<-data.frame(b1.as)
            else if(nrow(b1.as)==0) x<-data.frame(b1.ss)
            else x<-data.frame(rbind(b1.ss,b1.as))
            #x<-data.frame(cbind(b1$counts,sapply(as.list(b0$seq),toString)))
            colnames(x)<-c("counts","seq")
            
            x$seq<-substr(x$seq,1,ntlength)
            x1<-t(matrix(unlist(apply(matrix(x$seq,ncol=1),1,strsplit,"")),nrow=ntlength))
            result<-t(apply(x1,2,fun,as.numeric(x$counts)))
            
            ntorder<-result[,c(1:4)]
            freq<-result[,c(5:8)]
            for(i in c(1:ntlength)){
              freq[i,]<-freq[i,][order(ntorder[i,])]
            }
            freq<-matrix(as.numeric(freq),ncol=4)
            colnames(freq)<-ntorder[1,order(ntorder[1,])]
            if(toRNA) {
              colnames(freq)<-c("A","C","G","U")
            }
            return(as.data.frame(freq))    
          })

setGeneric("pingpong",
          function(bam_file){
            standardGeneric("pingpong")
          })

setMethod("pingpong",signature=c(bam_file="DataFrame"),
          definition=function(bam_file){
            # piRNA ping-pong analysis of complementary sequences
            # 1. selects for first base "U" (or "T" in this case)
            # 2. finds overlapping reads on different strands
            # 3. Returns overlap table
            # 4. Able to discriminate different rnames
            if(!("counts" %in% colnames(bam_file))) {
              message("No counts column provided in bam file. Running unique counts.\nOtherwise provide a count matrix or run getCountMatrix.")
              bam_file$counts<-1
            }
            
            position=c()  
            pos_reads<-bam_file[bam_file$strand=="+",]
            neg_reads<-bam_file[bam_file$strand=="-",]
            maxposlength<-max(pos_reads$qwidth)
            listrnames<-suppressWarnings(as.character(unique(bam_file$rname)))
            message("Checking rnames: ")
            message(listrnames)
            
            for(k in c(1:length(listrnames))){
              message(paste('\n\n',listrnames[k],'\n***************\n',sep=""))
              pr<-pos_reads[pos_reads$rname==listrnames[k],]
              nr<-neg_reads[neg_reads$rname==listrnames[k],]
              message(paste("Unique forward reads:",nrow(pr),', reverse reads:',nrow(nr),'\n'))
              if(nrow(pr))
                toiterate<-round_any(nrow(pr)/10,10,f=floor)
              if(toiterate==0) toiterate<-1
              for(i in 1:nrow(pr)){
                if(substr(sapply(as.list(pr[i,]$seq),toString),1,1)!="T") next
                if(i%%toiterate==1) message(paste("Forward read #",i,'out of',nrow(pr),"\n"))
                locationp<-pr[i,]$pos
                nr1<-nr[nr$pos<=(locationp) & nr$pos>=(locationp-maxposlength),]
                if(nrow(nr1)!=0){
                  for(j in 1:nrow(nr1)){
                    locationn<-nr1[j,]$pos+nr1[j,]$qwidth
                    difference<-locationn-locationp
                    position<-c(position,rep(difference,nr1[j,]$counts))
                  }
                }
              }    
            }
            
            ftable<-as.data.frame(table(position),stringAsFactors=FALSE)
            ftable$Freq<-ftable$Freq/sum(ftable$Freq)
            message("Done.")
            return(ftable)
          })

setGeneric("plotDistro",
           function(bamlist,type="qwidth",samplenames=NULL,unique=FALSE,ncounts=NULL,norm=FALSE,yname="Counts per million"){
             standardGeneric("plotDistro")
           })

setMethod("plotDistro",signature=c(bamlist="list"),
          definition=function(bamlist,type="qwidth",samplenames=NULL,unique=FALSE,
                              ncounts=NULL,norm=FALSE,yname="Counts per million"){
            value<-Condition<-NULL
            if(!(is.null(samplenames))){
              if(length(samplenames)!=length(bamlist)) {
                stop("Please provide number of samplesnames equivalent to number of bam files!")
              }
            } else {
              samplenames=seq(1,length(bamlist))
            }
            if(!(type %in% c("qwidth","rname","strand"))){
              stop("Invalid type. Should be one of \"qwidth\",\"rname\",\"strand\"")
            }
            
            if(!(is.null(ncounts))){
              if(length(ncounts)!=length(bamlist)) {
                stop("Please provide number of counts equivalent to number of bam files!")
              }
            }
            
            for(i in 1:length(bamlist)){
              bam_file<-bamlist[[i]]
              if(!("counts" %in% colnames(bam_file))) {
                message("No counts column provided in bam file. Running unique counts.\nOtherwise provide a count matrix or run getCountMatrix.")
                bam_file$counts<-1
              }
              b1<-unique(bam_file[c("counts",type,"qname")])
              if(unique) b1$counts<-1
              b11<-aggregate(b1$counts,by=list(b1[[type]]),sum)
              if(i>1) {
                counts<-merge(counts,b11,by=1,all=TRUE)    
              } else {
                counts<-b11
              }
              counts[is.na(counts)]<-0
            }
            rownames(counts)<-counts[,1]
            counts<-counts[-1]
            if(is.null(ncounts)){
              ncounts=colSums(counts)
            }
            counts<-t(t(counts)*1000000/ncounts)
            
            if(norm) {
              yname="Relative count level"
              counts<-counts/counts[,1]
            }
            colnames(counts)<-samplenames
            cplot<-melt(counts,id.vars=colnames(counts))
            colnames(cplot)<-c("length","Condition","value")
            cplot$Condition<-ordered(cplot$Condition,levels=samplenames) #proper sample ordering!
            
            switch(type,
                   qwidth={xlabname="Tag length"},
                   strand={xlabname="Strand"},
                   rname={xlabname="Location"}
            )
            
            theplot<-ggplot(cplot,aes(x=factor(length),y=value))+geom_bar(aes(fill=Condition),colour="black",stat="identity",position="dodge",width=0.7)+
              scale_fill_brewer(palette="Set1")+theme_bw(base_size=18)+scale_x_discrete(name=xlabname)+scale_y_continuous(name=yname)
            return(theplot)
          })

setGeneric("plotFreq",
           function(freqvector,percentage=TRUE){
             standardGeneric("plotFreq")
           })

setMethod("plotFreq",signature=c(freqvector="data.frame",percentage="logicalORmissing"),
          definition=function(freqvector,percentage=TRUE){
            value<-variable<-NULL
            y<-freqvector
            if(percentage){
              y<-y*100
              ylabname<-"Percentage"
            } else {ylabname="Frequency"}
            y$ID <- rownames(y)
            y.melt <- melt(y, id.var = 'ID')
            y.melt <- within(y.melt, ID <- factor(ID,rownames(y),ordered = TRUE))
            ggplot(y.melt, aes(x = ID, y = value, fill = variable)) +scale_fill_manual(values=c("#CC0000","#0000CC","#009900","#FFCC00"))+
              geom_bar(stat = 'identity',colour="black") + xlab("Nucleotide") +ylab(ylabname) + theme_bw(base_size=18)+theme(legend.title=element_blank())
            
          })

setGeneric("plotPP",
           function(pout,samplenames=NULL){
             standardGeneric("plotPP")
          })

setMethod("plotPP",signature=c(pout="list"),
          definition=function(pout,samplenames=NULL){
            position<-Condition<-value<-NULL
            if(!(is.null(samplenames))){
              if(length(samplenames)!=length(pout)) {
                stop("Please provide number of samplesnames equivalent to number of bam files!")
              }
            } else {
              samplenames=seq(1,length(pout))
            }
            
            pp_total=c()
            for(i in 1:length(pout)){
              pp1<-pout[[i]]
              tt<-as.numeric(levels(pp1$position))
              positions<-seq(min(tt),max(tt),1)
              pp1<-merge(cbind(positions,0),pp1,by=1,all=TRUE)[,c(1,3)]
              pp1[is.na(pp1)]<-0
              if(i>1) pp_total<-merge(pp_total,pp1,by=1,all=TRUE)    
              else pp_total<-pp1
              pp_total[is.na(pp_total)]<-0
            }
            pp_total<-pp_total[order(pp_total[,1]),]
            rownames(pp_total)=pp_total[,1]
            colnames(pp_total)<-c("position",samplenames)
            pp_total.m<-melt(pp_total,id="position")
            colnames(pp_total.m)<-c("position","Condition","value")
            ggplot(data=pp_total.m,aes(x=position,y=value,colour=Condition))+geom_line(aes(group=Condition),size=1.5)+
              xlab("Distance from 5' ends (nt)")+
              scale_color_brewer(palette="Set1")+theme_bw(base_size=18)+
              scale_y_continuous(name="Frequency")
          })

setGeneric("plotRegion",
           function(bamlist,region,howsmooth=2,ncounts=NULL,samplenames=NULL){
             standardGeneric("plotRegion")
           })

setMethod("plotRegion",signature=c(bamlist="list",region="character"),
          definition=function(bamlist,region="",howsmooth=2,ncounts=NULL,samplenames=NULL){
            position<-Condition<-value<-NULL
            #check bam contents, will not plot if not mapped to chromosome, or params wrong
            if(!(is.null(samplenames))){
              if(length(samplenames)!=length(bamlist)) {
                stop("Please provide number of samplesnames equivalent to number of bam files!")
              }
            } else {
              samplenames=seq(1,length(bamlist))
            }
            if(!(is.null(ncounts))){
              nFlag=FALSE
              if(length(ncounts)!=length(bamlist)) {
                stop("Please provide number of counts equivalent to number of bam files!")
              }
            } else {
              nFlag=TRUE
              ncounts=rep(1,length(bamlist))
            }
            #setting global vars
            maxdensity=0
            multiplier=0
            yseries=list()
            xseries=list()
            #begin looping for set of bam files
            for(ib in 1:length(bamlist)){
              bam_file<-bamlist[[ib]]
              
              if(substr(bam_file$rname[1],1,3)!="chr") {
                stop("Bam file not mapped to chromosome! Abort!")
              }
              
              if(!("counts" %in% colnames(bam_file))) {
                message("No counts column provided in bam file. Running unique counts.\nOtherwise provide a count matrix or run getCountMatrix.")
                bam_file$counts<-1
              }
              regionsplit<-unlist(strsplit(region,'[:-]'))
              theregion<-bam_file[bam_file$rname==regionsplit[1] & bam_file$pos>as.numeric(regionsplit[2]) & bam_file$pos<as.numeric(regionsplit[3])-bam_file$qwidth,]
              #limit bam file to region and plot reads
              cpos=c()
              total_counts=0
              message(paste("\nProcessing bam file #",ib))
              pb <- txtProgressBar(style=3)
              for(i in 1:nrow(theregion)){
                setTxtProgressBar(pb, i/nrow(theregion))
                #numreads<-as.numeric(strsplit(theregion[i,]$qname,'-')[[1]][2])
                numreads<-theregion[i,]$counts
                cpos<-c(cpos,rep(theregion[i,]$pos,numreads))
                total_counts<-total_counts+numreads
              }
              if(nFlag) ncounts[ib]<-total_counts
              message(paste('Num reads:',total_counts,' Unique positions:',length(unique(cpos))))
              histo<-hist(cpos,breaks=length(unique(cpos))*howsmooth,plot=FALSE)
              denso<-density(cpos,n=length(unique(cpos))*howsmooth)
              multiplier<-histo$counts[1]/histo$density[1]
              yseries[[ib]]<-denso$y*1000000/ncounts[ib]*multiplier
              xseries[[ib]]<-round(denso$x)
              
              if(ib>1){
                a<-testt
                testt<-merge(a,data.frame(cbind(xseries[[ib]],yseries[[ib]])),by=1,all=TRUE)
              } else {
                testt<-data.frame(cbind(xseries[[ib]],yseries[[ib]]))
              }
            } 
            
            #testt[is.na(testt)]<-0
            colnames(testt)<-c("position",samplenames)
            testt.m<-melt(testt,id="position")
            colnames(testt.m)<-c("position","Condition","value")
            gp<-ggplot(data=testt.m[!is.na(testt.m$value),],aes(x=position,y=value,color=Condition))+geom_line(size=2)+
              xlab(paste("Chromosome position (",regionsplit[1],")",sep=""))+ylab("Reads per million")+
              scale_color_brewer(palette="Set1")+theme_bw(base_size=18)
            print(gp)
            return(list(xseries,yseries))
          })