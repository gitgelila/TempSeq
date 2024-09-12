                        #### Temporal Sequencing of Documents ####

Data Source: https://cran.r-project.org/web/packages/sotu/sotu.pdf

                                        #R SOURCE CODE 


#IN THIS PROGRAM, WE SELECT 100 RANDOM SAMPLES WHERE EACH SAMPLE CONSISTS OF A SEQUENCE OF
#TEN ADDRESSES FROM STATE OF THE UNION ADDRESS (SOTU). THE DATES OF THE ADDRESSES ARE IN A
#CONSECUTIVE SEQUENTIAL ORDER 23 YEARS APART, RANGING FROM 1790 TO 2020. THERE ARE A TOTAL OF 240 ADDRESSES.
#WE RANDOMLY PERMUTE THE TEMPORAL ORDER OF THE SEQUENCE OF THE TEN
#ADDRESSES FROM EACH SAMPLE, WHICH ARE THE INPUTS, AND THE OUTPUTS ARE THE PREDICTIONS OF THE
#TRUE TEMPORAL ORDER OF EACH ONE OF THE PERMUTED SAMPLES.
#########################################################################################
##########################################################################################               

                              ######### PRE-PROCESSING OF TEXTUAL DATA ########
#Input: From the SOTU corpus, read all lines (each line being an Address) into an array
#SofUtext = read lines from file=SOTU into an array
SofUtext<-readLines("/home/SOTU", n=-1)
L<-length(SofUtext); SofU2gramtxt<-SofUtext

for(i in 1:L){
    s<-SofU2gramtxt[i]
    w<-strsplit(s, " ", fixed = TRUE)[[1L]]  ## Split into words:
    SofU2gramtxt[i]<-paste(vapply(ngrams(w, 1L), paste, "", collapse = "_"), sep=" ", collapse = " ")
}

SofUVCDocTerm<-as.matrix(DocumentTermMatrix(VCorpus(VectorSource(SofU2gramtxt))))
SofUVCDocTerm<-SofUVCDocTerm[, -(1:1254)]  # Remove corrupted or non-words from the word column list

indx<-which(colSums(SofUVCDocTerm)>=2)
SofUVCDocTerm<-SofUVCDocTerm[,indx]

FreqTableSofU<-t(SofUVCDocTerm)



##################### SAMPLE 100 TIMES A SEQUENCE OF 10 ADDRESSES, DATED TO BE 23 YEARS APART. ENSURE THAT THERE ARE NO
#DUPLICATED SAMPLES BY SAMPLING 200 TIMES, THEN REMOVE THE DUPLICATED TERMS, AND THEN SAVE THE FIRST 100 UNDUPLICATED
#SAMPLES. THE UNIQUE 100 SAMPLES ARE STORED IN MATRIX "SampleMat" #####################



colnames(FreqTableSofU)=1:dim(FreqTableSofU)[2]

lengthofstring<-10;  numofspeech<-dim(FreqTableSofU)[2]  #240 columns (Addresses)
numSample<- 200; SampleMat<- matrix(0,numSample,lengthofstring)
count=1

while(count<=numSample){
         initialsample<-sample(1:numofspeech,1)
         SampleSofU<-c(initialsample, initialsample+(1:(lengthofstring-1))*23)%%(numofspeech)
         indx<-which(SampleSofU==0)
         if(length(indx)>0){
             SampleSofU[indx]=240
         }
         SampleSofU<-sort(SampleSofU) 
         SampleMat[count,]<-SampleSofU
         count=count+1
}


numSampleNow<-100
SampleMat<-unique(SampleMat)[1:numSampleNow,] #Each row is a unique sequence


                                              ####### MAIN FUNCTION #######
######################  WITHIN EACH FOR LOOP:  EACH ROW OF THE SAMPLE MATRIX "SampleMat" CONSISTS OF A SET OF 10 RANDOM ADDRESSES SEQUENTIALLY ORDERED 23 YEARS APART. IN THE MATRIX "FreqMatrixFinal", THE COLUMNS ARE THE TEMPORALLY ORDERED SAMPLED ADDRESSES, AND THE ROWS ARE THE WORDS WHICH HAVE OCCURRED AT LEAST TWICE ACROSS THE SAMPLED ADDRESSES. THE ELEMENTS OF THE MATRIX ARE THE COUNTS OF WORDS (ROWS) IN THE ADDRESSES (COLUMNS). "Sint" IS A RANDOM PERMUTATION OF THE COLUMNS OF "FreqMatrixFinal". THE INPUT FOR THE SIMULATED ANNEALING ALGORITHM IS THE MATRIX "FreqMatrixFinal" ALONG WITH THE RANDOM COLUMN PERMUTATION "Sint". THE OUTPUT IN THE ROW OF THE MATRIX "SestHopt" IS THE PREDICTED ORDER OF THE 10 ADDRESSES AFTER THE PERMUTATION "Sint". ######################



SestHopt<-matrix(0,numSampleNow,lengthofstring+1) #columns is now 11; optimal bandwidth &length of 10 strings 
for(i in 1:numSampleNow){  
    S<-SampleMat[i,]
    FreqMatrixTemp<-FreqTableSofU[,S]
    IndxRow<- which(rowSums(FreqMatrixTemp)>=2)
    FreqMatrixFinal<- FreqMatrixTemp[IndxRow,]			
    M<-length(S)
    Sint=sample(1:M); print(c("Sint=", Sint))
    SestHopt[i,]<-SimAnnealHopt( FreqMatrixFinal, Sint)
    print(SestHopt[i,])
    x=c(i, SestHopt[i,])
    write(x,file="/home/SestHopt.txt",ncolumns =12, append=TRUE, sep = " ") #Write the final result "x" into a file "SestHopt.txt" . The first column is the index. The second column is the bandwith associated with the predicted temporal order (Section 5.1, equation (4)). The remaining 10 columns are the estimated order of the true temporal (the correct order is 1 to M sequentially, M=10)
}



########### SUB-FUNCTION: COMPUTES THE ELEMENTS REQUIRED TO COMPUTE BANDWIDTH (EQUATION 11 AND EQUATION 12 IN APPENDIX B) ###########


HoptInputFun<-function(){
    nu<-5
    xrand<-rt(10000,nu)   #1000 random variables from a t-distribution 
    IntgKsqrd<-mean(dt(xrand, nu)) #this is the last component of equation 11 in Appendix B (\int K^2(z)dz).
    secMomentSqr<- (nu/(nu-2))^2 #this is the first component of equation 12 in Appendix B ( (\int z^2K(z)dz)^2 )
    return(c(IntgKsqrd,secMomentSqr))
}
    



############# SUB-FUNCTION: COMPUTING THE MEDIAN AVERAGE BANDWIDTH VALUE ACROSS WORDS (SECTION 5.1, EQUATION 3) #############

HoptFun<-function( FreqMatrixFinal, Sint){

    hoptInput<-HoptInputFun()
    d1<- dim(FreqMatrixFinal)[1]
    IntgKsqrd<-hoptInput[1]
    secMomentSqr<-hoptInput[2]
    
    hoptArray<-numeric(d1)
    x<-1:length(Sint)
    N<-colSums(FreqMatrixFinal[,Sint])
    M<-length(Sint)
    
    for(i in 1:d1){
        x2<-x^2
        Trials<-cbind(FreqMatrixFinal[i,Sint],N-FreqMatrixFinal[i,Sint] )
        model<-glm(Trials~x+x2, family = binomial(link="logit")) 
        nu2ndderivArray<- 2*model$coefficients[3]
        Expnu2ndderivArray<-nu2ndderivArray^2
	invVar<-sum( ( model$fitted*(1-model$fitted) )^{-1} )
        hopt<-{ (IntgKsqrd*invVar)/(secMomentSqr*Expnu2ndderivArray*M)}^{1/5}
        hoptArray[i]<-hopt
    }
    return(median(hoptArray))
}


		
########## SUB-FUNCTION: CREATE A NEIGHBOURHOOD FROM THE CURRENT INPUT SEQUENCE BY REVERSING AND/OR MOVING A SUBSEQUENCE OF TERMS (SECTION 5.2) #############
         

CreateneighbourSolNew<-function(S){ 
    M<-length(S)
    anchorIndx<-sample(1:(M-3),1)
    indextake<-anchorIndx:(anchorIndx+3)

    dropIndx<-sample(1:(M-3),1)
    indexput<-dropIndx:(dropIndx+3)

    SNew<-numeric(M)
    SNew[indexput]<-S[indextake]
    SNew[-indexput]<-S[-indextake]

    flip<-sample(c(0,1),1); # to randomly reverse order
    if(flip==1){
        LenI<-length(indexput)
        SNew[indexput]<-SNew[indexput[LenI:1]]
    }
    return(SNew)
}    


########## SUB-FUNCTION: ANNEALING SCHEDULE ##########

CalculateTemp<-function(Constant,i){
	return(Constant^{-i})
}




########## SUB-FUNCTION: A CALL FROM THE SUB-FUNCTION "SIMULATED ANNEALING" ##########

EntropySeqHopt<-function(FreqMatrixFinal,S){ 
        hopt<- -(HoptFun(FreqMatrixFinal, S) ) 
        return(hopt)

}





############### SUB-FUNCTION: SIMULATED ANNEALING ALGORITHM. WHEN SEEN IN CONJUCTION WITH THE SUB-FUNCTION "EntropySeqHopt", THIS ALGORITHM FINDS A PERMUTAION SOLUTION THAT MAXIMIZES THE MEDIAN AVERAGE BANDWIDTH VALUE ACROSS WORDS (SECTION 5.1, EQUATION 4). ######################

SimAnnealHopt<-function(FreqMatrixFinal, Sint){
	Constant<- 1.011
	Scurrent<-Sint
	Sbest<-Scurrent
	i<-1
	while(i<=150){		
      		i<-i+0.5
      		Si<-CreateneighbourSolNew(Scurrent)
      		tempcurrent<-CalculateTemp(Constant, i)
      		EntropySeqSi<-EntropySeqHopt(FreqMatrixFinal, Si); EntropySeqScurrent<-EntropySeqHopt(FreqMatrixFinal, Scurrent)
      		deltaEntropy<- EntropySeqSi-EntropySeqScurrent
      		if(deltaEntropy < 0){
                    Scurrent<-Si
                    if(EntropySeqSi<=EntropySeqHopt(FreqMatrixFinal, Sbest)){
                        Sbest<-Si; 
                    }
      		}
      		else{
                    if( exp(-deltaEntropy/tempcurrent)> runif(1,0,1) ) {
                        Scurrent<-Si; 
                    }
     		}     
        }
        h<- -EntropySeqHopt(FreqMatrixFinal, Sbest)
        return(c(h, Sbest) )
}

