# Autor: Leonna Szangolies, Jonatan Siegmund, Potsdam Institute for Climate Impact Research, 2020
# This function is part of the R package CoinCalc().
# CC.binarize can be used to binarize a vector referring to a given threshold. This threshold can either be a percentile or an absolute value.

CC.binarize <- function(data,ev.def="percentile",thres=0.90,event="higher"){ 
        # ----------------------------------------------- TEST FOR ERRORS IN FUNCTION CALL  --------------------------------------------------------  
        # ---- #1: Time series not binary 
        if(!is.vector(data)){
            print("|-------------    ERROR #1    -------------|")
            print("|     given 'data' is not a vector         |")
            print("|------------------------------------------|")
            return()}
        # ----------------------------------------------- BINARISATION ------------------------------------------------------------------------------
        
        tmp=data
        if(ev.def=="percentile"){
            if(event=="higher"){
                data[tmp>quantile(tmp,thres,na.rm=TRUE)]=1;data[tmp<quantile(tmp,thres,na.rm=TRUE)]=0
            } # end if higher
            if(event=="lower"){
                data[tmp<quantile(tmp,thres,na.rm=TRUE)]=1;data[tmp>quantile(tmp,thres,na.rm=TRUE)]=0
            } # end if higher
        data[tmp==quantile(tmp,thres,na.rm=TRUE)]=1
        } # end if percentile

        if(ev.def=="absolute"){
            if(event=="higher"){
                data[tmp>thres]=1;data[tmp<thres]=0
            } # end if higher
            if(event=="lower"){
                data[tmp<thres]=1;data[tmp>thres]=0
            } # end if higher
            data[tmp==thres]=1
        } # end if absolute
        
        if(ev.def=="timeextrem"){
          period=thres
          rep=floor(length(data)/period)
          for(i in c(0:rep)){
            if((period*i+1)<=length(data)){
            if(event=="higher"){
              data[(period*i+1):min(length(data),(period*(i+1)))]=0
              data[i*period+which.max(tmp[(period*i+1):(period*(i+1))])]=1;
            } # end if higher
            if(event=="lower"){
              data[(period*i+1):min(length(data),(period*(i+1)))]=0
              data[i*period+which.min(tmp[(period*i+1):(period*(i+1))])]=1;
            } # end if higher
            }
          }
        } # end if absolute
        return(data)
        
    } # end function CC.binarize