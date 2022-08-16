#look at data from COVID19 using Marc's master files
# this script takes each time point, and compares it to each other time point for that volunteer
# and extracts all the TCRs which change between time points, by an amount greater than would be expected by chance.
#the script  plots this data for all volunteers and all time points. It also saves all the TCRs (as a combination of V, J and nucleeotide junction.

#plot COVID response
library(mgcv)

#on Linux
drive <- "/media/benny/data/"
#on Windows
drive<-"C:/users/Benny Chain/"
#this is the path to the master files on my computer
input_data<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\data_for_paper\\"
dir<-dir(paste0(drive,input_data))

#folder for collecting output
#folder<-"Dropbox/R_temp/01_09_2021/"
#folder for plots
folder_plots<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\Figs\\Pairwise\\"
#folder for data output
folder_data<-"Dropbox\\Temp (1)\\Papers\\COVIDsortium\\TCR_paper\\Data_for_paper\\"

#I run the script separately for alpha and beta chains; I have commented out beta chain at the moment
#alpha
load(paste0(drive,input_data,"data_tsv_alpha.Rdata"))
chain<-"alpha"
data<-data_tsv_alpha
rm(data_tsv_alpha)

#Beta
load(paste0(drive,input_data,"data_tsv_beta.Rdata"))
chain<-"beta"
data<-data_tsv_beta
rm(data_tsv_beta)

#analyse individual changes in TCRS
#list of ids
HCW_all<-unique(data$ID)

#convert time points to numeric and add a column to data
tp_week<-rep(100,dim(data)[1])
tp_week[which(data$timepoint=="BL")]<-0
tp_split<-strsplit(data$timepoint,"U")
i_spl<-which(sapply(tp_split,length)==2)
week<-as.numeric(sapply(tp_split[i_spl],function(x){x[2]}))
tp_week[i_spl]<-week
data$tp_week<-tp_week#

#square plots
par(mar=c(5, 10, 4, 2) + 0.1)

#####################################################################
#RESTART
#####################################################################

#collects all  TCRs which change at each time point for each volunteer, as a named list
TCR_change_HCW<-list()
#p is index of HCW names in alphabetical order
#v is number of volunteer
p<-2

#Note that PCR negative controls are 17, 42, 84,107,158,186
v<-17
v<-35
v<-107
############################################################
#if you want to run controls only
#for (v in c(17,42,84,107,158,186)){
###########################################################
#which(HCW_all==299)
#loop counter
n<-1
#loop for all volunteers
p<-1
#flag for whether to do plots
plot<-"TRUE"
#plot<-"FALSE"
for ( p in 1:length(HCW_all)){
  v<-HCW_all[p]
#select the data for this individual
data_i_v<-data[data$ID==v,]

#exclude any time point post week 5
data_i<-data_i_v[data_i_v$tp_week<5,]

#list of time points
#time<-unique(data_i_v$tp_week)
time<-unique(data_i$tp_week)
#print the volunteer number and the time points available for each week
#this is just useful so you know where you are in teh loop while its running
cat (v, "    ", time, "\n")
#only process volunteers with more thna one time point
if (length(time)>1){

#variables to collect up and down regulated
TCR_change_all<-c()

#cycle through combinations of time points
i<-1
j<-2
#note the use of a nested loop with the second counter staring at i+1, to compare all time points
for ( i in 1 : (length(time)-1)){
  for ( j in (i+1) : length(time)){
    #first time point
    data_i_1<-data_i[data_i$tp_week==time[i],]
    #list of TCRs, combinining V, J and junction
    TCR_names_1<-unique(paste(data_i_1$junction,data_i_1$v_call,data_i_1$j_call,data_i_1$decombinator_id, sep="_"))

    #second time point
    data_i_2<-data_i[data_i$tp_week==time[j],]
    #list of TCRs, combinining V, J and junction
    TCR_names_2<-unique(paste(data_i_2$junction,data_i_2$v_call,data_i_2$j_call,data_i_2$decombinator_id,sep="_"))

    ########################################
    #list of TCRs from both time points, combinining V, J and junction
    TCR_names<-unique(c(TCR_names_1,TCR_names_2))
    #set up matrix to receive all the TCRs from these two time points ;
    all_TCRs<-matrix(rep(0,2*length(TCR_names)),nrow=length(TCR_names))
    rownames(all_TCRs)<-TCR_names

    #load data into all_TCRs from these two time points
    all_TCRs[TCR_names_1,1]<-data_i_1$duplicate_count
    all_TCRs[TCR_names_2,2]<-data_i_2$duplicate_count

    #normalisation - convert absolute counts to TCRs/million in that sample
    norm_1<-1E6*all_TCRs[,1]/sum(all_TCRs[,1])
    norm_2<-1E6*all_TCRs[,2]/sum(all_TCRs[,2])

    sum_1<-sum(all_TCRs[,1])
    sum_2<-sum(all_TCRs[,2])
    sum_m<-min(sum_1,sum_2)


    min<-unique(sort(c(norm_1,norm_2)))[2]
    #replace zeros for plotting purposes with min/2
    i_1<-which(norm_1==0)
    norm_1[i_1]<-min/2

    i_2<-which(norm_2==0)
    norm_2[i_2]<-min/2

  x<-log2(norm_1)
  y<-log2(norm_2)

  #error margins - I am calculating the confidence intervals that, if you observe a TCR m times in a sample,
  #you will observe it n times in a second sample, if the only process at work is random Poisson sampling.
  #I multiply by a scalar because all the raw counts are normalised ot TCR counts/million.

  sum_m<-min(sum_1,sum_2)
  scalar_2<-1E6/length(TCR_names_1)
  scalar_1<-1E6/length(TCR_names_2)
  scalar<-max(scalar_1,scalar_2)
  poiss_l<- qpois(0.0001,c(1,2^(0:4)),lower.tail=FALSE)



#####################################################################################################
  #for plotting purposes only, not calculation or saving, subsample to 50,000 dots per plot
  if(plot==TRUE){
  if (length(x)>50000){i_sample<-sample(1:length(x),50000)} else {i_sample<-1:length(x)}
  x_ss<-x[i_sample]
  y_ss<-y[i_sample]

  #add a small ranodm noise so teh dots don;t overlie each other
  x_jit<-rnorm(length(x_ss),x_ss,abs(x_ss)/50)
  y_jit<-rnorm(length(y_ss),y_ss,abs(y_ss)/50)

  #plot the freqeuncy of each TCR at each of teh two time points
  plot(x_jit,y_jit,pch=19,cex=0.5,main=paste(v,time[i],time[j],chain),xlim= c(0,12),ylim=c(0,12),xaxt="n",yaxt="n",ylab="TCRs/million",xlab="TCRs/million",cex.lab=2)
  #plot(x,y,pch=19,cex=0.5,main=paste(v,time[i],time[j],chain),xlim= c(0,12),ylim=c(0,12),xaxt="n",yaxt="n",ylab="TCRs/million",xlab="TCRs/million")

  axis(1, at = c(min(x_ss),4,6,8,10),labels = c(0,2^c(4,6,8,10)),cex.axis=1.5)
  axis(2, at = c(min(y_ss),4,6,8,10),labels = c(0,2^c(4,6,8,10)),las =2 ,cex.axis=1.5)
  mx<-mean(x)
  my<-mean(y)

  #y_pois<-log2(scalar*poiss_l)
  #add the error limits as calcualted above using Poisson distribution
  #can correct the error margins by setting r = my/mx. But seems to work better setting r=1
  r<-1
  points(log2(c(min/4,(2^c(0:4,4))*1E6/sum_m)),log2(c(poiss_l*scalar,2^12)),col="blue",pch=19,type="l", lty = 2,lwd=2)
  points(log2(c(poiss_l*scalar,2^12)),log2(c(min/4,(2^c(0:4,4))*1E6/sum_m)),col="blue",pch=19,type="l", lty = 2,lwd=2)

  #plot a straight line, whose slope is the relative mean of the tiem time points
  abline(0,my/mx)

#save the plots using a homemade function called imageSave, which makes saving plots
#much more straightforward
  imageSave(file=paste0(drive,folder_plots,v,"_",time[i],"_",time[j],"_",chain,".png"))

  }
  #calculate which points are outside error limits
  # note that correct margin by mx/my
  #this sets the polygons  in which to count the TCRs
  
  bnd_up<-cbind(log2(c(min/4,(2^c(0:4,4))*1E6/sum_m,min/4)),log2(c(poiss_l*scalar,2^12,2^12)))
  bnd_down<-cbind(log2(c(poiss_l*scalar,2^12,2^12)),log2(c(min/4,(2^c(0:4,4))*1E6/sum_m,min/4)))


  #bnd_up<-cbind(c(1:6,6,1),c(y_pois*(r),12,12))
  #bnd_down<-cbind(c(y_pois*(1/r),12,12),c(1:6,6,1))

  #the function in.out from the package mgcv counts the number of points in a quadrant, given the
  # equation which defines the polygon
  i_up<-which(in.out(bnd_up,cbind(x,y)))
  i_down<-which(in.out(bnd_down,cbind(x,y)))

  #collect all TCRs up regulated; include flag which is +1 for up, -1 for down
  TCR_up<-cbind(x[i_up],y[i_up],rep(time[i],length(x[i_up])),rep(time[j],length(x[i_up])),rep(1,length(x[i_up])))
  TCR_down<-cbind(x[i_down],y[i_down],rep(time[i],length(x[i_down])),rep(time[j],length(x[i_down])),rep(-1,length(x[i_down])))
#TCR_change_all<-c()
  TCR_change_all<-rbind(TCR_change_all,TCR_up,TCR_down)
  colnames(TCR_change_all)<-c("Freq 1","Freq 2", "Time point 1", "Time point 2","Change" )
  #collect all the number of points beyond limits
  #t_up<-c(t_up,length(i_up))
  #t_down<-c(t_down,length(i_down))
  cat(paste(length(i_up),length(i_down)),"\n")
  }
    }#end of timepoint loops
      } else { #end of if statement
          TCR_change_all<-NA
            }

#once cycled though all combination of time points for that volunteer, add all TCRs to list
TCR_change_HCW[[n]]<-TCR_change_all
names(TCR_change_HCW)[[n]]<-v
#cat(v, "\n")
if(is.na(TCR_change_all)){cat( v ," has only one timepoint","\n")}
n<-n+1
}#end of HCW loop

#once the loop has gone through all individual volunteers save the output; either as alphe or beta
#depending on what was decided on lines 20-25.
if (chain=="alpha"){
  TCR_change_HCW_a<-TCR_change_HCW
  save(TCR_change_HCW_a, file = paste0(drive, folder_data,"TCR_change_all_alpha.RData"))
    }

if (chain=="beta"){
  TCR_change_HCW_b<-TCR_change_HCW
  save(TCR_change_HCW_b, file = paste0(drive, folder_data,"TCR_change_all_beta.RData"))
}

#load(file = paste0(drive, folder,"TCR_change_all_alpha.RData"))
#load(file = paste0(drive, folder,"TCR_change_all_beta.RData"))


