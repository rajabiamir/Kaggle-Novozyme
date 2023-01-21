##########libraries
library(tidyverse)


####to read the data 
dataset<- read.csv("train.csv")

correction<-read.csv("train_updates_20220929.csv")
View(correction)
####to clean the updated dataset and make it ready to be merged to the main dataset

correction%>%filter(protein_sequence=="")%>%as_tibble()->should_be_deleted
should_be_deleted

should_be_transposed<-anti_join(correction,should_be_deleted)
as_tibble(should_be_transposed)


##to delete the rows should be deleted plus the rows that needed to be updated
dataset<-as_tibble(anti_join(dataset,correction,by="seq_id"))


#to update the 25 row that needed to be modified
dataset<-as_tibble(full_join(dataset,should_be_transposed,"seq_id"))
View(dataset)

#to clean the columns
dataset$tm.y<-ifelse(is.na(dataset$tm.y),0,dataset$tm.y)
dataset$tm.x<-ifelse(is.na(dataset$tm.x),0,dataset$tm.x)
dataset%>%mutate(tm = tm.x+tm.y)%>%as_tibble()->dataset

dataset$pH.y<-ifelse(is.na(dataset$pH.y),0,dataset$pH.y)
dataset$pH.x<-ifelse(is.na(dataset$pH.x),0,dataset$pH.x)
dataset%>%mutate(pH = pH.x+pH.y)%>%as_tibble()->dataset
dataset%>%select(-tm.x,-tm.y,-pH.x,-pH.y)->dataset


dataset$protein_sequence.y[is.na(dataset$protein_sequence.y)]<-""
dataset$protein_sequence.x[is.na(dataset$protein_sequence.x)]<-""
dataset$protein_sequence<-str_c(dataset$protein_sequence.x,'',dataset$protein_sequence.y)

dataset$data_source.y[is.na(dataset$data_source.y)]<-""
dataset$data_source.x[is.na(dataset$data_source.x)]<-""
dataset$data_source<-str_c(dataset$data_source.x,'',dataset$data_source.y)

dataset%>%select(-protein_sequence.x,-protein_sequence.y,-data_source.x,-data_source.y)->dataset
dataset%>%select(seq_id,protein_sequence,pH,tm,data_source)->dataset
View(dataset)


#write.csv(dataset,"training_dataset.csv")
the_data <- read.csv("training_dataset.csv")


fireprot<- read.csv("fireprotdb_results.csv")
fireprot<-as_tibble(fireprot)

test_set<-read.csv("test.csv")
test_set<-as_tibble(test_set)

my_elu<-function(f){
  if (f<0.25 | f>0.75) {f<-0}
  else{f<-1}
  return(f)
}


#the wild type sequence of the hypothetical protein
wt<-("VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQRVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGTNAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKALGSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK")
wt<-str_split(wt,"")
##wt_vector<-unlist(wt)


#ihem<-wt_vector==prisoner2
#ihem

#the function that specifies at which amino acid index the point mutation has occured (only detects the first one)
#for(i in 1:length(ihem)) {
 # if (ihem[i]==TRUE){}
  #else {print(i)}
#}



#seqvec is the protein sequence of the test set splitted by ""
seqvec<-test_set$protein_sequence[1:2413]
seqvec<-str_split((seqvec),"")
seqvec

###function to find the point mutation from the wild type sequence 
#dualVector is the vector of T n F of the resulting wt[[]]==maVector

mutindex<-vector(mode="integer",length = 2413L)

for(i in 1:length(seqvec)){
    dualVector<- wt[[1]]==seqvec[[i]]
    for(j in 1:length(dualVector)){
  if (dualVector[j]==TRUE){}
  else {mutindex[i]<-j
  break}
      }
}
mutindex%>%as_tibble()%>%print(n=10000)
mutindex[1170]



#####function to find the single amino acid deletion from the wt---it gives 1 to any row that has deletion
deletion_ind_finder<-function(seqvec){
delindex<-vector(mode="logical",length = as.integer(length(seqvec)))
delindex
for (i in 1:length(seqvec)){
  Tvec<-(wt[[1]]==seqvec[[i]])
  if(sum(Tvec)>=220){}
  else{delindex[i]<-1}
}
}

deletion_ind_finder(seqvec)
delindex
delindex%>%as_tibble()%>%View()

mutindex
test_set%>%mutate(delindex)%>%mutate(mutindex)

##function to find the amino acid that is mutated TO

mutletter<-vector(mode="character",length=2413L)
mutletter<-c()
for(i in 1:length(mutindex)){
  if(mutindex[i]>0){mutletter[i]<-seqvec[[i]][mutindex[i]]}
  else{mutletter[i]<-""}
  
}

test_set%>%mutate(mutletter)%>%mutate(mutindex)%>%mutate(deletion=delindex)->test_set
View(test_set)

###to find the wt of the mutated

wt_letter<-c()

use_it_on_wt<-c()
use_it_on_wt<-mutindex-delindex
use_it_on_wt[1171]

for(i in 1:length(use_it_on_wt)){
  if(use_it_on_wt[i]>0){wt_letter[i]<-wt[[1]][use_it_on_wt[i]]}
  else{wt_letter[i]<-""}
  
}
wt_letter






##to find the location of the point mutation(0 means it's on 25% of the n or the c terminal)
mut_location<-c()
for(i in 1:length(mutindex)){
mut_location[i]<-mutindex[i]/length(seqvec[[i]])
mut_location[i]<-my_elu(mut_location[i])
}
mut_location


##to clean the mutletter of the deletion (DO IT AT THE LAST STEP)(Z is the single amino acid deletion)

for(i in 1:2413){
  if(delindex[i]==0){}
  else{
  mutletter[i]<-"Z"
  }
}

mutletter


test_set%>%select(-mutletter,-deletion)%>%mutate(mutletter)%>%mutate(wt_letter)%>%mutate(mut_location)->updated_test_set
View(test_set)

#write.csv(updated_test_set,"updated_test_set.csv")

####to clean up the fireprot dataset of the duplicate and ddG
distinct(fireprot,experiment_id)
fireprot%>%group_by(experiment_id)%>%summarise(count=n())%>%filter(count>1)%>%.$experiment_id->fireprot_duplicate
fireprot%>%filter(!(experiment_id%in%fireprot_duplicate))%>%View()
  
fireprot%>%distinct(experiment_id,.keep_all = TRUE)%>%filter(!is.na(dTm))%>%View()

####to find the same protein in the kaggle dataset
which_i_is_it<-vector(mode= "list",length = length(kaggle_data_splet))

which_i_is_it

kaggle_data_splet<-str_split(the_data$protein_sequence,"")

for (i in 1:length(kaggle_data_splet)){
  for (j in 1:length(kaggle_data_splet)){
    haha<-(!kaggle_data_splet[[i]]==kaggle_data_splet[[j]])
      if(sum(haha)>0 && sum(haha<2)){which_i_is_it[[i]]<-j}
        else{}
  }
}

!kaggle_data_splet[[2000]]==kaggle_data_splet[[2]]

which_i_is_it[2000]


the_data%>%filter(!data_source=="doi.org/10.1038/s41592-020-0801-4")->the_data_2
the_data_2%>%arrange(pH,data_source)->the_data_2

the_data_2_2<-str_split(the_data_2$protein_sequence,"")  
  

####to get as much data as possible for the single point mutation from the dataset provided by kaggle
ccc[-(which(ccc %in% ddd))]


ccc<-c(1:4456)
ccc[-(which(ccc %in% ddd))]

cccc<-c(1:4456)

this_wt<-c()
this_mt<-c()
this_dt<-c()

for (i in ccc){
  for(j in cccc){
    kj<-c()
    ccccc<-the_data_2$data_source[i]==the_data_2$data_source[j]
      if(ccccc==TRUE){kj<-(!the_data_2_2[[i]]==the_data_2_2[[j]])}else{}
    if(sum(kj)==1){
          for(k in 1:length(kj)){if(kj[k]==FALSE){} else{this_wt<-append(this_wt,the_data_2_2[[i]][k])}}} else{}
    if(sum(kj)==1){ccc<-ccc[-(which(ccc=="j"))]} else{}
}
}
this_wt
  this_mt

rm(i,j,k,ccc,cccc,ccccc,kj)




ccc<-c(1:4456)
cccc<-c(1:4456)

this_mt<-c()

this_mt_vt<-vector(mode="integer",length= 4456L)
this_wt_vt<-vector(mode="integer", length=4456L)
for (i in ccc){
  for(j in cccc){
    kj<-c()
    ccccc<-the_data_2$data_source[i]==the_data_2$data_source[j]
    if(ccccc==TRUE){kj<-(!the_data_2_2[[i]]==the_data_2_2[[j]])}
    else{}
    if(sum(kj)==1){
      for(k in 1:length(kj)){if(kj[k]==FALSE){}else{this_wt_vt[i]<-k}}
    }else{}
    if(sum(kj)==1){ccc<-ccc[-j]}else{}
  }
}


kj<-c()
this_wt_vt




to_probe_for_del<-str_split(the_data$protein_sequence,"")
length(to_probe_for_del[[1]])

for(p in 1:28981){
if(length(to_probe_for_del[[p]])==340){print(p)}
}          

to_probe_for_del[[1]]==to_probe_for_del[[499]]                    
