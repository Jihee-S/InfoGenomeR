##args <- commandArgs(TRUE)
args=c()
args[1]="SVs.CN_opt";
########
CN_MAX=200;
IP_output=read.table("copy_numbers.CN_opt", header=T,stringsAsFactors=F)
levels(IP_output[,1])[levels(IP_output[,1])=="23"]<-"X"
IP_output$modal_cn[is.na(IP_output$modal_cn)]=round(IP_output$raw_expected_cn[is.na(IP_output$modal_cn)])
IP_output$modal_cn[IP_output$modal_cn>CN_MAX]=CN_MAX;

segmean=read.table("copy_numbers",stringsAsFactors=F)
CN=IP_output[,18];
segmeanCN=cbind(segmean, CN);
write.table(segmeanCN, "copy_numbers.CN", quote=F, row.names=F, col.names=F, sep="\t")
t=segmeanCN
d=data.frame(chrom=t[,2],key=0,node=0, degree=0)
dindex=1
for (i in 1:nrow(t)){
	d[dindex,1]=t[i,2]
	d[dindex,2]=t[i,3]
	d[dindex,3]=dindex
	d[dindex,4]=t[i,13]
	dindex=dindex+1;
	d[dindex,1]=t[i,2]
	d[dindex,2]=t[i,4]
	d[dindex,3]=dindex
        d[dindex,4]=t[i,13]
	dindex=dindex+1;
 }

d$SV_CN=0;


SV=read.table(args[1],stringsAsFactors=F)
for(i in 1:nrow(SV)){
	
	d1=d[d$chrom==as.character(SV[i,2])&d$key==SV[i,3],3]
	d2=d[d$chrom==as.character(SV[i,4])&d$key==SV[i,5],3]
	if(SV[i,15]!=0){
	for(SV_CN in 1:SV[i,15]){
	d[d$chrom==as.character(SV[i,2])&d$key==SV[i,3],4]=d[d$chrom==as.character(SV[i,2])&d$key==SV[i,3],4]-1
        d[d$chrom==as.character(SV[i,2])&d$key==SV[i,3],5]=d[d$chrom==as.character(SV[i,2])&d$key==SV[i,3],5]+1

	d[d$chrom==as.character(SV[i,4])&d$key==SV[i,5],4]=d[d$chrom==as.character(SV[i,4])&d$key==SV[i,5],4]-1
        d[d$chrom==as.character(SV[i,4])&d$key==SV[i,5],5]=d[d$chrom==as.character(SV[i,4])&d$key==SV[i,5],5]+1

	}
	}
#	if(d[973,4]<0){print(i)};
}



r_edge_start=2

for(i in c(1:22,"X")){

	r_edge_end=r_edge_start+nrow(d[d$chrom==i,])-4

	if(r_edge_start<=r_edge_end){

		edge_index=r_edge_start;
		while(edge_index<=r_edge_end){

			if(min(d[edge_index,4],d[edge_index+1,4])!=0){
				for(j in 1:min(d[edge_index,4],d[edge_index+1,4])){
		
					d[edge_index,4]=d[edge_index,4]-1;
					d[edge_index+1,4]=d[edge_index+1,4]-1;
		
				}
			}
		edge_index=edge_index+2;
		}


	}
	r_edge_start=r_edge_end+4;
}

############### IF the number of odd edges are not even, error message will occur######################

heindex=1

d_hidden= d[d$degree!=0,]
##################################
hidden_node= max(d$node)+1
##################################
for(i in 1:(nrow(d_hidden))){
	while(d_hidden[i,4]>=0){
                        d_hidden[i,4]=d_hidden[i,4]-1;
                        heindex=heindex+1;
		}

}


################################################

############### IF the number of odd edges are not even, error message will occur######################

### 오일러 사이클을 위한 것. 오일러 경로를 뽑기 위해선 하나의 edge를 제거###


################################# for additional search of SVs ############
d_copy=d;

d=d[d$degree!=0,]
write.table(d, "unsatisfied.list", quote=F, sep="\t", row.names=F)
#################
##############
d_dup=d[0,];
for(i in c(1:22,"X")){
	d_tmp=d_copy[d_copy$chrom==i,];
	if(nrow(d_tmp)>=3){
		d_tmp=d_tmp[-c(1,nrow(d_tmp)),]
		d_dup=rbind(d_dup,d_tmp);
	}
}
d_dup=d_dup[d_dup$degree!=0,];

write.table(d_dup, "unsatisfied.list.no_telomere_ends", quote=F, sep="\t", row.names=F, col.names=F)
##############################
