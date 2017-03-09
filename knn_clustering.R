#############################################################################################################################	
# Description: k nearest neighbor iterative custering approach for microscope data
# Author: Tomas William Fitzgerald	
# Date: 08/04/2017
# Depends: FNN

# function to define how data points are combined
# mean of positions, error adjusted based on max errors, distance
# @ single collasped event
knn.collapse <- function(data, inds) {
		event = data[inds,][1,]
		event$xnm = mean(data$xnm[inds])
		event$ynm = mean(data$ynm[inds])
		event$znm = mean(data$znm[inds])
		event$npoints = sum(data$npoints[inds])
		event$locprecnm = max(c(data$locprecnm[inds], max(dist(data[inds,c("xnm", "ynm", "znm")]))))
		event$ids = paste(inds, collapse="-")
return(event)
}

# function to find pairs of nearest neighbors that fall within the error zones
# @ returns both singleton data points and paired data sets 
knn.pairs <- function(data) {
	k = knn(data[,c("xnm", "ynm", "znm")], data[,c("xnm", "ynm", "znm")], rep(1, nrow(data)), k=2)
	indices = attr(k, "nn.index"); dist = attr(k, "nn.dist")
	kds = cbind(1:nrow(dist), indices[,2], dist[,2]); 
	kidx = kds[kds[,3]<data$locprecnm,]
	if(is.null(dim(kidx))) { kidx = cbind(kidx[1], kidx[2], kidx[3]) }
	t = kidx[,1]; 
	kidx[kidx[,1]>kidx[,2],1] = kidx[kidx[,1]>kidx[,2],2]; 
	kidx[t>kidx[,2],2] = t[t>kidx[,2]]
	ids = paste(kidx[,1], kidx[,2], sep="-"); 
return(list("singles"=data[kds[kds[,3]>=data$locprecnm,1],], "dist.index"=kidx[!duplicated(ids),]))
}

# function to collaspe pairs
# all pairs with any overlapping ids are collapsed
# @ set of collapsed paired events
knn.collapse.pairs <- function(data, kps) {
	rc = list()
	kidx = kps$dist.index
	if(length(kidx)>0) {
		if(is.null(dim(kidx))) { kidx = cbind(kidx[1], kidx[2], kidx[3]) }
		pb <- txtProgressBar(min = 0, max = nrow(kidx), style = 3);
		for(x in 1:nrow(kidx)) {
			inds = c(kidx[x,1:2])
			event = knn.collapse(data, inds)
			rc[[x]] = event
			setTxtProgressBar(pb, x)
		}
		close(pb);
	}
	print("wait...")
return(do.call('rbind', rc))
}

# function to itratively collapse all paired events into clusters
# any overlapping ids across collasped pair events are found
# new events are created based all overlapping events from the original data
# @ set of fully collapsed events
knn.pairs.exhaust <- function(data, kcp) {
	nscs = list()
	if(!is.null(kcp)) {
		sids = table(unlist(strsplit(kcp$ids, "-")))
		search_ids = names(sids[sids>1])
		if(length(search_ids)>1) {
			pb <- txtProgressBar(min = 0, max = length(search_ids), style = 3);
			for(x in 1:length(search_ids)) {
				s1 = kcp[grep(paste("-", search_ids[x], "$", sep=""), kcp$ids),]
				s2 = kcp[grep(paste("^", search_ids[x], "-", sep=""), kcp$ids),]
				sk = rbind(s1, s2)
				inds = as.numeric(unique(unlist(strsplit(sk$ids, "-"))))
				event = knn.collapse(data, inds)
				nscs[[x]] = event
				setTxtProgressBar(pb, x)
			}
		close(pb);
		}
		print("wait...")
	}
return(do.call('rbind', nscs))
}

#### MAIN
 
require(FNN)

# read speific data set, sort and add some new ids
kdata = read.table("6_160620_NW_mEGFP-NCAPH_GFP-NB-sortase-AF647_1_driftc_sml.csv", header=T, sep=",")
kdata = kdata[order(kdata$xnm, kdata$ynm, kdata$znm),]
kdata$npoints = rep(1, nrow(kdata))
kdata$real_id = 1:nrow(kdata); 
kdata$ids = kdata$real_id

# Iteratively cluster events until full closure is achieved 
# NB. the remianing singletons are added back at each iterative - may now fall within a new events error
while(nrow(kdata)<n) {
	n = nrow(kdata)
	kps = knn.pairs(kdata)
	kcp = knn.collapse.pairs(kdata, kps)
	kpe = knn.pairs.exhaust(kdata, kcp)
	kdata = rbind(kps$singles, kpe)
	kdata = kdata[order(kdata$xnm, kdata$ynm, kdata$znm),]
	print(n)
}
