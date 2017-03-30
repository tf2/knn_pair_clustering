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
    #event$ids = paste(inds, collapse="-")
    inds1 = c(as.numeric(unique(unlist(strsplit(as.character(data$real_ids[inds]), "-")))), as.numeric(unique(unlist(strsplit(as.character(data$ids[inds]), "-")))))
    event$real_ids = paste(unique(inds1[order(inds1)]), collapse="-")
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
    d.i = kidx[!duplicated(ids),];
    if(is.null(dim(d.i))) { d.i = cbind(d.i[1], d.i[2], d.i[3]) }
    s.ids = data[kds[kds[,3]>=data$locprecnm,1],];
    s.ids = s.ids[!s.ids$ids%in%data$real_ids[d.i[,1]],];
    s.ids = s.ids[!s.ids$ids%in%data$real_ids[d.i[,2]],]
    return(list("singles"=s.ids, "dist.index"=d.i))
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
            rc[[x]] = inds
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
    nid = NULL
    oc = NULL
    if(!is.null(kcp)) {
        t = table(c(kcp))
        cpairs = names(t[t>1])
        inds = sapply(1:nrow(kcp), function(x) sum(kcp[x,]%in%cpairs))
        ic = kcp[inds>0,]
        oc = kcp[inds==0,]
        pin = 1
        # pb <- txtProgressBar(min = 0, max = nrow(ic), style = 3);
        while(length(ic)>0) {
            if(is.null(dim(ic))) { ic = cbind(ic[1], ic[2]) }
            w = which(ic[,1]==ic[1,1] | ic[,2]==ic[1,1] | ic[,1]==ic[1,2] | ic[,2]==ic[1,2] )
            while(length(w)!=length(which(ic[,1]%in%ic[w,1] | ic[,1]%in%ic[w,2] | ic[,2]%in%ic[w,1] | ic[,2]%in%ic[w,2]))) {
                w = which(ic[,1]%in%ic[w,1] | ic[,1]%in%ic[w,2] | ic[,2]%in%ic[w,1] | ic[,2]%in%ic[w,2])
            }
            ns = unique(c(ic[w,]))
            nscs[[pin]] = ns
            nid = c(nid, paste(ns, collapse="-"))
            ic = ic[-w,]
            pin = pin +1
            #setTxtProgressBar(pb, pin-nrow(ic))
        }
        #close(pb);
    }
    icls = do.call(rbind, lapply(nscs, function(x) knn.collapse(data, x)))
    ocls = NULL
    if(is.null(dim(oc))) { oc = cbind(oc[1], oc[2]) }
    if(length(oc)>0) {
        for(x in 1:nrow(oc)) {
            ocls = rbind(ocls, knn.collapse(data, oc[x,]))
        }
    }
    return(rbind(ocls, icls))
}

#### MAIN
require(FNN)

# read speific data set, sort and add some new ids
data = read.table("6_160620_NW_mEGFP-NCAPH_GFP-NB-sortase-AF647_1_driftc_sml.csv", header=T, sep=",")
kdata = data[1:100,]
kdata = kdata[order(kdata$xnm, kdata$ynm, kdata$znm),]
kdata$npoints = rep(1, nrow(kdata))
kdata$real_ids = 1:nrow(kdata);
kdata$ids = kdata$real_ids
n = nrow(kdata)+1

# Iteratively cluster events until full closure is achieved
# NB. the remianing singletons are added back at each iterative - may now fall within a new events error
while(nrow(kdata)<n) {
    n = nrow(kdata)
    kps = knn.pairs(kdata)
    kcp = knn.collapse.pairs(kdata, kps)
    kpe = knn.pairs.exhaust(kdata, kcp)
    kdata = rbind(kps$singles, kpe)
    kdata = kdata[order(kdata$xnm, kdata$ynm, kdata$znm),]
    kdata$npoints = unlist(lapply(strsplit(kdata$real_ids, "-"), length))
    print(n)
}




