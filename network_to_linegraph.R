rand_cog_graphs <- function(concepts,prob,n)
{
	N <- length(concepts)
	g<-rgraph(N,n,tprob=prob,mode='graph')
	g
}

mult.net.plot <-function(g,concepts,soc_ppar,do.small=FALSE)
{
	
	N <- length(concepts)
	n <- dim(g)[1]
	if (do.small){
	plot(as.network(soc_ppar$ADJ,directed=FALSE),
coord=soc_ppar$soc_coord,
vertex.cex=20,
vertex.col=NA,
xlim=soc_ppar$xlimy,
ylim=soc_ppar$ylimy,
jitter = FALSE,
pad=0,
label=soc_ppar$personnames,
	label.pos = 7,label.cex=5,label.col='grey')
}
if (do.small==FALSE){
	plot(as.network(soc_ppar$ADJ,directed=FALSE),
coord=soc_ppar$soc_coord,
vertex.cex=.5,
vertex.col=NA,
xlim=soc_ppar$xlimy,
ylim=soc_ppar$ylimy,
jitter = FALSE,
pad=0,
label=soc_ppar$personnames,
	label.pos = 7,label.cex=.25,label.col='grey')
}

	#browser()
	for (k in c(1:n))
	{
		
		Acoord <- gplot.layout.fruchtermanreingold(g[k,,],NULL)
		Acoord <- scale.coord(Acoord,soc_ppar$soc_coord[k,])
		#print(Acoord)
		if (do.small){
		plot(as.network(g[k,,],directed=FALSE),
		label=concepts,new=FALSE,coord=Acoord,
		xlim=soc_ppar$xlimy,ylim=soc_ppar$ylimy,
		jitter = FALSE,pad=0,
		label.pos = 7, boxed.labels=TRUE
		)#, label.pad =-1)
		}
		if (do.small==FALSE){
			plot(as.network(g[k,,],directed=FALSE),
		new=FALSE,coord=Acoord,vertex.cex=.25,vertex.col='grey',
		xlim=soc_ppar$xlimy,ylim=soc_ppar$ylimy,
		jitter = FALSE,pad=0,
		label.pos = 7		)#, label.pad =-1)
			}

		
		
	}
	
}

scale.coord <- function(Acoord,scale.macro)
{
	Acoord[,1]<- Acoord[,1] -min(Acoord[,1])
		Acoord[,2]<- Acoord[,2] -min(Acoord[,2])
		Acoord[,1]<- Acoord[,1]/max(Acoord[,1])
		Acoord[,2]<- Acoord[,2]/max(Acoord[,2])
		Acoord <- Acoord*8
		Acoord[,1]<-Acoord[,1]+(scale.macro[1]-4)
		Acoord[,2]<-Acoord[,2]+(scale.macro[2]-4)
		Acoord
}

socnet_lay_out <- function(Y,respnames)
{
	n <- length(respnames)
	if (n==2)
	{
		xcoord <-c(5,25)
		ycoord <- c(5,5)
		xlimy <- c(0,30)
		ylimy <- c(0,10) 
	}
	if (n==3)
	{
		xcoord <-c(5,15,25)
		ycoord <- c(5,25,5)
		xlimy <- c(0,30)
		ylimy <- c(0,30) 
	}
	if (n==4)
	{
		xcoord <-c(5,5,25,25)
		ycoord <- c(5,25,5,25)
		xlimy <- c(0,30)
		ylimy <- c(0,30) 
	}
	if (n >4)
	{
		require('network')
		coord <- plot(as.network(Y,directed=FALSE))
		xcoord <- coord[,1]
		ycoord <- coord[,2]
		xlimy <- range(xcoord)
		ylimy <- range(ycoord)
	}
	soc_coord <- matrix(0,n,2)
	soc_coord[,1] <- xcoord
	soc_coord[,2] <- ycoord
	soc_ppar <- list(soc_coord=soc_coord,xlimy=xlimy,ylimy,ADJ=Y,personnames=respnames)
	soc_ppar
}

get.union.net <- function(g,binary=TRUE)
{
	N <- dim(g)[3]
	n <- dim(g)[1]
	UniNet <- matrix(0,N,N)
	for (k in c(1:n))
	{
		UniNet <- UniNet+g[k,,]
	}
	if (binary){
	UniNet[UniNet>0] <- 1
	}
	UniNet
}

get.line.graph.edgelist <- function(UniNet,thresh=0,expert=NULL)
{
	upperUni <- UniNet
	upperUni[lower.tri(upperUni)] <- 0
	if (is.null(expert)){
	Eset <- which(upperUni>thresh, arr.ind = TRUE)
	}
	if (is.null(expert)==FALSE){
	upperexpert <- expert
	upperexpert[lower.tri(upperexpert)]	 <- 0
	Eset <- which(upperUni>thresh | upperexpert>0, arr.ind = TRUE)
	}
	Eset
}

get.kneser <- function(Eset,concepts)
{
	vert.set.comp <- unique(c(Eset))
	small.e <- dim(Eset)[1]
	Q <- matrix(0,small.e,small.e)
	node.name <- matrix(NA,small.e,1)
	for (i in c(1:(small.e-1)) )
	{
		node.name[i] <- paste(concepts[Eset[i,1]],'-',concepts[Eset[i,2]] )
		for (j in c((i+1):small.e))
		{
			# check if edges share node
			conc.nodes <- unique( c(Eset[i,],Eset[j,]) )
			if (length(conc.nodes)<4)
			{
				Q[i,j] <- 1
				Q[j,i] <- Q[i,j]
			}
		}
	}
	node.name[small.e] <- paste(concepts[Eset[small.e,1]],'-',concepts[Eset[small.e,2]] )
	row.names(Q) <- node.name
	Q
}

combine.in.block <- function(g,Q,Eset,concepts,Y,respnames,expert=NULL){
	N <- dim(g)[3]
	n <- dim(g)[1]
	vert.set.comp <- unique(c(Eset))
	small.e <- dim(Eset)[1]
	node.name <- row.names(Q)
	W <- matrix(0,n,small.e)
	for (i in c(1:n))
	{
		for (j in c(1:small.e))
		{
			sender <- Eset[j,1]
			receiver <- Eset[j,2]
			if (g[i,sender,receiver]==1)
			{
				W[i,j] <- 1
			}
		}
	}
	A <- matrix(NA,n+small.e,1)
	if (is.null(expert)==FALSE)
	{
		for (j in c(1:small.e))
		{
			sender <- Eset[j,1]
			receiver <- Eset[j,2]
			if (expert[sender,receiver]==1)
			{
				A[j] <- 'grey'
			}
			if (expert[sender,receiver]==0)
			{
				A[j] <- 'white'
			}
		}
		
	}
	C <- cbind(Q,t(W))
	C <- rbind(C,cbind(W,Y))
	row.names(C) <- A
	C
}

plot.expert <- function(expert,concepts)
{
	plot( as.network(expert,directed=FALSE),
	 boxed.labels=TRUE,
	 label=concepts,
	 label.pos = 5)
}

combine.plots <- function(){
	pdf('example_NoN.pdf',height=9,width=11)
	layout.matrix <- matrix(c(0,2,1, 2, 0,2,3, 3), nrow = 2, ncol = 4)

layout(mat = layout.matrix,
       heights = c(1, 3), # Heights of the two rows
       widths = c(1,2,1, 5)) # Widths of the two columns

#layout.show(3)
plot.expert(expert,concepts)
title(main='(a) Criterion graph')
mult.net.plot(g,concepts,soc_ppar)
title(main='(b) Network of networks')
# Multilevel graph
is.people <- matrix(FALSE,dim(C)[1],1)
is.people[(dim(Q)[1]+1):dim(C)[1]] <- TRUE
use.box<- matrix('grey',dim(C)[1],1)
use.box[is.people] <- NA
#lab.bg<- matrix('grey',dim(Q)[1],1)
#lab.bg[is.people] <- NA
kneser.col <- matrix(NA,dim(Q)[1],dim(Q)[2])
kneser.col[Q==1] <- 'red'
edge.colour <- matrix(NA,dim(C)[1],dim(C)[2])
edge.colour[C==1] <- 'black'
edge.colour[1:dim(Q)[1],1:dim(Q)[1]]<- kneser.col

C.net <- as.network(C,directed=FALSE)
#C.net %e% "myeval" <- edge.colour

lab.bg<- row.names(C)
lab.size <- matrix(1,dim(C)[1],1)
lab.size[is.people] <- 6
lab.col <- matrix('black',dim(C)[1],1)
lab.col[is.people] <- 'grey'
vert.size <- matrix(1,dim(C)[1],1)
vert.size [is.people] <- 7
plot( C.net, 
label=c(row.names(Q),respnames),
label.pos = 7,
boxed.labels=TRUE,
vertex.col=NA,
label.pad=1,
label.border=use.box,
label.bg = lab.bg,
label.col = lab.col,
vertex.cex=vert.size,
edge.col=edge.colour)
title(main='(c) Multielvel representation of NoN')
dev.off()
}

dont_run <- function(){
concepts <- c('risk','flood','management','barrier','pump','monitor')
g <- rand_cog_graphs(concepts,prob=.25,n=4)
Y <- matrix(0,4,4)
Y[1,2] <- 1
Y[2,3] <- 1
Y[3,4] <- 1
Y <- Y+t(Y)
respnames <- c('peter','sarah','ben','kim')
soc_ppar <- socnet_lay_out(Y,respnames)
mult.net.plot(g,concepts,soc_ppar)
UniNet <- get.union.net(g)
expert <- rgraph(length(concepts),1,tprob=.25,mode='graph')
Eset <- get.line.graph.edgelist(UniNet,expert=expert)
Q <- get.kneser(Eset,concepts)

C <- combine.in.block(g,Q,Eset,concepts,Y,respnames,expert)

plot.expert(expert,concepts)

is.people <- matrix(FALSE,dim(C)[1],1)
is.people[(dim(Q)[1]+1):dim(C)[1]] <- TRUE
use.box<- matrix('grey',dim(C)[1],1)
use.box[is.people] <- NA
#lab.bg<- matrix('grey',dim(Q)[1],1)
#lab.bg[is.people] <- NA
kneser.col <- matrix(NA,dim(Q)[1],dim(Q)[2])
kneser.col[Q==1] <- 'red'
edge.colour <- matrix(NA,dim(C)[1],dim(C)[2])
edge.colour[C==1] <- 'black'
edge.colour[1:dim(Q)[1],1:dim(Q)[1]]<- kneser.col

C.net <- as.network(C,directed=FALSE)
#C.net %e% "myeval" <- edge.colour

lab.bg<- row.names(C)
lab.size <- matrix(1,dim(C)[1],1)
lab.size[is.people] <- 6
lab.col <- matrix('black',dim(C)[1],1)
lab.col[is.people] <- 'grey'
vert.size <- matrix(1,dim(C)[1],1)
vert.size [is.people] <- 7
plot( C.net, 
label=c(row.names(Q),respnames),
label.pos = 7,
boxed.labels=TRUE,
vertex.col=NA,
label.pad=1,
label.border=use.box,
label.bg = lab.bg,
label.col = lab.col,
vertex.cex=vert.size,
edge.col=edge.colour)

use.box<- matrix('red',dim(Q)[1],1)
use.box[1:2] <- NA
lab.cox<- matrix('grey',dim(Q)[1],1)
lab.cox[1:2] <- NA
plot(as.network(Q,directed=FALSE),label=row.names(Q),label.pos = 7,boxed.labels=TRUE,vertex.co=NA,label.pad=1,
label.border=use.box,
label.bg = lab.cox,
edge.col=edge.colour)
}
