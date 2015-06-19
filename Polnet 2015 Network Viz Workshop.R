
##======================================================##
##                                                      ##
##    POLNET 2015 Network Visualization Workshop        ##
##    Porland OR, June 18, 2015                         ##
##                                                      ##
##    Katya Ognyanova, katya@ognyanova.net              ##
##    www.kateto.net/polnet2015                         ##
##                                                      ##
##======================================================##


# Download handouts and example data from www.kateto.net/polnet2015


# Key packages to install: 

install.packages("igraph") 
install.packages("network") 
install.packages("ndtv")

# Set the working directory to the folder containing the example data
setwd("C:/Data") 


# ================ Read the example data ================


# DATASET 1: edgelist 

nodes <- read.csv("Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)

# Examine the data:
head(nodes)
head(links)
nrow(nodes); length(unique(nodes$id))
nrow(links); nrow(unique(links[,c("from", "to")]))

# Collapse multiple links of the same type between the same two nodes
# by summing their weights, using aggregate() by "from", "to", & "type":
links <- aggregate(links[,3], links[,-3], sum)
links <- links[order(links$from, links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL


# DATASET 2: matrix 

nodes2 <- read.csv("Dataset2-Media-User-Example-NODES.csv", header=T, as.is=T)
links2 <- read.csv("Dataset2-Media-User-Example-EDGES.csv", header=T, row.names=1)

# Examine the data:
head(nodes2)
head(links2)

# links2 is an adjacency matrix for a two-mode network:
links2 <- as.matrix(links2)
dim(links2)
dim(nodes2)


# ================ Plotting networks with igraph ================
 
library(igraph)

# Converting the data to an igraph object:
net <- graph.data.frame(links, nodes, directed=T) 

# Examine the resulting object:
class(net)
net 

# It's easy to access nodes, edges, and their attributes:
E(net)
V(net)
E(net)$type
V(net)$media

# You can also manipulate the network matrix:
net[1,]
net[5,7]

# First attempt to plot the graph:
plot(net) # not pretty!

# Removing loops from the graph:
net <- simplify(net, remove.multiple = F, remove.loops = T) 

# Let's and reduce the arrow size and remove the labels:
plot(net, edge.arrow.size=.4,vertex.label=NA)


# ================ A brief detour on colors in R plots ================

# In most R functions, you can use named colors, hex, or rgb values:
# (In the simple base plot chart below x and y are point coordiantes, pch 
# is the point symbol shape, cex is the point size, and col is the color.
# to see the parameters for ploting in base R, check out ?par
plot(x=1:10, y=rep(5,10), pch=19, cex=5, col="dark red")
points(x=1:10, y=rep(6, 10), pch=19, cex=5, col="#557799")
points(x=1:10, y=rep(4, 10), pch=19, cex=5, col=rgb(.25, .5, .3))

# You may notice that rgb here ranges from 0 to 1. While this is the R default,
# you can also set it for the 0-255 range: 
rgb(10, 100, 100, maxColorValue=255) 

# We can also set the opacity/transparency using the parameter 'alpha' (range 0-1):
plot(x=1:5, y=rep(5,5), pch=19, cex=16, col=rgb(.25, .5, .3, alpha=.5), xlim=c(0,6))  

# If we have a hex color representation, we can set the transparency alpha 
# using 'adjustcolor' from package 'grDevices'. For fun, let's also set the
# the plot background to gray using the par() function for graphical parameters.
par(bg="black")
col.tr <- grDevices::adjustcolor("#557799", alpha=0.7)
plot(x=1:5, y=rep(5,5), pch=19, cex=20, col=col.tr, xlim=c(0,6)) 
par(bg="white")

# If you plan on using the built-in color names, here's what they are: 
colors()
grep("blue", colors(), value=T)

# In many cases, we need a number of contrasting colors, or multiple shades of a color.
# R comes with some predefined palette function that can generate those for us.
pal1 <- heat.colors(5, alpha=1)   # generate 5 colors from the heat palette, opaque
pal2 <- rainbow(5, alpha=.5)      # generate 5 colors from the heat palette, semi-transparent
plot(x=1:10, y=1:10, pch=19, cex=10, col=pal1)
plot(x=10:1, y=1:10, pch=19, cex=10, col=pal2)

# We can also generate our own gradients using colorRampPalette.
# Note that colorRampPalette returns a *function* that we can use 
# to generate as many colors from that palette as we need.

palf <- colorRampPalette(c("gray80", "dark red")) 
plot(x=10:1, y=1:10, pch=19, cex=10, col=palf(100)) 

# To add transparency to colorRampPalette, you need to add a parameter `alpha=TRUE`:
palf <- colorRampPalette(c(rgb(1,1,1, .2),rgb(.8,0,0, .7)), alpha=TRUE)
plot(x=10:1, y=1:10, pch=19, cex=10, col=palf(10)) 

# Finding good color combinations is a tough task - and the built-in R palettes
# are rather limited. Thankfully there are other available packages for this:

# install.packages("RColorBrewer")
library(RColorBrewer)

display.brewer.all()

# This package has one main function, called 'brewer.pal'.
# Using it, you just need to select the desired palette and a number of colors.
# Let's take a look at some of the RColorBrewer palettes:
display.brewer.pal(8, "Set3")
display.brewer.pal(8, "Spectral")
display.brewer.pal(8, "Blues")

pal3 <- brewer.pal(10, "Set3")
plot(x=10:1, y=10:1, pch=19, cex=6, col=pal3)
plot(x=10:1, y=10:1, pch=19, cex=6, col=rev(pal3)) # backwards

detach(package:RColorBrewer)

# ================ A brief detour on fonts in R plots ================

# Using different fonts for R plots may take a little bit of work.
# Especially for Windows - Mac & Linux users may not have to do this.
# First we'd use the 'extrafont' package to import the fonts from the OS into R:

# install.packages("extrafont")

library(extrafont)

# Import system fonts - may take a while, so DO NOT run this during the workshop.
# font_import() 
fonts() # See what font families are available to you now.
loadfonts(device = "win") # use device = "pdf" for pdf plot output. 

plot(net, vertex.size=30)

# Now you should be able to do  this:
plot(net, vertex.size=30, vertex.label.family="Arial Black" )

# To embed the fonts & use them in PDF files:
# First you may have to let R know where to find ghostscript
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.10/bin/gswin64c.exe")

# The command 'pdf' will send all the plots we output before dev.off() to a pdf file: 
pdf(file="ArialBlack.pdf")
plot(net, vertex.size=30, vertex.label.family="Arial Black" )
dev.off()

embed_fonts("ArialBlack.pdf", outfile="ArialBlack_embed.pdf")


detach(package:extrafont)


# ================ Back to network plots with igraph ================


#  ------->> Plot parameters in igraph --------

# Plotting with igraph: node options (starting with 'vertex.) and edge options
# (starting with 'edge.'). A list of options is included in your handout.
?igraph.plotting

# We can set the node & edge options in two ways - one is to specify
# them in the plot() function, as we are doing below.

# Plot with curved edges (edge.curved=.1) and reduce arrow size:
plot(net, edge.arrow.size=.4, edge.curved=.1)

# Set node color to orange and the border color to hex #555555
# Replace the vertex label with the node names stored in "media"
plot(net, edge.arrow.size=.4, edge.curved=0,
     vertex.color="orange", vertex.frame.color="#555555",
     vertex.label=V(net)$media, vertex.label.color="black",
     vertex.label.cex=.7) 


# The second way to set attributes is to add them to the igraph object.

# Generate colors base on media type:
colrs <- c("gray50", "tomato", "gold")
V(net)$color <- colrs[V(net)$media.type]

# Compute node degree (#links) and use it to set node size:
deg <- degree(net, mode="all")
V(net)$size <- deg*3
V(net)$size <- V(net)$audience.size*0.6

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net)$label.color <- "black"
V(net)$label <- NA

# Set edge width based on weight:
E(net)$width <- E(net)$weight/6

#change arrow size and edge color:
E(net)$arrow.size <- .2
E(net)$edge.color <- "gray80"

plot(net) 

# We can also override the attributes explicitly in the plot:
plot(net, edge.color="orange", vertex.color="gray50") 


# We can also add a legend explaining the meaning of the colors we used:
plot(net) 
legend(x=-1.1, y=-1.1, c("Newspaper","Television", "Online News"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2.5, bty="n", ncol=1)


# Sometimes, especially with semantic networks, we may be interested in 
# plotting only the labels of the nodes:

plot(net, vertex.shape="none", vertex.label=V(net)$media, 
     vertex.label.font=2, vertex.label.color="gray40",
     vertex.label.cex=1.2, edge.color="gray90")


# Let's color the edges of the graph based on their source node color.
# We'll get the starting node for each edge with "get.edges"
edge.start <- get.edges(net, 1:ecount(net))[,1]
edge.col <- V(net)$color[edge.start]

plot(net, edge.color=edge.col, edge.curved=.1)

#  ------->> Network Layouts --------

# Network layouts are algorithms that return coordinates for each
# node in a network.

# Let's generate a slightly larger 80-node graph.

net.bg <- barabasi.game(80) 
V(net.bg)$size <- 8
V(net.bg)$color <- "orange"
V(net.bg)$label <- "" 
E(net.bg)$arrow.mode <- 0
plot(net.bg)

# You can set the layout in the plot function:
plot(net.bg, layout=layout.random)

# Or calculate the vertex coordinates in advance:
l <- layout.circle(net.bg)
plot(net.bg, layout=l)

# l is simply a matrix of x,y coordinates (N x 2) for the N 
# nodes in the graph. You can generate your own:
l
l <- matrix(c(1:vcount(net.bg), c(1, vcount(net.bg):2)), vcount(net.bg), 2)
plot(net.bg, layout=l)

# This layout is just an example and not very helpful - thankfully
# igraph has a number of built-in layouts, including:

# Randomly placed vertices
l <- layout.random(net.bg)
plot(net.bg, layout=l)

# Circle layout
l <- layout.circle(net.bg)
plot(net.bg, layout=l)

# 3D sphere layout
l <- layout.sphere(net.bg)
plot(net.bg, layout=l)

# The Fruchterman-Reingold force-directed algorithm 
# Nice but slow, most often used in graphs smaller than ~1000 vertices. 

# Some parameters you can set are the area (default is the square of # nodes)
# and repulserad (cancelation radius for the repulsion - the area multiplied by # nodes). 
# Both parameters affect the spacing  of the plot - play with them until you like the results.
# The "weight" parameter increases the attraction among nodes connected by heavier edges.

l <- layout.fruchterman.reingold(net.bg, repulserad=vcount(net.bg)^3, 
                                      area=vcount(net.bg)^2.4)
par(mfrow=c(1,2)) # plot two figures - 1 row, 2 columns
plot(net.bg, layout=layout.fruchterman.reingold)
plot(net.bg, layout=l)
dev.off() # shut off the  graphic device to clear the two-figure configuration.

# You will also notice that the layout is not deterministic - different runs 
# will result in slightly different configurations. Saving the layout in l
# allows us to get the exact same result multiple times.
par(mfrow=c(2,2), mar=c(1,1,1,1))
plot(net.bg, layout=layout.fruchterman.reingold)
plot(net.bg, layout=layout.fruchterman.reingold)
plot(net.bg, layout=l)
plot(net.bg, layout=l)
dev.off()

# By default, the coordinates of the plots are rescaled to the [-1,1] interval
# for both x and y. You can change that with the parameter "rescale=FALSE"
# and rescale your plot manually by multiplying the coordinates by a scalar.
# You can use layout.norm to normalize the plot with the boundaries you want.

# Get the layout coordinates:
l <- layout.fruchterman.reingold(net.bg)
# Normalize them so that they are in the -1, 1 interval:
l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

par(mfrow=c(2,2))
plot(net.bg, rescale=F, layout=l*0.4)
plot(net.bg, rescale=F, layout=l*0.8)
plot(net.bg, rescale=F, layout=l*1.2)
plot(net.bg, rescale=F, layout=l*1.6)
dev.off()

# Another popular force-directed algorithm that produces nice results for
# connected graphs is Kamada Kawai. Like Fruchterman Reingold, it attempts to 
# minimize the energy in a spring system.
# Igraph also has a spring-embedded layout called layout.spring().

l <- layout.kamada.kawai(net.bg)
plot(net.bg, layout=l)

l <- layout.spring(net.bg, mass=.5)
plot(net.bg, layout=l)


# LGL algorithm for large connected graphs. Here you can specify a root - 
# the node that will be placed in the middle of the layout.
plot(net.bg, layout=layout.lgl)

# By default, igraph uses a layout called layout.auto which selects
# an appropriate layout algorithm based on the properties of the graph. 

# Check out all available layouts in igraph:
?igraph::layout

layouts <- grep("^layout\\.", ls("package:igraph"), value=TRUE) 
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama", layouts)]

par(mfrow=c(3,3), mar=c(1,1,1,1))

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(net)) 
  plot(net, edge.arrow.mode=0, layout=l, main=layout) }

dev.off()



# ------->> Highlighting aspects of the network --------

plot(net)

# Notice that this network plot is still not too helpful.
# We can identify the type and size of nodes, but cannot see
# much about the structure since the links we're examining are so dense.
# One way to approach this is to see if we can sparsify the network.

hist(links$weight)
mean(links$weight)
sd(links$weight)

# There are more sophisticated ways to extract the key edges,
# but for the purposes of this excercise we'll only keep ones
# that have weight higher than the mean for the network.

# We can delete edges using delete.edges(net, edges)
cut.off <- mean(links$weight) 
net.sp <- delete.edges(net, E(net)[weight<cut.off])
plot(net.sp) 

l <- layout.fruchterman.reingold(net.sp, repulserad=vcount(net.sp)^2.1)
plot(net.sp, layout=l) 


# Another way to think about this is to plot the two tie types 
# (hyperlik & mention) separately:

E(net)$width <- 2
plot(net, edge.color=c("dark red", "slategrey")[(E(net)$type=="hyperlink")+1],
      vertex.color="gray40", layout=layout.circle)

# Another way to delete edges:  
net.m <- net - E(net)[E(net)$type=="hyperlink"]
net.h <- net - E(net)[E(net)$type=="mention"]

# Plot the two links separately:
par(mfrow=c(1,2))

plot(net.h, vertex.color="orange", main="Tie: Hyperlink")
plot(net.m, vertex.color="lightsteelblue2", main="Tie: Mention")

dev.off()

# Make sure the nodes stay in place in both plots:
par(mfrow=c(1,2),mar=c(1,1,4,1))

l <- layout.fruchterman.reingold(net)
plot(net.h, vertex.color="orange", layout=l, main="Tie: Hyperlink")
plot(net.m, vertex.color="lightsteelblue2", layout=l, main="Tie: Mention")

dev.off()


# We can also try to make the network map more useful by
# showing the communities within it:
optimal.community(net)
V(net)$community <- optimal.community(net)$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(net, vertex.color=colrs[V(net)$community])




# ------->> Highlighting specific nodes or links --------


# Sometimes we want to focus the visualization on a particular node
# or a group of nodes. Let's represent distance from the NYT:
shortest.paths(net, algorithm="unweighted")
dist.from.NYT <- shortest.paths(net, algorithm="unweighted")[1,]

oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.NYT)+1)
col <- col[dist.from.NYT+1]

plot(net, vertex.color=col, vertex.label=dist.from.NYT, edge.arrow.size=.6, 
     vertex.label.color="white")


# Or we can show all the immediate neighbors of the WSJ:
col <- rep("grey40", vcount(net))
col[V(net)$media=="Wall Street Journal"] <- "#ff5100"

# The neighbors function finds all nodes one step out from the focal actor:
# (the corresponding function that finds all edges for a node is "incident")
neigh.nodes <- neighbors(net, V(net)[media=="Wall Street Journal"], mode="out")

col[neigh.nodes] <- "#ff9d00"
plot(net, vertex.color=col)


# Another way to draw attention to a group of nodes:
plot(net, mark.groups=c(1,4,5,8), mark.col="#C5E5E7", mark.border=NA)
# Mark multiple groups:
plot(net, mark.groups=list(c(1,4,5,8), c(15:17)), 
          mark.col=c("#C5E5E7","#ECD89A"), mark.border=NA)

# Highlight a path in the network
news.path <- get.shortest.paths(net, 
                                V(net)[media=="MSNBC"], 
                                V(net)[media=="New York Post"],
                                mode="all", output="both")

# Generate edge color variable:
ecol <- rep("gray80", ecount(net))
ecol[unlist(news.path$epath)] <- "orange"

# Generate edge width variable:
ew <- rep(2, ecount(net))
ew[unlist(news.path$epath)] <- 4

# Generate node color variable:
vcol <- rep("gray40", vcount(net))
vcol[unlist(news.path$vpath)] <- "gold"

plot(net, vertex.color=vcol, edge.color=ecol, 
     edge.width=ew, edge.arrow.mode=0)
      


# ------->> Interactive plotting with tkplot -------- 

# R and igraph offer interactive plotting capabilities
# (mostly helpful for small networks)

tkid <- tkplot(net) #tkid is the id of the tkplot
l <- tkplot.getcoords(tkid) # grab the coordinates from tkplot
plot(net, layout=l)



# ------->> Other ways to represent a network -------- 

# One reminder that there are other ways to represent a network:

# Heatmap of the network matrix:
netm <- get.adjacency(net, attr="weight", sparse=F)
colnames(netm) <- V(net)$media
rownames(netm) <- V(net)$media

palf <- colorRampPalette(c("gold", "dark orange")) 
heatmap(netm[,17:1], Rowv = NA, Colv = NA, col = palf(20), 
        scale="none", margins=c(10,10) )

# Degree distribution
dd <- degree.distribution(net, cumulative=T, mode="all")
plot(dd, pch=19, cex=2, col="orange", xlab="Degree", ylab="Cumulative Frequency")



# ================ Plotting two-mode networks with igraph ================

head(nodes2)
head(links2)

net2 <- graph.incidence(links2)
plot(net2)

# This time we will make nodes look different based on their type.
V(net2)$color <- c("steel blue", "orange")[V(net2)$type+1]
V(net2)$shape <- c("square", "circle")[V(net2)$type+1]
V(net2)$label <- ""
V(net2)$label[V(net2)$type==F] <- nodes2$media[V(net2)$type==F] 
V(net2)$label.cex=.6
V(net2)$label.font=2

plot(net2, vertex.label.color="white", vertex.size=(2-V(net2)$type)*8) 

plot(net2, vertex.label=NA, vertex.size=7, layout=layout.bipartite) 

 
# Using text as nodes:
par(mar=c(0,0,0,0))
plot(net2, vertex.shape="none", vertex.label=nodes2$media,
     vertex.label.color=V(net2)$color, vertex.label.font=2, 
     vertex.label.cex=.95, edge.color="gray70",  edge.width=2)

dev.off()


# Using images as nodes
# You will need the 'png' library to do this:

# install.packages("png")
library("png")
 
img.1 <- readPNG("./images/news.png")
img.2 <- readPNG("./images/user.png")

V(net2)$raster <- list(img.1, img.2)[V(net2)$type+1]

par(mar=c(3,3,3,3))

plot(net2, vertex.shape="raster", vertex.label=NA,
     vertex.size=16, vertex.size2=16, edge.width=2)

# By the way, you can also add any image you want to any plot.
# For example, many network graphs could be improved by a photo
# of a puppy carrying a basket full of kittens.
img.3 <- readPNG("./images/puppy.png")
rasterImage(img.3,  xleft=-1.7, xright=0, ybottom=-1.2, ytop=0)


# The numbers after the image are coordinates for the plot.
# The limits of your plotting area are given in par()$usr

dev.off()

detach(package:png) 
detach(package:igraph)


# ================ Quick example using the network package ================

# Plotting with the 'network' package is very similar to that with igraph -
# although the notation is slightly different (a whole new set of parameter names!)
# Here is a quick example using the (by now familiar) media network.

library(network)

# Convert the data into the network format used by the Statnet family.
# As in igraph, we can generate a 'network' object from an edgelist, 
# an adjacency matrix, or an incidence matrix. 

# Remember to set the ignore.eval to F for weighted networks.
net3 <- network(links, vertex.attr=nodes, matrix.type="edgelist", 
                loops=F, multiple=F, ignore.eval = F)
net3

# You can access the edges, vertices, and the network matrix using:
net3[,]
net3 %n% "net.name" <- "Media Network" #  network attribute
net3 %v% "media"    # Node attribute
net3 %e% "type"     # Node attribute
 
net3 %v% "col" <- c("gray70", "tomato", "gold")[net3 %v% "media.type"]

# plot the network:
plot(net3, vertex.cex=(net3 %v% "audience.size")/7, vertex.col="col")

# For a full list of parameters that you can use in this plot, 
# check out ?plot.network.
?plot.network

# Note that - as in igraph - the plot returns the node position coordinates.
l <- plot(net3, vertex.cex=(net3 %v% "audience.size")/7, vertex.col="col")
plot(net3, vertex.cex=(net3 %v% "audience.size")/7, vertex.col="col", coord=l)


detach(package:network)

# ================ Interactive D3 Networks ================

# There are a number of libraries like rcharts and htmlwidgets that can
# help you export interactive web charts from R. We'll take a quick look
# at networkD3 which exports networks from r to javascript. 

# install.packages("networkD3")

library(networkD3) 
 
# d3ForceNetwork expects node IDs that are numeric and start from 0
# so we have to transform our character node IDs:

el <- data.frame(from=as.numeric(factor(links$from))-1, 
                 to=as.numeric(factor(links$to))-1 )

# The nodes need to be in the same order as the "source" column in links:
nl <- cbind(idn=factor(nodes$media, levels=nodes$media), nodes) 

# The `Group` parameter is used to color the nodes.
# Nodesize is not (as you might think) the size of the node, but the
# number of the column in the node data that should be used for sizing.
# The `charge` parameter guides node repulsion (if negative) or 
# attraction (if positive).

forceNetwork(Links = el, Nodes = nl, Source="from", Target="to",
               NodeID = "idn", Group = "type.label",linkWidth = 1,
               linkColour = "#afafaf", fontSize=12, zoom=T, legend=T,
               Nodesize=6, opacity = 1, charge=-600, 
               width = 600, height = 600)
                 

detach(package: networkD3)

# ================ Simple Plot Animations in R ================

# If you have already installed "ndtv", you should also have
# a package used by it called "animation".

# install.packages('animation')
library(animation)
library(igraph)

# In order for this to work, you need not only the R package,
# but also an additional software called ImageMagick from imagemagick.org 

ani.options("convert") # Check that the package knows where to find ImageMagick
ani.options(convert="C:/Program Files/ImageMagick-6.8.8-Q16/convert.exe") 
 
# You can use this technique to create various (not necessarily network-related)
# animations in R by generating multiple plots and combining them in an animated GIF.

l <- layout.fruchterman.reingold(net)

saveGIF( {  col <- rep("grey40", vcount(net))
            plot(net, vertex.color=col, layout=l)
            
            step.1 <- V(net)[media=="Wall Street Journal"]
            col[step.1] <- "#ff5100"
            plot(net, vertex.color=col, layout=l)
            
            step.2 <- unlist(neighborhood(net, 1, step.1, mode="out"))
            col[setdiff(step.2, step.1)] <- "#ff9d00"
            plot(net, vertex.color=col, layout=l) 
            
            step.3 <- unlist(neighborhood(net, 2, step.1, mode="out"))
            col[setdiff(step.3, step.2)] <- "#FFDD1F"
            plot(net, vertex.color=col, layout=l)  },
          interval = .8, movie.name="network_animation.gif" )

 detach(package:igraph)
 detach(package:animation)



# ================ Interactive networks with ndtv-d3 ================



# ------->> Interactive network plots --------


# install.packages("ndtv", dependencies=T)

library(ndtv)

# You should not need additional software to produce web animations with D3 (below).
# If you want to save the animations as  video  files ( see ?saveVideo), you
# would have to install a video converter called FFmpeg (http://ffmpg.org) 
# To find out how to get the right installation for your OS, check out ?install.ffmpeg
# To use all available layouts, you would need to have Java installed on your machine.

net3 


render.d3movie(net3, usearrows = F, displaylabels = F, bg="#111111", 
               vertex.border="#ffffff", vertex.col =  net3 %v% "col",
               vertex.cex = (net3 %v% "audience.size")/8, 
               edge.lwd = (net3 %e% "weight")/3, edge.col = '#55555599',
               vertex.tooltip = paste("<b>Name:</b>", (net3 %v% 'media') , "<br>",
                                    "<b>Type:</b>", (net3 %v% 'type.label')),
               edge.tooltip = paste("<b>Edge type:</b>", (net3 %e% 'type'), "<br>", 
                                  "<b>Edge weight:</b>", (net3 %e% "weight" ) ),
               launchBrowser=T, filename="Media-Network.html" )  

    
# If you are going to embed this in a markdown document, 
# you would also need to use output.mode='inline' above.


#=================================#
# Network evolution animations
#=================================#

# In order to work with the network animations in ndtv, we need to understand the 
# dynamic network format used by Statnet packages, implemented in networkDynamic

# Let's look at one of the example datasets included in the package:
data(short.stergm.sim)
short.stergm.sim 
head(as.data.frame(short.stergm.sim))

# Plot the network ignoring time (all nodes & edges that were ever present):
plot(short.stergm.sim)  

# Plot the network at time 1
plot( network.extract(short.stergm.sim, at=1) )

# Plot nodes & vertices that were active from time 1 to time 5:
plot( network.extract(short.stergm.sim, onset=1, terminus=5, rule="all") )

# Plot all nodes and vertices that were active between time 1 & 10:
plot( network.extract(short.stergm.sim, onset=1, terminus=10, rule="any") ) 

# Let's make a quick d3 animation from the example network:
render.d3movie(short.stergm.sim,displaylabels=TRUE) 

# We are now ready to create and animate our own dynamic network.

# Dynamic network object can be generated in a number of ways: from 
# a set of networks/matrices representing different time points, or from
# data frames/matrices with node lists and edge lists indicating when each
# is active, or when they switch state. See ?networkDynamic for more information.

net3
plot(net3)

vs <- data.frame(onset=0, terminus=50, vertex.id=1:17)
es <- data.frame(onset=1:49, terminus=50, 
                 head=as.matrix(net3, matrix.type="edgelist")[,1],
                 tail=as.matrix(net3, matrix.type="edgelist")[,2])

net3.dyn <- networkDynamic(base.net=net3, edge.spells=es, vertex.spells=vs)


plot(net3.dyn, vertex.cex=(net3 %v% "audience.size")/7, vertex.col="col")


# We can pre-compute the animation coordinates (otherwise they get calculated when 
# you generate the animation). Here animation.mode is the layout algorithm - 
# one of "kamadakawai", "MDSJ", "Graphviz"and "useAttribute" (user-generated).

compute.animation(net3.dyn, animation.mode = "kamadakawai",
                  slice.par=list(start=0, end=49, interval=10, 
                         aggregate.dur=10, rule='any'))

# Show time evolution through static images at different time points:
 filmstrip(net3.dyn, displaylabels=F, mfrow=c(2, 3),
           slice.par=list(start=0, end=49, interval=10, 
                         aggregate.dur=10, rule='any'))

# Let's make an actial animation: 
 
compute.animation(net3.dyn, animation.mode = "kamadakawai",
                  slice.par=list(start=0, end=50, interval=1, 
                         aggregate.dur=1, rule='any'))


render.d3movie(net3.dyn, usearrows = F, 
               displaylabels = F, label=net3 %v% "media",
               bg="#ffffff", vertex.border="#333333",
               vertex.cex = degree(net3)/2,  
               vertex.col = net3.dyn %v% "col",
               edge.lwd = (net3.dyn %e% "weight")/3, 
               edge.col = '#55555599',
               vertex.tooltip = paste("<b>Name:</b>", (net3.dyn %v% "media") , "<br>",
                                    "<b>Type:</b>", (net3.dyn %v% "type.label")),
               edge.tooltip = paste("<b>Edge type:</b>", (net3.dyn %e% "type"), "<br>", 
                                  "<b>Edge weight:</b>", (net3.dyn %e% "weight" ) ),
               launchBrowser=T, filename="Media-Network-Dynamic.html",
               render.par=list(tween.frames = 30, show.time = F),
               plot.par=list(mar=c(0,0,0,0)) )


# In addition to dynamic nodes and edges, ndtv takes dynamic attributes.
# We could have added those to the es and vs data frames above.
# In addition, the plotting function can evaluate special parameters
# and generate dynamic arguments on the fly. For example,
# function(slice) { do some calculations with slice } will perform operations
# on the current time slice, allowing us to change parameters dynamically.
# See the node size below:

compute.animation(net3.dyn, animation.mode = "kamadakawai",
                  slice.par=list(start=0, end=50, interval=4, 
                         aggregate.dur=1, rule='any'))


render.d3movie(net3.dyn, usearrows = F, 
               displaylabels = F, label=net3 %v% "media",
               bg="#000000", vertex.border="#dddddd",
               vertex.cex = function(slice){ degree(slice)/2.5 },  
               vertex.col = net3.dyn %v% "col",
               edge.lwd = (net3.dyn %e% "weight")/3, 
               edge.col = '#55555599',
               vertex.tooltip = paste("<b>Name:</b>", (net3.dyn %v% "media") , "<br>",
                                    "<b>Type:</b>", (net3.dyn %v% "type.label")),
               edge.tooltip = paste("<b>Edge type:</b>", (net3.dyn %e% "type"), "<br>", 
                                  "<b>Edge weight:</b>", (net3.dyn %e% "weight" ) ),
               launchBrowser=T, filename="Media-Network-even-more-Dynamic.html",
               render.par=list(tween.frames = 25, show.time = F) )





# ================ || ================
 





