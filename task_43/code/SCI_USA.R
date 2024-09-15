
library(igraph) #basic Netowrk-handling/ graphical representation packages
library(ggplot2)
library(RColorBrewer)
library(rgl)

library(usmap)   # Contains FIPS codes for U.S. counties and geographical data
library(dplyr)   # For data manipulation
library(tigris)   #For retrieving counties' boundaries
library(sf)      #For handling geographic manipulations


#Our aim is to build a network for Facebook "friendships" between USA's counties


plot_usmap(regions = "counties") +
  labs(title = "US Counties") + 
  theme(panel.background = element_rect(color = "black", fill = "lightblue"))

#Network building

# Define a function to convert decimal degrees to DMS notation
decimal_to_dms <- function(decimal) {
  # Determine the direction based on the hemisphere
  direction <- ifelse(decimal >= 0, "N", "S")
  if (abs(decimal) > 180) {
    direction <- ifelse(decimal >= 0, "E", "W")
  }
  
  # Convert decimal degrees to absolute value for processing
  abs_decimal <- abs(decimal)
  
  # Extract degrees
  degrees <- floor(abs_decimal)
  
  # Extract minutes
  minutes <- floor((abs_decimal - degrees) * 60)
  
  # Extract seconds
  seconds <- ((abs_decimal - degrees) * 60 - minutes) * 60
  
  # Return DMS format as a string
  sprintf("%d° %d' %.2f\" %s", degrees, minutes, seconds, direction)
}



# Vectorized function to convert a vector of DMS notations to decimal degrees
dms_to_decimal <- function(dms_vector) {
  # Use sapply to apply the conversion function to each element of the vector
  sapply(dms_vector, function(dms) {
    # Scalar conversion (same as the scalar function)
    dms <- trimws(dms)
    direction <- substr(dms, nchar(dms), nchar(dms))
    dms <- gsub(paste0(" ", direction), "", dms)
    dms <- gsub('[°\'"]', '', dms)
    parts <- unlist(strsplit(dms, "[ ]+"))
    degrees <- as.numeric(parts[1])
    minutes <- as.numeric(parts[2])
    seconds <- as.numeric(parts[3])
    decimal <- degrees + minutes / 60 + seconds / 3600
    
    if (direction %in% c("S", "W")) {
      decimal <- -decimal
    }
    
    return(as.numeric(decimal))
  })
}



vec2pal <- function(v_centr, mypal){
     val_idxs <- 1 + floor((length(mypal)-1)*normalize_vec(v_centr))
     return(mypal[val_idxs])
}

normalize_vec <- function(v_centr){
     return((v_centr - min(v_centr))/( max(v_centr) - min(v_centr) ))
}

mypal <- colorRampPalette(rev(c(brewer.pal(9, "YlGnBu"))))(20)


#We now want to generate a database for labels which links each FIPS code to its corresponding county. Moreover it shows the corresponding latitude and longitude. We use the tigris package:


#Downloading counties shapefile
# Setting year=2018 to agree with the sci counties-counties dataset
counties <- counties(year=2018,cb = TRUE)


#Add the corresponding state as a new column


#dataset with fips codes and corresponding names
data("fips_codes")
# Ensure the fips_codes dataset only contains relevant columns
fips_codes <- fips_codes %>%
  select(state_code, state_name) %>%
  distinct()  # Make sure each state_code is unique

fips_codes$state_name[fips_codes$state_code==78] <- "78" # (U.S. Virgin Islands) this procedure will avoid conflicts with the tigris package later

counties <- counties %>% rename(state_code=STATEFP) %>%
  left_join(fips_codes, by = "state_code")



#Convert to sf object and transform to WGS84 (lat/lon in degrees)
county_sf <- st_as_sf(counties, coords = c("x", "y"), crs = st_crs(counties))  # Projected CRS (Albers)
county_sf <- st_transform(county_sf, crs = 4326)  # Transform to WGS84 (degrees)

county_nodes <- st_centroid(county_sf) %>% mutate(lat = st_coordinates(geometry)[, 2],
         lon = st_coordinates(geometry)[, 1]) %>% st_drop_geometry() %>%
  select(nodeID=GEOID,nodeLabel=NAME,latitude= lat,longitude= lon,state=state_name)
#  select(nodeID=GEOID,nodeLabel=NAME,latitude= lat,longitude= lon,state=STATE_NAME)

county_nodes <- county_nodes %>% mutate(latitude = sapply(latitude, decimal_to_dms),longitude=sapply(longitude, decimal_to_dms),nodeID=as.character(nodeID))


#We want to construct a node file for each state


nodes_by_state <- split(county_nodes, county_nodes$state)
nodes_by_state <- lapply(nodes_by_state,function(df){df <- df %>% select(-state)})


#Now we want to build the edges file


sci <- read.table("county_county.tsv", header=T) #social connectdness index dataset
sci <- sci %>% rename(nodeID_from=user_loc,nodeID_to=fr_loc) %>%
  mutate(
    nodeID_from = sprintf("%05d", as.numeric(nodeID_from)),  # Add leading zeros (in this way we have the same format for all the elements)
    nodeID_to = sprintf("%05d", as.numeric(nodeID_to))       # Add leading zeros
  )



# Extract the state code from the FIPS codes (first two digits)
sci$source_state <- substr(as.character(sci$nodeID_from), 1, 2)  # First two digits of the source FIPS
sci$target_state <- substr(as.character(sci$nodeID_to), 1, 2)  # First two digits of the target FIPS

# Filter out edges that connect counties from different states
same_state_edges <- sci[sci$source_state == sci$target_state, ]
#Clean columns
same_state_edges <- subset(same_state_edges,select=-target_state)
same_state_edges <- same_state_edges %>% rename(state_code=source_state) %>% mutate(state_code = as.character(state_code))

same_state_edges <- same_state_edges %>%
  left_join(fips_codes, by = "state_code")


same_state_edges <- same_state_edges %>% subset(select=-state_code) %>% rename(state=state_name,weight=scaled_sci)

#Delete self-loops and multiple edges in the file
# Remove self-loops
same_state_edges <- same_state_edges %>% 
  filter(nodeID_from != nodeID_to)

# Remove duplicate edges (ignoring direction)
same_state_edges <- same_state_edges %>%
  rowwise() %>%
  mutate(from = min(nodeID_from, nodeID_to), 
         to = max(nodeID_from, nodeID_to)) %>%
  ungroup() %>%
  distinct(from, to, .keep_all = TRUE) %>% select(index=-c(from,to))

# Split the data by state
edges_by_state <- split(same_state_edges, same_state_edges$state)
states <- names(edges_by_state) #state names



#don't consider counties with only one node
for (i in length(states):1) {
  if(dim(nodes_by_state[[i]])[1] <= 1) {
    states <- states[-i] 
  }
}

#correct label
state_label <- ifelse(states==78,"U.S. Virgin Islands",states)


#Example: California


state <- "California"

nodes <- nodes_by_state[[state]]
edges <- edges_by_state[[state]]

g <- graph.data.frame(edges, directed=F, vertices=as.character(nodes$nodeID))
if(all(E(simplify(g))==E(g))) cat(sprintf("The graph is simple"))

layout <- matrix(0, vcount(g), 2)
layout[,1] <- dms_to_decimal(nodes$longitude)
layout[,2] <- dms_to_decimal(nodes$latitude)
#temporary
names(layout[,1]) <- NULL
names(layout[,2]) <- NULL
#

#v_centr <- eigen_centrality(g, directed = FALSE, weights = E(g)$weight)$vector
#v_centr <- betweenness(g,v=V(g),normalized = F,weights=E(g)$weights)
#v_centr <- harmonic_centrality(g,normalized = F,weights=E(g)$weights)
v_centr <- strength(g, vids = V(g), weights = E(g)$weight) 
names(v_centr) <- NULL
#

colors <- vec2pal(v_centr, mypal)
sizes <- v_centr/max(v_centr)

plot(g, layout=layout, 
     vertex.label=nodes$nodeLabel, 
     vertex.color=colors,
     vertex.frame.color=NA,
     vertex.size=sizes,
     edge.color="gray80",
     edge.width=0.5)



# Prepare the nodes and edges as data.frames
nodes <- nodes %>% rename(lon=longitude,lat=latitude) %>% mutate(lon=layout[,1],lat=layout[,2])
edges <- edges %>% rename(from=nodeID_from,to=nodeID_to)

# Merge the nodes coordinates with edges data to get x and y for each endpoint
edges <- edges %>%
  left_join(nodes, by = c("from" = "nodeID")) %>%
  rename(x = lon, y = lat) %>%
  left_join(nodes, by = c("to" = "nodeID")) %>%
  rename(xend = lon, yend = lat)



png(filename = "network_California_strength.png", width = 2500, height = 2000, res = 300)

# Load the California map using tigris
california <- counties(state = state, cb = TRUE, year = 2018)

# Plot the map using ggplot2 and overlay the graph
ggplot(data = california) +
  geom_sf(fill = "white", color = "black") +
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
               color = "cornflowerblue", size = (2/max(edges$weight))*edges$weight) +
  geom_point(data = nodes, aes(x = lon, y = lat, color = v_centr),
             size = sizes*5) +
  geom_text(data = nodes, aes(x = lon, y = lat, label = nodeLabel), 
            size = 2) +
  theme_minimal() +
  scale_color_gradient(name = "Strength",
                       low = "green", high = "red") +
  labs(x = "Lon", y = "Lat") +
  ggtitle(paste("SCI Network", state, sep=" ")) +
  theme(plot.title = element_text(color="black", size=20, hjust=0.5, face="bold.italic"))

dev.off()

#Generalize the example to obtain a function which can plot the Network for a chosen state


plot_network <- function(state) {
  nodes <- nodes_by_state[[state]]
  edges <- edges_by_state[[state]]
  
  g <- graph.data.frame(edges, directed=F, vertices=as.character(nodes$nodeID))
  if((any(E(simplify(g))!=E(g)))) {stop("Some graph is not simple!")}
  
  layout <- matrix(0, vcount(g), 2)
  layout[,1] <- dms_to_decimal(nodes$longitude)
  layout[,2] <- dms_to_decimal(nodes$latitude)
  #temporary
  names(layout[,1]) <- NULL
  names(layout[,2]) <- NULL
  #
  
  v_centr <- strength(g, vids = V(g), weights = E(g)$weight) 
  #v_centr <- harmonic_centrality(g, normalize=F,weights = E(g)$weights)
  names(v_centr) <- NULL
  
  colors <- vec2pal(v_centr, mypal)
  sizes <- v_centr/max(v_centr)
  
  # Prepare the nodes and edges as data.frames
  nodes <- nodes %>% rename(lon=longitude,lat=latitude) %>% mutate(lon=layout[,1],lat=layout[,2])
  edges <- edges %>% rename(from=nodeID_from,to=nodeID_to)
  
  # Merge the nodes coordinates with edges data to get x and y for each endpoint
  edges <- edges %>%
    left_join(nodes, by = c("from" = "nodeID")) %>%
    rename(x = lon, y = lat) %>%
    left_join(nodes, by = c("to" = "nodeID")) %>%
    rename(xend = lon, yend = lat)
  
  # Load the state map using tigris
  state_data <- counties(state = state, cb = TRUE, year = 2018)
  
  
  # For representation's sake we exclude Aleutians West Census Area (FIPS: 02016) from Alaska
  if (state=="Alaska")
  {
    exclude_fips <- c("02016")
    
    state_data <- state_data[!(state_data$GEOID %in% exclude_fips), ]
  }
  
  #correct label
  state_label <- ifelse(state==78,"U.S. Virgin Islands",state)
  
  # Plot the map using ggplot2 and overlay the graph
  ggplot(data = state_data) +
    geom_sf(fill = "white", color = "black") +
    geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "cornflowerblue", size = (2/max(edges$weight))*edges$weight) +
    geom_point(data = nodes, aes(x = lon, y = lat, color = v_centr),
               size = sizes*5) +
    geom_text(data = nodes, aes(x = lon, y = lat, label = nodeLabel), 
              size = 2) +
    theme_minimal() +
    scale_color_gradient(name = "Strength",
                         low = "green", high = "red") +
    labs(x = "Lon", y = "Lat") +
    ggtitle(paste("SCI Network", state_label, sep=" ")) +
    theme(plot.title = element_text(color="black", size=20, hjust=0.5, face="bold.italic"))
}


#We produce the output files and plots


#pdf("network_plots_strength.pdf")
for (state in states) {
  png(filename = sprintf("Plot_strength/network_plots_strength_%s.png", state), 
      width = 5000, height = 4000, res = 300)
  N <- plot_network(state)
  print(N)
  dev.off()
}



##Data analysis

#Plot the empirical strength distribution and compute descriptors for each state

N_vec <- rep(0,length(states)) #number of nodes
E_vec <- rep(0,length(states)) #number of edges
hc_vec <- list() #harmonic centrality
st_vec <- list() #strength




pdf("strenght_plots.pdf")
for(i in 1:length(states)) {
  state <- states[i]
  nodes <- nodes_by_state[[state]]
  edges <- edges_by_state[[state]]
  
  g <- graph.data.frame(edges, directed=F, vertices=nodes$nodeID)
  g <- simplify(g)
  
  
  st <- strength(g, vids = V(g), weights = E(g)$weight) #strength
  
  #We store the descriptors for each network
  N_vec[i] <- length(V(g))
  E_vec[i] <- length(E(g))
  st_vec[[i]] <- st
  hc_vec[[i]] <- harmonic_centrality(g,normalized = F,weights=E(g)$weights) #harmonic centrality
  
  #We plot the strength distribution
  plot <- ggplot(data.frame(weighted_degree = st), aes(x = weighted_degree)) +
    geom_histogram(binwidth = 2*10e6, fill = "lightblue", color = "black") + 
    theme_minimal()+
    labs(title = paste("Strength Distribution", state_label[i], sep=" "), x = "Strength", y = "Frequency") +
    theme(plot.title = element_text(color="black", size=15, hjust=0.5, face="bold.italic"))
  
  print(plot)
}
dev.off()


#We plot the strength distributions for the largest networks

big_states <- states[N_vec>=100]
big_st_vec <- st_vec[N_vec>=100]

pdf("big_st_plots.pdf")
for(i in 1:length(big_states)) {
  state <- big_states[i]
  
  #We plot the strength distribution
  plot <- ggplot(data.frame(strength = big_st_vec[[i]]), aes(x = strength)) +
    geom_histogram(binwidth = 2*10e6, fill = "lightblue", color = "black") + 
    theme_minimal()+
    labs(title = paste("Strength Distribution", state, sep=" "), x = "Strength", y = "Frequency") +
    theme(plot.title = element_text(color="black", size=15, hjust=0.5, face="bold.italic"))
  
  print(plot)
}
dev.off()




#Some networks are not complete

idx <- which(N_vec*(N_vec-1)/2!=E_vec)
cat(sprintf("The non-complete Networks are relative to the regions of %s, %s, %s",states[idx[1]],states[idx[2]],states[idx[3]]))
cat("\n")
cat(sprintf("With respective number of nodes %d, %d, %d",N_vec[idx[1]],N_vec[idx[2]],N_vec[idx[3]]))
cat("\n")


#Alanalysis for the entire U.S.

tot_nodes <- county_nodes %>% subset(select=-state)

tot_edges <- sci
tot_edges <- tot_edges %>% rename(state_code=source_state) %>% mutate(state_code = as.character(state_code))

tot_edges <- tot_edges %>%
  left_join(fips_codes, by = "state_code")


tot_edges <- tot_edges %>% subset(select=c(nodeID_from,nodeID_to,scaled_sci,state_name)) %>% rename(weight=scaled_sci)



#Delete self-loops and multiple edges in the file
# Remove self-loops
tot_edges <- tot_edges %>% 
  filter(nodeID_from != nodeID_to)

# Remove duplicate edges (ignoring direction)
tot_edges  <- tot_edges  %>%
  rowwise() %>%
  mutate(from = min(nodeID_from, nodeID_to), 
         to = max(nodeID_from, nodeID_to)) %>%
  ungroup() %>%
  distinct(from, to, .keep_all = TRUE) %>% select(index=-c(from,to))



##Create nodes and edges files

# Create a directory if it does not exist
create_directory <- function(state) {
  dir_path <- paste("data/",state_label[which(state == states)], sep = "")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  return(dir_path)
}




#for all US
#edges
write.csv(tot_edges, file = "data/tot_edges.csv", row.names = FALSE)
#nodes
write.csv(tot_nodes, file = "data/tot_nodes.csv", row.names = FALSE)

#for each state
#edges
invisible(sapply(states, function(state) {
  dir_path <- create_directory(state)
  file_path <- file.path(dir_path, "edges.csv")
  write.csv(edges_by_state[[state]], file = file_path, row.names = FALSE)
}))

#nodes
invisible(sapply(states, function(state) {
  dir_path <- create_directory(state)
  file_path <- file.path(dir_path, "nodes.csv")
  write.csv(nodes_by_state[[state]], file = file_path, row.names = FALSE)
}))


#Overall Network

nodes <- tot_nodes
edges <- tot_edges

# For representation purposes filter nodes with bounds close to the U.S.
nodes <- nodes %>% mutate(latitude=dms_to_decimal(latitude),longitude=dms_to_decimal(longitude))
#temporary
names(nodes$latitude) <- NULL
names(nodes$longitude) <- NULL

nodes <- nodes %>% filter(longitude >= -125 & longitude <= -65 & latitude >= 24 & latitude <= 50)

# Keep only the edges that connect nodes within those nodes
edges <- edges %>%
  filter(nodeID_from %in% nodes$nodeID & nodeID_to %in% nodes$nodeID)

g <- graph.data.frame(edges, directed=F, vertices=nodes$nodeID)
g <- simplify(g)

layout <- matrix(0, vcount(g), 2)
layout[,1] <- nodes$longitude
layout[,2] <- nodes$latitude
#temporary
names(layout[,1]) <- NULL
names(layout[,2]) <- NULL


#v_centr <- harmonic_centrality(g,normalized = F,weights=E(g)$weights)
v_centr <- strength(g, vids = V(g), weights = E(g)$weight) 
names(v_centr) <- NULL

colors <- vec2pal(v_centr, mypal)
sizes <- v_centr/max(v_centr)

# Prepare the nodes and edges as data.frames
nodes <- nodes %>% rename(lon=longitude,lat=latitude) %>% mutate(lon=layout[,1],lat=layout[,2])
edges <- edges %>% rename(from=nodeID_from,to=nodeID_to)

# Merge the nodes coordinates with edges data to get x and y for each endpoint
edges <- edges %>%
  left_join(nodes, by = c("from" = "nodeID")) %>%
  rename(x = lon, y = lat) %>%
  left_join(nodes, by = c("to" = "nodeID")) %>%
  rename(xend = lon, yend = lat)

# Load the state map using tigris
data <- counties(cb = TRUE, year = 2018)

# We exclude the same counties from the map for better representation
keep_fips <- nodes$nodeID
data <- data[data$GEOID %in% keep_fips, ]




#png(filename = "US_network2.png", width = 5000, height = 4000, res = 300)
# Plot the map using ggplot2 and overlay the graph
ggplot(data = data) +
  geom_sf(fill = "white", color = "black") +
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
               color = "cornflowerblue", size = (10/max(edges$weight))*edges$weight) +
  geom_point(data = nodes, aes(x = lon, y = lat, color = v_centr),
             size = sizes) +
  #geom_text(data = nodes, aes(x = lon, y = lat, label = nodeLabel), 
  #         size = 2) +
  theme_minimal() +
  scale_color_gradient(name = "Strength",
                       low = "green", high = "red") +
  labs(x = "Lon", y = "Lat") +
  ggtitle("SCI Network U.S.") +
  theme(plot.title = element_text(color="black", size=20, hjust=0.5, face="bold.italic"))
#dev.off()



#Total strength distribution

st <- strength(g, vids = V(g), weights = E(g)$weight) #strength

#We plot the strength distribution
plot <- ggplot(data.frame(weighted_degree = st), aes(x = weighted_degree)) +
  geom_histogram(binwidth = 2*10e6, fill = "lightblue", color = "black") + 
  theme_minimal()+
  labs(title = "Strength Distribution U.S.", x = "Strength", y = "Frequency") +
  theme(plot.title = element_text(color="black", size=15, hjust=0.5, face="bold.italic"))

print(plot)

