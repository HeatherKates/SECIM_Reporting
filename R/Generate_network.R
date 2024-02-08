#Generate plot function
generate_network <- function(i, mydata,dir,dir_word) {
  mydata[[paste("vs",dir,sep=".")]][[i]] <- V(mydata[[paste(dir,"graph",sep=".")]][[i]]) #Get node list - this has names
  mydata[[paste("es",dir,sep=".")]][[i]]<- as.data.frame(get.edgelist(mydata[[paste(dir,"graph",sep=".")]][[i]])) # Get edgelist
  #This has the names and whether it was input
  mydata[[paste("node.data",dir,sep=".")]][[i]] <- get.data.frame(mydata[[paste(dir,"graph",sep=".")]][[i]],what="vertices") # Get node attribute dataframe
  #Add shape data
  mydata[[paste("node.data",dir,sep=".")]][[i]] <- mydata[[paste("node.data",dir,sep=".")]][[i]] %>%
    mutate(input = as.character(input))
  
  mydata[[paste("node.data",dir,sep=".")]][[i]] <- mydata[[paste("node.data",dir,sep=".")]][[i]] %>%
    mutate(com = case_when(
      endsWith(input, "TRUE") ~ "6",
      TRUE ~ as.character(com)
    ))
  
  mydata[[paste("node.data",dir,sep=".")]][[i]] <- mydata[[paste("node.data",dir,sep=".")]][[i]] %>%
    mutate(input = case_when(
      endsWith(input, "TRUE") ~ "Input",
      endsWith(input, "FALSE") ~ "",
    ))
  
  if(test_type=="nostats"){
    mydata[[paste("node.data",dir,sep=".")]][[i]] <- 
      mydata[[paste("node.data",dir,sep=".")]][[i]] %>%
      dplyr::mutate_at(vars(com), ~ stri_replace_all_regex(., pattern = as.character(c(1:6)), replacement = c("pathway", "module", "enzyme", "reaction", "compound", ">1.5 log2FC changed compound"), vectorize = FALSE))
  }else{
    mydata[[paste("node.data",dir,sep=".")]][[i]] <- 
      mydata[[paste("node.data",dir,sep=".")]][[i]] %>%
      dplyr::mutate_at(vars(com), ~ stri_replace_all_regex(., pattern = as.character(c(1:6)), replacement = c("pathway", "module", "enzyme", "reaction", "compound", "Significantly changed compound"), vectorize = FALSE))
  }
  
  mydata[[paste("Nv",dir,sep=".")]][[i]] <- length(mydata[[paste("vs",dir,sep=".")]][[i]]) #number of nodes
  mydata[[paste("Ne",dir,sep=".")]][[i]] <- length(mydata[[paste("es",dir,sep=".")]][[i]][[1]])# number of edges
  
  #Coordinates for nodes
  mydata[[paste("L",dir,sep=".")]][[i]] <- layout.fruchterman.reingold(mydata[[paste(dir,"graph",sep=".")]][[i]])
  mydata[[paste("Xn",dir,sep=".")]][[i]] <- mydata[[paste("L",dir,sep=".")]][[i]][,1]
  mydata[[paste("Yn",dir,sep=".")]][[i]] <- mydata[[paste("L",dir,sep=".")]][[i]][,2]
  
  #Add KEGG names as attribute to graph
  #Merge the vs$names with the KEGG_map df
  KEGG_map <- readRDS("../InputFiles/KEGG_map.df.RData")
  mydata[[paste("vs.names",dir,".")]][[i]] <- data.frame(names(mydata[[paste("vs",dir,sep=".")]][[i]]))
  
  mydata[[paste("vs.names",dir,".")]][[i]] <- mydata[[paste("vs.names",dir,".")]][[i]] %>%
    dplyr::mutate(num = as.numeric(row.names(mydata[[paste("vs.names",dir,".")]][[i]])))
  setnames(mydata[[paste("vs.names",dir,".")]][[i]],c("vs.name","num"))
  mydata[[paste("vs.names.KEGG",dir,".")]][[i]] <- merge(mydata[[paste("vs.names",dir,".")]][[i]],KEGG_map,by.x="vs.name",by.y="V1",all.x=TRUE)
  mydata[[paste("vs.names.KEGG",dir,".")]][[i]]  <- mydata[[paste("vs.names.KEGG",dir,".")]][[i]] [order(mydata[[paste("vs.names.KEGG",dir,".")]][[i]]$num), ]
  
  #Add confidence flag
  mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]]  <- 
    data.frame(Client_Data_Download[["report_results"]] %>% filter(contrast==contrasts[[i]]) %>% 
                 select("KEGG","Level") %>% filter(!KEGG=="") %>% distinct())
  
  mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]] <- 
    merge(mydata[[paste("vs.names.KEGG",dir,".")]][[i]],mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]],by.x="vs.name",by.y="KEGG")
  
  #Recode confidence
  mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]] <- mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]] %>% 
    mutate(Level=case_when(Level==1~"High Confidence ID",.default=as.character(Level)))
  
  mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]] <- mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]] %>% 
    mutate(Level=case_when(Level==3~"Low Confidence ID",.default=as.character(Level)))
  
  # Modify the column and assign it back to the same dataframe
  mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]] <- transform(mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]], report.results.Level = paste(V2, Level, sep = ";"))
  
  
  ####
  mydata[[paste("vs.names.KEGG.confidence.2",dir,".")]][[i]] <- mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]] %>% 
    filter(!vs.name %in% mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]][duplicated(mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]]$vs.name),]$vs.name) 
  
  mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]] <- rbind(mydata[[paste("vs.names.KEGG.confidence.2",dir,".")]][[i]],
                                                                    mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]] %>%
                                                                      filter(vs.name %in% mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]][duplicated(mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]]$vs.name),]$vs.name) %>% filter(Level=="High Confidence ID"))

  mydata[[paste("vs.names.KEGG",dir,".")]][[i]] <- merge(mydata[[paste("vs.names.KEGG",dir,".")]][[i]],mydata[[paste("vs.names.KEGG.confidence",dir,".")]][[i]],all.x=TRUE)
  
  mydata[[paste("vs.names.KEGG",dir,".")]][[i]] <- mydata[[paste("vs.names.KEGG",dir,".")]][[i]]  %>% 
    mutate(V2=case_when(!is.na(report.results.Level)~report.results.Level,.default=V2)) %>% dplyr::select(-report.results.Level)
  
  mydata[[paste("vs.names.KEGG",dir,".")]][[i]] <- mydata[[paste("vs.names.KEGG",dir,".")]][[i]][order(mydata[[paste("vs.names.KEGG",dir,".")]][[i]]$num), ]
  
  #add linebreaks for long descriptions (hovertext) 
  mydata[[paste("vs.names.KEGG",dir,".")]][[i]]$V2 <- gsub(";", "\n", mydata[[paste("vs.names.KEGG",dir,".")]][[i]]$V2 )
  
  #Creates the nodes (plots the points)
  mydata[[paste("node.data",dir,sep=".")]][[i]]$size <- ifelse(mydata[[paste("node.data",dir,sep=".")]][["element1"]]$input == "Input", 10, 5)

  network <- plot_ly(evaluate=TRUE)
  network <- add_trace(network,x = ~mydata[[paste("Xn",dir,sep=".")]][[i]], y = ~mydata[[paste("Yn",dir,sep=".")]][[i]], #Node points
                       mode = "markers", 
                       text = mydata[[paste("vs.names.KEGG",dir,".")]][[i]]$vs.name, 
                       hoverinfo = "text",
                       hovertext = mydata[[paste("vs.names.KEGG",dir,".")]][[i]]$V2,
                       color = as.factor(mydata[[paste("node.data",dir,sep=".")]][[i]]$com),
                       colors = c(brewer.pal(5,"Dark2"), "red"),
                       symbol = ~mydata[[paste("node.data",dir,sep=".")]][[i]]$input,
                       symbols = c("circle","square"),#size=10)
                      size = ~mydata[[paste("node.data",dir,sep=".")]][[i]]$size)
                      
  #Adds names to the nodes
  network <- network %>% add_text(x = ~mydata[[paste("Xn",dir,sep=".")]][[i]], y = ~mydata[[paste("Yn",dir,sep=".")]][[i]], text = mydata[[paste("vs.names.KEGG",dir,".")]][[i]]$vs.name,evaluate=TRUE)
  
  #Create edges
  #Initialize an empty list to store edge shapes
  mydata[[paste("edge_shapes",dir,sep=".")]][[i]] <- list()
  names(mydata[[paste("Xn",dir,sep=".")]][[i]])<- names(mydata[[paste("vs",dir,sep=".")]][[i]])
  names(mydata[[paste("Yn",dir,sep=".")]][[i]]) <-  names(mydata[[paste("vs",dir,sep=".")]][[i]])
  
  
  # Define the function to be applied for each j
  calculate_edge_shape <- function(j) {
    v0 <- as.character(mydata[[paste("es",dir,sep=".")]][[i]][j,]$V1)
    v1 <- as.character(mydata[[paste("es",dir,sep=".")]][[i]][j,]$V2)
    
    list(
      type = "line",
      line = list(color = "red", width = 0.3),
      x0 = mydata[[paste("Xn",dir,sep=".")]][[i]][v0],
      y0 = mydata[[paste("Yn",dir,sep=".")]][[i]][v0],
      x1 = mydata[[paste("Xn",dir,sep=".")]][[i]][v1],
      y1 = mydata[[paste("Yn",dir,sep=".")]][[i]][v1]
    )
  }
  
  # Use lapply to apply the function
  edge_shapes_list <- lapply(1:mydata[[paste("Ne",dir,sep=".")]][[i]], calculate_edge_shape)
  
  # Assign the result back to mydata
  mydata[["edge_shapes"]][[i]] <- edge_shapes_list
  
  
  mydata[[paste("edge_shapes",dir,sep=".")]][[i]] <-  mydata[["edge_shapes"]][[i]]  # Assign the list to your target variable
  mydata[[paste("axis",dir,sep=".")]][[i]] <- list(title = "", showgrid = FALSE, 
                                                   showticklabels = FALSE, zeroline = FALSE)
  
  if(test_type=="nostats"){
    network <- layout(
      network,
      title= paste("KEGG subnetwork for known metabolites", dir_word ,"abundant in\n",str_split(contrasts[[i]],pattern="-")[[1]][1],"relative to\n",
                   str_split(contrasts[[i]],pattern="-")[[1]][2]),
      shapes = force(mydata[[paste("edge_shapes",dir,sep=".")]][[i]]),
      xaxis = mydata[[paste("axis",dir,sep=".")]][[i]],
      yaxis = mydata[[paste("axis",dir,sep=".")]][[i]],
      width=800,height=800, #Added
      showlegend=TRUE,evaluate=TRUE)
  }else{
    network <- layout(
      network,
      title= paste("KEGG subnetwork for known metabolites significantly", dir_word ,"abundant in\n",str_split(contrasts[[i]],pattern="-")[[1]][1],"relative to\n",
                   str_split(contrasts[[i]],pattern="-")[[1]][2]),
      shapes = mydata[[paste("edge_shapes",dir,sep=".")]][[i]],
      xaxis = mydata[[paste("axis",dir,sep=".")]][[i]],
      yaxis = mydata[[paste("axis",dir,sep=".")]][[i]],
      width=800,height=800, #Added
      showlegend=TRUE,evaluate=TRUE)
    
  }
  return(network)
}