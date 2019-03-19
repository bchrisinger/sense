###########\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# R CODE TO PROCESS AND VISUALIZE DATA FROM DISCOVERY TOOL AND BIOSENSOR
# Citation: Chrisinger, B. W., & King, A. C. (2018). Stress experiences in neighborhood and social environments (SENSE): a pilot study to integrate the quantified self with citizen science to improve the built environment and health. International journal of health geographics, 17(1), 17.
# Paper available at: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5989430/
# Contact: Ben Chrisinger benjamin.chrisinger@spi.ox.ac.uk
# Last updated: 9 FEBRUARY, 2019
###########\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  

required_packages = c("leaflet","data.table","sp","rgdal","maptools","KernSmooth","XML","anytime","htmlwidgets","raster","ggplot2","formattable","lubridate","stringr")
install.packages(required_packages,dependencies=T)

#### SETUP INSTRUCTIONS:
# REQUIRED INPUTS:
# > inputs_folder = a single folder containing ONLY participant data folders, downloaded from the Empatica Connect web server.
# >>> NOTES: a valid .gpx file must be placed inside each participant data folder, and its name must match the folder name 
# >>> e.g., if a participant data folder is called "ParticipantXYZ", the .gpx file must be named "ParticipantXYZ.gpx"
# >>> folder names become ParticipantIDs
# > project_name = any character string that will be set as a project name
# > time_zone = a valid time zone; for formats, see: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones 

#OPTIONAL INPUTS:
# > trim_start = if data are going to be trimmed, the number of minutes to trim from beginning of EDA data; default is 5 min.
# > trim_end = the end time point if not the minimum overlapping time of GPS and Empatica E4 data; default is 0 min.
# > DTsummary_path = if Discovery Tool data are available, this will be a spreadsheet with summary data for all DT photos (e.g., lat/lon, transcription, etc.)
# > KML_path = if user wants to snap all GPS points to a path they've generated (like a walking/driving route from Google Maps)

#OUTPUT:
# > leaflet map(s) of participant data as .html file(s), saved to inputs_folder.
# >>> NOTE: this function will throw at least two Warning messages related to duplicated column names. This is ok.

sense_maps = function(inputs_folder,project_name,time_zone,trim_start=5,trim_end=0,DTsummary_path=NULL,KML_path=NULL) {
  #load required libraries
  library("leaflet")
  library("data.table")
  library("sp")
  library("rgdal")
  library("maptools")
  library("KernSmooth")
  library("XML")
  library("anytime")
  library("htmlwidgets")
  library("raster")
  library("ggplot2")
  library("formattable")
  library("lubridate")
  library("stringr")
  
  if (hasArg(DTsummary_path)) {
    DTsummary = read.csv(DTsummary_path,header=T, na.strings = "NA")
    DTsummary$Y_coords = DTsummary$Latitude
    DTsummary$X_coords = DTsummary$Longitude }
   
    #set DT photo marker colors by positive/negative rating
    getColor <- function() {
      sapply(DTsummary$GoodOrBad, function(GoodOrBad) {
        if(GoodOrBad == "Good") {
          "blue"
        } else if(GoodOrBad == "Bad") {
          "red"
        } else {
          "gray"
        } })
    }
  
    if (hasArg(KML_path)) {
      #KML_path KML file of walking path on sidewalks (e.g., straight lines)
      path = KML_path
      # create spatial lines from the kml data for sidewalks
      walk0 = as.data.frame(getKMLcoordinates(path,ignoreAltitude = T))
      names(walk0) = c("lon","lat")
      walk = Line(walk0)
      walkLine <- Line(walk0) %>% list() %>% Lines(ID='walk=line') %>% list() %>% SpatialLines() 
      walkLine = SpatialLines(walkLine@lines,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    }
  
  if (hasArg(KML_path) & hasArg(DTsummary_path)) {
    #snap DT data to sidewalk
      DTpoints = SpatialPointsDataFrame(data.frame(DTsummary$Longitude,DTsummary$Latitude), data = DTsummary)
      proj4string(DTpoints)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
      snapDT = snapPointsToLines(DTpoints, walkLine, maxDist = 10) 
      snapDTtab = as.data.table(cbind(snapDT,snapDT@coords)) 
      #remove rows with missing data
      snapDTtab = remove_missing(snapDTtab, na.rm = T, vars = c("lat", "lon"))
      #set snapped output as new inputs for mapping
      DTsummary$X_coords = snapDTtab$X.1
      DTsummary$Y_coords = snapDTtab$Y.1 
      }
  
  if (hasArg(DTsummary_path)) {
  #set up positive/negative density plots
  pospts=SpatialPointsDataFrame(data.frame(DTsummary[DTsummary$GoodOrBad=="Good",]$X_coords,
                                           DTsummary[DTsummary$GoodOrBad=="Good",]$Y_coords),
                                data=DTsummary[DTsummary$GoodOrBad=="Good",])
  negpts=SpatialPointsDataFrame(data.frame(DTsummary[DTsummary$GoodOrBad=="Bad",]$X_coords,
                                           DTsummary[DTsummary$GoodOrBad=="Bad",]$Y_coords),
                                data=DTsummary[DTsummary$GoodOrBad=="Bad",])
  #dat_shp=SpatialPointsDataFrame(data.frame(dat[,-14]$X_coords,dat[,-14]$Y_coords),data=dat[,-14])
  #proj4string(dat_shp)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  dtpolys = as.data.table(DTsummary)
  dtpolys = na.omit(dtpolys, c("X_coords","Y_coords"))
  #contour lines
  kde_pos <- bkde2D(dtpolys[dtpolys$GoodOrBad == "Good",][ , list(X_coords,Y_coords)],bandwidth=c(.00015, .00015), gridsize = c(800,800))
  CL_pos <- contourLines(kde_pos$x1 , kde_pos$x2 , kde_pos$fhat) 
  kde_neg <- bkde2D(dtpolys[dtpolys$GoodOrBad == "Bad",][ , list(X_coords,Y_coords)],bandwidth=c(.00015, .00015), gridsize = c(800,800))
  CL_neg <- contourLines(kde_neg$x1 , kde_neg$x2 , kde_neg$fhat) 
  #extract contour line levels
  LEVS_pos <- as.factor(sapply(CL_pos, `[[`, "level"))
  NLEV_pos <- length(levels(LEVS_pos))
  LEVS_neg <- as.factor(sapply(CL_neg, `[[`, "level"))
  NLEV_neg <- length(levels(LEVS_neg))
  ## CONVERT CONTOUR LINES TO POLYGONS
  pgons_pos <- lapply(1:length(CL_pos), function(i)
    Polygons(list(Polygon(cbind(CL_pos[[i]]$x, CL_pos[[i]]$y))), ID=i))
  spgons_pos = SpatialPolygons(pgons_pos)
  proj4string(spgons_pos)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  pgons_neg <- lapply(1:length(CL_neg), function(i)
    Polygons(list(Polygon(cbind(CL_neg[[i]]$x, CL_neg[[i]]$y))), ID=i))
  spgons_neg = SpatialPolygons(pgons_neg)
  proj4string(spgons_neg)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  }
  
  setwd(inputs_folder)
  partIDs = list.files()
  datalist = list()
  maplist = list()
  
for (i in seq(1:nrow(data.frame(partIDs)))) {
  
  PartID = partIDs[[i]]
  GPXfilename = paste(PartID,"/",PartID,".gpx",sep = "",collapse=NULL)
  
#set up helper function for GPX file reading  
  shift.vec <- function (vec, shift) {
    if(length(vec) <= abs(shift)) {
      rep(NA ,length(vec))
    }else{
      if (shift >= 0) {
        c(rep(NA, shift), vec[1:(length(vec)-shift)]) }
      else {
        c(vec[(abs(shift)+1):length(vec)], rep(NA, abs(shift))) } } }
  
  # Parse the GPX file
  pfile <- htmlTreeParse(GPXfilename,
                         error = function (...) {}, useInternalNodes = T)
  # Get all elevations, times and coordinates via the respective xpath
  elevations <- as.numeric(xpathSApply(pfile, path = "//trkpt/ele", xmlValue))
  times <- xpathSApply(pfile, path = "//trkpt/time", xmlValue)
  coords <- xpathSApply(pfile, path = "//trkpt", xmlAttrs)
  # Extract latitude and longitude from the coordinates
  lats <- as.numeric(coords["lat",])
  lons <- as.numeric(coords["lon",])
  # Put everything in a dataframe and get rid of old variables
  geodf <- data.frame(lat = lats, lon = lons, ele = elevations, time = times)
  rm(list=c("elevations", "lats", "lons", "pfile", "times", "coords"))
  
  # Shift vectors for lat and lon so that each row also contains the next position.
  geodf$lat.p1 <- shift.vec(geodf$lat, -1)
  geodf$lon.p1 <- shift.vec(geodf$lon, -1)
  # Calculate distances (in metres) using the function pointDistance from the ‘raster’ package.
  # Parameter ‘lonlat’ has to be TRUE!
  geodf$dist.to.prev <- apply(geodf, 1, FUN = function (row) {
    pointDistance(c(as.numeric(row["lat.p1"]),
                    as.numeric(row["lon.p1"])),
                  c(as.numeric(row["lat"]), as.numeric(row["lon"])),
                  lonlat = T)
  })
  
  # Shift the time vector, too.
  geodf$time.p1 <- shift.vec(geodf$time, -1)
  geodf$time <- as.POSIXct(geodf$time, format = "%Y-%m-%dT%H:%M:%OS")
  geodf$time <- force_tz(geodf$time, tzone=time_zone)

  EDA_FileName = paste(PartID,"/","EDA.csv",sep = "",collapse=NULL)
  ACC_FileName = paste(PartID,"/","ACC.csv",sep = "",collapse=NULL)
  BVP_FileName = paste(PartID,"/","BVP.csv",sep = "",collapse=NULL)
  Temp_FileName = paste(PartID,"/","TEMP.csv",sep = "",collapse=NULL)
  Tags_FileName = paste(PartID,"/","tags.csv",sep = "",collapse=NULL)
  
  data_output = list()
  data_list = list(EDA_FileName,ACC_FileName,BVP_FileName,Temp_FileName)

  for (x in seq(length(data_list))) {
    #set up eda dataframe with timestamps
    data_file = read.csv(data_list[[x]],header=F)
    start = data_file[1,1]
    interval = data_file[2,1]
    
    data = data.frame(data_file[3:nrow(data_file),])
    obvs = nrow(data)-1
    totsec = obvs/interval
    totmin = totsec/60
    end = start + totsec
    
    time = anytime(seq(start,end,1/interval),"UTC")
    time = data.frame(with_tz(time,time_zone))
    time = setNames(time,"time")

    #trim Empatica data based on user-defined start time
    datadf = data.frame("data"=data[(interval*60*trim_start):nrow(data),])
    time = data.frame("time"=time[(interval*60*trim_start):nrow(time),])
    datadf1 = cbind(time,datadf)
    
    #determine the minimum viable overlap between datasets
    startTime = max(min(datadf1$time,na.rm=T),min(geodf$time,na.rm=T))
    stopTime =  min(max(geodf$time),max(datadf1$time))

    datadf = datadf1[(datadf1$time>startTime & datadf1$time<stopTime),]

    
    normalize = function(data) {
      out = (data - min(data,na.rm=T))/(max(data)-min(data,na.rm=T))
      return(out)
    }
    
    #normalizing/smoothing data
    if (ncol(datadf)<3) {
    datadf$normal = normalize(datadf$data)
    datadf$lowess = lowess(datadf[,2], f = 0.2)$y
    datadf$scale = scale(datadf[,2])
      } else { 
        datadf$normal = NA 
        datadf$lowess = NA
        datadf$scale = NA}

   datadf$lowess = ifelse(ncol(datadf)<4,lowess(datadf[,2], f = 0.2)$y,NA)
   datadf$scale = ifelse(ncol(datadf)<5,scale(datadf[,2]),NA)

    data_output[[x]] = data.frame(datadf)
  }
  
  eda_df = data_output[[1]]
  acc_df = data_output[[2]]
  bvp_df = data_output[[3]]
  tmp_df = data_output[[4]]
  
  
  #set up tagged events  
  if (file.info(Tags_FileName)$size > 0) {
  tags = read.csv(Tags_FileName,header=F)
  tags$V1 = anytime(tags$V1,tz=time_zone)
  eda_df$tag =  ifelse(as.character(eda_df$time) %in% as.character(tags$V1),"1","0") }
  
  #merge geodf and biodf by time, discard incomplete cases
  merge1 = merge(x = geodf[(geodf$time>startTime & geodf$time<stopTime),],
                 y = eda_df[(eda_df$time>startTime & eda_df$time<stopTime),],
                 by.x = "time", by.y = "time", all.x = T)
  
  merge2 = merge(x = merge1[(merge1$time>startTime & merge1$time<stopTime),],
                 y = acc_df[(acc_df$time>startTime & acc_df$time<stopTime),],
                 by.x = "time", by.y = "time", all.x = T)
  
  merge3 = merge(x = merge2[(merge2$time>startTime & merge2$time<stopTime),],
                 y = bvp_df[(bvp_df$time>startTime & bvp_df$time<stopTime),],
                 by.x = "time", by.y = "time", all.x = T)
  
  merge4 = merge(x = merge3[(merge3$time>startTime & merge3$time<stopTime),],
                 y = tmp_df[(tmp_df$time>startTime & tmp_df$time<stopTime),],
                 by.x = "time", by.y = "time", all.x = T)
  
  merge4$PartID = PartID
  
  if (file.info(Tags_FileName)$size > 0) {
    alldf = setNames(merge4,c("time","lat","lon","ele","lat.p1","lon.p1","dist.to.prev","time.p1",
                            "dat_eda","norm_eda","lowess_eda","scale_eda","tag",
                            "dat_acc1","dat_acc2","dat_acc3","norm_acc","lowess_acc","scale_acc",
                            "dat_bvp","norm_bvp","lowess_bvp","scale_bvp",
                            "dat_tmp","norm_tmp","lowess_tmp","scale_tmp","PartID")) 
  } else {
    alldf = setNames(merge4,c("time","lat","lon","ele","lat.p1","lon.p1","dist.to.prev","time.p1",
                              "dat_eda","norm_eda","lowess_eda","scale_eda",
                              "dat_acc1","dat_acc2","dat_acc3","norm_acc","lowess_acc","scale_acc",
                              "dat_bvp","norm_bvp","lowess_bvp","scale_bvp",
                              "dat_tmp","norm_tmp","lowess_tmp","scale_tmp","PartID"))
  }
  
  alldf = alldf[!is.na(alldf$lat) & !is.na(alldf$lon),]
  alldf$counter = alldf$time- alldf$time[1]

    if (trim_end>0) {
      alldf = data.frame(alldf[1:(which(alldf$time == force_tz(anytime(trim_end),time_zone))),])
     } else { 
       alldf = data.frame(alldf)
    }
 
  # Use all available Empatica E4 data in a linear model to estimate EDA outcome 
  attach(alldf)
  
  lm0 = lm(norm_eda ~ time + ele + dat_acc1 + dat_acc2 + dat_acc3 + dat_tmp + dat_bvp, data=alldf)
  alldf = cbind(alldf,"residuals"=lm0$residuals)
  alldf = cbind(alldf,"fitted"=lm0$fitted.values)

  alldf$resid_sds = alldf$residuals / (sd(alldf$residuals))

  r = abs(alldf$residuals) / (sd(alldf$residuals))
  alldf$resid_brks = 0
 if (!is.na(table(alldf$residuals>0 & r>=2)[2])) {
  alldf[alldf$residuals>0 & r>=2,]$resid_brks = "f. 2+ Above"
  alldf[alldf$residuals>0 & r<2 & r>=1,]$resid_brks = "e. 1-2 Above"
  alldf[alldf$residuals>0 & r<1 & r>0,]$resid_brks = "d. <1 Above" 
 } else {
   alldf[alldf$residuals>0 & r<2 & r>=1,]$resid_brks = "e. 1-2 Above"
   alldf[alldf$residuals>0 & r<1 & r>0,]$resid_brks = "d. <1 Above" 
 }
 
 if (!is.na(table(alldf$residuals<0 & r>=2)[2])) {
  alldf[alldf$residuals<0 & r<1 & r>0,]$resid_brks = "c. <1 Below"  
  alldf[alldf$residuals<0 & r<2 & r>=1,]$resid_brks = "b. 1-2 Below"
  alldf[alldf$residuals<0 & r>=2,]$resid_brks = "a. 2+ Below"
 } else {
    alldf[alldf$residuals<0 & r<1 & r>0,]$resid_brks = "c. <1 Below"  
    alldf[alldf$residuals<0 & r<2 & r>=1,]$resid_brks = "b. 1-2 Below"
  }
  
  
  ###########\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  # MAP PARTICIPANT DATA
  # code below will map one participant's data. can be looped for multiple datasets.
  ###########\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  
  #create color palette for normalized EDA bothdfa
  qpal_eda = colorBin("plasma", unique(alldf[alldf$PartID == PartID,]$dat_eda),pretty=F,bins=10)
  qpal_scaleda = colorBin("plasma", unique(alldf[alldf$partID == PartID]$scale_eda),pretty=T,bins=50)
  qpal_normEDA = colorBin("plasma", unique(alldf[alldf$partID == PartID]$normEDA),pretty=T,bins=100)
  qpal_scale_norm_eda = colorBin("plasma", unique(alldf[alldf$partID == PartID]$scale_norm_eda),pretty=T,bins=100)
  qpal_lowessEDA = colorBin("plasma", unique(alldf[alldf$partID == PartID]$lowess.eda),pretty=F,bins=10)
  qpal_lowessNormEDA = colorBin("plasma", unique(alldf[alldf$partID == PartID]$lowess.norm.eda),pretty=F,bins=10)
  pal_legend = colorBin("plasma", seq(1:5),pretty=F,bins=5)
  resid_eda_sds = colorBin("plasma", unique(alldf$resid_sds),pretty=F,bins=6)
  fitted_eda = colorBin("plasma", alldf$fitted,pretty=F,bins=c(-2*sd(alldf$fitted),-1*sd(alldf$fitted),0,
                                                               1*sd(alldf$fitted),2*sd(alldf$fitted)))
  resid_eda = colorFactor("plasma",alldf$resid_brks)
                      
#compile leaflet map
  
  if (hasArg(DTsummary_path)) {
    
    #create pathnames for where photos will be stored. 
    #this example stores photos in a participant-specific online folder. 
    #it is possible to store photos locally. see this discussion on stackoverflow for more: https://stackoverflow.com/questions/36433899/image-in-r-leaflet-marker-popups
    photopath = paste("http://web.stanford.edu/~bchris/SENSE/",DTsummary[DTsummary$ParticipantID == PartID,]$Path,sep="")
    
    }
    
    ###########\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # SNAP GPS AND DISCOVERY TOOL DATA TO SPECIFIED WALKING PATH. 
    #NOTE: skip this section if you don't have walking path file for snapping GPS points.
    # must supply a KML file input
    ###########\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    if (hasArg(KML_path)) {
    #snap GPS data to sidewalk
    gpspoints = SpatialPointsDataFrame(data.frame(alldf$lon,alldf$lat),data = alldf)
    proj4string(gpspoints)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    snapgps = snapPointsToLines(gpspoints, walkLine, maxDist = 10)
    snapgpstab = as.data.table(cbind(snapgps,snapgps@coords))
    #remove rows with missing data
    snapgpstab = remove_missing(snapgpstab, na.rm = T, vars = c("lat", "lon"))
    #set snapped output as new inputs for mapping
    alldf$Y_coords = as.numeric(snapgpstab$Y.1)
    alldf$X_coords = as.numeric(snapgpstab$X.1) 
    }
    
m = leaflet(alldf) %>%  
  #setView(lng = setlng, lat = setlat, zoom = 16.5) %>% #Set zoom at start
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addTiles(urlTemplate = "http://mt0.google.com/vt/lyrs=m&hl=en&x={x}&y={y}&z={z}&s=Ga",
           group = "Google", attribution = 'Google') %>%
  addProviderTiles(providers$OpenStreetMap, group = "Open Street Map") %>%
  addProviderTiles(providers$Stamen.Toner, group = "Black & White") %>%
  addCircleMarkers(lat = alldf$lat, 
                   lng = alldf$lon, 
                   radius = 5, fillOpacity = 0.75,
                   color = ~resid_eda(alldf$resid_brks),
                   popup = paste("participant ID: ",alldf$PartID,"<br>",
                                 "time:",alldf$time,"<br>",
                                 "EDA:",round(alldf$norm_eda,4),"<br>",
                                 "Residual (std. devs):", round(alldf$resid_sds,4)),
                   group = "Electrodermal Activity (EDA)") %>%
    addMarkers(lat = alldf[alldf$tag == "1",]$lat,lng = alldf[alldf$tag == "1",]$lon,
             popup = paste("<b>Participant:</b>",alldf$PartID, "<br>",
                           "<b>Time:</b>", alldf[alldf$tag == "1",]$time),
             popupOptions = popupOptions(closeButton = T, closeOnClick = T,minWidth=50,maxWidth=300,
                                         maxHeight=300),
             group = "Tagged Events") %>%
  addMiniMap(tiles = providers$Stamen, #Add inset map
           toggleDisplay = TRUE, zoomLevelFixed = 11, width = 230,height = 230,position="bottomleft") %>%
  addLegend("bottomright", pal =resid_eda,values=~resid_brks,
            labFormat = labelFormat(digits=1),
            title = paste("EDA Residual","<br>"),
            opacity = 1) %>%
  addLegend("topright",colors = list(), labels = list(),
            title=paste("SENSE Participant Map:",PartID,sep=" ")) %>%
  addMeasure(position="topright","primaryLengthUnit"="meters",primaryAreaUnit="sqmeters") 
  

if (hasArg(DTsummary_path)) {
 m = m %>% 
    addAwesomeMarkers(lat = DTsummary[DTsummary$ParticipantID == PartID,]$Y_coords,
                    lng = DTsummary[DTsummary$ParticipantID == PartID,]$X_coords,
                    icon= awesomeIcons(
                      icon = "ion-android-camera",
                      iconColor = 'white',
                      library = 'ion',
                      markerColor = getColor(DTsummary[DTsummary$ParticipantID == PartID,])),
                    label = paste(DTsummary[DTsummary$ParticipantID == PartID,]$PhotoNum),
                    labelOptions = labelOptions(noHide = T, direction = "right", offset = c(12,-15),
                                                style = list(
                                                  "color" = "black",
                                                  "font-family" = "arial",
                                                  "font-style" = "bold",
                                                  "box-shadow" = "2px 2px rgba(0,0,0,0.25)",
                                                  "font-size" = "8.5px",
                                                  "border-color" = "rgba(0,0,0,0.5)"
                                                )),
                    popup = paste("<b>Walk:</b>",DTsummary[DTsummary$ParticipantID == PartID,]$ParticipantID, "<br>",
                                  "<b>Transcription:</b>", DTsummary[DTsummary$ParticipantID == PartID,]$Transcription, "<br>",
                                  "<b>Feeling:<b/>",DTsummary[DTsummary$ParticipantID == PartID,]$Feeling,"<br>",
                                  "<img src = ", photopath, ">"),
                    popupOptions = popupOptions(closeButton = T, closeOnClick = T,minWidth=300,maxWidth=300,
                                                maxHeight=300),
                    options = markerOptions(riseOnHover = T),
                    group = "Discovery Tool Data") %>%
   addLayersControl(position ="bottomleft",
                    baseGroups = c("Black & White","Open Street Map", "Google","Satellite","Walk Path"),
                    overlayGroups = c("Negative Clusters", "Positive Clusters",
                                      "Electrodermal Activity (EDA)","Tagged Events","Discovery Tool Data"),
                    options = layersControlOptions(collapsed = F)) %>%
   hideGroup(c("Negative Clusters", "Positive Clusters","Tagged Events")) 
 
 saveWidget(m,paste(PartID,".html",sep=""))  
 
 } else {
    m = m %>% addLayersControl(position ="bottomleft",
                     baseGroups = c("Black & White","Open Street Map", "Google","Satellite","Walk Path"),
                     overlayGroups = c("Electrodermal Activity (EDA)","Tagged Events"),
                     options = layersControlOptions(collapsed = F)) %>%
      hideGroup(c("Tagged Events"))
    
    saveWidget(m,paste(PartID,".html",sep=""))  }
  
  datalist[[i]] = alldf
  maplist[[i]] = m

}}




