library(ggplot2)
# install_github("easyGgplot2", "kassambara")
library(easyGgplot2)
library(shiny)
library(shinydashboard)
library(dplyr)
library(leaflet)
library(RColorBrewer)


load("data.Rdata")



# FONCTIONS

# Normalisation
normalize = function(df){
  df2 = df
  for (i in 3:5){
    df2[,i] = (df2[,i]-min(df2[,i])) / (max(df2[,i])-min(df2[,i]))
  }
  return(df2)
}

normalize_val = function(x,i){
  return((x-min(df_base[,i])) / (max(df_base[,i])-min(df_base[,i])) )
}

K = function(x){
  return(15/16*(1-x^2)^2)
}

Kh = function(x){
  tmp = x/h
  tmp[abs(tmp)>1]=1
  tmp = K(tmp)
  return(tmp/h)
}

indic_oz = function(N,champ){
  Y = matrix(rep(df$oz,length(champ)),ncol=length(champ))
  tmp = t(matrix(rep(champ,N),ncol=N))
  return(Y>tmp)
}

indic_tp = function(N,champ){
  Y = matrix(rep(df$tp,length(champ)),ncol=length(champ))
  tmp = t(matrix(rep(champ,N),ncol=N))
  return(Y>tmp)
}


qnc = function(alpha,Fn,champ){
  return(champ[min(which(Fn<=alpha))])
}

gam_ncH = function(champ,Fn){
  tau = 1:9
  Sdiv = sum(log(tau))
  Snum = c()
  for (t in alpha/tau){
    Snum = c(Snum,qnc(t,Fn,champ))
  }
  return(sum(log(Snum)-log(Snum[1]))/Sdiv)
}

qnW = function(Fn,beta,champ){
  res = qnc(alpha,Fn,champ)*(alpha/beta)^gam_ncH(champ,Fn)
  return(res)
}





# SCRIPT


# Réglages
df = df[df$lat<62,] # Suppression de l'Alaska
df = df[order(df$oz),] # Tri par température
df$tp = (df$tp-32)*5/9 # Passage en Celsius
df_base = df
df = normalize(df)


# Paramètres
h = 0.12
N = dim(df)[1]
alpha = 0.001 # ~= 11 * log(N)/N
Npt = 60
Cpt = as.data.frame(table(df$tp))
Cpt$Freq = Cpt$Freq/N
Cpt$Var1 = as.numeric(as.character(Cpt$Var1))
Cpt = Cpt[(Cpt$Freq>0.001 & Cpt$Var1>=10 & Cpt$Var1<50),]
champ_tp = c(Cpt$Var1,seq(44,70,4))
champ_oz = 40:160/1000
J = 1:Npt/Npt
D = df[,3:5]
Mat = Kh(t(matrix(rep(J,N),ncol=N))-D$day)

# Calcul de l'indicatrice dès le début (ne dépend pas de lat et lng)
ind_tp = indic_tp(N,champ_tp)
ind_oz = indic_oz(N,champ_oz)


DF = as.data.frame(dplyr::summarise(dplyr::group_by(df_base,lat), lng=mean(lng)))
DF = cbind(rownames(DF),DF)
colnames(DF)[1]="station"
DF$station = paste("Station",DF$station)

dataplot = as.data.frame(dplyr::summarise(dplyr::group_by(df_base,lat), lng=mean(lng), tp=mean(tp), oz=mean(oz)))
dataplot = cbind(rownames(dataplot),dataplot)
colnames(dataplot)[1]="station"
dataplot$station = paste("Station",dataplot$station)

H_lat = h*(max(DF$lat)-min(DF$lat))
H_lng = h*(max(DF$lng)-min(DF$lng))





# GRAPHIQUES

# ESTIMATIONS

carte = function(){   # lieu <- couple (lng,lat)
  leaflet() %>%
    addProviderTiles(providers$OpenStreetMap.Mapnik,
                     options = providerTileOptions(noWrap = TRUE)
    ) %>%
    addCircles(data = DF, lng = ~lng, lat = ~lat, label = ~station
    )
}

graphe_day = function(lieu,l){   # lieu <- couple (lng,lat)
  
    beta = l*alpha
    latn = normalize_val(lieu[2],3)
    lngn = normalize_val(lieu[1],4)
    v = Kh(latn-D$lat)*Kh(lngn-D$lng)
    
    if (length(which(v!=0))==0){
      ggplot() + geom_label(aes(0,0),label="Hors du secteur",size=10,col="coral3")
    }
    
    else {
      
    M = t(v*Mat)
    SK = apply(M,1,sum)
  
  
    # Températures
    Ki_tp = M%*%ind_tp
    Fn_tp = Ki_tp/SK
    QtpW = apply(Fn_tp,1,qnW,beta=beta,champ=champ_tp)
    QtpC = apply(Fn_tp,1,qnc,alpha=alpha,champ=champ_tp)
    q1 = min(QtpW,QtpC)
  
  
    # Ozone
    Ki_oz = M%*%ind_oz
    Fn_oz = Ki_oz/SK
    QozW = apply(Fn_oz,1,qnW,beta=beta,champ=champ_oz)
    QozC = apply(Fn_oz,1,qnc,alpha=alpha,champ=champ_oz)
    
    tableau = data.frame(Day=365*J,QTPW=QtpW,QTPC=QtpC,QOZW=QozW,QOZC=QozC)
  
    p = ggplot(tableau, aes(Day,QTPW,QTPC,QOZW,QOZC)) +
      theme(panel.background=element_rect(fill="#CEE1FA"),
          plot.title=element_text(color='darkblue',size=10,face="bold"),
          axis.line.x=element_line(color="darkblue"),
          axis.text=element_text(color="darkblue"))
  
    p1 = p + theme(axis.ticks.x=element_blank(),
                 axis.title.x=element_text(color="white"),
                 axis.title.y=element_text(color="coral3",face="bold")) +
      geom_ribbon(aes(x=Day,ymin=q1-1,ymax=QTPW),fill="coral",alpha=0.3) +
      geom_line(aes(x=Day,y=QTPW),colour='coral3',lwd=1) +
      geom_line(aes(x=Day,y=QTPC),colour='chocolate') +
      labs(title="Quantile de Weissman - Température",y="Températures (°C)")
    
  
    p2 = p + theme(axis.title.x=element_text(color="darkblue"),
                 axis.title.y=element_text(color="#7555DB",face="bold")) +
      geom_area(aes(x=Day,y=QOZW),fill='#7555DB',alpha=0.3) +
      geom_line(aes(x=Day,y=QOZW),colour='#7555DB',lwd=1) +
      geom_line(aes(x=Day,y=QOZC),colour='blue') +
      labs(title="Quantile de Weissman - Ozone",x="Jour",y="Concentration d'ozone (ppm)")
  
  
    ggplot2.multiplot(p1,p2,cols=1) # Joindre 2 graphes en 1
    
    }
}


influ = function(lieu){   # lieu <- couple (lng,lat)
  DF2 = DF[,2:3]
  DF2$lng = abs(DF2$lng-lieu[1])
  DF2$lat = abs(DF2$lat-lieu[2])
  DF2 = DF2[DF2$lng<H_lng,]
  DF2 = DF2[DF2$lat<H_lat,]
  
  if (0 %in% dim(DF2)) {
    ggplot() + geom_label(aes(0,0),label="Hors du secteur",size=10,col="coral3")
  }
  
  else {
  DF2 = as.data.frame(apply(K(h*DF2),1,prod))
  DF2 = cbind(rownames(DF2),DF2)
  colnames(DF2) = c("Station","Influence")
  
  ggplot(DF2) + 
    geom_bar(aes(Station,weight=Influence,fill=Influence)) + coord_polar() +
    scale_fill_gradient(low="yellow", high="red") +
    theme(panel.background=element_rect(fill="#CEE1FA"),
          plot.title=element_text(color='darkblue',size=10,face="bold"),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = "darkblue"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line(size =2,linetype ='solid',colour = "white"),
          panel.grid.major.y = element_line(size =1,linetype ='dashed',colour = "white")) +
    labs(title="Influence des stations sur la localisation")
  }
}

# SYNTHESE

moy_map_tp = function(){   # lieu <- couple (lng,lat)
  bins <- c(0, 10, 20, 30, 40, Inf)
  pal <- colorBin("YlOrRd", domain = dataplot$tp, bins = bins)
  leaflet() %>%
    addProviderTiles(providers$OpenStreetMap.Mapnik,
                     options = providerTileOptions(noWrap = TRUE)
    ) %>%
    addCircles(data = dataplot, lng = ~lng, lat = ~lat, color= ~pal(tp),
               radius = ~tp^2*50, label = ~station, fillColor= ~pal(tp)
    )
}

moy_map_oz = function(){   # lieu <- couple (lng,lat)
  bins <- c(0, 0.01, 0.03, 0.05, Inf)
  pal <- colorBin("RdPu", domain = dataplot$oz, bins = bins)
  leaflet() %>%
    addProviderTiles(providers$OpenStreetMap.Mapnik,
                     options = providerTileOptions(noWrap = TRUE)
    ) %>%
    addCircles(data = dataplot, lng = ~lng, lat = ~lat, color= ~pal(oz),
               radius = ~(oz*1000)^3/2, label = ~station,  fillColor= ~pal(oz)
    )
}



# STATIONS

evol_station_tp = function(i){
  if (i<1 | i>427) ggplot() + geom_label(aes(0,0),label="Hors du secteur",size=10,col="coral3")
  else{
    long = DF$lng[i]
    lati = DF$lat[i]
    st = df_base[(df_base$lng==long & df_base$lat==lati),]
    ggplot(st) + geom_line(aes(x=day,y=tp),colour="coral",lwd=0.5) +
      geom_smooth(aes(x=day,y=tp),colour="coral3",lwd=1.5,method = "loess") +
      labs(x="Jour",y="Température (°C)")
  }
}


evol_station_oz = function(i){
  if (i<1 | i>427) ggplot() + geom_label(aes(0,0),label="Hors du secteur",size=10,col="coral3")
  else{
    long = DF$lng[i]
    lati = DF$lat[i]
    st = df_base[(df_base$lng==long & df_base$lat==lati),]
    ggplot(st) + geom_line(aes(x=day,y=oz),colour="#7555DB",lwd=0.5) +
      geom_smooth(aes(x=day,y=oz),colour="#7555DB",lwd=1.5,method = "loess") +
      labs(x="Jour",y="Concentration d'ozone (ppm)")
  }
}


# JOURS

evol_jour_tp = function(i){
  bins <- c(-Inf,0, 10, 20, 30, 40, Inf)
  pal <- colorBin("YlOrRd", domain = dataplot$tp, bins = bins)
  leaflet() %>%
    addProviderTiles(providers$OpenStreetMap.Mapnik,
                     options = providerTileOptions(noWrap = TRUE)
    ) %>%
    addCircles(data = df_base[df_base$day==i,], lng = ~lng, lat = ~lat, color= ~pal(tp),
               radius = ~tp^2*50, fillColor= ~pal(tp)
    )
}


evol_jour_oz = function(i){
  bins <- c(0, 0.01, 0.03, 0.05, Inf)
  pal <- colorBin("RdPu", domain = dataplot$oz, bins = bins)
  leaflet() %>%
    addProviderTiles(providers$OpenStreetMap.Mapnik,
                     options = providerTileOptions(noWrap = TRUE)
    ) %>%
    addCircles(data = df_base[df_base$day==i,], lng = ~lng, lat = ~lat, color= ~pal(oz),
               radius = ~(oz*1000)^3/2,  fillColor= ~pal(oz)
    )
}

localise = function(x,y){
  leaflet() %>%
    addProviderTiles(providers$OpenStreetMap.Mapnik,
                     options = providerTileOptions(noWrap = TRUE)
    ) %>% fitBounds(-120, 25, -65, 50) %>%
    addMarkers(lng = x,lat = y,
               label = paste("Position :\n(",x,",",y,")",sep=""),
               icon = list(
                 iconUrl = "http://icons.iconarchive.com/icons/martz90/circle/512/maps-icon.png",
                 iconSize = c(20,20)))
}




# SHINY APP



ui <- dashboardPage(
  dashboardHeader(title='VisualData  ~  TER'),
  dashboardSidebar(
    sidebarMenu(
      br(),br(),
      HTML("<h4 align='center'>ÉTUDE DE DONNÉES<br>
           ATMOSPHÉRIQUES<br>
           EXTRÊMES<br>
           (USA)</h4>"),
      br(),br(),
      menuItem("Estimations", tabName = "visu", icon = icon("bar-chart-o"),badgeLabel = "►", badgeColor = "green"),
      menuItem("Synthèse", tabName = "moy", icon = icon("spinner"),badgeLabel = "►", badgeColor = "olive"),
      menuItem("Stations", tabName="stat", icon = icon("plus"),badgeLabel = "►", badgeColor = "olive"),
      menuItem("Jours", tabName="jour", icon = icon("plus"),badgeLabel = "►", badgeColor = "olive"),
      br(),br(),
      HTML("<p align='center'>Nous étudions les courbes de<br>
           température et de concentration<br>
           d'ozone au cours de l'année à partir<br>
           de données collectées dans 426<br>
           stations meteo aux États-Unis</p>"),
      br(),br(),
      HTML('<p align="center">Paramètres optimaux :</p>'),
      HTML('<table align="center" border=1 width=50%>
            <tr>
               <td align="center">α</td>
               <td align="center">0.001</td>
            </tr>
            <tr>
               <td align="center">h</td>
               <td align="center">0.12</td>
            </tr>
               </table>'),
      br(),br(),br(),
      div("Auteurs :"),
      div("Mathis CORDIER"),
      div("Amandine PEILLON")
    )
  ),
  dashboardBody(         
    tags$img(
      src = "https://www.tarmac-interim.com/wp-content/uploads/2017/06/background-gris.jpg",
      style = 'position:absolute',
      class = "topimg",
      height = "100%",
      width = 1400), 
    tags$style(".topimg {
                   margin-left:-15px;
                   margin-right:-15px;
                   margin-top:-15px;}"),
    tags$style(".leaflet-container {
                   cursor: auto !important;}"),
    tabItems(
      tabItem(tabName="visu",
                box(width=3,htmlOutput('consignes'),sliderInput("l","Valeur du coefficient d'extrapolation (β/α)",min=0.1,max=0.9,0.5,step=0.1),
                  textOutput('long'),textOutput('lat')),
                box(width=9,title="Carte et analyse de la localisation",background="light-blue",leafletOutput('map')),
                box(plotOutput('graphe')),
                box(plotOutput('influ'))
              ),
      tabItem(tabName="moy",
                box(title="Tableau des moyennes observées",width=4,background="light-blue",
                  tableOutput("table"),
                  HTML("<p align='center'>...</p>")),
                box(background="light-blue",width=8,title="Carte des moyennes observées",
                  box(title="Moyenne des températures par station sur l'année",width="100%",
                      leafletOutput('moy_tp')),
                  box(title="Moyenne des valeurs d'ozone par station sur l'année",width="100%",
                      leafletOutput('moy_oz'))
                )
              ),
      tabItem(tabName="stat",
            box(width=6,
              box(title="Sélection de la station",width="100%",
                  sliderInput("n_stat","",min=1,max=426,100,pre="Station ")),
              box(title="Moyenne pour cette station sur l'année",background="light-blue",width="100%",
                  tableOutput("tab_stat"))
            ),
            box(title="Localisation",width=6,leafletOutput('loc')
            ),
            fluidRow(
              box(title="Évolution des températures pour la station sur l'année",width=6,background="light-blue",
                  plotOutput('evol_stat_tp')),
              box(title="Évolution des concentrations d'ozone pour la station sur l'année",width=6,background="light-blue",
                  plotOutput('evol_stat_oz'))
                )
              ),
      tabItem(tabName="jour",
              fluidRow(
                box(title="Sélection du jour",width=12,
                    sliderInput("n_jour","",min=1,max=365,100,pre="Jour ",step=5,animate=animationOptions(interval = 1000)))
              ),
              fluidRow(
                box(title="Carte des températures",background="light-blue",
                    leafletOutput('evol_day_tp')),
                box(title="Carte des concentrations d'ozone",background="light-blue",
                    leafletOutput('evol_day_oz'))
              )
      )
    )
  )
)

server = function(input, output, session) {

  # STATIONS
  output$loc <- renderLeaflet({
    localise(DF$lng[input$n_stat],DF$lat[input$n_stat])
  })
  
  output$tab_stat <- renderTable({
    dataplot[input$n_stat,]
  })
  
  output$evol_stat_tp <- renderPlot({
    evol_station_tp(input$n_stat)
  })
  
  output$evol_stat_oz <- renderPlot({
    evol_station_oz(input$n_stat)
  })
  
  # JOURS
  output$evol_day_tp <- renderLeaflet({
    evol_jour_tp(input$n_jour)
  })
  
  output$evol_day_oz <- renderLeaflet({
    evol_jour_oz(input$n_jour)
  })
  
  # SYNTHESE
  output$table <- renderTable(head(dataplot,20))
  
  output$moy_tp <- renderLeaflet({
    moy_map_tp()
  })
  
  output$moy_oz <- renderLeaflet({
    moy_map_oz()
  })
  
  # ESTIMATIONS
  
  output$map <- renderLeaflet({
    carte()
  })
  
  output$consignes <- renderUI(HTML("<font color='#3D508A'><b>CLIQUEZ SUR LA CARTE POUR LANCER</b></font><br><br><p align='center'>
                                    <img src='https://cdn0.iconfinder.com/data/icons/office-237/64/graph-graphic-stat-business-chart-512.png' height=100></img>
                                    </p><br><br>"))
  
  observeEvent(input$map_click, {
    
    click <- input$map_click
    lieu <- c(click$lng,click$lat)
    
    leafletProxy('map') %>% clearMarkers() %>%
      addMarkers(lng = lieu[1],lat = lieu[2],
                              label = paste("Position :\n(",lieu[1],",",lieu[2],")",sep=""),
                              icon = list(iconUrl = "http://icons.iconarchive.com/icons/martz90/circle/512/maps-icon.png",iconSize = c(20,20)))
    
    output$long <- renderText(paste("Longitude  :  ",lieu[1]))
    output$lat <- renderText(paste("Latitude    :  ",lieu[2]))
    
    output$influ <- renderPlot(influ(lieu))
    
    output$graphe <- renderPlot({
      graphe_day(lieu,input$l)
    })
    
  })
}



shinyApp(ui,server)
