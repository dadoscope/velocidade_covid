# Libraries
library(tidyverse)
library(distcrete)
library(epitrix)
library(incidence)
library(projections)
library(EpiEstim)
library(jsonlite)
library(RcppRoll)
library(forecast)
library(lubridate)
library(data.table)
library(tibble)


theme_set(theme_bw() +  theme_bw(base_size = 18) +  theme(panel.border = element_rect(size = 1.2)))
plasma_pal <- c("red", viridis::plasma(n = 6))


# Get R effective estimates at state and country levels 
estimate_effective_R <- function(data) {
  # Estimates R effective in state level and Brazil
  # of COVID-19 Cases.
  # Data must follow the following convention:
  #   colnames: "UF", "DataBanco", "CasosNovos" 
  
  list_of_r_estimates <- list()
  
  # Preprocess data
  data$date_f <- as.Date(data$DataBanco, format = '%Y-%m-%d')
  data$CasosNovos <- as.numeric(data$CasosNovos)
  
  data$CasosNovos <- ifelse(data$CasosNovos <= 0, 0, data$CasosNovos)
  
  data <- 
    data %>% 
    dplyr::select(UF, date_f, CasosNovos) %>% 
    ungroup()
  names(data) <- c("state", "date", "newCases")
  
  br <- 
    data %>% 
    group_by(date) %>%  
    summarize(newCases = sum(newCases))
  br$state <- "BR"
  br <- br[, c(3, 1, 2)]
  data <- rbind(data, br)
  
  state_vector <- unique(data$state)
  
  for (i in 1:length(state_vector)) {
    
    data_incidence_function_data <- 
      data %>%
      dplyr::filter(`state` == state_vector[i]) %>%
      dplyr::select(`date`, `newCases`) %>%
      uncount(newCases)
    
    data_incidence_function_data$date <- as.Date(data_incidence_function_data$date)
    seq_dates <- data.frame(date = seq.Date(from = as.Date(min(data_incidence_function_data$date)), 
                                            to = as.Date(max(data_incidence_function_data$date)),
                                            by = 'day'))
    
    data_incidence_function_data <- left_join(seq_dates, data_incidence_function_data, by = "date")
    data_incidence_object        <- incidence(data_incidence_function_data$date)
    
    data_res_parametric_si <-  estimate_R(data_incidence_object, 
                                          method = "parametric_si",
                                          config = make_config(
                                            list(mean_si = 7.5, std_si = 3.4)
                                          ))
    # Combine R effective and dates
    data_r_effective <- data_res_parametric_si$R
    data_dates <- data.frame(date = data_incidence_object$dates[8:length(data_incidence_object$dates)])
    
    # Create data frame
    data_r_effective <- bind_cols(data_dates, data_r_effective)
    data_r_effective$UF <- state_vector[i]
    # To know the last record
    is_last <- rep(FALSE,nrow(data_r_effective))
    is_last[length(is_last)] <- TRUE
    data_r_effective$is_last <- is_last
    
    list_of_r_estimates[[i]] <- data_r_effective
    
  }
  
  data_r_effective <- bind_rows(list_of_r_estimates)

  return(data_r_effective)
}



setwd("/Users/isiscosta/RScript/moving_average_covid_deaths")
link = "https://data.brasil.io/dataset/covid19/caso.csv.gz"
dat = fread(link, stringsAsFactors = FALSE)
regioes <- c("Norte","Nordeste","Sudeste","Sul","Centro-Oeste")

dat <- dat %>% 
  mutate(date = ymd(date)) %>%
  filter(place_type == "city") %>% 
  mutate(cod_regiao = trunc(city_ibge_code / 10),
         regiao = regioes[cod_regiao])

#filtrando cidades com mais de 1000 casos
cities_mais_1000_casos <- dat %>%
  filter(is_last) %>%
  filter(confirmed >= 1000) %>%
  select(city_ibge_code)
  
daily <- dat %>% 
  filter(city_ibge_code %in% cities_mais_1000_casos$city_ibge_code) %>%
  filter(city_ibge_code != 733) %>% #733: city = "Importados/Indefinidos", state = "PB
  group_by(city, state) %>%
  arrange(date) %>%
  mutate(confirmed = c(0,diff(confirmed))) %>%
  #mutate(confirmed = forecast::ma(confirmed,order=7)) %>% 
  select(state, date, confirmed)

daily <- daily %>% 
  ungroup() %>% 
  mutate(UF = paste(city,state,sep="_"), DataBanco = date, CasosNovos = confirmed) %>% 
  select(UF, DataBanco, CasosNovos)

all_rt <- data.frame()
for (uf in unique(daily$UF)){
  df <- daily %>% filter(UF == uf)
  all_rt <- rbind(all_rt,estimate_effective_R(df))
}

#number of cities
summ_sta <- all_rt %>% ungroup()%>% filter(UF != "BR") %>%
  filter(date == "2020-07-09") %>% 
  group_by(date) %>%
  summarise(cities = n(),
            verdes = sum(`Mean(R)` < 0.8)/cities*100,
            vermelhos = sum(`Mean(R)` > 1.2)/cities*100,
            laranjas = sum(`Mean(R)` >= 0.8 & `Mean(R)` <= 1.2)/cities*100,
            mediart = mean(`Mean(R)`),
            maxrt = max(`Mean(R)`),
            minrt = min(`Mean(R)`),
            minUF = UF[which(`Mean(R)`==minrt)],
            maxUF = UF[which(`Mean(R)`==maxrt)])
  
p1 <- all_rt %>% filter(UF == summ_sta$minUF) %>%
  ggplot(aes(x = date, y = `Mean(R)`)) +
  geom_line(col = "green3", cex = 1.5) +
  labs(title = "Taxa de expansão da epidemia",
       subtitle = summ_sta$minUF)

p2 <- all_rt %>% filter(UF == summ_sta$maxUF) %>%
  ggplot(aes(x = date, y = `Mean(R)`)) +
  geom_line(col = "red3", cex = 1.5) +
  labs(title = "Taxa de expansão da epidemia",
       subtitle = summ_sta$maxUF)

png("minimo_rt.png",width=3200,height=1800,res=300)
print(p1)
dev.off()

png("maximo_rt.png",width=3200,height=1800,res=300)
print(p2)
dev.off()



