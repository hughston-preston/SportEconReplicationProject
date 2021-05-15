library(tidyverse)
library(baseballr)
library(dplyr)
library(dbplyr)
library(ggplot2)
library(stargazer)
#library(R2jags)
#library(R2WinBUGS)
setwd("~/Desktop/Year 3/Semester 2/SAL 366/Data")

##2015
april2015_1 <- scrape_statcast_savant_pitcher_all("2015-04-05", "2015-04-07")
april2015_2 <- scrape_statcast_savant_pitcher_all("2015-04-08", "2015-04-14")
april2015_3 <- scrape_statcast_savant_pitcher_all("2015-04-15", "2015-04-22")
april2015_4 <- scrape_statcast_savant_pitcher_all("2015-04-23", "2015-04-30")

April2015 <- rbind(april2015_1, april2015_2, april2015_3, april2015_4)

may2015_1 <- scrape_statcast_savant_pitcher_all("2015-05-01", "2015-05-07")
may2015_2 <- scrape_statcast_savant_pitcher_all("2015-05-08", "2015-05-14")
may2015_3 <- scrape_statcast_savant_pitcher_all("2015-05-15", "2015-05-22")
may2015_4 <- scrape_statcast_savant_pitcher_all("2015-05-23", "2015-05-31")

May2015 <- rbind(may2015_1, may2015_2, may2015_3, may2015_4)

june2015_1 <- scrape_statcast_savant_pitcher_all("2015-06-01", "2015-06-07")
june2015_2 <- scrape_statcast_savant_pitcher_all("2015-06-08", "2015-06-14")
june2015_3 <- scrape_statcast_savant_pitcher_all("2015-06-15", "2015-06-22")
june2015_4 <- scrape_statcast_savant_pitcher_all("2015-06-23", "2015-06-30")

June2015 <- rbind(june2015_1, june2015_2, june2015_3, june2015_4)

july2015_1 <- scrape_statcast_savant_pitcher_all("2015-07-01", "2015-07-07")
july2015_2 <- scrape_statcast_savant_pitcher_all("2015-07-08", "2015-07-14")
july2015_3 <- scrape_statcast_savant_pitcher_all("2015-07-15", "2015-07-22")
july2015_4 <- scrape_statcast_savant_pitcher_all("2015-07-23", "2015-07-31")

July2015 <- rbind(july2015_1, july2015_2, july2015_3, july2015_4)

august2015_1 <- scrape_statcast_savant_pitcher_all("2015-08-01", "2015-08-07")
august2015_2 <- scrape_statcast_savant_pitcher_all("2015-08-08", "2015-08-14")
august2015_3 <- scrape_statcast_savant_pitcher_all("2015-08-15", "2015-08-22")
august2015_4 <- scrape_statcast_savant_pitcher_all("2015-08-23", "2015-08-31")

August2015 <- rbind(august2015_1, august2015_2, august2015_3, august2015_4)

september2015_1 <- scrape_statcast_savant_pitcher_all("2015-09-01", "2015-09-07")
september2015_2 <- scrape_statcast_savant_pitcher_all("2015-09-08", "2015-09-14")
september2015_3 <- scrape_statcast_savant_pitcher_all("2015-09-15", "2015-09-22")
september2015_4 <- scrape_statcast_savant_pitcher_all("2015-09-23", "2015-09-30")
october <- scrape_statcast_savant_pitcher_all("2015-10-01", "2015-10-04")

September2015 <- rbind(september2015_1, september2015_2, september2015_3, september2015_4, october)

x2015data <- rbind(April2015, May2015, July2015, September2015, June2015, August2015)

##2016
april2016_1 <- scrape_statcast_savant_pitcher_all("2016-04-03", "2016-04-07")
april2016_2 <- scrape_statcast_savant_pitcher_all("2016-04-08", "2016-04-14")
april2016_3 <- scrape_statcast_savant_pitcher_all("2016-04-15", "2016-04-22")
april2016_4 <- scrape_statcast_savant_pitcher_all("2016-04-23", "2016-04-30")

April2016 <- rbind(april2016_1, april2016_2, april2016_3, april2016_4)

may2016_1 <- scrape_statcast_savant_pitcher_all("2016-05-01", "2016-05-07")
may2016_2 <- scrape_statcast_savant_pitcher_all("2016-05-08", "2016-05-14")
may2016_3 <- scrape_statcast_savant_pitcher_all("2016-05-15", "2016-05-22")
may2016_4 <- scrape_statcast_savant_pitcher_all("2016-05-23", "2016-05-31")

May2016 <- rbind(may2016_1, may2016_2, may2016_3, may2016_4)

june2016_1 <- scrape_statcast_savant_pitcher_all("2016-06-01", "2016-06-07")
june2016_2 <- scrape_statcast_savant_pitcher_all("2016-06-08", "2016-06-14")
june2016_3 <- scrape_statcast_savant_pitcher_all("2016-06-15", "2016-06-22")
june2016_4 <- scrape_statcast_savant_pitcher_all("2016-06-23", "2016-06-30")

June2016 <- rbind(june2016_1, june2016_2, june2016_3, june2016_4)

july2016_1 <- scrape_statcast_savant_pitcher_all("2016-07-01", "2016-07-07")
july2016_2 <- scrape_statcast_savant_pitcher_all("2016-07-08", "2016-07-14")
july2016_3 <- scrape_statcast_savant_pitcher_all("2016-07-15", "2016-07-22")
july2016_4 <- scrape_statcast_savant_pitcher_all("2016-07-23", "2016-07-31")

July2016 <- rbind(july2016_1, july2016_2, july2016_3, july2016_4)

august2016_1 <- scrape_statcast_savant_pitcher_all("2016-08-01", "2016-08-07")
august2016_2 <- scrape_statcast_savant_pitcher_all("2016-08-08", "2016-08-14")
august2016_3 <- scrape_statcast_savant_pitcher_all("2016-08-15", "2016-08-22")
august2016_4 <- scrape_statcast_savant_pitcher_all("2016-08-23", "2016-08-31")

August2016 <- rbind(august2016_1, august2016_2, august2016_3, august2016_4)

september2016_1 <- scrape_statcast_savant_pitcher_all("2016-09-01", "2016-09-07")
september2016_2 <- scrape_statcast_savant_pitcher_all("2016-09-08", "2016-09-14")
september2016_3 <- scrape_statcast_savant_pitcher_all("2016-09-15", "2016-09-22")
september2016_4 <- scrape_statcast_savant_pitcher_all("2016-09-23", "2016-09-30")
october <- scrape_statcast_savant_pitcher_all("2016-10-01", "2016-10-02")

September2016 <- rbind(september2016_1, september2016_2, september2016_3, september2016_4, october)

x2016data <- rbind(April2016, May2016, July2016, September2016, June2016, August2016)

##2017

april2017_1 <- scrape_statcast_savant_pitcher_all("2017-04-02", "2017-04-07")
april2017_2 <- scrape_statcast_savant_pitcher_all("2017-04-08", "2017-04-14")
april2017_3 <- scrape_statcast_savant_pitcher_all("2017-04-15", "2017-04-22")
april2017_4 <- scrape_statcast_savant_pitcher_all("2017-04-23", "2017-04-30")

April2017 <- rbind(april2017_1, april2017_2, april2017_3, april2017_4)

may2017_1 <- scrape_statcast_savant_pitcher_all("2017-05-01", "2017-05-07")
may2017_2 <- scrape_statcast_savant_pitcher_all("2017-05-08", "2017-05-14")
may2017_3 <- scrape_statcast_savant_pitcher_all("2017-05-15", "2017-05-22")
may2017_4 <- scrape_statcast_savant_pitcher_all("2017-05-23", "2017-05-31")

May2017 <- rbind(may2017_1, may2017_2, may2017_3, may2017_4)

june2017_1 <- scrape_statcast_savant_pitcher_all("2017-06-01", "2017-06-07")
june2017_2 <- scrape_statcast_savant_pitcher_all("2017-06-08", "2017-06-14")
june2017_3 <- scrape_statcast_savant_pitcher_all("2017-06-15", "2017-06-22")
june2017_4 <- scrape_statcast_savant_pitcher_all("2017-06-23", "2017-06-30")

June2017 <- rbind(june2017_1, june2017_2, june2017_3, june2017_4)

july2017_1 <- scrape_statcast_savant_pitcher_all("2017-07-01", "2017-07-07")
july2017_2 <- scrape_statcast_savant_pitcher_all("2017-07-08", "2017-07-14")
july2017_3 <- scrape_statcast_savant_pitcher_all("2017-07-15", "2017-07-22")
july2017_4 <- scrape_statcast_savant_pitcher_all("2017-07-23", "2017-07-31")

July2017 <- rbind(july2017_1, july2017_2, july2017_3, july2017_4)

august2017_1 <- scrape_statcast_savant_pitcher_all("2017-08-01", "2017-08-07")
august2017_2 <- scrape_statcast_savant_pitcher_all("2017-08-08", "2017-08-14")
august2017_3 <- scrape_statcast_savant_pitcher_all("2017-08-15", "2017-08-22")
august2017_4 <- scrape_statcast_savant_pitcher_all("2017-08-23", "2017-08-31")

August2017 <- rbind(august2017_1, august2017_2, august2017_3, august2017_4)

september2017_1 <- scrape_statcast_savant_pitcher_all("2017-09-01", "2017-09-07")
september2017_2 <- scrape_statcast_savant_pitcher_all("2017-09-08", "2017-09-14")
september2017_3 <- scrape_statcast_savant_pitcher_all("2017-09-15", "2017-09-22")
september2017_4 <- scrape_statcast_savant_pitcher_all("2017-09-23", "2017-10-01")

September2017 <- rbind(september2017_1, september2017_2, september2017_3, september2017_4)

x2017data <- rbind(April2017, May2017, July2017, September2017, June2017, August2017)

replication_data <- rbind(x2015data, x2016data, x2017data)
write_csv(replication_data, "replication_data.csv")
replication_data <- read_csv("replication_data.csv")

replication_data$game_date <- as.Date(replication_data$game_date, format = "%m/%d/%y")
replication_data$year <- format(replication_data$game_date, "%Y")
replication_data1 <- replication_data %>% filter(!is.na(release_speed))

omit <- replication_data1 %>% 
  filter(inning == "1") %>% 
  group_by(pitcher...8) %>% 
  summarise(player_name = player_name) %>% 
  unique()
notstarters <- replication_data1[!replication_data1$pitcher...8 %in% omit$pitcher...8, ]

#Piecewise Linear Function Knots: 25th, 50th, 75th, 90th percentiles
percentiles <- c(0.25, 0.5, 0.75, 0.9)
knots <- notstarters %>% group_by(pitcher...8, player_name, game_date) %>% summarise(pitches = n())
quantile(knots$pitches, percentiles)

pitchcounts <- notstarters %>% 
  group_by(pitcher...8, year) %>% 
  summarise(player_name = player_name, 
            pitches = n()) %>% 
  unique() %>% 
  filter(pitches > 200)

check324Sample <- pitchcounts %>% 
  group_by(pitcher...8) %>% 
  summarise(player_name = player_name) %>% 
  unique() 
#402 pitchers is close to their claimed sample of 324 pitchers. Final trim may come when we filter sample to only include outings within 10 days of each other next

#10% of pitches thrown, which is hardest?
total_pitchcounts <- pitchcounts %>% 
  group_by(pitcher...8) %>% 
  summarise(player_name = player_name, 
            TotalPitches = sum(pitches)) %>% 
  unique()

omit1 <- notstarters %>% 
  #filter(!is.na(release_speed)) %>% 
  group_by(pitcher...8, pitch_type) %>% 
  summarise(player_name = player_name, 
            velo = mean(release_speed), 
            pitches = n()) %>% 
  unique() %>% 
  filter(pitch_type == "FF"| 
           pitch_type == "FT" | 
           pitch_type == "FC" | 
           pitch_type == "SL" | 
           pitch_type == "FS" | 
           pitch_type == "SI")

PitchesThrown <- total_pitchcounts %>% 
  left_join(omit1, by = c("pitcher...8" = "pitcher...8", 
                          "player_name" = "player_name")) %>% 
  unique() %>% 
  mutate(PCT_Thrown = pitches / TotalPitches) %>% 
  filter(PCT_Thrown >= 0.1) %>%
  group_by(pitcher...8) %>% 
  top_n(1, velo)

#Games within 10 days of eachother 
Trim1Sample <- replication_data1 %>% 
  right_join(PitchesThrown, by = c("pitcher...8" = "pitcher...8", 
                                   "player_name" = "player_name", 
                                   "pitch_type" = "pitch_type")) %>%
  group_by(game_date, 
           pitcher...8, 
           player_name) %>% 
  summarise(FBs_Thrown = n(), 
            pitch_type = pitch_type) %>%
  unique() %>%
  filter(FBs_Thrown >= 3) %>%
  arrange(player_name, 
        game_date) %>%
  ungroup()

Trim1Sample$year <- format(Trim1Sample$game_date, "%Y")

Trim2Sample <- 
  Trim1Sample %>%
  group_by(pitcher...8, player_name, year) %>%
  mutate(lag.date = lag(game_date, n = 1, default = NA)) %>% #This is a hell of a bit of lag variable code I found. Saving for later for personal use
  mutate(DOR = difftime(game_date ,lag.date , units = c("days"))) %>%
  filter(DOR <= 10) %>%
  ungroup()

ReadyData <- Trim2Sample %>% 
  left_join(replication_data1, by = c("pitcher...8" = "pitcher...8", 
                                     "player_name" = "player_name", 
                                     "game_date" = "game_date", 
                                     "year" = "year", 
                                     "pitch_type" = "pitch_type")) %>% 
  unique()

#Piecewise Linear Function Knots: 25th, 50th, 75th, 90th percentiles
percentiles <- c(0.25, 0.5, 0.75, 0.9)
knots1 <- ReadyData %>% group_by(pitcher...8, player_name, game_date) %>% summarise(pitches = n())
quantile(knots$pitches, percentiles)

#after testing the knots for both the all relievers and our sample relievers we have found that both showed percentile values of 25th:10, 50th:15, 75th:20, and 90th:27
#our 25th and 75th percentiles are 1 off of those claimed in the research paper
kSecondsInDay = 24 * 3600
pgame <- ReadyData %>%
  group_by(game_date, 
           pitcher...8, 
           player_name) %>%
  summarise(pitches = n(), 
            DOR = DOR, 
            Velo = mean(release_speed), 
            SDVelo = sd(release_speed), 
            year = year) %>%
  mutate(num.date = as.numeric(strptime(game_date, 
                                        format = "%Y-%m-%d")) / kSecondsInDay) %>%
  unique() %>%
  arrange(pitcher...8, 
          num.date) %>%
  group_by(pitcher...8, 
           player_name) %>%
  mutate(PreviousVelo = lag(Velo)) %>%
  mutate(Velo_Difference = Velo - PreviousVelo)

pgameJoin <- ReadyData %>% 
  group_by(pitcher...8, year) %>% 
  summarize(SeasonAvgVelo = mean(release_speed))

pgame <- pgame %>% 
  left_join(pgameJoin, by = c("pitcher...8" = "pitcher...8", "year" = "year")) %>%
  mutate(SeasonVeloDiff = Velo - SeasonAvgVelo)

totalpitches <- replication_data %>% 
  group_by(pitcher...8, player_name, game_date) %>% 
  summarise(PitchCount = n())

pgame <- pgame %>% left_join(totalpitches, by = c("pitcher...8" = "pitcher...8", 
                                                         "game_date" = "game_date", 
                                                         "player_name" = "player_name"))

##########
Pitches_Thrown = matrix(0, nrow = nrow(pgame), ncol = 6)
colnames(Pitches_Thrown) = paste0("Pitches_Thrown_Day",1:6)
for(i in 1:nrow(pgame)){
  for(j in 1:6){
    index = ifelse((i-j)>0, i-j, NA) 
    if(!is.na(index) &&
       (pgame$player_name[i] == pgame$player_name[(i-j)]) && 
       (pgame$num.date[i] - pgame$num.date[(i-j)] <= 6)){ 
      d = (pgame$num.date[i] - pgame$num.date[i-j]) 
      Pitches_Thrown[i,d] = pgame$PitchCount[i-j]
    }
  }
}

games <- as.data.frame(cbind(as.matrix(pgame),Pitches_Thrown)) %>%
  select(-num.date) %>%
  ungroup() %>%
  mutate(player_name = factor(player_name)) %>%
  unique()

#Linear Mixed Effects Model
games$pitcher...8 <- as.character(games$pitcher...8)
games$Velo_Difference <- as.numeric(as.character(games$Velo_Difference))
games$SeasonVeloDiff <- as.numeric(as.character(games$SeasonVeloDiff))
games$PitchCount <- as.numeric(as.character(games$PitchCount))
games$Pitches_Thrown_Day1 <- as.numeric(as.character(games$Pitches_Thrown_Day1))
games$Pitches_Thrown_Day2 <- as.numeric(as.character(games$Pitches_Thrown_Day2))
games$Pitches_Thrown_Day3 <- as.numeric(as.character(games$Pitches_Thrown_Day3))
games$Pitches_Thrown_Day4 <- as.numeric(as.character(games$Pitches_Thrown_Day4))
games$Pitches_Thrown_Day5 <- as.numeric(as.character(games$Pitches_Thrown_Day5))
games$Pitches_Thrown_Day6 <- as.numeric(as.character(games$Pitches_Thrown_Day6))

games <- games %>% 
  mutate(TotalPitchesLast6Days = Pitches_Thrown_Day1 + 
           Pitches_Thrown_Day2 + 
           Pitches_Thrown_Day3 + 
           Pitches_Thrown_Day4 + 
           Pitches_Thrown_Day5 + 
           Pitches_Thrown_Day6)

library(lme4)
library(jtools)
library(lmerTest)
lmer1 <- lmer(SeasonVeloDiff ~ (Pitches_Thrown_Day1|pitcher...8) + 
                (Pitches_Thrown_Day2|pitcher...8) + 
                (Pitches_Thrown_Day3|pitcher...8) + 
                (Pitches_Thrown_Day4|pitcher...8) + 
                (Pitches_Thrown_Day5|pitcher...8) + 
                (Pitches_Thrown_Day6|pitcher...8) +
                Pitches_Thrown_Day1 + 
                Pitches_Thrown_Day2 + 
                Pitches_Thrown_Day3 + 
                Pitches_Thrown_Day4 + 
                Pitches_Thrown_Day5 + 
                Pitches_Thrown_Day6, 
              data = games)
#Summary
summary(lmer1)
summ(lmer1)
ranova(lmer1)
ranef(lmer1)
PitcherRandomEffects1 <- ranef(lmer1)[[1]] %>% 
  unique()
plot(ranef(lmer1))
qqnorm(resid(lmer1))
qqline(resid(lmer1))


#Getting Betas for Spline
PitcherRandomEffects1$pitcher...8 = rownames(PitcherRandomEffects1)
names(PitcherRandomEffects1)[1] <- "Day1Intercept"
names(PitcherRandomEffects1)[2] <- "Day1Coef"
names(PitcherRandomEffects1)[3] <- "Day2Intercept"
names(PitcherRandomEffects1)[4] <- "Day2Coef"
names(PitcherRandomEffects1)[5] <- "Day3Intercept"
names(PitcherRandomEffects1)[6] <- "Day3Coef"
names(PitcherRandomEffects1)[7] <- "Day4Intercept"
names(PitcherRandomEffects1)[8] <- "Day4Coef"
names(PitcherRandomEffects1)[9] <- "Day5Intercept"
names(PitcherRandomEffects1)[10] <- "Day5Coef"
names(PitcherRandomEffects1)[11] <- "Day6Intercept"
names(PitcherRandomEffects1)[12] <- "Day6Coef"

percentiles <- c(0.25, 0.5, 0.75, 0.9)
knots1 <- ReadyData %>% 
  group_by(pitcher...8, player_name, game_date) %>% 
  summarise(PitchCount = n())
quantile(knots1$PitchCount, percentiles)


Games <- games %>% 
  right_join(PitcherRandomEffects1, by = c("pitcher...8" = "pitcher...8"))

#Betas
b1 <- Games %>% filter(quantile(PitchCount, 0.25) >= PitchCount)
b2 <- Games %>% filter(quantile(PitchCount, 0.25) < PitchCount & quantile(PitchCount, 0.5) >= PitchCount)
b3 <- Games %>% filter(quantile(PitchCount, 0.5) < PitchCount & quantile(PitchCount, 0.75) >= PitchCount)
b4 <- Games %>% filter(quantile(PitchCount, 0.75) < PitchCount & quantile(PitchCount, 0.9) >= PitchCount)
b5 <- Games %>% filter(quantile(PitchCount, 0.9) < PitchCount)
sampscale <- 24024 * (5/6)
B1 <- ((((nrow(b1)/nrow(Games)) * sampscale)) / nrow(Games))
B2 <- ((((nrow(b2)/nrow(Games)) * sampscale)) / nrow(Games)) + B1
B3 <- ((((nrow(b3)/nrow(Games)) * sampscale)) / nrow(Games)) + B2
B4 <- ((((nrow(b4)/nrow(Games)) * sampscale)) / nrow(Games)) + B3
B5 <- ((((nrow(b5)/nrow(Games)) * sampscale)) / nrow(Games)) + B4

Games <- Games %>%
  mutate(Dose = ifelse(PitchCount <= 10, B1 * PitchCount, 
                       ifelse(PitchCount > 10 & PitchCount <= 15, (10 * B1) + B2 * (PitchCount - 10), 
                              ifelse(PitchCount > 15 & PitchCount <= 20, (10 * B1) + (5 * B2) + B3 * (PitchCount - 15),
                                     ifelse(PitchCount > 20 & PitchCount <= 27, (10 * B1) + (5 * B2) + (5 * B3) + B4 * (PitchCount - 20),
                                            ifelse(PitchCount > 27, (10 * B1) + (5 * B2) + (5 * B3) + (7 * B4) + B5 * (PitchCount - 27),
                                                   NA))))),
         PC_Group = ifelse(PitchCount <= 10, 1, 
                           ifelse(PitchCount > 10 & PitchCount <= 15, 2, 
                                  ifelse(PitchCount > 15 & PitchCount <= 20, 3,
                                         ifelse(PitchCount > 20 & PitchCount <= 27, 4,
                                                ifelse(PitchCount > 27, 5,
                                                       NA))))))



ggplot(Games, aes(x = PitchCount, y = Dose)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)
Games <- Games %>% mutate(num.date = as.numeric(strptime(game_date, 
                                      format = "%Y-%m-%d")) / kSecondsInDay)
PitchingDosage = matrix(0, nrow = nrow(Games), ncol = 6)
colnames(PitchingDosage) = paste0("PitchingDosageDay",1:6)
for(i in 1:nrow(Games)){
  for(j in 1:6){
    index = ifelse((i-j)>0, i-j, NA) 
    if(!is.na(index) &&
       (Games$player_name[i] == Games$player_name[(i-j)]) && 
       (Games$num.date[i] - Games$num.date[(i-j)] <= 6)){ 
      d = (Games$num.date[i] - Games$num.date[i-j]) 
      PitchingDosage[i,d] = Games$Dose[i-j]
    }
  }
}

DosageGames <- as.data.frame(cbind(as.matrix(Games),PitchingDosage)) %>%
  select(-num.date) %>%
  ungroup() %>%
  mutate(player_name = factor(player_name)) %>%
  unique()  

DosageGames <- as.data.frame(cbind(as.matrix(DosageGames))) %>%
  #select(-num.date) %>%
  ungroup() %>%
  mutate(player_name = factor(player_name)) %>%
  unique()  

DosageGames$PitchingDosageDay1 <- as.numeric(as.character(DosageGames$PitchingDosageDay1))
DosageGames$PitchingDosageDay2 <- as.numeric(as.character(DosageGames$PitchingDosageDay2))
DosageGames$PitchingDosageDay3 <- as.numeric(as.character(DosageGames$PitchingDosageDay3))
DosageGames$PitchingDosageDay4 <- as.numeric(as.character(DosageGames$PitchingDosageDay4))
DosageGames$PitchingDosageDay5 <- as.numeric(as.character(DosageGames$PitchingDosageDay5))
DosageGames$PitchingDosageDay6 <- as.numeric(as.character(DosageGames$PitchingDosageDay6))
DosageGames$SeasonVeloDiff <- as.numeric(as.character(DosageGames$SeasonVeloDiff))
DosageGames$Velo_Difference <- as.numeric(as.character(DosageGames$Velo_Difference))


DosageGames <- DosageGames %>% 
  mutate(TotalDose = PitchingDosageDay1 + 
           PitchingDosageDay2 + 
           PitchingDosageDay3 + 
           PitchingDosageDay4 + 
           PitchingDosageDay5 + 
           PitchingDosageDay6)
summary(lm(SeasonVeloDiff ~ TotalDose, DosageGames))

ggplot(DosageGames, aes(x = TotalDose, y = SeasonVeloDiff)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)

#proceeding with the Assumption that the elimination constant ePhi is 0.55 because that was found in the paper. I have not found a clear method of inferring this through a frequentist model
ePhi <- 0.55
Vjt <- lmer(SeasonVeloDiff ~ (PitchingDosageDay1|pitcher...8) + 
              (PitchingDosageDay2|pitcher...8) + 
              (PitchingDosageDay3|pitcher...8) + 
              (PitchingDosageDay4|pitcher...8) + 
              (PitchingDosageDay5|pitcher...8) + 
              (PitchingDosageDay6|pitcher...8) +
              PitchingDosageDay1 +
              PitchingDosageDay2 +
              PitchingDosageDay3 +
              PitchingDosageDay4 +
              PitchingDosageDay5 +
              PitchingDosageDay6,
             data = DosageGames)
summary(Vjt)
Vjt.df <- as.data.frame(ranef(Vjt)[[1]])
Vjt.df$pitcher...8 = rownames(Vjt.df)
names(Vjt.df)[1] <- "DosageDay1Intercept"
names(Vjt.df)[2] <- "DosageDay1Coef"
names(Vjt.df)[3] <- "DosageDay2Intercept"
names(Vjt.df)[4] <- "DosageDay2Coef"
names(Vjt.df)[5] <- "DosageDay3Intercept"
names(Vjt.df)[6] <- "DosageDay3Coef"
names(Vjt.df)[7] <- "DosageDay4Intercept"
names(Vjt.df)[8] <- "DosageDay4Coef"
names(Vjt.df)[9] <- "DosageDay5Intercept"
names(Vjt.df)[10] <- "DosageDay5Coef"
names(Vjt.df)[11] <- "DosageDay6Intercept"
names(Vjt.df)[12] <- "DosageDay6Coef"

fix.vjt <- as.data.frame(fixef(Vjt))
Vjt.df$DosageDay1Intercept <- Vjt.df$DosageDay1Intercept + fix.vjt[1,]
Vjt.df$DosageDay2Intercept <- Vjt.df$DosageDay2Intercept + fix.vjt[1,]
Vjt.df$DosageDay3Intercept <- Vjt.df$DosageDay3Intercept + fix.vjt[1,]
Vjt.df$DosageDay4Intercept <- Vjt.df$DosageDay4Intercept + fix.vjt[1,]
Vjt.df$DosageDay5Intercept <- Vjt.df$DosageDay5Intercept + fix.vjt[1,]
Vjt.df$DosageDay6Intercept <- Vjt.df$DosageDay6Intercept + fix.vjt[1,]
Vjt.df$DosageDay1Coef <- Vjt.df$DosageDay1Coef + fix.vjt[2,]
Vjt.df$DosageDay2Coef <- Vjt.df$DosageDay2Coef + fix.vjt[3,]
Vjt.df$DosageDay3Coef <- Vjt.df$DosageDay3Coef + fix.vjt[4,]
Vjt.df$DosageDay4Coef <- Vjt.df$DosageDay4Coef + fix.vjt[5,]
Vjt.df$DosageDay5Coef <- Vjt.df$DosageDay5Coef + fix.vjt[6,]
Vjt.df$DosageDay6Coef <- Vjt.df$DosageDay6Coef + fix.vjt[7,]


DosageGames <- DosageGames %>% 
  right_join(Vjt.df, 
            by = c("pitcher...8" = "pitcher...8")) 

DosageGames1 <- DosageGames %>% mutate(concentration = (0.45 * 
                                                         (PitchingDosageDay1 + 
                                                            (0.45 * 
                                                               (PitchingDosageDay2 + 
                                                                  (0.45 * 
                                                                     (PitchingDosageDay3 + 
                                                                        (0.45 * 
                                                                           (PitchingDosageDay4 + 
                                                                              (0.45 * 
                                                                                 (PitchingDosageDay5 + 
                                                                                    (0.45 * PitchingDosageDay6))))))))))))

#Trying first with Velo Difference from season average Velo
ConcentrationEffect <- lmer(SeasonVeloDiff ~ (concentration|pitcher...8) + 
                              concentration, 
                            DosageGames)
#Summary
summary(ConcentrationEffect)
summ(ConcentrationEffect)
ranova(ConcentrationEffect)
ranef(ConcentrationEffect)
fixef(ConcentrationEffect)
plot(ranef(ConcentrationEffect))
qqnorm(resid(ConcentrationEffect))
qqline(resid(ConcentrationEffect))
class(ConcentrationEffect) <- "lmerMod"
stargazer(ConcentrationEffect, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
ran.CE.df <- as.data.frame(ranef(ConcentrationEffect)[[1]])
fix.CE.df <- as.data.frame(fixef(ConcentrationEffect))
ran.CE.df$pitcher...8 <- rownames(ran.CE.df)
names(ran.CE.df)[1] <- "ConcentrationEffectIntercept"
ran.CE.df$ConcentrationEffectIntercept <- ran.CE.df$ConcentrationEffectIntercept + fix.CE.df[1,]
names(ran.CE.df)[2] <- "ConcentrationEffectCoef"
ran.CE.df$ConcentrationEffectCoef <- ran.CE.df$ConcentrationEffectCoef + fix.CE.df[2,]



CE.df2 <- ran.CE.df %>% left_join(DosageGames, 
                             by = c("pitcher...8" = "pitcher...8")) %>% 
  mutate(xVeloLoss_byCE = (concentration * ConcentrationEffectCoef) + ConcentrationEffectIntercept) %>%
  select(concentration, xVeloLoss_byCE, ConcentrationEffectIntercept:DosageDay6Coef) %>%
  group_by(pitcher...8, player_name) %>%
  summarise(#concentration = concentration, 
            #xVeloLoss_byCE = xVeloLoss_byCE, 
            ConcentrationEffectCoef = ConcentrationEffectCoef,
            ConcentrationEffectIntercept = ConcentrationEffectIntercept) %>% unique()
  arrange(xVeloLoss_byCE) 

CE.df1 <- CE.df %>%
  group_by(pitcher...8, player_name) %>% 
  summarise(xVeloLoss_byCESUM = sum(xVeloLoss_byCE)) %>%
  unique()

#as briefly mentioned in paper, I don't know if normal distribution of velocity centered at 93mph is a responsible conclusion. this graph seem slightly left skewed and centered closer to 94mph
ggplot(CE.df, aes(x=Velo)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white") + 
  geom_density(alpha=.2, fill="#FF6666")
library(lattice)
dotplot(head(CE.df$xVeloLoss_byCE, 10),labels=row.names(CE.df$player_name),cex=.7,
         main="xVeloLoss Per Player", 
         xlab="xVeloLoss")

Plot_Head <- CE.df %>% 
  filter(ConcentrationEffectCoef > -0.01) %>% 
  filter(!is.na(Velo_Difference))
Plot_Tail <- CE.df %>% 
  filter(ConcentrationEffectCoef < -0.09) %>% 
  filter(!is.na(Velo_Difference))

ggplot(Plot_Head, aes(xVeloLoss_byCE, player_name)) + 
  geom_point()

ggplot(
  Plot_Head, 
  aes(x = `xVeloLoss_byCE`, y = `player_name`)
) +
  geom_density_ridges_gradient(
    aes(fill = year), scale = 1, size = 0.3
  ) +
  labs(title = 'Least Expected Velo Lost') 

ggplot(
  Plot_Tail, 
  aes(x = `xVeloLoss_byCE`, y = `player_name`)
) +
  geom_density_ridges_gradient(
    aes(fill = year), scale = 1, size = 0.3
  ) +
  labs(title = 'Most Expected Velo Lost') 

#Trying it with game to game Velo differences
ConcentrationEffect2 <- lmer(Velo_Difference ~ (concentration|pitcher...8), DosageGames)
summary(ConcentrationEffect2)
CE.df2 <- as.data.frame(ranef(ConcentrationEffect2)[[1]])
CE.df2$pitcher...8 <- rownames(CE.df2)
names(CE.df2)[1] <- "ConcentrationEffectIntercept2"
names(CE.df2)[2] <- "ConcentrationEffectCoef2"

CE.df2 <- CE.df2 %>% left_join(DosageGames, 
                             by = c("pitcher...8" = "pitcher...8")) %>% 
  group_by(pitcher...8, player_name) %>% 
  summarise(ConcentrationEffectIntercept2 = ConcentrationEffectIntercept2,
            ConcentrationEffectCoef2 = ConcentrationEffectCoef2) %>%
  unique()

#based off of inspection of the two data frames it seems that our velo difference effects from season average had most overlap with Paper's findings in terms of best and worst fatigue performers so we will continue with that basis
VeloFatigueEffectImplimentation <- DosageGames %>% 
  left_join(CE.df, by = c("pitcher...8" = "pitcher...8", 
                          "player_name" = "player_name")) %>%
  mutate(xVeloChange = ConcentrationEffectIntercept + (ConcentrationEffectCoef * concentration)) %>%
  arrange(-xVeloChange)











Games <- games %>% right_join(PitcherRandomEffects, 
                             by = c("pitcher...8" = "pitcher...8")) %>%
  mutate(Dose = ifelse(PitchCount <= 10, Day1Coef * PitchCount, 
                       ifelse(PitchCount > 10 & PitchCount <= 15, (10 * Day1Coef) + Day2Coef * (PitchCount - 10), 
                              ifelse(PitchCount > 15 & PitchCount <= 20, (10 * Day1Coef) + (4*Day2Coef) + Day3Coef * (PitchCount - 15),
                                     ifelse(PitchCount > 20 & PitchCount <= 27, (10 * Day1Coef) + (5 * Day2Coef) + (5 * Day3Coef) + Day4Coef * (PitchCount - 20),
                                            ifelse(PitchCount > 27, (10 * Day1Coef) + (5 * Day2Coef) + (5 * Day3Coef) + (7 * Day4Coef) + Day5Coef * (PitchCount - 20),
                                                   NA))))))


library(quantreg)
CE.df$SeasonAvgVelo <- as.numeric(CE.df$SeasonAvgVelo)
Velo_QuantReg <- rq(SeasonAvgVelo ~ ConcentrationEffectCoef, tau = c(0.25, 0.50, 0.75, 0.9), CE.df)
summary(Velo_QuantReg, se = "boot")
class(Velo_QuantReg) <- "rqsMod"
stargazer(Velo_QuantReg) #, type="html",out="QuantReg1.html")
stargazer(Velo_QuantReg, rq.se = "iid", type = "text", title="Regression Results", initial.zero = F,single.row=TRUE)
stargazer(Velo_QuantReg, type = "text", title="Regression Results", initial.zero = F,single.row=TRUE)


#Output
class(lmer1) <- "lmerMod"
stargazer(lmer1, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
stargazer(lmer1,type="html",out="LinearMix1.html")

class(Vjt) <- "lmerMod"
stargazer(Vjt, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
stargazer(Vjt,type="html",out="LinearMix2.html")

class(ConcentrationEffect) <- "lmerMod"
stargazer(ConcentrationEffect, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
stargazer(ConcentrationEffect,type="html",out="LinearMix3.html")
