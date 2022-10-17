library(ggplot2)
library(data.table)

RawData <- fread("PreProcessedData_Final_new.csv", na.strings = "nil", check.names = TRUE)

RawData$Month <- lubridate::month(RawData$Date)
RawData$Hour <- as.numeric(substring(RawData$Time, 1, 2))
RawData$Minute <- as.numeric(substring(RawData$Time, 4, 5))
RawData$HOD <- RawData$Hour + RawData$Minute/60

RawData <- RawData[!is.na(RawData$Temperature),]
RawData <- RawData[!is.na(RawData$cloud),]

PigPres <- t(sapply(1:nrow(RawData), function(j) sapply(1:20, function(i)
  sum(RawData[j, 3:12]==paste0("Pig", i)))))
colnames(PigPres) <- paste0("Pig", 1:20)

RawData <- cbind(RawData, PigPres)

ggplot(RawData, aes(x = Temperature, y = Count)) + geom_smooth()
ggsave("TempFuggesEgyben.pdf", width = 16, height = 9)

ggplot(melt(RawData[ , c("Temperature", paste0("Pig", 1:20))], id.vars = "Temperature"),
       aes(x = Temperature, y = value, group = variable, color = variable)) +
  geom_smooth()
ggsave("TempFuggesEgyesevel.pdf", width = 16, height = 9)

goodpigs <- as.character(melt(RawData[ ,paste0("Pig", 1:20)])[
  , .(sum(value)), .(variable)][V1>200]$variable)

cl <- parallel::makeCluster(parallel::detectCores()-1)
parallel::clusterExport(cl, "RawData")

fits <- parallel::parLapply(cl, goodpigs, function(pig)
  mgcv::gam(as.formula(paste0(pig, " ~ s(Temperature) + s(Pressure) + s(Humidity) +
                              s(sr...solar.radiation) + s(rau24...rain.fall) +
                              s(f...wind) + s(fx...wind) + s(cloud) +
                              s(HOD, bs = 'cc', k = 12) + s(Month, k = 7)")),
            knots = list(HOD = seq(0, 24, length.out = 12)),
            data = RawData[RawData$Month!=6,], family = binomial(link = "logit")))

parallel::stopCluster(cl)

predgrid <- CJ(Temperature = median(RawData$Temperature),
               Pressure = median(RawData$Pressure),
               Humidity = median(RawData$Humidity),
               sr...solar.radiation = median(RawData$sr...solar.radiation),
               rau24...rain.fall = median(RawData$rau24...rain.fall),
               f...wind = median(RawData$f...wind),
               fx...wind = median(RawData$fx...wind),
               cloud = median(RawData$cloud),
               Month = c(8+4/30, 10+21/30, 1+9/30),
               HOD = seq(0, 24, length.out = 500))

for(i in 1:length(fits)) {
  res <- gratia::fitted_values(fits[[i]], predgrid)
  res <- merge(res, data.frame(Month = c(8+4/30, 10+21/30, 1+9/30),
                               Date = c("August 5", "October 22", "January 10")))
  ggplot(res, aes(x = HOD, y = fitted*100, ymin = lower*100, ymax = upper*100,
                  color = Date, fill = Date, group = Date)) +  geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    labs(x = "Hour of day", y = "Probability of appearance [%]",
         title = paste0(goodpigs[i], ", p",
                        preport(summary(fits[[i]])$s.table["s(HOD)", "p-value"],
                                Sign = TRUE)),
         subtitle = paste0("Adjusted to: ",
                           "Temperature: ", round(median(RawData$Temperature), 1), " oC, ",
                           "Pressure: ", round(median(RawData$Pressure), 1), " hPa, ",
                           "Humidity: ", round(median(RawData$Humidity), 1), "%, ",
                           "Solar radiation: ", round(median(RawData$sr...solar.radiation), 1), ",\n",
                           "Rainfall: ", round(median(RawData$rau24...rain.fall), 1), ", ",
                           "f-Wind: ", round(median(RawData$f...wind), 1), ", ",
                           "fx-Wind: ", round(median(RawData$fx...wind), 1), ", ",
                           "Cloud coverage: ", round(median(RawData$cloud), 1), ". "))
  ggsave(paste0("./Results/Result_2_", goodpigs[i], ".png"), width = 16/2, height = 9/2,
         units = "in", dpi = 300)
}

for(i in 1:length(fits)) {
  print(gratia::draw(fits[[i]], constant = coef(fits[[i]])[1], fun = function(x) plogis(x)*100) +
          patchwork::plot_annotation(title = goodpigs[i]))
  ggsave(paste0("./Results/Result_3_", goodpigs[i], ".png"), width = 16, height = 9,
         units = "in", dpi = 300)
}