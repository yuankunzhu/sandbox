housing <- read.csv('dataSets/landdata-states.csv')
hp2001Q1 <- subset(housing, Date == 2001.25)
hp2001Q1$pred.SC <- predict(lm(Structure.Cost ~ log(Land.Value), data = hp2001Q1))
ggplot(
  hp2001Q1,
  aes(
    y = Structure.Cost,
    x = log(Land.Value))
) + geom_point(aes(color = Home.Value)) + geom_smooth() + geom_line(aes(y = pred.SC))
p1 <- ggplot(hp2001Q1, aes(x = log(Land.Value), y = Structure.Cost))
p1 + geom_point(color = rainbow(nrow(hp2001Q1)), size = 20, alpha = 4/10)

dat <- read.csv("dataSets/EconomistData.csv")
# Exercise I
# 1. Create a scatter plot with CPI on the x axis and HDI on the y axis.
ggplot(dat, aes(x=CPI, y=HDI)) + geom_point() + geom_smooth(method = 'lm')
# 5. Map the size of the points to HDI.Rank
ggplot(dat, aes(x=CPI, y=HDI)) + geom_point(aes(color=Region, size=HDI.Rank), alpha=.6)

ggplot(housing, aes(x = State, y = Home.Value)) + geom_col()
ggplot(mpg, aes(class)) + geom_bar()

# Exercise II
# 2. Overlay a smoothing line on top of the scatter plot using geom_smooth.
ggplot(dat, aes(x=CPI, y=HDI)) + geom_point(aes(col=Region)) + geom_smooth()
# 3. Overlay a smoothing line on top of the scatter plot using geom_smooth, but use a linear model for the predictions
ggplot(dat, aes(x=CPI, y=HDI)) + geom_point() + geom_smooth(method = 'lm')
# 4. Overlay a smoothing line on top of the scatter plot using geom_line.
ggplot(dat, aes(x=CPI, y=HDI)) + geom_point() + geom_line(aes(y = predict(lm(HDI ~ CPI, data = dat))))

p3 <- ggplot(housing,
             aes(x = State,
                 y = Home.Price.Index)) + 
  theme(legend.position="top", axis.text=element_text(size = 6))
(p4 <- p3 + geom_point(aes(color = Date),
                       alpha = 0.5,
                       size = 1.5,
                       position = position_jitter(width = 0.1, height = 0)))

p4 + scale_x_discrete(name="State Abbreviation") +
  scale_y_discrete(name="HPI") +
  scale_color_continuous(name="Date",
                         breaks = c(1976, 1994, 2013),
                         labels = c("'76", "'94", "'13"),
                         low = "blue", high = "red")
