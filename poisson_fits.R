setwd("~/DSI/capstone/ProteinFragmenter/")

library(tidyverse)
library(ggplot2)

# Try to fit a Poisson model to cap lengths

lens <- read.csv(file = "caps_lengths.csv")

pmodel <- glm(len_cnt ~ lengthCAP, family = "poisson", data = lens)
summary(pmodel)

fitted <- fitted(pmodel)
qplot(x = lens$lengthCAP, fitted)
ggplot() +
  geom_point(mapping = aes(x = lens$lengthCAP, y = lens$len_cnt)) +
  geom_path(mapping = aes(x = lens$lengthCAP, y = fitted)) +
  ggtitle("Helix Length vs Poisson Fitted Counts")

ggplot() +
  geom_point(mapping = aes(x = lens$lengthCAP, y = resid(pmodel))) +
  geom_hline(yintercept = 0) +
  ggtitle("Cap Length vs Poisson Deviance Residuals")

# It doesn't seem to be a good fit



# Try to fit a Poisson model to helix lengths - this should be a good fit

hlens <- read.csv(file = "helix_lengths.csv")

pmodel2 <- glm(len_cnt ~ lengthHELIX, family = "poisson", data = hlens)
summary(pmodel2)

fitted2 <- fitted(pmodel2)
qplot(x = hlens$lengthHELIX, fitted2)
ggplot() +
  geom_point(mapping = aes(x = hlens$lengthHELIX, y = hlens$len_cnt)) +
  geom_path(mapping = aes(x = hlens$lengthHELIX, y = fitted2)) +
  ggtitle("Helix Length vs Poisson Fitted Counts")

ggplot() +
  geom_point(mapping = aes(x = hlens$lengthHELIX, y = resid(pmodel2))) +
  geom_hline(yintercept = 0) +
  ggtitle("Helix Length vs Poisson Deviance Residuals")

# It doesn't seem to be a good fit at small numbers of residues, possibly because of an outlier at group lengthHELIX=4

# Try removing outlier and refitting

hlens2 <- hlens[hlens$lengthHELIX != 4,]

pmodel3 <- glm(len_cnt ~ lengthHELIX, family = "poisson", data = hlens2)
summary(pmodel3)

fitted2 <- fitted(pmodel3)
qplot(x = hlens2$lengthHELIX, fitted2)
ggplot() +
  geom_point(mapping = aes(x = hlens2$lengthHELIX, y = hlens2$len_cnt)) +
  geom_path(mapping = aes(x = hlens2$lengthHELIX, y = fitted2)) +
  ggtitle("Helix Length vs Poisson Fitted Counts")

ggplot() +
  geom_point(mapping = aes(x = hlens2$lengthHELIX, y = resid(pmodel3))) +
  geom_hline(yintercept = 0) +
  ggtitle("Helix Length vs Poisson Deviance Residuals")

# It still doesn't seem to be a good fit

