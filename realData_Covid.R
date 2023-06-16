


response <- "median_R_mean"# "casegrowth"

data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency)$data
