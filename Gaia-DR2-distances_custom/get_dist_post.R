args <- commandArgs(TRUE)
fname <- args[1]
w <- as.double(args[2])
wsd <- as.double(args[3])
glon <- as.double(args[4])
glat <- as.double(args[5])
rlen <- NA

# approx distance
r <- 1e3 / w
e_r <- r * (wsd/w)
rlo <- r - 10*e_r
rhi <- r + 10*e_r
rplotlo <- rlo
rplothi <- rhi
source("distest_single.R")
