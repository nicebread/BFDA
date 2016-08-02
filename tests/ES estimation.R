load("../finalSims/BFDA.0.5.prior.RData")

BFDA.analysis(BFDA.0.5.prior, boundary=6)

BFDA.analysis(BFDA.0.5.prior, boundary=10, n.max=200)



n.min=20
n.max=200
boundary=c(1/6, 6)

sim <- BFDA.0.5.prior$sim

# reduce simulation to relevant data
sim <- sim %>% filter(n >= n.min, n <= n.max)

logBoundary <- log(boundary)

# For the densities: Data frames of stopping times / stopping BFs
n.max.hit <- sim %>% group_by(id) %>% filter(n == n.max, max(logBF) <= logBoundary[2] & min(logBF) >= logBoundary[1])

# reduce to *first* break of a boundary
boundary.hit <- sim %>% group_by(id) %>%
	filter(logBF>=logBoundary[2] | logBF<=logBoundary[1]) %>%
	filter(row_number()==1) %>% ungroup()	

# endpoint stores all hits (boundary or n.max)
endpoint <- bind_rows(n.max.hit, boundary.hit)

# compute variance of Cohen's d
endpoint <- endpoint %>% mutate(vi=((n+n) / n^2) + d.emp^2/(2*n+n))

mean(sim$d) # 0.5
weighted.mean(endpoint$d.emp, w=1/endpoint$vi)

plot(endpoint$n, endpoint$d.emp)