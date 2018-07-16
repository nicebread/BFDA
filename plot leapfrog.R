library(tidyverse)
load(file="fullRes.RData")
res <- bind_rows(fullRes)
res$unique <- factor(paste0(res$batch, ":", res$iteration, ":", res$upper.boundary, ":", res$lower.boundary, ":", res$n.min, ":", res$ n.max))

res <- res %>% group_by(unique) %>% 
	mutate(treatmentsTested = n()
)


final <- res %>% group_by(unique) %>%
	filter(row_number() == n())
	
summ <- final %>% 	
	group_by(upper.boundary, lower.boundary, n.min, n.max) %>% 
	summarise(
		MEAN = mean(standard.ES),
		LQ = quantile(standard.ES, prob=.25),
		UQ = quantile(standard.ES, prob=.75),
		treatmentsTested = mean(treatmentsTested)
	)

ggplot(summ, aes(x=factor(upper.boundary), shape=factor(lower.boundary))) + geom_pointrange(aes(y=MEAN, ymin=LQ, ymax=UQ), position=position_dodge(width=.4)) + facet_grid(n.min~n.max) + theme_bw() + ylab("Mean final effect size of new treatment (+ quartiles)")


ggplot(summ, aes(x=factor(upper.boundary), y=treatmentsTested, shape=factor(lower.boundary))) + geom_point(position=position_dodge(width=.4)) + facet_grid(n.min~n.max) + theme_bw()

# ---------------------------------------------------------------------
# Plot a single run

run1 <- res %>% 
	filter(lower.boundary == 1/3, upper.boundary == 10, n.min == 10, n.max == 200, batch==1, iteration==1) %>% 
	#filter(unique == "3:5:3:0.333333333333333:10:200") %>% 
	mutate(
			n.cum = cumsum(n.used),
			arm = 1:n(),
			decision = lead(decisionQuality)
		)

ggplot(run1, aes(x=arm, y=standard.ES)) + 	
	geom_line(aes(group=1)) + 
	geom_point(aes(shape=decision, color=decision), size=2) + 
	theme_bw()