

library(CTS)

progressr::handlers(global=T)
progressr::handlers("progress")


#simulated data
s1 <- gen_Sample()
res <- calc_be(s1)
ext_res <- extractBE(res)
make_report(sim=s1,be=ext_res)


#Buproprion

s2 <- get_data("~/LAPK/Pmetrics/Buproprion/src/all_data.csv")
res2 <- calc_be(s2)
ext_res2 <- extractBE(res2)
make_report(sim=s2,be=ext_res2)
