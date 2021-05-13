

library(CTS)

progressr::handlers(global=T)
progressr::handlers("progress")


s1 <- gen_Sample()

plotQ(s1)

res <- calc_be(s1)

ext_res <- extractBE(res)

plotBE(ext_res$auc_par, type="ratio")
plotBE(ext_res$cmax_par, type="ratio")
plotBE(ext_res$auc_cross, type="ratio")
plotBE(ext_res$cmax_cross, type="ratio")
plotBE(ext_res$auc_rep, type="ratio")
plotBE(ext_res$cmax_rep, type="ratio")

plotBE(ext_res$auc_par, type="pBE")
plotBE(ext_res$cmax_par, type="pBE")
plotBE(ext_res$auc_cross, type="pBE")
plotBE(ext_res$cmax_cross, type="pBE")
plotBE(ext_res$auc_rep, type="pBE")
plotBE(ext_res$cmax_rep, type="pBE")


