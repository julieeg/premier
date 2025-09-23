# PID of current job: 1988117
mSet<-InitDataObjects("conc", "msetora", FALSE, 150)
cmpd.vec<-c("HMDB0011516","HMDB0010403","HMDB0000631","HMDB0005060","HMDB0000138","HMDB0000439","HMDB0000271","HMDB0011517")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "hmdb");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 150, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_1_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_1_", "png", 150, width=NA)
mSet<-SaveTransformedData(mSet)
