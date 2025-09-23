# PID of current job: 1990372
mSet<-InitDataObjects("conc", "msetora", FALSE, 150)
cmpd.vec<-c("HMDB0010393","HMDB0029723","HMDB0011516","HMDB0011506","HMDB0010404","HMDB0042063","HMDB0011507","HMDB0011477","HMDB0042301","HMDB0000247","HMDB0010386","HMDB0011487","HMDB0042068","HMDB0010403","HMDB0000482","HMDB0042099","HMDB0042093","HMDB0007170","HMDB0005356","HMDB0042278","HMDB0042751","HMDB0011526","HMDB0099325","HMDB0010401","HMDB0008993","HMDB0007012","HMDB0005365","HMDB0072780","HMDB0042098")
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
