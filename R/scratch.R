## Extract column names for metadata:
write.csv(a[8:length(a)], file = "16S_metadata.csv", row.names = FALSE)

#The following statement goest through a data frame and returns all 
#rows such that the value !=0
row_sub = apply(otuTable, 1, function(row) all(row != 0))
all <- otuTable[row_sub,]