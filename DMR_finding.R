# conduct the LLR test###

# make all variables as numeric 
#the sample information 
samples <- geograbi.get.samples("GSE59065")
vars <- geograbi.extract.characteristics(samples)
rownames(vars) <- rownames(samples)
vars$`age (yr)` <- as.numeric(vars$`age (yr)`)
#rename the age column since this dataset the age is binary distribution
names(vars)[names(vars) == "age (yr)"] <- "age"

# make the age as categorical data
for (i in 1:nrow(vars)) {
        if(vars$age[i] > 50){
                vars$age[i]=1
        } else {
                vars$age[i]=0
        }
}
vars$age <- as.numeric(vars$age)
for (i in 1:nrow(vars)) {
        if(vars$gender[i] == "male"){
                vars$gender[i]=1
        } else {
                vars$gender[i]=0
        }
}
vars$gender <- as.numeric(vars$gender)

n=nrow(t.beta)

#number.positive.cpg.wald <- rep(NA,length(L.vec))
number.positive.cpg.LR <- rep(NA,length(L.vec))
count.L = 0
index1 <- NA
positi.record.LR1 <- NA
start.loc.B1 <- NA
end.loc.B1 <- NA
for (L in L.vec) {
        count.L = count.L +1
        
        positi.record.LR=rep(NA,length(which(end.loc.B-start.loc.B+1==L)))
        count=0
        for (i in which(end.loc.B-start.loc.B+1==L)){
                count=count+1
                index=start.loc.B[i]:end.loc.B[i]
                
                X = as.data.frame(cbind(vars$age,vars$gender))
                X= as.matrix(X,n,2) # two variables
                
                res=ht(t.beta[,index],X)
                positi.record.LR[count]=res$p.lr
                index2 <- index  
                index1 <- c(index1, index2)
                positi.record.LR1[index[1]] <- positi.record.LR[count]
                print(count)
        }
        number.positive.cpg.LR[count.L] <- L * length(positi.record.LR[positi.record.LR < 0.05/total_compare_LR])
}

# combine the position information and p-value
positi.record.LR2 <- positi.record.LR1[!is.na(positi.record.LR1)]
length(positi.record.LR2)
result_data_LR <- cbind(start.loc.B,end.loc.B,positi.record.LR2)
result_data_LR <- as.data.frame(result_data_LR)
colnames(result_data_LR) <- c("start","end","p_value")

save(result_data_LR,file = "result_data_LR.Rdata")

total_compare_LR
save(total_compare_LR,file = "total_compare_LR.Rdata")
# find the final DMRs
DMR_LR <- result_data_LR[result_data_LR$p_value < 0.05/total_compare_LR,]
DMR_LR 
DMR_LR$n <- DMR_LR$end - DMR_LR$start + 1
total_cpg <- sum(DMR_LR$n)
total_cpg
save(DMR_LR,file = "DMR_LR.Rdata")

# 
posit <- NA
for (i in 1:nrow(DMR_LR)) {
        posit1 <- DMR_LR$start[i]:DMR_LR$end[i]
        posit <- c(posit,posit1)
        
}
posit <- posit[-1]


#locate the cpg sites

sig_cpg <- beta[posit,]
cpg_name_LR <- c(rownames(sig_cpg))
cpg_name_LR


end_LR <- Sys.time()

time_LR <- end_LR - start_LR
time_LR
