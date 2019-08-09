#R source code


library(glmnet)
library(pROC)
#read the miRNA, RNA and clinical information respectively
miRNA=read.csv(file.choose(),header=T)
RNA=read.csv(file.choose(),header=T)
y=read.csv(file.choose(),header=T)

RNA=RNA[,6:dim(RNA)[2]]
rownames(miRNA)=miRNA[,1]
rownames(RNA)=RNA[,1]
miRNA=miRNA[,2:dim(miRNA)[2]]
RNA=RNA[,2:dim(RNA)[2]]

y_sub=y
miRNA_sub=t(miRNA)
RNA_sub=t(RNA)

#log transform
RNA_sub=log2(RNA_sub+1e-16)
miRNA_sub=log2(miRNA_sub+1e-16)

#############################
######1. Pre-selection#######
#############################
#mRNA pre-selection using paired t-test
p_value=rep(NA,dim(RNA_sub)[2])
for(i in 1:dim(RNA_sub)[2])
{
  if(class(try(t.test(RNA_sub[,i]~y_sub$labels,paired = TRUE)$p.value))!="try-error")
  {p_value[i]=t.test(RNA_sub[,i]~y_sub$labels,paired = TRUE)$p.value}
  
}
fdr=p.adjust(p_value,method="fdr")
fc=rep(NA,dim(RNA_sub)[2])
for(i in 1:dim(RNA_sub)[2])
{
  fc[i]=median(2^(RNA_sub[y_sub$labels=="nHCC",i]))/median(2^(RNA_sub[y_sub$labels=="Adj.Norm",i]))
}
fc[which(fc<1)]=-1/fc[which(fc<1)]
table=cbind(colnames(RNA_sub),p_value,fdr,fc)
#To check p-values, FDR and Fold Change of all variables can use the command
#write.csv(table,"RNA.csv")

#We only selected the RNAs which FDR<0.05
RNA_sub_selected=RNA_sub[,which(fdr<0.05)]


#miRNA pre-selection using t-test
p_value=rep(NA,dim(miRNA_sub)[2])
for(i in 1:dim(miRNA_sub)[2])
{
  if(class(try(t.test(miRNA_sub[,i]~y_sub$labels,paired = TRUE)$p.value))!="try-error")
  {p_value[i]=t.test(miRNA_sub[,i]~y_sub$labels,paired = TRUE)$p.value}
  
}
fdr=p.adjust(p_value,method="fdr")
fc=rep(NA,dim(miRNA_sub)[2])
for(i in 1:dim(miRNA_sub)[2])
{
  fc[i]=median(2^(miRNA_sub[y_sub$labels=="nHCC",i]))/median(2^(miRNA_sub[y_sub$labels=="Adj.Norm",i]))
}
fc[which(fc<1)]=-1/fc[which(fc<1)]
table=cbind(colnames(miRNA_sub),p_value,fdr,fc)
#To check p-values, FDR and Fold Change of all variables can use the command
#write.csv(table,"miRNA.csv")

#We only selected the miRNAs which FDR<0.05
miRNA_sub_selected=miRNA_sub[,which(fdr<0.05)]


####################################
######2. Mechanistic submodel#######
####################################
pairs=matrix(ncol=3)
for(g in colnames(RNA_sub_selected))
{
  x_tmp=miRNA_sub_selected
  y_tmp=as.matrix(RNA_sub_selected[,colnames(RNA_sub_selected)==g])
  if(sum(y_tmp==log2(0+1e-16))!=58)
  {
    colnames(y_tmp)=g
    fit=glmnet(x=x_tmp,y=y_tmp,alpha=1)
    cv.fit=try(cv.glmnet(x=x_tmp,y=y_tmp,alpha=1,nfolds=10))
    if(class(cv.fit)!="try-error")
    {
      coef=coef(fit,s=cv.fit$lambda.1se)
      nom=rownames(coef)[-1]
      coef=coef[-1]
      if(sum(coef!=0)>=1)
      {
        ana=lm(y_tmp~x_tmp[,which(coef!=0)])
        table=rep(g,sum(coef!=0))
        table=cbind(table,nom[which(coef!=0)])
        table=cbind(table,ana$coefficients[-1])
        pairs=rbind(pairs,table)
      }
    }
  }
}
pairs=pairs[-1,]
colnames(pairs)=c("mRNA","miRNA","coefficient")
pairs=as.data.frame(pairs)
rownames(pairs)=NULL
#To check miRNA-mRNA associations can use the command
#write.csv(pairs,"pairs.csv")


#################################
######3. clinical submodel#######
#################################

######### Model lasso ##########
#unselected miRNAs and mRNAs
miRNA_unselected=setdiff(colnames(miRNA_sub_selected),pairs$miRNA)
RNA_unselected=setdiff(colnames(RNA_sub_selected),pairs$mRNA)

A=matrix(nrow=58)
B=matrix(nrow=58)
for(g in unique(pairs$mRNA))
{
  y_tmp=RNA_sub_selected[,colnames(RNA_sub_selected)==g]
  x=pairs$miRNA[pairs$mRNA==g]
  if(length(x)==1){
    x_tmp=miRNA_sub[,colnames(miRNA_sub)==x]
  }else{
    x_tmp=miRNA_sub[,colnames(miRNA_sub)%in%x]}
  ana=lm(y_tmp~x_tmp)
  A=cbind(A,predict(ana))
  colnames(A)[dim(A)[2]]=paste(g,"miRNA",sep = "-")
  B=cbind(B,residuals(ana))
  colnames(B)[dim(B)[2]]=paste(g,"other",sep = "-")
}
A=A[,-1]
B=B[,-1]
if(length(miRNA_unselected)==1)
{
  C=miRNA_sub_selected[,colnames(miRNA_sub_selected)==miRNA_unselected]
}else{
  C=miRNA_sub_selected[,colnames(miRNA_sub_selected)%in%miRNA_unselected]
}
D=RNA_sub_selected[,colnames(RNA_sub_selected)%in%RNA_unselected]
X=cbind(A,B,C,D)

fit=glmnet(x=X,y=y_sub$labels,family = "binomial")
cv.fit=cv.glmnet(x=X,y=y_sub$labels,family = "binomial")
coef=coef(fit,s=cv.fit$lambda.1se)
coef=coef[-1]
coef_notzero=which(coef!=0)
g=unique(pairs$mRNA)[coef_notzero[coef_notzero<dim(A)[2]]]
pairs_Modellasso=pairs[pairs$mRNA%in%g,1:2]
rownames(pairs_Modellasso)=NULL
if(max(coef_notzero)>dim(A)[2])
{
  table=matrix(c(rep("Y",sum(coef_notzero>dim(A)[2])),colnames(X)[coef_notzero[coef_notzero>dim(A)[2]]]),ncol=2)
  colnames(table)=colnames(pairs_Modellasso)
  pairs_Modellasso=rbind(pairs_Modellasso,table)
}
print("miRNA-mRNA pairs selected by Model lasso is")
print(pairs_Modellasso)
