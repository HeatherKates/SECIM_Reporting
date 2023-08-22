#######################
####Feature-wise#######
#######################
#Features are in cols##
#######################
Norm_Plots <- function(proc.data,norm.data){

GetRandomSubsetIndex<-function(total, sub.num = 50){
  if(total < sub.num){
    1:total
  }else{
    sample(1:total, sub.num)
  }
}

pre.inx<-GetRandomSubsetIndex(ncol(proc.data), sub.num=50);
namesVec <- colnames(proc.data[,pre.inx, drop=FALSE]);

# only get common ones
nm.inx <- namesVec %in% colnames(norm.data)
namesVec <- namesVec[nm.inx]
pre.inx <- pre.inx[nm.inx]

norm.inx<-match(namesVec, colnames(norm.data))
namesVec <- substr(namesVec, 1, 12)

rangex.pre <- range(proc.data[, pre.inx, drop=FALSE], na.rm=T)
rangex.norm <- range(norm.data[, norm.inx, drop=FALSE], na.rm=T)


# 1/4 Density plot of features before normalization
prenorm.density.data.feature <- density(apply(proc.data, 2, mean, na.rm=TRUE))
p1=
  ggplot(data.frame(cbind(prenorm.density.data.feature[[1]],prenorm.density.data.feature[[2]])),
         aes(X1,X2))+geom_line()+labs(x="",y="Density")+
  ggtitle("Before Normalization")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(size = 8, face = "bold"))

# 2/4 Boxplot of features before normalization
p2=
ggplot(stack(data.frame(proc.data[,norm.inx , drop=FALSE])),aes(values,ind))+geom_boxplot()+
  labs(x="Raw Peak Intensity",y="Up to 50 peaks")+
  theme(axis.text.y=element_blank())#+ylim(rangex.norm)

# 3/4 Density plot of features after normalization
postnorm.density.data.feature <- density(apply(norm.data, 2, mean, na.rm=TRUE))
p3=
ggplot(data.frame(cbind(postnorm.density.data.feature[[1]],postnorm.density.data.feature[[2]])),
       aes(X1,X2))+geom_line()+labs(x="",y="")+
  ggtitle("After Normalization")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(size = 8, face = "bold"))

# 4/4 Boxplot of features post normalization
p4=
ggplot(stack(data.frame(norm.data[,norm.inx , drop=FALSE])),aes(values,ind))+geom_boxplot()+
  labs(x="Normalized Peak Intensity",y="")+
  theme(axis.text.y=element_blank())#+ylim(rangex.norm)


#######################
####Sample-wise########
#######################
#samples are in rows###
#######################
pre.inx<-GetRandomSubsetIndex(nrow(proc.data), sub.num=50);
namesVec <- rownames(proc.data[pre.inx, , drop=FALSE]);

# only get common ones
nm.inx <- namesVec %in% rownames(norm.data)
namesVec <- namesVec[nm.inx]
pre.inx <- pre.inx[nm.inx]

norm.inx<-match(namesVec, rownames(norm.data))
namesVec <- gsub("QE2_jm_040_","",namesVec)
rangex.pre <- range(proc.data[pre.inx, , drop=FALSE], na.rm=T)
rangex.norm <- range(norm.data[norm.inx, , drop=FALSE], na.rm=T)

# 1/4 Density plot of samples before normalization
prenorm.density.data.sample <- density(apply(proc.data, 1, mean, na.rm=TRUE))
p5=
ggplot(data.frame(cbind(prenorm.density.data.sample[[1]],prenorm.density.data.sample[[2]])),
                           aes(X1,X2))+geom_line()+labs(x="",y="Density")+
  ggtitle("Before Normalization")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(size = 8, face = "bold"))

# 2/4 Boxplot of samples before normalization
p6=
ggplot(stack(data.frame(t(proc.data[norm.inx, , drop=FALSE]))),aes(values,ind))+geom_boxplot()+
  labs(x="Raw Peak Intensity",y="Samples")+
  theme(axis.text.y=element_blank())#+ylim(rangex.norm)
# 3/4 Density plot of samples after normalization
postnorm.density.data.sample <- density(apply(norm.data, 1, mean, na.rm=TRUE))
p7=
  ggplot(data.frame(cbind(postnorm.density.data.sample[[1]],postnorm.density.data.sample[[2]])),
          aes(X1,X2))+geom_line()+labs(y="",x="")+
  ggtitle("After Normalization")+
  theme(plot.title = element_text(size = 8, face = "bold"))

# 4/4 Boxplot of samples after normalization THIS ONE I FIXED, DO THE REST
p8=
  ggplot(stack(data.frame(t(norm.data[norm.inx, , drop=FALSE]))),aes(values,ind))+geom_boxplot()+
  labs(x="Normalized Peak Intensity",y="")+
  theme(axis.text.y=element_blank())#+ylim(rangex.norm)
return(p1)
return(p2)
return(p3)
return(p4)
return(p5)
return(p6)
return(p7)
return(p8)
}

