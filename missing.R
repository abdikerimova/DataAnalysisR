library(leaps);
library(doMC);
library(dplyr);
library(purrr);
library(caret);
library(corrplot);
library(knitr);

data <- read.csv("AmesHousing.csv")
orig_data <- data;
all_x <- data

plot_missing_values <- function( data ) {
	nrows <- nrow(data)
	missing_data <- sort(map_dbl(data, function(x) sum(is.na(x))/nrows), decreasing=TRUE)
	# names_missing <- names(missing_data[missing_data > 0 ])
	missing_data <- missing_data[missing_data > 0.05];
	# print(length(missing_data));
	# print(names(missing_data));
	require(grDevices);
	png(filename="~/Desktop/missingValues.png");
	barplot(missing_data,yaxt='n',ylab="missing percentages",las=3);
	axis(1,labels(length(missing_data)),las=3);
	perc <- as.numeric(unname(missing_data));
	axis(2,at=pretty(perc),lab=paste(pretty(perc)*100,"%"),las=TRUE);
	dev.off();
}

plot_missing_values(all_x);

all_x$Pool.QC <- NULL
all_x$Misc.Feature <- NULL
all_x$Alley <- NULL
all_x$Fence <- NULL

m <- all_x$Lot.Frontage
x <- m[!is.na(m)]
med <- median(x)
m[is.na(m)] <- med
all_x$Lot.Frontage <- m

# Checking/removing near zero variance predictors from the dataset
nzv_cols<-nearZeroVar(data)
if ( length(nzv_cols) > 0 ) 
	data <- data[,-nzv_cols]

# Fireplace.Qu has >= 40% NAs, and less than 80%, that's why we introduce new level "SA"
m <- all_x$Fireplace.Qu
library(forcats)
all_x$Fireplace.Qu <- fct_explicit_na(m,"SA")

nrows <- nrow(all_x)
missing_data <- sort(map_dbl(all_x, function(x) sum(is.na(x))/nrows), decreasing=TRUE)
names_missing <- names(missing_data[missing_data > 0 ])
head(missing_data)
data <- all_x

# boilerplate code to replace
# missing values accordingly for each type of variable
replace_missing_values <- function( names_missing, data ) {
	m <- length(names_missing);
	res <- data;
	for ( i in c(1:m) ) {
		cat(names_missing[i]);
		s <- names_missing[i];
		if ( is.factor(data[,s]) ) {
			cat(" is categorical, replacing its N/As with ");
			mycol <- data[,s];
			x <- table(mycol[!is.na(mycol)]);
			moda <- names(which.max(x))
			cat(paste(moda,"\n",sep=""));
			mycol[is.na(mycol)] <- moda;
			res[,s] <- mycol;
		}
		else if ( is.integer(data[,s]) ) {
			cat(" is integer, replacing its N/As with ");
			mycol <- data[,s];
			x <- table(mycol[!is.na(mycol)]);
			moda <- as.integer(names(which.max(x)));
			cat(paste(moda,"\n",sep=""));
			mycol[is.na(mycol)] <- moda;
			res[,s] <- mycol;
		}
		else if ( is.logical(data[,s]) ) {
			cat(" is logical, replacing its N/As with ");
			mycol <- data[,s];
			x <- table(mycol[!is.na(mycol)]);
			moda <- names(which.max(x))
			cat(paste(moda,"\n",sep=""));
			mycol[is.na(mycol)] <- moda;
			res[,s] <- mycol;
		}
		else if ( is.numeric(data[,s]) ) {
			cat(" is numeric, replacing its N/As with ");
			mycol <- data[,s];
			x <- mycol[!is.na(mycol)];
			med <- median(x);
			cat(paste(med,"\n",sep=""));
			mycol[is.na(mycol)] <- med;
			res[,s] <- mycol;
		}
	}
	return(res);
};

mydata <- replace_missing_values(names_missing,data);
all_x <- mydata;
nrows <- nrow(all_x)
missing_data <- sort(map_dbl(all_x, function(x) sum(is.na(x))/nrows), decreasing=TRUE)
names_missing <- names(missing_data[missing_data > 0 ])
head(missing_data)

glue_together <- function( toExclude, prefix, vec, suffix, delim ) {
	if ( length(vec) == 0 ) {
		return(paste(suffix,prefix,sep=delim));
	}
	res <- "";
	len <- length(vec);
	for ( i in c(1:len) ) {
		# if ( vec[i] == "SalePrice" ) {
		if ( vec[i] == toExclude ) {
			# print("Excluding");
		}
		else if ( i > 1 ) {
			res <- paste(res,vec[i],sep=delim);
		}
		else {
			res <- paste(res,vec[i],sep="");
		}
	}
	res <- paste(prefix,res,sep="");
	res <- paste(res,suffix,sep="");
	return(res);
}

# stepwise selection of variables
iterative_anova <- function( mydata, threshold ) {
	model <- lm( SalePrice ~ ., data = mydata);
	prevLen <- 0;
	iterations <- 0;
	firstRoundExcluded <- "";
	piValues <- "";
	# anova(model)
	while ( TRUE ) {
		an <- anova(model);
		goodVariables <- row.names(an[an$"Pr(>F)" <= threshold,]);
		goodVariables <- head(goodVariables,-1);
		
		# print(goodVariables);
		if ( length(goodVariables) == prevLen ) {
			break;
		}
		myarg <- glue_together("SalePrice","SalePrice ~ ",goodVariables,"","+");
		var.out <- names(mydata)[!names(mydata) %in% c(goodVariables,"SalePrice")];
		if ( iterations == 0 ) {
			firstRoundExcluded <- head(var.out,-1);
			piValues <- head(an[an$"Pr(>F)" > threshold,]$"Pr(>F)",-1);
		}
		# print(var.out);
		#for ( i in c(1:length(var.out)) ) {
			# noquote(paste("mydata$",var.out[i],"<-NULL",sep=""));
		#}
		mydata <- mydata[,!names(mydata) %in% var.out, drop = F]
		# print(names(mydata));
		model <- lm(as.formula(myarg),data = mydata);
		prevLen <- length(goodVariables);
		iterations <- iterations+1;
	}
	return(list("model"=model,"data"=mydata,"excluded"=firstRoundExcluded,"pivalues"=piValues));
}

get_collinear <- function( mydata, thre = 0.75, notToTouch ) {
  # Expects data dataframe
  x <- mydata;
  x$"SalePrice" <- NULL;
  num_cols <- ncol(x)
  collinear_vec <- vector("character", num_cols) 
  index <- 1;
  vec <- c();
  for (i in seq(1:num_cols)) {
    corMat <- cor(x);
	corMat[is.na(corMat)] <- 0;
    diag(corMat) <- 0;
    df_cols <- names(x);
    AB <- which(abs(corMat) == max(abs(corMat)), arr.ind = TRUE);
    if ( abs(corMat[AB][[1]]) > thre ) {
      names_AB <- rownames(AB)
      if (sum(abs(corMat[names_AB[1], ])) > sum(abs(corMat[names_AB[2], ]))) {
        collinear_vec[index] = names_AB[1]
        index <- index + 1;
		vec <- c(vec,names_AB[1],names_AB[2]);
      } else {
		collinear_vec[index] = names_AB[2]
        index <- index + 1;
		vec <- c(vec,names_AB[2],names_AB[1]);
	  }
      x <- select(x, one_of(setdiff(df_cols, collinear_vec[index - 1])))
    } else {
		break
	}
  }
  return(list("collinear_vec"=collinear_vec[collinear_vec != ""],"pairs"=vec));
}

num_data <- select_if(mydata, is.numeric);
collinear_data <- get_collinear(num_data,0.75,"SalePrice");
pl <- collinear_data$collinear_vec;
vec <- collinear_data$pairs;
newdata <- mydata[!names(mydata) %in% pl];
i = 1;
retained <- c();
excluded <- c();
while ( i <= length(vec) ) {
	cat(paste("retained: ",vec[i],", exluded: ",vec[i+1],"\n"));
	retained <- c(retained,vec[i]);
	excluded <- c(excluded,vec[i+1]);
	i <- i+2;
}
df <- data.frame("retained"=retained,"excluded"=excluded,"correlation"=mapply(function(x,y) cor(x,y), all_x[,retained],all_x[,excluded]));
write.csv(df,file="~/Desktop/retainedExcluded.csv");
mydata <- newdata;
modelAndData <- iterative_anova(mydata,0.05);
finaldata <- modelAndData$data;
# write.csv(modelAndData$model,"~/Desktop/model.csv");
write.csv(finaldata,file = "~/Desktop/houseDataCleanedOfMissingData.csv");
vec <- modelAndData$pivalues;
names(vec) <- modelAndData$excluded;
write.csv(vec,file="~/Desktop/pivaluesOfExcluded.csv");
# plot missing data versus column name
# an <- modelAndData$model;
# github.com
#threshold = 0.05;
#model <- lm( SalePrice ~ ., data = mydata);
#goodVariables <- row.names(an[an$"Pr(>F" <= threshold,]);
# n <- dim(modelAndData$data)[1];
#model_1=regsubsets(SalePrice  ~ ., data = modelAndData$data, direction = "forward", nvmax = 4);
#rs=summary(model_1)
#msize=1:4
#Aic = n*log(rs$rss/n) + 2*(msize+1)
#res <- which.min(Aic);
#ress <- msize[which.min(Aic)]
