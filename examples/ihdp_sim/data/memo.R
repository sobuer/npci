# source("examples/ihdp_sim/data.R")

# dataFile <- file.path("examples/ihdp_sim/data", "ihdp.RData")
# load(dataFile)
# ihdp <- subset(ihdp, treat != 1 | momwhite != 0)

# covariateNames <- c("bw", "b.head", "preterm", "birth.o", "nnhealth", "momage",
#                     "sex", "twin", "b.marr", "mom.lths", "mom.hs", "mom.scoll",
#                     "cig", "first", "booze", "drugs", "work.dur", "prenatal",
#                     "ark", "ein", "har", "mia", "pen", "tex", "was")

# x <- ihdp[,covariateNames]          #データを共変量だけに
# trans <- npci:::getTransformations(x)
# x <- npci:::transform(x, trans$standardize)     #連続値を標準化
# z <- ihdp$treat     #処置の有無
# w <- 0.5

# # callingEnv$x <- x
# # callingEnv$z <- z

# # ## places mu.0, mu.1, y.0, y.1, and y into calling env, as well as x.r and z.r
# # generateDataForIterInCurrentEnvironment(iter, x, z, w, overlap = TRUE, covariates = "select", setting = "A")
# iter <- 1L
# overlap = TRUE

# terms <- character()        #？？？よく分からん
# isBinary <- sapply(seq_len(ncol(x)), function(j) length(unique(x[,j])) == 2)        #離散値確認？
# for (i in seq_len(ncol(x))) {
#     if (!isBinary[i]) terms[length(terms) + 1L] <- paste0("I(", colnames(x)[i], "^2)")
    
#     if (i < ncol(x)) for (j in seq(i + 1L, ncol(x)))
#         terms[length(terms) + 1L] <- paste0(colnames(x)[i], ":", colnames(x)[j])
# }

# # callingEnv <- parent.frame(1L)
# set.seed(iter * 5L + if (iter <= 500L) 565L else 7565L)

# #共変量のデータフレーム
# x.m <- if (is.data.frame(x)) dbarts::makeModelMatrixFromDataFrame(x) else x
# # callingEnv$x.r <- x.m

# n <- nrow(x)
# p <- ncol(x) ## main effects

# #x.m <- cbind(1.0, x.m)
# sigma.y <- 1.0
# tau <- 4.0
# # w.full <- matrix(c(0.0, rep_len(w, ncol(x.m) - 1)), n, ncol(x.m), byrow = TRUE)

# # vals <- seq(0.0, 0.4, 0.1)      #0.0~0.4まで0.1ずつ（B:ランダムにサンプリングされた回帰係数）
# # probs <- max(1.0 - 2 / sqrt(p), 0.2)        #0.6
# # probs <- c(probs, rep(0.25 * (1.0 - probs), 4L))        #β_B
# # beta <- c(sample(seq(-1, 1, 0.25), 1), sample(vals, p, replace = TRUE, prob = probs))

# vals <- seq(0, 4, 1)            #βの係数：0,1,2,3,4
# probs <- c(0.5, 0.2, 0.15, 0.1, 0.05)           #β_A
# beta <- sample(vals, p, replace = TRUE, prob = probs)

# # mu.0 <- exp((x.m + w.full) %*% beta)
# # mu.1 <- x.m %*% beta
# mu.0 <- x.m %*% beta
# mu.1 <- x.m %*% beta + 4

# # callingEnv$z.r <- z

# # omega <- if (overlap == TRUE){
# #     mean(mu.1[z == 1] - mu.0[z == 1]) - tau
# # }else{
# #     mean(mu.1[z == 0] - mu.0[z == 0]) - tau
# # }

# # mu.1 <- mu.1 - omega

# y.0 <- rnorm(n, mu.0, sigma.y)
# y.1 <- rnorm(n, mu.1, sigma.y)
# y <- ifelse(z == 1, y.1, y.0)

# head(mu.0)
# head(mu.1)
# head(y.0)   #非処置時の結果
# head(y.1)   #処置時の結果
head(y)



####################################################################################################

source("examples/ihdp_sim/data.R")
dataFile <- file.path("examples/ihdp_sim/data", "ihdp.RData")
load(dataFile)

ihdp <- subset(ihdp, treat != 1 | momwhite != 0)
covariateNames <- c("bw", "b.head", "preterm", "birth.o", "nnhealth", "momage",
                    "sex", "twin", "b.marr", "mom.lths", "mom.hs", "mom.scoll",
                    "cig", "first", "booze", "drugs", "work.dur", "prenatal",
                    "ark", "ein", "har", "mia", "pen", "tex", "was")

x <- ihdp[,covariateNames]          #データを共変量だけに
trans <- npci:::getTransformations(x)
x <- npci:::transform(x, trans$standardize)     #連続値を標準化
z <- ihdp$treat     #処置の有無
w <- 0.5
overlap = TRUE

terms <- character()        #？？？よく分からん
isBinary <- sapply(seq_len(ncol(x)), function(j) length(unique(x[,j])) == 2)        #離散値確認？
for (i in seq_len(ncol(x))) {
    if (!isBinary[i]) terms[length(terms) + 1L] <- paste0("I(", colnames(x)[i], "^2)")
    
    if (i < ncol(x)) for (j in seq(i + 1L, ncol(x)))
        terms[length(terms) + 1L] <- paste0(colnames(x)[i], ":", colnames(x)[j])
}


#1000回分のシミュレーション（A）
for(iter in seq(2)){
    set.seed(iter * 5L + if (iter <= 500L) 565L else 7565L)
    #set.seed(iter)

    #共変量のデータフレーム
    x.m <- if (is.data.frame(x)) dbarts::makeModelMatrixFromDataFrame(x) else x
    n <- nrow(x)
    p <- ncol(x) ## main effects
    
    ##x.m <- cbind(1.0, x.m)
    sigma.y <- 1.0
    tau <- 4.0
    vals <- seq(0, 4, 1)            #βの係数：0,1,2,3,4
    probs <- c(0.5, 0.2, 0.15, 0.1, 0.05)           #β_A


    beta <- sample(vals, p, replace = TRUE, prob = probs)
    print(beta)
    mu.0 <- x.m %*% beta
    mu.1 <- x.m %*% beta + 4

    y.0 <- rnorm(n, mu.0, sigma.y)
    y.1 <- rnorm(n, mu.1, sigma.y)
    y <- ifelse(z == 1, y.1, y.0)

    results <- cbind(x.m, y, y.1, y.0)
    print(head(results))
    ##csvで保存
    # fileName <- paste0("ihdp_", iter, ".csv")
    # save(results, file = file.path("examples/ihdp_sim/data", fileName))
}



####################################################################################################
##シミュレーション（B）
## これでやろう！
source("examples/ihdp_sim/data.R")
dataFile <- file.path("examples/ihdp_sim/data", "ihdp.RData")
load(dataFile)

ihdp <- subset(ihdp, treat != 1 | momwhite != 0)
covariateNames <- c("bw", "b.head", "preterm", "birth.o", "nnhealth", "momage",
                    "sex", "twin", "b.marr", "mom.lths", "mom.hs", "mom.scoll",
                    "cig", "first", "booze", "drugs", "work.dur", "prenatal",
                    "ark", "ein", "har", "mia", "pen", "tex", "was")

x <- ihdp[,covariateNames]          #データを共変量だけに
trans <- npci:::getTransformations(x)
x <- npci:::transform(x, trans$standardize)     #連続値を標準化
z <- ihdp$treat     #処置の有無
w <- 0.5
overlap = TRUE

terms <- character()        #？？？よく分からん
isBinary <- sapply(seq_len(ncol(x)), function(j) length(unique(x[,j])) == 2)        #離散値確認？
for (i in seq_len(ncol(x))) {
    if (!isBinary[i]) terms[length(terms) + 1L] <- paste0("I(", colnames(x)[i], "^2)")
    
    if (i < ncol(x)) for (j in seq(i + 1L, ncol(x)))
        terms[length(terms) + 1L] <- paste0(colnames(x)[i], ":", colnames(x)[j])
}

#1000回分のシミュレーション（B）
for(iter in seq(1000)){
    set.seed(iter * 5L + if (iter <= 500L) 565L else 7565L)
    #set.seed(iter)

    #共変量のデータフレーム
    x.m <- if (is.data.frame(x)) dbarts::makeModelMatrixFromDataFrame(x) else x
    n <- nrow(x)
    p <- ncol(x) ## main effects

    x.m <- cbind(1.0, x.m)
    sigma.y <- 1.0
    tau <- 4.0
    w.full <- matrix(c(0.0, rep_len(w, ncol(x.m) - 1)), n, ncol(x.m), byrow = TRUE)

    vals <- seq(0.0, 0.4, 0.1)      #0.0~0.4まで0.1ずつ（B:ランダムにサンプリングされた回帰係数）
    probs <- max(1.0 - 2 / sqrt(p), 0.2)        #0.6
    probs <- c(probs, rep(0.25 * (1.0 - probs), 4L))        #β_B
    beta <- c(sample(seq(-1, 1, 0.25), 1), sample(vals, p, replace = TRUE, prob = probs))

    #print(beta)

    mu.0 <- exp((x.m + w.full) %*% beta)
    mu.1 <- x.m %*% beta

    omega <- if (overlap == TRUE){
        mean(mu.1[z == 1] - mu.0[z == 1]) - tau
    }else{
        mean(mu.1[z == 0] - mu.0[z == 0]) - tau
    }

    mu.1 <- mu.1 - omega


    y.0 <- rnorm(n, mu.0, sigma.y)
    y.1 <- rnorm(n, mu.1, sigma.y)
    y <- ifelse(z == 1, y.1, y.0)
    yc <- ifelse(z == 0, y.1, y.0)
    colnames(mu.1) <- "mu.1"
    colnames(mu.0) <- "mu.0"


    results <- cbind(x.m[,2:ncol(x.m)], y, yc, y.1, y.0, mu.1, mu.0)
    #print(head(results))
    # print(mean(mu.1))
    # print(mean(mu.0))
    # print(mean(mu.1)-mean(mu.0))
    # print(mean(y.1)-mean(y.0))
    # print(head(y))
    # print(head(yc))
    # print("#####")
    #print(mu.1[z == 1] - mu.0[z == 1])
    # print(mean(y.1[z == 1] -  y.0[z == 1]))         #答え
    # print(mean(y.1 -  y.0))
    # print(mean(mu.1[z == 1] -  mu.0[z == 1]))       #4

    ##csvで保存
    fileName <- paste0("ihdp_", iter, ".csv")
    folderPath <- "IHDP"
    filePath <- file.path(folderPath, fileName)
    write.csv(results, file = filePath, row.names = FALSE)
}

# data_read <- read.csv("IHDP/ihdp_2.csv")
# head(data_read)



####################################################################################################
## npciパッケージ（B）
## CEVAEのデータ
source("examples/ihdp_sim/data.R")
dataFile <- file.path("examples/ihdp_sim/data", "ihdp.RData")
load(dataFile)
ihdp <- subset(ihdp, treat != 1 | momwhite != 0)
covariateNames <- c("bw", "b.head", "preterm", "birth.o", "nnhealth", "momage",
                    "sex", "twin", "b.marr", "mom.lths", "mom.hs", "mom.scoll",
                    "cig", "first", "booze", "drugs", "work.dur", "prenatal",
                    "ark", "ein", "har", "mia", "pen", "tex", "was")

x <- ihdp[,covariateNames]          #データを共変量だけに
trans <- npci:::getTransformations(x)
x <- npci:::transform(x, trans$standardize)     #連続値を標準化
z <- ihdp$treat     #処置の有無
w <- 0.5
overlap = TRUE
covariates <- "full"

terms <- character()        #？？？よく分からん
isBinary <- sapply(seq_len(ncol(x)), function(j) length(unique(x[,j])) == 2)        #離散値確認？
for (i in seq_len(ncol(x))) {
    if (!isBinary[i]) terms[length(terms) + 1L] <- paste0("I(", colnames(x)[i], "^2)")
    
    if (i < ncol(x)) for (j in seq(i + 1L, ncol(x)))
        terms[length(terms) + 1L] <- paste0(colnames(x)[i], ":", colnames(x)[j])
}

getQuadraticTerms <- function(x){
    terms <- character()
    isBinary <- sapply(seq_len(ncol(x)), function(j) length(unique(x[,j])) == 2)
    for (i in seq_len(ncol(x))) {
      if (!isBinary[i]) terms[length(terms) + 1L] <- paste0("I(", colnames(x)[i], "^2)")
      
      if (i < ncol(x)) for (j in seq(i + 1L, ncol(x)))
        terms[length(terms) + 1L] <- paste0(colnames(x)[i], ":", colnames(x)[j])
    }
    terms
}



for (iter in c(501,1000)) { 
  message("  running iter ", iter, "\n", sep = "")
  set.seed(iter * 5L + if (iter <= 500L) 565L else 7565L)

  if (is.numeric(covariates)) {
    x <- x[,sample(ncol(x), covariates)]
  }

  x.m <- if (is.data.frame(x)) dbarts::makeModelMatrixFromDataFrame(x) else x
  if (is.integer(x.m))
    x.m <- matrix(as.double(x.m), nrow(x.m), dimnames = dimnames(x.m))

  n <- nrow(x)
  p <- ncol(x)

  mainEffects <- colnames(x)
  quadEffects <- getQuadraticTerms(x)
  probs.q.0 <- max(1.0 - 5^(1 - 1 / 8) / (p^(15/16)), 0.4)
  quadEffects <- quadEffects[rbinom(length(quadEffects), 1, 1 - probs.q.0) == 1]
  formulaString <- paste0("y ~ ", paste0(mainEffects, collapse = " + "),
                          " + ", paste0(quadEffects, collapse = " + "))
  temp <- x
  temp$y <- numeric(nrow(x))
  mod <- glm(formulaString, data = temp, x = TRUE)
  coefs <- mod$coef[-1L]
  x.m <- mod$x[,-1L]
  x.m <- x.m[,!is.na(coefs)]

  x.m <- cbind(1.0, x.m)
  sigma.y <- 1.0
  tau <- 4.0
  w.full <- matrix(c(0.0, rep_len(w, ncol(x.m) - 1)), n, ncol(x.m), byrow = TRUE)


  vals.m  <- c(0.0, 1.0, 2.0)
  probs.m <- max(1.0 - 2.0 / sqrt(p), 0.2) ## hits 0.6 w/25 covariates, 0.8 w/108
  probs.m <- c(probs.m, 0.75 * (1.0 - probs.m), 0.25 * (1.0 - probs.m))

  vals.q  <- c(0.5, 1.0)
  probs.q <- c(0.75, 0.25)

  beta.m0 <- sample(vals.m, p + 1, replace = TRUE, prob = probs.m)
  beta.m1 <- sample(vals.m, p + 1, replace = TRUE, prob = probs.m)
  beta.q0 <- sample(vals.q, ncol(x.m) - p - 1, replace = TRUE, prob = probs.q)
  beta.q1 <- sample(vals.q, ncol(x.m) - p - 1, replace = TRUE, prob = probs.q)

  mu.0 <- x.m %*% c(beta.m0, beta.q0)
  mu.1 <- x.m %*% c(beta.m1, beta.q1)

  omega <- if (overlap == TRUE)
    mean(mu.1[z == 1] - mu.0[z == 1]) - tau
  else
    mean(mu.1[z == 0] - mu.0[z == 0]) - tau

  mu.1 <- mu.1 - omega

  y.0 <- rnorm(n, mu.0, sigma.y)
  y.1 <- rnorm(n, mu.1, sigma.y)
  y <- ifelse(z == 1, y.1, y.0)
  yc <- ifelse(z == 0, y.1, y.0)

  results <- cbind(x.m, y, y.1, y.0)
  #print(head(results))
  print(mean(mu.0))
  print(mean(mu.1))
  # print(head(y.0))
  # print(head(y.1))
  # print(head(y))
  # print(head(yc))
  # print("#####")
  # print(mu.1[z == 1] - mu.0[z == 1])
#   print(mean(y.1[z == 1] -  y.0[z == 1]))         #答え
#   print(mean(y.1 -  y.0))
  # print(mean(mu.1[z == 1] -  mu.0[z == 1]))       #4




}
