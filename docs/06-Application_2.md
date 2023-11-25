# Aplikace na reálných datech 2

V této části se podíváme podrobněji na klasifikaci pomocí postupu, při kterém nejprve data projektujeme na Reproducing Kernel Hilbert Space definovaný jádrem $K$ a následně pomocí koeficientů této projekce klasifikujeme data pomocí metody SVM. K získání koeficientů projekce se také využívá metoda SVM, avšak její analogie pro řešení regresních preblémů (SVM for regression). Tento postup je podrobněji popsán ve článku[^k], jehož autory jsou Muñoz a González.

[^k]: Muñoz, A. and González, J. (2010) *Representing functional data using support vector machines*, Pattern Recognition Letters, 31(6), pp. 511--516. [doi:10.1016/j.patrec.2009.07.014](https://www.sciencedirect.com/science/article/pii/S0167865509001913).

Autoři ve svém článku popisují mimo jiné i chování koeficientů z vyjádření pomocí jádrových funkcí (kernel expansion) a zlepšení stability koeficientů pomocí RKHS. Tuto vlastnost si nyní ilustrujme na analogickém příkladě, jako uvádějí autoři.


```r
library(fda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ddalpha)
library(patchwork)
library(e1071)
library(caret)
library(kernlab)
```

## *handwrit* data

Nejprve si načteme data, která představují souřadnice ručně psaných slov *fda*. Pro naše účely, kdy chceme porovnat dvě různá pozorování, si načtěme pouze první dva záznamy.


```r
# nacteni dat 
samples <- c(4, 5) + 10
hand <- handwrit[, samples, ]
```

Pro lepší orientaci si data nejprve vykresleme. Nejprve celkový záznam nápisu *fda*, následně jednotlivé souřadnice $x$ a $y$ proti času.



Nyní se podívejme, jak dopadnou vektory koeficientů $\boldsymbol \alpha_1$, $\boldsymbol \alpha_2$ z jádrového vyjádření získaného pomocí SVM pro regresi a také jak vypadají $\boldsymbol \lambda_1$, $\boldsymbol \lambda_2$ pro RKHS reprezentaci. Vektory $\boldsymbol \alpha_1$ a $\boldsymbol \alpha_2$ definují reprezentaci křivek

$$
c_l^*(\boldsymbol x) = \sum_{i = 1}^m \alpha_{il} K(\boldsymbol x_i, \boldsymbol x), \quad \forall \boldsymbol x \in \mathcal X,
$$

kde $\alpha_{il} \in \mathbb R, l= 1, 2.$. Podobně vektory $\boldsymbol \lambda_1$ a $\boldsymbol \lambda_2$ definují vyjádření křivek v bázi RKHS $\mathcal H_K$

$$
c_l^*(\boldsymbol x) = \sum_{j = 1}^d \lambda_{jl}^* \phi_j(\boldsymbol x),
$$

které můžeme odhadnout z dat pomocí 

$$
\hat\lambda_{jl}^* = \hat\lambda_{jl} \sum_{i = 1}^m \alpha_{il}\hat\phi_{ji}, \quad j = 1, 2, \dots, \hat d.
$$

Uvažme Gaussovské jádro s parametrem $\gamma = 0.5$ (analogicky jako ve výše zmíněném článku).


```r
gamma <- 0.5

# jadro a jadrova matice ... Gaussovske s parametrem gamma
Gauss.kernel <- function(x, y, gamma) {
  return(exp(-gamma * norm(c(x - y) |> t(), type = 'F')^2))
}

Kernel.RKHS <- function(x, gamma) {
  K <- matrix(NA, ncol = length(x), nrow = length(x))
  for(i in 1:nrow(K)) {
    for(j in 1:ncol(K)) {
      K[i, j] <- Gauss.kernel(x = x[i], y = x[j], gamma = gamma)
    }
  }
  return(K)
}
```

Nyní si spočítejme pomocí funkce `svm()` z knihovny `e1071` vektory koeficientů $\boldsymbol \alpha_1$ a $\boldsymbol \alpha_2$. Uvědomme si, že funkce `svm()` vrací koeficienty pouze pro podpůrné vektory, pro ostatní jsou tyto koeficienty nulové.


```r
# urceni koeficientu alpha z SVM
alpha <- matrix(0, nrow = length(handwritTime),
                ncol = length(samples)) # prazdny objekt

# model
for(i in 1:length(samples)) {
  df.svm <- data.frame(x = handwritTime,
                       y = hand[, i, 1])
  svm.RKHS <- svm(y ~ x, data = df.svm, 
                  kernel = 'radial',
                  type = 'eps-regression',
                  epsilon = 0.03,
                  gamma = gamma,
                  cost = 1, 
                  tolerance = 0.001,
                  shrinking = TRUE,
                  scale = TRUE)
  # alpha[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
  svm.RKHS <- train(y ~ x, data = df.svm,
                    method = 'svmRadial',
                    metric = "RMSE", 
                    preProcess = c('center', 'scale'),
                    # trControl = trainControl(
                    #   method = "repeatedcv", 
                    #   number = 5, 
                    #   repeats = 2, 
                    #   verboseIter = FALSE
                    #   )
                    trControl = trainControl(method = "none"),   
                    # Telling caret not to re-tune
                    tuneGrid = data.frame(sigma = 19, C = 1000)   
                    # Specifying the parameters
                    )
  alpha[svm.RKHS$finalModel@alphaindex, i] <- svm.RKHS$finalModel@alpha *
    svm.RKHS$finalModel@scaling$y.scale$`scaled:scale`
}

  # urceni alpha
  # b <- - svm.RKHS$rho # rho ... the negative intercept
  # betas <- svm.RKHS$coefs
  # alphas <- betas + b #/ Gauss.kernel()
```


```r
plot(handwritTime, svm.RKHS$finalModel@fitted * svm.RKHS$finalModel@scaling$y.scale$`scaled:scale`, type = 'l')
lines(handwritTime, hand[, i, 1], lwd = 3, col = 3)
```

<img src="06-Application_2_files/figure-html/unnamed-chunk-6-1.png" width="672" />



```r
# plot(alpha[, 1], type = 'n', ylab = expression(alpha[i]))
# segments(x0 = 1:1401, y0 = 0, y1 = alpha[, 1])

data.frame(alpha = c(alpha[, 1], alpha[, 2]),
           Time = c(handwritTime, handwritTime),
           Sample = rep(c('A', 'B'), each = length(handwritTime))) |>
  ggplot(aes(xend = Time, x = Time, y = 0, yend = alpha, colour = Sample)) + 
  geom_segment() + 
  theme_bw() + 
  labs(x = 'Time',
       y = expression(bold(alpha)[i])) 
```

<img src="06-Application_2_files/figure-html/unnamed-chunk-7-1.png" width="672" />

```r
#plot(solve(KKK) %*% svm.RKHS$fitted)
#plot(alpha[, 1] + solve(KKK) %*% rep(b, 1401))
```



```r
# spocitame matici K
K <- Kernel.RKHS(handwritTime, gamma = gamma)

# urcime vlastni cisla a vektory
Eig <- eigen(K)
eig.vals <- Eig$values
eig.vectors <- Eig$vectors
```


```r
# d
d.RKHS <- rankMM(K)

# urceni vektoru lambda
Lambda.RKHS <- matrix(NA, 
                      ncol = dim(hand)[2], 
                      nrow = d.RKHS) # vytvoreni prazdneho objektu

# vypocet reprezentace
Lambda.RKHS[, 1] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha[, 1]) *
  eig.vals[1:d.RKHS]
Lambda.RKHS[, 2] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha[, 2]) *
  eig.vals[1:d.RKHS]
```


```r
plot(Lambda.RKHS[, 2])
points(Lambda.RKHS[, 1], col = 'red')
```

<img src="06-Application_2_files/figure-html/unnamed-chunk-10-1.png" width="672" />


```r
plot(Lambda.RKHS[, 1] - Lambda.RKHS[, 2])#, ylim = c(-0.25, 0.05))
```

<img src="06-Application_2_files/figure-html/unnamed-chunk-11-1.png" width="672" />

## `growth` data

Podívejme se nyní ještě na tutéž analýzu, ale z pohledu dat `dataf.growth`.



```r
# nacteni dat 
samples <- c(54, 55) - 4
dataf <- dataf.growth()
t1 <- dataf$dataf[[samples[1]]]$args
t2 <- dataf$dataf[[samples[2]]]$args
y1 <- dataf$dataf[[samples[1]]]$vals
y2 <- dataf$dataf[[samples[2]]]$vals
data.gr <- data.frame(t = t1, H1 = y1, H2 = y2)
```

<img src="06-Application_2_files/figure-html/unnamed-chunk-13-1.png" width="672" />


```r
gamma <- 0.1

# jadro a jadrova matice ... Gaussovske s parametrem gamma
Gauss.kernel <- function(x, y, gamma) {
  return(exp(-gamma * norm(c(x - y) |> t(), type = 'F')^2))
}

Kernel.RKHS <- function(x, gamma) {
  K <- matrix(NA, ncol = length(x), nrow = length(x))
  for(i in 1:nrow(K)) {
    for(j in 1:ncol(K)) {
      K[i, j] <- Gauss.kernel(x = x[i], y = x[j], gamma = gamma)
    }
  }
  return(K)
}
```

Nyní si spočítejme pomocí funkce `svm()` z knihovny `e1071` vektory koeficientů $\boldsymbol \alpha_1$ a $\boldsymbol \alpha_2$. Uvědomme si, že funkce `svm()` vrací koeficienty pouze pro podpůrné vektory, pro ostatní jsou tyto koeficienty nulové.


```r
# urceni koeficientu alpha z SVM
alpha <- matrix(0, nrow = length(data.gr$t),
                ncol = 2) # prazdny objekt

# model
for(i in 1:2) {
  df.svm <- data.frame(x = data.gr$t,
                       y = data.gr[, i + 1])
  svm.RKHS <- svm(y ~ x, data = df.svm, 
                  kernel = 'radial',
                  type = 'eps-regression',
                  epsilon = 0.5,
                  gamma = gamma,
                  cost = 100000, 
                  tolerance = 0.001,
                  shrinking = FALSE,
                  scale = TRUE)
  # svm.RKHS <- train(y ~ x, data = df.svm,
  #                   method = 'svmLinear')
  # urceni alpha
  # b <- - svm.RKHS$rho # rho ... the negative intercept
  # betas <- svm.RKHS$coefs
  # alphas <- betas + b #/ Gauss.kernel()
  alpha[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
}
```


```r
# plot(alpha[, 1], type = 'n', ylab = expression(alpha[i]))
# segments(x0 = 1:1401, y0 = 0, y1 = alpha[, 1])

data.frame(alpha = c(alpha[, 1], alpha[, 2]),
           Time = c(data.gr$t, data.gr$t),
           Sample = rep(c('A', 'B'), each = length(data.gr$t))) |>
  ggplot(aes(xend = Time, x = Time, y = 0, yend = alpha, colour = Sample)) + 
  geom_segment() + 
  theme_bw() + 
  labs(x = 'Time',
       y = expression(bold(alpha)[i])) 
```

<img src="06-Application_2_files/figure-html/unnamed-chunk-16-1.png" width="672" />

```r
#plot(solve(KKK) %*% svm.RKHS$fitted)
#plot(alpha[, 1] + solve(KKK) %*% rep(b, 1401))
```



```r
# spocitame matici K
K <- Kernel.RKHS(data.gr$t, gamma = gamma)

# urcime vlastni cisla a vektory
Eig <- eigen(K)
eig.vals <- Eig$values
eig.vectors <- Eig$vectors
```


```r
# d
d.RKHS <- rankMM(K)

# urceni vektoru lambda
Lambda.RKHS <- matrix(NA, 
                      ncol = 2, 
                      nrow = d.RKHS) # vytvoreni prazdneho objektu

# vypocet reprezentace
Lambda.RKHS[, 1] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha[, 1]) *
  eig.vals[1:d.RKHS]
Lambda.RKHS[, 2] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha[, 2]) *
  eig.vals[1:d.RKHS]
```


```r
plot(Lambda.RKHS[, 2])
points(Lambda.RKHS[, 1], col = 'red')
```

<img src="06-Application_2_files/figure-html/unnamed-chunk-19-1.png" width="672" />


```r
plot(Lambda.RKHS[, 1] - Lambda.RKHS[, 2], ylim = c(-5, 3))
```

<img src="06-Application_2_files/figure-html/unnamed-chunk-20-1.png" width="672" />
