# Podpůrné materiály pro diplomovou práci

V této poslední kapitole jsou uvedeny zdrojové kódy pro vygenerování grafů a dalších případných materiálů, které jsou použity v diplomové práci. Jedná se především o ilustrativní grafy určitých vlastností a fenoménů spojených s funkcionálními daty. 

Kapitola je členěna do sekcí, které odpovídají jednotlivým kapitolám v diplomové práci. Všechny grafy jsou vytvořeny pomocí balíčku `ggplot2`, který poskytuje celou řadu grafických funkcionalit, pomocí kterých jsme (alespoň subjektivně) schopni dosáhnout podstatně lépe a profesionálněji vypadajících grafických výstupů v porovnání s klasickou grafikou v `R`.

Všechny grafy jsou uloženy pomocí funkce `ggsave()` ve formátu **pdf** nebo **tikz**, který umožňuje lepší kombinaci grafiky a symbolů v $\LaTeX$u.


```r
# nacteme potrebne balicky 
library(fda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ddalpha)
library(tidyverse)
library(patchwork)
library(tikzDevice)

set.seed(42)
options(tz = "UTC")

# nacteni dat
data <- read.delim2('phoneme.txt', header = T, sep = ',')

# zmenime dve promenne na typ factor
data <- data |> 
  mutate(g = factor(g),
         speaker = factor(speaker))

# numericke promenne prevedeme opravdu na numericke
data[, 2:257] <- as.numeric(data[, 2:257] |> as.matrix())

tr_vs_test <- str_split(data$speaker, '\\.') |> unlist()
tr_vs_test <- tr_vs_test[seq(1, length(tr_vs_test), by = 4)]
data$train <- ifelse(tr_vs_test == 'train', TRUE, FALSE)

# vybrane fonemy ke klasifikaci
phoneme_subset <- c('aa', 'ao')

# testovaci a trenovaci data
data_train <- data |> filter(train) |> filter(g %in% phoneme_subset)
data_test <- data |> filter(!train) |> filter(g %in% phoneme_subset)

# odstranime sloupce, ktere nenesou informaci o frekvenci a 
# transponujeme tak, aby ve sloupcich byly jednotlive zaznamy
X_train <- data_train[, -c(1, 258, 259, 260)] |> t()
X_test <- data_test[, -c(1, 258, 259, 260)] |> t()

# prejmenujeme radky a sloupce
rownames(X_train) <- 1:256
colnames(X_train) <- paste0('train', data_train$row.names)
rownames(X_test) <- 1:256
colnames(X_test) <- paste0('test', data_test$row.names)

# definujeme vektor fonemu
y_train <- data_train[, 258] |> factor(levels = phoneme_subset)
y_test <- data_test[, 258] |> factor(levels = phoneme_subset)
y <- c(y_train, y_test)
```


## Materiály pro Kapitolu 1

V této sekci uvedeme podpůrné grafy pro první kapitolu diplomové práce.

### Funkcionální průměr

Pro data `phoneme` spočítáme průměrný průběh log-periodogramů.


```r
t <- 1:256
rangeval <- range(t)
breaks <- t
norder <- 4

bbasis <- create.bspline.basis(rangeval = rangeval, 
                               norder = norder, 
                               breaks = breaks)

curv.Lfd <- int2Lfd(2) # penalizujeme 2. derivaci

# spojeni pozorovani do jedne matice
XX <- cbind(X_train, X_test) |> as.matrix()
XXaa <- XX[, y == phoneme_subset[1]]

lambda.vect <- 10^seq(from = 1, to = 3, length.out = 35) # vektor lambd
gcv <- rep(NA, length = length(lambda.vect)) # prazdny vektor pro ulozebi GCV

for(index in 1:length(lambda.vect)) {
  curv.Fdpar <- fdPar(bbasis, curv.Lfd, lambda.vect[index])
  BSmooth <- smooth.basis(t, XX, curv.Fdpar) # vyhlazeni
  gcv[index] <- mean(BSmooth$gcv) # prumer pres vsechny pozorovane krivky
}

GCV <- data.frame(
  lambda = round(log10(lambda.vect), 3),
  GCV = gcv
)

# najdeme hodnotu minima
lambda.opt <- lambda.vect[which.min(gcv)]

curv.fdPar <- fdPar(bbasis, curv.Lfd, lambda.opt)
BSmooth <- smooth.basis(t, XX, curv.fdPar)
XXfd <- BSmooth$fd

fdobjSmootheval <- eval.fd(fdobj = XXfd, evalarg = t)

## pouze pro aa
lambda.vect <- 10^seq(from = 1, to = 3, length.out = 35) # vektor lambd
gcv <- rep(NA, length = length(lambda.vect)) # prazdny vektor pro ulozebi GCV

for(index in 1:length(lambda.vect)) {
  curv.Fdpar <- fdPar(bbasis, curv.Lfd, lambda.vect[index])
  BSmoothaa <- smooth.basis(t, XXaa, curv.Fdpar) # vyhlazeni
  gcv[index] <- mean(BSmoothaa$gcv) # prumer pres vsechny pozorovane krivky
}

GCV <- data.frame(
  lambda = round(log10(lambda.vect), 3),
  GCV = gcv
)

# najdeme hodnotu minima
lambda.opt <- lambda.vect[which.min(gcv)]

curv.fdPar <- fdPar(bbasis, curv.Lfd, lambda.opt)
BSmoothaa <- smooth.basis(t, XXaa, curv.fdPar)
XXfdaa <- BSmoothaa$fd

fdobjSmoothevalaa <- eval.fd(fdobj = XXfdaa, evalarg = t)

# prumer 
meanfd <- mean.fd(XXfdaa)
fdmean <- eval.fd(fdobj = meanfd, evalarg = t)
```



```r
n <- dim(XX)[2]
DFsmooth <- data.frame(
  t = rep(t, n),
  time = factor(rep(1:n, each = length(t))),
  Smooth = c(fdobjSmootheval),
  Phoneme = rep(y, each = length(t))) |>
  filter(Phoneme == 'aa')

DFmean <- data.frame(
  t = rep(t, 2),
  Mean = c(fdmean, fdmean),
  Phoneme = factor(rep(phoneme_subset, each = length(t)),
                 levels = levels(y))
) |> filter(Phoneme == 'aa')

# tikz(file = "figures/DP_kap1_mean.tex", width = 4.2, height = 3.5)

DFsmooth |> 
  filter(time %in% as.character(1:100)) |>
  ggplot(aes(x = t, y = Smooth)) + 
  geom_line(aes(group = time), linewidth = 0.2, colour = 'deepskyblue2',
            alpha = 0.6) +
  theme_bw() +
  labs(x = 'Frekvence',
       y = 'Log-periodogram',
       colour = 'Phoneme') +
  scale_colour_discrete(labels = phoneme_subset) + 
  geom_line(data = DFmean, aes(x = t, y = Mean, 
                               group = Phoneme), 
            linewidth = 1, linetype = 'solid', colour = 'grey2') + 
  theme(legend.position = 'none')
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-3-1.png" alt="Vykreslení prvních 100 vyhlazených pozorovaných křivek. Černou čarou je zakreslen průměr." width="672" />
<p class="caption">(\#fig:unnamed-chunk-3)Vykreslení prvních 100 vyhlazených pozorovaných křivek. Černou čarou je zakreslen průměr.</p>
</div>

```r
# dev.off()
# ggsave("figures/DP_kap1_mean.pdf", width = 6, height = 5)
```

### Variance

Pro data `phoneme` spočítáme průběh varianční funkce.


```r
varfd <- var.fd(XXfdaa)
fdvar <- eval.bifd(t, t, varfd)
```


```r
dfs <- data.frame(
  time = t, 
  value = c(fdobjSmoothevalaa))
df <- data.frame(dfs, fdmean = fdmean, fdvar = diag(fdvar))

# tikz(file = "figures/DP_kap1_variance.tex", width = 6, height = 5)
# df <- df[seq(1, length(df$time), length = 1001), ]

ggplot(data = df, aes(x = time, y = fdvar)) +
  geom_line(color = 'deepskyblue2', linewidth = 0.8) +
  labs(x = 'Frekvence',
       y = 'Variance',
       colour = 'Phoneme') +
  theme_bw()
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-5-1.png" alt="Varianční funkce." width="672" />
<p class="caption">(\#fig:unnamed-chunk-5)Varianční funkce.</p>
</div>

```r
 # dev.off()
# ggsave("figures/DP_kap1_variance.tex", width = 4.2, height = 3.5, device = tikz)
```

### Kovariance a Korelace

Pro data `phoneme` spočítáme kovarianční a korelační funkce.


```r
fdcor <- cor.fd(t, XXfdaa) 

df <- merge(t, t)
df <- data.frame(df, fdcov = c(fdvar), fdcor = c(fdcor))
df <- df[seq(1, length(df$x), length = 68001), ]

# tikz(file = "figures/DP_kap1_cov.tex", width = 4, height = 4)

p1 <- ggplot(data = df, aes (x, y, z = fdcov)) +
  geom_raster(aes(fill = fdcov)) +
  geom_contour(colour = "white") +
  labs(x = 'Frekvence',
       y = 'Frekvence',
       fill = 'Kovariance') +
  coord_fixed(ratio = 1) + 
  theme_bw() + 
  theme(legend.position = 'bottom') + 
  scale_y_continuous(expand = c(0,0) + 0.01) + 
  scale_x_continuous(expand = c(0,0) + 0.01)
p1
```

<img src="14-Application_5_files/figure-html/unnamed-chunk-6-1.png" width="672" />

```r
# dev.off()
# ggsave("figures/DP_kap1_cov.tex", width = 6, height = 6, device = tikz)
# tikz(file = "figures/DP_kap1_cor.tex", width = 4, height = 4)
p2 <- ggplot(data = df, aes (x, y, z = fdcor)) +
  geom_raster(aes(fill = fdcor)) +
  geom_contour(colour = "white") +
  labs(x = 'Frekvence',
       y = 'Frekvence',
       fill = 'Korelace') +
  coord_fixed(ratio = 1) + 
  theme_bw() + 
  theme(legend.position = 'bottom') + 
  scale_y_continuous(expand = c(0,0) + 0.01) + 
  scale_x_continuous(expand = c(0,0) + 0.01)
p2
```

<img src="14-Application_5_files/figure-html/unnamed-chunk-6-2.png" width="672" />

```r
# dev.off()
# ggsave("figures/DP_kap1_cor.tex", width = 6, height = 6, device = tikz)
```


### B-splinová báze

Podívejme se na princip, jak se pomocí splinové báze dostaneme od diskrétních naměřených hodnot k funkcionálním datům.

Uvažujme pro přehlednost opět data `phoneme` a pouze malý počet bázových funkcí. Uvedeme tři obrázky, jeden se znázorněnými bázovými funkcemi, druhý s bázovými funkcemi přenásobenými vypočtenou hodnotou parametru a třetí výslednou křivku poskládanou sečtením jednotlivých přeškálovaných bázových funkcí.


```r
# definice barev 
cols7 <- c("#12DEE8", "#4ECBF3", "#127DE8", "#4C3CD3", "#4E65F3", "#4E9EF3", "#081D58")
cols5 <- c("#12DEE8",  "#1D89BC", "#4E9EF3", "#4C3CD3", "#081D58")
```


#### pro `norder = 2`


```r
df <- data.frame(x = 1:256, y = data[5, 2:257] |> c() |> unlist())
breaks <- quantile(df$x, probs = seq(0.1, 0.9, by = 0.2))
norder <- 2
rangeval <- range(df$x)

bbasis <- create.bspline.basis(rangeval, norder = norder, breaks = breaks)
BSmooth <- smooth.basis(df$x, df$y, bbasis)
```



```r
fdBSmootheval <- eval.fd(fdobj = BSmooth$fd, evalarg = df$x)
fdB <- eval.basis(basisobj = bbasis, evalarg = df$x)

basisdf1 <-  data.frame(bs = c(fdB), 
                        x = df$x, 
                        basis = rep(colnames(fdB), each = length(df$x)))
ebchan <- fdB * matrix(rep(BSmooth$fd$coefs, each = length(df$x)),
                       nrow = length(df$x))
basisdf2 <-  data.frame(bs = c(ebchan), 
                        x = df$x, 
                        basis = rep(colnames(fdB), each = length(df$x)))

library(RColorBrewer)

# tikz(file = "figures/DP_kap1_Bbasis_norder2.tex", width = 9, height = 3)

# samotna baze
p1 <- ggplot(data = basisdf1, aes(x = x, y = bs * 10, colour = basis)) +
  geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey2') +
  geom_line() +
  labs(x = 'Frekvence',
       y = 'B-splajnová báze',
       colour = 'Foném') +
  theme_classic() + 
  guides(colour = FALSE) +
  ylim(c(0, 22)) + 
  # scale_color_brewer(palette = 'Blues') +
  # scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[c(6,8,9,10,12)])
  scale_color_manual(values = cols5)

# prenasobena koeficienty
p2 <- ggplot(data = basisdf2, aes(x = x, y = bs, colour = basis)) +
  geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey2') +
  geom_line() + 
  labs(x = 'Frekvence',
       y = 'B-splajnová báze (škálovaná)',
       colour = 'Foném') +
  theme_classic() + 
  guides(colour = FALSE) + 
  ylim(c(0, 22)) + 
  # scale_color_brewer()
  # scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[c(6,8,9,10,12)])
  scale_color_manual(values = cols5)

# vyhlazena data
p3 <- ggplot(data = df, aes(x, y)) + 
  geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey2') +
  geom_point(colour = 'deepskyblue2', size = 0.8, alpha = 0.75) +
  geom_line(aes(y = fdBSmootheval)) +
  theme_classic() + 
  #guides (colour = FALSE) + 
  ylim(c(0, 22)) + 
  labs(x = 'Frekvence',
       y = 'Log-periodogram',
       colour = 'Foném')

(p1 | p2 | p3)
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-9-1.png" alt="B-spliny." width="1152" />
<p class="caption">(\#fig:unnamed-chunk-9)B-spliny.</p>
</div>

```r
# dev.off()
# ggsave('figures/DP_kap1_Bbasis_norder2.tex', device = tikz)
```


#### pro `norder = 4`


```r
df <- data.frame(x = 1:256, y = data[5, 2:257] |> c() |> unlist())
breaks <- quantile(df$x, probs = seq(0.1, 0.9, by = 0.2))
norder <- 4
rangeval <- range(df$x)

bbasis <- create.bspline.basis (rangeval, norder = norder, breaks = breaks)
BSmooth <- smooth.basis(df$x, df$y, bbasis)
```



```r
fdBSmootheval <- eval.fd(fdobj = BSmooth$fd, evalarg = df$x)
fdB <- eval.basis(basisobj = bbasis, evalarg = df$x)

basisdf1 <-  data.frame(bs = c(fdB), 
                        x = df$x, 
                        basis = rep(colnames(fdB), each = length(df$x)))
ebchan <- fdB * matrix(rep(BSmooth$fd$coefs, each = length(df$x)),
                       nrow = length(df$x))
basisdf2 <-  data.frame(bs = c(ebchan), 
                        x = df$x, 
                        basis = rep(colnames(fdB), each = length(df$x)))

# tikz(file = "figures/DP_kap1_Bbasis_norder4.tex", width = 9, height = 3)

# samotna baze
p1 <- ggplot(data = basisdf1, aes(x = x, y = bs * 10, colour = basis)) +
  geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey2') +
  geom_line() +
  labs(x = 'Frekvence',
       y = 'B-splajnová báze',
       colour = 'Foném') +
  theme_classic() + 
  guides(colour = FALSE) +
  ylim(c(0, 22)) + 
  # scale_color_brewer(palette = 'Blues') + 
  # scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[c(6,8,9,10,12)])
  scale_color_manual(values = cols7)

# prenasobena koeficienty
p2 <- ggplot(data = basisdf2, aes(x = x, y = bs, colour = basis)) +
  geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey2') +
  geom_line() + 
  labs(x = 'Frekvence',
       y = 'B-splajnová báze (škálovaná)',
       colour = 'Foném') +
  theme_classic() + 
  guides(colour = FALSE) + 
  ylim(c(0, 22)) + 
  # scale_color_brewer()
  # scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[c(6,8,9,10,12)])
  scale_color_manual(values = cols7)
  

# vyhlazena data
p3 <- ggplot(data = df, aes(x, y)) + 
  geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey2') +
  geom_point(colour = 'deepskyblue2', size = 0.8, alpha = 0.75) +
  geom_line(aes(y = fdBSmootheval)) +
  theme_classic() + 
  #guides (colour = FALSE) + 
  ylim(c(0, 22)) + 
  labs(x = 'Frekvence',
       y = 'Log-periodogram',
       colour = 'Foném')

(p1 | p2 | p3)
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-11-1.png" alt="B-spliny." width="1152" />
<p class="caption">(\#fig:unnamed-chunk-11)B-spliny.</p>
</div>

```r
# dev.off()
#ggsave('figures/DP_kap1_Bbasis_norder4.pdf')
```

### Fourierova báze

Podívejme se na princip, jak se pomocí Fourierovské báze dostaneme od diskrétních naměřených hodnot k funkcionálním datům.

Uvažujme pro přehlednost opět data `phoneme` a pouze malý počet bázových funkcí. Uvedeme tři obrázky, jeden se znázorněnými bázovými funkcemi, druhý s bázovými funkcemi přenásobenými vypočtenou hodnotou parametru a třetí výslednou křivku poskládanou sečtením jednotlivých přeškálovaných bázových funkcí.

#### pro `nbasis = 5`


```r
df <- data.frame(x = 1:256, y = data[5, 2:257] |> c() |> unlist())
nbasis <- 5
rangeval <- range(df$x)

fbasis <- create.fourier.basis(rangeval, nbasis = nbasis, period = 256)
FSmooth <- smooth.basis(df$x, df$y, fbasis)
```



```r
fdBSmootheval <- eval.fd(fdobj = FSmooth$fd, evalarg = df$x)
fdF <- eval.basis(basisobj = fbasis, evalarg = df$x)

basisdf1 <-  data.frame(bs = c(fdF), 
                        x = df$x, 
                        basis = rep(colnames(fdF), each = length(df$x)))
ebchan <- fdF * matrix(rep(FSmooth$fd$coefs, each = length(df$x)),
                       nrow = length(df$x))
basisdf2 <-  data.frame(bs = c(ebchan), 
                        x = df$x, 
                        basis = rep(colnames(fdF), each = length(df$x)))

# tikz(file = "figures/DP_kap1_Fbasis_nbasis5.tex", width = 9, height = 3)

# samotna baze
p1 <- ggplot(data = basisdf1, aes(x = x, y = bs, colour = basis)) +
  # geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey') +
  geom_line() +
  labs(x = 'Frekvence',
       y = 'Fourierova báze',
       colour = 'Foném') +
  theme_classic() + 
  guides(colour = FALSE) +
  #ylim(c(0, 22)) + 
  # scale_color_brewer(palette = 'Blues')
  scale_color_manual(values = cols5)

# prenasobena koeficienty
p2 <- ggplot(data = basisdf2, aes(x = x, y = bs, colour = basis)) +
  # geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey') +
  geom_line() + 
  labs(x = 'Frekvence',
       y = 'Fourierova báze (škálovaná)',
       colour = 'Foném') +
  theme_classic() + 
  guides(colour = FALSE) + 
  #ylim(c(0, 22)) + 
  # scale_color_brewer()
  scale_color_manual(values = cols5)

# vyhlazena data
p3 <- ggplot(data = df, aes(x, y)) + 
  # geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey') +
  geom_point(colour = 'deepskyblue2', size = 0.8, alpha = 0.75) +
  geom_line(aes(y = fdBSmootheval)) +
  theme_classic() + 
  #guides (colour = FALSE) + 
  ylim(c(0, 22)) + 
  labs(x = 'Frekvence',
       y = 'Log-periodogram',
       colour = 'Foném')

(p1 | p2 | p3)
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-13-1.png" alt="Fourierova baze." width="1152" />
<p class="caption">(\#fig:unnamed-chunk-13)Fourierova baze.</p>
</div>

```r
# dev.off()
#ggsave('figures/DP_kap1_Fbasis_nbasis5.pdf')
```



#### pro `nbasis = 7`


```r
df <- data.frame(x = 1:256, y = data[5, 2:257] |> c() |> unlist())
nbasis <- 7
rangeval <- range(df$x)

fbasis <- create.fourier.basis(rangeval, nbasis = nbasis, period = 256)
FSmooth <- smooth.basis(df$x, df$y, fbasis)
```



```r
fdBSmootheval <- eval.fd(fdobj = FSmooth$fd, evalarg = df$x)
fdF <- eval.basis(basisobj = fbasis, evalarg = df$x)

basisdf1 <-  data.frame(bs = c(fdF), 
                        x = df$x, 
                        basis = rep(colnames(fdF), each = length(df$x)))
ebchan <- fdF * matrix(rep(FSmooth$fd$coefs, each = length(df$x)),
                       nrow = length(df$x))
basisdf2 <-  data.frame(bs = c(ebchan), 
                        x = df$x, 
                        basis = rep(colnames(fdF), each = length(df$x)))

# tikz(file = "figures/DP_kap1_Fbasis_nbasis7.tex", width = 12, height = 4)

# samotna baze
p1 <- ggplot(data = basisdf1, aes(x = x, y = bs, colour = basis)) +
  # geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey') +
  geom_line() +
  labs(x = 'Frekvence [Hz]',
       y = 'Fourierova báze',
       colour = 'Foném') +
  theme_classic() + 
  guides(colour = FALSE) +
  #ylim(c(0, 22)) + 
  # scale_color_brewer(palette = 'Blues')
  scale_color_manual(values = cols7)

# prenasobena koeficienty
p2 <- ggplot(data = basisdf2, aes(x = x, y = bs, colour = basis)) +
  # geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey') +
  geom_line() + 
  labs(x = 'Frekvence [Hz]',
       y = 'Fourierova báze (škálovaná)',
       colour = 'Foném') +
  theme_classic() + 
  guides(colour = FALSE) + 
  #ylim(c(0, 22)) + 
  # scale_color_brewer()
  scale_color_manual(values = cols7)

# vyhlazena data
p3 <- ggplot(data = df, aes(x, y)) + 
  # geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey') +
  geom_point(colour = 'deepskyblue2', size = 0.8, alpha = 0.75) +
  geom_line(aes(y = fdBSmootheval)) +
  theme_classic() + 
  #guides (colour = FALSE) + 
  ylim(c(0, 22)) + 
  labs(x = 'Frekvence [Hz]',
       y = 'Log-periodogram',
       colour = 'Foném')

(p1 | p2 | p3)
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-15-1.png" alt="Fourierova baze." width="1152" />
<p class="caption">(\#fig:unnamed-chunk-15)Fourierova baze.</p>
</div>

```r
# dev.off()
#ggsave('figures/DP_kap1_Fbasis_nbasis7.pdf')
```

## Materiály pro Kapitolu 2

Ve druhé sekci se podíváme na materiály pro Kapitolu 2 diplomové práce.

Bude nás zajímat vliv vyhlazovacího parametru $\lambda$ na výslednou odhadnutou křivku z diskrétních dat. Dále se podíváme na funkcionální analýzu hlavních komponent. Nejprve se ale podívejme na vliv počtu bázových funkcí na výsledný odhad funkce.

### Počet bázových funkcí a výsledný odhad



```r
df <- data.frame(x = 1:256, y = data[5, 2:257] |> c() |> unlist())
norder <- 4
rangeval <- range(df$x)

bbasis1 <- create.bspline.basis (rangeval, norder = norder, nbasis = 5)
BSmooth1 <- smooth.basis(df$x, df$y, bbasis1)
fdBSmootheval1 <- eval.fd(fdobj = BSmooth1$fd, evalarg = df$x)

bbasis2 <- create.bspline.basis (rangeval, norder = norder, nbasis = 15)
BSmooth2 <- smooth.basis(df$x, df$y, bbasis2)
fdBSmootheval2 <- eval.fd(fdobj = BSmooth2$fd, evalarg = df$x)

bbasis3 <- create.bspline.basis (rangeval, norder = norder, nbasis = 25)
BSmooth3 <- smooth.basis(df$x, df$y, bbasis3)
fdBSmootheval3 <- eval.fd(fdobj = BSmooth3$fd, evalarg = df$x)

# 10 bazovych funkci
p1 <- ggplot(data = df, aes(x, y)) + 
  # geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey') +
  geom_point(colour = 'deepskyblue2', size = 0.8, alpha = 0.75) +
  geom_line(aes(y = fdBSmootheval1), linewidth = 0.7) +
  theme_classic() + 
  ylim(c(0, 22)) + 
  labs(x = 'Frekvence',
       y = 'Log-periodogram',
       colour = 'Foném')

# 10 bazovych funkci
p2 <- ggplot(data = df, aes(x, y)) + 
  # geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey') +
  geom_point(colour = 'deepskyblue2', size = 0.8, alpha = 0.75) +
  geom_line(aes(y = fdBSmootheval2), linewidth = 0.7) +
  theme_classic() + 
  ylim(c(0, 22)) + 
  labs(x = 'Frekvence',
       y = 'Log-periodogram',
       colour = 'Foném')

# 10 bazovych funkci
p3 <- ggplot(data = df, aes(x, y)) + 
  # geom_vline(xintercept = breaks, linetype = "dotted", linewidth = 0.1, colour = 'grey') +
  geom_point(colour = 'deepskyblue2', size = 0.8, alpha = 0.75) +
  geom_line(aes(y = fdBSmootheval3), linewidth = 0.7) +
  theme_classic() + 
  ylim(c(0, 22)) + 
  labs(x = 'Frekvence',
       y = 'Log-periodogram',
       colour = 'Foném')

# tikz(file = "figures/DP_kap2_differentnbasis.tex", width = 12, height = 4)

(p1 | p2 | p3)
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-16-1.png" alt="Počet bázových funkcí a výsledný odhad." width="1152" />
<p class="caption">(\#fig:unnamed-chunk-16)Počet bázových funkcí a výsledný odhad.</p>
</div>

```r
# dev.off()
ggsave('figures/DP_kap2_differentnbasis.tex', width = 9, height = 3, device = tikz)
```



### Volba $\lambda$

Začněme volbou vyhlazovacího parametru $\lambda > 0$. S rostoucí hodnotou $\lambda$ dáváme v penalizované sumě čtverců 
$$
SS_{pen} = (\boldsymbol y - \boldsymbol B \boldsymbol c)^\top (\boldsymbol y - \boldsymbol B \boldsymbol c) + \lambda \boldsymbol c^\top \boldsymbol R \boldsymbol c
$$
větší váhu penalizačnímu členu, tedy dostaneme více penalizované, více hladké křivky blížící se lineární funkci.
Vykreslíme si obrázky, ve kterých bude zřejmé, jak se s měnící se hodnotou $\lambda$ mění výsledná vyhlazená křivka.

Ke znázornění tohoto chování použijeme data `phoneme` z jedné z předchozích kapitol. Vybereme jedno zajímavé pozorování a ukážeme na něm toto chování.

Za uzly bereme celý vektor frekvencí (1 až 256 Hz), standardně uvažujeme kubické spliny, proto volíme (implicitní volba v `R`) `norder = 4`.
Budeme penalizovat druhou derivaci funkcí.


```r
t <- 1:256
rangeval <- range(t)
breaks <- t
norder <- 4

bbasis <- create.bspline.basis(rangeval = rangeval, 
                               norder = norder, 
                               breaks = breaks)

curv.Lfd <- int2Lfd(2) # penalizujeme 2. derivaci

# spojeni pozorovani do jedne matice
XX <- cbind(X_train, X_test) |> as.matrix()
```

Zvolme nyní nějakých 6 hodnot pro vyhlazovací parametr $\lambda$ a spočítejme vyhlazené křivky pro jeden vybraný záznam.


```r
lambdas <- c(0.01, 0.1, 50, 500, 10000, 1000000) # vektor lambd

tt <- seq(min(t), max(t), length = 1001)

# objekt, do ktereho ulozime hodnoty
res_plot <- matrix(NA, ncol = length(lambdas), nrow = length(tt))

for(i in 1:length(lambdas)) {
  curv.fdPar <- fdPar(bbasis, curv.Lfd, lambdas[i])
  BSmooth <- smooth.basis(t, XX, curv.fdPar)
  XXfd <- BSmooth$fd
  
  fdobjSmootheval <- eval.fd(fdobj = XXfd, evalarg = tt)[, 1]
  res_plot[, i] <- fdobjSmootheval
}
```


```r
options(scipen = 999)
library(scales)

lam_labs <- paste0('$\\lambda = ', lambdas, "$")
names(lam_labs) <- lambdas

# tikz(file = "figures/DP_kap2_lambdas.tex", width = 9, height = 6)

data.frame(time = rep(tt, length(lambdas)),
           value = c(res_plot),
           lambda = rep(lambdas, each = length(tt))) |>
  # mutate(lambda = factor(lambda)) |>
  ggplot(aes(x = time, y = value)) + 
  geom_point(data = data.frame(time = rep(t, length(lambdas)),
                               value = rep(c(data[5, 2:257]) |> unlist(), length(lambdas)),
                               lambda = rep(lambdas, each = length(t))) ,
             alpha = 0.5, size = 0.75, colour = "deepskyblue2") + 
  geom_line(linewidth = 0.7, colour = "grey2") + 
  facet_wrap(~lambda, ncol = 3, nrow = 2, labeller = labeller(lambda = lam_labs)) +
  theme_bw() + 
  theme(legend.position = 'none') + 
  labs(x = 'Frekvence',
       y = 'Log-periodogram',
       colour = 'Foném') 
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-19-1.png" alt="Log-periodogram vybraného fonému pro zvolené hodnoty vyhlazovacího parametru." width="864" />
<p class="caption">(\#fig:unnamed-chunk-19)Log-periodogram vybraného fonému pro zvolené hodnoty vyhlazovacího parametru.</p>
</div>

```r
  # scale_color_brewer()
# dev.off()
# ggsave('figures/DP_kap2_lambdas.pdf')
```

### Vyhlazení s optimální $\lambda$

V Kapitole Aplikace na reálných datech 2 jsme zjistili optimální hodnotu vyhlazovacího parametru. Tu nyní použijeme.



```r
lambdas <- c(175.75)

tt <- seq(min(t), max(t), length = 1001)

# objekt, do ktereho ulozime hodnoty
res_plot <- matrix(NA, ncol = length(lambdas), nrow = length(tt))

for(i in 1:length(lambdas)) {
  curv.fdPar <- fdPar(bbasis, curv.Lfd, lambdas[i])
  BSmooth <- smooth.basis(t, XX, curv.fdPar)
  XXfd <- BSmooth$fd
  
  fdobjSmootheval <- eval.fd(fdobj = XXfd, evalarg = tt)[, 1]
  res_plot[, i] <- fdobjSmootheval
}
```


```r
options(scipen = 999)
library(scales)

lam_labs <- paste0('$\\lambda = ', lambdas, "$")
names(lam_labs) <- lambdas

# tikz(file = "figures/DP_kap2_optimal_lambda.tex", width = 6, height = 4)

data.frame(time = rep(tt, length(lambdas)),
           value = c(res_plot),
           lambda = rep(lambdas, each = length(tt))) |>
  # mutate(lambda = factor(lambda)) |>
  ggplot(aes(x = time, y = value)) + 
  geom_point(data = data.frame(time = rep(t, length(lambdas)),
                               value = rep(c(data[5, 2:257]) |> unlist(), length(lambdas)),
                               lambda = rep(lambdas, each = length(t))) ,
             alpha = 0.5, size = 0.75, colour = "deepskyblue2") + 
  geom_line(linewidth = 0.7, colour = "grey2") + 
  # facet_wrap(~lambda, ncol = 3, nrow = 2, labeller = labeller(lambda = lam_labs)) +
  theme_bw() + 
  theme(legend.position = 'none') + 
  labs(x = 'Frekvence',
       y = 'Log-periodogram',
       colour = 'Foném') 
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-21-1.png" alt="Log-periodogram vybraného fonému pro zvolené hodnoty vyhlazovacího parametru." width="864" />
<p class="caption">(\#fig:unnamed-chunk-21)Log-periodogram vybraného fonému pro zvolené hodnoty vyhlazovacího parametru.</p>
</div>

```r
  # scale_color_brewer()
# dev.off()
# ggsave('figures/DP_kap2_lambdas.pdf')
```


### Funkcionální PCA

Najdeme vhodnou hodnotu vyhlazovacího parametru $\lambda > 0$ pomocí $GCV(\lambda)$, tedy pomocí zobecněné cross--validace.


```r
t <- 1:256
rangeval <- range(t)
breaks <- t
norder <- 4

bbasis <- create.bspline.basis(rangeval = rangeval, 
                               norder = norder, 
                               breaks = breaks)

curv.Lfd <- int2Lfd(2) # penalizujeme 2. derivaci
```



```r
# spojeni pozorovani do jedne matice
XX <- cbind(X_train, X_test) |> as.matrix()

lambda.vect <- 10^seq(from = 1, to = 3, length.out = 35) # vektor lambd
gcv <- rep(NA, length = length(lambda.vect)) # prazdny vektor pro ulozebi GCV

for(index in 1:length(lambda.vect)) {
  curv.Fdpar <- fdPar(bbasis, curv.Lfd, lambda.vect[index])
  BSmooth <- smooth.basis(t, XX, curv.Fdpar) # vyhlazeni
  gcv[index] <- mean(BSmooth$gcv) # prumer pres vsechny pozorovane krivky
}

GCV <- data.frame(
  lambda = round(log10(lambda.vect), 3),
  GCV = gcv
)

# najdeme hodnotu minima
lambda.opt <- lambda.vect[which.min(gcv)]
```

S touto optimální volbou vyhlazovacího parametru $\lambda$ nyní vyhladíme všechny funkce.


```r
curv.fdPar <- fdPar(bbasis, curv.Lfd, lambda.opt)
BSmooth <- smooth.basis(t, XX, curv.fdPar)
XXfd <- BSmooth$fd

fdobjSmootheval <- eval.fd(fdobj = XXfd, evalarg = t)
```

Proveďme tedy nejprve funkcionální analýzu hlavních komponent a určeme počet $p$.


```r
# analyza hlavnich komponent
data.PCA <- pca.fd(XXfd, nharm = 20) # nharm - maximalni pocet HK
nharm <- which(cumsum(data.PCA$varprop) >= 0.9)[1] # urceni p

# data.PCA <- pca.fd(XXfd, nharm = 20) 
data.PCA.train <- as.data.frame(data.PCA$scores) # skore prvnich p HK
data.PCA.train$Y <- factor(y) # prislusnost do trid
```

V tomto konkrétním případě jsme za počet hlavních komponent vzali $p$ = 9, které dohromady vysvětlují 90.47 % variability v datech.
První hlavní komponenta potom vysvětluje 44.79 % a druhá 13.37 % variability.
Graficky si můžeme zobrazit hodnoty skórů prvních dvou hlavních komponent, barevně odlišených podle příslušnosti do klasifikační třídy.


```r
p1 <- data.PCA.train |> ggplot(aes(x = V1, y = V2, colour = Y)) +
  geom_point(size = 1.2, alpha = 0.75) + 
  labs(x = paste('1. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[1], 2), '\\%)'),
       y = paste('2. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[2], 2), '\\%)'),
       colour = 'Foném') +
  # scale_colour_discrete(labels = phoneme_subset) +
  theme_bw() + 
  theme(legend.position = 'none') + 
  lims(x = c(-70, 62), y = c(-70, 62)) +
  scale_colour_manual(values = c('tomato', 'deepskyblue2'))

p2 <- data.PCA.train |> ggplot(aes(x = V1, y = V3, colour = Y)) +
  geom_point(size = 1.2, alpha = 0.75) + 
  labs(x = paste('1. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[1], 2), '\\%)'),
       y = paste('3. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[3], 2), '\\%)'),
       colour = 'Foném') +
  # scale_colour_discrete(labels = phoneme_subset) +
  theme_bw() + 
  theme(legend.position = 'none') + 
  lims(x = c(-70, 62), y = c(-70, 62)) +
  scale_colour_manual(values = c('tomato', 'deepskyblue2'))

p3 <- data.PCA.train |> ggplot(aes(x = V2, y = V3, colour = Y)) +
  geom_point(size = 1.2, alpha = 0.75) + 
  labs(x = paste('2. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[2], 2), '\\%)'),
       y = paste('3. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[3], 2), '\\%)'),
       colour = 'Foném') +
  # scale_colour_discrete(labels = phoneme_subset) +
  theme_bw() + 
  lims(x = c(-70, 62), y = c(-70, 62)) +
  theme(legend.position = c(0.84, 0.84)) +
  scale_colour_manual(values = c('tomato', 'deepskyblue2'))

(p1|p2|p3)
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-26-1.png" alt="Hodnoty skórů prvních tří hlavních komponent pro trénovací data. Barevně jsou odlišeny body podle příslušnosti do klasifikační třídy." width="1152" />
<p class="caption">(\#fig:unnamed-chunk-26)Hodnoty skórů prvních tří hlavních komponent pro trénovací data. Barevně jsou odlišeny body podle příslušnosti do klasifikační třídy.</p>
</div>

```r
# tikz(file = "figures/kap2_PCA_scores1.tex", width = 3, height = 3)
# p1
# dev.off()
# tikz(file = "figures/kap2_PCA_scores2.tex", width = 3, height = 3)
# p2
# dev.off()
# tikz(file = "figures/kap2_PCA_scores3.tex", width = 3, height = 3)
# p3
# dev.off()
# ggsave("figures/kap2_PCA_scores.tex", device = tikz, width = 12, height = 4)
```

Podívejme se ještě na 3D graf skórů prvních třech hlavních komponent.


```r
# 3D plot 
library(plotly)
plot_ly(data = data.PCA.train,
        x = ~V1, y = ~V2, z = ~V3,
        color = ~Y, colors = c('tomato', 'deepskyblue2'),
        type="scatter3d", mode="markers", marker = list(size = 2.5)) %>% 
  layout(scene = list(xaxis = list(title = '1. hlavní komponenta'),
                     yaxis = list(title = '2. hlavní komponenta'),
                     zaxis = list(title = '3. hlavní komponenta')))
```

<div class="figure">

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-30e4409849e39179307f" style="width:672px;height:480px;"></div>
<script type="application/json" data-for="htmlwidget-30e4409849e39179307f">{"x":{"visdat":{"5dec61f2109f":["function () ","plotlyVisDat"]},"cur_data":"5dec61f2109f","attrs":{"5dec61f2109f":{"x":{},"y":{},"z":{},"mode":"markers","marker":{"size":2.5},"color":{},"colors":["tomato","deepskyblue2"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter3d"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"scene":{"xaxis":{"title":"1. hlavní komponenta"},"yaxis":{"title":"2. hlavní komponenta"},"zaxis":{"title":"3. hlavní komponenta"}},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-17.851886079165403,-13.415144453970905,24.023982470826997,21.234230199128238,10.846849109852061,4.4430911109400739,4.4007922081057593,55.772742658764791,30.89011054184958,57.547087174462817,43.246166716602445,2.8614801053436505,23.301513892159754,-3.6866963725526629,-5.9973030828747378,5.1377342202016818,31.458941015505406,20.572079704033037,0.22827344879490499,29.198568754638259,17.579998914256123,23.406128701065651,-2.2447301235355752,-25.191351368142314,-15.875851115779094,59.853463686174869,21.861305597334223,-22.736813923585174,14.77336785744723,34.756152047878636,44.940252538740012,0.81500189777827559,41.655878686102533,28.551303296595002,13.166146847876755,23.85381506824444,25.180516218760815,6.2170221963516061,26.588975536190194,-4.1740904285215015,16.475889763223826,13.483547987152507,-3.6488669556747046,33.495228886847293,3.0016185635150694,22.514947415774635,15.986862049380383,-19.202646668331589,-20.701187117902553,53.312709599233607,28.479757956454165,24.10131949679981,10.74727639928423,-15.191769057060393,0.87072635949440635,28.458742775669986,1.8332556455110811,9.7666267813445593,15.865195204926225,11.054686534278559,-3.1240448325716463,28.7054270841011,30.559520803076722,34.941715034893164,22.229026288252037,6.2020226370099767,7.9637418569058696,24.586499983612594,2.7734152683892863,33.468540334411188,49.524173400093112,11.506946998911754,8.1820989770734602,28.628554708185348,7.2326102654684092,28.171617307274548,14.554278245977757,7.1978270744596538,-19.942766579402921,46.913037368132613,34.679285056833059,16.624144442030083,12.142781508286534,12.557425182070345,-15.181462756278659,15.117166454466842,33.072405064488805,10.146256461165427,22.406954705565639,9.1890010146749042,33.114674900869502,23.154567833153273,-8.2354338945202006,-8.1452613691940883,35.919623388923569,2.3529427065521276,13.495736151643763,10.884352860026615,23.595222080079182,10.181324298545562,13.160595839312505,13.593076994267502,2.4371317254412035,-10.172045623886239,-4.8525084007341555,-28.463311820279841,23.726707286401677,21.980879747824137,37.777337581425961,16.499780605136674,25.311727939043344,41.124967294611416,7.1956864286078925,25.755120586820262,23.450988026622952,13.67792788293756,2.1636312652895024,51.598061543764828,22.166386355652559,49.33001756013492,47.230400555016793,1.2186516421140541,11.276859819478807,13.214297323882008,8.1240739897042982,15.740780132961721,21.304310328493436,-17.642325617893881,-27.168974980534188,16.169955787434848,37.236833679914625,26.407221540547614,49.342107410721042,10.454318190560668,11.258969716556214,-33.030739558992003,13.027910413557802,-1.1545621170370117,22.620264810212081,5.1352112409691077,17.821623790389811,4.2334558083392464,-2.1120430421938452,22.293774097538755,4.8131033703562842,10.325790474086848,-8.6558843649010928,-2.8718218167807068,-31.03747913097946,13.25565961169589,-9.998712752961497,20.324887873236641,32.174902257295614,44.630321349414636,38.792760182243995,-11.125429419380845,-27.624883594601428,47.071148229962795,21.679274076367307,-5.7604654408918634,45.030089651725703,12.925659389257797,1.561748299289534,-28.606593590797392,19.802282376426405,-13.238320821070388,-9.316027026592165,29.961001310145381,-6.7694902591450399,20.806427972171218,-0.14151398926454198,5.5876835729464389,14.295300547329507,5.5337604225287063,3.3370724543823935,11.340188502353696,17.575982698254109,31.047797042604277,32.662838642465772,21.606384861602574,-27.304619460049963,-0.029835524676423138,-21.7725863927521,23.734451402458557,16.655567608642986,10.976990777765971,21.151273993122892,29.614446864358921,-0.88778933505031599,-4.7516991198918674,-4.8489458583831304,41.880774028445408,48.912043009669389,17.600897400668657,19.8852687470192,-12.9265371113648,0.11054612743181136,17.583937510132625,5.1845714159771932,4.2154115205617515,-9.9125399891049479,-6.7315730931728739,-11.061169725221811,-0.11288681760998608,14.611232487363239,32.361754533426911,-2.5820496526936392,-2.4679230447179594,9.5726894543353165,41.911695827050956,17.098282607028807,-4.8227365513646667,3.0552107378105573,17.302672879457113,-15.413725827311028,-2.9346655309262673,25.359422376855569,18.994384673077182,28.027726963141887,-18.165190692216306,21.471369863391605,8.0502382134576145,-1.0949457716622124,8.4273153238072691,-7.1254550674694768,14.534784850558157,8.2951984919461506,-9.5676593704513913,-3.3880897002819377,-3.742668796281563,14.431766759482388,-29.925258038082344,-36.668701000832677,-43.132981985700077,20.134133232315683,-2.6235051970877836,60.120691254503903,2.8556596871902125,24.234874638846712,-9.6721220944332433,47.204821505469639,19.439961347636942,19.671441948146924,21.602852530149057,21.046152551334078,11.053237120116883,6.1081579395842542,15.525148000876177,-0.043023337744299561,40.180388760146883,0.59405893555051248,-13.030311038559153,50.356382743111048,-2.0588813178974732,2.7513445014774431,26.411302725425394,13.081838994168114,8.5478393220756246,-9.6573096614979743,12.576025415340037,17.471790968717166,18.550207809616907,28.000368755679862,8.1966760270303141,13.852813072075072,15.025809805152573,3.3914823487407562,15.590139486174445,11.208947616062293,10.299829477308613,-0.98234636574317979,-7.2206486889177928,-11.580902270883401,22.865471198688191,17.619710095220682,31.990148153330388,13.240312218306824,1.9051126412153991,27.272953128931206,43.913442893581617,28.158479707263968,-19.38286120886022,21.364184245763511,17.574690919668754,25.590689786154766,21.31613126906058,1.4386434006377529,0.76496987744981171,6.6114617460164551,6.3197285867720536,25.601844827682431,13.725236992229764,-2.7017537944588104,24.244067114599453,18.269576834123164,17.521266254912941,24.690179155623905,11.198335220820004,10.482651161127309,10.665929404914328,29.565959581236065,14.807580832345089,-12.725731774929148,30.882741386213326,22.965040846517159,2.6296961449721703,-19.103261309356188,23.208069587520633,28.485641096304466,-10.166620571295717,-16.278914417371304,4.5652537475632418,5.2003515118058958,45.64431901757419,36.996665898818819,42.387322286350972,10.577876327245914,8.3421054032896045,26.422671543266599,23.450801362839524,38.144188718381173,22.702020779616738,16.089395338904438,1.0862808759131073,17.544630677938038,0.75830723374730058,31.333610920183709,-5.1599804787170225,18.778165004800261,10.294124108511729,2.2154982709223936,15.972248005043372,21.22919570066675,23.262170685416908,27.44741413692946,17.87305366932204,0.10163500450764157,-9.3292329885391521,34.939381091307759,26.691281136849938,3.5595824952138018,23.77246360293141,43.88242803076546,-2.0051648550816714,37.941654613582799,20.791139680109293,26.847002209973294,35.799481699658571,29.636192093268459,43.878081211856212,12.794776090251817,9.737451968365832,13.854433202648934,-27.532465110495796,9.780961426571924,-16.68186094263864,-11.824616644299319,26.303040204374206,27.37170653872678,30.912026133991535,2.4500239793609468,-14.732037088872165,10.024068065841544,3.0015721773897366,34.690058664959693,-7.375776926619527,4.9368249255160874,38.104541631778403,11.961599491500982,31.442336410549423,14.744396097784056,3.2910249347242826,23.6728111491459,19.703549029828871,8.2322028617912437,4.9961377393889288,23.98451878726544,5.0108459703971198,9.7455824807546581,12.184889357773054,36.722489773271199,57.722435975119957,-10.948298792995931,-21.768633942915862,20.624015223759962,-0.78134362667869817,-5.9750176327689744,19.267436065551482,31.629554275171611,22.792804526103346,7.5279836711237467,-5.6720378009035377,-11.64293084569737,-25.833446032145314,-29.446813316898204,21.8052390943202,11.22393176224324,-5.946185361459202,-12.062028513581717,-29.895690393293012,31.110022080686914,17.750568160852342,11.712286269161309,9.2703184702040069,29.774220375722663,41.39257299606021,-16.739436870485608,39.123032519902686,8.9753020436552227,20.028598490868944,26.36521144295541,8.4712227490967784,1.7041816025761969,51.790992438499636,40.963699156550938,24.329172024912854,16.812301565406916,-12.094300990310723,45.622855148209119,25.189301665905056,21.002239590714851,39.647191126304563,-1.5736382499595094,27.662475351637664,-6.3065281456492794,29.611936116975571,13.089187189005687,11.691899389957332,16.166020880472171,-13.112871492173953,32.829016289226949,-28.200494317852527,-12.635663569849028,-37.1600558335155,-18.51559701751059,-12.926017905063627,-19.356881558746608,14.72848975470602,-1.3562835016133943,29.629005112429279,-0.88105282448103206,-18.842579216565511,-26.318663112671953,29.753766409809735,2.2313540203643711,-4.9207521613985179,16.443547942067273,25.419067969650758,8.7668126955663723,-4.7045530539378184,8.732800652611127,1.4752839662737383,-8.1717931919773417,-12.145528135432595,-0.25663159355685783,15.669654828458814,26.387972014984019,-6.3592067811082513,-13.679569635941766,-13.976781808259426,-27.061533300419281,25.317601872580731,14.091571310012164,13.92518396431934,-0.29504758687116645,7.6223575719135566,16.960328446075316,20.80516082630453,-2.2851261049992742,-19.199092700754239,-28.260158775343093,-46.883767317640931,3.3251896045648088,-3.9653842470460061,3.9415316226628727,-10.235422498980681,33.238614714790472,35.13235345054585,23.35380462588186,28.396146586870742,-15.41426700528149,-15.663593097641369,9.041859366861777,-17.667257295372092,5.6326549551125122,-11.951880988826092,25.385918667764788,11.856238695427466,28.741191722678298,15.361493664422595,28.143438650986425,14.48580222583681,0.44967755357839057,-3.9500015203683736,4.428612002746843,2.0268071415776525,-25.490041558531061,-11.952352760189596,17.980222390088397,-15.235735735536815,-17.910506263791063,-23.864398777038069,-34.207170304225038,3.7831102624806627,-15.712935290347243,-15.348626299168902,-19.606674854234811,11.087021473789699,36.849226951657542,43.11059219940767,26.253958734085735,-0.71804078991959708,-4.1666154671864266,20.259295316456317,3.2519248979488253,37.040344231485257,-1.7000681131809494,15.267887271339816,-5.8991296040053633,-3.4052420622648842,-30.766799212114826,42.671595215722867,43.97601508403946,18.697424400307337,21.424286722172997,11.932943147047638,-5.0491762824402082,20.738726872176375,1.2924510346767357,53.496221803170286,14.612769486533137,13.401474324754638,53.941670654841424,32.967631801958177,35.562222238362224,1.1877420771612091,-3.4539535918258673,19.272967424232917,23.543103669485937,28.422969870596393,-7.3958705844275654,29.131167576635811,6.7373321410485314,36.76793056196496,6.1436924078373636,-12.46955866097484,33.112327304163351,13.073015887307202,32.454272764622232,38.006212324840718,42.612954472036392,4.9454489771511732,20.519617071756088,8.5380584570687255,42.442106636013627,60.518946399425396,18.772074452134571,-7.9520841654287153,-10.865259460715027,-15.025607722049067,1.1927546105274609,1.2042672403869905,13.909952694200918,37.576610629553919,27.944609257817238,-2.8507978280323161,20.22699379670394,7.5973764727416322,7.1283146157636077,16.228392737377867,5.7235995855318267,-1.7036242117734257,-30.980240681813218,-34.349605510776513,-37.49445037064833,3.6227932642009049,-4.9913320253299602,-26.297988101403458,11.236092441912964,23.576435656970645,-16.442537154203844,18.154111090722864,36.284011492550164,5.0482876763107409,-25.352493531452172,15.097684743879935,17.61494039006465,-0.96753642967147879,10.742177707553049,10.348460020020735,5.5357056677125822,-9.6415289314428012,-1.0412800910619413,-8.9095877983256742,36.296181305447085,-9.1659140955698017,15.83148037033806,18.258163186786426,17.185998544901302,20.187162677893955,1.4814413993207836,15.926421720467276,15.269470765406332,41.510511447040997,15.876543620337554,1.2125979520943224,13.820167581927716,2.4095275359402786,-14.113505909684454,19.712382785375063,31.242156602513795,-4.397445910219262,-5.7479138959044649,-11.55214631336122,-6.7790773558772068,-13.728877814273337,-25.026206176545884,-2.4365281514120189,-6.8793749597808631,-1.7753866330936705,30.621641695878164,36.524910816602258,35.370419192113275,-6.6537428471209461,27.803087951418586,4.4754867883255516,3.0464054437489625,18.023925529037587,-0.91437716051662643,-8.306135922615077,-12.080657824194766,21.37399969449524,15.426517528124181,-21.664209640079516,-8.280409444632495,-16.041772404851606,-13.667797539794046,-9.4324259094314407,40.161089058450571,-9.8473709737689585,20.213533689269365,6.1664387754678849,47.023841122739775,32.360324758040996,32.707335495228072,18.876470936778961,8.5454948985765924,12.447875437989829,26.003519956590445,26.665279167051089,21.954412805190159,13.460665429823985,28.106840767170876,37.840913675822335,-31.045611854064482,-36.030558297138917,-41.51673020650815,21.13555615356665,-2.1887305557709578,10.400901716762469,1.0372745046345979,2.8005268880599266,10.754472030573828,9.8026193143474369,-7.5371182124490552,-11.298376103977889,-21.76092164311309,2.0138421035054357,-4.2036489466123399,19.286334869866565,-8.3534801161815952,6.3266024146932889,-12.021199494155672,-36.214286058332306,8.5597850006661496,3.1300255513726993,3.3830191829234599,-1.8849945406219102,11.544949402689969,16.89773359329077,18.974026330552078,13.326941386422643,2.0065209469304013,-18.139040346806564,30.144280154281208,43.642142956732684,-20.199212941069202,-24.378247724003575,-25.458430111516726,-13.289967538672332,24.096782249036345,31.361547423543335,12.147177985862561,-28.508590084248841,-22.08921712729224],"y":[-9.7470919088419787,2.3913408117606281,17.28041332721352,2.6526736863790656,-27.210479819881861,-21.482304226677208,-18.49794308865523,-10.497094054321682,2.0515798290093969,19.001792908995782,19.643058992982713,-8.0225459833103354,2.1335773655261328,-10.510039484237067,4.3146646920272547,-8.40452750275022,-0.52048292139669172,-0.19095604665192437,3.2159687171666662,1.7520388834298757,6.0771767851096694,-1.2784354939654448,-0.69110703351867231,-9.2478214834237296,-8.2157567128106415,-1.5411781215089644,-15.997238737177961,-16.849259781892719,-6.0305493562267634,22.999054613811708,-6.7408152485135302,-9.0384397668422167,0.099436634640915164,-7.2350991477799029,-4.7427470731653401,-0.1885832871084491,-5.8288605691694748,-10.574934323213213,-16.318814034387927,4.0311657339069109,-5.0271916811344841,7.2622628683094996,9.8152606374042008,-8.7484097013764366,-14.976190345727385,-16.371989114644997,-3.8764117373098821,-4.8760975915676905,4.5546280612476773,-9.5634641360572665,7.2501777905740656,-5.2505303393529124,-3.4605888850554156,13.542817767824499,15.015443902773566,-12.173134876835272,-19.891808265572738,-2.6418104393555111,0.71955771781576428,-10.560519237079086,-9.6481679912887977,-24.099331756734163,-15.358618849253721,-14.874621841963645,-15.477220865062039,-12.838552545435878,10.344614144173271,0.19424886933297253,3.0545941804595818,1.0629180978782642,-7.4057513642710928,-15.64509702476804,-16.900150165914454,1.4084217604778564,0.16338549365715882,7.5871810026459858,-3.8828434843732951,-10.388811730484946,-11.304911157941161,-10.085852179604547,-1.0496318956010371,-11.620688432882384,-4.2261146088634245,-2.8907515338698779,-6.3253738029079196,3.8067226570656674,-3.5725801389046765,2.0465332518050992,5.443681278432936,-12.632670698214978,1.5431571741785217,5.6602726595976911,8.5744846160455559,1.6996859482001336,1.1672213661506943,3.8167222209305867,-17.124751469815184,-6.1176080486711424,-8.9982221439472223,-16.213909033957773,14.66390528175244,15.791873529182284,-2.3953674696543299,-5.5184260115744772,-11.054774044185148,-9.3316463769024853,1.42683833210666,-3.6625588486246419,-22.255038349684398,-3.093474997937776,-9.341934062606839,-11.272562591672957,-9.886549800649151,5.7958019347368399,-5.348497446759291,8.5685120475487668,-24.345836794852268,-5.6248232010383097,-9.4874039199099247,5.2502292209735257,16.495157135166959,-11.037593185145012,-9.4170050522369522,-5.2789322654580291,-3.7862463774711195,-1.6446990614368096,-16.657139180929825,-20.407997085269912,-16.146173255448556,-10.98305002637505,-8.2195886537520959,3.1289392024655833,18.222066282808189,-6.3145270317757074,-3.6368899429856905,-13.245281162364353,5.4170699520901717,1.6058476316781223,3.8548386246203403,5.3295142387177341,7.5741660557916104,1.9486986421394297,1.4953792081343245,-15.669089676921686,-5.8926693830045638,-19.45540530308536,-18.166910360572672,-10.430514410993446,-11.571820012020476,9.7058094367614789,1.1364964742744668,-8.7210880919981069,-7.959481770944401,0.15837387798714109,3.3691394309062956,-3.843794843141715,7.0131589846334643,-4.3429148285389596,-10.962978016818257,3.9953218667032329,-7.0348175734552472,-15.34080746194042,-16.829090422181537,-1.4654900345183088,2.7021231028164139,-9.6482081890636398,-12.712656371383895,-20.086933117321465,4.7823800446723546,-25.898256076765701,2.4453303395776289,7.5146896674890389,-0.3079915312934115,7.7318589534362232,2.0966105517929892,15.338892567660993,-10.00442769489517,-7.1101519742617691,3.0931454487259105,-8.9191022816638075,-8.7419910078005483,-13.96247297765775,-7.2929596474001617,-13.409318681937386,-16.085718745267862,3.1120535718728899,-8.6493236735342069,14.92373299661101,-1.6226570371425275,-6.4526490649863506,-9.1584449744456968,-2.3635789749237968,10.751595471147214,1.5637950683755446,2.2619101215202151,-17.12021679462854,-8.6723470810347543,5.0255685845391964,2.8605526430876775,12.935274493518341,2.9745954221149287,-2.3595592442826057,-6.5901837929288858,-6.0428571450747075,-14.180064418416816,4.8192923997780994,-1.7032190554390363,8.492726009683194,-4.1393117120069061,8.4143669741916174,-0.96495939046581769,-8.5630628204778869,4.0836113823923492,12.008367088055088,11.11322742567781,3.2220985409037155,7.004364736204777,2.6673318998392817,-10.606448685020801,-7.3287352797571943,4.2471938564606599,2.5094051242853044,1.9364558671020178,13.074116457969858,1.7379682834427648,-14.034974014636271,13.236662228080005,-2.8088374747310128,-9.5754914467281065,-10.514185657264816,-16.669841270233938,4.355367727088411,4.1671003029517255,5.1320851628368604,-12.655459092487479,-1.2381831983028193,10.266422490787248,-20.975621000226479,-2.4970930099874225,-12.856550229949404,-13.839303326157779,-15.894738126861984,0.34695384141392605,-5.7972472605426431,-21.166480818017607,-10.719342790530796,-4.4179121823413787,-10.465826077374407,-13.453867918342851,18.404916412977421,6.6671256036783326,4.8975162559654901,21.692617391172003,-6.6276564967867904,-2.2287869250074714,13.927344728747959,-7.100666445267847,-23.716403194465517,9.435680162017313,17.37866544999865,-0.54613211951586471,-13.135258515851252,-3.1336885144349362,-12.577839981623036,-12.165008510606052,-28.81172885616699,-14.28050509761519,-5.4129421201125298,-2.6894396880549909,-7.2778748850754074,-9.543257989806472,-3.8422440224232601,10.05064762294559,-5.0690014320114338,2.2220761788729066,-2.170309943443181,-0.68956955686382049,12.746536985834149,1.7944120225300357,-4.386389133979872,-8.731194904918949,8.2788162562381125,17.334792565550579,16.566422286738153,9.7324060846462572,8.9849364795554454,-18.229832603487633,-18.1451841005928,-6.5116373249412209,2.6346192415270324,12.155935120879237,-15.830644852299427,-0.61864356551396349,-13.503482324490427,-5.5213882432368422,-1.7374376104877252,11.185233832657834,14.07618288771511,-4.6022121830700362,-3.534492749786486,18.864745982586985,17.743244052506586,-5.2150956321729112,2.6866798623504131,9.3228846935447915,-13.88239454843808,-20.669047288477397,6.0052145571731606,9.0398278024683734,-0.53516854887593135,7.6982842436671115,-2.981273025051645,-5.9652524821445754,-7.0370217028882882,-9.8497791764198972,-0.24967256543667626,-5.6192362491637562,14.916300848360041,0.77271901942335108,-10.61101037617224,-11.36873825703808,-4.4560295376178871,-10.611505714851376,-3.8958224595335631,10.439664693885335,2.663794091918918,17.656372226979535,-14.630095275868355,12.704400551324724,-6.1143769271106443,-13.11318828443528,-16.552311664472729,-10.380563120007505,2.090383754661552,-4.3052332838263201,-2.7850884739734365,8.0944886686937316,7.0330492807625156,7.7014364895825596,8.1811151129090547,-6.1631651935376874,7.5516595887548519,9.041491852245521,-9.6067400422755629,7.9407006923089565,4.3131470908679601,3.5387603074749725,8.7353937706250857,-27.314325100635241,4.617527826781707,9.7264778376975407,0.65463859127942758,11.976723409849306,1.7973314155937909,4.3286486364441146,-11.359022460456764,8.1310837330793522,-5.5843145620757122,-4.8834815489722754,-3.48260595376694,-15.25792804822869,-8.9355065129059899,17.491572518860988,24.267006815609257,-17.693288662978397,-22.735906435976172,-9.8794939115339151,-10.314146923262369,-34.397625438509252,-16.030303781702173,-13.584462186095108,8.2959904370140336,9.1736582685581141,24.798111624050755,-0.66234017587249761,2.3187838267043461,17.545222785754408,-24.063843690800859,-4.2744936120474133,-10.903588834214135,-11.539826160174583,-16.214887536142864,34.279209623538648,26.701206187565415,-16.595267991036533,-6.0931933495681596,-5.9396062902628888,-4.7678099415790562,8.7829603726636147,0.088794038729723704,-2.4898681581875062,-22.423997952154789,-31.667298431107152,-20.315030170444064,-17.935101468455318,-5.2643322335129117,-0.503936704764819,1.3476765539877322,-2.0377510683288445,9.3551971826960152,-9.915128703034684,0.4214377166792605,2.3893892462269473,-1.6825147811807561,-15.54689518419238,16.366724588611188,-16.972622963059706,-3.405720796638831,-33.030966260289624,-3.2750356886124363,-1.1080183442522962,9.4585869512564003,-4.2673443263756834,-7.0987998429936718,-9.0440054945748649,-3.0464127163655954,-6.9853753040479152,-17.189987025835251,1.7967650422248693,13.259146028842311,14.960009771097003,5.070194420380842,-14.133588313927088,17.915398989545849,3.2887670661315731,-6.0897656194096319,3.006810065953045,-9.5758802347016037,-13.943525966076789,-9.9683111048513968,-14.166059673577656,-3.4849072436652029,-1.2137304047753792,8.4557139178902716,-10.992535923481096,-9.2989839115829103,5.6124179000146137,-6.6935024739094304,-15.168999419775133,-23.630506352893445,-13.281838712995446,-29.141014918309896,-18.278533687957776,-6.0145117497910059,-8.0893858298587205,8.6568420259339014,-3.1834375626746145,-8.1672999602812872,-7.9855259924577915,-4.4478751305639301,-9.6687721926876016,-11.845163104622854,-7.4429010123449757,-4.1955045200258088,9.7530726479989323,13.037888745545317,22.324434664800982,-4.5631246823594251,1.1068807319701317,-9.415502831057287,-0.065124970090080628,-12.577307826186559,-7.1162797606632182,-2.6088862210747275,-6.696032272423917,-2.9271833594254626,-31.814630982444285,-9.456783357351167,-15.161640534528525,-17.987226187316686,-19.022944236005966,-10.071985215608052,-5.9929888041718051,-6.5703063191494655,8.3944536187523244,-5.9451888485977884,-10.294099606607872,-0.024918171368729328,-4.454547534038845,10.480946157244219,-2.5623967095711011,6.7913447589586307,-14.682778524653806,-0.43790396169692419,-17.802487586840567,-4.3259546998847433,-9.1563848660105407,1.4655891668294099,-16.181680761248511,-17.227183605202612,-10.352670436444546,-2.6428548292505094,-8.7686389742890469,-9.0002659589770531,-9.2174140504286548,-8.1654195548595858,-8.7220981914298807,-9.0855102600604258,-12.13462855167815,-25.388810304972573,-24.634440421734134,-4.8489166157196806,-9.1790187885765846,-15.772946497342533,1.0943392993699719,-16.258300584252197,2.4959939497802552,-8.4640557091549464,-9.8184006256605194,-0.9410013430226597,3.5350855200192464,8.1684643210473542,12.934736472819067,-5.0914753270509525,-23.423416901867565,-3.1577887259371589,4.3910249121561371,-3.3160980590167926,-16.222376130507438,-10.921870208476957,-2.5934606385601704,2.7606464602157974,10.724811321217585,-8.3253611757021631,-3.1775819605658726,-7.1538693621442375,-16.802137981997458,1.5608381421109696,10.4584652010001,-23.33385765635666,0.81403651371464847,-3.8266808006679027,2.0907517796438477,4.2454458172752201,-16.241900904525608,7.3007148900367032,2.6872228511590111,2.0235559216066674,3.7464651956534483,-17.91534883281134,-5.135276225816817,-13.894239042340155,-20.88014776447465,-15.906330619952392,3.3372161354998719,-3.2014705474032836,-28.800030359855871,-22.24745594475565,4.1053930114608503,-2.5411241743444668,14.979759979376832,6.7499541359379087,15.969563380989589,-20.159120147640689,-2.1639775130681249,-6.4286250898230977,-4.0281384998618037,-3.327323623054312,10.690690993576087,10.946809825380624,30.047984024175005,-6.9476855716757022,-0.97355994147550551,0.78162508426552224,0.8344022752897039,-2.6758512334531086,1.9659393456324847,-2.9912721631391279,-1.9198333999419517,1.1466012507789058,-5.0511116685806963,-5.1939675641628167,-13.600144016757014,-14.957425252573097,-4.7805157834117198,-21.137558402833015,-12.702965542536463,-22.865728260655903,2.6550122519486874,-2.1659971097200952,9.5703281180266035,12.155980234381708,27.766991920302022,-1.6666086965343494,-1.3519640081264868,6.9886830608852293,-2.1451087898181043,-5.714560757667277,16.534833268616719,18.244927004005245,7.2453232395091671,-7.6949957453119486,0.02608368403625104,4.3681436049505882,8.2579889949634548,3.8699668253183668,-0.32488517947066375,-2.1955636939221943,6.8166830156458387,0.77111672264706876,-3.3090931606767584,0.16518130946973628,-1.4191113926962569,0.34226526051298672,20.874785152767888,-19.738851492183777,0.89710810789760187,-6.1684566361902498,-22.767938957909564,-2.8966764629764099,18.924684862565773,10.864527206590967,-17.313591456043461,6.6422489115253329,-1.2008206182621854,-4.9806222622768654,-6.0184425758406377,-13.306783825205109,-4.239575698296421,-27.365522462702593,-19.793491842556154,-3.7042926353631365,-7.5120046194315604,-4.5033442456680239,-5.4355618822015819,1.4852719596823363,-4.6546991721955067,-2.967349926879197,5.0139743515211732,9.6068684729092748,6.7955428753283398,10.150722880190589,-1.8634756245511579,-5.0972720689832283,-12.945226875043366,-11.783933700823193,10.199795709825526,-19.238587328956442,5.2019538035928354,-28.004593765151565,0.16852729603694364,8.791759043150142,19.29343607590307,20.872859057404341,-7.4343372180657914,-12.844345213492627,3.8194822751918491,4.2770729496547535,9.9024293182734322,13.738281633568757,2.1707480356717879,-7.6750362907375029,-2.5955672882625889,-1.7491005698761888,-0.44686620411864952,6.4377530852607618,2.8870745304778045,-8.1641786755130905,-23.728099931813411,2.3486248303775539,-4.5078656832061217,2.8023804307641691,4.2339429113125107,-2.8145127489128048,-0.19600866334640521,-0.15590454380268662,1.7456771288088226,-28.532488906179399,1.8050688436834303,-6.760727860204776,9.2573959125873806,-2.56353524973549,-13.386210408600743,-3.7595958014422277,-8.3933386525077101,6.8079382840736979,-12.963594300450609,-16.097117864906856,-11.364067826364966,3.7255961085627232,6.331129733180469,-4.7114133395682165,4.8836185758855928,19.234671006072027,-11.429026970335212,-9.2636807759033086,-16.959775359766532,-16.588766118623024,2.4604932906175874,-13.105963552927271,-3.9216685336898949,12.174496793230894,-9.3581852356447559,-16.05274840601119,-21.487095890050561],"z":[-15.443307527216424,0.0097927242631325444,-0.87444044964232115,-10.348383436973625,7.3160227681541956,16.837096264799428,-4.8737085884351057,2.3010604417187932,-9.0797526658556471,-3.3943493783442662,6.2514163427807157,-21.941366005267646,-9.7079454962286746,-11.10125099619018,-11.844644335908173,-10.354948117924408,-3.3838756603546987,-7.2495661281405868,-0.88777912525551883,-9.3512787822053109,-4.2498637749948269,-7.8091809342947887,-9.2576276450834065,-9.2824351791289335,-8.7316451983439549,1.3771517859677243,-12.130226611558832,2.9639939397539212,-11.731267806308667,7.5152259172512554,12.84916511199039,-8.1842260260769333,-3.9127681827651086,14.459458588012122,8.4091538229638552,-6.1494399995320856,3.6902872783905498,3.5613134624846663,-10.940010027872702,-8.6639632869234315,-4.8350553458288017,6.2808924895003715,0.70437104170120457,-6.5492890918731277,3.415863757413883,8.3327335051167086,15.131551691644708,-6.7669630707524924,-1.8551846865289061,-15.475575100283013,-17.086982558359434,-14.098203494117964,-14.701396200678014,-2.1218318400246532,4.2545265431260582,-6.3905396916835588,-13.54076290790583,-13.701484146950241,-9.94140946635512,-9.5541144478598703,-10.754119909464292,-4.880335575195649,10.722006919177529,13.796025273562119,1.101048710732073,11.068901544403511,-15.617993808335239,-13.003550170448367,-6.1300904451409908,-6.1545301381122171,-6.1295210513101939,-7.4571305093869178,-0.38567048648953084,-10.193527236245899,1.8370233620716281,-19.005162705547136,3.2937204682688788,-7.6639156875172869,-2.9265476616418291,3.5238039082957182,4.1003048714856947,-7.7815201626307902,-12.282522252002011,-11.267916990564034,-13.640811919617054,-7.436167583181839,-3.8650982716244435,-10.19348039228889,-16.816359384668583,0.79400074961684619,-10.774051804690163,-12.427037257096956,-12.109265287646688,-3.0300537800101255,-8.0345310406820225,-0.3556518084548605,-9.2846436231431877,-5.1690400467379582,0.90196516462710885,-1.8608783588983622,-9.4277238816138951,-0.98575167694488575,-1.1768583616474213,16.671965975306371,-13.804162182384486,-16.392028728359001,-2.0649170356918662,7.571832199206713,1.2630830147769687,6.9262761302542337,-10.20239474167596,10.798091638638613,-0.76972368694827908,-10.66896681071416,-5.1385150910602526,-6.0468642074198655,-9.8992688283956518,6.3937663841615269,-2.3029267533706852,-12.39735391623689,-4.1397662744004791,-22.631458324218585,-0.20466476334353403,-0.83555892557580214,11.177814172029002,-8.684993449502814,-4.1873401905164478,-2.5837036955533139,-5.558852090720281,1.0856548711190892,11.634762819756439,-13.092364371150332,-8.4889895609894506,-11.628875095783277,-6.4381108617141587,-13.647229910584947,-13.091245643229851,-8.4173695895842862,-13.478196105637881,-10.25921247711786,-14.385795266321958,-10.212238216264669,-4.3193379164988652,5.331358217766617,16.865646964816335,-12.874579227334522,-12.738158464989167,-2.9406474061123609,-7.8261034849325526,-4.5245748145274263,-0.22022566270720478,-17.345921634120685,-5.1108835700284185,5.5911261902220382,-4.4529001080155233,-5.2087434155529921,0.97556874746509503,0.92845037226307048,-12.715721818465211,7.9212540929063699,1.747375474750714,-16.549016868004873,-3.699908241222122,14.555425507396082,-7.8896366156708435,-0.89750531589533888,8.6852532883342857,-4.9114656832057983,-25.186987186163712,-13.009457999473785,-2.1540239838832447,6.3196982911623616,0.59454024225392732,-25.976893183230331,-12.913132820438785,-11.602611662996887,-10.966601013709584,0.56234402680383444,5.1957063595383177,-3.0212855257834468,-20.305073566126698,-12.839946754037765,-8.7471065122134473,-1.1232444104058235,0.049906889001088153,-5.6841147428586725,-7.1202520620160694,-5.6042279798801049,-13.220537826095335,-14.041066201330477,-10.941248871280887,-4.5160759877323651,2.2725999022968759,-6.8236988555761631,7.7893440989827791,-10.266366960029117,-8.4564424138403318,0.64766143807878518,7.788701939568945,0.54600473041549458,1.7176275280612379,2.1755669200755223,3.924150645335911,-17.419268872056637,-1.9151591738359142,-9.6765219426719806,-18.176822891416201,-4.4518819593847114,-0.84850556879075933,-3.0360378434738382,-16.093517898458924,-13.299927979257737,-1.3674706680273701,-7.1158514633371732,-12.728854875862625,-11.65392277241213,-1.9824102743108669,7.9967517451897239,-5.2316267054576908,5.1115465231760258,-13.102144472386168,-2.5972937339176485,-3.2040994544304318,2.0452444306194084,-8.6752071216586284,-9.0597422936470924,-5.793381063035854,-11.148388189410278,-1.9265639889249697,-17.673978353095709,12.37652108326763,-23.105956424899695,-14.969733626585908,-14.936362227282279,-3.1655124217372568,-4.7708001225869161,-14.357348902285272,1.8474457326258173,-4.789189245275387,-3.6838545793129929,3.7072683140938398,10.453243361425393,-16.311522797726731,-6.6861746111355185,7.2762524263567538,-17.420609987294046,-10.838110427195886,-4.9528678679026124,10.656284822067596,-8.8565834619808328,-0.4562434762483949,15.985306052041214,-7.0699754662136201,-12.507676858174325,-21.930050821714904,-27.539864433740899,-8.1250547554747374,2.5940436298037812,-11.505280379813097,-20.379419926330421,6.7276630926499834,-15.154740835318172,-3.7060162907894973,3.1120993108988242,0.91455064932018493,7.4513496473280973,0.96488090104680768,-4.2612283118987646,2.4885547220906723,-1.3920555914299002,3.0049081604324859,-6.409194689727296,-2.2920985105573046,-11.937072151424777,-5.5941491482004135,1.5931438137222873,11.947428596728766,-10.10008269565658,-6.3042129448399828,-4.7630847622241266,-2.1589827467560876,-12.19043606219328,7.3369578318773438,1.5997187332210656,0.5715590646327261,-11.445980202076749,-7.6475259980272554,-7.4230595178094596,4.4570155473864324,5.8310553396069151,-7.7374370497780145,9.9669368413724513,-6.3854946906406678,7.0240185441946963,-1.7956160399674965,2.6696259816323162,-8.8707861967752493,-9.5673298130150215,-1.8803209159053522,-6.5811206501893675,6.5992422474329544,7.9309973976570349,-11.184096493511316,-4.4873419441720586,-0.15864273800876577,-8.5884590593153991,9.3424202829815499,4.31255859934332,4.1366973064891717,6.692517873315075,-9.1613599413287883,-19.835337406660116,-16.642272271287037,8.1815698728990487,3.290079806319731,16.025852146481181,-9.098783743764459,-6.4304009720844828,-6.3377502055466088,5.8583248040664433,0.64794504479604231,-9.2937951517647122,-2.3521422580691294,-17.584199375127074,-2.4842288492305102,1.6387202553401798,3.2300915413972215,-1.1471958771352941,-8.717546154526838,1.2669998691376723,-5.5238331076776781,6.850619776490066,-0.85351755527406625,7.5711974033637679,-6.0406708332481225,-5.5442653394309049,-8.8544121828806741,-0.16431024304549091,1.6168768125885649,1.9215948743735796,0.81479374196844578,-3.8312619444231575,12.274313938668744,1.677629200049485,10.559477180299202,13.629914793592352,-11.740467397923229,-9.2642130144478809,2.9562501856376846,-8.4952203738312821,-3.3368025964410846,-11.234162344578499,0.66637196195161774,-7.2769623589205947,-4.1148902197968704,2.5638423589547346,12.321815903480985,3.3613851283251179,-9.4779162596656157,-2.0132966312299136,-10.965662017354262,0.55366175390304118,-9.3656529731276201,-11.41609183834797,-2.9348634576024404,0.93201143610178261,1.6139186180261933,-8.0107595206265341,3.8506178006245775,-1.7089210933986612,9.3734834160525029,-10.201224619799463,4.5165042940966611,9.2719723646439416,-14.41007748349778,-8.2094867786765153,-8.4854262647472254,2.2483114116809984,-1.6657824120983362,0.74650370042146141,-6.6345171800387863,13.468616032114149,-2.0691981452162671,-4.1130956787383175,-5.2433011802930878,-9.7873663889259586,4.3778025297576342,-20.432523819478696,11.497632379558581,-19.391147168814829,-9.9213650042467982,-3.4840309132790788,-1.7989760453609276,-0.02191253067623224,2.9397168419067978,-12.873896031341522,-6.2005793674720042,3.387223027021371,-24.573752423609196,-10.718936138921691,13.525081789920835,13.963269614544041,-9.3264637036266365,-14.315651660987527,-8.7422994005046686,-20.298928724551537,-7.5374498922440889,4.7595199378330637,-4.0240456648412222,-7.2904354860307654,-2.3298071027631551,-8.5241461764939697,-7.5965818286484419,6.1379882970826651,10.995007055464855,2.8143097957651015,9.5919585557875173,11.563910477633362,-14.080864522315817,-10.550951207163886,19.789243847254884,-2.2523461584699564,4.2369617341375312,1.9735229431835324,-13.922945235511634,6.2448585295964234,3.773854691338475,6.4714157234674854,14.145551717762192,15.864841006933794,2.4358712018551483,-13.722567412417698,3.6973002378421254,-4.9048963725887829,3.9280051681679371,16.469513699576769,7.1362286491221338,-8.7126193187700345,3.4408382143476777,3.6895005845335551,-11.099503514412534,-15.364454070211757,-0.59852957967229414,-8.9296412517616659,4.3102513454575133,6.5797407727078188,-11.370425065111251,-10.198990061277708,-7.2581427411531365,1.8565351707717246,-17.889073998846612,-7.6628242898461405,-2.4456810560086919,-15.298143307270605,-3.7272943654601005,-6.92959618562712,7.7177569014564957,-6.8693175071670804,-9.3295275220930023,-14.243178649508666,1.9894133112607044,9.8046781768732618,17.683264701613563,-5.9502992018292735,2.9024186604298099,1.7729423353276739,8.2113770304206994,1.2828816272013046,-12.061882225447178,-7.738340106402819,-6.2954622950161641,1.3534278267482893,-3.3451275214680414,-6.4902546422346541,2.733960769661794,-1.7967956286810804,-6.0871338191498747,-10.188565802477758,-1.0649723465255563,-10.068400931165524,-27.694934029383507,-8.031309646889877,8.2711846064640362,14.573615758465122,0.10500423716416095,9.1026801559966763,-7.1043705987762698,8.4899076275912027,-10.966515363430224,-1.8546052771379697,-5.4091419796166118,-1.4754085931626351,-11.045978620071971,4.0988152303296914,-12.678731486718593,-15.391587942075784,-1.5402202499256366,-5.7722761125597533,8.9055641729528521,0.044225361296012421,1.1942793706192312,5.3995411546449388,14.58846003439149,-18.371584403234881,-3.5065644154995916,7.0062009736397961,16.10866519833473,-20.124689770424212,6.4196674675379271,13.950178722879247,5.9166352279171797,-28.580110995718371,-23.561361677290211,-5.6432036474915481,4.1772030035683407,0.81267478259179637,-7.3952590771260267,-6.7777678374409698,2.9667986804565794,-1.7906772709783201,-0.54084038784598465,-8.2352413845402506,-0.99794479514652745,-12.093910221509327,-7.7286173148815642,-1.7609456323380277,-10.00343531898371,-11.01768741261642,-5.5583745671504339,15.995220618640433,-6.2551281853895357,-19.778252264379073,-3.1983078855151477,-12.315021490453971,0.56106433746666662,-19.660973386460707,-12.23672378706959,-16.031635781644855,-4.1833501483419226,2.9395305246037222,-6.1937041757323019,-5.2768340587958997,6.6698199384433705,12.882409924424893,-19.266573803046189,-13.806460853660617,6.3800610634813975,17.039632572718215,-12.242909370717166,-9.7899768281404445,-0.10258815386342142,-15.352905796407036,-2.7176797787455844,-11.758660293332452,-7.7164674538119931,-7.7048421146441814,-1.0427521421207926,11.665069031886516,-3.0509240547841761,9.7194073171304431,15.724229757766848,-18.446026835875536,-1.6364657600062766,-11.992544292451418,-3.1944111730502271,-11.064613023303204,-13.156104588320956,-16.777554204853548,-11.623497370787364,7.253301024350792,-1.6305163111472738,4.1189324812701775,-8.1286701826618142,-11.628649903762955,-4.1348248517198609,-10.239526077793526,-15.063507330030882,-10.177224372504803,-7.041651626592289,-2.3029284607988654,-22.522219568489163,2.9995830595265343,11.434539966217764,-4.2570720487803246,-22.839782944651837,4.2650405052993943,-4.2842946181557995,3.9854398974794134,-13.205681643724125,5.1843797528410738,-4.8604458024299539,-12.323088702898412,-16.007504534158855,-13.339689011146614,-16.157438000533865,-7.2378955166940493,-2.1478606767024702,-7.7633201906483951,-10.484284076453523,-6.6292125792142427,0.51065760652541448,-8.6993754872486093,-2.1190094955281236,-4.3280643468974764,-2.0838940694811621,0.074459964409456875,13.450781785088262,-14.975297378793877,-11.397517546382828,-7.4904923803552395,-26.623488449745206,-29.082845534925969,3.3839338057898805,-12.355619152119761,-15.568179154442413,-13.835455069061741,-1.4618533169975263,3.1629103396221296,12.035433739658664,-8.5117985198874901,-4.4592820069880235,-0.42277682823711121,14.30024901913483,2.4739053276504848,-4.7368509701486818,-12.284150445950202,0.16426090677602961,-12.917726927264775,-9.8878296733613649,10.602514267735792,12.931599478257557,2.7984173849534995,-12.802354058063322,-4.255996623677853,-2.679967995708902,0.35634873355127078,-6.6866457008397369,-0.26397232591257352,-10.884493107541514,-19.825303316923364,-9.1681833333152341,-8.3344744504122872,10.302280999510202,21.880279446722803,-11.740371252059793,-2.6178130092983851,-19.659812381638993,-3.2875763051466325,-5.0129138695073312,0.98298862131218823,-13.263003892180494,11.329161010067285,-10.571993946054986,6.824625027032388,-14.396156655733403,-10.897160749250929,-1.8515187901165688,-13.620446170080912,4.6033671576928175,-5.354965213721905,-3.119915667110845,7.7733999062047943,3.737995857054758,-15.937364209048985,3.7958409122010726,14.896955380987352,-14.479658524130206,-13.818756023037379,-2.2430236272176298,-3.1191848773307309,-12.646462069369644,-9.1192425846594034,-11.353538701417722,-9.9665879241599988,-3.7229947540462174,-0.063011203750294856,-11.95860663750099,-4.4566116825571793,0.078865595228002211,-7.922965801918755,3.4127848660745368,-11.383380304994022,-3.8521089988414294,-13.263110576733592,-10.145092721041765,-6.833678192070872,-4.3722675347329396,-1.480153089755978,-6.1419718994500032,1.0834867967392849,-7.8197794065590225,-0.21764072180899444,-1.3591265310516749,-19.46578748317102,-1.7351988240452956],"mode":"markers","marker":{"color":"rgba(255,99,71,1)","size":2.5,"line":{"color":"rgba(255,99,71,1)"}},"type":"scatter3d","name":"aa","textfont":{"color":"rgba(255,99,71,1)"},"error_y":{"color":"rgba(255,99,71,1)"},"error_x":{"color":"rgba(255,99,71,1)"},"line":{"color":"rgba(255,99,71,1)"},"frame":null},{"x":[-28.983241281140938,-33.242102564361829,-6.277162728288876,-12.312850981217297,-22.380293805487806,-3.1397729807943726,-31.95300235349309,-23.06148581434876,6.8088452457998816,-16.705644169891045,-30.237892512838012,-29.319120613840891,-29.239743633379852,5.2080885633527494,-13.232410836108556,23.559857715304737,20.730068499247615,2.7029015096361095,-8.623134668702793,-15.606346650326904,-17.254812302717934,3.3093690042854575,10.009760181119425,-17.121451741356896,4.4785511813247432,1.557772620195268,19.466239920819532,6.6228950925361669,-10.845805052791341,-4.0372024607179444,-10.01021743367407,-13.699028702383258,-0.70145118203582912,-21.825033051359561,-16.423250522472191,11.125115398526649,5.6755144851589296,39.781443760468221,13.795411170277989,4.6439929219265581,-14.111554326527834,-2.3902860149274074,-16.093518392286025,-23.119487644273129,-28.39715870852897,4.3658840498116325,2.612135307834182,-0.98140204268658315,20.964023859954629,7.5857220419851803,-6.648015237401169,14.76427811384711,-2.396417894800885,11.367142977396647,-22.716276325388893,-25.744805772284991,-20.138978413391051,-6.7058843504641619,-31.892076442616155,4.3554777462826504,-13.232458272694332,1.0985200929146832,-2.7234554195825846,-19.009397582774479,-3.6965693449464756,1.3201897906807887,-25.505685189888769,-50.666478337671172,-4.4391742284817948,-33.423491297659254,-6.8828441355296137,-27.603378683012277,-35.510579891140807,-20.463601496686938,-34.428354583440559,-44.281038794110167,36.831860853978355,15.200971875288491,6.2841959800572207,-13.02364413246174,-24.1198039275627,-36.439208887574466,8.2884126875446853,-0.0030602906972291575,30.824925566963547,9.9532038128983142,23.399295187773522,35.672961046453821,21.983025175407775,5.4639335515495926,-15.009210989143364,-25.354755617963363,-29.799942550419654,-32.685216583739802,10.618051170090814,9.6267438914511114,-13.400798958778591,-28.18403338106469,5.0634562776155487,6.1958314871571565,3.1153232960100383,-15.081026106602595,-16.796486254090944,18.285864112647321,31.410417484270397,8.2037712829843557,-4.3252968764863633,2.6221868146663674,0.44297229318332937,-8.1660138710631607,-26.263350443280277,44.240538778598918,8.6405791270519696,10.300120110642256,8.2425802954727292,-2.7181116375789061,12.758852409640721,-8.182125438087505,2.1416427086858034,-27.995123747711595,-26.714830379219503,-2.6635222023557112,-22.370005026233763,3.1036858703876917,3.8806820311598491,-10.703009145846524,-8.7474249750825486,9.0684140308835683,2.655658432888039,12.062969916117188,3.1173058139387013,-12.542846252193804,-11.463896213172815,-35.022850117819871,-26.874698476474563,-10.791298750971373,8.1547034485523682,7.5027767126145655,5.2891761323593496,16.957282886991727,1.4602969151013758,2.51353304788464,25.342210308681683,37.040845892023341,12.530462561961459,-25.247338514701507,-30.705355019777301,-14.733775321104641,-19.90401319332727,-22.296866762706472,-14.935776608058523,-19.934696140412669,21.171657737070799,-29.843439364767001,-14.387219719695612,31.65168328669132,12.541640203992502,27.649889454853103,-1.7738400416438909,-1.3251569719565535,10.511019144046777,19.675271275383917,-0.099079953939341478,-17.344853719523673,-24.096217125203772,-14.705356956237015,-30.961035089254754,-17.113327974782241,5.6290046764543398,-5.0895469206085613,-37.317716906992864,18.30340730718574,13.734031165680197,-7.225581155942173,-14.258195634684634,16.873479926507031,-0.83499537682259328,-24.484628517848414,-20.631091456966871,-34.756349809241343,17.830143155329175,13.645833198291422,-12.305324115736578,-14.048577666795731,-19.318078441799951,6.6294787569622189,-2.6297107349863387,-26.148544316925317,-24.282912068534095,-35.747321621086222,-33.365250117716926,-3.8625691619872158,-29.821598282294858,-0.36068534045957856,-20.654728040926965,-28.499495263180254,-19.636014182298215,-1.3663687457630977,-10.626681354980462,-23.372935537190813,10.827537517396294,11.387942459602451,9.5712204672404475,-17.818792360143561,-37.46865439924489,-29.841789234855163,-39.95040394284996,-22.10872343575852,-43.98230817824281,4.9764564021874538,1.2719795688229103,9.6580216688999698,7.3375706869862611,-23.289312814674826,24.201709565030839,-0.67761960885309791,-5.350362204501244,15.072354623031456,-10.761099531564723,-20.564254031211675,-35.553542108685349,21.286903083871245,13.376491999867257,3.9152909339427207,16.075062238672963,-6.8386507780846175,41.258247618552033,-0.75492851801449434,24.539656963294554,-3.4458374359204997,26.110808421931491,16.907083851000312,16.855636078255603,-23.976201293730831,-26.672382297903141,-12.476236093576103,-0.22305549083745377,-6.6220712346648822,-44.631723622545159,-19.001657511107414,-3.477077265575923,-2.6604897408224693,-12.75306572725129,-4.3020791313566695,2.1939064244014168,8.1338784729696876,11.952514028026256,1.0815153680796035,13.83247863721563,2.6350169565328283,8.3929742404714727,-8.4620668532214403,-5.0925405037395972,-33.711266830918099,-19.875009726675145,-34.969209297599967,-17.94086641434124,-29.223480246361451,-33.054951233442338,-16.798722295603536,17.114250606162788,13.955812642338413,10.518621004404778,19.631447762119933,-9.1736593025376862,1.6293528081576765,38.869112940473514,21.891037154464929,0.70781332447540046,19.281094638258097,5.2584535729421633,5.2795546786082141,2.5977923456491165,-11.907560876031933,-35.050765789238412,17.772927341215397,-18.956763082338338,-28.348259249829503,9.4867572354361585,12.319084831616784,5.0459241705812241,-3.7909910860702225,-9.0408644653989469,17.098512942125915,1.4959493877135608,-2.2702190130079205,-16.426039776580303,-16.813428593415278,-26.746818157708184,5.8573358163355502,-25.166803267869504,-14.503667396534087,-19.639847144357024,-14.645001286598255,-34.732382580827405,-10.914105881001666,-7.1029198361935295,-3.8402382034161486,-1.367558732604764,-22.18279318416846,-46.574663195555658,-38.914554850577147,14.385524260492371,22.243782170408437,20.724918745851586,-14.005438739482626,-32.693696379308726,-0.023309418666515598,-25.970603842275604,-20.110000395070244,-20.100044284026573,-2.2074344182865668,-17.529820071700058,-16.137348685111274,-25.879689177428851,-3.7468608964581289,-16.231117345236566,-30.62853242560476,-36.309849524021779,-6.705719154376836,-17.304709881278871,-37.714138905758837,-36.439861163143625,2.2276429391825094,8.0101450512264911,36.967194792672267,-3.1396514609154909,-2.993185454153874,-14.960931920987537,-33.809673732789214,-20.458064395751201,5.8975925785073153,-0.66550184071021834,-14.683105042578928,-38.371386357045907,-1.8606925756829726,-26.879335134441742,-38.621863376048623,20.925513177006049,4.4490312309192239,-7.7584008016000325,-18.909170465199747,-35.879374653984101,2.6277782143770607,-5.2406258061146147,-6.8659615071921873,-8.1804401998577969,-27.744856388213247,-13.592594933205717,-6.9239879126739661,22.58942911813228,29.249614218305112,24.260526002754958,20.066588190030654,-51.366179873263377,-32.150337690322907,-12.600182106875717,-7.0867154046852443,-1.3279397514601592,-17.503801334618615,-20.56609238779415,-25.528381247556638,-38.801364281047206,-3.9136738457174389,-7.3039435225391109,-0.8240971612915462,-34.629695372641457,2.4144478774363929,-2.3629768293666276,-1.7574019372118121,-8.5405236990373616,-32.635870220515386,34.17926699482048,30.701420317397027,21.668718733287232,-19.406986944406238,-12.513100409914351,-30.853503007541608,-0.61981128708034106,-7.0369794001001678,-0.44483874800956946,4.6334063566746488,-32.109908567133679,0.20473904605548685,-1.7885165124562679,1.465670040391376,7.9987568988760591,4.4136465100133222,9.745237691802469,17.656574199305762,19.907251159905311,-20.515654779625383,28.11438660027563,3.6014891677079577,5.3480737316765463,1.645408285046591,6.958829246864676,-8.692118448752824,-10.355313887908238,-0.74717063093349079,-26.53644094329719,3.9083017705122831,-5.5192845589640562,29.019154667454831,33.502041450923436,43.943956521346188,-20.625720852386888,-31.073106936980732,-11.173202227737928,-13.291945967662389,-13.993165930612992,-23.693990056081578,-24.804294959724245,-16.117621599343316,8.3600763325933674,7.9711472019393739,-19.586748131681333,6.9338156724760074,5.8824412579906653,-10.415002142212355,-21.192112412800547,-13.509039117646031,-1.5675657699092413,-34.893592070581995,-24.742804319798477,9.4972448311083699,1.1401166926128854,-0.0069906256669206076,3.5715921006603688,13.023825772148669,-4.8850584661929348,-2.2541394983637364,3.2840662951947897,9.2764312788506444,-9.4717091384902492,-4.8420099350268169,-21.115254628084546,-1.7011840822501469,-26.027852790601404,-28.72086685449948,13.572892104822424,-7.0370285617958492,-16.808699884727989,-5.3822232111210209,-41.039619297384952,-9.7511946591208059,-8.0676927451003309,-18.494619896061369,3.7290621248541753,3.3731620422838486,19.75109999543351,-2.79000253424935,-22.03898120624968,-2.8170446857234959,-2.4187594313353089,30.420131735294181,-9.5000281384910235,6.8404694649774145,17.692944775595674,-19.180696086122218,-36.462278778209267,-38.193620452171231,-14.174816715783926,3.0418916638508304,-20.823953719821287,-8.2609422665792032,-23.385781885730765,-33.150502955534101,-5.54190879367163,-12.344113889443484,-1.6381073353786522,1.6024079038688244,5.7881850033398168,14.051564550862148,22.061867100185619,24.496364169037129,-8.4909502884037273,-2.9900365667480768,-16.142884101174772,-0.10100997626467192,-22.711315032081089,-8.5340893592194913,-1.5820024992921118,8.0313585255666844,-15.618149844983332,-24.153193770084943,-8.6081977592491707,-1.1576689724026885,16.917342694578196,-22.567131977147913,-33.101919514090717,-17.828022856654339,-3.1607882986018456,-23.263159707557584,-28.90675864340075,-9.4366386377412983,-29.921091830020558,7.040360146085594,-0.92898844352299181,-24.597917940446802,-12.268093204898319,-25.714165982939942,-22.440823017581859,30.208714045330151,14.619553774193884,17.712467167198028,7.0087162623566952,-6.0305350371745927,0.30607831054294449,14.065019686956639,15.430364267968343,10.794382207767013,-5.2560675452265944,10.228299438560155,0.24229625467686944,-7.3304133264853446,-33.3114709176498,-36.326038486597412,-18.614356437202169,-22.539323047494332,14.727287045288501,11.309081513287312,-10.321355525594003,-8.4234378704397237,-4.6835816321596333,-13.258454249283202,-19.795111605793487,-10.481141306059055,15.293403326564922,0.36028461199632655,-11.019163058108189,-25.227011683198434,1.9181700496655305,1.6390274719776912,12.620860342638723,-0.63481762602379788,-27.16292103037199,-24.946434777028216,3.9037019242459436,-2.5149307555324718,-14.766295819426729,-24.731973847257009,-22.570266909873229,-2.4329330726826059,-21.000838195589711,-21.988541141773535,-27.12708768530732,-36.401416763600359,25.117149962623802,22.436564939553477,-10.227080863799676,31.228858577174041,15.276277898688296,29.619943889200172,-9.6954361409061409,-19.308963728228075,0.69303650529103267,7.0602568289748548,15.857105116597523,1.6398287204365225,4.3848210846134581,-7.7205743024705145,-17.118483514144362,4.4332464224553831,-17.25209961017028,-20.713502202271382,-20.228088176936097,-28.624542458671403,-20.77733766067179,15.834806753130795,-11.680627787186429,-22.401338057423679,-11.955538525027833,-7.5846878643763995,-51.026823632388606,-20.831208530854418,0.040899821634755786,26.934986762827315,10.259048469854012,14.774493512818637,-28.205738957328681,-11.458030252088079,17.203497326786085,-1.0469468648865494,-10.367211330919048,21.099376814243968,-3.7042094738880489,-2.0790068746812924,5.174259167507997,-5.3654750235285622,-18.732787151831236,18.36223862537539,7.0057575997011066,-5.7281691471745484,13.511245360798359,-14.835035274782911,7.0124153580410828,-7.3860823077685751,-29.883053364994598,40.408422822237924,24.146859834698439,-16.247664661192761,-43.563097324047959,-0.86703821760901212,-13.757792402448889,-25.055569210070651,-32.020112696348029,12.040046150510701,2.0333246141166921,3.0601569640671706,-17.532126779676638,-22.290822143624194,-5.9691087457371959,-13.621140398760838,-6.6197578511966126,-2.1415154785451453,16.423448615448041,-8.2957777926935563,-8.176615098842726,-26.306508072652143,-26.740051093838346,-15.852844242568111,-28.249826918920974,1.9721100303586585,10.629666899449722,7.071307625911059,1.6214161748436073,-7.3829632044974733,11.124082791692542,-10.239712745064793,-6.4801568074745415,14.029824826488854,1.5115221649853128,0.28694665949249387,28.843432940246863,39.920215549686134,22.188680243009422,-20.962577666163053,-30.095964794286786,-31.43697535412019,9.3630075858278765,12.365037499382765,4.4116561718652285,10.761436729440343,-22.16924090144866,-17.469277815090724,20.183513444384971,4.7303896899538218,0.93712524171056655,5.4510042070655409,20.06353162162458,-13.586374488746522,-24.270751687376215,11.764893566051517,-13.338883802172459,0.612343178353558,4.0280589996097174,-5.1354662178270605,11.151444049690056,4.2378205884981757,9.1696010532969385,-2.0210332652245402,-22.125659000134121,-4.4400829081200053,-15.381430289019196,-47.082388005989941,4.9017045634486474,9.3830512422951777,-15.907375469057394,-6.024713489420126,4.3479569671712213,-18.221593606517537,-41.160656975071568,-37.757688246660706,-7.3953724313569609,-4.4094180189879948,-1.3186134305692172,7.1981215079841512,23.379747738506783,1.3826499605779776,-2.1784993264913091,-21.987242788736708,-46.480890063806449,-45.47401204256343,21.977544052843758,-8.4467119090017508,0.3884745475427353,-16.481960741314825,-15.130489357320604,1.2133335940529497,-23.142921911982697,-20.365873082572087,-4.3481306580839565,-2.5323530973804744,-35.839296451539354,-33.804597837023366,-37.499245020120227,-22.273689699718961,3.76482243581661,-10.838479327952339,9.1390381979454816,0.44020581247471602,-23.60402314840379,-33.4664001510139,-43.850739395825791,-24.288512528427908,-2.7078395182632686,-12.670737344057285,-2.2161346495894678,27.936499047552722,7.769791425983219,-16.09038504402389,17.782632549544552,1.7237286476155884,14.735982610976274,-5.2982824555521795,-16.746063651141444,-18.706652646555792,-17.654189499881209,-22.76211315427258,-33.330890029822093,-22.106915538466684,-7.6701807961613397,10.344522413859458,2.3131200458464005,23.547042124601088,-11.041465893321833,-16.947873639019306,-25.516585809033526,-0.29955955460946365,2.4206147666760285,-16.479696783698163,-22.307044431085377,-22.749473366760064,-21.342732269448881,-12.708855255959362,15.839892997585205,-10.526068428720814,-19.655632217216539,-13.62008354813422,-30.45288000046353,-43.956639303714887,-40.179787866418366,13.203954406661639,-3.3187058975914203,9.6500975432401077,1.2776646284140654,15.461295561560034,-11.921443164414713,-17.791815692418073,4.327716833046475,-6.3904499110048816,-14.093714718332215,-13.093956118593077,-5.6893488930958886,17.534368627722003,-20.774400595864218,-7.3273590063109948,-18.556922612432277,-43.60112572973572,36.944063573521873,37.751047779415636,-19.911146224746389,0.37289755254849644,-27.984734830195588,-28.227397751308384,-19.006611460114986,-12.002532976589615,-48.53748872560292,-6.9766824926413564,-11.171920471420997,14.752791774029211,14.961399660154903,18.9919078614572,16.304430191786093,0.97577259800442118,-28.57019204249432,12.107256764457031,12.984180441784712,12.09410486837457,26.285661187360127,14.042244071610918,4.4109237429023507,-4.9324244745670587,-6.4563426215105455,-17.707239422898414,-15.741449852048317,-26.418128022346469,8.4072973201781629,22.559180378365429,-9.3369889691503225,3.9967321047868194,26.870900428554037,9.7671997529518872,4.3175637371493565,-9.9221084080863626,-23.805943969761028,-0.43193365446865739,8.7979217641808649,27.260228547465061,28.105179244084127,6.094859654191163,30.791212491383096,10.652266010199465,-0.91164381502022307,-19.679000343685594,3.5077237902701732,-17.427709740545627,-18.229448320217134,9.7855643379133195,-14.587269156602577,-31.2530550364866,41.546007121185035,37.767976307511319,20.037599965624995,-14.850573703050115,-8.6944671632703852,-30.008346130330615,-5.5564159414193206,-19.722118818051431,33.984192354651022,21.557773531403413,10.24540364251013,0.72401088234780142,-21.974356178683827,-35.36640416092775,-10.417200964405939,-7.2825002344853855,-35.253178792104997,-10.77599365724152,-1.4411547117695676,-26.985527775478193,2.5744780265123124,-8.4894001458658312,10.342820302930759,12.925761442853448,24.802873643116264,5.7086928851573786,-17.443420098692599,-24.1982403024844,-5.097350683821583,-14.862431803877001,-42.33710337711782,-43.811190871289995,-48.157694480697337,5.1157334483608334,17.348826829239748,-11.124560818262397,-20.139787437181447,-18.293575481486158,8.7466953607365916,24.555707653828577,-6.4487180901270058,2.0026760560519761,-2.1047162650718882,-28.144741718294515,-4.4331155050197149,-41.044606611542676,-31.013470283515264,-12.176964774546777,3.1315483167071454,-11.898432049064283,-35.495676209323278,-1.2837974922678514,12.566351078031014,-8.309236313695866,-8.7028970233247627,-23.31015074739927,-9.9629067880531998,-35.560822860981233,-3.8489383472723562,-16.799395614899847,-26.364591840764938,-17.596221831900589,-63.700825015319637,22.949317149261343,-1.9902983344867198,-6.4708725623330574,-14.019147743433841,-14.96274189132768,-34.848422533627087,-17.35324981506675,-18.499383434438137,-4.9477979114119961,1.4033874992633462,1.4136802967464854,-13.643553191816601,-12.044911287592692,7.5132081556861676,4.4560179467179051,16.461919299650713,-21.379125728390644,-23.273768794743184,2.8185097506216139,21.660186120961132,23.775363144555861,9.6701773017541157,13.889840216663485,-23.933399540172395,15.425951643650787,17.633522050600572,-3.291554453328053,1.1649805321005211,4.1073861559696985,13.14993042102814,38.669461532104037,-41.937606235215583,-19.111254680338909,-26.456274330062886,-9.3354200650809851,-18.52987296442258,-18.167576860673915,-35.366263237949596,-41.611602716492328,-46.144376580933006,-30.07933760725038,-32.107134602485857,17.725055559859236,4.1317899989321294,17.694285136156157,-17.713875656065497,-21.391291121590839,-5.4709163461282779,-21.717703514230941,-26.559539172901609,-7.6948230159298125,-13.061636892646728,-7.347967745812876,-40.96710690659463,-19.754690419805559,-16.290722491222123,29.448567827743823,-2.0506674338911735,4.5017361456481861,-6.2960897192363072,-7.2640151187499331,-0.73126130148087942,-33.512764190881612,-32.166531756949439,11.519199686674446,-2.0267236798692019,10.030357612694957,17.223237190047758,-0.13998069579113659,-2.9048705904713539,-9.6284233122694136,-14.894513233661867,-17.813527867141993,-28.161360847270554,-30.647220665459059,-10.793748648521646,-13.790238327206573,-3.2246526522338232,-16.839699423929989,-23.629373545188979,-16.809658153723205,-29.144079583292672,-2.2323749735449714,-3.5644758942892238,-17.982063768335447,-16.964594053185085,18.236794234543986,-2.8432646500177592,-14.312375476345192,-38.252359179243612,-43.103283279502065,-8.4386581298253152,0.89363190512847002,5.5709643148334562,-17.830428024786382,5.1624970879493262,-69.361979376204971,-9.8579656832259577,-0.33947162299696909,-20.968255231616375,-16.738749694598631,-27.879809166464156,-25.407326009723384,-19.85213108270441,-12.579389269279574,-10.992498647883179,-11.188292536098558,-14.58606541710442,-6.1456446976937249,-20.1313207694716,-10.218733832248855,-2.1462210939432262,-5.0174925851730787,-15.309070520644479,-1.8010982011039391,1.7114681436279606,-46.905612764335473,-25.378198245548724,-40.290904563054106,1.7152162102433819,7.4838190999878016,-18.46705085492772,-36.786114803159357,-42.400249094582755,-11.973806356233062,-23.062756316893402,-13.500271815739216,-6.0309837602962162,-26.474972589185946,-1.1592522546979327,-28.767124860726263,0.02541687486258375,1.6311160897463133,-28.324946407254824,9.7779119976273989,20.435502341343749,-31.188019436140259,-36.415702956590572,-31.374511455342592,-52.306155215065438,16.547883530060869,-7.9091006357199456,-10.228571213843438,-3.5960507887013136,0.79701562249944091,-44.106295597692366,-53.593602728706081],"y":[-0.36318288021553574,7.5610928015450822,3.1158056644190761,-1.2406925490892635,-3.0993430464014144,6.7958472296888948,-4.1853599657672778,13.530538134707081,-7.6329632903707108,-1.6641424303810264,-15.716700289322812,0.43165433380504092,-15.186500850079259,-3.9190651518879136,5.0438904056364438,13.984699109010004,30.446037088696858,0.055044865188616655,3.6433628901943185,-5.0501106025074689,-0.20792264868638236,-7.8217616520247315,18.81353892757912,-0.21272828673505564,0.70846830935470506,4.4186844322389662,12.441873741680817,11.512923217804341,12.016300784136606,14.392548771098522,12.375229654136639,1.0384425851544674,3.8266203402486711,2.4275544534565232,1.9161162035731241,11.155531708337596,11.827557397471525,18.384633422739174,3.9942382242733196,6.4213951482811575,0.075749721015056115,-4.4774753745666294,-5.1204234548459908,3.7767057897828153,-3.3127888903324925,2.9450612337145969,6.0989083639324981,16.014387047697706,13.462127020067484,16.744974180253966,29.536790850345724,-2.883932187863727,-5.9024712052652841,4.8574284923438418,1.922022323639015,-6.2363735004136567,-10.400053459461354,-17.928222978585456,-11.056510840392434,5.4908037874808979,-1.3035304987003562,5.4216908070601209,8.1258721968599055,18.156519362118896,-0.66617123179691595,8.6161429043449278,7.1620340792895973,15.643964958619215,0.73054859187036469,32.481210400847488,-8.8270149806710378,-1.7366213827114154,-4.4345962290986858,0.89799075290939456,-2.844603014411939,2.3580962557051932,11.427170771248209,3.8953397119979956,-10.800123841947972,-0.43654589436431135,12.481395976536012,1.340494277713228,10.619761651464904,2.4698006582761245,3.0073975261315291,1.0209563088441569,5.8719989731306663,6.6937721401906316,8.2280521180258237,11.063380723363384,-2.3562666436649407,8.3832331224514007,12.531151476716866,15.25906586398448,8.1184486014036192,8.8285605596102581,-3.7278107894774939,-11.204574903012389,-13.385473570856544,-3.2899258277393053,-12.180798092073433,-6.5861105238970667,-9.678427282470965,25.495620489750941,18.666425924918904,15.366571909052112,-12.597567676442475,-3.2594451036772707,3.5286120219367314,-1.0908838716036899,-3.4587333831559572,7.5465590347464051,9.3306985745785891,3.2233744566205025,16.420030849827182,5.9405077698261985,-22.295945841853278,-12.050086229659367,-7.2158385036519981,-5.5247413215661529,-0.48806886722124637,13.606367551154435,11.594694673128608,-4.4905657999170945,-5.1607386305719674,-12.382325784596715,1.6047632884039911,5.5587436813304345,8.459545368205827,24.19675481515462,17.909782991002086,4.9951008354792119,11.023953038465821,-11.81406161660302,-13.442253146182567,-6.9947322956306381,-9.6172217580197152,-7.8766001009899291,-7.13981371211697,6.7915916531207001,-2.4010707012759056,14.046243299448269,-7.1682229158211008,0.46161796853982945,-4.3395114680774061,-11.600171829351682,0.59272619136723215,24.948206353335785,35.456425526149637,39.71848512613812,-0.99301301486097315,12.182657900897107,6.5857161715728427,-15.402155975731159,-10.448624187788011,-17.653585832534262,-4.4384224892596853,-7.5526326353794735,7.8779233254518282,4.0926154313714553,12.547637597975262,14.033171069500392,8.1300613319829456,-11.198481280315653,-10.837920399347054,-8.1056202345117594,3.9952599639503772,-6.590340526079701,-6.3274537396470754,3.4711008665144911,3.8816143358984916,1.737276460481495,7.4867087885040569,7.1732618377457431,-8.6565615969831384,-7.8896839473871383,-8.5238031836951826,-2.1312594212270568,-8.1397138952042543,-7.1091042977921246,18.25355468666239,16.730571414229988,-7.1586087465533925,-1.8091532195718216,6.8876044540879047,-15.057392353895544,-0.86623171834369261,0.89063946109685288,-14.424861007095936,-9.4820293776370956,-3.6286959723299108,0.46432114239124656,-8.3990787923661276,18.901694311764672,9.3095031111216038,0.52377844335020418,16.117827716611764,6.6533710385078795,13.988542883852078,-4.2763217314113637,-6.1080947830590251,5.0021485653454052,-18.424179364918874,-6.7838504262086774,-10.732629139098824,-7.4447029240568785,-9.7880235370632427,-7.3998992191165307,-0.71343790249190475,5.6826042883238141,2.8422538355931883,5.4862697235405209,21.565773905472415,15.542002689726639,1.396391403286859,11.708201021425994,-3.1557390043009499,-10.11195897947224,-14.348169165864409,-3.1264774630129426,7.4099415389624017,-27.955557773220189,-3.6877329088322535,-10.099482167165178,-8.2710049463839699,-15.388655698635802,-16.082868848166232,-1.1682371347757967,11.14051639510823,-6.147883700855262,13.406915638622982,17.400532758506429,8.5588620027054905,-9.1334785071314855,2.1383850365601456,-21.647297880219043,9.7601436157766521,6.4912704688042115,-4.3512802733428098,3.4128008432440384,-3.8712295184742138,-23.505095287814978,-20.274320460760553,12.87307208678846,24.810921014948676,-2.665618032290904,-5.4234218749340339,-3.9561583381161691,-4.1563231431175316,3.6941118453028774,-6.3876105905874754,-6.000599428306745,-3.3337350254470555,1.9707518562568371,-10.242850122252829,0.76239678988438042,-8.6285807189862744,0.30078721850710283,3.8198755070297317,1.6454056937466095,8.1501770540907401,7.7381737388765819,11.554562369204405,-10.549231347120807,-13.459730342019308,-4.1839344880473437,12.929745124663601,14.482106427294232,30.899889657885332,15.61975172632568,11.134906017961985,5.527158933129102,-1.0543984991457032,3.7079343567654943,9.8508305495677089,8.3299816951003258,-5.9954537021083407,3.4831873769115949,4.8716135904179909,33.104167521253601,15.339772905769989,19.576957182540525,14.497545384962608,15.597027859665545,18.904627831944932,16.455759782969324,-1.5554608047882152,1.650461781406531,3.7450822645938997,0.28636906908405352,-11.520559989404493,5.2307126060489582,-0.027892120985151581,6.8354178700274506,4.6886807599943072,17.117453694719078,12.258977356298063,11.815377534154237,20.036374456434718,-15.033211646461329,-8.9903174399452297,2.8605093250890463,18.88301704263597,19.791719169902223,22.292644329586555,3.242809862245136,4.542804955423799,29.307816333464753,18.475221228918905,-1.6308938497030288,1.5126358653587335,-11.851181349981164,-1.7457745850653263,8.1763249102261977,5.894499927230326,7.407623979286087,7.8711248460095176,10.324818442278266,13.253550535583619,29.005896408129974,-13.31148788124316,-1.0286938174847142,-3.1181904287496502,-1.3446842847265619,6.9535288219072307,4.7890222349544596,-2.6634674303397556,-8.4912639024642491,-3.2546886412676277,4.3113307099437144,6.1101341057988003,6.5680166905896389,14.439760907869626,8.423762444472926,17.260005900736058,-11.62193454335142,-1.0649323978580076,15.420605190677223,5.340758643524957,4.0026149099896378,15.216119895968644,5.846333726727015,10.004233531341493,5.1952621746412486,-3.8827520462855154,17.261034821389288,9.4068790953679233,10.896159241732869,15.545340991536911,18.551058122955627,15.156746534976877,6.9818914801862348,9.3096457629588105,9.8711405512118091,8.9203947550814515,0.31298370261367375,-5.7713983342024386,-22.66466583890503,3.1772415021898199,5.6807623851178324,-5.0970944823581981,1.9119146623302599,-4.12502956679793,-3.5803289924539148,-9.9788413904819748,1.4313871402002381,-4.0020384653814354,2.8249827546457338,-1.6066598031037314,-6.6285968308498537,-6.164158011363841,-5.3725527279254228,11.545832458111626,17.129777933882,17.142252189985033,-0.55433071106704246,9.3436215635573703,-5.0054027545136952,-2.5374938615241494,-13.769916073042282,-6.4664084239915285,-16.283908326298636,-15.314822258003931,-3.0427276717866172,-21.424319079507359,-5.0929300165005005,3.1557748318004784,6.5421738639099321,24.343313944279398,-20.521912009830007,-10.613629662865,6.2995908482945611,-24.145937268065715,-16.110394272614442,-1.3647170511191244,3.4863319406522861,14.410524075455177,6.601641397595257,-5.9377985837995304,-0.32471214596830089,9.7760116759657709,-14.654633648262736,2.2286392980976677,-0.33306747105470397,14.389459351053191,21.222659058889114,3.7196881587238293,2.3800249864108651,11.903516680375986,23.206939336470025,19.676791306173886,-2.508238898877428,1.9284077251626779,-0.46247591218174688,19.070162508951412,8.670029505950744,-1.4002771382431802,8.3617878310617222,10.64381765075596,-1.3039010481560305,-15.960074211283224,-9.8368262973655405,-9.3421535622386767,-3.2448619899631517,10.403972708897461,14.522674336733825,21.126028918837424,19.050430102642583,-9.4748402865166472,-3.6308671355836459,-0.4360104447915662,-6.3936003935759267,-9.6309829526740831,-15.483496247597518,-9.2763683248428812,-18.363202881074134,5.3147307169330791,-16.76886234540542,-11.866542792856306,-12.917994556532523,-3.1811798876497264,-16.189940853490391,-4.3200658238064609,1.6062933938330086,-2.9395149816212527,-14.097477497714557,9.5557294875805958,6.2368615754140837,-2.5499926421603751,-7.1306711272467718,4.1404629675206124,-9.8852474131681074,11.71332427796022,20.552752696203317,12.152978579407645,15.205423221576075,12.624123968731022,3.0440635841394341,2.4195971763942783,6.562190648143039,-2.4118799851790316,-1.3777580341086972,13.980391087561472,15.318532750824367,19.976088357111184,15.770020836320331,-8.4327180479788506,-15.927793125498431,10.703190399821111,2.5657470512629352,13.611604350617117,13.959045516082101,20.175357603325388,-6.3361150660286425,-10.688403450338244,-1.4170545323149162,-3.7066929905922992,-3.3028815528488726,-5.38506127302057,-6.6339019109108737,-16.853000099483982,-14.425955462929608,-5.6959885151122798,8.7505804951015662,20.818811201153725,23.582200650430813,-3.2589344128029514,5.9277783980076206,21.148532501449054,8.203024402745994,4.256344928272223,4.6855096265224994,-4.0788094645946646,3.7163815504676632,-7.4499674595870795,-6.7577299213991449,-18.393016268206473,3.1967955524341742,12.962420811161552,3.6266179903245273,13.012277497280156,3.3763393091647531,-0.32436795973237448,-3.0176864086962225,-6.7507593581583603,-1.1774775167299976,-15.918926891844475,-8.8373119666412716,16.767998157877635,23.609386376986151,8.9143086001921237,26.709155242799735,17.091843315234481,-0.91712558126828025,5.723043136818303,2.4255402044977799,0.67272627682181074,10.709853063367886,13.384656711887926,-6.0446869175747722,5.9007728148524894,2.2033758380577608,13.366719333807636,-16.196097855812287,-12.513760577577687,-0.21223813627463903,8.2097116038783167,8.3232267356332414,16.36680789810157,-16.993498400292257,-5.6682471847807729,-3.7050960644646835,2.1528492213628518,-18.911947159728488,6.4300381113654703,3.2448856169527258,4.517768976906158,4.0303092570029158,-9.2024478875789235,3.4151679058775759,10.025569347115223,3.2340089895768567,2.1461365459460433,13.73868980640937,15.59319370015861,0.29338634638953043,0.72137382283559404,6.9632499709446405,2.9445852212697781,18.778319235248283,19.944724954052528,5.7149583206667627,11.003632050262066,21.530894242656945,-7.8876127295549319,6.374734162771901,12.742457583149742,2.7423654087553406,13.999887515916114,17.321842485341627,17.16759979185408,-9.5354334914301848,15.819974867251096,-0.26013761462372142,18.690364365883404,13.035456183822495,12.320258894883295,13.777906799201888,8.7756509261161302,0.66109626924459208,-8.99619267205847,-12.626570218824657,-3.7223120843683573,-7.2168966181193843,3.1122947910702372,28.071303685758263,24.046013433064626,-14.086846332861231,1.7425369330316089,-1.5137386136458217,-10.540763295940799,6.3600534932021624,-12.489086923694629,-4.9590284344929225,1.5643512142518177,-7.464679413969824,-1.9951115166780984,0.83755093683742665,21.051071124542784,5.5730417173109208,7.0525621124042033,14.642584746528714,12.58406812557616,12.956750859115976,20.814394895996777,9.2443741059855356,-5.4108042593645704,-0.98122688895602661,-0.90553152888226585,1.1815459925617497,0.99588546717843973,19.998926801567752,26.975631948921951,20.771889631358171,24.685305696234185,25.384956206293182,27.16733916519901,-17.577088855446341,-8.5980226697058217,-7.412709520048173,2.3413056436793034,8.4083896716792932,1.3237674605757992,2.1050927039460698,6.1747070300154512,7.7034094733052543,-15.47953136450637,-16.285762944817662,-0.1864762082867242,4.8319641985859372,-1.3391471130437533,17.802147511841245,15.639612813386035,3.3747596619371176,6.2187044799249884,1.6013035674140419,13.20050126312012,8.7735280100344095,-4.1911830216579311,-0.69420925975490388,2.2463310215374226,-14.385448230750518,-12.884900130534302,-14.952516919881806,30.446244737671353,29.89236645944283,32.869214507668651,-13.466035533961602,-9.3453273239644918,-2.3223658220732557,-1.0732208631845959,16.421554269280197,14.658913768953054,-10.874026621735634,-28.77392920666982,-19.50783251448922,-4.0408846520406421,-6.5793827154695421,-17.610362394860257,15.169214066334924,-7.0024856294262507,-10.781302803083053,2.0048190407054025,-4.4205582667547754,-12.043681155323647,9.3124237826991116,11.23897635203592,15.558742022023431,-8.091948565874878,-0.31768533879236271,20.557906571914522,24.65089195618738,20.558033162243916,-2.9114708976967631,-2.0845523580270253,-4.353290338195138,9.6360696156004373,14.082585900825789,-0.55351085181858162,-12.870866808225312,2.1934032270504558,2.7252165114459008,15.263877725520238,-2.1080005385184202,-12.926078773219167,-1.5084069684975931,-22.95014970217871,-28.821447465960478,-3.540284875778458,2.8838995880061504,6.0644407017545499,4.046836503760999,-9.1502114328821058,-8.8691758007329948,0.0065413435025442872,-2.7342328487988565,3.3510566334929401,-2.9778876620554513,-4.0480422040140862,-4.0840262288803864,13.681370634192257,14.390544704587942,23.707000092734891,30.478553467714956,-1.2429935121076323,8.1229433413539809,3.3280434108176067,14.511309144518391,8.506029392119272,10.822792132859362,14.25105495527764,-5.5118614355764501,-4.1746781536991735,-1.6211453107651077,-1.625899148994703,-5.5668815309530766,8.0994589128872878,11.696937094425735,6.5804476633640379,-4.3309017936226359,-20.151995354910866,-11.529138435115449,5.7732988258460729,9.7518911915862176,5.2675356114060561,4.8342874886694993,4.190879273855896,-0.2530893009218228,12.109440069779208,-15.409984952151,-5.1877987211761765,3.4856250827472959,6.7547589570181819,-15.861169544651389,-14.469935004362783,-14.879181573766273,-5.1564752402229308,-3.6660862528470091,-3.5537180292631141,-3.0464578236811426,-3.8689972628057867,-2.3545230062553535,-4.5107869883010885,-11.721804197432458,-0.90949808337804849,-1.6354873192239519,-8.5331980209892269,-13.024009851136796,1.1686381063470019,3.3224179245507797,6.7297787981651789,-3.6146678180772045,5.3644176894051547,-11.021523427022583,-9.9471744746623063,4.8927537784671484,-11.083982004321621,0.96841853776413012,8.4932570023685869,18.000407697749168,4.7234810594556382,7.5474907039023762,-1.776532349791212,-18.742099968048596,-2.2852303772972165,0.33999506896330911,0.64131150944594983,-16.777084214367282,-7.4323285921002773,10.220237256695324,11.511697064440796,15.306406847212251,-25.019596611588671,5.7478364128900665,-3.709825126045835,-9.3555785902055462,5.4333157103593699,-0.44100292001072,6.4105777685519554,28.75762868170219,16.487127386169274,-15.147071898056007,-15.90512027865948,-3.9507713549507484,0.68133061950905149,-13.798573321728123,-4.9120964390786952,3.9921721909472718,8.329789679788016,-8.3413404418962589,3.022917865131725,-5.2557245044070457,2.748462196473255,-6.2547506066655609,-11.349021337406947,3.465618109637735,-0.1422497593870122,1.4710878469823863,3.9406473320008635,7.4229231733368,2.3476439625864076,-1.9589748338579758,-5.3600035298459581,-15.814395635710968,10.602509857402593,2.345149189821361,4.9014068887787614,-17.848056515195911,-6.0226111279692143,10.435233656297873,11.37635800998352,10.612901003704545,14.458040808703448,28.730717195675936,4.2934161760231229,9.7410248597672169,-9.6279518548017897,-12.517931469082184,3.6156335445724852,-15.215850861806084,-20.033752484034032,-3.7913757466917102,-12.963075642608899,-6.0304261997569952,-8.2316523771857746,-1.9370032073852579,8.9938305510891343,24.93560269508205,-0.1337672860404191,2.2361604076005577,-6.7400794819271033,-11.858180771995391,-4.7954856680735949,6.9565974275522855,12.02923854487514,16.126958909342356,-11.118332475544033,13.835180428747535,0.51942970769432939,-0.07235391160909864,13.584332030510318,10.334827781126165,-3.2095202896068873,-3.3751044482755987,-4.5742138015691376,5.7444299408526476,2.080178213278451,8.9898354651358297,-2.1595792370392295,-1.5236464214204268,-13.562636957069516,-7.6696745656424197,-5.7155065907049334,-11.472449725988735,2.5612182000783985,-7.7875416337694015,7.9480694099608726,12.09662688175146,23.562440897440094,8.2932059497775406,11.907613514043403,27.079007896803496,3.6068577860650999,6.0894475131637673,9.1897994616463734,7.2749910731039931,21.345999107381022,9.3710162246103579,6.9286233275653988,13.463036923375366,25.335708374516088,-5.8231955341710107,11.644658568584894,4.9201879924120897,24.775499983920682,-27.555026720915215,-19.420879936065052,-13.350807071659339,-12.982760093118575,2.8279824461354486,7.4543400875644652,-1.3803577104912708,15.980825952081419,-13.203852136110267,6.7553975592384869,11.13469529734256,0.8429510850616071,5.3626348512720963,6.9163984269188417,12.931662628218968,11.041130447013167,-4.6957699307190603,-6.2267989241559301,-0.13227909839315666,0.098804790573835941,0.67489676451243907,9.5353517446029628,2.3513479353629916,-29.056552856433598,-19.105309147339007,5.2885238635768577,-4.7390318427767655,-8.9520083241244919,-5.429949507055464,13.238774605591376,18.184687161567293,4.9319138237636047,5.3489530288437175,-4.7050336526040235,-4.2503195559539888,-11.794445300172796,6.8208121040778176,-6.7994464880595231,23.202927140022375,16.326063308849665,12.982780617922781,8.998759473177298,10.909854323922813,15.206729427536395,-18.798767262226601,-14.437407157662161,-9.7892008711085996,-5.461327646751224,5.6898003769937349,-3.8623219267031681,2.4759814648207277,1.192015678355167,-19.321287771384966,0.63587827937569441,2.7309059352645844,-1.8332563091949121,-0.0028826114491992644,-7.3645754021527008,-4.9800785892981692,-6.6525313451597494,-8.4532141598473807,11.330204598527619,12.354298316218717,17.330804653531793,9.2612345922036834,14.796532001514814,2.7322042596047544,10.26706953809602,-2.3780626470342514,11.15665446306838,-2.038826746166436,9.2799914564877408,9.6662937540107059,17.961510870355628,12.248139909104719,13.996365048720444,20.457326829784222,10.84310902018961,6.6251702216043213,5.7938651846597757,-20.632816732315653,-5.6154557461825707,-9.9849137227047926,-1.4449135316693265,-7.1958879543272589,13.373584039166108,12.82246143475435,11.596720482109864,2.3245172219891517,-23.092964450557442,-16.155238989536024,-18.87101756050933,0.50805083668564022,-6.3318657771053761,17.626493902649162,13.576464674161626,6.8086090403830504,-5.6324636036060962,-12.78923151383969,7.3675939665114774,8.7377211151596583,4.1593572052468355,10.774446168085381,-2.5690375746415222,4.1800587780665053,6.5364500182964118,13.614224572357886,-9.7635074470203023,4.7176943651267687,-14.73317107731639,-5.6507931169352164,-13.831376351715079,2.0430830660212722,-11.51022513790136,8.1554078184722023,9.5152573361357309,7.8508227225196503,13.349224408767961,-16.782471535581728,-19.29485609974401,-8.0284263976425638,14.771887112312893,22.075179797885983,29.46330751311589,-7.1194227434618922,0.25852421889758531,3.5451924435867666,17.624529165863141,18.500293928203668,-8.5886981454061999,2.4897064283697081,2.5275180763692515,-3.8304723636391778,-4.6448447276074347,8.1234805605136167,7.1388485542920037,2.8991796858118088,17.137981572357482,13.14957566196548,-12.94165847036221,-5.7147282557825418,1.7733206985352246,4.9637874583491142,19.177644490932853,-9.8186126323151495,-7.0081968040164409,-8.0232009840094687,5.0604200158286439,-0.30647809969573347,-4.9057215859706114,12.472109557747757,9.3825841328771631,-7.6708441473344964,-9.1297563626510119,0.91747624012515427,-7.4887465208426329,6.1766789811538541],"z":[-7.0561319531423194,-8.0762771190413964,3.9391305503706908,3.5653495130143207,0.72352179755051271,-0.81697261081760308,-1.6180889365172806,-1.721154220613472,10.678881236276341,24.149655776921939,-5.8743520880737581,-3.6232474535811434,-3.3753396167931324,-11.393977348587022,-16.840951162690633,11.088879147807978,9.7557469352240869,-7.3266445926394921,-7.398589404691517,-2.3473444838323796,-3.7492125377826051,-4.3613431592264273,2.5456739068582328,6.2403201401104944,-2.5427229136030456,2.4733835277928904,2.8119677680090791,-4.9118001306843064,1.3241931826083952,-6.5241554516075944,-16.173685011432095,-9.3238141330159987,-1.3802534346611781,3.1648149858993109,4.80705128003195,-10.643217221699656,-3.0523921716607032,0.099239337576000411,6.5841423020883525,-10.777859934859904,-17.529093609939327,4.8999833824087222,15.597061263652181,9.1689634011896111,21.906282507979668,-3.7703326077027404,-4.3130136555327292,7.8780312955700182,5.4959654453052176,6.90244663350399,19.93141566352433,16.789616525957765,11.361825755656803,19.888139104388969,-3.7605944395421709,13.756059012592731,20.86980346915913,-0.81972053883704499,16.496404445682071,-0.24537677860213919,-5.8615484182290336,9.9182407919671682,-0.43706449289350974,4.4694777991974783,16.548707535963104,22.889327052336846,0.8843168370298371,7.5720208959191329,-5.6334804598369308,10.918169757354134,-2.7462710650615092,-9.3527519976386557,-1.8458191664021235,-1.4034307811032922,-1.8101317532208314,2.0799726523247513,14.833086606607406,11.631894964207561,23.037764358819615,-10.971636731903605,-9.9566861605505039,-1.9399457422509392,-3.6627466920016976,4.5461429320722404,10.323383634371121,3.539613208207157,13.077336170180697,7.8867421269301152,-0.36453882544611754,-0.30434427407651377,-1.7817424816003302,-3.8279900285119628,0.050049265615693485,1.2782793060870834,0.93936191087347232,20.336581965617825,0.72716067109205051,3.3335625363496852,3.1918429592955042,17.122956017699124,7.9852009310803584,5.615386859003471,16.398437811334048,2.9130271200788633,3.238008300692687,9.0670295342054033,1.7903421330681539,5.853780965125801,17.248994112214749,-7.343761415388423,-1.3595603886218162,5.3083973393265893,15.961428523846283,1.2025125732439237,-5.0651620151799186,2.9177893528302965,6.8993895831073511,-3.9746982540366282,10.255091023624232,-2.8608737387759717,-0.71386479237895994,6.9879570873044186,6.7398767701979354,6.8904964603864007,10.378511556800698,-11.225529775694172,-5.7231860326408341,10.428870097651892,8.8404016138693819,-0.2618554160649198,-0.051294106634502214,7.9820703224354528,18.352295023117282,-2.4465778727173868,3.5234818161338293,16.693056692362386,19.250403647355583,3.1837992270898678,11.341577627870594,-1.5282144106171758,1.6217938677624806,13.316676428187304,12.998187035687529,5.1597860747774122,27.456698917987644,7.6220565898927308,15.162534898610641,-11.388564189399379,-13.638411112102355,-13.149841739898406,6.7004794903504452,11.228577583244816,-5.8263790285191561,-6.2459155990253459,-4.0354130484820177,12.049515595317844,10.101253441468462,17.900558563838668,-19.321992005201842,-4.0110442429633109,-0.38459504782144399,-1.2089400917914508,0.76018188177241852,-7.3637151552203006,1.4750943819947795,5.6409961604051837,9.7612316474400735,-9.438202780178587,-2.499168192244412,6.5995456661239302,0.97053891903063305,10.847559150524374,6.8875827297828973,25.396594595905135,10.66389545943751,23.137315900405092,33.739432482486599,-13.522341528891468,0.75620992763941508,2.7205930926255326,-8.5574705125688908,-1.6040830662683829,-11.40207640817864,-7.668963465338833,3.6033336103000333,9.2985275079439766,-1.9288454735884437,-0.32742330729468699,-3.0991003969267998,-5.5140353467899166,-1.6660477644952809,-6.0000170134354045,-4.3910816966046982,3.8697189752471393,0.28388855516225342,0.57213100546076934,-7.5936414957508331,7.8500897400926153,10.664082712970032,3.2491931425028793,6.1513408935890999,11.922697986259859,3.5158992142036221,13.64735215042491,-6.6046615499132839,3.4342674744230997,-5.7839935372972295,6.1749775315649797,8.3190029124271643,-0.76719216830191506,4.1374447771134699,1.5454558669596989,-1.9246046800765186,-0.55905938019136325,-0.97556966454927196,7.5202092072415869,4.4212376923493837,-1.8559732012470522,7.6188943498855553,5.1463585896410011,-3.5004948794231563,9.4214370382887989,10.879276033903208,24.101584085387262,-3.7406378475927449,-4.0907323279003114,12.248245773010945,12.869380839050569,17.57138492653948,3.2082722804295214,-0.90863929578681957,-2.1148047133202672,9.5840769802162509,0.47631462855367335,11.099525489622616,5.7030788533439312,11.264939597141399,-7.8726888856877064,-9.1791640537184023,1.4999275511428267,-9.3608089969215786,-2.7622835359784181,10.218388169714203,9.5017211437951374,-7.2205234135818239,0.062782831029703723,13.828120828824703,20.668586920838639,12.589749975051605,13.694450186616953,-3.6261543448285005,13.904904470727036,16.338540932181118,1.9005778452971995,4.4607322055612579,14.929251847216356,-0.017940809627157325,-3.0462208579617212,-14.019328458136783,-5.0259760134866145,10.932831169207461,10.441535547842133,13.700796967037176,4.9351412708266675,5.4548292838470642,11.411340088194926,8.3956429095206744,10.273322172213192,22.714935306121877,0.99608953649527898,4.753144063103016,12.90592862107383,-5.8078554454323568,-8.1221623192299237,-3.8589315394378683,-3.6714045017726291,-12.797511484059035,-7.0283972736262408,-1.6691108143462057,3.9050617648695454,-4.1388628568257957,-3.5774736866027741,-2.1170730881012751,-1.5862849670719454,1.5830473758908807,2.5991282287188739,16.658982371371554,8.614399109304351,9.9240254757093744,0.033716077058618982,7.9040181512661363,2.8084447917251825,2.749529133250125,-3.0718451809170442,3.8291471424214394,0.94689005027752438,0.48606186731571788,-4.9388590001781978,6.1298202853930377,-4.9400722803249408,-21.401783736768078,-5.5378150116626479,9.3301793039209109,7.2253877444905061,2.4579929646620733,-21.66872147602826,-11.354196199882596,4.3955913411848986,4.8462425011112744,-20.431613086802987,-4.7694873876927559,-2.6126029932245372,0.92456629537041357,5.5757449164756343,-0.42703446664868033,3.156347611851551,-1.2642281820977301,-2.0114079140383749,5.7714374045666563,3.7930132318491618,-0.6947224472257868,-8.8335505277437942,0.83706316521199509,-10.983799485758519,13.961470055314926,-2.1651220423545636,-2.2272753810194241,-9.5411136987123228,9.2490544049455021,5.6399082420896027,0.45661563779465875,1.3259111161715951,3.0416008879168226,-4.4188482447355506,3.3109034480960577,-2.8722630458328635,4.4863601481116433,4.9520690153594877,10.076572599888728,3.4295236678313148,7.6810707511277467,-2.2059588599943143,-9.6182879270739825,1.5546753855131625,-6.5393352184517441,-1.5719311685489965,2.2573851856788654,1.6561459574031132,1.1534519609400962,15.892466444057542,8.9416573872483553,4.7320686994224861,2.8000836456370024,1.0728467409934546,-10.436971129574088,-9.2000361048704686,1.9603967953432986,8.116115106174421,6.3047785424041658,4.5246939478991184,-3.5218904357890168,-5.7466425015257947,-5.7407716378902229,7.0009731587156576,7.4250159839953103,12.568490443263727,-15.190495321844246,6.3766722734807608,6.8327108947223278,18.566623147384675,0.24072406435167529,-6.7607280890908577,-2.5302446678690318,0.99701714181533352,3.5972038743343635,-6.3407798788245797,3.5332721402356695,10.597982006194542,8.6053364065487745,8.7674315451031344,6.8959417512768217,-5.7785056816723452,-4.9196509306638552,4.2513901898409081,3.9835063052192381,16.799449968111652,-11.62494698482312,-15.740693762819978,12.282950310011495,16.078438418397948,6.7926783589237258,13.875758473865595,9.2853878265976224,3.738591477363145,2.394533561660277,-6.3212350737127583,-9.2157011683952899,-2.6997599207531424,-7.9175554506200108,4.1570106001112892,1.4193076907917335,9.2386688527709104,15.885838522406749,-10.377296191268906,-1.7256049016322541,8.3853859866784131,-1.0251217190062092,11.948803909609644,-7.9526667531948503,-3.9985127661800202,-3.2402927300043061,-1.4364248558816075,-7.5201624635480186,-5.7539058569754147,2.0048664530865015,-4.7994054509473782,8.047757012323018,-6.5259902915568455,7.0262985341978919,7.437282816527115,5.8500574426335055,6.6801216927182763,-4.714176601392742,-5.6862592091571242,-1.6951420324598916,-8.1659159598262789,-3.6284841947116355,-1.6257903095523654,0.70558228889510666,9.3277538806296914,9.0647932494930483,7.7014709216488235,-15.87206342974728,5.7996070666027331,-2.0060476689081979,4.7456223202231742,2.9453566821391184,1.7339040267729435,-8.2119111173290094,3.9542478853514669,16.325155729851453,17.36165759275509,4.9184732908435,8.3871093832455159,3.9597610232488574,4.1391412557172966,0.19471348393987409,-2.7025032793077632,-5.6878835210531413,1.9541916288507044,6.0990629832242025,33.534507645097634,-5.6031572556465949,-0.94563954630547919,-0.32893984509432228,7.0688188998485879,6.4929269418833808,-6.3942698581870046,1.6551678384679425,-8.6076117487366197,5.2768500129991933,6.7893538987712239,-1.9488009471296375,2.2055860580646405,-6.3132269907168546,18.697966120029957,-1.5697404625231408,0.60561478974790994,17.645127248183336,9.6833149304306847,-1.1687208979603887,-3.4180445307442433,-7.2101503895111287,2.0798099016736535,6.4668656726092442,12.352015207533722,6.1116701056783338,4.6319920339682561,17.06352993403107,4.5128571748150863,9.6849898430751811,6.6032576208470415,4.4913261918427221,16.150992163869343,-4.1298986968133997,3.565350699955145,-5.4369268325558213,-1.3810805471841507,19.694614018616541,3.6278042663054801,6.3617505294502763,-3.0145106234018755,4.3279464177573486,-1.9820571112634617,-9.2463495528463326,1.8009931780314257,0.17707876954364199,10.832464391209397,4.9832729295279812,2.0373537984034167,6.5785953564836852,-12.972449613765775,-21.964870977670309,5.933882521754791,1.2593223248907581,9.4827171246452959,20.450147419060499,15.751096213460634,-13.006113456046341,7.7321933982284916,-6.1113896566109851,-1.3150699668259531,-1.1040503316194712,3.3752661695748771,0.22565399345469675,2.4372074263131327,-11.837197052199333,-12.873325549900683,-2.094331288742219,0.14805070485722457,2.8718969431953347,-1.1674070361081965,0.031751144336378802,2.6476788943691427,8.5443278139307424,14.792006515588591,-0.90287902589126512,-13.295619549168995,-3.02531277662312,11.080249844722131,12.665435681571291,-2.4414371226357381,15.186561515007739,18.531401761751134,2.2669352497001252,2.6361147994456542,9.9326327110631834,-1.2240096499162687,1.461020458439837,2.0027881562362455,-1.317426781843049,7.9019827009374977,5.6081879631622469,-0.38805833720710564,3.0546297744180237,8.324227559576574,-2.8486236869998569,3.4781466267787851,-0.52600180374079197,12.264005578677015,-0.88234267454011628,17.806043325611441,11.498763318914634,20.0584908082031,3.6395765592020739,-3.8656513031534661,9.7667037012860405,7.8593296658571816,18.32082889074707,-17.503329563226636,-11.963147817839936,5.2741968051645038,7.8785379314016044,11.847671667119293,6.3643775719515334,14.786732334893349,-1.5135529828102392,-1.6364356294433455,-3.4571162991807318,6.1484543680249306,1.7438800838989712,18.101707932569916,10.332663631166715,12.83314915721162,6.0338421117219649,-10.935545474282332,-8.8857669765847138,14.747898659775201,10.789743550281848,25.110973177212667,7.7740873856676265,10.048740905377022,2.8886495653177686,-10.854355209668684,13.310422378824169,7.7918430929232168,19.696808399271539,-0.48454089458296934,8.8004326433671611,-2.1957093353560095,1.2454981173052013,-1.4508276593393572,3.5544491926951101,6.0357648279441163,3.2567454306543566,13.482659151123164,5.9161195195518177,11.845533727975328,16.246248521538735,14.845926957472066,18.673139115481828,-7.9171458502572856,6.1824893644516914,-19.27127529638479,-10.552549572180482,3.7270866265328877,4.300168872624603,0.5625256381526037,-3.9837647535116658,9.1239463583190332,-0.11551810052117614,8.5856723087254814,1.0113816144462344,15.237493246725212,-15.646836663121007,-12.441365515432446,-1.9767300259264815,-4.3646165343062995,2.9037641494944175,1.8112823860476115,6.3776203972836409,-7.2200123626472319,2.429048765426669,17.185321315491226,27.361402686642911,-3.3406410373668214,-2.0198912656508132,11.386347840556461,-2.7366105740990059,-3.9337658952795578,15.123515196889604,-11.542301380758504,-12.934576067716657,-5.4933919000740064,-5.1196666454554745,-7.9881280344773469,-0.11136279650263362,6.6725446038444352,-0.79439923548244107,6.2716433044196185,14.896140148313956,2.310490641742224,16.36761952317644,-8.3681931426907603,-12.375254695147225,-5.5557632145506997,0.18590028138176898,3.8002465459544803,20.316994115046132,-6.8904343268304196,17.02470514338555,19.61821491194533,4.1562408170260063,-7.9718446013896198,-15.216637123127898,-1.23171754095036,10.452840040729864,23.521655079109618,4.794965367365581,8.276289145999165,12.566383358162623,-1.1671859042668082,9.3499437453120837,2.1173947611015689,1.4584189924611299,-4.0872325245340058,5.9705109393623434,8.840893325688846,11.240621309980471,3.3142386989635755,-1.0457810250290136,8.199219212923758,2.4745495814082004,11.832974168679614,4.6727336993046062,7.5665215988218844,13.486741283195856,7.7131637343639632,6.9727619319123981,-10.208777604353237,1.5512658752329314,-2.0161618679795872,5.7917931966695706,10.308070902641113,-2.2515558196912719,-7.6133701022177656,1.2168870444191517,5.6022304354406103,5.5513634341355012,7.4006472408055926,-6.1585249822674397,-0.41566003821715525,-6.2391250506547937,6.5670470036312638,-9.5592545823809338,-0.83615375414421089,-0.28284485208169441,2.0643852523743371,17.454522566082865,20.175858533378914,5.8253522473684569,7.240601696025843,-4.2005449933873811,-2.2314333868926308,-3.5422171773222075,18.095290476725278,-6.0379241210451173,3.372058076808309,9.1238219259146955,1.0973251615585993,8.1707799953826363,16.883799594221362,-7.6443410687695961,2.7467888674043714,4.1514166283619556,0.68888305630388391,16.532570853495873,-10.170361670691818,-2.7865101892081334,1.1488537127227498,7.2425304995626973,8.3340254795890694,10.56601912566228,-1.4662189721687717,-6.2438492867235835,2.1356892608172791,4.3809771089114111,0.33055702359224015,9.393581876826687,-10.923011544210722,0.29733005345706065,1.566069167715612,7.6954867918413887,10.82858958451415,-3.0827120627204558,3.349620409624392,13.733516594628743,9.4150157992808321,23.667896577491696,-6.7654495756098356,9.2787277921956939,14.979506399321002,19.099301550295689,-5.4283254380946975,-4.4943749162001234,0.68139388701122461,-8.8112078032636845,15.468559514397185,-12.301852656214086,-15.942445496459754,-2.7129109970896312,-4.8535256868398857,-6.7655984152826933,15.279076597332399,17.478177574942716,0.33172794418301293,4.8666175310395658,4.0814076178173462,8.1395474004205841,1.4056485983425038,-11.723160839918629,-0.85720740873462697,-18.481550104147146,-0.57015539459194287,-0.11485965292369238,3.3326936414005361,9.2933688343119254,-6.3397482883249907,3.4303779525939828,-15.477086906317114,15.305523347660095,14.636321072193084,18.22920229048583,11.067707050431899,-1.216632836012344,3.398489032660895,2.6026178393884267,-0.07657829123781594,-5.9550925681313283,4.596957261424893,5.8465495351944208,2.1097709117496035,-2.9108203376009056,3.9109584649348177,17.659804409474273,-12.221497901060477,-13.194353395063208,-6.4307811799114161,0.19672126511371321,0.68257170165429326,-1.4483932863921964,8.3046058337956143,2.6488388348648551,30.622233988237063,10.04598582870854,-8.0324231066558713,-4.1055051622192416,15.431923174113544,4.226212576543503,-5.6250425981650585,0.162047224118025,7.5825025202060061,2.660752092853282,1.1161933207522412,-5.39036565738569,9.0308561203369013,0.097221113181623298,2.053943251813914,4.8582448590119105,-4.3674652789624648,-0.63598599688661994,5.1839412320825105,8.2416275926813327,7.7980484940962196,19.124091883698714,-3.1643508609146287,11.950782340443943,18.824684487173855,-0.65648574979324414,-3.6006237981928648,-4.9389183082625685,5.5470752415381073,15.977052748327305,-4.6657478904922121,-2.0969959491377232,12.182143119447289,-8.4272614897489326,-5.8707687893284595,6.3827635634374937,6.2744442008780261,3.086008033948838,1.2156561526315,-6.3549500728109285,-7.3290433506370372,3.8620746798822996,-3.1919199069832138,19.410394387326022,12.423127277197853,15.085613549468746,20.704483964707908,-12.621132811926799,-9.0343475444287495,-6.9062857212106419,-7.9209864969788653,0.55580259096752593,-1.5422693615847385,5.3199722398714524,-31.422487402919401,-28.503537705038028,-19.071045741802397,8.055692196732311,8.2905141303904966,5.06774028620426,5.3901981918590005,13.24345463104285,-8.4160319345213725,-8.5913455061735284,4.051671103572974,-2.3185727846678583,4.653415460141213,19.660454137302139,5.56568354577532,-1.1346319660187705,-4.420474288133116,8.7837464299360093,15.205321039174198,11.562527812285419,18.023871751590221,5.4194299472621843,19.555594043494942,-3.3911493751252402,-0.24841678696532865,2.3958359298037792,-7.98598794582922,-4.6547165795041447,-4.8854342132864099,-1.2360708989073994,7.0204104062585007,-5.5450389137996829,-1.5997752997755408,2.1604974756004465,-9.3583292405276062,-0.56479094590016576,21.17621960659968,3.354778242336574,4.3046838816453867,-2.6766013745890858,2.1487271386666618,-3.4325467170453861,9.5592314111096783,0.85137144628591832,0.57797783318239737,-2.7385428434289709,9.4534584811534845,-1.160360853507979,1.0983545728950215,12.173425532744266,13.945966291914377,-14.588211118269953,-7.538552891185736,0.26335355858413256,4.3164524052787341,5.4668901829786751,-6.6860439134454719,-14.71890004126357,-9.8774815181984064,-19.891884532365228,-15.844801053947656,-7.2243649118645363,5.640989842292484,4.4347302491860248,8.394136497023462,11.975724407617463,-7.9123601483758783,-5.031412667458385,-9.6378199085235448,-9.8903038295567853,5.9898557342082546,0.68010998352218965,12.542234378604499,-1.8062720446830971,14.608808691852442,18.276738060087016,13.043830055595221,15.776065255274547,4.565515909385133,-5.4758539224208649,-5.4452609909277534,-0.73308848522084025,-6.5782456212283291,0.20560291289498919,-9.5984740621396671,3.5022599725074413,7.1433452888094067,-2.1678684805596786,-4.4663904317602583,-1.7547025734744977,14.996608748390184,9.498276738176866,2.289752819078203,-5.23440131754018,3.0861607126207855,10.505059934356096,20.350557743286394,15.258849946861993,6.0625835156633032,14.131549908821377,10.909787959426975,5.1017289458523765,3.1309717858646744,10.916498370615235,-1.4020859148315314,-2.0218277788356267,1.3606772845123376,-0.28554744532727377,9.8009743569400207,-4.3409674327321897,-11.206833192153287,-14.146880247948436,8.758462168041973,16.232104688358017,37.556777383990159,-4.5903780982688822,4.2143167766380376,-14.031886067927386,2.3973146289815919,-3.9410798603462514,4.6830163915574827,9.8315110694448524,13.715533995484238,-1.9675107437565174,8.5737884478124791,0.31679474245956035,11.7632489296514,-7.0574859180339784,-4.6733842390453901,0.28715885825716836,11.475412236717959,1.3157500582542427,5.9744659008771288,-11.903908946075767,-5.9676663337429883,-3.4660028128475688,-10.311202731839826,-13.733922423509854,-7.7807746341687745,5.1703012331697256,9.811756074372413,18.183883223474933,9.0022018578071084,21.598298603323258,-3.2121169032539982,5.2142920934350068,11.046770981881599,-6.5903358102588196,-9.2693764786847073,5.3043131030801751,-11.291063648463496,-14.182074757984589,-3.1342343425568764,-2.0278999204160413,0.52582913140846665,-7.2948011846350278,-2.6928430096232243,0.36513904507794837,5.6072180062143167,5.399341060803482,8.0337394298965208,-2.5752630114048745,-12.683398991726261,1.6874693505589788,-4.18473019333255,2.2102365994571369,9.9832735836022515,4.2323471019999026,4.6098885299387149,6.1836335287756166,5.3891460798111632,7.3926181705275127,20.665276857226509,-4.1559634735286393,-4.2897253833578679],"mode":"markers","marker":{"color":"rgba(0,178,238,1)","size":2.5,"line":{"color":"rgba(0,178,238,1)"}},"type":"scatter3d","name":"ao","textfont":{"color":"rgba(0,178,238,1)"},"error_y":{"color":"rgba(0,178,238,1)"},"error_x":{"color":"rgba(0,178,238,1)"},"line":{"color":"rgba(0,178,238,1)"},"frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

<p class="caption">(\#fig:unnamed-chunk-27)3D graf skórů prvních třech hlavních komponent.</p>
</div>

Vykresleme si i průběh kumulativní vysvětlené variability. 


```r
data.frame(x = 1:15,
           y = pca.fd(XXfd, nharm = 15)$varprop |> cumsum() * 100) |>
  ggplot(aes(x, y)) + 
  geom_point(col = 'deepskyblue2') + 
  geom_line(col = 'deepskyblue2') + 
  theme_bw() + 
  labs(x = 'Počet hlavních komponent',
       y = "Kumulativní vysvětlená variabilita [v \\%]") + 
  geom_hline(aes(yintercept = 90), linetype = 'dashed', col = 'grey2') + 
  scale_x_continuous(breaks = 1:15) + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-28-1.png" alt="Kumulativní vysvětlená variabilita [v %] proti počtu hlavních komponent." width="576" />
<p class="caption">(\#fig:unnamed-chunk-28)Kumulativní vysvětlená variabilita [v %] proti počtu hlavních komponent.</p>
</div>

```r
# ggsave("figures/kap2_PCA_nharm.tex", device = tikz, width = 5, height = 3)
```

Také nás zajímá průběh prvních třech funkcionálních hlavních komponent.


```r
## Looking at the principal components:

fdobjPCAeval <- eval.fd(fdobj = data.PCA$harmonics[1:3], evalarg = t)
df.comp <- data.frame(
  time = t, 
  harmonics = c(fdobjPCAeval), 
  component = factor(rep(1:3, each = 256))
  )

df.comp |> ggplot(aes(x = time, y = harmonics, color = component)) + 
  geom_line(linewidth = 0.7) +
  labs(x = "Frekvence", 
       y = "Hlavní komponenty",
       colour = '') +
  theme_bw() + 
  scale_color_manual(values=c("#127DE8", "#4C3CD3", "#12DEE8"))
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-29-1.png" alt="Hlavní komponenty." width="576" />
<p class="caption">(\#fig:unnamed-chunk-29)Hlavní komponenty.</p>
</div>

```r
# ggsave("figures/kap2_PCA_components.tex", device = tikz, width = 5, height = 3)
```

Můžeme se také podívat na vliv prvních třech hlavních komponent na průměrnou křivku. Vždy je od průměru přičten nebo odečten dvojnásobek hlavní komponenty (škálovaný rozptylem). 


```r
library(plyr)
freq3 <- seq(1,256)
fdobjPCAeval3 <- eval.fd(fdobj = data.PCA$harmonics[1:3], evalarg = freq3)

df.comp3 <- data.frame(
  freq = freq3, 
  harmonics = c(fdobjPCAeval3), 
  component = factor(rep(1:3, each = length(freq3)))
  )

fdm3 <- c(eval.fd(fdobj = meanfd, evalarg = freq3))

df.3 <- data.frame(df.comp3, m = fdm3)
df.pv3 <- ddply(df.3, .(component), mutate, 
                m1 = m + 2*sqrt(data.PCA$values[component])*harmonics,
                m2 = m - 2*sqrt(data.PCA$values[component])*harmonics,
                pov = paste0("Komponenta ", 
                             component,", Vysvětlená variabilita = ",
                             round(100*data.PCA$varprop[component], 1), 
                             ' \\%'))

df.pv3 |> ggplot(aes (x = freq, y = m)) +
  geom_line() +
  geom_line(aes(y = m1), linetype = 'solid', color = 'deepskyblue2') +
  geom_line(aes(y = m2), linetype = 'dashed', color = 'deepskyblue2') +
  labs(x = "Frekvence", 
       y = "Log-periodogram",
       colour = 'Komponenta') +
  theme_bw() +
  theme(legend.position = 'none') + 
  facet_wrap(~ pov, nrow = 1)
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-30-1.png" alt="Vliv komponent." width="1152" />
<p class="caption">(\#fig:unnamed-chunk-30)Vliv komponent.</p>
</div>

```r
# ggsave("figures/kap2_PCA_impactofcomponents.tex", device = tikz, width = 9, height = 3)
```

Podívejme se ještě na rekonstrukci původního log-periodogramu pomocí hlavních komponent. Nejprve uvažujme pouze průměr. 


```r
meanfd <- mean.fd(XXfd)
fdm <- eval.fd(fdobj = meanfd, evalarg = t)
colnames(fdm) <- NULL
scores <- data.PCA$scores
PCs <- eval.fd(fdobj = data.PCA$harmonics, evalarg = t) # vyhodnoceni

df <- data.frame(dfs[1:256, ], reconstruction =  fdm, estimate = "mean")
p0 <- ggplot(data = df, aes (x = time, y = value)) +
  geom_line(color = "grey2", linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = reconstruction), colour = "deepskyblue2", linewidth = 0.6) +
  labs(x = "Frekvence", y = "Log-periodogram") +
  theme_bw()
print(p0)
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-31-1.png" alt="Původní křivka a průměr." width="576" />
<p class="caption">(\#fig:unnamed-chunk-31)Původní křivka a průměr.</p>
</div>

Konečně se podívejme na vybraný počet hlavních komponent a příslušnou rekonstrukci.


```r
for (k in 1:20){
 df1 <- data.frame(dfs[1:256, ], 
                   reconstruction = fdm + c(PCs[, 1:k] %*% 
                                              t(scores[, 1:k]))[1:256],
                   estimate = paste0("comp", k))
 df <- rbind(df, df1)
 p1 <- ggplot(data = df1, aes (x = time, y = value)) +
  geom_line(color = "grey2", linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = reconstruction), colour = "deepskyblue2", linewidth = 0.6) +
  labs(x = "Frekvence", y = "Log-periodogram") +
  theme_bw()
 # print(p1)
}

df |> mutate(estimate = factor(estimate)) |>
  filter(estimate %in% c('mean', 'comp1', 'comp2', 'comp3', 'comp9', 'comp20')) |> 
  mutate(estimate = factor(estimate, levels = c('mean', 'comp1', 'comp2', 'comp3', 'comp9', 'comp20'))) |> 
  ggplot(aes (x = time, y = value)) +
  geom_line(color = "grey2", linewidth = 0.5, alpha = 0.5) +
  geom_line(aes(y = reconstruction), colour = "deepskyblue2", linewidth = 0.7) +
  labs(x = "Frekvence", y = "Log-periodogram") +
  theme_bw() + 
  facet_wrap(~estimate, ncol = 3, nrow = 2)
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-32-1.png" alt="Původní křivka a její rekonstrukce." width="864" />
<p class="caption">(\#fig:unnamed-chunk-32)Původní křivka a její rekonstrukce.</p>
</div>

```r
# ggsave("figures/kap2_PCA_reconstruction.tex", device = tikz, width = 9, height = 6)
```

## Materiály pro Kapitolu 3

Tyto materiály jsou převzaty z Kapitoly \@ref(aplikace2). 

## Materiály pro Kapitolu 4

V této sekci uvedeme podpůrné grafy pro čtvrtou kapitolu diplomové práce.

### Maximal margin classifier

Nejprve simulujeme data ze dvou klasifikačních tříd, které budou lineárně separabilní. 


```r
library(MASS)
library(dplyr)
library(ggplot2)

set.seed(21)

# simulace dat
n_0 <- 40
n_1 <- 40
mu_0 <- c(0, 0)
mu_1 <- c(3, 4.5)
Sigma_0 <- matrix(c(1.3, -0.7, -0.7, 1.3), ncol = 2)
Sigma_1 <- matrix(c(1.5, -0.25, -0.25, 1.5), ncol = 2)

df_MMC <- rbind(
  mvrnorm(n = n_0, mu = mu_0, Sigma = Sigma_0),
  mvrnorm(n = n_1, mu = mu_1, Sigma = Sigma_1)) |>
  as.data.frame() |>
  mutate(Y = rep(c('-1', '1'), c(n_0, n_1)))
colnames(df_MMC) <- c('x1', 'x2', 'Y')
```

Nyní vykreslíme data.


```r
p1 <- ggplot(data = df_MMC,
             aes(x = x1, y = x2, colour = Y)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title.align = 0.5) + 
  scale_x_continuous(breaks = seq(-2, 6, by = 2), 
                     limits = c(-3.5, 6.5)) + 
  scale_y_continuous(breaks = seq(-4, 8, by = 2),
                     limits = c(-2.5, 7)) + 
  scale_colour_manual(values = c('tomato', 'deepskyblue2')) + 
  labs(x = '$X_1$', y = '$X_2$', colour = 'Klasifikační\n třída')

p1
```

<img src="14-Application_5_files/figure-html/unnamed-chunk-34-1.png" width="672" />

Natrénujeme klasifikátor a vykreslíme dělicí nadrovinu společně s podpůrnými vektory.


```r
library(e1071)

clf <- svm(Y ~ ., data = df_MMC,
                 type = 'C-classification',
                 scale = FALSE,
                 kernel = 'linear')
```


```r
df_SV <- df_MMC[clf$index, ]
p2 <- p1 + geom_point(data = df_SV, col = 'grey2', alpha = 0.7,
                      size = 2)
p2
```

<img src="14-Application_5_files/figure-html/unnamed-chunk-36-1.png" width="672" />

Dokreslíme dělicí nadrovinu.


```r
# vektor koeficientů
w <- t(clf$coefs) %*% clf$SV

slope <- - w[1] / w[2]
intercept <- clf$rho / w[2]

p3 <- p2 + 
  geom_abline(slope = slope, 
              intercept = intercept,
              col = 'grey2', linewidth = 0.7, alpha = 0.8) + 
  geom_abline(slope = slope, 
              intercept = intercept - 1 / w[2],
              col = 'grey2', linewidth = 0.5, alpha = 0.8,
              linetype = 'dashed') + 
  geom_abline(slope = slope, 
              intercept = intercept + 1 / w[2],
              col = 'grey2', linewidth = 0.5, alpha = 0.8,
              linetype = 'dashed') + 
  geom_point() + 
  geom_point(data = df_SV, col = 'grey2', alpha = 0.4,
                      size = 2)

p3
```

<img src="14-Application_5_files/figure-html/unnamed-chunk-37-1.png" width="672" />

```r
# ggsave("figures/kap4_MMC.tex", device = tikz, width = 6, height = 4)
```


### Support vector classifier

Nejprve simulujeme data ze dvou klasifikačních tříd, které budou lineárně neseparabilní. 


```r
set.seed(42)

# simulace dat
n_0 <- 50
n_1 <- 50
mu_0 <- c(0, 0)
mu_1 <- c(3, 4.5)
Sigma_0 <- matrix(c(2, -0.55, -0.55, 2), ncol = 2)
Sigma_1 <- matrix(c(2.75, -0.3, -0.3, 2.75), ncol = 2)

df_MMC <- rbind(
  mvrnorm(n = n_0, mu = mu_0, Sigma = Sigma_0),
  mvrnorm(n = n_1, mu = mu_1, Sigma = Sigma_1)) |>
  as.data.frame() |>
  mutate(Y = rep(c('-1', '1'), c(n_0, n_1)))
colnames(df_MMC) <- c('x1', 'x2', 'Y')
```

Nyní vykreslíme data.


```r
p1 <- ggplot(data = df_MMC,
             aes(x = x1, y = x2, colour = Y)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title.align = 0.5) + 
  scale_x_continuous(breaks = seq(-2, 6, by = 2),
                     limits = c(-2.75, 6)) +
  scale_y_continuous(breaks = seq(-4, 8, by = 2),
                     limits = c(-4.5, 8)) +
  scale_colour_manual(values = c('tomato', 'deepskyblue2')) + 
  labs(x = '$X_1$', y = '$X_2$', colour = 'Klasifikační\n třída')

p1
```

<img src="14-Application_5_files/figure-html/unnamed-chunk-39-1.png" width="672" />

Natrénujeme klasifikátor a vykreslíme dělicí nadrovinu společně s podpůrnými vektory.


```r
clf <- svm(Y ~ ., data = df_MMC,
                 type = 'C-classification',
                 scale = FALSE,
                 kernel = 'linear')
```


```r
df_SV <- df_MMC[clf$index, ]
p2 <- p1 + geom_point(data = df_SV, col = 'grey2', alpha = 0.7,
                      size = 2)
p2
```

<img src="14-Application_5_files/figure-html/unnamed-chunk-41-1.png" width="672" />

Dokreslíme dělicí nadrovinu.


```r
# vektor koeficientů
w <- t(clf$coefs) %*% clf$SV

slope <- - w[1] / w[2]
intercept <- clf$rho / w[2]

p3 <- p2 + 
  geom_abline(slope = slope, 
              intercept = intercept,
              col = 'grey2', linewidth = 0.7, alpha = 0.8) + 
  geom_abline(slope = slope, 
              intercept = intercept - 1 / w[2],
              col = 'grey2', linewidth = 0.5, alpha = 0.8,
              linetype = 'dashed') + 
  geom_abline(slope = slope, 
              intercept = intercept + 1 / w[2],
              col = 'grey2', linewidth = 0.5, alpha = 0.8,
              linetype = 'dashed') + 
  geom_point() + 
  geom_point(data = df_SV, col = 'grey2', alpha = 0.4,
                      size = 2)

p3
```

<img src="14-Application_5_files/figure-html/unnamed-chunk-42-1.png" width="672" />

```r
# ggsave("figures/kap4_SVC.tex", device = tikz, width = 6, height = 4)
```

Nakonec přidáme popisky k podpůrným vektorům.


```r
df_SVlab <- cbind(df_SV, data.frame(label = 1:dim(df_SV)[1]))
p4 <- p3 + 
  geom_text(data = df_SVlab, aes(label = label), colour = 'grey2',
            check_overlap = T, size = 3, vjust = -0.55, hjust = 1)

p4
```

<img src="14-Application_5_files/figure-html/unnamed-chunk-43-1.png" width="672" />

```r
# ggsave("figures/kap4_SVC.tex", device = tikz, width = 6, height = 4)
```


#### Změna šířky tolerančního pásma při změně hyperparametru $C$

Podívejme se ještě na změnu tolerančního pásma v závislosti na hyperparametru $C$.


```r
C <- c(0.005, 0.1, 100)

clf1 <- svm(Y ~ ., data = df_MMC,
                 type = 'C-classification',
                 scale = FALSE,
                 kernel = 'linear',
           cost = C[1])

clf2 <- svm(Y ~ ., data = df_MMC,
                 type = 'C-classification',
                 scale = FALSE,
                 kernel = 'linear',
           cost = C[2])

clf3 <- svm(Y ~ ., data = df_MMC,
                 type = 'C-classification',
                 scale = FALSE,
                 kernel = 'linear',
           cost = C[3])

df_SV <- rbind(df_MMC[clf1$index, ] |> mutate(cost = C[1]),
               df_MMC[clf2$index, ] |> mutate(cost = C[2]),
               df_MMC[clf3$index, ] |> mutate(cost = C[3]))
p2 <- p1 + geom_point(data = df_SV, col = 'grey2', alpha = 0.7,
                      size = 2) + 
  facet_wrap(~cost)

# vektor koeficientů
w <- t(clf1$coefs) %*% clf1$SV
slope <- - w[1] / w[2]
intercept <- clf1$rho / w[2]

df_lines <- data.frame(slope = slope,
                       intercept = intercept,
                       lb = intercept - 1 / w[2],
                       rb = intercept + 1 / w[2],
                       cost = C[1])

# pro clf2
w <- t(clf2$coefs) %*% clf2$SV
slope <- - w[1] / w[2]
intercept <- clf2$rho / w[2]

df_lines <- rbind(df_lines,
                  data.frame(slope = slope,
                       intercept = intercept,
                       lb = intercept - 1 / w[2],
                       rb = intercept + 1 / w[2],
                       cost = C[2])
                  )

# pro clf3
w <- t(clf3$coefs) %*% clf3$SV
slope <- - w[1] / w[2]
intercept <- clf3$rho / w[2]

df_lines <- rbind(df_lines,
                  data.frame(slope = slope,
                       intercept = intercept,
                       lb = intercept - 1 / w[2],
                       rb = intercept + 1 / w[2],
                       cost = C[3])
                  )

p3 <- p2 + 
  geom_abline(data = df_lines,
              aes(slope = slope, 
              intercept = intercept),
              col = 'grey2', linewidth = 0.7, alpha = 0.8) + 
  geom_abline(data = df_lines, 
              aes(slope = slope, 
              intercept = lb),
              col = 'grey2', linewidth = 0.5, alpha = 0.8,
              linetype = 'dashed') + 
  geom_abline(data = df_lines, 
              aes(slope = slope, 
              intercept = rb),
              col = 'grey2', linewidth = 0.5, alpha = 0.8,
              linetype = 'dashed') + 
  geom_point() + 
  geom_point(data = df_SV, col = 'grey2', alpha = 0.4,
                      size = 2) + 
  theme(legend.position = 'none')

p3
```

<img src="14-Application_5_files/figure-html/unnamed-chunk-44-1.png" width="672" />

```r
# ggsave("figures/kap4_SVC_C_dependency.tex", device = tikz, width = 8, height = 3)
```

### Support vector machines

Pro ilustraci této metody využijeme data `tecator`, kterým se podrobně věnujeme v Kapitole 11. 


```r
# nacteni dat 
library(fda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ddalpha)

data <- ddalpha::dataf.tecator()

data.gr <- data$dataf[[1]]$vals
for(i in 2:length(data$labels)) {
  data.gr <- rbind(data.gr, data$dataf[[i]]$vals)
  }
data.gr <- cbind(data.frame(wave = data$dataf[[1]]$args),
                 t(data.gr))

# vektor trid
labels <- data$labels |> unlist()
# prejmenovani podle tridy
colnames(data.gr) <- c('wavelength',
                       paste0(labels, 1:length(data$labels)))
```



```r
t <- data.gr$wavelength
rangeval <- range(t)
breaks <- t
norder <- 6

bbasis <- create.bspline.basis(rangeval = rangeval, 
                               norder = norder, 
                               breaks = breaks)

curv.Lfd <- int2Lfd(4) # penalizujeme 4. derivaci

# spojeni pozorovani do jedne matice
XX <- data.gr[, -1] |> as.matrix()

lambda.vect <- 10^seq(from = -2, to = 1, length.out = 50) # vektor lambd
gcv <- rep(NA, length = length(lambda.vect)) # prazdny vektor pro ulozebi GCV

for(index in 1:length(lambda.vect)) {
  curv.Fdpar <- fdPar(bbasis, curv.Lfd, lambda.vect[index])
  BSmooth <- smooth.basis(t, XX, curv.Fdpar) # vyhlazeni
  gcv[index] <- mean(BSmooth$gcv) # prumer pres vsechny pozorovane krivky
}

GCV <- data.frame(
  lambda = round(log10(lambda.vect), 3),
  GCV = gcv
)

# najdeme hodnotu minima
lambda.opt <- lambda.vect[which.min(gcv)]

curv.fdPar <- fdPar(bbasis, curv.Lfd, lambda.opt)
BSmooth <- smooth.basis(t, XX, curv.fdPar)
XXfd <- BSmooth$fd

fdobjSmootheval <- eval.fd(fdobj = XXfd, evalarg = t)

# rozdeleni na testovaci a trenovaci cast
set.seed(42)
library(caTools)
split <- sample.split(XXfd$fdnames$reps, SplitRatio = 0.7)

# vytvoreni vektoru 0 a 1, 0 pro < 20 a 1 pro > 20 
Y <- ifelse(labels == 'large', 1, 0)

X.train <- subset(XXfd, split == TRUE)
X.test <- subset(XXfd, split == FALSE)

Y.train <- subset(Y, split == TRUE)
Y.test <- subset(Y, split == FALSE)

# vytvoreni vektoru 0 a 1, 0 pro < 20 a 1 pro > 20 
# Y <- ifelse(labels == 'large', 1, 0)
# X.train <- XXfd
# Y.train <- Y

# table(Y.train)
# 
# # relativni zastoupeni
# table(Y.train) / sum(table(Y.train))
```


```r
# analyza hlavnich komponent
data.PCA <- pca.fd(X.train, nharm = 10) # nharm - maximalni pocet HK
nharm <- which(cumsum(data.PCA$varprop) >= 0.9)[1] # urceni p
if(nharm == 1) nharm <- 2 # aby bylo mozne vykreslovat grafy,
# potrebujeme alespon 2 HK

data.PCA <- pca.fd(X.train, nharm = nharm) 
data.PCA.train <- as.data.frame(data.PCA$scores) # skore prvnich p HK
data.PCA.train$Y <- factor(Y.train) # prislusnost do trid
```

U všech třech jader projdeme hodnoty hyperparametru $C$ v intervalu $[10^{-3}, 10^{3}]$, přičemž u jádra polynomiálního zafixujeme hyperparametr $p$ na hodnotě 3, neboť pro jiné celočíselné hodnoty metoda nedává zdaleka tak dobré výsledky. Naopak pro radiální jádro využijeme k volbě optimální hodnoty hyperparametru $\gamma$ opět 10-násobnou CV, přičemž uvažujeme hodnoty v intervalu $[10^{-3}, 10^{2}]$. Zvolíme `coef0` $= 1$. 


```r
set.seed(42)
k_cv <- 10
# rozdelime trenovaci data na k casti
library(caret)
folds <- createMultiFolds(1:length(Y.train), k = k_cv, time = 1)

# ktere hodnoty gamma chceme uvazovat
gamma.cv <- 10^seq(-3, 2, length = 15)
C.cv <- 10^seq(-3, 3, length = 20)
p.cv <- c(2, 3, 4, 5)
coef0 <- 1

# list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro dane
# v radcich budou hodnoty pro danou gamma a vrstvy odpovidaji folds
CV.results <- list(
  SVM.l = array(NA, dim = c(length(C.cv), k_cv)),
  SVM.p = array(NA, dim = c(length(C.cv), length(p.cv), k_cv)),
  SVM.r = array(NA, dim = c(length(C.cv), length(gamma.cv), k_cv))
)

# nejprve projdeme hodnoty C
for (C in C.cv) {
  # projdeme jednotlive folds
  for (index_cv in 1:k_cv) {
    # definice testovaci a trenovaci casti pro CV
    fold <- folds[[index_cv]]
    cv_sample <- 1:dim(data.PCA.train)[1] %in% fold
    
    data.PCA.train.cv <- as.data.frame(data.PCA.train[cv_sample, ])
    data.PCA.test.cv <- as.data.frame(data.PCA.train[!cv_sample, ])
    
    ## LINEARNI JADRO
    # sestrojeni modelu
    clf.SVM.l <- svm(Y ~ ., data = data.PCA.train.cv,
                     type = 'C-classification',
                     scale = TRUE,
                     cost = C,
                     kernel = 'linear')
    
    # presnost na validacnich datech
    predictions.test.l <- predict(clf.SVM.l, newdata = data.PCA.test.cv)
    presnost.test.l <- table(data.PCA.test.cv$Y, predictions.test.l) |>
      prop.table() |> diag() |> sum()
    
    # presnosti vlozime na pozice pro dane C a fold
    CV.results$SVM.l[(1:length(C.cv))[C.cv == C], 
                     index_cv] <- presnost.test.l
    
    ## POLYNOMIALNI JADRO
    for (p in p.cv) {
      # sestrojeni modelu
      clf.SVM.p <- svm(Y ~ ., data = data.PCA.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       coef0 = coef0,
                       degree = p,
                       kernel = 'polynomial')
      
      # presnost na validacnich datech
      predictions.test.p <- predict(clf.SVM.p, newdata = data.PCA.test.cv)
      presnost.test.p <- table(data.PCA.test.cv$Y, predictions.test.p) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane C, p a fold
      CV.results$SVM.p[(1:length(C.cv))[C.cv == C], 
                       (1:length(p.cv))[p.cv == p],
                       index_cv] <- presnost.test.p
    }
        
    ## RADIALNI JADRO
    for (gamma in gamma.cv) {
      # sestrojeni modelu
      clf.SVM.r <- svm(Y ~ ., data = data.PCA.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       gamma = gamma,
                       kernel = 'radial')
      
      # presnost na validacnich datech
      predictions.test.r <- predict(clf.SVM.r, newdata = data.PCA.test.cv)
      presnost.test.r <- table(data.PCA.test.cv$Y, predictions.test.r) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane C, gamma a fold
      CV.results$SVM.r[(1:length(C.cv))[C.cv == C], 
                       (1:length(gamma.cv))[gamma.cv == gamma],
                       index_cv] <- presnost.test.r
    }
  }
}
```

Nyní zprůměrujeme výsledky 10-násobné CV tak, abychom pro jednu hodnotu hyperparametru (případně jednu kombinaci hodnot) měli jeden odhad validační chybovosti. Přitom určíme i optimální hodnoty jednotlivých hyperparametrů.


```r
# spocitame prumerne presnosti pro jednotliva C pres folds
## Linearni jadro
CV.results$SVM.l <- apply(CV.results$SVM.l, 1, mean)
## Polynomialni jadro
CV.results$SVM.p <- apply(CV.results$SVM.p, c(1, 2), mean)
## Radialni jadro
CV.results$SVM.r <- apply(CV.results$SVM.r, c(1, 2), mean)

C.opt <- c(which.max(CV.results$SVM.l), 
           which.max(CV.results$SVM.p) %% length(C.cv), 
           which.max(CV.results$SVM.r) %% length(C.cv))
C.opt[C.opt == 0] <- length(C.cv)
C.opt <- C.cv[C.opt]

gamma.opt <- which.max(t(CV.results$SVM.r)) %% length(gamma.cv)
gamma.opt[gamma.opt == 0] <- length(gamma.cv)
gamma.opt <- gamma.cv[gamma.opt]

p.opt <- which.max(t(CV.results$SVM.p)) %% length(p.cv)
p.opt[p.opt == 0] <- length(p.cv)
p.opt <- p.cv[p.opt]

presnost.opt.cv <- c(max(CV.results$SVM.l), 
                     max(CV.results$SVM.p),
                     max(CV.results$SVM.r))
```




```r
# sestrojeni modelu
clf.SVM.l.PCA <- svm(Y ~ ., data = data.PCA.train,
                     type = 'C-classification',
                     scale = TRUE,
                     cost = C.opt[1],
                     kernel = 'linear')

clf.SVM.p.PCA <- svm(Y ~ ., data = data.PCA.train,
                     type = 'C-classification',
                     scale = TRUE,
                     cost = C.opt[2],
                     degree = p.opt,
                     coef0 = coef0,
                     kernel = 'polynomial')

clf.SVM.r.PCA <- svm(Y ~ ., data = data.PCA.train,
                     type = 'C-classification',
                     scale = TRUE,
                     cost = C.opt[3],
                     gamma = gamma.opt,
                     kernel = 'radial')
```




```r
# pridame diskriminacni hranici
np <- 1001 # pocet bodu site
# x-ova osa ... 1. HK
nd.x <- seq(from = min(data.PCA.train$V1) - 5, 
            to = max(data.PCA.train$V1) + 5, length.out = np)
# y-ova osa ... 2. HK
nd.y <- seq(from = min(data.PCA.train$V2) - 5, 
            to = max(data.PCA.train$V2) + 5, length.out = np)
# pripad pro 2 HK ... p = 2
nd <- expand.grid(V1 = nd.x, V2 = nd.y)

nd <- rbind(nd, nd, nd) |> mutate(
   prd = c(as.numeric(predict(clf.SVM.l.PCA, newdata = nd, type = 'response')),
           as.numeric(predict(clf.SVM.p.PCA, newdata = nd, type = 'response')),
           as.numeric(predict(clf.SVM.r.PCA, newdata = nd, type = 'response'))),
   kernel = rep(c('Lineární', 'Polynomiální', 'Radiální'),
                each = length(as.numeric(predict(clf.SVM.l.PCA, 
                                                 newdata = nd,
                                                 type = 'response')))) |>
     as.factor())

df_SV <- rbind(data.PCA.train[clf.SVM.l.PCA$index, ] |>
                 mutate(kernel = 'Lineární'),
               data.PCA.train[clf.SVM.p.PCA$index, ] |>
                 mutate(kernel = 'Polynomiální'),
               data.PCA.train[clf.SVM.r.PCA$index, ] |> 
                 mutate(kernel = 'Radiální'))

pSVM <- ggplot(data = data.PCA.train,
       aes(x = V1, y = V2, colour = Y)) +
  geom_contour(data = nd, aes(x = V1, y = V2, z = prd), 
              colour = 'grey2', size = 0.25) + 
  labs(x = '$X_1$',
       y = '$X_2$',
       colour = 'Klasifikační\n třída',
       fill = 'none') +
  scale_colour_manual(values = c('tomato', 'deepskyblue2')) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title.align = 0.5,
        legend.position = 'none') + 
  scale_y_continuous(expand = c(-0.02, -0.02),
                     limits = c(-3.5, 2.5)) + 
  scale_x_continuous(expand = c(-0.02, -0.02),
                     limits = c(-15, 26)) + 
  facet_wrap(~kernel) + 
  geom_contour_filled(data = nd, aes(x = V1, y = V2, z = prd, colour = prd),
                      breaks = c(1, 2, 3), alpha = 0.1,
                      show.legend = F) +
  scale_fill_manual(values = c('tomato', 'deepskyblue2')) +
  geom_point(data = df_SV, col = 'grey2', alpha = 0.7,
                      size = 1.5) +
  geom_point(size = 1.2) + 
  geom_point(data = df_SV, col = 'grey2', alpha = 0.4,
                      size = 1.5)

pSVM
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-51-1.png" alt="Skóre prvních dvou hlavních komponent, barevně odlišené podle příslušnosti do klasifikační třídy. Černě je vyznačena dělící hranice (přímka, resp. křivky v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí metody SVM." width="672" />
<p class="caption">(\#fig:unnamed-chunk-51)Skóre prvních dvou hlavních komponent, barevně odlišené podle příslušnosti do klasifikační třídy. Černě je vyznačena dělící hranice (přímka, resp. křivky v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí metody SVM.</p>
</div>

```r
# ggsave("figures/kap4_SVM.tex", device = tikz, width = 8, height = 3)
```

## Materiály pro Kapitolu 5

V této sekci uvedeme podpůrné grafy pro pátou kapitolu diplomové práce.

### Diskretizace intervalu

Chtěli bychom se podívat na hodnoty skalárních součinů funkcí, které jsou blízko u sebe a naopak které se tvarem velmi liší.

#### `tecator` data

Podívejme se také na data `tecator`.



```r
data <- ddalpha::dataf.tecator()

data.gr <- data$dataf[[1]]$vals
for(i in 2:length(data$labels)) {
  data.gr <- rbind(data.gr, data$dataf[[i]]$vals)
  }
data.gr <- cbind(data.frame(wave = data$dataf[[1]]$args),
                 t(data.gr))

# vektor trid
labels <- data$labels |> unlist()
# prejmenovani podle tridy
colnames(data.gr) <- c('wavelength',
                       paste0(labels, 1:length(data$labels)))
```


```r
library(fda.usc)
t <- data.gr$wavelength
rangeval <- range(t)
breaks <- t
norder <- 6

bbasis <- create.bspline.basis(rangeval = rangeval, 
                               norder = norder, 
                               breaks = breaks)

curv.Lfd <- int2Lfd(4) # penalizujeme 4. derivaci

# spojeni pozorovani do jedne matice
XX <- data.gr[, -1] |> as.matrix()

lambda.vect <- 10^seq(from = -2, to = 1, length.out = 50) # vektor lambd
gcv <- rep(NA, length = length(lambda.vect)) # prazdny vektor pro ulozebi GCV

for(index in 1:length(lambda.vect)) {
  curv.Fdpar <- fdPar(bbasis, curv.Lfd, lambda.vect[index])
  BSmooth <- smooth.basis(t, XX, curv.Fdpar) # vyhlazeni
  gcv[index] <- mean(BSmooth$gcv) # prumer pres vsechny pozorovane krivky
}

GCV <- data.frame(
  lambda = round(log10(lambda.vect), 3),
  GCV = gcv
)

# najdeme hodnotu minima
lambda.opt <- lambda.vect[which.min(gcv)]

curv.fdPar <- fdPar(bbasis, curv.Lfd, lambda.opt)
BSmooth <- smooth.basis(t, XX, curv.fdPar)
XXfd <- BSmooth$fd # * as.numeric(1 / norm.fd(BSmooth$fd[1]))
# set norm equal to one
norms <- c()
for (i in 1:215) {norms <- c(norms, as.numeric(1 / norm.fd(BSmooth$fd[i])))}
XXfd_norm <- XXfd 
XXfd_norm$coefs <- XXfd_norm$coefs * matrix(norms, ncol = 215, nrow = 104, byrow = T)

fdobjSmootheval <- eval.fd(fdobj = XXfd_norm, evalarg = t)

# rozdeleni na testovaci a trenovaci cast
set.seed(42)
library(caTools)
split <- sample.split(XXfd$fdnames$reps, SplitRatio = 0.7)

# vytvoreni vektoru 0 a 1, 0 pro < 20 a 1 pro > 20 
Y <- ifelse(labels == 'large', 1, 0)

X.train <- subset(XXfd, split == TRUE)
X.test <- subset(XXfd, split == FALSE)

Y.train <- subset(Y, split == TRUE)
Y.test <- subset(Y, split == FALSE)

# vytvoreni vektoru 0 a 1, 0 pro < 20 a 1 pro > 20 
# Y <- ifelse(labels == 'large', 1, 0)
# X.train <- XXfd
# Y.train <- Y
```

Spočítáme skalární součiny prvního s ostatními.






```r
n <- dim(XX)[2]
abs.labs <- c("$< 20 \\%$", "$> 20 \\%$")
# abs.labs <- c("$Y = {-1}$", "$Y = 1$")
names(abs.labs) <- c('small', 'large')

DFsmooth <- data.frame(
  t = rep(t, n),
  time = factor(rep(1:n, each = length(t))),
  Smooth = c(fdobjSmootheval),
  Fat = factor(rep(labels, each = length(t)), levels = c('small', 'large'))
)

DFmean <- data.frame(
  t = rep(t, 2),
  Mean = c(apply(fdobjSmootheval[ , labels == 'small'], 1, mean), 
            apply(fdobjSmootheval[ , labels == 'large'], 1, mean)),
  Fat = factor(rep(c('small', 'large'), each = length(t)),
                 levels = c('small', 'large'))
)

DFsmooth |> filter(time %in% as.character(c(nn))) |>
  ggplot(aes(x = t, y = Smooth, color = Fat)) +
  geom_line(linewidth = 1.1, aes(group = time, linetype = 'apozor x1')) +
  geom_line(data = DFsmooth |> 
              filter(time %in% as.character(order(Inprod_vect)[1:4])),
            aes(group = time, linetype = 'nejmensi'), linewidth = 0.6) +
  geom_line(data = DFsmooth |> 
              filter(time %in% as.character(rev(order(Inprod_vect))[1:5])),
            aes(group = time, linetype = 'nejvetsi'), linewidth = 0.6) +
  geom_line(linewidth = 1.1, aes(group = time, linetype = 'apozor x1')) +
  theme_bw() +
  # facet_wrap(~Fat,
  #            labeller = labeller(Fat = abs.labs)) + 
  labs(x = "Vlnová délka [v nm]",
       y = "Absorbance",
       colour = 'Obsah tuku',#"Klasifikační\n      třída",
       linetype = 'Styl čáry') + 
  # scale_color_discrete(guide="none") +
  scale_colour_manual(values = c('tomato', 'deepskyblue2'), labels = abs.labs) + 
  # guides(color = guide_legend(position = 'none')) +
  # scale_color_discrete(labels = abs.labs) +
  scale_linetype_manual(values=c('solid', "dotted", "longdash")) +
  theme(legend.position = c(0.15, 0.75),
        #legend.box="vertical",
        panel.grid = element_blank())
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-55-1.png" alt="Vykreslení čtyř pozorování s minimální a maximální hodnotou skalárního součinu." width="672" />
<p class="caption">(\#fig:unnamed-chunk-55)Vykreslení čtyř pozorování s minimální a maximální hodnotou skalárního součinu.</p>
</div>

```r
# ggsave("figures/kap5_discretization_tecator.tex", device = tikz, width = 4.5, height = 4.5)
```


#### `phoneme` data

Použijeme data `phoneme`.


```r
library(fda.usc)
# nacteni dat
data <- read.delim2('phoneme.txt', header = T, sep = ',')

# zmenime dve promenne na typ factor
data <- data |> 
  mutate(g = factor(g),
         speaker = factor(speaker))

# numericke promenne prevedeme opravdu na numericke
data[, 2:257] <- as.numeric(data[, 2:257] |> as.matrix())

tr_vs_test <- str_split(data$speaker, '\\.') |> unlist()
tr_vs_test <- tr_vs_test[seq(1, length(tr_vs_test), by = 4)]
data$train <- ifelse(tr_vs_test == 'train', TRUE, FALSE)

# vybrane fonemy ke klasifikaci
phoneme_subset <- c('aa', 'ao')

# testovaci a trenovaci data
data_train <- data |> filter(train) |> filter(g %in% phoneme_subset)
data_test <- data |> filter(!train) |> filter(g %in% phoneme_subset)

# odstranime sloupce, ktere nenesou informaci o frekvenci a 
# transponujeme tak, aby ve sloupcich byly jednotlive zaznamy
X_train <- data_train[, -c(1, 258, 259, 260)] |> t()
X_test <- data_test[, -c(1, 258, 259, 260)] |> t()

# prejmenujeme radky a sloupce
rownames(X_train) <- 1:256
colnames(X_train) <- paste0('train', data_train$row.names)
rownames(X_test) <- 1:256
colnames(X_test) <- paste0('test', data_test$row.names)

# definujeme vektor fonemu
y_train <- data_train[, 258] |> factor(levels = phoneme_subset)
y_test <- data_test[, 258] |> factor(levels = phoneme_subset)
y <- c(y_train, y_test)
```



```r
t <- 1:256
rangeval <- range(t)
breaks <- t
norder <- 4

bbasis <- create.bspline.basis(rangeval = rangeval, 
                               norder = norder, 
                               breaks = breaks)

curv.Lfd <- int2Lfd(2) # penalizujeme 2. derivaci

# spojeni pozorovani do jedne matice
XX <- cbind(X_train, X_test) |> as.matrix()

lambda.vect <- 10^seq(from = 1, to = 3, length.out = 35) # vektor lambd
gcv <- rep(NA, length = length(lambda.vect)) # prazdny vektor pro ulozebi GCV

for(index in 1:length(lambda.vect)) {
  curv.Fdpar <- fdPar(bbasis, curv.Lfd, lambda.vect[index])
  BSmooth <- smooth.basis(t, XX, curv.Fdpar) # vyhlazeni
  gcv[index] <- mean(BSmooth$gcv) # prumer pres vsechny pozorovane krivky
}

GCV <- data.frame(
  lambda = round(log10(lambda.vect), 3),
  GCV = gcv
)

# najdeme hodnotu minima
lambda.opt <- lambda.vect[which.min(gcv)]

curv.fdPar <- fdPar(bbasis, curv.Lfd, lambda.opt)
BSmooth <- smooth.basis(t, XX, curv.fdPar)
XXfd <- BSmooth$fd

# set norm equal to one
norms <- c()
for (i in 1:dim(XXfd$coefs)[2]) {
  norms <- c(norms, as.numeric(1 / norm.fd(BSmooth$fd[i])))
  }
XXfd_norm <- XXfd 
XXfd_norm$coefs <- XXfd_norm$coefs * matrix(norms, 
                                            ncol = dim(XXfd$coefs)[2],
                                            nrow = dim(XXfd$coefs)[1],
                                            byrow = T)

fdobjSmootheval <- eval.fd(fdobj = XXfd_norm, evalarg = t)
```

Spočítáme skalární součiny prvního log-periodogramu s ostatními.






```r
n <- dim(XX)[2]
y <- c(y_train, y_test)
DFsmooth <- data.frame(
  t = rep(t, n),
  time = factor(rep(1:n, each = length(t))),
  Smooth = c(fdobjSmootheval),
  Phoneme = rep(y, each = length(t)))

DFmean <- data.frame(
  t = rep(t, 2),
  Mean = c(apply(fdobjSmootheval[ , y == phoneme_subset[1]], 1, mean),
            apply(fdobjSmootheval[ , y == phoneme_subset[2]], 1, mean)),
  Phoneme = factor(rep(phoneme_subset, each = length(t)),
                 levels = levels(y))
)

DFsmooth |> 
  filter(time %in% as.character(nn)) |>
  ggplot(aes(x = t, y = Smooth, color = Phoneme)) + 
  geom_line(linewidth = 1.1, aes(group = time, linetype = 'apozor x1')) +
  geom_line(data = DFsmooth |> 
              filter(time %in% as.character(order(Inprod_vect)[1:4])),
            aes(group = time, linetype = 'nejmensi'), linewidth = 0.6) +
  geom_line(data = DFsmooth |> 
              filter(time %in% as.character(rev(order(Inprod_vect))[1:5])),
            aes(group = time, linetype = 'nejvetsi'), linewidth = 0.6) +
  geom_line(linewidth = 1.1, aes(group = time, linetype = 'apozor x1')) +
  theme_bw() +
  labs(x = 'Frekvence',
       y = 'Log-periodogram',
       colour = 'Foném',
       linetype = 'Styl čáry') +
  # scale_colour_discrete(labels = phoneme_subset) +
  scale_colour_manual(values = c('tomato', 'deepskyblue2'),
                      labels = phoneme_subset) +
  scale_linetype_manual(values=c('solid', "dotted", "longdash")) +
  theme(legend.position = c(0.8, 0.75),
        panel.grid = element_blank())
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-59-1.png" alt="Vykreslení čtyř pozorování s minimální a maximální hodnotou skalárního součinu." width="672" />
<p class="caption">(\#fig:unnamed-chunk-59)Vykreslení čtyř pozorování s minimální a maximální hodnotou skalárního součinu.</p>
</div>

```r
# ggsave("figures/kap5_discretization_phoneme.tex", device = tikz, width = 4.5, height = 4.5)
```

### Support vector regression (SVR)

Ukázka metody SVR na obou datových souborech.

#### `tecator` data


```r
t <- data.gr$wavelength
rangeval <- range(t)
breaks <- t
norder <- 6

bbasis <- create.bspline.basis(rangeval = rangeval, 
                               norder = norder, 
                               breaks = breaks)

curv.Lfd <- int2Lfd(4) # penalizujeme 4. derivaci

# spojeni pozorovani do jedne matice
XX <- data.gr[, -1] |> as.matrix()

lambda.vect <- 10^seq(from = -2, to = 1, length.out = 50) # vektor lambd
gcv <- rep(NA, length = length(lambda.vect)) # prazdny vektor pro ulozebi GCV

for(index in 1:length(lambda.vect)) {
  curv.Fdpar <- fdPar(bbasis, curv.Lfd, lambda.vect[index])
  BSmooth <- smooth.basis(t, XX, curv.Fdpar) # vyhlazeni
  gcv[index] <- mean(BSmooth$gcv) # prumer pres vsechny pozorovane krivky
}

GCV <- data.frame(
  lambda = round(log10(lambda.vect), 3),
  GCV = gcv
)

# najdeme hodnotu minima
lambda.opt <- lambda.vect[which.min(gcv)]

curv.fdPar <- fdPar(bbasis, curv.Lfd, lambda.opt)
BSmooth <- smooth.basis(t, XX, curv.fdPar)
XXfd <- BSmooth$fd 

fdobjSmootheval <- eval.fd(fdobj = XXfd, evalarg = t)
```



```r
library(e1071)
library(caret)

df_plot <- data.frame()

# model
for(i in 1:5) {
  df.svm <- data.frame(x = t,
                       y = fdobjSmootheval[, i])
  svm.RKHS <- svm(y ~ x, data = df.svm, 
                  kernel = 'radial',
                  type = 'eps-regression',
                  epsilon = 0.03,
                  gamma = 0.5,
                  cost = 1, 
                  tolerance = 0.001,
                  shrinking = TRUE,
                  scale = TRUE)
  
  svm.RKHS <- train(y ~ x, data = df.svm,
                    method = 'svmRadial',
                    metric = "RMSE",
                    preProcess = c('center', 'scale'),
                    trControl = trainControl(
                      method = "repeatedcv",
                      number = 10,
                      repeats = 10,
                      verboseIter = FALSE
                    )
                    # trControl = trainControl(method = "none"),
                    # Telling caret not to re-tune
                    # tuneGrid = data.frame(sigma = 1000, C = 1000)
                    # Specifying the parameters
                    )
  df_plot <- rbind(
    df_plot, 
    data.frame(
      x = t,
      y = svm.RKHS$finalModel@fitted * 
             svm.RKHS$finalModel@scaling$y.scale$`scaled:scale` +
             svm.RKHS$finalModel@scaling$y.scale$`scaled:center`,
      line = 'estimate',
      curve = as.character(i)) |>
      rbind(data.frame(
        x = t,
        y = fdobjSmootheval[, i],
        line = 'sample',
        curve = as.character(i)
  )))
}
```


Vykresleme si pro lepší představu odhad křivky (červeně) společně s pozorovanou křivkou (modře).


```r
df_plot |> filter(curve %in% c('1', '2', '3')) |> 
  ggplot(aes(x, y, col = line, linetype = curve)) +
  geom_line(linewidth = 0.8) + 
  theme_bw() +
  labs(x = "Vlnová délka [v nm]",
       y = "Absorbance",
       colour = 'Křivka') +
  # scale_colour_discrete(labels = phoneme_subset) +
  scale_colour_manual(values = c('tomato', 'deepskyblue2'),
                      labels = c("odhadnutá", "pozorovaná")) +
  scale_linetype_manual(values = c('solid', "dotted", "dashed")) +
  theme(legend.position = c(0.17, 0.85),
        panel.grid = element_blank()) + 
  guides(linetype = 'none')
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-62-1.png" alt="Porovnání pozorované a odhadnuté křivky." width="672" />
<p class="caption">(\#fig:unnamed-chunk-62)Porovnání pozorované a odhadnuté křivky.</p>
</div>

```r
# ggsave("figures/kap5_SVR_tecator.tex", device = tikz, width = 4.5, height = 4.5)
```

#### `phoneme` data


```r
data <- read.delim2('phoneme.txt', header = T, sep = ',')

# zmenime dve promenne na typ factor
data <- data |> 
  mutate(g = factor(g),
         speaker = factor(speaker))

# numericke promenne prevedeme opravdu na numericke
data[, 2:257] <- as.numeric(data[, 2:257] |> as.matrix())

tr_vs_test <- str_split(data$speaker, '\\.') |> unlist()
tr_vs_test <- tr_vs_test[seq(1, length(tr_vs_test), by = 4)]
data$train <- ifelse(tr_vs_test == 'train', TRUE, FALSE)

# vybrane fonemy ke klasifikaci
phoneme_subset <- c('aa', 'ao')

# testovaci a trenovaci data
data_train <- data |> filter(train) |> filter(g %in% phoneme_subset)
data_test <- data |> filter(!train) |> filter(g %in% phoneme_subset)

# odstranime sloupce, ktere nenesou informaci o frekvenci a 
# transponujeme tak, aby ve sloupcich byly jednotlive zaznamy
X_train <- data_train[, -c(1, 258, 259, 260)] |> t()
X_test <- data_test[, -c(1, 258, 259, 260)] |> t()

# prejmenujeme radky a sloupce
rownames(X_train) <- 1:256
colnames(X_train) <- paste0('train', data_train$row.names)
rownames(X_test) <- 1:256
colnames(X_test) <- paste0('test', data_test$row.names)

# definujeme vektor fonemu
y_train <- data_train[, 258] |> factor(levels = phoneme_subset)
y_test <- data_test[, 258] |> factor(levels = phoneme_subset)
y <- c(y_train, y_test)
```


```r
t <- 1:256
rangeval <- range(t)
breaks <- t
norder <- 4

bbasis <- create.bspline.basis(rangeval = rangeval, 
                               norder = norder, 
                               breaks = breaks)

curv.Lfd <- int2Lfd(2) # penalizujeme 2. derivaci

# spojeni pozorovani do jedne matice
XX <- cbind(X_train, X_test) |> as.matrix()

lambda.vect <- 10^seq(from = 1, to = 3, length.out = 35) # vektor lambd
gcv <- rep(NA, length = length(lambda.vect)) # prazdny vektor pro ulozebi GCV

for(index in 1:length(lambda.vect)) {
  curv.Fdpar <- fdPar(bbasis, curv.Lfd, lambda.vect[index])
  BSmooth <- smooth.basis(t, XX, curv.Fdpar) # vyhlazeni
  gcv[index] <- mean(BSmooth$gcv) # prumer pres vsechny pozorovane krivky
}

GCV <- data.frame(
  lambda = round(log10(lambda.vect), 3),
  GCV = gcv
)

# najdeme hodnotu minima
lambda.opt <- lambda.vect[which.min(gcv)]

curv.fdPar <- fdPar(bbasis, curv.Lfd, lambda.opt)
BSmooth <- smooth.basis(t, XX, curv.fdPar)
XXfd <- BSmooth$fd

fdobjSmootheval <- eval.fd(fdobj = XXfd, evalarg = t)
```




```r
df_plot <- data.frame()
# model
for(i in 1:5) {
  df.svm <- data.frame(x = t,
                       y = fdobjSmootheval[, i])
  svm.RKHS <- svm(y ~ x, data = df.svm, 
                  kernel = 'radial',
                  type = 'eps-regression',
                  epsilon = 0.03,
                  gamma = 0.5,
                  cost = 1, 
                  tolerance = 0.001,
                  shrinking = TRUE,
                  scale = TRUE)
  
  svm.RKHS <- train(y ~ x, data = df.svm,
                    method = 'svmRadial',
                    metric = "RMSE",
                    preProcess = c('center', 'scale'),
                    trControl = trainControl(
                      method = "repeatedcv",
                      number = 10,
                      repeats = 10,
                      verboseIter = FALSE
                    )
                    # trControl = trainControl(method = "none"),
                    # Telling caret not to re-tune
                    # tuneGrid = data.frame(sigma = 1000, C = 1000)
                    # Specifying the parameters
                    )
  df_plot <- rbind(
    df_plot, 
    data.frame(
      x = t,
      y = svm.RKHS$finalModel@fitted * 
             svm.RKHS$finalModel@scaling$y.scale$`scaled:scale` +
             svm.RKHS$finalModel@scaling$y.scale$`scaled:center`,
      line = 'estimate',
      curve = as.character(i)) |>
      rbind(data.frame(
        x = t,
        y = fdobjSmootheval[, i],
        line = 'sample',
        curve = as.character(i)
  )))
}
```


Vykresleme si pro lepší představu odhad křivky (červeně) společně s pozorovanou křivkou (modře).


```r
df_plot |> filter(curve %in% c('1', '5', '3')) |> 
  ggplot(aes(x, y, col = line, linetype = curve)) +
  geom_line(linewidth = 0.8) + 
  theme_bw() +
  labs(x = 'Frekvence',
       y = 'Log-periodogram',
       colour = 'Křivka') +
  # scale_colour_discrete(labels = phoneme_subset) +
  scale_colour_manual(values = c('tomato', 'deepskyblue2'),
                      labels = c("odhadnutá", "pozorovaná")) +
  scale_linetype_manual(values = c('solid', "dotted", "dashed")) +
  theme(legend.position = c(0.8, 0.85),
        panel.grid = element_blank()) + 
  guides(linetype = 'none')
```

<div class="figure">
<img src="14-Application_5_files/figure-html/unnamed-chunk-66-1.png" alt="Porovnání pozorované a odhadnuté křivky." width="672" />
<p class="caption">(\#fig:unnamed-chunk-66)Porovnání pozorované a odhadnuté křivky.</p>
</div>

```r
# ggsave("figures/kap5_SVR_phoneme.tex", device = tikz, width = 4.5, height = 4.5)
```

## Materiály pro Kapitolu 6

Veškeré grafické podklady i číselné výstupy prezentované v Diplomové práci v Kapitole 6 jsou k dispozici v Kapitole \@ref(simulace3) (případně také v Kapitolách \@ref(simulace3sigma), \@ref(simulace3shift) a \@ref(simulace3diskr)) a v Kapitole \@ref(simulace4). Krabicové diagramy testovacích chybovostí jsme vizuálně upravili pro potřeby diplomové práce (barevnost, změna měřítka, popisky), kód použitý k jejich vygenerování je k vidění v příslušných buňkách (zapoznámkovaný).


## Materiály pro Kapitolu 7

Veškeré grafické podklady i číselné výstupy prezentované v Diplomové práci v Kapitole 7 jsou k dispozici v Kapitole \@ref(aplikace2), která je věnována datovému souboru `phoneme`, a Kapitole \@ref(aplikace3) věnované datovému souboru `tecator`. Krabicové diagramy testovacích chybovostí jsme vizuálně upravili pro potřeby diplomové práce (barevnost, změna měřítka, popisky), kód použitý k jejich vygenerování je k vidění v příslušných buňkách (zapoznámkovaný).

