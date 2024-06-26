# Aplikace na reálných datech 3

Stejně jako v předchozí kapitole se budeme v této části dokumentu zabývat aplikací dříve popsaných metod (pro více podrobností viz například sekci \@ref(simulace1)) na reálná data `tecator`, která jsou dostupná například v balíčku `ddalpha`. Podrobný popis dat pak můžeme nalézt [zde](https://search.r-project.org/CRAN/refmans/ddalpha/html/dataf.tecator.html). Jedná se o datový soubor obsahující spektrometrické křivky (absorbanční křivky měřené ve 100 vlnových délkách). 

Pro každý kus jemně nasekaného masa pozorujeme jednu spektrometrickou křivku, která odpovídá absorbanci naměřené při 100 vlnových délkách. Kusy jsou rozděleny podle [Ferratyho a Vieu (2006)](https://link.springer.com/book/10.1007/0-387-36620-2) do dvou tříd: s malým ($< 20\,\%$) a velkým ($\geq 20\,\%$) obsahem tuku získaným analytickým chemickým zpracováním. Naším cílem bude klasifikovat spektrometrické křivky na intervalu $I = [850 \text{ nm}, 1050 \text{ nm}]$ na základě obsahu tuku. Jak uvidíme z výsledků v části \@ref(klasA3deriv), je výhodné uvažovat druhou derivaci křivek. 

Začněme nejprve s načtením a vykreslením dat. Data jsou uložena poněkud složitě, proto pro lepší budou práci s nimi si je uložíme do praktičtějšího formátu. Pojmenujeme si také příslušné sloupce podle toho, zda obsah tuku je malý (`small`) nebo velký (`large`).


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

Spektrometrické křivky si vykresleme podle skupiny.


```r
abs.labs <- c("Obsah tuku < 20 %", "Obsah tuku > 20 %")
names(abs.labs) <- c('small', 'large')

pivot_longer(data.gr, cols = large1:large215, names_to = 'sample',
                        values_to = 'absorbance', cols_vary = 'slowest') |>
  mutate(sample = as.factor(sample),
         Abs = factor(rep(labels, each = length(data.gr$wavelength)), 
                      levels = c('small', 'large'))) |>
  ggplot(aes(x = wavelength, y = absorbance, colour = Abs, group = sample)) + 
  geom_line(linewidth = 0.5) + 
  theme_bw() +
  facet_wrap(~Abs,
             labeller = labeller(Abs = abs.labs)) + 
  labs(x = "Vlnová délka [v nm]",
       y = "Absorbance",
       colour = "Obsah tuku") + 
  theme(legend.position = 'none') +
  scale_color_discrete(labels = abs.labs)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-2-1.png" alt="Absorpční křivky podle skupiny." width="672" />
<p class="caption">(\#fig:unnamed-chunk-2)Absorpční křivky podle skupiny.</p>
</div>

## Vyhlazení pozorovaných křivek

Nyní převedeme pozorované diskrétní hodnoty (vektory hodnot) na funkcionální objekty, se kterými budeme následně pracovat.
Jelikož se nejedná o periodické křivky na intervalu $I = [850, 1050]$, využijeme k vyhlazení B-sline bázi.

Za uzly bereme celý vektor `wavelength`, standardně bychom uvažovali kubické spliny, protože ale budeme chtít pracovat s druhou derivací, volíme `norder = 6`. Ze stejného důvodu budeme penalizovat čtvrtou derivaci funkcí.


```r
t <- data.gr$wavelength
rangeval <- range(t)
breaks <- t
norder <- 6

bbasis <- create.bspline.basis(rangeval = rangeval, 
                               norder = norder, 
                               breaks = breaks)

curv.Lfd <- int2Lfd(4) # penalizujeme 4. derivaci
```

Najdeme vhodnou hodnotu vyhlazovacího parametru $\lambda > 0$ pomocí $GCV(\lambda)$, tedy pomocí zobecněné cross--validace.
Hodnotu $\lambda$ budeme uvažovat pro obě pohlaví stejnou, neboť pro testovací pozorování bychom dopředu nevěděli, kterou hodnotu $\lambda$ máme v případě rozdílné volby pro každou třídu volit.


```r
# spojeni pozorovani do jedne matice
XX <- data.gr[, -1] |> as.matrix()

lambda.vect <- 10^seq(from = -1, to = 0.5, length.out = 50) # vektor lambd
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

Pro lepší znázornění si vykreslíme průběh $GCV(\lambda)$.


```r
GCV |> ggplot(aes(x = lambda, y = GCV)) + 
  geom_line(linetype = 'solid', linewidth = 0.6) + 
  geom_point(size = 1.7) + 
  theme_bw() + 
  labs(x = bquote(paste(log[10](lambda), ' ;   ', 
                        lambda[optimal] == .(round(lambda.opt, 4)))),
       y = expression(GCV(lambda))) + 
  geom_point(aes(x = log10(lambda.opt), y = min(gcv)), colour = 'red', size = 3)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-5-1.png" alt="Průběh $GCV(\lambda)$ pro zvolený vektor $\boldsymbol\lambda$. Na ose $x$ jsou hodnoty vyneseny v logaritmické škále při základu 10. Červeně je znázorněna optimální hodnota vyhlazovacího parametru $\lambda_{optimal}$." width="672" />
<p class="caption">(\#fig:unnamed-chunk-5)Průběh $GCV(\lambda)$ pro zvolený vektor $\boldsymbol\lambda$. Na ose $x$ jsou hodnoty vyneseny v logaritmické škále při základu 10. Červeně je znázorněna optimální hodnota vyhlazovacího parametru $\lambda_{optimal}$.</p>
</div>

S touto optimální volbou vyhlazovacího parametru $\lambda$ nyní vyhladíme všechny funkce.


```r
curv.fdPar <- fdPar(bbasis, curv.Lfd, lambda.opt)
BSmooth <- smooth.basis(t, XX, curv.fdPar)
XXfd <- BSmooth$fd

fdobjSmootheval <- eval.fd(fdobj = XXfd, evalarg = t)
```




Ještě znázorněme všechny křivky včetně průměru zvlášť pro každou třídu.


```r
library(tikzDevice)
n <- dim(XX)[2]

DFsmooth <- data.frame(
  t = rep(t, n),
  time = factor(rep(1:n, each = length(t))),
  Smooth = c(fdobjSmootheval),
  Fat = factor(rep(labels, each = length(t)), levels = c('small', 'large'))
)

DFmean <- data.frame(
  t = rep(t, 2),
  Mean = c(eval.fd(fdobj = mean.fd(XXfd[labels == 'small']), evalarg = t),
           eval.fd(fdobj = mean.fd(XXfd[labels == 'large']), evalarg = t)),
  # c(apply(fdobjSmootheval[ , labels == 'small'], 1, mean), 
  #           apply(fdobjSmootheval[ , labels == 'large'], 1, mean)),
  Fat = factor(rep(c('small', 'large'), each = length(t)),
                 levels = c('small', 'large'))
)

DFsmooth |> ggplot(aes(x = t, y = Smooth, color = Fat)) + 
  geom_line(linewidth = 0.05, aes(group = time), alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Fat,
             labeller = labeller(Fat = abs.labs)
             ) + 
  labs(x = "Vlnová délka",
       y = "Absorbance",
       colour = "Obsah tuku") + 
  theme(legend.position = 'none') +
  geom_line(data = DFmean, aes(x = t, y = Mean), 
            colour = 'grey2', linewidth = 1.25, linetype = 'solid') +
  scale_colour_manual(values = c('tomato', 'deepskyblue2'))
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-7-1.png" alt="Vykreslení všech vyhlazených pozorovaných křivek, barevně jsou odlišeny křivky podle pohlaví. Černou čerchovanou čarou je zakreslen průměr pro každou třídu." width="672" />
<p class="caption">(\#fig:unnamed-chunk-7)Vykreslení všech vyhlazených pozorovaných křivek, barevně jsou odlišeny křivky podle pohlaví. Černou čerchovanou čarou je zakreslen průměr pro každou třídu.</p>
</div>

```r
# ggsave("figures/kap7_tecator_curves_mean.tex", device = tikz, width = 9, height = 4.5)
```

Vidíme, že křivky pro obě skupiny (podle obsahu tuku) jsou poměrně podobné, černou čerchovanou čarou je znázorněn průměr. Křivky se liší zejména uprostřed intervalu, kde u tučnějších vzorků nastává o jeden lokální extrém více, naopak u méně tučných vzorků vypadají křivky jednodušeji pouze s jedním globálním extrémem.

## Výpočet derivací

Jak jsme již zmínili výše, bude výhodné klasifikovat křivky na základě druhé derivace. K výpočtu derivace pro funkcionální objekt využijeme v `R` funkci `deriv.fd()` z balíčku `fda`. Jelikož chceme klasifikovat na základě druhé derivace, volíme argument `Lfdobj = 2`. Využití těchto dat bude ukázáno v sekci \@ref(klasA3deriv).


```r
XXder <- deriv.fd(XXfd, 2)
ttt <- seq(min(t), max(t), length = 501)
fdobjSmootheval_der2 <- eval.fd(fdobj = XXder, 
                                evalarg = ttt)
```

Ještě znázorněme všechny křivky včetně průměru zvlášť pro každou třídu.


```r
DFsmooth <- data.frame(
  t = rep(ttt, n),
  time = factor(rep(1:n, each = length(ttt))),
  Smooth = c(fdobjSmootheval_der2),
  Fat = factor(rep(labels, each = length(ttt)), levels = c('small', 'large'))
)

DFmean <- data.frame(
  t = rep(ttt, 2),
  Mean = c(eval.fd(fdobj = mean.fd(XXder[labels == 'small']), evalarg = ttt),
           eval.fd(fdobj = mean.fd(XXder[labels == 'large']), evalarg = ttt)),
  Fat = factor(rep(c('small', 'large'), each = length(ttt)),
                 levels = c('small', 'large'))
)

DFsmooth |> ggplot(aes(x = t, y = Smooth, color = Fat)) + 
  geom_line(linewidth = 0.05, aes(group = time), alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Fat#,
             #labeller = labeller(Fat = abs.labs)
             ) + 
  labs(x = "Vlnová délka",
       y = "Absorbance",
       colour = "Obsah tuku") + 
  theme(legend.position = 'none') +
  geom_line(data = DFmean, aes(x = t, y = Mean), 
            colour = 'grey2', linewidth = 1.25, linetype = 'solid') +
  scale_colour_manual(values = c('tomato', 'deepskyblue2'))
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-9-1.png" alt="Vykreslení všech vyhlazených pozorovaných křivek, barevně jsou odlišeny křivky podle příslušnosti do klasifikační třídy. Černou čerchovanou čarou je zakreslen průměr pro každou třídu." width="672" />
<p class="caption">(\#fig:unnamed-chunk-9)Vykreslení všech vyhlazených pozorovaných křivek, barevně jsou odlišeny křivky podle příslušnosti do klasifikační třídy. Černou čerchovanou čarou je zakreslen průměr pro každou třídu.</p>
</div>

```r
# ggsave("figures/kap7_tecator_curves_derivatives.tex", device = tikz, width = 9, height = 4.5)
```

Vidíme z obrázku výše, že nyní se průměrné křivky mezi oběma skupinami vzorků liší mnohem výrazněji než v případě původních nederivovaných křivek.

## Klasifikace křivek

Nejprve načteme potřebné knihovny pro klasifikaci.


```r
library(caTools) # pro rozdeleni na testovaci a trenovaci
library(caret) # pro k-fold CV
library(fda.usc) # pro KNN, fLR
library(MASS) # pro LDA
library(fdapace)
library(pracma)
library(refund) # pro LR na skorech
library(nnet) # pro LR na skorech
library(caret)
library(rpart) # stromy
library(rattle) # grafika
library(e1071)
library(randomForest) # nahodny les

set.seed(42)
```

Rozdělíme data v poměru 30:70 na testovací a trénovací část, abychom mohli stanovit úspěšnost klasifikace jednotlivých metod. Trénovací část použijeme při konstrukci klasifikátoru a testovací na výpočet chyby klasifikace a případně dalších charakteristik našeho modelu. Výsledné klasifikátory podle těchto spočtených charakteristik můžeme následně porovnat mezi sebou z pohledu jejich úspěšnosti klasifikace.


```r
# rozdeleni na testovaci a trenovaci cast
set.seed(42)
split <- sample.split(XXfd$fdnames$reps, SplitRatio = 0.7)

# vytvoreni vektoru 0 a 1, 0 pro < 20 a 1 pro > 20 
Y <- ifelse(labels == 'large', 1, 0)

X.train <- subset(XXfd, split == TRUE)
X.test <- subset(XXfd, split == FALSE)

Y.train <- subset(Y, split == TRUE)
Y.test <- subset(Y, split == FALSE)
```

Ještě se podíváme na zastoupení jednotlivých skupin v testovací a trénovací části dat.


```r
# absolutni zastoupeni
table(Y.train)
```

```
## Y.train
##  0  1 
## 91 59
```

```r
table(Y.test)
```

```
## Y.test
##  0  1 
## 47 18
```

```r
# relativni zastoupeni
table(Y.train) / sum(table(Y.train))
```

```
## Y.train
##         0         1 
## 0.6066667 0.3933333
```

```r
table(Y.test) / sum(table(Y.test))
```

```
## Y.test
##         0         1 
## 0.7230769 0.2769231
```

### $K$ nejbližších sousedů

Začněme neparametrickou klasifikační metodou, a to metodou $K$ nejbližších sousedů.
Nejprve si vytvoříme potřebné objekty tak, abychom s nimi mohli pomocí funkce `classif.knn()` z knihovny `fda.usc` dále pracovat.


```r
x.train <- fdata(X.train)
y.train <- as.numeric(factor(Y.train))
```

Nyní můžeme definovat model a podívat se na jeho úspěšnost klasifikace.
Poslední otázkou však zůstává, jak volit optimální počet sousedů $K$.
Mohli bychom tento počet volit jako takové $K$, při kterém nastává minimální chybovost na trénovacích datech.
To by ale mohlo vést k přeučení modelu, proto využijeme cross-validaci.
Vzhledem k výpočetní náročnosti a rozsahu souboru zvolíme $k$-násobnou CV, my zvolíme například hodnotu $k = {10}$.


```r
# model pro vsechna trenovaci data pro K = 1, 2, ..., sqrt(n_train)
neighb.model <- classif.knn(group = y.train, 
                            fdataobj = x.train, 
                            knn = c(1:round(sqrt(length(y.train))))) 

neighb.model$max.prob # maximalni presnost
```

```
## [1] 0.8533333
```

```r
(K.opt <- neighb.model$h.opt) # optimalni hodnota K
```

```
## [1] 1
```

Proveďme předchozí postup pro trénovací data, která rozdělíme na $k$ částí a tedy zopakujeme tuto část kódu $k$-krát.


```r
k_cv <- 10 # k-fold CV
neighbours <- c(1:(2 * ceiling(sqrt(length(y.train))))) # pocet sousedu 

# rozdelime trenovaci data na k casti
folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)

# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro danou cast trenovaci mnoziny
# v radcich budou hodnoty pro danou hodnotu K sousedu
CV.results <- matrix(NA, nrow = length(neighbours), ncol = k_cv)

for (index in 1:k_cv) {
  # definujeme danou indexovou mnozinu
  fold <- folds[[index]]
    
  x.train.cv <- subset(X.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
    fdata()
  y.train.cv <- subset(Y.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
    factor() |> as.numeric()
  
  x.test.cv <- subset(X.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
    fdata()
  y.test.cv <- subset(Y.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
    factor() |> as.numeric()
  
  # projdeme kazdou cast ... k-krat zopakujeme
  for(neighbour in neighbours) {
    # model pro konkretni volbu K
    neighb.model <- classif.knn(group = y.train.cv, 
                              fdataobj = x.train.cv, 
                              knn = neighbour) 
    # predikce na validacni casti
    model.neighb.predict <- predict(neighb.model, 
                                    new.fdataobj = x.test.cv)
    # presnost na validacni casti
    presnost <- table(y.test.cv, model.neighb.predict) |> 
      prop.table() |> diag() |> sum()
    
    # presnost vlozime na pozici pro dane K a fold
    CV.results[neighbour, index] <- presnost
  }
}

# spocitame prumerne presnosti pro jednotliva K pres folds
CV.results <- apply(CV.results, 1, mean)
K.opt <- which.max(CV.results)
presnost.opt.cv <- max(CV.results)
CV.results <- data.frame(K = neighbours, CV = CV.results)
```

Vidíme, že nejlépe vychází hodnota parametru $K$ jako 1 s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.8522.
Pro přehlednost si ještě vykresleme průběh validační chybovosti v závislosti na počtu sousedů $K$.


```r
CV.results |> ggplot(aes(x = K, y = 1 - CV)) + 
  geom_line(linetype = 'dashed', colour = 'grey') + 
  geom_point(size = 1.5) + 
  geom_point(aes(x = K.opt, y = 1 - presnost.opt.cv), colour = 'red', size = 2) +
  theme_bw() + 
  labs(x = bquote(paste(K, ' ;   ', 
                        K[optimal] == .(K.opt))),
       y = 'Validační chybovost') + 
  scale_x_continuous(breaks = neighbours)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-16-1.png" alt="Závislost validační chybovosti na hodnotě $K$, tedy na počtu sousedů." width="672" />
<p class="caption">(\#fig:unnamed-chunk-16)Závislost validační chybovosti na hodnotě $K$, tedy na počtu sousedů.</p>
</div>

Nyní známe optimální hodnotu parametru $K$ a tudíž můžeme sestavit finální model.


```r
neighb.model <- classif.knn(group = y.train, fdataobj = x.train, knn = K.opt)

# predikce
model.neighb.predict <- predict(neighb.model, 
                                new.fdataobj = fdata(X.test))

# presnost na testovacich datech
presnost <- table(as.numeric(factor(Y.test)), model.neighb.predict) |>
  prop.table() |>
  diag() |>
  sum()
```

Vidíme tedy, že přesnost modelu sestrojeného pomocí metody $K$ nejbližších sousedů s optimální volbou $K_{optimal}$ rovnou 1, kterou jsme určili cross-validací, je na trénovacích datech rovna 0.1467 a na testovacích datech 0.1692.

K porovnání jendotlivých modelů můžeme použít oba typy chybovostí, pro přehlednost si je budeme ukládat do tabulky.


```r
RESULTS <- data.frame(model = 'KNN', 
                      Err.train = 1 - neighb.model$max.prob,
                      Err.test = 1 - presnost)
```

### Lineární diskriminační analýza

Jako druhou metodu pro sestrojení klasifikátoru budeme uvažovat lineární diskriminační analýzu (LDA).
Jelikož tato metoda nelze aplikovat na funkcionální data, musíme je nejprve diskretizovat, což provedeme pomocí funkcionální analýzy hlavních komponent.
Klasifikační algoritmus následně provedeme na skórech prvních $p$ hlavních komponent.
Počet komponent $p$ zvolíme tak, aby prvních $p$ hlavních komponent dohromady vysvětlovalo alespoň 90 % variability v datech.

Proveďme tedy nejprve funkcionální analýzu hlavních komponent a určeme počet $p$.


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

V tomto konkrétním případě jsme za počet hlavních komponent vzali $p=$ 2, které dohromady vysvětlují 99.57 $\%$ variability v datech.
První hlavní komponenta potom vysvětluje 98.47 % a druhá 1.09 $\%$ variability.
Graficky si můžeme zobrazit hodnoty skórů prvních dvou hlavních komponent, barevně odlišených podle příslušnosti do klasifikační třídy.


```r
data.PCA.train |> ggplot(aes(x = V1, y = V2, colour = Y)) +
  geom_point(size = 1.5) + 
  labs(x = paste('1. hlavní komponenta (vysvětlená variabilita', 
                 round(100 * data.PCA$varprop[1], 2), '%)'),
       y = paste('2. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[2], 2), '%)'),
       colour = 'Obsah tuku') +
  scale_color_discrete(labels = c("malý", "velký")) + 
  theme_bw()
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-20-1.png" alt="Hodnoty skórů prvních dvou hlavních komponent pro trénovací data. Barevně jsou odlišeny body podle příslušnosti do klasifikační třídy." width="672" />
<p class="caption">(\#fig:unnamed-chunk-20)Hodnoty skórů prvních dvou hlavních komponent pro trénovací data. Barevně jsou odlišeny body podle příslušnosti do klasifikační třídy.</p>
</div>

Abychom mohli určit přesnost klasifikace na testovacích datech, potřebujeme spočítat skóre pro první 2 hlavní komponenty pro testovací data.
Tato skóre určíme pomocí vzorce:

$$
\xi_{i, j} = \int \left( X_i(t) - \mu(t)\right) \cdot \rho_j(t)\text dt,
$$ 
kde $\mu(t)$ je střední hodnota (průměrná funkce) a $\rho_j(t)$ vlastní fukce (funkcionální hlavní komponenty).


```r
# vypocet skoru testovacich funkci
scores <- matrix(NA, ncol = nharm, nrow = length(Y.test)) # prazdna matice

for(k in 1:dim(scores)[1]) {
  xfd = X.test[k] - data.PCA$meanfd[1] # k-te pozorovani - prumerna funkce
  scores[k, ] = inprod(xfd, data.PCA$harmonics) 
  # skalarni soucin rezidua a vlastnich funkci rho (funkcionalni hlavni komponenty)
}

data.PCA.test <- as.data.frame(scores)
data.PCA.test$Y <- factor(Y.test)
colnames(data.PCA.test) <- colnames(data.PCA.train) 
```

Nyní již můžeme sestrojit klasifikátor na trénovací části dat.


```r
# model
clf.LDA <- lda(Y ~ ., data = data.PCA.train)

# presnost na trenovacich datech
predictions.train <- predict(clf.LDA, newdata = data.PCA.train)
presnost.train <- table(data.PCA.train$Y, predictions.train$class) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.LDA, newdata = data.PCA.test)
presnost.test <- table(data.PCA.test$Y, predictions.test$class) |>
  prop.table() |> diag() |> sum()
```

Spočítali jsme jednak přesnost klasifikátoru na trénovacích (68 %), tak i na testovacích datech (70.77 %).

Pro grafické znázornění metody můžeme zaznačit dělící hranici do grafu skórů prvních dvou hlavních komponent.
Tuto hranici spočítáme na husté síti bodů a zobrazíme ji pomocí funkce `geom_contour()`.


```r
# pridame diskriminacni hranici
np <- 1001 # pocet bodu site
# x-ova osa ... 1. HK
nd.x <- seq(from = min(data.PCA.train$V1), 
            to = max(data.PCA.train$V1), length.out = np)
# y-ova osa ... 2. HK
nd.y <- seq(from = min(data.PCA.train$V2), 
            to = max(data.PCA.train$V2), length.out = np)
# pripad pro 2 HK ... p = 2
nd <- expand.grid(V1 = nd.x, V2 = nd.y)
# pokud p = 3
if(dim(data.PCA.train)[2] == 4) {
  nd <- expand.grid(V1 = nd.x, V2 = nd.y, V3 = data.PCA.train$V3[1])}
# pokud p = 4
if(dim(data.PCA.train)[2] == 5) {
  nd <- expand.grid(V1 = nd.x, V2 = nd.y, V3 = data.PCA.train$V3[1],
                    V4 = data.PCA.train$V4[1])}
# pokud p = 5
if(dim(data.PCA.train)[2] == 6) {
  nd <- expand.grid(V1 = nd.x, V2 = nd.y, V3 = data.PCA.train$V3[1],
                    V4 = data.PCA.train$V4[1], V5 = data.PCA.train$V5[1])}

# pridame Y = 0, 1
nd <- nd |> mutate(prd = as.numeric(predict(clf.LDA, newdata = nd)$class))

data.PCA.train |> ggplot(aes(x = V1, y = V2, colour = Y)) +
  geom_point(size = 1.5) + 
  labs(x = paste('1. hlavní komponenta (vysvětlená variabilita', 
                 round(100 * data.PCA$varprop[1], 2), '%)'),
       y = paste('2. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[2], 2), '%)'),
       colour = 'Obsah tuku') +
  scale_color_discrete(labels = c("malý", "velký")) + 
  theme_bw() +
  geom_contour(data = nd, aes(x = V1, y = V2, z = prd), colour = 'black')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-23-1.png" alt="Skóre prvních dvou hlavních komponent, barevně odlišené podle pohlaví. Černě je vyznačena dělící hranice (přímka v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí LDA." width="672" />
<p class="caption">(\#fig:unnamed-chunk-23)Skóre prvních dvou hlavních komponent, barevně odlišené podle pohlaví. Černě je vyznačena dělící hranice (přímka v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí LDA.</p>
</div>

Vidíme, že dělící hranicí je přímka, lineární funkce v prostoru 2D, což jsme ostatně od LDA čekali.
Nakonec přidáme chybovosti do souhrnné tabulky.


```r
Res <- data.frame(model = 'LDA', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

### Kvadratická diskriminační analýza

Jako další sestrojme klasifikátor pomocí kvadratické diskriminační analýzy (QDA).
Jedná se o analogický případ jako LDA s tím rozdílem, že nyní připouštíme pro každou ze tříd rozdílnou kovarianční matici normálního rozdělení, ze kterého pocházejí příslušné skóry.
Tento vypuštěný předpoklad o rovnosti kovariančních matic vede ke kvadratické hranici mezi třídami.

V `R` se provede QDA analogicky jako LDA v předchozí části, tedy opět bychom pomocí funkcionální analýzy hlavních komponent spočítali skóre pro trénovací i testovací funkce, sestrojili klasifikátor na skórech prvních $p$ hlavních komponent a pomocí něj predikovali příslušnost testovacích křivek do třídy $Y^* \in \{0, 1\}$.

Funkcionální PCA provádět nemusíme, využijeme výsledků z části LDA.





Můžeme tedy rovnou přistoupit k sestrojení klasifikátoru, což provedeme pomocí funkce `qda()`.
Následně spočítáme přesnost klasifikátoru na testovacích a trénovacích datech.


```r
# model
clf.QDA <- qda(Y ~ ., data = data.PCA.train)

# presnost na trenovacich datech
predictions.train <- predict(clf.QDA, newdata = data.PCA.train)
presnost.train <- table(data.PCA.train$Y, predictions.train$class) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.QDA, newdata = data.PCA.test)
presnost.test <- table(data.PCA.test$Y, predictions.test$class) |>
  prop.table() |> diag() |> sum()
```

Spočítali jsme tedy jednak přesnost klasifikátoru na trénovacích (68 %), tak i na testovacích datech (69.23 %).

Pro grafické znázornění metody můžeme zaznačit dělící hranici do grafu skórů prvních dvou hlavních komponent.
Tuto hranici spočítáme na husté síti bodů a zobrazíme ji pomocí funkce `geom_contour()` stejně jako v případě LDA.




```r
nd <- nd |> mutate(prd = as.numeric(predict(clf.QDA, newdata = nd)$class))

data.PCA.train |> ggplot(aes(x = V1, y = V2, colour = Y)) +
  geom_point(size = 1.5) + 
  labs(x = paste('1. hlavní komponenta (vysvětlená variabilita', 
                 round(100 * data.PCA$varprop[1], 2), '%)'),
       y = paste('2. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[2], 2), '%)'),
       colour = 'Obsah tuku') +
  scale_color_discrete(labels = c("malý", "velký")) + 
  theme_bw() +
  geom_contour(data = nd, aes(x = V1, y = V2, z = prd), colour = 'black')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-29-1.png" alt="Skóre prvních dvou hlavních komponent, barevně odlišené podle pohlaví. Černě je vyznačena dělící hranice (parabola v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí QDA." width="672" />
<p class="caption">(\#fig:unnamed-chunk-29)Skóre prvních dvou hlavních komponent, barevně odlišené podle pohlaví. Černě je vyznačena dělící hranice (parabola v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí QDA.</p>
</div>

Všimněme si, že dělící hranicí mezi klasifikačními třídami je nyní parabola, avšak se jen (alespoň opticky) velmi málo liší od přímky.

Nakonec ještě doplníme chybovosti do souhrnné tabulky.


```r
Res <- data.frame(model = 'QDA', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

### Logistická regrese

Logistickou regresi můžeme provést dvěma způsoby.
Jednak použít funkcionální obdobu klasické logistické regrese, druhak klasickou mnohorozměrnou logistickou regresi, kterou provedeme na skórech prvních $p$ hlavních komponent.

#### Funkcionální logistická regrese

Analogicky jako v případě konečné dimenze vstupních dat uvažujeme logistický model ve tvaru:

$$
g\left(\mathbb E [Y|X = x]\right) = \eta (x) = g(\pi(x)) = \alpha + \int \beta(t)\cdot x(t) \text d t,
$$ 
kde $\eta(x)$ je lineární prediktor nabývající hodnot z intervalu $(-\infty, \infty)$, $g(\cdot)$ je *linková funkce*, v případě logistické regrese se jedná o logitovou funkci $g: (0,1) \rightarrow \mathbb R,\ g(p) = \ln\frac{p}{1-p}$ a $\pi(x)$ podmíněná pravděpodobnost

$$
\pi(x) = \text{Pr}(Y = 1 | X = x) = g^{-1}(\eta(x)) = \frac{\text e^{\alpha + \int \beta(t)\cdot x(t) \text d t}}{1 + \text e^{\alpha + \int \beta(t)\cdot x(t) \text d t}},
$$

přičemž $\alpha$ je konstanta a $\beta(t) \in L^2[a, b]$ je parametrická funkce.
Naším cílem je odhadnout tuto parametrickou funkci.

Pro funkcionální logistickou regresi použijeme funkci `fregre.glm()` z balíčku `fda.usc`.
Nejprve si vytvoříme vhodné objekty pro konstrukci klasifikátoru.


```r
# vytvorime vhodne objekty
x.train <- fdata(X.train)
y.train <- as.numeric(Y.train) 

# body, ve kterych jsou funkce vyhodnoceny
tt <- x.train[["argvals"]]

dataf <- as.data.frame(y.train) 
colnames(dataf) <- "Y"

nbasis.x <- 7

# B-spline baze 
basis1 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.x)
```

Abychom mohli odhadnout parametrickou funkci $\beta(t)$, potřebujeme ji vyjádřit v nějaké bazické reprezentaci, v našem případě B-splinové bázi.
K tomu však potřebujeme najít vhodný počet bázových funkcí.
To bychom mohli určit na základě chybovosti na trénovacích datech, avšak tato data budou upřenostňovat výběr velkého počtu bází a bude docházet k přeučení modelu.

Ilustrujme si to na následujícím případě.
Pro každý z počtu bází $n_{basis} \in \{4, 5, \dots, 30\}$ natrénujeme model na trénovacích datech, určíme na nich chybovost a také spočítáme chybovost na testovacích datech.
Připomeňme, že k výběru vhodného počtu bází nemůžeme využít stejná data jako pro odhad testovací chybovosti, neboť bychom tuto chybovost podcenili.


```r
n.basis.max <- 30
n.basis <- 4:n.basis.max
pred.baz <- matrix(NA, nrow = length(n.basis), ncol = 2, 
                   dimnames = list(n.basis, c('Err.train', 'Err.test')))

for (i in n.basis) {
  # baze pro bety
  basis2 <- create.bspline.basis(rangeval = range(tt), nbasis = i)
  # vztah
  f <- Y ~ x
  # baze pro x a bety
  basis.x <- list("x" = basis1) # vyhlazene data
  basis.b <- list("x" = basis2)
  # vstupni data do modelu
  ldata <- list("df" = dataf, "x" = x.train)
  # binomicky model ... model logisticke regrese
  model.glm <- fregre.glm(f, family = binomial(), data = ldata,
                          basis.x = basis.x, basis.b = basis.b)
  
  # presnost na trenovacich datech
  predictions.train <- predict(model.glm, newx = ldata)
  predictions.train <- data.frame(Y.pred = ifelse(predictions.train < 1/2, 0, 1))
  presnost.train <- table(Y.train, predictions.train$Y.pred) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  newldata = list("df" = as.data.frame(Y.test), "x" = fdata(X.test))
  predictions.test <- predict(model.glm, newx = newldata)
  predictions.test <- data.frame(Y.pred = ifelse(predictions.test < 1/2, 0, 1))
  presnost.test <- table(Y.test, predictions.test$Y.pred) |>
    prop.table() |> diag() |> sum()
  
  # vlozime do matice
  pred.baz[as.character(i), ] <- 1 - c(presnost.train, presnost.test)
} 

pred.baz <- as.data.frame(pred.baz)
pred.baz$n.basis <- n.basis
```

Znázorněme si průběh obou typů chybovostí v grafu v závislosti na počtu bazických funkcí.


```r
n.basis.beta.opt <- pred.baz$n.basis[which.min(pred.baz$Err.test)]

pred.baz |> ggplot(aes(x = n.basis, y = Err.test)) + 
  geom_line(linetype = 'dashed', colour = 'black') + 
  geom_line(aes(x = n.basis, y = Err.train), colour = 'deepskyblue3', 
            linetype = 'dashed', linewidth = 0.5) + 
  geom_point(size = 1.5) + 
  geom_point(aes(x = n.basis, y = Err.train), colour = 'deepskyblue3', 
             size = 1.5) + 
  geom_point(aes(x = n.basis.beta.opt, y = min(pred.baz$Err.test)),
             colour = 'red', size = 2) +
  theme_bw() + 
  labs(x = bquote(paste(n[basis], ' ;   ', 
                        n[optimal] == .(n.basis.beta.opt))),
        y = 'Chybovost')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-33-1.png" alt="Závislost testovací a trénovací chybovosti na počtu bázových funkcí pro $\beta$. Červeným bodem je znázorněn optimální počet $n_{optimal}$ zvolený jako minimum testovací chybovosti, černou čarou je vykreslena testovací a modrou přerušovanou čarou je vykreslen průběh trénovací chybovosti." width="672" />
<p class="caption">(\#fig:unnamed-chunk-33)Závislost testovací a trénovací chybovosti na počtu bázových funkcí pro $\beta$. Červeným bodem je znázorněn optimální počet $n_{optimal}$ zvolený jako minimum testovací chybovosti, černou čarou je vykreslena testovací a modrou přerušovanou čarou je vykreslen průběh trénovací chybovosti.</p>
</div>

Vidíme, že s rostoucím počtem bází pro $\beta(t)$ má trénovací chybovost (modrá čára) tendenci klesat a tedy bychom na jejím základě volili velké hodnoty $n_{basis}$.
Naopak optimální volbou na základě testovací chybovosti je $n$ rovno 6, tedy výrazně menší hodnota než 30.
Naopak s rostoucím $n$ roste testovací chyvost, což ukazuje na přeučení modelu.

Z výše uvedených důvodů pro určení optimálního počtu bazických funkcí pro $\beta(t)$ využijeme 10-ti násobnou cross-validaci.
Jako maximální počet uvažovaných bazických funkcí bereme 25, neboť jak jsme viděli výše, nad touto hodnotou dochází již k přeučení modelu.


```r
### 10-fold cross-validation
n.basis.max <- 25
n.basis <- 4:n.basis.max
k_cv <- 10 # k-fold CV
# rozdelime trenovaci data na k casti
folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)
## prvky, ktere se behem cyklu nemeni
# body, ve kterych jsou funkce vyhodnoceny
tt <- x.train[["argvals"]]
rangeval <- range(tt)
# B-spline baze 
basis1 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.x)
# vztah
f <- Y ~ x
# baze pro x
basis.x <- list("x" = basis1)
# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro danou cast trenovaci mnoziny
# v radcich budou hodnoty pro dany pocet bazi
CV.results <- matrix(NA, nrow = length(n.basis), ncol = k_cv, 
                     dimnames = list(n.basis, 1:k_cv))
```

Nyní již máme vše připravené pro spočítání chybovosti na každé z deseti podmnožin trénovací množiny.
Následně určíme průměr a jako optimální $n$ vezmeme argument minima validační chybovosti.


```r
for (index in 1:k_cv) {
  # definujeme danou indexovou mnozinu
  fold <- folds[[index]]
    
  x.train.cv <- subset(X.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
    fdata()
  y.train.cv <- subset(Y.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
    as.numeric()
  
  x.test.cv <- subset(X.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
    fdata()
  y.test.cv <- subset(Y.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
    as.numeric()
  
  dataf <- as.data.frame(y.train.cv) 
  colnames(dataf) <- "Y"
  
  for (i in n.basis) {
    # baze pro bety
    basis2 <- create.bspline.basis(rangeval = rangeval, nbasis = i)
    
    basis.b <- list("x" = basis2)
    # vstupni data do modelu
    ldata <- list("df" = dataf, "x" = x.train.cv)
    # binomicky model ... model logisticke regrese
    model.glm <- fregre.glm(f, family = binomial(), data = ldata,
                            basis.x = basis.x, basis.b = basis.b)
    
    # presnost na validacni casti 
    newldata = list("df" = as.data.frame(y.test.cv), "x" = x.test.cv)
    predictions.valid <- predict(model.glm, newx = newldata)
    predictions.valid <- data.frame(Y.pred = ifelse(predictions.valid < 1/2, 0, 1))
    presnost.valid <- table(y.test.cv, predictions.valid$Y.pred) |>
      prop.table() |> diag() |> sum()
    
    # vlozime do matice
    CV.results[as.character(i), as.character(index)] <- presnost.valid
  } 
}

# spocitame prumerne presnosti pro jednotliva n pres folds
CV.results <- apply(CV.results, 1, mean)
n.basis.opt <- n.basis[which.max(CV.results)]
presnost.opt.cv <- max(CV.results)
```

Vykresleme si ještě průběh validační chybovosti i se zvýrazněnou optimální hodnotou $n_{optimal}$ rovnou 24 s validační chybovostí 0.0488.


```r
CV.results <- data.frame(n.basis = n.basis, CV = CV.results)
CV.results |> ggplot(aes(x = n.basis, y = 1 - CV)) + 
  geom_line(linetype = 'dashed', colour = 'grey') + 
  geom_point(size = 1.5) + 
  geom_point(aes(x = n.basis.opt, y = 1 - presnost.opt.cv), colour = 'red', size = 2) +
  theme_bw() + 
  labs(x = bquote(paste(n[basis], ' ;   ', 
                        n[optimal] == .(n.basis.opt))),
       y = 'Validační chybovost') + 
  scale_x_continuous(breaks = n.basis) + 
  theme(panel.grid.minor = element_blank())
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-36-1.png" alt="Závislost validační chybovosti na hodnotě $n_{basis}$, tedy na počtu bází." width="672" />
<p class="caption">(\#fig:unnamed-chunk-36)Závislost validační chybovosti na hodnotě $n_{basis}$, tedy na počtu bází.</p>
</div>

Nyní již tedy můžeme definovat finální model pomocí funkcionální logistické regrese, přičemž bázi pro $\beta(t)$ volíme B-splinovou bázi s 24 bázemi.


```r
# optimalni model
basis2 <- create.bspline.basis(rangeval = range(tt), nbasis = n.basis.opt)
f <- Y ~ x
# baze pro x a bety
basis.x <- list("x" = basis1) 
basis.b <- list("x" = basis2)
# vstupni data do modelu
dataf <- as.data.frame(y.train) 
colnames(dataf) <- "Y"
ldata <- list("df" = dataf, "x" = x.train)
# binomicky model ... model logisticke regrese
model.glm <- fregre.glm(f, family = binomial(), data = ldata,
                        basis.x = basis.x, basis.b = basis.b,
                        maxit = 1000, epsilon = 1e-2)

# presnost na trenovacich datech
predictions.train <- predict(model.glm, newx = ldata)
predictions.train <- data.frame(Y.pred = ifelse(predictions.train < 1/2, 0, 1))
presnost.train <- table(Y.train, predictions.train$Y.pred) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
newldata = list("df" = as.data.frame(Y.test), "x" = fdata(X.test))
predictions.test <- predict(model.glm, newx = newldata)
predictions.test <- data.frame(Y.pred = ifelse(predictions.test < 1/2, 0, 1))
presnost.test <- table(Y.test, predictions.test$Y.pred) |>
  prop.table() |> diag() |> sum()
```

Spočítali jsme trénovací chybovost (rovna 0 %) i testovací chybovost (rovna 6.15 %).
Pro lepší představu si ještě můžeme vykreslit hodnoty odhadnutých pravděpodobností příslušnosti do klasifikační třídy $Y = 1$ na trénovacích datech v závislosti na hodnotách lineárního prediktoru.


```r
data.frame(
  linear.predictor = model.glm$linear.predictors,
  response = model.glm$fitted.values,
  Y = factor(y.train)
) |> ggplot(aes(x = linear.predictor, y = response, colour = Y)) + 
  geom_point(size = 1.5) + 
  scale_color_discrete(labels = c("malý", "velký")) + 
  geom_abline(aes(slope = 0, intercept = 0.5), linetype = 'dashed') + 
  theme_bw() + 
  labs(x = 'Lineární prediktor',
       y = 'Odhadnuté pravděpodobnosti Pr(Y = 1|X = x)',
       colour = 'Obsah tuku') 
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-38-1.png" alt="Závoslost odhadnutých pravděpodobností na hodnotách lineárního prediktoru. Barevně jsou odlišeny body podle příslušnosti do klasifikační třídy." width="672" />
<p class="caption">(\#fig:unnamed-chunk-38)Závoslost odhadnutých pravděpodobností na hodnotách lineárního prediktoru. Barevně jsou odlišeny body podle příslušnosti do klasifikační třídy.</p>
</div>

Můžeme si ještě pro informaci zobrazit průběh odhadnuté parametrické funkce $\beta(t)$.


```r
t.seq <- seq(min(t), max(t), length = 1001)
beta.seq <- eval.fd(evalarg = t.seq, fdobj = model.glm$beta.l$x)

data.frame(t = t.seq, beta = beta.seq) |> 
  ggplot(aes(t, beta)) + 
  geom_abline(aes(slope = 0, intercept = 0), linetype = 'dashed', 
              linewidth = 0.5, colour = 'grey') +
  geom_line() + 
  theme_bw() +
  labs(x = expression(x[1]),
       y = expression(widehat(beta)(t))) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-39-1.png" alt="Průběh odhadu parametrické funkce $\beta(t), t \in [850, 1050]$." width="672" />
<p class="caption">(\#fig:unnamed-chunk-39)Průběh odhadu parametrické funkce $\beta(t), t \in [850, 1050]$.</p>
</div>

Vidíme, že hodnoty funkce $\hat\beta(t)$ se drží kolem nuly pro časy $t$ z prostředka  a konce intervalu $[850, 1050]$, zatímco pro dřívější časy jsou hodnoty vyšší.
To implikuje rozdílnost funkcí z klasifikačních tříd na začátku  intervalu, zatímco uprostřed a na konci intervalu jsou funkce velmi podobné.

Výsledky opět přidáme do souhrnné tabulky.


```r
Res <- data.frame(model = 'LR functional', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

#### Logistická regrese s analýzou hlavních komponent

Abychom mohli sestrojit tento klasifikátor, potřebujeme provést funkcionální analýzu hlavních komponent, určit vhodný počet komponent a spočítat hodnoty skórů pro testovací data.
To jsme již provedli v části u lineární diskriminační analýzy, proto využijeme tyto výsledky v následující části.





Můžeme tedy rovnou sestrojit model logistické regrese pomocí funkce `glm(, family = binomial)`.


```r
# model
clf.LR <- glm(Y ~  ., data = data.PCA.train, family = binomial)

# presnost na trenovacich datech
predictions.train <- predict(clf.LR, newdata = data.PCA.train, type = 'response')
predictions.train <- ifelse(predictions.train > 0.5, 1, 0)
presnost.train <- table(data.PCA.train$Y, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.LR, newdata = data.PCA.test, type = 'response')
predictions.test <- ifelse(predictions.test > 0.5, 1, 0)
presnost.test <- table(data.PCA.test$Y, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Spočítali jsme tedy přesnost klasifikátoru na trénovacích (69.33 %) i na testovacích datech (70.77 %).

Pro grafické znázornění metody můžeme zaznačit dělící hranici do grafu skórů prvních dvou hlavních komponent.
Tuto hranici spočítáme na husté síti bodů a zobrazíme ji pomocí funkce `geom_contour()` stejně jako v případě LDA i QDA.




```r
nd <- nd |> mutate(prd = as.numeric(predict(clf.LR, newdata = nd,
                                            type = 'response')))
nd$prd <- ifelse(nd$prd > 0.5, 1, 0)

data.PCA.train |> ggplot(aes(x = V1, y = V2, colour = Y)) +
  geom_point(size = 1.5) + 
  labs(x = paste('1. hlavní komponenta (vysvětlená variabilita', 
                 round(100 * data.PCA$varprop[1], 2), '%)'),
       y = paste('2. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[2], 2), '%)'),
       colour = 'Obsah tuku') +
  scale_colour_discrete(labels = c("malý", "velký")) +
  theme_bw() +
  geom_contour(data = nd, aes(x = V1, y = V2, z = prd), colour = 'black')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-45-1.png" alt="Skóre prvních dvou hlavních komponent, barevně odlišené podle příslušnosti do klasifikační třídy. Černě je vyznačena dělící hranice (přímka v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí logistické regrese." width="672" />
<p class="caption">(\#fig:unnamed-chunk-45)Skóre prvních dvou hlavních komponent, barevně odlišené podle příslušnosti do klasifikační třídy. Černě je vyznačena dělící hranice (přímka v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí logistické regrese.</p>
</div>

Všimněme si, že dělící hranicí mezi klasifikačními třídami je nyní přímka jako v případě LDA.

Nakonec ještě doplníme chybovosti do souhrnné tabulky.


```r
Res <- data.frame(model = 'LR score', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

### Rozhodovací stromy

V této části se podíváme na velmi odlišný přístup k sestrojení klasifikátoru, než byly například LDA či logistická regrese.
Rozhodovací stromy jsou velmi oblíbeným nástrojem ke klasifikaci, avšak jako v případě některých předchozích metod nejsou přímo určeny pro funkcionální data.
Existují však postupy, jak funkcionální objekty převést na mnohorozměrné a následně na ně aplikovat algoritmus rozhodovacích stromů.
Můžeme uvažovat následující postupy:

-   algoritmus sestrojený na bázových koeficientech,

-   využití skórů hlavních komponent,

-   použít diskretizaci intervalu a vyhodnotit funkci jen na nějaké konečné síti bodů.

My se nejprve zaměříme na diskretizaci intervalu a následně porovnáme výsledky se zbylými dvěma přístupy k sestrojení rozhodovacího stromu.

#### Diskretizace intervalu

Nejprve si musíme definovat body z intervalu $I = [850, 1050]$, ve kterých funkce vyhodnotíme.
Následně vytvoříme objekt, ve kterém budou řádky představovat jednotlivé (diskretizované) funkce a sloupce časy.
Nakonec připojíme sloupec $Y$ s informací o příslušnosti do klasifikační třídy a totéž zopakujeme i pro testovací data.


```r
# posloupnost bodu, ve kterych funkce vyhodnotime
t.seq <- seq(min(t), max(t), length = 101)
   
grid.data <- eval.fd(fdobj = X.train, evalarg = t.seq)
grid.data <- as.data.frame(t(grid.data)) # transpozice kvuli funkcim v radku
grid.data$Y <- Y.train |> factor()

grid.data.test <- eval.fd(fdobj = X.test, evalarg = t.seq)
grid.data.test <- as.data.frame(t(grid.data.test))
grid.data.test$Y <- Y.test |> factor()
```

Nyní mážeme sestrojit rozhodovací strom, ve kterém budou jakožto prediktory vystupovat všechny časy z vektoru `t.seq`.
Tato klasifikační není náchylná na multikolinearitu, tudíž se jí nemusíme zabývat.
Jako metriku zvolíme přesnost.


```r
# sestrojeni modelu
clf.tree <- train(Y ~ ., data = grid.data, 
                 method = "rpart", 
                 trControl = trainControl(method = "CV", number = 10),
                 metric = "Accuracy")

# presnost na trenovacich datech
predictions.train <- predict(clf.tree, newdata = grid.data)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.tree, newdata = grid.data.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost klasifikátoru na testovacích datech je tedy 64.62 % a na trénovacích datech 70.67 %.

Graficky si rozhodovací strom můžeme vykreslit pomocí funkce `fancyRpartPlot()`.
Nastavíme barvy uzlů tak, aby reflektovaly předchozí barevné odlišení.
Jedná se o neprořezaný strom.


```r
colnames(grid.data) <- c(paste0('time:', t.seq), 'Y')
fancyRpartPlot(rpart(Y ~ ., data = grid.data, method = "class"),
               sub = '', palettes = c('Reds', 'Blues')) 
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-49-1.png" alt="Grafické znázornění neprořezaného rozhodovacího stromu. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0." width="672" />
<p class="caption">(\#fig:unnamed-chunk-49)Grafické znázornění neprořezaného rozhodovacího stromu. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0.</p>
</div>

Můžeme si také vykreslit již prořezaný finální rozhodovací strom.


```r
rpart.plot::rpart.plot(clf.tree$finalModel, # finalni model ... prorezany strom
                       extra = 104, # zobrazeni pozadovanych informaci
                       box.palette = c('Reds', 'Blues'),
                       branch.lty = 3, # dotted branch lines
                       shadow.col = 0, # shadows under the node boxes
                       nn = FALSE, 
                       under = FALSE, 
                       digits = 2)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-50-1.png" alt="Finální prořezaný rozhodovací strom." width="672" />
<p class="caption">(\#fig:unnamed-chunk-50)Finální prořezaný rozhodovací strom.</p>
</div>

Nakonec opět přidejme trénovací a testovací chybovost do souhrnné tabulky.


```r
Res <- data.frame(model = 'Tree - diskr.', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

#### Skóre hlavních komponent

Další možností pro sestrojení rozhodovacího stromu je použít skóre hlavních komponent.
Jelikož jsme již skóre počítali pro předchozí klasifikační metody, využijeme těchto poznatků a sestrojíme rozhodovací strom na skórech prvních 2 hlavních komponent.


```r
# sestrojeni modelu
clf.tree.PCA <- train(Y ~ ., data = data.PCA.train, 
                 method = "rpart", 
                 trControl = trainControl(method = "CV", number = 10),
                 metric = "Accuracy")

# presnost na trenovacich datech
predictions.train <- predict(clf.tree.PCA, newdata = data.PCA.train)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.tree.PCA, newdata = data.PCA.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost rozhodovacího stromu na testovacích datech je tedy 61.54 % a na trénovacích datech 67.33 %.

Graficky si rozhodovací strom sestrojený na skórech hlavních komponent můžeme vykreslit pomocí funkce `fancyRpartPlot()`.
Nastavíme barvy uzlů tak, aby reflektovaly předchozí barevné odlišení.
Jedná se o neprořezaný strom.


```r
fancyRpartPlot(rpart(Y ~ ., data = data.PCA.train, method = "class"),
               sub = '', palettes = c('Reds', 'Blues')) 
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-53-1.png" alt="Grafické znázornění neprořezaného rozhodovacího stromu sestrojeného na skórech hlavních komponent. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0." width="672" />
<p class="caption">(\#fig:unnamed-chunk-53)Grafické znázornění neprořezaného rozhodovacího stromu sestrojeného na skórech hlavních komponent. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0.</p>
</div>

Můžeme si také vykreslit již prořezaný finální rozhodovací strom.


```r
rpart.plot::rpart.plot(clf.tree.PCA$finalModel, # finalni model 
                       extra = 104, # zobrazeni pozadovanych informaci
                       box.palette = c('Reds', 'Blues'),
                       branch.lty = 3, # dotted branch lines
                       shadow.col = 0, # shadows under the node boxes
                       nn = FALSE, under = FALSE, digits = 2)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-54-1.png" alt="Finální prořezaný rozhodovací strom." width="672" />
<p class="caption">(\#fig:unnamed-chunk-54)Finální prořezaný rozhodovací strom.</p>
</div>

Nakonec opět přidejme trénovací a testovací chybovost do souhrnné tabulky.


```r
Res <- data.frame(model = 'Tree - score', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

#### Bázové koeficienty

Poslední možností, kterou využijeme pro sestrojení rozhodovacího stromu, je použití koeficientů ve vyjádření funkcí v B-splinové bázi.

Nejprve si definujme potřebné datové soubory s koeficienty.


```r
# trenovaci dataset
data.Bbasis.train <- t(X.train$coefs) |> as.data.frame()
data.Bbasis.train$Y <- factor(Y.train)

# testovaci dataset
data.Bbasis.test <- t(X.test$coefs) |> as.data.frame()
data.Bbasis.test$Y <- factor(Y.test)
```

Nyní již můžeme sestrojit klasifikátor.


```r
# sestrojeni modelu
clf.tree.Bbasis <- train(Y ~ ., data = data.Bbasis.train, 
                 method = "rpart", 
                 trControl = trainControl(method = "CV", number = 10),
                 metric = "Accuracy")

# presnost na trenovacich datech
predictions.train <- predict(clf.tree.Bbasis, newdata = data.Bbasis.train)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.tree.Bbasis, newdata = data.Bbasis.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost rozhodovacího stromu na trénovacích datech je tedy 77.33 % a na testovacích datech 75.38 %.

Graficky si rozhodovací strom sestrojený na koeficientech B-splinového vyjádření můžeme vykreslit pomocí funkce `fancyRpartPlot()`.
Nastavíme barvy uzlů tak, aby reflektovaly předchozí barevné odlišení.
Jedná se o neprořezaný strom.


```r
fancyRpartPlot(rpart(Y ~ ., data = data.Bbasis.train, method = "class"),
               sub = '', palettes = c('Reds', 'Blues')) 
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-58-1.png" alt="Grafické znázornění neprořezaného rozhodovacího stromu sestrojeného na bázových koeficientech. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0." width="672" />
<p class="caption">(\#fig:unnamed-chunk-58)Grafické znázornění neprořezaného rozhodovacího stromu sestrojeného na bázových koeficientech. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0.</p>
</div>

Můžeme si také vykreslit již prořezaný finální rozhodovací strom.


```r
rpart.plot::rpart.plot(clf.tree.Bbasis$finalModel, # finalni model 
                       extra = 104, # zobrazeni pozadovanych informaci
                       box.palette = c('Reds', 'Blues'),
                       branch.lty = 3, # dotted branch lines
                       shadow.col = 0, # shadows under the node boxes
                       nn = FALSE, under = FALSE, digits = 2)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-59-1.png" alt="Finální prořezaný rozhodovací strom." width="672" />
<p class="caption">(\#fig:unnamed-chunk-59)Finální prořezaný rozhodovací strom.</p>
</div>

Nakonec opět přidejme trénovací a testovací chybovost do souhrnné tabulky.


```r
Res <- data.frame(model = 'Tree - Bbasis', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

### Náhodné lesy

Klasifikátor sestrojený pomocí metody náhodných lesů spočívá v sestrojení několika jednotlivých rozhodovacích stromů, které se následně zkombinují a vytvoří společný klasifikátor (společným "hlasováním").

Tak jako v případě rozhodovacích stromů máme několik možností na to, jaká data (konečně-rozměrná) použijeme pro sestrojení modelu.
Budeme opět uvažovat výše diskutované tři přístupy.
Datové soubory s příslušnými veličinami pro všechny tři přístupy již máme připravené z minulé sekce, proto můžeme přímo sestrojit dané modely, spočítat charakteristiky daného klasifikátoru a přidat výsledky do souhrnné tabulky.

#### Diskretizace intervalu

V prvním případě využíváme vyhodnocení funkcí na dané síti bodů intervalu $I = [850, 1050]$.




```r
# sestrojeni modelu
clf.RF <- randomForest(Y ~ ., data = grid.data, 
                       ntree = 500, # pocet stromu
                       importance = TRUE,
                       nodesize = 5)

# presnost na trenovacich datech
predictions.train <- predict(clf.RF, newdata = grid.data)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.RF, newdata = grid.data.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost náhodného lesu na trénovacích datech je tedy 98 % a na testovacích datech 87.69 %.


```r
Res <- data.frame(model = 'RForest - diskr', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

#### Skóre hlavních komponent

V tomto případě využijeme skóre prvních $p =$ 2 hlavních komponent.


```r
# sestrojeni modelu
clf.RF.PCA <- randomForest(Y ~ ., data = data.PCA.train, 
                           ntree = 500, # pocet stromu
                           importance = TRUE,
                           nodesize = 5)

# presnost na trenovacich datech
predictions.train <- predict(clf.RF.PCA, newdata = data.PCA.train)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.RF.PCA, newdata = data.PCA.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost rozhodovacího stromu na trénovacích datech je tedy 95.33 % a na testovacích datech 69.23 %.


```r
Res <- data.frame(model = 'RForest - score', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

#### Bázové koeficienty

Nakonec použijeme vyjádření funkcí pomocí B-splinové báze.


```r
# sestrojeni modelu
clf.RF.Bbasis <- randomForest(Y ~ ., data = data.Bbasis.train, 
                              ntree = 500, # pocet stromu
                              importance = TRUE,
                              nodesize = 5)

# presnost na trenovacich datech
predictions.train <- predict(clf.RF.Bbasis, newdata = data.Bbasis.train)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.RF.Bbasis, newdata = data.Bbasis.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost tohoto klasifikátoru na trénovacích datech je 98.67 % a na testovacích datech 87.69 %.


```r
Res <- data.frame(model = 'RForest - Bbasis', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

### Support Vector Machines

Nyní se podívejme na klasifikaci našich nasimulovaných křivek pomocí metody podpůrných vektorů (ang. Support Vector Machines, SVM).
Výhodou této klasifikační metody je její výpočetní nenáročnost, neboť pro definici hraniční křivky mezi třídami využívá pouze několik (často málo) pozorování.

Hlavní výhodou SVM je použití tzv.
*jádrového triku* (kernel trick), pomocí kterého nahradíme obyčejný skalární součin jiným skalárním součinem transformovaných dat, aniž bychom tuto transformaci museli přímo definovat.
Tím dostaneme obecně nelineární dělící hranici mezi klasifikačními třídami.
*Jádro* (jádrová funkce, ang. kernel, kernel function) $K$ je taková funkce, která splňuje

$$
K(x_i, x_j) = \langle \phi(x_i), \phi(x_j) \rangle_{\mathcal H}, 
$$
kde $\phi$ je nějaká (neznámá) transformace (ang. feature map), $\mathcal H$ je Hilbertův prostor a $\langle \cdot, \cdot \rangle_{\mathcal H}$ je nějaký skalární součin na tomto Hilbertově prostoru.

Nejčastěji se v praxi volí tři typy jádrových funkcí:

-   lineární jádro -- $K(x_i, x_j) = \langle x_i, x_j \rangle$,
-   polynomiální jádro -- $K(x_i, x_j) = \big(\alpha_0 + \gamma \langle x_i, x_j \rangle \big)^d$,
-   radiální (gaussovské) jádro -- $\displaystyle{K(x_i, x_j) = \text e^{-\gamma \|x_i - x_j \|^2}}$.

U všech výše zmíněných jader musíme zvolit konstantu $C > 0$, která udává míru penalizace za překročení dělící hranice mezi třídami (ang. inverse regularization parameter).
S rostoucí hodnotou $C$ bude metoda více penalizovat špatně klasifikovaná data a méně tvar hranice, naopak pro malé hodnoty $C$ metoda nedává takový význam špatně klasifikovaným datům, ale zaměřuje se více na penalizaci tvaru hranice.
Tato konstanta $C$ se defaultně volí rovna 1, můžeme ji určit i přímo například pomocí cross-validace.

Využitím cross-validace můžeme také určit optimální hodnoty ostatních hyperparametrů, které nyní závisí na naší volbě jádrové funkce.
V případě lineárního jádra nevolíme žádný další parametr kromě konstanty $C$, u polynomiálního jádra musíme určit hodnoty hyperparametrů $\alpha_0, \gamma \text{ a } d$, jejichž defaultní hodnoty v `R` jsou postupně $\alpha_0^{default} = 0, \gamma^{default} = \frac{1}{dim(\texttt{data})} \text{ a } d^{default} = 3$.
Při volbě radiálního jádra máme pouze jeden další hyperparametr $\gamma$, jehož defaultní hodnota v `R` je totožná jako u polynomiálního jádra.
Opět bychom mohli hodnoty hyperparametrů určit jako optimální pro naše data, avšak vzhledem k relativní výpočetní náročnosti necháme hodnoty příslušných hyperparametrů na jejich defaultních hodnotách.

V případě funkcionálních dat máme několik možností, jak použít metodu SVM.
Nejjednodušší variantou je použít tuto klasifikační metodu přímo na diskretizovanou funkci (sekce \@ref(diskrA3)).
Další možností je opět využít skóre hlavních komponent a klasifikovat křivky pomocí jejich reprezentace \@ref(PCA-SVMA3).
Další přímočarou variantou je využít vyjádření křivek pomocí B-splinové báze a klasifikovat křivky na základě koeficientů jejich vyjádření v této bázi (sekce \@ref(basis-SVMA3)).

Složitější úvahou můžeme dospět k několika dalším možnostem, které využívají funkcionální podstatu dat.
Jednak můžeme místo klasifikace původní křivky využít její derivaci (případně druhou derivaci, třetí, ...), druhak můžeme využít projekce funkcí na podprostor generovaný, např. B-splinovými, funkcemi (sekce \@ref(projection-SVMA3)). Poslední metoda, kterou použijeme pro klasifikaci funkcionálních dat, spočívá v kombinaci projekce na určitý podprostor generovaný funkcemi (Reproducing Kernel Hilbert Space, RKHS) a klasifikace příslušné reprezentace. Tato metoda využívá kromě klasického SVM i SVM pro regresi, více uvádíme v sekci RKHS + SVM \@ref(RKHS-SVMA3).

#### Diskretizace intervalu {#diskrA3}

Začněme nejprve aplikací metody podpůrných vektorů přímo na diskretizovaná data (vyhodnocení funkce na dané síti bodů na intervalu $I = [850, 1050]$), přičemž budeme uvažovat všech tři výše zmíněné jádrové funkce.


```r
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

# rozdeleni na testovaci a trenovaci cast
X.train_norm <- subset(XXfd_norm, split == TRUE)
X.test_norm <- subset(XXfd_norm, split == FALSE)

Y.train_norm <- subset(Y, split == TRUE)
Y.test_norm <- subset(Y, split == FALSE)

grid.data <- eval.fd(fdobj = X.train_norm, evalarg = t.seq)
grid.data <- as.data.frame(t(grid.data)) 
grid.data$Y <- Y.train_norm |> factor()

grid.data.test <- eval.fd(fdobj = X.test_norm, evalarg = t.seq)
grid.data.test <- as.data.frame(t(grid.data.test))
grid.data.test$Y <- Y.test_norm |> factor()
```


Nyní se pokusme, na rozdíl od postupu v předchozích kapitolách, hyperparametry klasifikátorů odhadnout z dat pomocí 10-násobné cross-validace. Vzhledem k tomu, že každé jádro má ve své definici jiné hyperparametry, budeme ke každé jádrové funkci přistupovat zvlášť. Nicméně hyperparametr $C$ vystupuje u všech jádrových funkcí, přičemž ale připouštíme, že se může jeho optimální hodnota mezi jádry lišit.

U všech třech jader projdeme hodnoty hyperparametru $C$ v intervalu $[10^{-3}, 10^{3}]$, přičemž u jádra polynomiálního zafixujeme hyperparametr $p$ na hodnotě 3, neboť pro jiné celočíselné hodnoty metoda nedává zdaleka tak dobré výsledky. Naopak pro radiální jádro využijeme k volbě optimální hodnoty hyperparametru $\gamma$ opět 10-násobnou CV, přičemž uvažujeme hodnoty v intervalu $[10^{-3}, 10^{2}]$. Zvolíme `coef0` $= 1$. 


```r
set.seed(42)

k_cv <- 10 #  k-fold CV

# rozdelime trenovaci data na k casti
folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)

# ktere hodnoty gamma chceme uvazovat
gamma.cv <- 10^seq(-3, 2, length = 15)
C.cv <- 10^seq(-3, 3, length = 20)
p.cv <- 3
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
    cv_sample <- 1:dim(grid.data)[1] %in% fold
    
    data.grid.train.cv <- as.data.frame(grid.data[cv_sample, ])
    data.grid.test.cv <- as.data.frame(grid.data[!cv_sample, ])
    
    ## LINEARNI JADRO
    # sestrojeni modelu
    clf.SVM.l <- svm(Y ~ ., data = data.grid.train.cv,
                     type = 'C-classification',
                     scale = TRUE,
                     cost = C,
                     kernel = 'linear')
    
    # presnost na validacnich datech
    predictions.test.l <- predict(clf.SVM.l, newdata = data.grid.test.cv)
    presnost.test.l <- table(data.grid.test.cv$Y, predictions.test.l) |>
      prop.table() |> diag() |> sum()
    
    # presnosti vlozime na pozice pro dane C a fold
    CV.results$SVM.l[(1:length(C.cv))[C.cv == C], 
                     index_cv] <- presnost.test.l
    
    ## POLYNOMIALNI JADRO
    for (p in p.cv) {
      # sestrojeni modelu
      clf.SVM.p <- svm(Y ~ ., data = data.grid.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       coef0 = coef0,
                       degree = p,
                       kernel = 'polynomial')
      
      # presnost na validacnich datech
      predictions.test.p <- predict(clf.SVM.p, newdata = data.grid.test.cv)
      presnost.test.p <- table(data.grid.test.cv$Y, predictions.test.p) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane C, p a fold
      CV.results$SVM.p[(1:length(C.cv))[C.cv == C], 
                       (1:length(p.cv))[p.cv == p],
                       index_cv] <- presnost.test.p
    }
        
    ## RADIALNI JADRO
    for (gamma in gamma.cv) {
      # sestrojeni modelu
      clf.SVM.r <- svm(Y ~ ., data = data.grid.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       gamma = gamma,
                       kernel = 'radial')
      
      # presnost na validacnich datech
      predictions.test.r <- predict(clf.SVM.r, newdata = data.grid.test.cv)
      presnost.test.r <- table(data.grid.test.cv$Y, predictions.test.r) |>
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

Podívejme se, jak dopadly optimální hodnoty. Pro *lineární jádro* máme optimální hodnotu $C$ rovnu 0.0379, pro *polynomiální jádro* je $C$ rovno 1.4384 a pro *radiální jádro* máme dvě optimální hodnoty, pro $C$ je optimální hodnota 12.7427 a pro $\gamma$ je to 0.0052. Validační přesnosti jsou postupně 0.9933333 pro lineární, 0.9804167 pro polynomiální a 0.9866667 pro radiální jádro.

Konečně můžeme sestrojit finální klasifikátory na celých trénovacích datech s hodnotami hyperparametrů určenými pomocí 10-násobné CV. Určíme také chybovosti na testovacích a také na trénovacích datech.


```r
# sestrojeni modelu
clf.SVM.l <- svm(Y ~ ., data = grid.data,
                 type = 'C-classification',
                 scale = TRUE,
                 cost = C.opt[1],
                 kernel = 'linear')

clf.SVM.p <- svm(Y ~ ., data = grid.data,
                 type = 'C-classification',
                 scale = TRUE,
                 cost = C.opt[2],
                 degree = p.opt,
                 coef0 = coef0,
                 kernel = 'polynomial')

clf.SVM.r <- svm(Y ~ ., data = grid.data,
                 type = 'C-classification',
                 scale = TRUE, 
                 cost = C.opt[3],
                 gamma = gamma.opt,
                 kernel = 'radial')

# presnost na trenovacich datech
predictions.train.l <- predict(clf.SVM.l, newdata = grid.data)
presnost.train.l <- table(Y.train, predictions.train.l) |>
  prop.table() |> diag() |> sum()

predictions.train.p <- predict(clf.SVM.p, newdata = grid.data)
presnost.train.p <- table(Y.train, predictions.train.p) |>
  prop.table() |> diag() |> sum()

predictions.train.r <- predict(clf.SVM.r, newdata = grid.data)
presnost.train.r <- table(Y.train, predictions.train.r) |>
  prop.table() |> diag() |> sum()

# presnost na testovacich datech
predictions.test.l <- predict(clf.SVM.l, newdata = grid.data.test)
presnost.test.l <- table(Y.test, predictions.test.l) |>
  prop.table() |> diag() |> sum()

predictions.test.p <- predict(clf.SVM.p, newdata = grid.data.test)
presnost.test.p <- table(Y.test, predictions.test.p) |>
  prop.table() |> diag() |> sum()

predictions.test.r <- predict(clf.SVM.r, newdata = grid.data.test)
presnost.test.r <- table(Y.test, predictions.test.r) |>
  prop.table() |> diag() |> sum()
```

Přesnost metody SVM na trénovacích datech je tedy 99.3333 % pro lineární jádro, 98.6667 % pro polynomiální jádro a 98.6667 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 92.3077 % pro lineární jádro, 92.3077 % pro polynomiální jádro a 95.3846 % pro radiální jádro.


```r
Res <- data.frame(model = c('SVM linear - diskr', 
                            'SVM poly - diskr', 
                            'SVM rbf - diskr'), 
                  Err.train = 1 - c(presnost.train.l, presnost.train.p, presnost.train.r),
                  Err.test = 1 - c(presnost.test.l, presnost.test.p, presnost.test.r))

RESULTS <- rbind(RESULTS, Res)
```

#### Skóre hlavních komponent {#PCA-SVMA3}

V tomto případě využijeme skóre prvních $p =$ 2 hlavních komponent.

Nyní se pokusme, na rozdíl od postupu v předchozích kapitolách, hyperparametry klasifikátorů odhadnout z dat pomocí 10-násobné cross-validace. Vzhledem k tomu, že každé jádro má ve své definici jiné hyperparametry, budeme ke každé jádrové funkci přistupovat zvlášť. Nicméně hyperparametr $C$ vystupuje u všech jádrových funkcí, přičemž ale připouštíme, že se může jeho optimální hodnota mezi jádry lišit.

U všech třech jader projdeme hodnoty hyperparametru $C$ v intervalu $[10^{-3}, 10^{3}]$, přičemž u jádra polynomiálního zafixujeme hyperparametr $p$ na hodnotě 3, neboť pro jiné celočíselné hodnoty metoda nedává zdaleka tak dobré výsledky. Naopak pro radiální jádro využijeme k volbě optimální hodnoty hyperparametru $\gamma$ opět 10-násobnou CV, přičemž uvažujeme hodnoty v intervalu $[10^{-3}, 10^{2}]$. Zvolíme `coef0` $= 1$. 


```r
set.seed(42)

# ktere hodnoty gamma chceme uvazovat
gamma.cv <- 10^seq(-3, 2, length = 15)
C.cv <- 10^seq(-3, 3, length = 20)
p.cv <- 3
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

Podívejme se, jak dopadly optimální hodnoty. Pro *lineární jádro* máme optimální hodnotu $C$ rovnu 0.1624, pro *polynomiální jádro* je $C$ rovno 0.0785 a pro *radiální jádro* máme dvě optimální hodnoty, pro $C$ je optimální hodnota 26.3665 a pro $\gamma$ je to 1.6379. Validační přesnosti jsou postupně 0.671131 pro lineární, 0.671131 pro polynomiální a 0.7040476 pro radiální jádro.

Konečně můžeme sestrojit finální klasifikátory na celých trénovacích datech s hodnotami hyperparametrů určenými pomocí 10-násobné CV. Určíme také chybovosti na testovacích a také na trénovacích datech.


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

# presnost na trenovacich datech
predictions.train.l <- predict(clf.SVM.l.PCA, newdata = data.PCA.train)
presnost.train.l <- table(data.PCA.train$Y, predictions.train.l) |>
  prop.table() |> diag() |> sum()

predictions.train.p <- predict(clf.SVM.p.PCA, newdata = data.PCA.train)
presnost.train.p <- table(data.PCA.train$Y, predictions.train.p) |>
  prop.table() |> diag() |> sum()

predictions.train.r <- predict(clf.SVM.r.PCA, newdata = data.PCA.train)
presnost.train.r <- table(data.PCA.train$Y, predictions.train.r) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test.l <- predict(clf.SVM.l.PCA, newdata = data.PCA.test)
presnost.test.l <- table(data.PCA.test$Y, predictions.test.l) |>
  prop.table() |> diag() |> sum()

predictions.test.p <- predict(clf.SVM.p.PCA, newdata = data.PCA.test)
presnost.test.p <- table(data.PCA.test$Y, predictions.test.p) |>
  prop.table() |> diag() |> sum()

predictions.test.r <- predict(clf.SVM.r.PCA, newdata = data.PCA.test)
presnost.test.r <- table(data.PCA.test$Y, predictions.test.r) |>
  prop.table() |> diag() |> sum()
```

Přesnost metody SVM aplikované na skóre hlavních komponent na trénovacích datech je tedy 68 % pro lineární jádro, 66.67 % pro polynomiální jádro a 84 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 70.7692 % pro lineární jádro, 72.3077 % pro polynomiální jádro a 63.0769 % pro radiální jádro.

Pro grafické znázornění metody můžeme zaznačit dělící hranici do grafu skórů prvních dvou hlavních komponent.
Tuto hranici spočítáme na husté síti bodů a zobrazíme ji pomocí funkce `geom_contour()` stejně jako v předchozích případech, kdy jsme také vykreslovali klasifikační hranici.


```r
nd <- rbind(nd, nd, nd) |> mutate(
   prd = c(as.numeric(predict(clf.SVM.l.PCA, newdata = nd, type = 'response')),
           as.numeric(predict(clf.SVM.p.PCA, newdata = nd, type = 'response')),
           as.numeric(predict(clf.SVM.r.PCA, newdata = nd, type = 'response'))),
   kernel = rep(c('lineární', 'polynomiální', 'radiální'),
                each = length(as.numeric(predict(clf.SVM.l.PCA, 
                                                 newdata = nd,
                                                 type = 'response')))) |>
     as.factor())

data.PCA.train |> ggplot(aes(x = V1, y = V2, colour = Y)) +
 geom_point(size = 1.5) + 
 labs(x = paste('1. hlavní komponenta (vysvětlená variabilita', 
                round(100 * data.PCA$varprop[1], 2), '%)'),
      y = paste('2. hlavní komponenta (', 
                round(100 * data.PCA$varprop[2], 2), '%)'),
      colour = 'Obsah tuku', 
      linetype = 'Jádro') +
 scale_colour_discrete(labels = c("malý", "velký")) +
 theme_bw() +
 geom_contour(data = nd, aes(x = V1, y = V2, z = prd, linetype = kernel), 
              colour = 'black') 
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-76-1.png" alt="Skóre prvních dvou hlavních komponent, barevně odlišené podle příslušnosti do klasifikační třídy. Černě je vyznačena dělící hranice (přímka, resp. křivky v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí metody SVM." width="672" />
<p class="caption">(\#fig:unnamed-chunk-76)Skóre prvních dvou hlavních komponent, barevně odlišené podle příslušnosti do klasifikační třídy. Černě je vyznačena dělící hranice (přímka, resp. křivky v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí metody SVM.</p>
</div>


```r
Res <- data.frame(model = c('SVM linear - PCA', 
                            'SVM poly - PCA', 
                            'SVM rbf - PCA'), 
                  Err.train = 1 - c(presnost.train.l, presnost.train.p, presnost.train.r),
                  Err.test = 1 - c(presnost.test.l, presnost.test.p, presnost.test.r))

RESULTS <- rbind(RESULTS, Res)
```

#### Bázové koeficienty {#basis-SVMA3}

Nakonec použijeme vyjádření funkcí pomocí B-splinové báze.
U všech třech jader projdeme hodnoty hyperparametru $C$ v intervalu $[10^{-3}, 10^{3}]$, přičemž u jádra polynomiálního zafixujeme hyperparametr $p$ na hodnotě 3, neboť pro jiné celočíselné hodnoty metoda nedává zdaleka tak dobré výsledky. Naopak pro radiální jádro využijeme k volbě optimální hodnoty hyperparametru $\gamma$ opět 10-násobnou CV, přičemž uvažujeme hodnoty v intervalu $[10^{-3}, 10^{2}]$. Zvolíme `coef0` $= 1$.


```r
set.seed(42)

# ktere hodnoty gamma chceme uvazovat
gamma.cv <- 10^seq(-3, 2, length = 15)
C.cv <- 10^seq(-3, 3, length = 20)
p.cv <- 3
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
    cv_sample <- 1:dim(grid.data)[1] %in% fold
    
    data.Bbasis.train.cv <- as.data.frame(data.Bbasis.train[cv_sample, ])
    data.Bbasis.test.cv <- as.data.frame(data.Bbasis.train[!cv_sample, ])
    
    ## LINEARNI JADRO
    # sestrojeni modelu
    clf.SVM.l <- svm(Y ~ ., data = data.Bbasis.train.cv,
                     type = 'C-classification',
                     scale = TRUE,
                     cost = C,
                     kernel = 'linear')
    
    # presnost na validacnich datech
    predictions.test.l <- predict(clf.SVM.l, newdata = data.Bbasis.test.cv)
    presnost.test.l <- table(data.Bbasis.test.cv$Y, predictions.test.l) |>
      prop.table() |> diag() |> sum()
    
    # presnosti vlozime na pozice pro dane C a fold
    CV.results$SVM.l[(1:length(C.cv))[C.cv == C], 
                     index_cv] <- presnost.test.l
    
    ## POLYNOMIALNI JADRO
    for (p in p.cv) {
      # sestrojeni modelu
      clf.SVM.p <- svm(Y ~ ., data = data.Bbasis.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       coef0 = coef0,
                       degree = p,
                       kernel = 'polynomial')
      
      # presnost na validacnich datech
      predictions.test.p <- predict(clf.SVM.p, 
                                    newdata = data.Bbasis.test.cv)
      presnost.test.p <- table(data.Bbasis.test.cv$Y, predictions.test.p) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane C, p a fold
      CV.results$SVM.p[(1:length(C.cv))[C.cv == C], 
                       (1:length(p.cv))[p.cv == p],
                       index_cv] <- presnost.test.p
    }
        
    ## RADIALNI JADRO
    for (gamma in gamma.cv) {
      # sestrojeni modelu
      clf.SVM.r <- svm(Y ~ ., data = data.Bbasis.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       gamma = gamma,
                       kernel = 'radial')
      
      # presnost na validacnich datech
      predictions.test.r <- predict(clf.SVM.r, 
                                    newdata = data.Bbasis.test.cv)
      presnost.test.r <- table(data.Bbasis.test.cv$Y, predictions.test.r) |>
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

Podívejme se, jak dopadly optimální hodnoty. Pro *lineární jádro* máme optimální hodnotu $C$ rovnu 6.1585, pro *polynomiální jádro* je $C$ rovno 54.5559 a pro *radiální jádro* máme dvě optimální hodnoty, pro $C$ je optimální hodnota 26.3665 a pro $\gamma$ je to 0.0118. Validační přesnosti jsou postupně 0.9804167 pro lineární, 0.9675 pro polynomiální a 0.9741667 pro radiální jádro.

Konečně můžeme sestrojit finální klasifikátory na celých trénovacích datech s hodnotami hyperparametrů určenými pomocí 10-násobné CV. Určíme také chybovosti na testovacích a také na trénovacích datech.


```r
# sestrojeni modelu
clf.SVM.l.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                        type = 'C-classification',
                        scale = TRUE,
                        cost = C.opt[1],
                        kernel = 'linear')

clf.SVM.p.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                        type = 'C-classification',
                        scale = TRUE,
                        cost = C.opt[2],
                        degree = p.opt,
                        coef0 = coef0,
                        kernel = 'polynomial')

clf.SVM.r.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                        type = 'C-classification',
                        scale = TRUE,
                        cost = C.opt[3],
                        gamma = gamma.opt,
                        kernel = 'radial')

# presnost na trenovacich datech
predictions.train.l <- predict(clf.SVM.l.Bbasis, newdata = data.Bbasis.train)
presnost.train.l <- table(Y.train, predictions.train.l) |>
  prop.table() |> diag() |> sum()

predictions.train.p <- predict(clf.SVM.p.Bbasis, newdata = data.Bbasis.train)
presnost.train.p <- table(Y.train, predictions.train.p) |>
  prop.table() |> diag() |> sum()

predictions.train.r <- predict(clf.SVM.r.Bbasis, newdata = data.Bbasis.train)
presnost.train.r <- table(Y.train, predictions.train.r) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test.l <- predict(clf.SVM.l.Bbasis, newdata = data.Bbasis.test)
presnost.test.l <- table(Y.test, predictions.test.l) |>
  prop.table() |> diag() |> sum()

predictions.test.p <- predict(clf.SVM.p.Bbasis, newdata = data.Bbasis.test)
presnost.test.p <- table(Y.test, predictions.test.p) |>
  prop.table() |> diag() |> sum()

predictions.test.r <- predict(clf.SVM.r.Bbasis, newdata = data.Bbasis.test)
presnost.test.r <- table(Y.test, predictions.test.r) |>
  prop.table() |> diag() |> sum()
```

Přesnost metody SVM aplikované na bázové koeficienty na trénovacích datech je tedy 99.33 % pro lineární jádro, 99.33 % pro polynomiální jádro a 98.67 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 93.8462 % pro lineární jádro, 90.7692 % pro polynomiální jádro a 93.8462 % pro radiální jádro.


```r
Res <- data.frame(model = c('SVM linear - Bbasis', 
                            'SVM poly - Bbasis', 
                            'SVM rbf - Bbasis'), 
                  Err.train = 1 - c(presnost.train.l, presnost.train.p, presnost.train.r),
                  Err.test = 1 - c(presnost.test.l, presnost.test.p, presnost.test.r))

RESULTS <- rbind(RESULTS, Res)
```

#### Projekce na B-splinovou bázi {#projection-SVMA3}

Další možností, jak použít klasickou metodu SVM pro funkcionální data, je projektovat původní data na nějaký $d$-dimenzionální podprostor našeho Hilbertova prostoru $\mathcal H$, označme jej $V_d$.
Předpokládejme, že tento podprostor $V_d$ má ortonormální bázi $\{\Psi_j\}_{j = 1, \dots, d}$.
Definujeme transformaci $P_{V_d}$ jakožto ortogonální projekci na podprostor $V_d$, tedy můžeme psát

$$
P_{V_d} (x) = \sum_{j = 1}^d \langle x, \Psi_j \rangle \Psi_j.
$$

Nyní můžeme pro klasifikaci použít koeficienty z ortogonální projekce, tedy aplikujeme standardní SVM na vektory $\left( \langle x, \Psi_1 \rangle, \dots, \langle x, \Psi_d \rangle\right)^\top$.
Využitím této transformace jsme tedy definovali nové, tzv.
adaptované jádro, které je složené z ortogonální projekce $P_{V_d}$ a jádrové funkce standardní metody podpůrných vektorů.
Máme tedy (adaptované) jádro $Q(x_i, x_j) = K(P_{V_d}(x_i), P_{V_d}(x_j))$.
Jde tedy o metodu redukce dimenze, kterou můžeme nazvat *filtrace*.

Pro samotnou projekci použijeme v `R` funkci `project.basis()` z knihovny `fda`.
Na jejím vstupu bude matice původních diskrétních (nevyhlazených) dat, hodnoty, ve kterých měříme hodnoty v matici původních dat a bázový objekt, na který chceme data projektovat.
My zvolíme projekci na B-splinovou bázi, protože využití Fourierovy báze není pro naše neperiodická data vhodné.

Dimenzi $d$ volíme buď z nějaké předchozí expertní znalosti, nebo pomocí cross-validace.
V našem případě určíme optimální dimenzi podprostoru $V_d$ pomocí $k$-násobné cross-validace (volíme $k \ll n$ kvůli výpočetní náročnosti metody, často se volí $k = 5$ nebo $k = 10$).
Požadujeme B-spliny řádu 4, pro počet bázových funkcí potom platí vztah

$$
n_{basis} = n_{breaks} + n_{order} - 2,
$$

kde $n_{breaks}$ je počet uzlů a $n_{order} = 4$.
Minimální dimenzi tedy (pro $n_{breaks} = 1$) volíme $n_{basis} = 3$ a maximální (pro $n_{breaks} = 51$ odpovídající počtu původních diskrétních dat) $n_{basis} = 53$.
V `R` však hodnota $n_{basis}$ musí být alespoň $n_{order} = 4$ a pro velké hodnoty $n_{basis}$ již dochází k přefitování modelu, tudíž volíme za maximální $n_{basis}$ menší číslo, řekněme 20.


```r
k_cv <- 10 # k-fold CV

# hodnoty pro B-splinovou bazi
rangeval <- range(t)
norder <- 4
n_basis_min <- norder
n_basis_max <- 20

dimensions <- n_basis_min:n_basis_max # vsechny dimenze, ktere chceme vyzkouset

# rozdelime trenovaci data na k casti
folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)

# list se tremi slozkami ... maticemi pro jednotlive jadra -> linear, poly, radial
# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro danou cast trenovaci mnoziny
# v radcich budou hodnoty pro danou hodnotu dimenze
CV.results <- list(SVM.l = matrix(NA, nrow = length(dimensions), ncol = k_cv),
                   SVM.p = matrix(NA, nrow = length(dimensions), ncol = k_cv),
                   SVM.r = matrix(NA, nrow = length(dimensions), ncol = k_cv))

for (d in dimensions) {
  # bazovy objekt
  bbasis <- create.bspline.basis(rangeval = rangeval, 
                                 nbasis = d)
  
  # projekce diskretnich dat na B-splinovou bazi o dimenzi d
  Projection <- project.basis(y = XX, # matice diskretnich dat
                              argvals = t, # vektor argumentu
                              basisobj = bbasis) # bazovy objekt
  
  # rozdeleni na trenovaci a testovaci data v ramci CV
  XX.train <- subset(t(Projection), split == TRUE)
  
  for (index_cv in 1:k_cv) {
    # definice testovaci a trenovaci casti pro CV
    fold <- folds[[index_cv]]
    cv_sample <- 1:dim(XX.train)[1] %in% fold
    
    data.projection.train.cv <- as.data.frame(XX.train[cv_sample, ])
    data.projection.train.cv$Y <- factor(Y.train[cv_sample])
    
    data.projection.test.cv <- as.data.frame(XX.train[!cv_sample, ])
    Y.test.cv <- Y.train[!cv_sample]
    data.projection.test.cv$Y <- factor(Y.test.cv)
  
    # sestrojeni modelu
    clf.SVM.l.projection <- svm(Y ~ ., data = data.projection.train.cv,
                            type = 'C-classification',
                            scale = TRUE,
                            kernel = 'linear')
    
    clf.SVM.p.projection <- svm(Y ~ ., data = data.projection.train.cv,
                            type = 'C-classification',
                            scale = TRUE,
                            coef0 = coef0,
                            kernel = 'polynomial')
    
    clf.SVM.r.projection <- svm(Y ~ ., data = data.projection.train.cv,
                            type = 'C-classification',
                            scale = TRUE,
                            kernel = 'radial')
      
    # presnost na validacnich datech
    ## linear kernel
    predictions.test.l <- predict(clf.SVM.l.projection,
                                  newdata = data.projection.test.cv)
    presnost.test.l <- table(Y.test.cv, predictions.test.l) |>
      prop.table() |> diag() |> sum()
    ## polynomial kernel
    predictions.test.p <- predict(clf.SVM.p.projection, 
                                  newdata = data.projection.test.cv)
    presnost.test.p <- table(Y.test.cv, predictions.test.p) |>
      prop.table() |> diag() |> sum()
    ## radial kernel
    predictions.test.r <- predict(clf.SVM.r.projection,
                                  newdata = data.projection.test.cv)
    presnost.test.r <- table(Y.test.cv, predictions.test.r) |>
      prop.table() |> diag() |> sum()
    
    # presnosti vlozime na pozice pro dane d a fold
    CV.results$SVM.l[d - min(dimensions) + 1, index_cv] <- presnost.test.l
    CV.results$SVM.p[d - min(dimensions) + 1, index_cv] <- presnost.test.p
    CV.results$SVM.r[d - min(dimensions) + 1, index_cv] <- presnost.test.r
  }
}
  
# spocitame prumerne presnosti pro jednotliva d pres folds
for (n_method in 1:length(CV.results)) {
  CV.results[[n_method]] <- apply(CV.results[[n_method]], 1, mean)
}

d.opt <- c(which.max(CV.results$SVM.l) + n_basis_min - 1, 
           which.max(CV.results$SVM.p) + n_basis_min - 1, 
           which.max(CV.results$SVM.r) + n_basis_min - 1)
presnost.opt.cv <- c(max(CV.results$SVM.l), 
                     max(CV.results$SVM.p),
                     max(CV.results$SVM.r))
data.frame(d_opt = d.opt, ERR = 1 - presnost.opt.cv,
           row.names = c('linear', 'poly', 'radial'))
```

```
##        d_opt        ERR
## linear     9 0.01958333
## poly       7 0.03297619
## radial     6 0.14125000
```

Vidíme, že nejlépe vychází hodnota parametru $d$ jako 9 pro lineární jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.9804, 7 pro polynomiální jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.967 a 6 pro radiální jádro s hodnotou přesnosti 0.8588.
Pro přehlednost si ještě vykresleme průběh validačních chybovostí v závislosti na dimenzi $d$.


```r
CV.results <- data.frame(d = dimensions |> rep(3), 
                         CV = c(CV.results$SVM.l, 
                                CV.results$SVM.p, 
                                CV.results$SVM.r),
                         Kernel = rep(c('lineární', 'polynomiální', 'radiální'), 
                                      each = length(dimensions)) |> factor())
CV.results |> ggplot(aes(x = d, y = 1 - CV, colour = Kernel)) + 
  geom_line(linetype = 'dashed') + 
  geom_point(size = 1.5) + 
  geom_point(data = data.frame(d.opt,
                               presnost.opt.cv),
             aes(x = d.opt, y = 1 - presnost.opt.cv), colour = 'black', size = 2) +
  theme_bw() + 
  labs(x = bquote(paste(d)),
       y = 'Validační chybovost',
       colour = 'Jádro') + 
  theme(legend.position = "bottom") + 
  scale_x_continuous(breaks = dimensions)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-83-1.png" alt="Závislost validační chybovosti na dimenzi podprostoru $V_d$, zvlášť pro všechna tři uvažovaná jádra v metodě SVM. Černými body jsou vyznačeny optimální hodnoty dimenze $V_d$ pro jednotlivé jádrové funkce." width="672" />
<p class="caption">(\#fig:unnamed-chunk-83)Závislost validační chybovosti na dimenzi podprostoru $V_d$, zvlášť pro všechna tři uvažovaná jádra v metodě SVM. Černými body jsou vyznačeny optimální hodnoty dimenze $V_d$ pro jednotlivé jádrové funkce.</p>
</div>

Nyní již můžeme natrénovat jednotlivé klasifikátory na všech trénovacích datech a podívat se na jejich úspěšnost na testovacích datech.
Pro každou jádrovou funkci volíme dimenzi podprostoru, na který projektujeme, podle výsledků cross-validace.

V proměnné `Projection` máme uloženou matici koeficientů ortogonální projekce, tedy

$$
\texttt{Projection} = \begin{pmatrix}
\langle x_1, \Psi_1 \rangle & \langle x_2, \Psi_1 \rangle & \cdots & \langle x_n, \Psi_1 \rangle\\
\langle x_1, \Psi_2 \rangle & \langle x_2, \Psi_2 \rangle & \cdots & \langle x_n, \Psi_2 \rangle\\
\vdots & \vdots & \ddots & \vdots \\
\langle x_1, \Psi_d \rangle & \langle x_2, \Psi_d \rangle & \dots & \langle x_n, \Psi_d \rangle
\end{pmatrix}_{d \times n}.
$$


```r
# pripravime si datovou tabulku pro ulozeni vysledku
Res <- data.frame(model = c('SVM linear - projection', 
                            'SVM poly - projection', 
                            'SVM rbf - projection'), 
                  Err.train = NA,
                  Err.test = NA)

# projdeme jednotliva jadra
for (kernel_number in 1:3) {
  kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
  # bazovy objekt
  bbasis <- create.bspline.basis(rangeval = rangeval, 
                                 nbasis = d.opt[kernel_number])
  
  # projekce diskretnich dat na B-splinovou bazi
  Projection <- project.basis(y = XX, # matice diskretnich dat
                              argvals = t, # vektor argumentu
                              basisobj = bbasis) # bazovy objekt
  
  # rozdeleni na trenovaci a testovaci data
  XX.train <- subset(t(Projection), split == TRUE)
  XX.test <- subset(t(Projection), split == FALSE)
  
  data.projection.train <- as.data.frame(XX.train)
  data.projection.train$Y <- factor(Y.train)
  
  data.projection.test <- as.data.frame(XX.test)
  data.projection.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.SVM.projection <- svm(Y ~ ., data = data.projection.train,
                            type = 'C-classification',
                            scale = TRUE,
                            coef0 = coef0,
                            kernel = kernel_type)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.SVM.projection, newdata = data.projection.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test <- predict(clf.SVM.projection, newdata = data.projection.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  # ulozeni vysledku
  Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
}
```

Přesnost metody SVM aplikované na bázové koeficienty na trénovacích datech je tedy 2 % pro lineární jádro, 2.67 % pro polynomiální jádro a 9.33 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 6.15 % pro lineární jádro, 6.15 % pro polynomiální jádro a 10.77 % pro radiální jádro.


```r
RESULTS <- rbind(RESULTS, Res)
```

#### RKHS + SVM {#RKHS-SVMA3} 

V této sekci se podíváme na další možnost, jak využít metodu podpůrných vektorů pro klasifikaci funkcionálních dat.
V tomto případě půjde opět o již nám známý princip, kdy nejprve funkcionální data vyjádříme jakožto nějaké konečně-rozměrné objekty a na tyto objekty následně aplikujeme klasickou metodu SVM.

Nyní však metodu SVM použijeme i pro samotnou reprezentaci funkcionálních dat pomocí určitého konečně-rozměrného objektu.
Jak již název napovídá, půjde o kombinaci dvou konceptů -- jednak metody podpůrných vektorů a druhak prostoru, který se nazývá v anglické literatuře *Reproducing Kernel Hilbert Space*.
Pro tento prostor je klíčovým pojmem *jádro* -- *kernel*.

::: {.definition #kernel name="Jádro"} 
Jádro je taková funkce $K : \mathcal X \times \mathcal X \rightarrow \mathbb R$, že pro každou dvojici $\boldsymbol x, \tilde{\boldsymbol x} \in \mathcal X$ platí 
\begin{equation*}
K(\boldsymbol x, \tilde{\boldsymbol x}) = \big\langle \boldsymbol\phi(\boldsymbol x), \boldsymbol\phi(\tilde{\boldsymbol x}) \big\rangle_{\mathcal H},
\end{equation*}
kde $\boldsymbol\phi : \mathcal X \rightarrow \mathcal H$ je zobrazení z prostoru $\mathcal X$ do prostoru $\mathcal H$.
:::

Aby funkce byla jádrem, musí splňovat určité podmínky.

::: {.lemma} 
Nechť $\mathcal X$ je nějaký Hilbertův prostor. Potom symetrická funkce $K : \mathcal X \times \mathcal X \rightarrow \mathbb R$ je jádrem, pokud $\forall k \geq 1, \boldsymbol x_1, \dots, \boldsymbol x_k \in \mathcal X$ a $c_1, \dots, c_k \in \mathbb R$ platí
\begin{equation*}
\sum_{i, j = 1}^k c_ic_j K(\boldsymbol x_i, \boldsymbol x_j) \geq 0.
\end{equation*}
:::

Vlastnost výše se nazývá pozitivní semidefinitnost.
Platí také následující tvrzení.

::: {.theorem} 
Funkce $K: \mathcal X \times \mathcal X \rightarrow \mathbb R$ je jádrem právě tehdy, když existuje Hilbertův prostor $\mathcal H$ a zobrazení $\boldsymbol\phi : \mathcal X \rightarrow \mathcal H$ takové, že
\begin{equation*}
K(\boldsymbol x, \tilde{\boldsymbol x}) = \big\langle \boldsymbol\phi(\boldsymbol x), \boldsymbol\phi(\tilde{\boldsymbol x}) \big\rangle_{\mathcal H} \quad \forall \boldsymbol x, \tilde{\boldsymbol x}\in \mathcal X.
\end{equation*}
:::

Nyní již máme připravenou půdu pro zavedení pojmu *Reproducing Kernel Hilbert Space*.

##### Reproducing Kernel Hilbert Space (RKHS)

Uvažujme Hilbertův prostor $\mathcal H$ jakožto prostor funkcí.
Naším cílem je definovat prostor $\mathcal H$ a zobrazení $\phi$ takové, že $\phi(x) \in \mathcal H, \ \forall x \in \mathcal X$.
Označme $\phi(x) = k_x$.
Každé funkci $x \in \mathcal X$ tedy přiřadíme funkci $x \mapsto k_x \in \mathcal H, k_x := K(x, \cdot), k_x: \mathcal X \rightarrow \mathbb R$.
Potom $\phi: \mathcal X \rightarrow \mathbb R^{\mathcal X}$, můžeme tedy souhrnně napsat

$$
x \in \mathcal X \mapsto \phi(x) = k_x = K(x, \cdot) \in \mathcal H,
$$

Bod (funkce) $x \in \mathcal X$ je zobrazen na funkci $k_x: \mathcal X \rightarrow \mathbb R, k_x(y) = K(x, y)$.

Uvažujme množinu všech obrazů $\{k_x | x \in \mathcal X\}$ a definujme lineární obal této množiny vektorů jakožto

$$
\mathcal G := \text{span}\{k_x | x \in \mathcal X\} = \left\{\sum_{i = 1}^r\alpha_i K(x_i, \cdot)\ \Big|\ \alpha_i \in \mathbb R, r \in \mathbb N, x_i \in \mathcal X\right\}.
$$

Potom skalární součin

$$
\langle k_x, k_y \rangle = \langle K(x, \cdot), K(y, \cdot) \rangle = K(x, y),\quad x, y \in \mathcal X
$$

a obecně

$$
f, g \in \mathcal G, f = \sum_i \alpha_i K(x_i, \cdot), g = \sum_j \beta_j K(y_j, \cdot), \\
\langle f, g \rangle_{\mathcal G} = \Big\langle \sum_i \alpha_i K(x_i, \cdot), \sum_j \beta_j K(y_j, \cdot) \Big\rangle = \sum_i\sum_j\alpha_i\beta_j \langle K(x_i, \cdot), K(y_j, \cdot) \rangle = \sum_i\sum_j\alpha_i\beta_j K(x_i, y_j).
$$

Prostor $\mathcal H := \overline{\mathcal G}$, který je zúplněním prostoru $\mathcal G$, nazýváme *Reproducing Kernel Hilbert Space* (RKHS).
Významnou vlastností tohoto prostoru je

$$
K(x, y) = \Big\langle \phi(x), \phi(y) \Big\rangle_{\mathcal H}.
$$

*Poznámka:* Jméno Reproducing vychází z následujícího faktu.
Mějme libovolnou funkci $f = \sum_i \alpha_i K(x_i, \cdot)$.
Potom

```{=tex}
\begin{align*}
\langle K(x, \cdot), f\rangle &= \langle K(x, \cdot), \sum_i \alpha_i K(x_i, \cdot) \rangle =\\
&= \sum_i \alpha_i \langle K(x, \cdot), K(x_i, \cdot) \rangle = \sum_i \alpha_i K(x_i, x) = \\
&= f(x)
\end{align*}
```

*Vlastnosti:*

-   nechť $\mathcal H$ je Hilbertův prostor funkcí $g: \mathcal X \rightarrow \mathbb R$.
    Potom $\mathcal H$ je RKHS $\Leftrightarrow$ všechny funkcionály (evaluation functionals) $\delta_x: \mathcal H \rightarrow \mathbb R, g \mapsto g(x)$ jsou spojité,

-   pro dané jádro $K$ existuje právě jeden prostor RKHS (až na isometrickou izomofrii),

-   pro daný RKHS je jádro $K$ určeno jednoznačně,

-   funkce v RKHS jsou bodově korektně definovány,

-   RKHS je obecně nekonečně-rozměrný vektorový prostor, v praxi však pracujeme pouze s jeho konečně-rozměrným podprostorem.

Na konec této sekce si uveďme jedno důležité tvrzení.

::: {.theorem #representer name="The representer theorem"} 
Nechť $K$ je jádro a $\mathcal H$ je příslušný RKHS s normou a skalárním součinem $\|\cdot\|_{\mathcal H}$ a $\langle \cdot, \cdot \rangle_{\mathcal H}$. Předpokládejme, že chceme zjistit lineární funkci $f: \mathcal H \rightarrow \mathbb R$ na Hilbertově prostoru $\mathcal H$ definovaného jádrem $K$. Funkce $f$ má tvar $f(x) = \langle \omega, x \rangle_{\mathcal H}$ pro nějaké $\omega \in \mathcal H$. Uvažujme regularizovaný minimalizační problém 
\begin{equation}
\min_{\omega \in \mathcal H} R_n(\omega) + \lambda \Omega(\|\omega\|_{\mathcal H}),
(\#eq:repr-thm)
\end{equation}
kde $\Omega: [0, \infty) \rightarrow \mathbb R$ je striktně monotonně rostoucí funkce (regularizer), $R_n(\cdot)$ je empirická ztráta (empirical risk) klasifikátoru vzhledem ke ztrátové funkci $\ell$. Potom optimalizační úloha \@ref(eq:repr-thm) má vždy optimální řešení a to je tvaru 
\begin{equation}
\omega^* = \sum_{i = 1}^n \alpha_i K(x_i, \cdot),
(\#eq:optimal-thm)
\end{equation}
kde $(x_i, y_i)_{i = 1, 2, \dots, n} \in \mathcal X \times \mathcal Y$ je množina trénovacích hodnot.
:::

$\mathcal H$ je obecně nekočně-rozměrný prostor, ale pro konečný datový soubor velikosti $n$ má $\mathcal H$ dimenzi nejvýše $n$.
Každý $n$-dimenzionální podprostor Hilbertova prostoru je navíc izometrický s $\mathbb R^n$, tudíž můžeme předpokládat, že zobrazení (feature map) zobrazuje právě do $\mathbb R^n$.

Jádro $K$ je *univerzální* pokud RKHS $\mathcal H$ je hustá množina v $\mathcal C(\mathcal X)$ (množina spojitých funkcí).
Navíc platí následující poznatky:

-   univerzální jádra jsou dobrá pro aproximaci,
-   Gaussovo jádro s pevnou hodnotou $\sigma$ je univerzální,
-   univerzalita je nutnou podmínkou pro konzistenci.

##### Klasifikace pomocí RKHS

Základní myšlenkou je projekce původních dat na podprostor prostoru RKHS, označme jej $\mathcal H_K$ (index ${}_K$ odkazuje na fakt, že tento prostor je definován jádrem $K$).
Cílem je tedy transformovat křivku (pozorovaný objekt, funkce) na bod v RKHS.
Označme $\{\hat c_1, \dots, \hat c_n\}$ množinu pozorovaných křivek, přičemž každá křivka $\hat c_l$ je definována daty $\{(\boldsymbol x_i, \boldsymbol y_{il}) \in \mathcal X \times \mathcal Y\}_{i = 1}^m$, kde $\mathcal X$ je prostor vstupních proměnných a nejčastěji $\mathcal Y = \mathbb R$.
Předpokládejme, že pro každou funkci $\hat c_l$ existuje spojitá funkce $c_l:\mathcal X \rightarrow \mathcal Y, \mathbb E[y_l|\boldsymbol x] = c_l(\boldsymbol x)$.
Předpokládejme také, že $\boldsymbol x_i$ jsou společné pro všechny křivky.

Muñoz a González ve svém článku[^1] navrhují následující postup.
Křivku $c_l^*$ můžeme napsat ve tvaru

[^1]: Muñoz, A.
    and González, J.
    (2010) *Representing functional data using support vector machines*, Pattern Recognition Letters, 31(6), pp. 511--516.
    [doi:10.1016/j.patrec.2009.07.014](https://www.sciencedirect.com/science/article/pii/S0167865509001913).

$$
c_l^*(\boldsymbol x) = \sum_{i = 1}^m \alpha_{il} K(\boldsymbol x_i, \boldsymbol x), \quad \forall \boldsymbol x \in \mathcal X,
$$

kde $\alpha_{il} \in \mathbb R$.
Tyto koeficienty získáme v praxi řešením optimalizačního problému $$
\text{argmin}_{c \in \mathcal H_K} \frac{1}{m} \sum_{i = 1}^m \big[|c(\boldsymbol x_i) - y_i| - \varepsilon\big]_+ + \gamma \|c\|_{K}^2, \gamma > 0, \varepsilon \geq 0,
$$ tedy právě například pomocí metody SVM.
Díky známé vlastnosti této metody pak bude mnoho koeficientů $\alpha_{il} = 0$.
Minimalizací výše uvedeného výrazu získáme funkce $c_1^*, \dots, c_n^*$ odpovídající původním křivkám $\hat c_1, \dots, \hat c_n$.
Metoda SVM tedy dává smysluplnou reprezentaci původních křivek pomocí vektoru koeficientů $\boldsymbol \alpha_l = (\alpha_{1l}, \dots, \alpha_{ml})^\top$ pro $\hat c_l$.
Tato reprezentace je však velmi nestabilní, neboť i při malé změně původních hodnot může dojít ke změně v množině podpůrných vektorů pro danou funkci, a tedy dojde k výrazné změně celé reprezentace této křivky (reprezentace není spojitá ve vstupních hodnotách).
Definujeme proto novou reprezentaci původních křivek, která již nebude trpět tímto nedostatkem.

::: {.theorem #MaG name="Muñoz-González"} 
Nechť $c$ je funkce, jejíž pozorovaná verze je $\hat c = \{(\boldsymbol x_i, y_{i}) \in \mathcal X \times \mathcal Y\}_{i = 1}^m$ a $K$ je jádro s vlastními funkcemi $\{\phi_1, \dots, \phi_d, \dots\}$ (báze $\mathcal H_K$). Potom funkce $c^*(\boldsymbol x)$ může být vyjádřena ve tvaru
\begin{equation*}
c^*(\boldsymbol x) = \sum_{j = 1}^d \lambda_j^* \phi_j(\boldsymbol x),
\end{equation*}
kde $\lambda_j^*$ jsou váhy projekce $c^*(\boldsymbol x)$ na prostor funkcí generovaný vlastními funkcemi jádra $K$ a $d$ je dimenze prostoru $\mathcal H$. V praxi, kdy máme k dispozici pouze konečně mnoho pozorování, $\lambda_j^*$ mohou být odhadnuty pomocí
\begin{equation*}
\hat\lambda_j^* = \hat\lambda_j \sum_{i = 1}^m \alpha_i\hat\phi_{ji}, \quad j = 1, 2, \dots, \hat d,
\end{equation*}
kde $\hat\lambda_j$ je $j$-té vlastní číslo příslušné $j$-tému vlastnímu vektoru $\hat\phi_j$ matice $K_S = \big(K(\boldsymbol x_i, \boldsymbol x_j)\big)_{i, j = 1}^m, \hat d = \text{rank}(K_S)$ a $\alpha_i$ jsou řešením optimalizačního problému.
:::

##### Implementace metody v `R`

Z poslední části Tvrzení \@ref(thm:MaG) vyplývá, jak máme spočítat v praxi reprezentace křivek.
Budeme pracovat s diskretizovanými daty po vyhlazení křivek.
Nejprve si definujeme jádro pro prostor RKHS.
Využijeme Gaussovské jádro s parametrem $\gamma$.
Hodnota tohoto hyperparametru výrazně ovlivňuje chování a tedy i úspěšnost metody, proto jeho volbě musíme věnovat zvláštní pozornost (volíme pomocí cross-validace).

Jako dobrá volba hyperparametrů se po vyzkoušení zdají být hodnoty $\varepsilon = 0.01$ a $C = 1$. Vzhledem k výpočetní náročnosti nebudeme tyto hyperparametry odhadovat pomocí CV.


```r
eps <- 0.01
C <- 1 
```

###### Gaussovké jádro


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())

# jadro a jadrova matice ... Gaussovske s parametrem gamma
Gauss.kernel <- function(x, y, gamma) {
  return(exp(-gamma * norm(c(x - y) |> t(), type = 'F')))
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

Spočítejme nyní matici $K_S$ a její vlastní čísla a příslušné vlastní vektory.


```r
# spocitame matici K
gamma <- 0.1 # pevna hodnota gamma, optimalni urcime pomoci CV
K <- Kernel.RKHS(t.seq, gamma = gamma)

# urcime vlastni cisla a vektory
Eig <- eigen(K)
eig.vals <- Eig$values
eig.vectors <- Eig$vectors
```

K výpočtu koeficientů v reprezentaci křivek, tedy výpočtu vektorů $\hat{\boldsymbol \lambda}_l^* = \left( \hat\lambda_{1l}^*, \dots, \hat\lambda_{\hat dl}^*\right)^\top, l = 1, 2, \dots, n$, potřebujeme ještě koeficienty z SVM.
Narozdíl od klasifikačního problému nyní řešíme problém regrese, neboť se snažíme vyjádřit naše pozorované křivky v nějaké (námi zvolené pomocí jádra $K$) bázi.
Proto využijeme metodu *Support Vector Regression*, z níž následně získáme koeficienty $\alpha_{il}$.


```r
# urceni koeficientu alpha z SVM
alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                     ncol = dim(data.RKHS)[2]) # prazdny objekt

# model
for(i in 1:dim(data.RKHS)[2]) {
  df.svm <- data.frame(x = t.seq,
                       y = data.RKHS[, i])
  svm.RKHS <- svm(y ~ x, data = df.svm, 
                  kernel = 'radial',
                  type = 'eps-regression',
                  epsilon = eps,
                  cost = C,
                  gamma = gamma)
  # urceni alpha
  alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
}
```

Nyní již můžeme spočítat reprezentace jednotlivých křivek.
Nejprve zvolme za $\hat d$ celou dimenzi, tedy $\hat d = m ={}$ 101, následně určíme optimální $\hat d$ pomocí cross-validace.


```r
# d
d.RKHS <- dim(alpha.RKHS)[1]

# urceni vektoru lambda
Lambda.RKHS <- matrix(NA, 
                      ncol = dim(data.RKHS)[2], 
                      nrow = d.RKHS) # vytvoreni prazdneho objektu

# vypocet reprezentace
for(l in 1:dim(data.RKHS)[2]) {
  Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
}
```

Nyní máme v matici `Lambda.RKHS` uloženy ve sloupcích vektory $\hat{\boldsymbol \lambda}_l^*, l = 1, 2, \dots, n$ pro každou křivku.
Tyto vektory nyní využijeme jakožto reprezentaci daných křivek a klasifikujeme data podle této diskretizace.


```r
# rozdeleni na trenovaci a testovaci data
XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]

# pripravime si datovou tabulku pro ulozeni vysledku
Res <- data.frame(model = c('SVM linear - RKHS', 
                            'SVM poly - RKHS', 
                            'SVM rbf - RKHS'), 
                  Err.train = NA,
                  Err.test = NA)

# projdeme jednotliva jadra
for (kernel_number in 1:3) {
  kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]

  data.RKHS.train <- as.data.frame(t(XX.train))
  data.RKHS.train$Y <- factor(Y.train)
  
  data.RKHS.test <- as.data.frame(t(XX.test))
  data.RKHS.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                      type = 'C-classification',
                      cost = C,
                      coef0 = coef0,
                      scale = TRUE,
                      kernel = kernel_type)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  # ulozeni vysledku
  Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
}
```


Table: (\#tab:unnamed-chunk-92)Souhrnné výsledky metody SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

Model                                $\widehat{Err}_{train}\quad\quad\quad\quad\quad$         $\widehat{Err}_{test}\quad\quad\quad\quad\quad$       
----------------------------------  -------------------------------------------------------  -------------------------------------------------------
SVM linear - RKHS                                                                    0.0133                                                   0.1231
SVM poly - RKHS                                                                      0.0267                                                   0.0308
SVM rbf - RKHS                                                                       0.0400                                                   0.0308

Vidíme, že model u všech třech jader velmi dobře klasifikuje trénovací data, zatímco jeho úspěšnost na testovacích datech není vůbec dobrá.
Je zřejmé, že došlo k overfittingu, proto využijeme cross-validaci, abychom určili optimální hodnoty $\gamma$ a $d$.


```r
# rozdelime trenovaci data na k casti
folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()

# hodnoty hyperparametru, ktere budeme prochazet
dimensions <- 3:30 # rozumny rozsah hodnot d
gamma.cv <- 10^seq(-2, 2, length = 15)

# list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro dane
# v radcich budou hodnoty pro danou gamma a vrstvy odpovidaji folds
dim.names <- list(gamma = paste0('gamma:', round(gamma.cv, 3)),
                  d = paste0('d:', dimensions),
                  CV = paste0('cv:', 1:k_cv))

CV.results <- list(
  SVM.l = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.p = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.r = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                dimnames = dim.names))
```


```r
# samotna CV
for (gamma in gamma.cv) {
  K <- Kernel.RKHS(t.seq, gamma = gamma)
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'radial',
                    type = 'eps-regression',
                    epsilon = eps,
                    cost = C,
                    gamma = gamma)
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
  }
  
  # projdeme dimenze
  for(d.RKHS in dimensions) {
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) 
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                             alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    # projdeme folds
    for (index_cv in 1:k_cv) {
      # definice testovaci a trenovaci casti pro CV
      fold <- folds[[index_cv]]
      # rozdeleni na trenovaci a validacni data
      XX.train <- Lambda.RKHS[, fold]
      XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
      # pripravime si datovou tabulku pro ulozeni vysledku
      Res <- data.frame(model = c('SVM linear - RKHS', 
                                  'SVM poly - RKHS', 
                                  'SVM rbf - RKHS'), 
                        Err.test = NA)
      # projdeme jednotliva jadra
      for (kernel_number in 1:3) {
        kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
      
        data.RKHS.train <- as.data.frame(t(XX.train))
        data.RKHS.train$Y <- factor(Y.train[fold])
        
        data.RKHS.test <- as.data.frame(t(XX.test))
        data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
        
        # sestrojeni modelu
        clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                            type = 'C-classification',
                            scale = TRUE,
                            cost = C,
                            coef0 = coef0,
                            kernel = kernel_type)
        
        # presnost na validacnich datech
        predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
        presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
          prop.table() |> diag() |> sum()
        
        # ulozeni vysledku
        Res[kernel_number, 2] <- 1 - presnost.test
      }
      # presnosti vlozime na pozice pro dane d, gamma a fold
      CV.results$SVM.l[paste0('gamma:', round(gamma, 3)), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[1, 2]
      CV.results$SVM.p[paste0('gamma:', round(gamma, 3)), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[2, 2]
      CV.results$SVM.r[paste0('gamma:', round(gamma, 3)), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[3, 2]
    }
  }
}
```


```r
# spocitame prumerne presnosti pro jednotliva d pres folds
for (n_method in 1:length(CV.results)) {
  CV.results[[n_method]] <- apply(CV.results[[n_method]], c(1, 2), mean)
}

gamma.opt <- c(which.min(CV.results$SVM.l) %% length(gamma.cv), 
               which.min(CV.results$SVM.p) %% length(gamma.cv), 
               which.min(CV.results$SVM.r) %% length(gamma.cv))
gamma.opt[gamma.opt == 0] <- length(gamma.cv)
gamma.opt <- gamma.cv[gamma.opt]

d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
           which.min(t(CV.results$SVM.p)) %% length(dimensions), 
           which.min(t(CV.results$SVM.r)) %% length(dimensions))
d.opt[d.opt == 0] <- length(dimensions)
d.opt <- dimensions[d.opt]

err.opt.cv <- c(min(CV.results$SVM.l), 
                     min(CV.results$SVM.p),
                     min(CV.results$SVM.r))
df.RKHS.res <- data.frame(d = d.opt, gamma = gamma.opt, CV = err.opt.cv,
           Kernel = c('linear', 'polynomial', 'radial') |> factor(),
           row.names = c('linear', 'poly', 'radial'))
```


Table: (\#tab:unnamed-chunk-95)Souhrnné výsledky cross-validace pro metodu SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

          $\quad\quad\quad\quad\quad d$   $\quad\quad\quad\quad\quad\gamma$   $\widehat{Err}_{cross\_validace}$  Model                             
-------  ------------------------------  ----------------------------------  ----------------------------------  ----------------------------------
linear                               11                              1.0000                                   0  linear                            
poly                                 27                              1.0000                                   0  polynomial                        
radial                               30                              3.7276                                   0  radial                            

Vidíme, že nejlépe vychází hodnota parametru $d={}$ 11 a $\gamma={}$ 1 pro lineární jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 1, $d={}$ 27 a $\gamma={}$ 1 pro polynomiální jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 1 a $d={}$ 30 a $\gamma={}$ 3.7276 pro radiální jádro s hodnotou přesnosti 1.
Pro zajímavost si ještě vykresleme funkci validační chybovosti v závislosti na dimenzi $d$ a hodnotě hyperparametru $\gamma$.


```r
CV.results.plot <- data.frame(d = rep(dimensions |> rep(3), each = length(gamma.cv)), 
                              gamma = rep(gamma.cv, length(dimensions)) |> rep(3),
                               CV = c(c(CV.results$SVM.l), 
                                      c(CV.results$SVM.p), 
                                      c(CV.results$SVM.r)),
                               Kernel = rep(c('linear', 'polynomial', 'radial'), 
                                            each = length(dimensions) * 
                                              length(gamma.cv)) |> factor())
CV.results.plot |> 
  ggplot(aes(x = d, y = gamma, z = CV)) + 
  geom_contour_filled() +
  scale_y_continuous(trans='log10') +
  facet_wrap(~Kernel) +
  theme_bw() + 
  labs(x = expression(d),
       y = expression(gamma)) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_point(data = df.RKHS.res, aes(x = d, y = gamma),
             size = 5, pch = '+')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-96-1.png" alt="Závislost validační chybovosti na volbě hyperparametrů $d$ a $\gamma$, zvlášť pro všechna tři uvažovaná jádra v metodě SVM." width="672" />
<p class="caption">(\#fig:unnamed-chunk-96)Závislost validační chybovosti na volbě hyperparametrů $d$ a $\gamma$, zvlášť pro všechna tři uvažovaná jádra v metodě SVM.</p>
</div>

Na grafech výše vidíme, jak se měnila validační chybovost v závislosti na hodnotách hyperparametrů $d$ a $\gamma$.
Všimněme si zejména, že ve všech třech grafech pro jednotlivá jádra jsou patrné výrazné horizontální útvary.
Z toho můžeme usoudit významné teoretické i praktické zjištění -- uvažovaná klasifikační metoda (projekce na RKHS pomocí SVM + klasifikace SVM) je robustní na volbu hyperparametru $d$ (tj. při malé změně v hodnotě tohoto parametru nedojde k výraznému zhoršení validační chybovosti), zatímco při volbě hyperparametru $\gamma$ musíme být velmi obezřetní (i malá změna v jeho hodnotě může vést k velké změně validační chybovosti).
Toto chování je nejlépe patrné u Gaussova jádra.

Jelikož již máme nalezeny optimální hodnoty hyperparametrů, můžeme zkounstruovat finální modely a určit jejich úspěšnost klasifikace na testovacích datech.


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
```


```r
# pripravime si datovou tabulku pro ulozeni vysledku
Res <- data.frame(model = c('SVM linear - RKHS - radial', 
                            'SVM poly - RKHS - radial', 
                            'SVM rbf - RKHS - radial'), 
                  Err.train = NA,
                  Err.test = NA)

# projdeme jednotliva jadra
for (kernel_number in 1:3) {
  # spocitame matici K
  gamma <- gamma.opt[kernel_number] # hodnota gamma pomoci CV
  K <- Kernel.RKHS(t.seq, gamma = gamma)
  
  # urcime vlastni cisla a vektory
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  # urceni koeficientu alpha z SVM
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                       ncol = dim(data.RKHS)[2]) # prazdny objekt
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'radial',
                    type = 'eps-regression',
                    epsilon = eps,
                    cost = C,
                    gamma = gamma)
    # urceni alpha
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
  }
  # d
  d.RKHS <- d.opt[kernel_number]
  
  # urceni vektoru lambda
  Lambda.RKHS <- matrix(NA, 
                        ncol = dim(data.RKHS)[2], 
                        nrow = d.RKHS) # vytvoreni prazdneho objektu
  
  # vypocet reprezentace
  for(l in 1:dim(data.RKHS)[2]) {
    Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
  }
  
  # rozdeleni na trenovaci a testovaci data
  XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
  XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]

  kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]

  data.RKHS.train <- as.data.frame(t(XX.train))
  data.RKHS.train$Y <- factor(Y.train)
  
  data.RKHS.test <- as.data.frame(t(XX.test))
  data.RKHS.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                      type = 'C-classification',
                      scale = TRUE,
                      cost = C,
                      coef0 = coef0,
                      kernel = kernel_type)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  # ulozeni vysledku
  Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
}
```


Table: (\#tab:unnamed-chunk-99)Souhrnné výsledky metody SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

Model                                $\widehat{Err}_{train}\quad\quad\quad\quad\quad$         $\widehat{Err}_{test}\quad\quad\quad\quad\quad$       
----------------------------------  -------------------------------------------------------  -------------------------------------------------------
SVM linear - RKHS - radial                                                                0                                                   0.0308
SVM poly - RKHS - radial                                                                  0                                                   0.0308
SVM rbf - RKHS - radial                                                                   0                                                   0.0154

Přesnost metody SVM v kombinaci s projekcí na Reproducing Kernel Hilbert Space je tedy na trénovacích datech rovna 0 % pro lineární jádro, 0 % pro polynomiální jádro a 0 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 3.08 % pro lineární jádro, 3.08 % pro polynomiální jádro a 1.54 % pro radiální jádro.


```r
RESULTS <- rbind(RESULTS, Res)
```

###### Polynomiální jádro


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())

# jadro a jadrova matice ... polynomialni s parametrem p
Poly.kernel <- function(x, y, p) {
  return((1 + x * y)^p)
}

Kernel.RKHS <- function(x, p) {
  K <- matrix(NA, ncol = length(x), nrow = length(x))
  for(i in 1:nrow(K)) {
    for(j in 1:ncol(K)) {
      K[i, j] <- Poly.kernel(x = x[i], y = x[j], p)
    }
  }
  return(K)
}
```


```r
# rozdelime trenovaci data na k casti
folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()

# hodnoty hyperparametru, ktere budeme prochazet
dimensions <- 2:30 # rozumny rozsah hodnot d
poly.cv <- 2:5

# list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro dane
# v radcich budou hodnoty pro dane p a vrstvy odpovidaji folds
dim.names <- list(p = paste0('p:', poly.cv),
                  d = paste0('d:', dimensions),
                  CV = paste0('cv:', 1:k_cv))

CV.results <- list(
  SVM.l = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.p = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.r = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                dimnames = dim.names))
```


```r
# samotna CV
for (p in poly.cv) {
  K <- Kernel.RKHS(t.seq, p = p)
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'polynomial',
                    type = 'eps-regression',
                    coef0 = 1,
                    cost = C,
                    epsilon = eps,
                    degree = p)
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
  }
  
  # projdeme dimenze
  for(d.RKHS in dimensions) {
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) 
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                             alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    # projdeme folds
    for (index_cv in 1:k_cv) {
      # definice testovaci a trenovaci casti pro CV
      fold <- folds[[index_cv]]
      # rozdeleni na trenovaci a validacni data
      XX.train <- Lambda.RKHS[, fold]
      XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
      # pripravime si datovou tabulku pro ulozeni vysledku
      Res <- data.frame(model = c('SVM linear - RKHS', 
                                  'SVM poly - RKHS', 
                                  'SVM rbf - RKHS'), 
                        Err.test = NA)
      # projdeme jednotliva jadra
      for (kernel_number in 1:3) {
        kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
      
        data.RKHS.train <- as.data.frame(t(XX.train))
        data.RKHS.train$Y <- factor(Y.train[fold])
        
        data.RKHS.test <- as.data.frame(t(XX.test))
        data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
        
        # sestrojeni modelu
        clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                            type = 'C-classification',
                            scale = TRUE,
                            cost = C,                    
                            coef0 = 1,
                            gamma = 1,
                            kernel = kernel_type)
        
        # presnost na validacnich datech
        predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
        presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
          prop.table() |> diag() |> sum()
        
        # ulozeni vysledku
        Res[kernel_number, 2] <- 1 - presnost.test
      }
      # presnosti vlozime na pozice pro dane d, gamma a fold
      CV.results$SVM.l[paste0('p:', p), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[1, 2]
      CV.results$SVM.p[paste0('p:', p), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[2, 2]
      CV.results$SVM.r[paste0('p:', p), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[3, 2]
    }
  }
}
```


```r
# spocitame prumerne presnosti pro jednotliva d pres folds
for (n_method in 1:length(CV.results)) {
  CV.results[[n_method]] <- apply(CV.results[[n_method]], c(1, 2), mean)
}

poly.opt <- c(which.min(CV.results$SVM.l) %% length(poly.cv), 
               which.min(CV.results$SVM.p) %% length(poly.cv), 
               which.min(CV.results$SVM.r) %% length(poly.cv))
poly.opt[poly.opt == 0] <- length(poly.cv)
poly.opt <- poly.cv[poly.opt]

d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
           which.min(t(CV.results$SVM.p)) %% length(dimensions), 
           which.min(t(CV.results$SVM.r)) %% length(dimensions))
d.opt[d.opt == 0] <- length(dimensions)
d.opt <- dimensions[d.opt]

err.opt.cv <- c(min(CV.results$SVM.l), 
                     min(CV.results$SVM.p),
                     min(CV.results$SVM.r))
df.RKHS.res <- data.frame(d = d.opt, p = poly.opt, CV = err.opt.cv,
           Kernel = c('linear', 'polynomial', 'radial') |> factor(),
           row.names = c('linear', 'poly', 'radial'))
```


Table: (\#tab:unnamed-chunk-104)Souhrnné výsledky cross-validace pro metodu SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

          $\quad\quad\quad\quad\quad d$   $\quad\quad\quad\quad\quad p$   $\widehat{Err}_{cross\_validace}$  Model                             
-------  ------------------------------  ------------------------------  ----------------------------------  ----------------------------------
linear                               20                               5                              0.0474  linear                            
poly                                 10                               3                              0.0461  polynomial                        
radial                                7                               5                              0.0403  radial                            

Vidíme, že nejlépe vychází hodnota parametru $d={}$ 20 a $p={}$ 5 pro lineární jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.9526, $d={}$ 10 a $p={}$ 3 pro polynomiální jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.9539 a $d={}$ 7 a $p={}$ 5 pro radiální jádro s hodnotou přesnosti 0.9597.

Jelikož již máme nalezeny optimální hodnoty hyperparametrů, můžeme zkounstruovat finální modely a určit jejich úspěšnost klasifikace na testovacích datech.


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
```


```r
# pripravime si datovou tabulku pro ulozeni vysledku
Res <- data.frame(model = c('SVM linear - RKHS - poly', 
                            'SVM poly - RKHS - poly', 
                            'SVM rbf - RKHS - poly'), 
                  Err.train = NA,
                  Err.test = NA)

# projdeme jednotliva jadra
for (kernel_number in 1:3) {
  # spocitame matici K
  p <- poly.opt[kernel_number] # hodnota gamma pomoci CV
  K <- Kernel.RKHS(t.seq, p = p)
  
  # urcime vlastni cisla a vektory
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  # urceni koeficientu alpha z SVM
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                       ncol = dim(data.RKHS)[2]) # prazdny objekt
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'polynomial',
                    type = 'eps-regression',
                    epsilon = eps,
                    coef0 = 1,
                    cost = C,
                    gamma = 1,
                    degree = p)
    # urceni alpha
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
  }
  # d
  d.RKHS <- d.opt[kernel_number]
  
  # urceni vektoru lambda
  Lambda.RKHS <- matrix(NA, 
                        ncol = dim(data.RKHS)[2], 
                        nrow = d.RKHS) # vytvoreni prazdneho objektu
  
  # vypocet reprezentace
  for(l in 1:dim(data.RKHS)[2]) {
    Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
  }
  
  # rozdeleni na trenovaci a testovaci data
  XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
  XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]

  kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]

  data.RKHS.train <- as.data.frame(t(XX.train))
  data.RKHS.train$Y <- factor(Y.train)
  
  data.RKHS.test <- as.data.frame(t(XX.test))
  data.RKHS.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                      type = 'C-classification',
                      scale = TRUE,
                      coef0 = 1,
                      cost = C,
                      gamma = 1,
                      kernel = kernel_type)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  # ulozeni vysledku
  Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
}
```


Table: (\#tab:unnamed-chunk-107)Souhrnné výsledky metody SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

Model                                $\widehat{Err}_{train}\quad\quad\quad\quad\quad$         $\widehat{Err}_{test}\quad\quad\quad\quad\quad$       
----------------------------------  -------------------------------------------------------  -------------------------------------------------------
SVM linear - RKHS - poly                                                             0.0333                                                   0.0615
SVM poly - RKHS - poly                                                               0.0267                                                   0.1077
SVM rbf - RKHS - poly                                                                0.0333                                                   0.1077

Přesnost metody SVM v kombinaci s projekcí na Reproducing Kernel Hilbert Space je tedy na trénovacích datech rovna 3.33 % pro lineární jádro, 2.67 % pro polynomiální jádro a 3.33 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 6.15 % pro lineární jádro, 10.77 % pro polynomiální jádro a 10.77 % pro radiální jádro.


```r
RESULTS <- rbind(RESULTS, Res)
```

###### Lineární jádro


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())

# jadro a jadrova matice ... polynomialni s parametrem p
Linear.kernel <- function(x, y) {
  return(x * y)
}

Kernel.RKHS <- function(x) {
  K <- matrix(NA, ncol = length(x), nrow = length(x))
  for(i in 1:nrow(K)) {
    for(j in 1:ncol(K)) {
      K[i, j] <- Linear.kernel(x = x[i], y = x[j])
    }
  }
  return(K)
}
```


```r
# rozdelime trenovaci data na k casti
folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()

# hodnoty hyperparametru, ktere budeme prochazet
dimensions <- 2:40 # rozumny rozsah hodnot d

# list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro dane d
# v radcich budou hodnoty pro vrstvy odpovidaji folds
dim.names <- list(d = paste0('d:', dimensions),
                  CV = paste0('cv:', 1:k_cv))

CV.results <- list(
  SVM.l = array(NA, dim = c(length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.p = array(NA, dim = c(length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.r = array(NA, dim = c(length(dimensions), k_cv),
                dimnames = dim.names))
```


```r
# samotna CV
K <- Kernel.RKHS(t.seq)
Eig <- eigen(K)
eig.vals <- Eig$values
eig.vectors <- Eig$vectors
alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 

# model
for(i in 1:dim(data.RKHS)[2]) {
  df.svm <- data.frame(x = t.seq,
                       y = data.RKHS[, i])
  svm.RKHS <- svm(y ~ x, data = df.svm, 
                  kernel = 'linear',
                  type = 'eps-regression',
                  epsilon = eps,                   
                  coef0 = 1,
                  gamma = 1,
                  cost = C)
  alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
}

# projdeme dimenze
for(d.RKHS in dimensions) {
  Lambda.RKHS <- matrix(NA, 
                        ncol = dim(data.RKHS)[2], 
                        nrow = d.RKHS) 
  # vypocet reprezentace
  for(l in 1:dim(data.RKHS)[2]) {
    Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                           alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
  }
  # projdeme folds
  for (index_cv in 1:k_cv) {
    # definice testovaci a trenovaci casti pro CV
    fold <- folds[[index_cv]]
    # rozdeleni na trenovaci a validacni data
    XX.train <- Lambda.RKHS[, fold]
    XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
    # pripravime si datovou tabulku pro ulozeni vysledku
    Res <- data.frame(model = c('SVM linear - RKHS', 
                                'SVM poly - RKHS', 
                                'SVM rbf - RKHS'), 
                      Err.test = NA)
    # projdeme jednotliva jadra
    for (kernel_number in 1:3) {
      kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
    
      data.RKHS.train <- as.data.frame(t(XX.train))
      data.RKHS.train$Y <- factor(Y.train[fold])
      
      data.RKHS.test <- as.data.frame(t(XX.test))
      data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
      
      # sestrojeni modelu
      clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                          type = 'C-classification',
                          scale = TRUE,
                          kernel = kernel_type,
                          epsilon = eps,                   
                          coef0 = 1,
                          gamma = 1,
                          cost = C)
      
      # presnost na validacnich datech
      predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
      presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
        prop.table() |> diag() |> sum()
      
      # ulozeni vysledku
      Res[kernel_number, 2] <- 1 - presnost.test
    }
    # presnosti vlozime na pozice pro dane d, gamma a fold
    CV.results$SVM.l[d.RKHS - min(dimensions) + 1, 
                     index_cv] <- Res[1, 2]
    CV.results$SVM.p[d.RKHS - min(dimensions) + 1, 
                     index_cv] <- Res[2, 2]
    CV.results$SVM.r[d.RKHS - min(dimensions) + 1, 
                     index_cv] <- Res[3, 2]
  }
}
```


```r
# spocitame prumerne presnosti pro jednotliva d pres folds
for (n_method in 1:length(CV.results)) {
  CV.results[[n_method]] <- apply(CV.results[[n_method]], 1, mean)
}

d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
           which.min(t(CV.results$SVM.p)) %% length(dimensions), 
           which.min(t(CV.results$SVM.r)) %% length(dimensions))
d.opt[d.opt == 0] <- length(dimensions)
d.opt <- dimensions[d.opt]

err.opt.cv <- c(min(CV.results$SVM.l), 
                     min(CV.results$SVM.p),
                     min(CV.results$SVM.r))
df.RKHS.res <- data.frame(d = d.opt, CV = err.opt.cv,
           Kernel = c('linear', 'polynomial', 'radial') |> factor(),
           row.names = c('linear', 'poly', 'radial'))
```


Table: (\#tab:unnamed-chunk-112)Souhrnné výsledky cross-validace pro metodu SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

          $\quad\quad\quad\quad\quad d$   $\widehat{Err}_{cross\_validace}$  Model                             
-------  ------------------------------  ----------------------------------  ----------------------------------
linear                               15                              0.0667  linear                            
poly                                 16                              0.0454  polynomial                        
radial                               25                              0.0526  radial                            

Vidíme, že nejlépe vychází hodnota parametru $d={}$ 15 pro lineární jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.9333, $d={}$ 16 pro polynomiální jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.9546 a $d={}$ 25 pro radiální jádro s hodnotou přesnosti 0.9474.

Jelikož již máme nalezeny optimální hodnoty hyperparametrů, můžeme zkounstruovat finální modely a určit jejich úspěšnost klasifikace na testovacích datech.


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
```


```r
# pripravime si datovou tabulku pro ulozeni vysledku
Res <- data.frame(model = c('SVM linear - RKHS - linear', 
                            'SVM poly - RKHS - linear', 
                            'SVM rbf - RKHS - linear'), 
                  Err.train = NA,
                  Err.test = NA)

# projdeme jednotliva jadra
for (kernel_number in 1:3) {
  # spocitame matici K
  K <- Kernel.RKHS(t.seq)
  
  # urcime vlastni cisla a vektory
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  # urceni koeficientu alpha z SVM
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                       ncol = dim(data.RKHS)[2]) # prazdny objekt
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'linear',
                    type = 'eps-regression',
                    epsilon = eps,                   
                    coef0 = 1,
                    gamma = 1,
                    cost = C)
    # urceni alpha
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
  }
  # d
  d.RKHS <- d.opt[kernel_number]
  
  # urceni vektoru lambda
  Lambda.RKHS <- matrix(NA, 
                        ncol = dim(data.RKHS)[2], 
                        nrow = d.RKHS) # vytvoreni prazdneho objektu
  
  # vypocet reprezentace
  for(l in 1:dim(data.RKHS)[2]) {
    Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
  }
  
  # rozdeleni na trenovaci a testovaci data
  XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
  XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]

  kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]

  data.RKHS.train <- as.data.frame(t(XX.train))
  data.RKHS.train$Y <- factor(Y.train)
  
  data.RKHS.test <- as.data.frame(t(XX.test))
  data.RKHS.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                      type = 'C-classification',
                      scale = TRUE,
                      kernel = kernel_type,
                      epsilon = eps,                   
                      coef0 = 1,
                      gamma = 1,
                      cost = C)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  # ulozeni vysledku
  Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
}
```


Table: (\#tab:unnamed-chunk-115)Souhrnné výsledky metody SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

Model                                $\widehat{Err}_{train}\quad\quad\quad\quad\quad$         $\widehat{Err}_{test}\quad\quad\quad\quad\quad$       
----------------------------------  -------------------------------------------------------  -------------------------------------------------------
SVM linear - RKHS - linear                                                           0.0400                                                   0.0923
SVM poly - RKHS - linear                                                             0.0133                                                   0.0308
SVM rbf - RKHS - linear                                                              0.0200                                                   0.0308

Přesnost metody SVM v kombinaci s projekcí na Reproducing Kernel Hilbert Space je tedy na trénovacích datech rovna 4 % pro lineární jádro, 1.33 % pro polynomiální jádro a 2 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 9.23 % pro lineární jádro, 3.08 % pro polynomiální jádro a 3.08 % pro radiální jádro.


```r
RESULTS <- rbind(RESULTS, Res)
```

## Tabulka výsledků

Vidíme z tabulky níže, že jednotlivé klasifikační metody mají mezi sebou výrazné rozdíly co do úspěšnosti klasifikace. Zejména klasické metody, jako je KNN, LDA nebo QDA si vedou velmi bídně. Můžeme si všimnout, že všechny metody postavené na funkcionální analýze hlavních komponent nedosahují zdaleka podobných výsledků jako některé jiné metody.

Naopak nyní se vymyká svou dobrou klasifikační schopností metoda RKHS společně s SVM. Poznamenejme, že také klasická SVM s lineárním jádrem si vede velmi obstojně. Obecně je lineární jádro dobrou volbou (jak jsme se ostatně mohli přesvědčit již dříve), neboť pro dostatečně hustou síť bodů dobře aproximuje určitý integrál na uvažovaném intervalu $I$.


Table: (\#tab:unnamed-chunk-117)Souhrnné výsledky použitých metod na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

Model                                $\widehat{Err}_{train}\quad\quad\quad\quad\quad$         $\widehat{Err}_{test}\quad\quad\quad\quad\quad$       
----------------------------------  -------------------------------------------------------  -------------------------------------------------------
KNN                                                                                  0.1467                                                   0.1692
LDA                                                                                  0.3200                                                   0.2923
QDA                                                                                  0.3200                                                   0.3077
LR functional                                                                        0.0000                                                   0.0615
LR score                                                                             0.3067                                                   0.2923
Tree - diskr.                                                                        0.2933                                                   0.3538
Tree - score                                                                         0.3267                                                   0.3846
Tree - Bbasis                                                                        0.2267                                                   0.2462
RForest - diskr                                                                      0.0200                                                   0.1231
RForest - score                                                                      0.0467                                                   0.3077
RForest - Bbasis                                                                     0.0133                                                   0.1231
SVM linear - diskr                                                                   0.0067                                                   0.0769
SVM poly - diskr                                                                     0.0133                                                   0.0769
SVM rbf - diskr                                                                      0.0133                                                   0.0462
SVM linear - PCA                                                                     0.3200                                                   0.2923
SVM poly - PCA                                                                       0.3333                                                   0.2769
SVM rbf - PCA                                                                        0.1600                                                   0.3692
SVM linear - Bbasis                                                                  0.0067                                                   0.0615
SVM poly - Bbasis                                                                    0.0067                                                   0.0923
SVM rbf - Bbasis                                                                     0.0133                                                   0.0615
SVM linear - projection                                                              0.0200                                                   0.0615
SVM poly - projection                                                                0.0267                                                   0.0615
SVM rbf - projection                                                                 0.0933                                                   0.1077
SVM linear - RKHS - radial                                                           0.0000                                                   0.0308
SVM poly - RKHS - radial                                                             0.0000                                                   0.0308
SVM rbf - RKHS - radial                                                              0.0000                                                   0.0154
SVM linear - RKHS - poly                                                             0.0333                                                   0.0615
SVM poly - RKHS - poly                                                               0.0267                                                   0.1077
SVM rbf - RKHS - poly                                                                0.0333                                                   0.1077
SVM linear - RKHS - linear                                                           0.0400                                                   0.0923
SVM poly - RKHS - linear                                                             0.0133                                                   0.0308
SVM rbf - RKHS - linear                                                              0.0200                                                   0.0308

## Klasifikace pomocí druhé derivace {#klasA3deriv}

Jak jsme již avizovali dříve, pro tato data je vhodné ke klasifikaci uvažovat jejich druhou derivaci. Tu jsme si již spočetli výše, proto už se nyní můžeme pustit rovnou do konstrukce modelů.

Proveďme obdobnou analýzu jako v situaci výše, následně (jelikož data náhodně rozdělujeme na testovací a trénovací část), provedeme simulační studii, pomocí které budeme jednotlivé klasifikační metody schopni lépe a s mnohem větší silou porovnat.


```r
# rozdeleni na testovaci a trenovaci cast
set.seed(42)
split <- sample.split(XXder$fdnames$reps, SplitRatio = 0.7)

# vytvoreni vektoru 0 a 1, 0 pro < 20 a 1 pro > 20 
Y <- ifelse(labels == 'large', 1, 0)

X.train <- subset(XXder, split == TRUE)
X.test <- subset(XXder, split == FALSE)

Y.train <- subset(Y, split == TRUE)
Y.test <- subset(Y, split == FALSE)
```

Ještě se podíváme na zastoupení jednotlivých skupin v testovací a trénovací části dat.


```r
# absolutni zastoupeni
table(Y.train)
```

```
## Y.train
##  0  1 
## 91 59
```

```r
table(Y.test)
```

```
## Y.test
##  0  1 
## 47 18
```

```r
# relativni zastoupeni
table(Y.train) / sum(table(Y.train))
```

```
## Y.train
##         0         1 
## 0.6066667 0.3933333
```

```r
table(Y.test) / sum(table(Y.test))
```

```
## Y.test
##         0         1 
## 0.7230769 0.2769231
```

### $K$ nejbližších sousedů

Začněme neparametrickou klasifikační metodou, a to metodou $K$ nejbližších sousedů.
Nejprve si vytvoříme potřebné objekty tak, abychom s nimi mohli pomocí funkce `classif.knn()` z knihovny `fda.usc` dále pracovat.


```r
x.train <- fdata(X.train)
y.train <- as.numeric(factor(Y.train))
```

Nyní můžeme definovat model a podívat se na jeho úspěšnost klasifikace.
Poslední otázkou však zůstává, jak volit optimální počet sousedů $K$.
Mohli bychom tento počet volit jako takové $K$, při kterém nastává minimální chybovost na trénovacích datech.
To by ale mohlo vést k přeučení modelu, proto využijeme cross-validaci.
Vzhledem k výpočetní náročnosti a rozsahu souboru zvolíme $k$-násobnou CV, my zvolíme například hodnotu $k = {10}$.


```r
# model pro vsechna trenovaci data pro K = 1, 2, ..., sqrt(n_train)
neighb.model <- classif.knn(group = y.train, 
                            fdataobj = x.train, 
                            knn = c(1:round(sqrt(length(y.train))))) 

neighb.model$max.prob # maximalni presnost
```

```
## [1] 0.9866667
```

```r
(K.opt <- neighb.model$h.opt) # optimalni hodnota K
```

```
## [1] 3
```

Proveďme předchozí postup pro trénovací data, která rozdělíme na $k$ částí a tedy zopakujeme tuto část kódu $k$-krát.


```r
k_cv <- 10 # k-fold CV
neighbours <- c(1:(2 * ceiling(sqrt(length(y.train))))) # pocet sousedu 

# rozdelime trenovaci data na k casti
folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)

# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro danou cast trenovaci mnoziny
# v radcich budou hodnoty pro danou hodnotu K sousedu
CV.results <- matrix(NA, nrow = length(neighbours), ncol = k_cv)

for (index in 1:k_cv) {
  # definujeme danou indexovou mnozinu
  fold <- folds[[index]]
    
  x.train.cv <- subset(X.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
    fdata()
  y.train.cv <- subset(Y.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
    factor() |> as.numeric()
  
  x.test.cv <- subset(X.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
    fdata()
  y.test.cv <- subset(Y.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
    factor() |> as.numeric()
  
  # projdeme kazdou cast ... k-krat zopakujeme
  for(neighbour in neighbours) {
    # model pro konkretni volbu K
    neighb.model <- classif.knn(group = y.train.cv, 
                              fdataobj = x.train.cv, 
                              knn = neighbour) 
    # predikce na validacni casti
    model.neighb.predict <- predict(neighb.model, 
                                    new.fdataobj = x.test.cv)
    # presnost na validacni casti
    presnost <- table(y.test.cv, model.neighb.predict) |> 
      prop.table() |> diag() |> sum()
    
    # presnost vlozime na pozici pro dane K a fold
    CV.results[neighbour, index] <- presnost
  }
}

# spocitame prumerne presnosti pro jednotliva K pres folds
CV.results <- apply(CV.results, 1, mean)
K.opt <- which.max(CV.results)
presnost.opt.cv <- max(CV.results)
CV.results <- data.frame(K = neighbours, CV = CV.results)
```

Vidíme, že nejlépe vychází hodnota parametru $K$ jako 3 s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.9897.
Pro přehlednost si ještě vykresleme průběh validační chybovosti v závislosti na počtu sousedů $K$.


```r
CV.results |> ggplot(aes(x = K, y = 1 - CV)) + 
  geom_line(linetype = 'dashed', colour = 'grey') + 
  geom_point(size = 1.5) + 
  geom_point(aes(x = K.opt, y = 1 - presnost.opt.cv), colour = 'red', size = 2) +
  theme_bw() + 
  labs(x = bquote(paste(K, ' ;   ', 
                        K[optimal] == .(K.opt))),
       y = 'Validační chybovost') + 
  scale_x_continuous(breaks = neighbours)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-123-1.png" alt="Závislost validační chybovosti na hodnotě $K$, tedy na počtu sousedů." width="672" />
<p class="caption">(\#fig:unnamed-chunk-123)Závislost validační chybovosti na hodnotě $K$, tedy na počtu sousedů.</p>
</div>

Nyní známe optimální hodnotu parametru $K$ a tudíž můžeme sestavit finální model.


```r
neighb.model <- classif.knn(group = y.train, fdataobj = x.train, knn = K.opt)

# predikce
model.neighb.predict <- predict(neighb.model, 
                                new.fdataobj = fdata(X.test))

# presnost na testovacich datech
presnost <- table(as.numeric(factor(Y.test)), model.neighb.predict) |>
  prop.table() |>
  diag() |>
  sum()
```

Vidíme tedy, že přesnost modelu sestrojeného pomocí metody $K$ nejbližších sousedů s optimální volbou $K_{optimal}$ rovnou 3, kterou jsme určili cross-validací, je na trénovacích datech rovna 0.0133 a na testovacích datech 0.0769.

K porovnání jendotlivých modelů můžeme použít oba typy chybovostí, pro přehlednost si je budeme ukládat do tabulky.


```r
RESULTS <- data.frame(model = 'KNN', 
                      Err.train = 1 - neighb.model$max.prob,
                      Err.test = 1 - presnost)
```

### Lineární diskriminační analýza

Jako druhou metodu pro sestrojení klasifikátoru budeme uvažovat lineární diskriminační analýzu (LDA).
Jelikož tato metoda nelze aplikovat na funkcionální data, musíme je nejprve diskretizovat, což provedeme pomocí funkcionální analýzy hlavních komponent.
Klasifikační algoritmus následně provedeme na skórech prvních $p$ hlavních komponent.
Počet komponent $p$ zvolíme tak, aby prvních $p$ hlavních komponent dohromady vysvětlovalo alespoň 90 % variability v datech.

Proveďme tedy nejprve funkcionální analýzu hlavních komponent a určeme počet $p$.


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

V tomto konkrétním případě jsme za počet hlavních komponent vzali $p=$ 2, které dohromady vysvětlují 93.12 $\%$ variability v datech.
První hlavní komponenta potom vysvětluje 77.7 % a druhá 15.42 $\%$ variability.
Graficky si můžeme zobrazit hodnoty skórů prvních dvou hlavních komponent, barevně odlišených podle příslušnosti do klasifikační třídy.


```r
data.PCA.train |> ggplot(aes(x = V1, y = V2, colour = Y)) +
  geom_point(size = 1.5) + 
  labs(x = paste('1. hlavní komponenta (vysvětlená variabilita', 
                 round(100 * data.PCA$varprop[1], 2), '%)'),
       y = paste('2. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[2], 2), '%)'),
       colour = 'Obsah tuku') +
  scale_color_discrete(labels = c("malý", "velký")) + 
  theme_bw()
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-127-1.png" alt="Hodnoty skórů prvních dvou hlavních komponent pro trénovací data. Barevně jsou odlišeny body podle příslušnosti do klasifikační třídy." width="672" />
<p class="caption">(\#fig:unnamed-chunk-127)Hodnoty skórů prvních dvou hlavních komponent pro trénovací data. Barevně jsou odlišeny body podle příslušnosti do klasifikační třídy.</p>
</div>

Abychom mohli určit přesnost klasifikace na testovacích datech, potřebujeme spočítat skóre pro první 2 hlavní komponenty pro testovací data.
Tato skóre určíme pomocí vzorce:

$$
\xi_{i, j} = \int \left( X_i(t) - \mu(t)\right) \cdot \rho_j(t)\text dt,
$$ 
kde $\mu(t)$ je střední hodnota (průměrná funkce) a $\rho_j(t)$ vlastní fukce (funkcionální hlavní komponenty).


```r
# vypocet skoru testovacich funkci
scores <- matrix(NA, ncol = nharm, nrow = length(Y.test)) # prazdna matice

for(k in 1:dim(scores)[1]) {
  xfd = X.test[k] - data.PCA$meanfd[1] # k-te pozorovani - prumerna funkce
  scores[k, ] = inprod(xfd, data.PCA$harmonics) 
  # skalarni soucin rezidua a vlastnich funkci rho (funkcionalni hlavni komponenty)
}

data.PCA.test <- as.data.frame(scores)
data.PCA.test$Y <- factor(Y.test)
colnames(data.PCA.test) <- colnames(data.PCA.train) 
```

Nyní již můžeme sestrojit klasifikátor na trénovací části dat.


```r
# model
clf.LDA <- lda(Y ~ ., data = data.PCA.train)

# presnost na trenovacich datech
predictions.train <- predict(clf.LDA, newdata = data.PCA.train)
presnost.train <- table(data.PCA.train$Y, predictions.train$class) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.LDA, newdata = data.PCA.test)
presnost.test <- table(data.PCA.test$Y, predictions.test$class) |>
  prop.table() |> diag() |> sum()
```

Spočítali jsme jednak přesnost klasifikátoru na trénovacích (96 %), tak i na testovacích datech (90.77 %).

Pro grafické znázornění metody můžeme zaznačit dělící hranici do grafu skórů prvních dvou hlavních komponent.
Tuto hranici spočítáme na husté síti bodů a zobrazíme ji pomocí funkce `geom_contour()`.


```r
# pridame diskriminacni hranici
np <- 1001 # pocet bodu site
# x-ova osa ... 1. HK
nd.x <- seq(from = min(data.PCA.train$V1), 
            to = max(data.PCA.train$V1), length.out = np)
# y-ova osa ... 2. HK
nd.y <- seq(from = min(data.PCA.train$V2), 
            to = max(data.PCA.train$V2), length.out = np)
# pripad pro 2 HK ... p = 2
nd <- expand.grid(V1 = nd.x, V2 = nd.y)
# pokud p = 3
if(dim(data.PCA.train)[2] == 4) {
  nd <- expand.grid(V1 = nd.x, V2 = nd.y, V3 = data.PCA.train$V3[1])}
# pokud p = 4
if(dim(data.PCA.train)[2] == 5) {
  nd <- expand.grid(V1 = nd.x, V2 = nd.y, V3 = data.PCA.train$V3[1],
                    V4 = data.PCA.train$V4[1])}
# pokud p = 5
if(dim(data.PCA.train)[2] == 6) {
  nd <- expand.grid(V1 = nd.x, V2 = nd.y, V3 = data.PCA.train$V3[1],
                    V4 = data.PCA.train$V4[1], V5 = data.PCA.train$V5[1])}

# pridame Y = 0, 1
nd <- nd |> mutate(prd = as.numeric(predict(clf.LDA, newdata = nd)$class))

data.PCA.train |> ggplot(aes(x = V1, y = V2, colour = Y)) +
  geom_point(size = 1.5) + 
  labs(x = paste('1. hlavní komponenta (vysvětlená variabilita', 
                 round(100 * data.PCA$varprop[1], 2), '%)'),
       y = paste('2. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[2], 2), '%)'),
       colour = 'Obsah tuku') +
  scale_color_discrete(labels = c("malý", "velký")) + 
  theme_bw() +
  geom_contour(data = nd, aes(x = V1, y = V2, z = prd), colour = 'black')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-130-1.png" alt="Skóre prvních dvou hlavních komponent, barevně odlišené podle pohlaví. Černě je vyznačena dělící hranice (přímka v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí LDA." width="672" />
<p class="caption">(\#fig:unnamed-chunk-130)Skóre prvních dvou hlavních komponent, barevně odlišené podle pohlaví. Černě je vyznačena dělící hranice (přímka v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí LDA.</p>
</div>

Vidíme, že dělící hranicí je přímka, lineární funkce v prostoru 2D, což jsme ostatně od LDA čekali.
Nakonec přidáme chybovosti do souhrnné tabulky.


```r
Res <- data.frame(model = 'LDA', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

### Kvadratická diskriminační analýza

Jako další sestrojme klasifikátor pomocí kvadratické diskriminační analýzy (QDA).
Jedná se o analogický případ jako LDA s tím rozdílem, že nyní připouštíme pro každou ze tříd rozdílnou kovarianční matici normálního rozdělení, ze kterého pocházejí příslušné skóry.
Tento vypuštěný předpoklad o rovnosti kovariančních matic vede ke kvadratické hranici mezi třídami.

V `R` se provede QDA analogicky jako LDA v předchozí části, tedy opět bychom pomocí funkcionální analýzy hlavních komponent spočítali skóre pro trénovací i testovací funkce, sestrojili klasifikátor na skórech prvních $p$ hlavních komponent a pomocí něj predikovali příslušnost testovacích křivek do třídy $Y^* \in \{0, 1\}$.

Funkcionální PCA provádět nemusíme, využijeme výsledků z části LDA.





Můžeme tedy rovnou přistoupit k sestrojení klasifikátoru, což provedeme pomocí funkce `qda()`.
Následně spočítáme přesnost klasifikátoru na testovacích a trénovacích datech.


```r
# model
clf.QDA <- qda(Y ~ ., data = data.PCA.train)

# presnost na trenovacich datech
predictions.train <- predict(clf.QDA, newdata = data.PCA.train)
presnost.train <- table(data.PCA.train$Y, predictions.train$class) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.QDA, newdata = data.PCA.test)
presnost.test <- table(data.PCA.test$Y, predictions.test$class) |>
  prop.table() |> diag() |> sum()
```

Spočítali jsme tedy jednak přesnost klasifikátoru na trénovacích (99.33 %), tak i na testovacích datech (98.46 %).

Pro grafické znázornění metody můžeme zaznačit dělící hranici do grafu skórů prvních dvou hlavních komponent.
Tuto hranici spočítáme na husté síti bodů a zobrazíme ji pomocí funkce `geom_contour()` stejně jako v případě LDA.




```r
nd <- nd |> mutate(prd = as.numeric(predict(clf.QDA, newdata = nd)$class))

data.PCA.train |> ggplot(aes(x = V1, y = V2, colour = Y)) +
  geom_point(size = 1.5) + 
  labs(x = paste('1. hlavní komponenta (vysvětlená variabilita', 
                 round(100 * data.PCA$varprop[1], 2), '%)'),
       y = paste('2. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[2], 2), '%)'),
       colour = 'Obsah tuku') +
  scale_color_discrete(labels = c("malý", "velký")) + 
  theme_bw() +
  geom_contour(data = nd, aes(x = V1, y = V2, z = prd), colour = 'black')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-136-1.png" alt="Skóre prvních dvou hlavních komponent, barevně odlišené podle pohlaví. Černě je vyznačena dělící hranice (parabola v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí QDA." width="672" />
<p class="caption">(\#fig:unnamed-chunk-136)Skóre prvních dvou hlavních komponent, barevně odlišené podle pohlaví. Černě je vyznačena dělící hranice (parabola v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí QDA.</p>
</div>

Všimněme si, že dělící hranicí mezi klasifikačními třídami je nyní parabola, avšak se jen (alespoň opticky) velmi málo liší od přímky.

Nakonec ještě doplníme chybovosti do souhrnné tabulky.


```r
Res <- data.frame(model = 'QDA', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

### Logistická regrese

Logistickou regresi můžeme provést dvěma způsoby.
Jednak použít funkcionální obdobu klasické logistické regrese, druhak klasickou mnohorozměrnou logistickou regresi, kterou provedeme na skórech prvních $p$ hlavních komponent.

#### Funkcionální logistická regrese

Analogicky jako v případě konečné dimenze vstupních dat uvažujeme logistický model ve tvaru:

$$
g\left(\mathbb E [Y|X = x]\right) = \eta (x) = g(\pi(x)) = \alpha + \int \beta(t)\cdot x(t) \text d t,
$$ 
kde $\eta(x)$ je lineární prediktor nabývající hodnot z intervalu $(-\infty, \infty)$, $g(\cdot)$ je *linková funkce*, v případě logistické regrese se jedná o logitovou funkci $g: (0,1) \rightarrow \mathbb R,\ g(p) = \ln\frac{p}{1-p}$ a $\pi(x)$ podmíněná pravděpodobnost

$$
\pi(x) = \text{Pr}(Y = 1 | X = x) = g^{-1}(\eta(x)) = \frac{\text e^{\alpha + \int \beta(t)\cdot x(t) \text d t}}{1 + \text e^{\alpha + \int \beta(t)\cdot x(t) \text d t}},
$$

přičemž $\alpha$ je konstanta a $\beta(t) \in L^2[a, b]$ je parametrická funkce.
Naším cílem je odhadnout tuto parametrickou funkci.

Pro funkcionální logistickou regresi použijeme funkci `fregre.glm()` z balíčku `fda.usc`.
Nejprve si vytvoříme vhodné objekty pro konstrukci klasifikátoru.


```r
# vytvorime vhodne objekty
x.train <- fdata(X.train)
y.train <- as.numeric(Y.train) 

# body, ve kterych jsou funkce vyhodnoceny
tt <- x.train[["argvals"]]

dataf <- as.data.frame(y.train) 
colnames(dataf) <- "Y"

nbasis.x <- 7

# B-spline baze 
basis1 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.x)
```

Abychom mohli odhadnout parametrickou funkci $\beta(t)$, potřebujeme ji vyjádřit v nějaké bazické reprezentaci, v našem případě B-splinové bázi.
K tomu však potřebujeme najít vhodný počet bázových funkcí.
To bychom mohli určit na základě chybovosti na trénovacích datech, avšak tato data budou upřenostňovat výběr velkého počtu bází a bude docházet k přeučení modelu.

Ilustrujme si to na následujícím případě.
Pro každý z počtu bází $n_{basis} \in \{4, 5, \dots, 30\}$ natrénujeme model na trénovacích datech, určíme na nich chybovost a také spočítáme chybovost na testovacích datech.
Připomeňme, že k výběru vhodného počtu bází nemůžeme využít stejná data jako pro odhad testovací chybovosti, neboť bychom tuto chybovost podcenili.


```r
n.basis.max <- 30
n.basis <- 4:n.basis.max
pred.baz <- matrix(NA, nrow = length(n.basis), ncol = 2, 
                   dimnames = list(n.basis, c('Err.train', 'Err.test')))

for (i in n.basis) {
  # baze pro bety
  basis2 <- create.bspline.basis(rangeval = range(tt), nbasis = i)
  # vztah
  f <- Y ~ x
  # baze pro x a bety
  basis.x <- list("x" = basis1) # vyhlazene data
  basis.b <- list("x" = basis2)
  # vstupni data do modelu
  ldata <- list("df" = dataf, "x" = x.train)
  # binomicky model ... model logisticke regrese
  model.glm <- fregre.glm(f, family = binomial(), data = ldata,
                          basis.x = basis.x, basis.b = basis.b)
  
  # presnost na trenovacich datech
  predictions.train <- predict(model.glm, newx = ldata)
  predictions.train <- data.frame(Y.pred = ifelse(predictions.train < 1/2, 0, 1))
  presnost.train <- table(Y.train, predictions.train$Y.pred) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  newldata = list("df" = as.data.frame(Y.test), "x" = fdata(X.test))
  predictions.test <- predict(model.glm, newx = newldata)
  predictions.test <- data.frame(Y.pred = ifelse(predictions.test < 1/2, 0, 1))
  presnost.test <- table(Y.test, predictions.test$Y.pred) |>
    prop.table() |> diag() |> sum()
  
  # vlozime do matice
  pred.baz[as.character(i), ] <- 1 - c(presnost.train, presnost.test)
} 

pred.baz <- as.data.frame(pred.baz)
pred.baz$n.basis <- n.basis
```

Znázorněme si průběh obou typů chybovostí v grafu v závislosti na počtu bazických funkcí.


```r
n.basis.beta.opt <- pred.baz$n.basis[which.min(pred.baz$Err.test)]

pred.baz |> ggplot(aes(x = n.basis, y = Err.test)) + 
  geom_line(linetype = 'dashed', colour = 'black') + 
  geom_line(aes(x = n.basis, y = Err.train), colour = 'deepskyblue3', 
            linetype = 'dashed', linewidth = 0.5) + 
  geom_point(size = 1.5) + 
  geom_point(aes(x = n.basis, y = Err.train), colour = 'deepskyblue3', 
             size = 1.5) + 
  geom_point(aes(x = n.basis.beta.opt, y = min(pred.baz$Err.test)),
             colour = 'red', size = 2) +
  theme_bw() + 
  labs(x = bquote(paste(n[basis], ' ;   ', 
                        n[optimal] == .(n.basis.beta.opt))),
        y = 'Chybovost')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-140-1.png" alt="Závislost testovací a trénovací chybovosti na počtu bázových funkcí pro $\beta$. Červeným bodem je znázorněn optimální počet $n_{optimal}$ zvolený jako minimum testovací chybovosti, černou čarou je vykreslena testovací a modrou přerušovanou čarou je vykreslen průběh trénovací chybovosti." width="672" />
<p class="caption">(\#fig:unnamed-chunk-140)Závislost testovací a trénovací chybovosti na počtu bázových funkcí pro $\beta$. Červeným bodem je znázorněn optimální počet $n_{optimal}$ zvolený jako minimum testovací chybovosti, černou čarou je vykreslena testovací a modrou přerušovanou čarou je vykreslen průběh trénovací chybovosti.</p>
</div>

Vidíme, že s rostoucím počtem bází pro $\beta(t)$ má trénovací chybovost (modrá čára) tendenci klesat a tedy bychom na jejím základě volili velké hodnoty $n_{basis}$.
Naopak optimální volbou na základě testovací chybovosti je $n$ rovno 5, tedy výrazně menší hodnota než 30.
Naopak s rostoucím $n$ roste testovací chyvost, což ukazuje na přeučení modelu.

Z výše uvedených důvodů pro určení optimálního počtu bazických funkcí pro $\beta(t)$ využijeme 10-ti násobnou cross-validaci.
Jako maximální počet uvažovaných bazických funkcí bereme 25, neboť jak jsme viděli výše, nad touto hodnotou dochází již k přeučení modelu.


```r
### 10-fold cross-validation
n.basis.max <- 25
n.basis <- 4:n.basis.max
k_cv <- 10 # k-fold CV
# rozdelime trenovaci data na k casti
folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)
## prvky, ktere se behem cyklu nemeni
# body, ve kterych jsou funkce vyhodnoceny
tt <- x.train[["argvals"]]
rangeval <- range(tt)
# B-spline baze 
basis1 <- create.bspline.basis(rangeval = range(tt), nbasis = nbasis.x)
# vztah
f <- Y ~ x
# baze pro x
basis.x <- list("x" = basis1)
# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro danou cast trenovaci mnoziny
# v radcich budou hodnoty pro dany pocet bazi
CV.results <- matrix(NA, nrow = length(n.basis), ncol = k_cv, 
                     dimnames = list(n.basis, 1:k_cv))
```

Nyní již máme vše připravené pro spočítání chybovosti na každé z deseti podmnožin trénovací množiny.
Následně určíme průměr a jako optimální $n$ vezmeme argument minima validační chybovosti.


```r
for (index in 1:k_cv) {
  # definujeme danou indexovou mnozinu
  fold <- folds[[index]]
    
  x.train.cv <- subset(X.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
    fdata()
  y.train.cv <- subset(Y.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
    as.numeric()
  
  x.test.cv <- subset(X.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
    fdata()
  y.test.cv <- subset(Y.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
    as.numeric()
  
  dataf <- as.data.frame(y.train.cv) 
  colnames(dataf) <- "Y"
  
  for (i in n.basis) {
    # baze pro bety
    basis2 <- create.bspline.basis(rangeval = rangeval, nbasis = i)
    
    basis.b <- list("x" = basis2)
    # vstupni data do modelu
    ldata <- list("df" = dataf, "x" = x.train.cv)
    # binomicky model ... model logisticke regrese
    model.glm <- fregre.glm(f, family = binomial(), data = ldata,
                            basis.x = basis.x, basis.b = basis.b)
    
    # presnost na validacni casti 
    newldata = list("df" = as.data.frame(y.test.cv), "x" = x.test.cv)
    predictions.valid <- predict(model.glm, newx = newldata)
    predictions.valid <- data.frame(Y.pred = ifelse(predictions.valid < 1/2, 0, 1))
    presnost.valid <- table(y.test.cv, predictions.valid$Y.pred) |>
      prop.table() |> diag() |> sum()
    
    # vlozime do matice
    CV.results[as.character(i), as.character(index)] <- presnost.valid
  } 
}

# spocitame prumerne presnosti pro jednotliva n pres folds
CV.results <- apply(CV.results, 1, mean)
n.basis.opt <- n.basis[which.max(CV.results)]
presnost.opt.cv <- max(CV.results)
```

Vykresleme si ještě průběh validační chybovosti i se zvýrazněnou optimální hodnotou $n_{optimal}$ rovnou 7 s validační chybovostí 0.0574.


```r
CV.results <- data.frame(n.basis = n.basis, CV = CV.results)
CV.results |> ggplot(aes(x = n.basis, y = 1 - CV)) + 
  geom_line(linetype = 'dashed', colour = 'grey') + 
  geom_point(size = 1.5) + 
  geom_point(aes(x = n.basis.opt, y = 1 - presnost.opt.cv), colour = 'red', size = 2) +
  theme_bw() + 
  labs(x = bquote(paste(n[basis], ' ;   ', 
                        n[optimal] == .(n.basis.opt))),
       y = 'Validační chybovost') + 
  scale_x_continuous(breaks = n.basis) + 
  theme(panel.grid.minor = element_blank())
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-143-1.png" alt="Závislost validační chybovosti na hodnotě $n_{basis}$, tedy na počtu bází." width="672" />
<p class="caption">(\#fig:unnamed-chunk-143)Závislost validační chybovosti na hodnotě $n_{basis}$, tedy na počtu bází.</p>
</div>

Nyní již tedy můžeme definovat finální model pomocí funkcionální logistické regrese, přičemž bázi pro $\beta(t)$ volíme B-splinovou bázi s 7 bázemi.


```r
# optimalni model
basis2 <- create.bspline.basis(rangeval = range(tt), nbasis = n.basis.opt)
f <- Y ~ x
# baze pro x a bety
basis.x <- list("x" = basis1) 
basis.b <- list("x" = basis2)
# vstupni data do modelu
dataf <- as.data.frame(y.train) 
colnames(dataf) <- "Y"
ldata <- list("df" = dataf, "x" = x.train)
# binomicky model ... model logisticke regrese
model.glm <- fregre.glm(f, family = binomial(), data = ldata,
                        basis.x = basis.x, basis.b = basis.b,
                        maxit = 1000, epsilon = 1e-2)
```

```
## Warning: glm.fit: algorithm did not converge
```

```
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```r
# presnost na trenovacich datech
predictions.train <- predict(model.glm, newx = ldata)
predictions.train <- data.frame(Y.pred = ifelse(predictions.train < 1/2, 0, 1))
presnost.train <- table(Y.train, predictions.train$Y.pred) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
newldata = list("df" = as.data.frame(Y.test), "x" = fdata(X.test))
predictions.test <- predict(model.glm, newx = newldata)
predictions.test <- data.frame(Y.pred = ifelse(predictions.test < 1/2, 0, 1))
presnost.test <- table(Y.test, predictions.test$Y.pred) |>
  prop.table() |> diag() |> sum()
```

Spočítali jsme trénovací chybovost (rovna 0 %) i testovací chybovost (rovna 7.69 %).
Pro lepší představu si ještě můžeme vykreslit hodnoty odhadnutých pravděpodobností příslušnosti do klasifikační třídy $Y = 1$ na trénovacích datech v závislosti na hodnotách lineárního prediktoru.


```r
data.frame(
  linear.predictor = model.glm$linear.predictors,
  response = model.glm$fitted.values,
  Y = factor(y.train)
) |> ggplot(aes(x = linear.predictor, y = response, colour = Y)) + 
  geom_point(size = 1.5) + 
  scale_color_discrete(labels = c("malý", "velký")) + 
  geom_abline(aes(slope = 0, intercept = 0.5), linetype = 'dashed') + 
  theme_bw() + 
  labs(x = 'Lineární prediktor',
       y = 'Odhadnuté pravděpodobnosti Pr(Y = 1|X = x)',
       colour = 'Obsah tuku') 
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-145-1.png" alt="Závoslost odhadnutých pravděpodobností na hodnotách lineárního prediktoru. Barevně jsou odlišeny body podle příslušnosti do klasifikační třídy." width="672" />
<p class="caption">(\#fig:unnamed-chunk-145)Závoslost odhadnutých pravděpodobností na hodnotách lineárního prediktoru. Barevně jsou odlišeny body podle příslušnosti do klasifikační třídy.</p>
</div>

Můžeme si ještě pro informaci zobrazit průběh odhadnuté parametrické funkce $\beta(t)$.


```r
t.seq <- seq(min(t), max(t), length = 1001)
beta.seq <- eval.fd(evalarg = t.seq, fdobj = model.glm$beta.l$x)

data.frame(t = t.seq, beta = beta.seq) |> 
  ggplot(aes(t, beta)) + 
  geom_abline(aes(slope = 0, intercept = 0), linetype = 'dashed', 
              linewidth = 0.5, colour = 'grey') +
  geom_line() + 
  theme_bw() +
  labs(x = expression(x[1]),
       y = expression(widehat(beta)(t))) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-146-1.png" alt="Průběh odhadu parametrické funkce $\beta(t), t \in [850, 1050]$." width="672" />
<p class="caption">(\#fig:unnamed-chunk-146)Průběh odhadu parametrické funkce $\beta(t), t \in [850, 1050]$.</p>
</div>

Vidíme, že hodnoty funkce $\hat\beta(t)$ se drží kolem nuly pro časy $t$ z prostředka  a začátku intervalu $[850, 1050]$, zatímco pro pozdějsí časy jsou hodnoty vyšší.
To implikuje rozdílnost funkcí z klasifikačních tříd na začátku a konci intervalu, zatímco uprostřed intervalu jsou funkce velmi podobné.

Výsledky opět přidáme do souhrnné tabulky.


```r
Res <- data.frame(model = 'LR functional', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

#### Logistická regrese s analýzou hlavních komponent

Abychom mohli sesrojit tento klasifikátor, potřebujeme provést funkcionální analýzu hlavních komponent, určit vhodný počet komponent a spočítat hodnoty skórů pro testovací data.
To jsme již provedli v části u lineární diskriminační analýzy, proto využijeme tyto výsledky v následující části.





Můžeme tedy rovnou sestrojit model logistické regrese pomocí funkce `glm(, family = binomial)`.


```r
# model
clf.LR <- glm(Y ~  ., data = data.PCA.train, family = binomial)
```

```
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```r
# presnost na trenovacich datech
predictions.train <- predict(clf.LR, newdata = data.PCA.train, type = 'response')
predictions.train <- ifelse(predictions.train > 0.5, 1, 0)
presnost.train <- table(data.PCA.train$Y, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.LR, newdata = data.PCA.test, type = 'response')
predictions.test <- ifelse(predictions.test > 0.5, 1, 0)
presnost.test <- table(data.PCA.test$Y, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Spočítali jsme tedy přesnost klasifikátoru na trénovacích (99.33 %) i na testovacích datech (95.38 %).

Pro grafické znázornění metody můžeme zaznačit dělící hranici do grafu skórů prvních dvou hlavních komponent.
Tuto hranici spočítáme na husté síti bodů a zobrazíme ji pomocí funkce `geom_contour()` stejně jako v případě LDA i QDA.




```r
nd <- nd |> mutate(prd = as.numeric(predict(clf.LR, newdata = nd,
                                            type = 'response')))
nd$prd <- ifelse(nd$prd > 0.5, 1, 0)

data.PCA.train |> ggplot(aes(x = V1, y = V2, colour = Y)) +
  geom_point(size = 1.5) + 
  labs(x = paste('1. hlavní komponenta (vysvětlená variabilita', 
                 round(100 * data.PCA$varprop[1], 2), '%)'),
       y = paste('2. hlavní komponenta (', 
                 round(100 * data.PCA$varprop[2], 2), '%)'),
       colour = 'Obsah tuku') +
  scale_colour_discrete(labels = c("malý", "velký")) +
  theme_bw() +
  geom_contour(data = nd, aes(x = V1, y = V2, z = prd), colour = 'black')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-152-1.png" alt="Skóre prvních dvou hlavních komponent, barevně odlišené podle příslušnosti do klasifikační třídy. Černě je vyznačena dělící hranice (přímka v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí logistické regrese." width="672" />
<p class="caption">(\#fig:unnamed-chunk-152)Skóre prvních dvou hlavních komponent, barevně odlišené podle příslušnosti do klasifikační třídy. Černě je vyznačena dělící hranice (přímka v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí logistické regrese.</p>
</div>

Všimněme si, že dělící hranicí mezi klasifikačními třídami je nyní přímka jako v případě LDA.

Nakonec ještě doplníme chybovosti do souhrnné tabulky.


```r
Res <- data.frame(model = 'LR score', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

### Rozhodovací stromy

V této části se podíváme na velmi odlišný přístup k sestrojení klasifikátoru, než byly například LDA či logistická regrese.
Rozhodovací stromy jsou velmi oblíbeným nástrojem ke klasifikaci, avšak jako v případě některých předchozích metod nejsou přímo určeny pro funkcionální data.
Existují však postupy, jak funkcionální objekty převést na mnohorozměrné a následně na ně aplikovat algoritmus rozhodovacích stromů.
Můžeme uvažovat následující postupy:

-   algoritmus sestrojený na bázových koeficientech,

-   využití skórů hlavních komponent,

-   použít diskretizaci intervalu a vyhodnotit funkci jen na nějaké konečné síti bodů.

My se nejprve zaměříme na diskretizaci intervalu a následně porovnáme výsledky se zbylými dvěma přístupy k sestrojení rozhodovacího stromu.

#### Diskretizace intervalu

Nejprve si musíme definovat body z intervalu $I = [850, 1050]$, ve kterých funkce vyhodnotíme.
Následně vytvoříme objekt, ve kterém budou řádky představovat jednotlivé (diskretizované) funkce a sloupce časy.
Nakonec připojíme sloupec $Y$ s informací o příslušnosti do klasifikační třídy a totéž zopakujeme i pro testovací data.


```r
# posloupnost bodu, ve kterych funkce vyhodnotime
t.seq <- seq(min(t), max(t), length = 101)
   
grid.data <- eval.fd(fdobj = X.train, evalarg = t.seq)
grid.data <- as.data.frame(t(grid.data)) # transpozice kvuli funkcim v radku
grid.data$Y <- Y.train |> factor()

grid.data.test <- eval.fd(fdobj = X.test, evalarg = t.seq)
grid.data.test <- as.data.frame(t(grid.data.test))
grid.data.test$Y <- Y.test |> factor()
```

Nyní mážeme sestrojit rozhodovací strom, ve kterém budou jakožto prediktory vystupovat všechny časy z vektoru `t.seq`.
Tato klasifikační není náchylná na multikolinearitu, tudíž se jí nemusíme zabývat.
Jako metriku zvolíme přesnost.


```r
# sestrojeni modelu
clf.tree <- train(Y ~ ., data = grid.data, 
                 method = "rpart", 
                 trControl = trainControl(method = "CV", number = 10),
                 metric = "Accuracy")

# presnost na trenovacich datech
predictions.train <- predict(clf.tree, newdata = grid.data)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.tree, newdata = grid.data.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost klasifikátoru na testovacích datech je tedy 98.46 % a na trénovacích datech 99.33 %.

Graficky si rozhodovací strom můžeme vykreslit pomocí funkce `fancyRpartPlot()`.
Nastavíme barvy uzlů tak, aby reflektovaly předchozí barevné odlišení.
Jedná se o neprořezaný strom.


```r
colnames(grid.data) <- c(paste0('time:', t.seq), 'Y')
fancyRpartPlot(rpart(Y ~ ., data = grid.data, method = "class"),
               sub = '', palettes = c('Reds', 'Blues')) 
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-156-1.png" alt="Grafické znázornění neprořezaného rozhodovacího stromu. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0." width="672" />
<p class="caption">(\#fig:unnamed-chunk-156)Grafické znázornění neprořezaného rozhodovacího stromu. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0.</p>
</div>

Můžeme si také vykreslit již prořezaný finální rozhodovací strom.


```r
rpart.plot::rpart.plot(clf.tree$finalModel, # finalni model ... prorezany strom
                       extra = 104, # zobrazeni pozadovanych informaci
                       box.palette = c('Reds', 'Blues'),
                       branch.lty = 3, # dotted branch lines
                       shadow.col = 0, # shadows under the node boxes
                       nn = FALSE, 
                       under = FALSE, 
                       digits = 2)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-157-1.png" alt="Finální prořezaný rozhodovací strom." width="672" />
<p class="caption">(\#fig:unnamed-chunk-157)Finální prořezaný rozhodovací strom.</p>
</div>

Nakonec opět přidejme trénovací a testovací chybovost do souhrnné tabulky.


```r
Res <- data.frame(model = 'Tree - diskr.', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

#### Skóre hlavních komponent

Další možností pro sestrojení rozhodovacího stromu je použít skóre hlavních komponent.
Jelikož jsme již skóre počítali pro předchozí klasifikační metody, využijeme těchto poznatků a sestrojíme rozhodovací strom na skórech prvních 2 hlavních komponent.


```r
# sestrojeni modelu
clf.tree.PCA <- train(Y ~ ., data = data.PCA.train, 
                 method = "rpart", 
                 trControl = trainControl(method = "CV", number = 10),
                 metric = "Accuracy")

# presnost na trenovacich datech
predictions.train <- predict(clf.tree.PCA, newdata = data.PCA.train)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.tree.PCA, newdata = data.PCA.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost rozhodovacího stromu na testovacích datech je tedy 93.85 % a na trénovacích datech 99.33 %.

Graficky si rozhodovací strom sestrojený na skórech hlavních komponent můžeme vykreslit pomocí funkce `fancyRpartPlot()`.
Nastavíme barvy uzlů tak, aby reflektovaly předchozí barevné odlišení.
Jedná se o neprořezaný strom.


```r
fancyRpartPlot(rpart(Y ~ ., data = data.PCA.train, method = "class"),
               sub = '', palettes = c('Reds', 'Blues')) 
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-160-1.png" alt="Grafické znázornění neprořezaného rozhodovacího stromu sestrojeného na skórech hlavních komponent. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0." width="672" />
<p class="caption">(\#fig:unnamed-chunk-160)Grafické znázornění neprořezaného rozhodovacího stromu sestrojeného na skórech hlavních komponent. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0.</p>
</div>

Můžeme si také vykreslit již prořezaný finální rozhodovací strom.


```r
rpart.plot::rpart.plot(clf.tree.PCA$finalModel, # finalni model 
                       extra = 104, # zobrazeni pozadovanych informaci
                       box.palette = c('Reds', 'Blues'),
                       branch.lty = 3, # dotted branch lines
                       shadow.col = 0, # shadows under the node boxes
                       nn = FALSE, under = FALSE, digits = 2)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-161-1.png" alt="Finální prořezaný rozhodovací strom." width="672" />
<p class="caption">(\#fig:unnamed-chunk-161)Finální prořezaný rozhodovací strom.</p>
</div>

Nakonec opět přidejme trénovací a testovací chybovost do souhrnné tabulky.


```r
Res <- data.frame(model = 'Tree - score', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

#### Bázové koeficienty

Poslední možností, kterou využijeme pro sestrojení rozhodovacího stromu, je použití koeficientů ve vyjádření funkcí v B-splinové bázi.

Nejprve si definujme potřebné datové soubory s koeficienty.


```r
# trenovaci dataset
data.Bbasis.train <- t(X.train$coefs) |> as.data.frame()
data.Bbasis.train$Y <- factor(Y.train)

# testovaci dataset
data.Bbasis.test <- t(X.test$coefs) |> as.data.frame()
data.Bbasis.test$Y <- factor(Y.test)
```

Nyní již můžeme sestrojit klasifikátor.


```r
# sestrojeni modelu
clf.tree.Bbasis <- train(Y ~ ., data = data.Bbasis.train, 
                 method = "rpart", 
                 trControl = trainControl(method = "CV", number = 10),
                 metric = "Accuracy")

# presnost na trenovacich datech
predictions.train <- predict(clf.tree.Bbasis, newdata = data.Bbasis.train)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.tree.Bbasis, newdata = data.Bbasis.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost rozhodovacího stromu na trénovacích datech je tedy 99.33 % a na testovacích datech 98.46 %.

Graficky si rozhodovací strom sestrojený na koeficientech B-splinového vyjádření můžeme vykreslit pomocí funkce `fancyRpartPlot()`.
Nastavíme barvy uzlů tak, aby reflektovaly předchozí barevné odlišení.
Jedná se o neprořezaný strom.


```r
fancyRpartPlot(rpart(Y ~ ., data = data.Bbasis.train, method = "class"),
               sub = '', palettes = c('Reds', 'Blues')) 
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-165-1.png" alt="Grafické znázornění neprořezaného rozhodovacího stromu sestrojeného na bázových koeficientech. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0." width="672" />
<p class="caption">(\#fig:unnamed-chunk-165)Grafické znázornění neprořezaného rozhodovacího stromu sestrojeného na bázových koeficientech. Modrými odstíny jsou vykresleny uzly příslušející klasifikační třídě 1 a červenými odstíny třídě 0.</p>
</div>

Můžeme si také vykreslit již prořezaný finální rozhodovací strom.


```r
rpart.plot::rpart.plot(clf.tree.Bbasis$finalModel, # finalni model 
                       extra = 104, # zobrazeni pozadovanych informaci
                       box.palette = c('Reds', 'Blues'),
                       branch.lty = 3, # dotted branch lines
                       shadow.col = 0, # shadows under the node boxes
                       nn = FALSE, under = FALSE, digits = 2)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-166-1.png" alt="Finální prořezaný rozhodovací strom." width="672" />
<p class="caption">(\#fig:unnamed-chunk-166)Finální prořezaný rozhodovací strom.</p>
</div>

Nakonec opět přidejme trénovací a testovací chybovost do souhrnné tabulky.


```r
Res <- data.frame(model = 'Tree - Bbasis', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

### Náhodné lesy

Klasifikátor sestrojený pomocí metody náhodných lesů spočívá v sestrojení několika jednotlivých rozhodovacích stromů, které se následně zkombinují a vytvoří společný klasifikátor (společným "hlasováním").

Tak jako v případě rozhodovacích stromů máme několik možností na to, jaká data (konečně-rozměrná) použijeme pro sestrojení modelu.
Budeme opět uvažovat výše diskutované tři přístupy.
Datové soubory s příslušnými veličinami pro všechny tři přístupy již máme připravené z minulé sekce, proto můžeme přímo sestrojit dané modely, spočítat charakteristiky daného klasifikátoru a přidat výsledky do souhrnné tabulky.

#### Diskretizace intervalu

V prvním případě využíváme vyhodnocení funkcí na dané síti bodů intervalu $I = [850, 1050]$.




```r
# sestrojeni modelu
clf.RF <- randomForest(Y ~ ., data = grid.data, 
                       ntree = 500, # pocet stromu
                       importance = TRUE,
                       nodesize = 5)

# presnost na trenovacich datech
predictions.train <- predict(clf.RF, newdata = grid.data)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.RF, newdata = grid.data.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost náhodného lesu na trénovacích datech je tedy 100 % a na testovacích datech 95.38 %.


```r
Res <- data.frame(model = 'RForest - diskr', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

#### Skóre hlavních komponent

V tomto případě využijeme skóre prvních $p =$ 2 hlavních komponent.


```r
# sestrojeni modelu
clf.RF.PCA <- randomForest(Y ~ ., data = data.PCA.train, 
                           ntree = 500, # pocet stromu
                           importance = TRUE,
                           nodesize = 5)

# presnost na trenovacich datech
predictions.train <- predict(clf.RF.PCA, newdata = data.PCA.train)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.RF.PCA, newdata = data.PCA.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost rozhodovacího stromu na trénovacích datech je tedy 99.33 % a na testovacích datech 95.38 %.


```r
Res <- data.frame(model = 'RForest - score', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

#### Bázové koeficienty

Nakonec použijeme vyjádření funkcí pomocí B-splinové báze.


```r
# sestrojeni modelu
clf.RF.Bbasis <- randomForest(Y ~ ., data = data.Bbasis.train, 
                              ntree = 500, # pocet stromu
                              importance = TRUE,
                              nodesize = 5)

# presnost na trenovacich datech
predictions.train <- predict(clf.RF.Bbasis, newdata = data.Bbasis.train)
presnost.train <- table(Y.train, predictions.train) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test <- predict(clf.RF.Bbasis, newdata = data.Bbasis.test)
presnost.test <- table(Y.test, predictions.test) |>
  prop.table() |> diag() |> sum()
```

Přesnost tohoto klasifikátoru na trénovacích datech je 100 % a na testovacích datech 96.92 %.


```r
Res <- data.frame(model = 'RForest - Bbasis', 
                  Err.train = 1 - presnost.train,
                  Err.test = 1 - presnost.test)

RESULTS <- rbind(RESULTS, Res)
```

### Support Vector Machines

Nyní se podívejme na klasifikaci našich nasimulovaných křivek pomocí metody podpůrných vektorů (ang. Support Vector Machines, SVM).
Výhodou této klasifikační metody je její výpočetní nenáročnost, neboť pro definici hraniční křivky mezi třídami využívá pouze několik (často málo) pozorování.

Hlavní výhodou SVM je použití tzv.
*jádrového triku* (kernel trick), pomocí kterého nahradíme obyčejný skalární součin jiným skalárním součinem transformovaných dat, aniž bychom tuto transformaci museli přímo definovat.
Tím dostaneme obecně nelineární dělící hranici mezi klasifikačními třídami.
*Jádro* (jádrová funkce, ang. kernel, kernel function) $K$ je taková funkce, která splňuje

$$
K(x_i, x_j) = \langle \phi(x_i), \phi(x_j) \rangle_{\mathcal H}, 
$$
kde $\phi$ je nějaká (neznámá) transformace (ang. feature map), $\mathcal H$ je Hilbertův prostor a $\langle \cdot, \cdot \rangle_{\mathcal H}$ je nějaký skalární součin na tomto Hilbertově prostoru.

Nejčastěji se v praxi volí tři typy jádrových funkcí:

-   lineární jádro -- $K(x_i, x_j) = \langle x_i, x_j \rangle$,
-   polynomiální jádro -- $K(x_i, x_j) = \big(\alpha_0 + \gamma \langle x_i, x_j \rangle \big)^d$,
-   radiální (gaussovské) jádro -- $\displaystyle{K(x_i, x_j) = \text e^{-\gamma \|x_i - x_j \|^2}}$.

U všech výše zmíněných jader musíme zvolit konstantu $C > 0$, která udává míru penalizace za překročení dělící hranice mezi třídami (ang. inverse regularization parameter).
S rostoucí hodnotou $C$ bude metoda více penalizovat špatně klasifikovaná data a méně tvar hranice, naopak pro malé hodnoty $C$ metoda nedává takový význam špatně klasifikovaným datům, ale zaměřuje se více na penalizaci tvaru hranice.
Tato konstanta $C$ se defaultně volí rovna 1, můžeme ji určit i přímo například pomocí cross-validace.

Využitím cross-validace můžeme také určit optimální hodnoty ostatních hyperparametrů, které nyní závisí na naší volbě jádrové funkce.
V případě lineárního jádra nevolíme žádný další parametr kromě konstanty $C$, u polynomiálního jádra musíme určit hodnoty hyperparametrů $\alpha_0, \gamma \text{ a } d$, jejichž defaultní hodnoty v `R` jsou postupně $\alpha_0^{default} = 0, \gamma^{default} = \frac{1}{dim(\texttt{data})} \text{ a } d^{default} = 3$.
Při volbě radiálního jádra máme pouze jeden další hyperparametr $\gamma$, jehož defaultní hodnota v `R` je totožná jako u polynomiálního jádra.
Opět bychom mohli hodnoty hyperparametrů určit jako optimální pro naše data, avšak vzhledem k relativní výpočetní náročnosti necháme hodnoty příslušných hyperparametrů na jejich defaultních hodnotách.

V případě funkcionálních dat máme několik možností, jak použít metodu SVM.
Nejjednodušší variantou je použít tuto klasifikační metodu přímo na diskretizovanou funkci (sekce \@ref(diskrA3b)).
Další možností je opět využít skóre hlavních komponent a klasifikovat křivky pomocí jejich reprezentace \@ref(PCA-SVMA3b).
Další přímočarou variantou je využít vyjádření křivek pomocí B-splinové báze a klasifikovat křivky na základě koeficientů jejich vyjádření v této bázi (sekce \@ref(basis-SVMA3b)).

Složitější úvahou můžeme dospět k několika dalším možnostem, které využívají funkcionální podstatu dat.
Jednak můžeme místo klasifikace původní křivky využít její derivaci (případně druhou derivaci, třetí, ...), druhak můžeme využít projekce funkcí na podprostor generovaný, např. B-splinovými, funkcemi (sekce \@ref(projection-SVMA3b)). Poslední metoda, kterou použijeme pro klasifikaci funkcionálních dat, spočívá v kombinaci projekce na určitý podprostor generovaný funkcemi (Reproducing Kernel Hilbert Space, RKHS) a klasifikace příslušné reprezentace. Tato metoda využívá kromě klasického SVM i SVM pro regresi, více uvádíme v sekci RKHS + SVM \@ref(RKHS-SVMA3b).

#### Diskretizace intervalu {#diskrA3b}

Začněme nejprve aplikací metody podpůrných vektorů přímo na diskretizovaná data (vyhodnocení funkce na dané síti bodů na intervalu $I = [850, 1050]$), přičemž budeme uvažovat všech tři výše zmíněné jádrové funkce.


```r
# set norm equal to one
norms <- c()
for (i in 1:dim(XXder$coefs)[2]) {
  norms <- c(norms, as.numeric(1 / norm.fd(XXder[i])))
  }
XXfd_norm_der <- XXder 
XXfd_norm_der$coefs <- XXfd_norm_der$coefs * matrix(norms, 
                                            ncol = dim(XXder$coefs)[2],
                                            nrow = dim(XXder$coefs)[1],
                                            byrow = T)

# rozdeleni na testovaci a trenovaci cast
X.train_norm <- subset(XXfd_norm_der, split == TRUE)
X.test_norm <- subset(XXfd_norm_der, split == FALSE)

Y.train_norm <- subset(Y, split == TRUE)
Y.test_norm <- subset(Y, split == FALSE)

grid.data <- eval.fd(fdobj = X.train_norm, evalarg = t.seq)
grid.data <- as.data.frame(t(grid.data)) 
grid.data$Y <- Y.train_norm |> factor()

grid.data.test <- eval.fd(fdobj = X.test_norm, evalarg = t.seq)
grid.data.test <- as.data.frame(t(grid.data.test))
grid.data.test$Y <- Y.test_norm |> factor()
```

Nyní se pokusme, na rozdíl od postupu v předchozích kapitolách, hyperparametry klasifikátorů odhadnout z dat pomocí 10-násobné cross-validace. Vzhledem k tomu, že každé jádro má ve své definici jiné hyperparametry, budeme ke každé jádrové funkci přistupovat zvlášť. Nicméně hyperparametr $C$ vystupuje u všech jádrových funkcí, přičemž ale připouštíme, že se může jeho optimální hodnota mezi jádry lišit.

U všech třech jader projdeme hodnoty hyperparametru $C$ v intervalu $[10^{-3}, 10^{3}]$, přičemž u jádra polynomiálního zafixujeme hyperparametr $p$ na hodnotě 3, neboť pro jiné celočíselné hodnoty metoda nedává zdaleka tak dobré výsledky. Naopak pro radiální jádro využijeme k volbě optimální hodnoty hyperparametru $\gamma$ opět 10-násobnou CV, přičemž uvažujeme hodnoty v intervalu $[10^{-3}, 10^{2}]$. Zvolíme `coef0` $= 1$. 


```r
set.seed(42)

k_cv <- 10 #  k-fold CV

# rozdelime trenovaci data na k casti
folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)

# ktere hodnoty gamma chceme uvazovat
gamma.cv <- 10^seq(-3, 2, length = 15)
C.cv <- 10^seq(-3, 3, length = 20)
p.cv <- 3
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
    cv_sample <- 1:dim(grid.data)[1] %in% fold
    
    data.grid.train.cv <- as.data.frame(grid.data[cv_sample, ])
    data.grid.test.cv <- as.data.frame(grid.data[!cv_sample, ])
    
    ## LINEARNI JADRO
    # sestrojeni modelu
    clf.SVM.l <- svm(Y ~ ., data = data.grid.train.cv,
                     type = 'C-classification',
                     scale = TRUE,
                     cost = C,
                     kernel = 'linear')
    
    # presnost na validacnich datech
    predictions.test.l <- predict(clf.SVM.l, newdata = data.grid.test.cv)
    presnost.test.l <- table(data.grid.test.cv$Y, predictions.test.l) |>
      prop.table() |> diag() |> sum()
    
    # presnosti vlozime na pozice pro dane C a fold
    CV.results$SVM.l[(1:length(C.cv))[C.cv == C], 
                     index_cv] <- presnost.test.l
    
    ## POLYNOMIALNI JADRO
    for (p in p.cv) {
      # sestrojeni modelu
      clf.SVM.p <- svm(Y ~ ., data = data.grid.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       coef0 = coef0,
                       degree = p,
                       kernel = 'polynomial')
      
      # presnost na validacnich datech
      predictions.test.p <- predict(clf.SVM.p, newdata = data.grid.test.cv)
      presnost.test.p <- table(data.grid.test.cv$Y, predictions.test.p) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane C, p a fold
      CV.results$SVM.p[(1:length(C.cv))[C.cv == C], 
                       (1:length(p.cv))[p.cv == p],
                       index_cv] <- presnost.test.p
    }
        
    ## RADIALNI JADRO
    for (gamma in gamma.cv) {
      # sestrojeni modelu
      clf.SVM.r <- svm(Y ~ ., data = data.grid.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       gamma = gamma,
                       kernel = 'radial')
      
      # presnost na validacnich datech
      predictions.test.r <- predict(clf.SVM.r, newdata = data.grid.test.cv)
      presnost.test.r <- table(data.grid.test.cv$Y, predictions.test.r) |>
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

Podívejme se, jak dopadly optimální hodnoty. Pro *lineární jádro* máme optimální hodnotu $C$ rovnu 0.001, pro *polynomiální jádro* je $C$ rovno 0.0785 a pro *radiální jádro* máme dvě optimální hodnoty, pro $C$ je optimální hodnota 0.336 a pro $\gamma$ je to 0.0052. Validační přesnosti jsou postupně 0.9933333 pro lineární, 0.9933333 pro polynomiální a 0.9933333 pro radiální jádro.

Konečně můžeme sestrojit finální klasifikátory na celých trénovacích datech s hodnotami hyperparametrů určenými pomocí 10-násobné CV. Určíme také chybovosti na testovacích a také na trénovacích datech.


```r
# sestrojeni modelu
clf.SVM.l <- svm(Y ~ ., data = grid.data,
                 type = 'C-classification',
                 scale = TRUE,
                 cost = C.opt[1],
                 kernel = 'linear')

clf.SVM.p <- svm(Y ~ ., data = grid.data,
                 type = 'C-classification',
                 scale = TRUE,
                 cost = C.opt[2],
                 degree = p.opt,
                 coef0 = coef0,
                 kernel = 'polynomial')

clf.SVM.r <- svm(Y ~ ., data = grid.data,
                 type = 'C-classification',
                 scale = TRUE, 
                 cost = C.opt[3],
                 gamma = gamma.opt,
                 kernel = 'radial')

# presnost na trenovacich datech
predictions.train.l <- predict(clf.SVM.l, newdata = grid.data)
presnost.train.l <- table(Y.train, predictions.train.l) |>
  prop.table() |> diag() |> sum()

predictions.train.p <- predict(clf.SVM.p, newdata = grid.data)
presnost.train.p <- table(Y.train, predictions.train.p) |>
  prop.table() |> diag() |> sum()

predictions.train.r <- predict(clf.SVM.r, newdata = grid.data)
presnost.train.r <- table(Y.train, predictions.train.r) |>
  prop.table() |> diag() |> sum()

# presnost na testovacich datech
predictions.test.l <- predict(clf.SVM.l, newdata = grid.data.test)
presnost.test.l <- table(Y.test, predictions.test.l) |>
  prop.table() |> diag() |> sum()

predictions.test.p <- predict(clf.SVM.p, newdata = grid.data.test)
presnost.test.p <- table(Y.test, predictions.test.p) |>
  prop.table() |> diag() |> sum()

predictions.test.r <- predict(clf.SVM.r, newdata = grid.data.test)
presnost.test.r <- table(Y.test, predictions.test.r) |>
  prop.table() |> diag() |> sum()
```

Přesnost metody SVM na trénovacích datech je tedy 99.3333 % pro lineární jádro, 99.3333 % pro polynomiální jádro a 99.3333 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 98.4615 % pro lineární jádro, 93.8462 % pro polynomiální jádro a 98.4615 % pro radiální jádro.


```r
Res <- data.frame(model = c('SVM linear - diskr', 
                            'SVM poly - diskr', 
                            'SVM rbf - diskr'), 
                  Err.train = 1 - c(presnost.train.l, presnost.train.p, presnost.train.r),
                  Err.test = 1 - c(presnost.test.l, presnost.test.p, presnost.test.r))

RESULTS <- rbind(RESULTS, Res)
```

#### Skóre hlavních komponent {#PCA-SVMA3b}

V tomto případě využijeme skóre prvních $p =$ 2 hlavních komponent.

Nyní se pokusme, na rozdíl od postupu v předchozích kapitolách, hyperparametry klasifikátorů odhadnout z dat pomocí 10-násobné cross-validace. Vzhledem k tomu, že každé jádro má ve své definici jiné hyperparametry, budeme ke každé jádrové funkci přistupovat zvlášť. Nicméně hyperparametr $C$ vystupuje u všech jádrových funkcí, přičemž ale připouštíme, že se může jeho optimální hodnota mezi jádry lišit.

U všech třech jader projdeme hodnoty hyperparametru $C$ v intervalu $[10^{-3}, 10^{3}]$, přičemž u jádra polynomiálního zafixujeme hyperparametr $p$ na hodnotě 3, neboť pro jiné celočíselné hodnoty metoda nedává zdaleka tak dobré výsledky. Naopak pro radiální jádro využijeme k volbě optimální hodnoty hyperparametru $\gamma$ opět 10-násobnou CV, přičemž uvažujeme hodnoty v intervalu $[10^{-3}, 10^{2}]$. Zvolíme `coef0` $= 1$. 


```r
set.seed(42)

# ktere hodnoty gamma chceme uvazovat
gamma.cv <- 10^seq(-3, 2, length = 15)
C.cv <- 10^seq(-3, 3, length = 20)
p.cv <- 3
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

Podívejme se, jak dopadly optimální hodnoty. Pro *lineární jádro* máme optimální hodnotu $C$ rovnu 2.9764, pro *polynomiální jádro* je $C$ rovno 0.6952 a pro *radiální jádro* máme dvě optimální hodnoty, pro $C$ je optimální hodnota 1000 a pro $\gamma$ je to 1.6379. Validační přesnosti jsou postupně 0.9933333 pro lineární, 0.9933333 pro polynomiální a 0.9933333 pro radiální jádro.

Konečně můžeme sestrojit finální klasifikátory na celých trénovacích datech s hodnotami hyperparametrů určenými pomocí 10-násobné CV. Určíme také chybovosti na testovacích a také na trénovacích datech.


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

# presnost na trenovacich datech
predictions.train.l <- predict(clf.SVM.l.PCA, newdata = data.PCA.train)
presnost.train.l <- table(data.PCA.train$Y, predictions.train.l) |>
  prop.table() |> diag() |> sum()

predictions.train.p <- predict(clf.SVM.p.PCA, newdata = data.PCA.train)
presnost.train.p <- table(data.PCA.train$Y, predictions.train.p) |>
  prop.table() |> diag() |> sum()

predictions.train.r <- predict(clf.SVM.r.PCA, newdata = data.PCA.train)
presnost.train.r <- table(data.PCA.train$Y, predictions.train.r) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test.l <- predict(clf.SVM.l.PCA, newdata = data.PCA.test)
presnost.test.l <- table(data.PCA.test$Y, predictions.test.l) |>
  prop.table() |> diag() |> sum()

predictions.test.p <- predict(clf.SVM.p.PCA, newdata = data.PCA.test)
presnost.test.p <- table(data.PCA.test$Y, predictions.test.p) |>
  prop.table() |> diag() |> sum()

predictions.test.r <- predict(clf.SVM.r.PCA, newdata = data.PCA.test)
presnost.test.r <- table(data.PCA.test$Y, predictions.test.r) |>
  prop.table() |> diag() |> sum()
```

Přesnost metody SVM aplikované na skóre hlavních komponent na trénovacích datech je tedy 99.33 % pro lineární jádro, 99.33 % pro polynomiální jádro a 99.33 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 95.3846 % pro lineární jádro, 98.4615 % pro polynomiální jádro a 96.9231 % pro radiální jádro.

Pro grafické znázornění metody můžeme zaznačit dělící hranici do grafu skórů prvních dvou hlavních komponent.
Tuto hranici spočítáme na husté síti bodů a zobrazíme ji pomocí funkce `geom_contour()` stejně jako v předchozích případech, kdy jsme také vykreslovali klasifikační hranici.


```r
nd <- rbind(nd, nd, nd) |> mutate(
   prd = c(as.numeric(predict(clf.SVM.l.PCA, newdata = nd, type = 'response')),
           as.numeric(predict(clf.SVM.p.PCA, newdata = nd, type = 'response')),
           as.numeric(predict(clf.SVM.r.PCA, newdata = nd, type = 'response'))),
   kernel = rep(c('lineární', 'polynomiální', 'radiální'),
                each = length(as.numeric(predict(clf.SVM.l.PCA, 
                                                 newdata = nd,
                                                 type = 'response')))) |>
     as.factor())

data.PCA.train |> ggplot(aes(x = V1, y = V2, colour = Y)) +
 geom_point(size = 1.5) + 
 labs(x = paste('1. hlavní komponenta (vysvětlená variabilita', 
                round(100 * data.PCA$varprop[1], 2), '%)'),
      y = paste('2. hlavní komponenta (', 
                round(100 * data.PCA$varprop[2], 2), '%)'),
      colour = 'Obsah tuku', 
      linetype = 'Jádro') +
 scale_colour_discrete(labels = c("malý", "velký")) +
 theme_bw() +
 geom_contour(data = nd, aes(x = V1, y = V2, z = prd, linetype = kernel), 
              colour = 'black') 
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-183-1.png" alt="Skóre prvních dvou hlavních komponent, barevně odlišené podle příslušnosti do klasifikační třídy. Černě je vyznačena dělící hranice (přímka, resp. křivky v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí metody SVM." width="672" />
<p class="caption">(\#fig:unnamed-chunk-183)Skóre prvních dvou hlavních komponent, barevně odlišené podle příslušnosti do klasifikační třídy. Černě je vyznačena dělící hranice (přímka, resp. křivky v rovině prvních dvou hlavních komponent) mezi třídami sestrojená pomocí metody SVM.</p>
</div>


```r
Res <- data.frame(model = c('SVM linear - PCA', 
                            'SVM poly - PCA', 
                            'SVM rbf - PCA'), 
                  Err.train = 1 - c(presnost.train.l, presnost.train.p, presnost.train.r),
                  Err.test = 1 - c(presnost.test.l, presnost.test.p, presnost.test.r))

RESULTS <- rbind(RESULTS, Res)
```

#### Bázové koeficienty {#basis-SVMA3b}

Nakonec použijeme vyjádření funkcí pomocí B-splinové báze.
U všech třech jader projdeme hodnoty hyperparametru $C$ v intervalu $[10^{-3}, 10^{3}]$, přičemž u jádra polynomiálního zafixujeme hyperparametr $p$ na hodnotě 3, neboť pro jiné celočíselné hodnoty metoda nedává zdaleka tak dobré výsledky. Naopak pro radiální jádro využijeme k volbě optimální hodnoty hyperparametru $\gamma$ opět 10-násobnou CV, přičemž uvažujeme hodnoty v intervalu $[10^{-3}, 10^{2}]$. Zvolíme `coef0` $= 1$.


```r
set.seed(42)

# ktere hodnoty gamma chceme uvazovat
gamma.cv <- 10^seq(-3, 2, length = 15)
C.cv <- 10^seq(-3, 3, length = 20)
p.cv <- 3
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
    cv_sample <- 1:dim(grid.data)[1] %in% fold
    
    data.Bbasis.train.cv <- as.data.frame(data.Bbasis.train[cv_sample, ])
    data.Bbasis.test.cv <- as.data.frame(data.Bbasis.train[!cv_sample, ])
    
    ## LINEARNI JADRO
    # sestrojeni modelu
    clf.SVM.l <- svm(Y ~ ., data = data.Bbasis.train.cv,
                     type = 'C-classification',
                     scale = TRUE,
                     cost = C,
                     kernel = 'linear')
    
    # presnost na validacnich datech
    predictions.test.l <- predict(clf.SVM.l, newdata = data.Bbasis.test.cv)
    presnost.test.l <- table(data.Bbasis.test.cv$Y, predictions.test.l) |>
      prop.table() |> diag() |> sum()
    
    # presnosti vlozime na pozice pro dane C a fold
    CV.results$SVM.l[(1:length(C.cv))[C.cv == C], 
                     index_cv] <- presnost.test.l
    
    ## POLYNOMIALNI JADRO
    for (p in p.cv) {
      # sestrojeni modelu
      clf.SVM.p <- svm(Y ~ ., data = data.Bbasis.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       coef0 = coef0,
                       degree = p,
                       kernel = 'polynomial')
      
      # presnost na validacnich datech
      predictions.test.p <- predict(clf.SVM.p, 
                                    newdata = data.Bbasis.test.cv)
      presnost.test.p <- table(data.Bbasis.test.cv$Y, predictions.test.p) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane C, p a fold
      CV.results$SVM.p[(1:length(C.cv))[C.cv == C], 
                       (1:length(p.cv))[p.cv == p],
                       index_cv] <- presnost.test.p
    }
        
    ## RADIALNI JADRO
    for (gamma in gamma.cv) {
      # sestrojeni modelu
      clf.SVM.r <- svm(Y ~ ., data = data.Bbasis.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       gamma = gamma,
                       kernel = 'radial')
      
      # presnost na validacnich datech
      predictions.test.r <- predict(clf.SVM.r, 
                                    newdata = data.Bbasis.test.cv)
      presnost.test.r <- table(data.Bbasis.test.cv$Y, predictions.test.r) |>
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

Podívejme se, jak dopadly optimální hodnoty. Pro *lineární jádro* máme optimální hodnotu $C$ rovnu 0.0089, pro *polynomiální jádro* je $C$ rovno 0.1624 a pro *radiální jádro* máme dvě optimální hodnoty, pro $C$ je optimální hodnota 2.9764 a pro $\gamma$ je to 0.0118. Validační přesnosti jsou postupně 0.9933333 pro lineární, 0.9933333 pro polynomiální a 0.9933333 pro radiální jádro.

Konečně můžeme sestrojit finální klasifikátory na celých trénovacích datech s hodnotami hyperparametrů určenými pomocí 10-násobné CV. Určíme také chybovosti na testovacích a také na trénovacích datech.


```r
# sestrojeni modelu
clf.SVM.l.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                        type = 'C-classification',
                        scale = TRUE,
                        cost = C.opt[1],
                        kernel = 'linear')

clf.SVM.p.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                        type = 'C-classification',
                        scale = TRUE,
                        cost = C.opt[2],
                        degree = p.opt,
                        coef0 = coef0,
                        kernel = 'polynomial')

clf.SVM.r.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                        type = 'C-classification',
                        scale = TRUE,
                        cost = C.opt[3],
                        gamma = gamma.opt,
                        kernel = 'radial')

# presnost na trenovacich datech
predictions.train.l <- predict(clf.SVM.l.Bbasis, newdata = data.Bbasis.train)
presnost.train.l <- table(Y.train, predictions.train.l) |>
  prop.table() |> diag() |> sum()

predictions.train.p <- predict(clf.SVM.p.Bbasis, newdata = data.Bbasis.train)
presnost.train.p <- table(Y.train, predictions.train.p) |>
  prop.table() |> diag() |> sum()

predictions.train.r <- predict(clf.SVM.r.Bbasis, newdata = data.Bbasis.train)
presnost.train.r <- table(Y.train, predictions.train.r) |>
  prop.table() |> diag() |> sum()
  
# presnost na testovacich datech
predictions.test.l <- predict(clf.SVM.l.Bbasis, newdata = data.Bbasis.test)
presnost.test.l <- table(Y.test, predictions.test.l) |>
  prop.table() |> diag() |> sum()

predictions.test.p <- predict(clf.SVM.p.Bbasis, newdata = data.Bbasis.test)
presnost.test.p <- table(Y.test, predictions.test.p) |>
  prop.table() |> diag() |> sum()

predictions.test.r <- predict(clf.SVM.r.Bbasis, newdata = data.Bbasis.test)
presnost.test.r <- table(Y.test, predictions.test.r) |>
  prop.table() |> diag() |> sum()
```

Přesnost metody SVM aplikované na bázové koeficienty na trénovacích datech je tedy 99.33 % pro lineární jádro, 99.33 % pro polynomiální jádro a 99.33 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 92.3077 % pro lineární jádro, 92.3077 % pro polynomiální jádro a 90.7692 % pro radiální jádro.


```r
Res <- data.frame(model = c('SVM linear - Bbasis', 
                            'SVM poly - Bbasis', 
                            'SVM rbf - Bbasis'), 
                  Err.train = 1 - c(presnost.train.l, presnost.train.p, presnost.train.r),
                  Err.test = 1 - c(presnost.test.l, presnost.test.p, presnost.test.r))

RESULTS <- rbind(RESULTS, Res)
```

#### Projekce na B-splinovou bázi {#projection-SVMA3b}

Další možností, jak použít klasickou metodu SVM pro funkcionální data, je projektovat původní data na nějaký $d$-dimenzionální podprostor našeho Hilbertova prostoru $\mathcal H$, označme jej $V_d$.
Předpokládejme, že tento podprostor $V_d$ má ortonormální bázi $\{\Psi_j\}_{j = 1, \dots, d}$.
Definujeme transformaci $P_{V_d}$ jakožto ortogonální projekci na podprostor $V_d$, tedy můžeme psát

$$
P_{V_d} (x) = \sum_{j = 1}^d \langle x, \Psi_j \rangle \Psi_j.
$$

Nyní můžeme pro klasifikaci použít koeficienty z ortogonální projekce, tedy aplikujeme standardní SVM na vektory $\left( \langle x, \Psi_1 \rangle, \dots, \langle x, \Psi_d \rangle\right)^\top$.
Využitím této transformace jsme tedy definovali nové, tzv.
adaptované jádro, které je složené z ortogonální projekce $P_{V_d}$ a jádrové funkce standardní metody podpůrných vektorů.
Máme tedy (adaptované) jádro $Q(x_i, x_j) = K(P_{V_d}(x_i), P_{V_d}(x_j))$.
Jde tedy o metodu redukce dimenze, kterou můžeme nazvat *filtrace*.

Pro samotnou projekci použijeme v `R` funkci `project.basis()` z knihovny `fda`.
Na jejím vstupu bude matice původních diskrétních (nevyhlazených) dat, hodnoty, ve kterých měříme hodnoty v matici původních dat a bázový objekt, na který chceme data projektovat.
My zvolíme projekci na B-splinovou bázi, protože využití Fourierovy báze není pro naše neperiodická data vhodné.

Dimenzi $d$ volíme buď z nějaké předchozí expertní znalosti, nebo pomocí cross-validace.
V našem případě určíme optimální dimenzi podprostoru $V_d$ pomocí $k$-násobné cross-validace (volíme $k \ll n$ kvůli výpočetní náročnosti metody, často se volí $k = 5$ nebo $k = 10$).
Požadujeme B-spliny řádu 4, pro počet bázových funkcí potom platí vztah

$$
n_{basis} = n_{breaks} + n_{order} - 2,
$$

kde $n_{breaks}$ je počet uzlů a $n_{order} = 4$.
Minimální dimenzi tedy (pro $n_{breaks} = 1$) volíme $n_{basis} = 3$ a maximální (pro $n_{breaks} = 51$ odpovídající počtu původních diskrétních dat) $n_{basis} = 53$.
V `R` však hodnota $n_{basis}$ musí být alespoň $n_{order} = 4$ a pro velké hodnoty $n_{basis}$ již dochází k přefitování modelu, tudíž volíme za maximální $n_{basis}$ menší číslo, řekněme 20.


```r
k_cv <- 10 # k-fold CV

# hodnoty pro B-splinovou bazi
rangeval <- range(t)
norder <- 4
n_basis_min <- norder
n_basis_max <- 20

dimensions <- n_basis_min:n_basis_max # vsechny dimenze, ktere chceme vyzkouset

# rozdelime trenovaci data na k casti
folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)

# list se tremi slozkami ... maticemi pro jednotlive jadra -> linear, poly, radial
# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro danou cast trenovaci mnoziny
# v radcich budou hodnoty pro danou hodnotu dimenze
CV.results <- list(SVM.l = matrix(NA, nrow = length(dimensions), ncol = k_cv),
                   SVM.p = matrix(NA, nrow = length(dimensions), ncol = k_cv),
                   SVM.r = matrix(NA, nrow = length(dimensions), ncol = k_cv))

for (d in dimensions) {
  # bazovy objekt
  bbasis <- create.bspline.basis(rangeval = rangeval, 
                                 nbasis = d)
  
  # projekce diskretnich dat na B-splinovou bazi o dimenzi d
  Projection <- project.basis(y = XX, # matice diskretnich dat
                              argvals = t, # vektor argumentu
                              basisobj = bbasis) # bazovy objekt
  
  # rozdeleni na trenovaci a testovaci data v ramci CV
  XX.train <- subset(t(Projection), split == TRUE)
  
  for (index_cv in 1:k_cv) {
    # definice testovaci a trenovaci casti pro CV
    fold <- folds[[index_cv]]
    cv_sample <- 1:dim(XX.train)[1] %in% fold
    
    data.projection.train.cv <- as.data.frame(XX.train[cv_sample, ])
    data.projection.train.cv$Y <- factor(Y.train[cv_sample])
    
    data.projection.test.cv <- as.data.frame(XX.train[!cv_sample, ])
    Y.test.cv <- Y.train[!cv_sample]
    data.projection.test.cv$Y <- factor(Y.test.cv)
  
    # sestrojeni modelu
    clf.SVM.l.projection <- svm(Y ~ ., data = data.projection.train.cv,
                            type = 'C-classification',
                            scale = TRUE,
                            kernel = 'linear')
    
    clf.SVM.p.projection <- svm(Y ~ ., data = data.projection.train.cv,
                            type = 'C-classification',
                            scale = TRUE,
                            coef0 = coef0,
                            kernel = 'polynomial')
    
    clf.SVM.r.projection <- svm(Y ~ ., data = data.projection.train.cv,
                            type = 'C-classification',
                            scale = TRUE,
                            kernel = 'radial')
      
    # presnost na validacnich datech
    ## linear kernel
    predictions.test.l <- predict(clf.SVM.l.projection,
                                  newdata = data.projection.test.cv)
    presnost.test.l <- table(Y.test.cv, predictions.test.l) |>
      prop.table() |> diag() |> sum()
    ## polynomial kernel
    predictions.test.p <- predict(clf.SVM.p.projection, 
                                  newdata = data.projection.test.cv)
    presnost.test.p <- table(Y.test.cv, predictions.test.p) |>
      prop.table() |> diag() |> sum()
    ## radial kernel
    predictions.test.r <- predict(clf.SVM.r.projection,
                                  newdata = data.projection.test.cv)
    presnost.test.r <- table(Y.test.cv, predictions.test.r) |>
      prop.table() |> diag() |> sum()
    
    # presnosti vlozime na pozice pro dane d a fold
    CV.results$SVM.l[d - min(dimensions) + 1, index_cv] <- presnost.test.l
    CV.results$SVM.p[d - min(dimensions) + 1, index_cv] <- presnost.test.p
    CV.results$SVM.r[d - min(dimensions) + 1, index_cv] <- presnost.test.r
  }
}
  
# spocitame prumerne presnosti pro jednotliva d pres folds
for (n_method in 1:length(CV.results)) {
  CV.results[[n_method]] <- apply(CV.results[[n_method]], 1, mean)
}

d.opt <- c(which.max(CV.results$SVM.l) + n_basis_min - 1, 
           which.max(CV.results$SVM.p) + n_basis_min - 1, 
           which.max(CV.results$SVM.r) + n_basis_min - 1)
presnost.opt.cv <- c(max(CV.results$SVM.l), 
                     max(CV.results$SVM.p),
                     max(CV.results$SVM.r))
data.frame(d_opt = d.opt, ERR = 1 - presnost.opt.cv,
           row.names = c('linear', 'poly', 'radial'))
```

```
##        d_opt        ERR
## linear     9 0.01958333
## poly       7 0.03297619
## radial     6 0.14125000
```

Vidíme, že nejlépe vychází hodnota parametru $d$ jako 9 pro lineární jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.9804, 7 pro polynomiální jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.967 a 6 pro radiální jádro s hodnotou přesnosti 0.8588.
Pro přehlednost si ještě vykresleme průběh validačních chybovostí v závislosti na dimenzi $d$.


```r
CV.results <- data.frame(d = dimensions |> rep(3), 
                         CV = c(CV.results$SVM.l, 
                                CV.results$SVM.p, 
                                CV.results$SVM.r),
                         Kernel = rep(c('lineární', 'polynomiální', 'radiální'), 
                                      each = length(dimensions)) |> factor())
CV.results |> ggplot(aes(x = d, y = 1 - CV, colour = Kernel)) + 
  geom_line(linetype = 'dashed') + 
  geom_point(size = 1.5) + 
  geom_point(data = data.frame(d.opt,
                               presnost.opt.cv),
             aes(x = d.opt, y = 1 - presnost.opt.cv), colour = 'black', size = 2) +
  theme_bw() + 
  labs(x = bquote(paste(d)),
       y = 'Validační chybovost',
       colour = 'Jádro') + 
  theme(legend.position = "bottom") + 
  scale_x_continuous(breaks = dimensions)
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-190-1.png" alt="Závislost validační chybovosti na dimenzi podprostoru $V_d$, zvlášť pro všechna tři uvažovaná jádra v metodě SVM. Černými body jsou vyznačeny optimální hodnoty dimenze $V_d$ pro jednotlivé jádrové funkce." width="672" />
<p class="caption">(\#fig:unnamed-chunk-190)Závislost validační chybovosti na dimenzi podprostoru $V_d$, zvlášť pro všechna tři uvažovaná jádra v metodě SVM. Černými body jsou vyznačeny optimální hodnoty dimenze $V_d$ pro jednotlivé jádrové funkce.</p>
</div>

Nyní již můžeme natrénovat jednotlivé klasifikátory na všech trénovacích datech a podívat se na jejich úspěšnost na testovacích datech.
Pro každou jádrovou funkci volíme dimenzi podprostoru, na který projektujeme, podle výsledků cross-validace.

V proměnné `Projection` máme uloženou matici koeficientů ortogonální projekce, tedy

$$
\texttt{Projection} = \begin{pmatrix}
\langle x_1, \Psi_1 \rangle & \langle x_2, \Psi_1 \rangle & \cdots & \langle x_n, \Psi_1 \rangle\\
\langle x_1, \Psi_2 \rangle & \langle x_2, \Psi_2 \rangle & \cdots & \langle x_n, \Psi_2 \rangle\\
\vdots & \vdots & \ddots & \vdots \\
\langle x_1, \Psi_d \rangle & \langle x_2, \Psi_d \rangle & \dots & \langle x_n, \Psi_d \rangle
\end{pmatrix}_{d \times n}.
$$


```r
# pripravime si datovou tabulku pro ulozeni vysledku
Res <- data.frame(model = c('SVM linear - projection', 
                            'SVM poly - projection', 
                            'SVM rbf - projection'), 
                  Err.train = NA,
                  Err.test = NA)

# projdeme jednotliva jadra
for (kernel_number in 1:3) {
  kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
  # bazovy objekt
  bbasis <- create.bspline.basis(rangeval = rangeval, 
                                 nbasis = d.opt[kernel_number])
  
  # projekce diskretnich dat na B-splinovou bazi
  Projection <- project.basis(y = XX, # matice diskretnich dat
                              argvals = t, # vektor argumentu
                              basisobj = bbasis) # bazovy objekt
  
  # rozdeleni na trenovaci a testovaci data
  XX.train <- subset(t(Projection), split == TRUE)
  XX.test <- subset(t(Projection), split == FALSE)
  
  data.projection.train <- as.data.frame(XX.train)
  data.projection.train$Y <- factor(Y.train)
  
  data.projection.test <- as.data.frame(XX.test)
  data.projection.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.SVM.projection <- svm(Y ~ ., data = data.projection.train,
                            type = 'C-classification',
                            scale = TRUE,
                            coef0 = coef0,
                            kernel = kernel_type)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.SVM.projection, newdata = data.projection.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test <- predict(clf.SVM.projection, newdata = data.projection.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  # ulozeni vysledku
  Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
}
```

Přesnost metody SVM aplikované na bázové koeficienty na trénovacích datech je tedy 2 % pro lineární jádro, 2.67 % pro polynomiální jádro a 9.33 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 6.15 % pro lineární jádro, 6.15 % pro polynomiální jádro a 10.77 % pro radiální jádro.


```r
RESULTS <- rbind(RESULTS, Res)
```

#### RKHS + SVM {#RKHS-SVMA3b} 

V této sekci se podíváme na další možnost, jak využít metodu podpůrných vektorů pro klasifikaci funkcionálních dat.
V tomto případě půjde opět o již nám známý princip, kdy nejprve funkcionální data vyjádříme jakožto nějaké konečně-rozměrné objekty a na tyto objekty následně aplikujeme klasickou metodu SVM.

Z poslední části Tvrzení \@ref(thm:MaG) vyplývá, jak máme spočítat v praxi reprezentace křivek.
Budeme pracovat s diskretizovanými daty po vyhlazení křivek.
Nejprve si definujeme jádro pro prostor RKHS.
Využijeme Gaussovské jádro s parametrem $\gamma$.
Hodnota tohoto hyperparametru výrazně ovlivňuje chování a tedy i úspěšnost metody, proto jeho volbě musíme věnovat zvláštní pozornost (volíme pomocí cross-validace).


```r
# hodnoty hypermarametru stejne jako v minule casti
eps <- 0.01
C <- 1 
```

###### Gaussovké jádro


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())

# jadro a jadrova matice ... Gaussovske s parametrem gamma
Gauss.kernel <- function(x, y, gamma) {
  return(exp(-gamma * norm(c(x - y) |> t(), type = 'F')))
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

Spočítejme nyní matici $K_S$ a její vlastní čísla a příslušné vlastní vektory.


```r
# spocitame matici K
gamma <- 0.1 # pevna hodnota gamma, optimalni urcime pomoci CV
K <- Kernel.RKHS(t.seq, gamma = gamma)

# urcime vlastni cisla a vektory
Eig <- eigen(K)
eig.vals <- Eig$values
eig.vectors <- Eig$vectors
```

K výpočtu koeficientů v reprezentaci křivek, tedy výpočtu vektorů $\hat{\boldsymbol \lambda}_l^* = \left( \hat\lambda_{1l}^*, \dots, \hat\lambda_{\hat dl}^*\right)^\top, l = 1, 2, \dots, n$, potřebujeme ještě koeficienty z SVM.
Narozdíl od klasifikačního problému nyní řešíme problém regrese, neboť se snažíme vyjádřit naše pozorované křivky v nějaké (námi zvolené pomocí jádra $K$) bázi.
Proto využijeme metodu *Support Vector Regression*, z níž následně získáme koeficienty $\alpha_{il}$.


```r
# urceni koeficientu alpha z SVM
alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                     ncol = dim(data.RKHS)[2]) # prazdny objekt

# model
for(i in 1:dim(data.RKHS)[2]) {
  df.svm <- data.frame(x = t.seq,
                       y = data.RKHS[, i])
  svm.RKHS <- svm(y ~ x, data = df.svm, 
                  kernel = 'radial',
                  type = 'eps-regression',
                  epsilon = eps,
                  cost = C,
                  gamma = gamma)
  # urceni alpha
  alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
}
```

Nyní již můžeme spočítat reprezentace jednotlivých křivek.
Nejprve zvolme za $\hat d$ celou dimenzi, tedy $\hat d = m ={}$ 101, následně určíme optimální $\hat d$ pomocí cross-validace.


```r
# d
d.RKHS <- dim(alpha.RKHS)[1]

# urceni vektoru lambda
Lambda.RKHS <- matrix(NA, 
                      ncol = dim(data.RKHS)[2], 
                      nrow = d.RKHS) # vytvoreni prazdneho objektu

# vypocet reprezentace
for(l in 1:dim(data.RKHS)[2]) {
  Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
}
```

Nyní máme v matici `Lambda.RKHS` uloženy ve sloupcích vektory $\hat{\boldsymbol \lambda}_l^*, l = 1, 2, \dots, n$ pro každou křivku.
Tyto vektory nyní využijeme jakožto reprezentaci daných křivek a klasifikujeme data podle této diskretizace.


```r
# rozdeleni na trenovaci a testovaci data
XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]

# pripravime si datovou tabulku pro ulozeni vysledku
Res <- data.frame(model = c('SVM linear - RKHS', 
                            'SVM poly - RKHS', 
                            'SVM rbf - RKHS'), 
                  Err.train = NA,
                  Err.test = NA)

# projdeme jednotliva jadra
for (kernel_number in 1:3) {
  kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]

  data.RKHS.train <- as.data.frame(t(XX.train))
  data.RKHS.train$Y <- factor(Y.train)
  
  data.RKHS.test <- as.data.frame(t(XX.test))
  data.RKHS.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                      type = 'C-classification',
                      cost = C,
                      coef0 = coef0,
                      scale = TRUE,
                      kernel = kernel_type)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  # ulozeni vysledku
  Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
}
```


Table: (\#tab:unnamed-chunk-199)Souhrnné výsledky metody SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

Model                                $\widehat{Err}_{train}\quad\quad\quad\quad\quad$         $\widehat{Err}_{test}\quad\quad\quad\quad\quad$       
----------------------------------  -------------------------------------------------------  -------------------------------------------------------
SVM linear - RKHS                                                                         0                                                        0
SVM poly - RKHS                                                                           0                                                        0
SVM rbf - RKHS                                                                            0                                                        0

Vidíme, že model u všech třech jader velmi dobře klasifikuje trénovací data, zatímco jeho úspěšnost na testovacích datech není vůbec dobrá.
Je zřejmé, že došlo k overfittingu, proto využijeme cross-validaci, abychom určili optimální hodnoty $\gamma$ a $d$.


```r
# rozdelime trenovaci data na k casti
folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()

# hodnoty hyperparametru, ktere budeme prochazet
dimensions <- 3:30 # rozumny rozsah hodnot d
gamma.cv <- 10^seq(-2, 2, length = 15)

# list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro dane
# v radcich budou hodnoty pro danou gamma a vrstvy odpovidaji folds
dim.names <- list(gamma = paste0('gamma:', round(gamma.cv, 3)),
                  d = paste0('d:', dimensions),
                  CV = paste0('cv:', 1:k_cv))

CV.results <- list(
  SVM.l = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.p = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.r = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                dimnames = dim.names))
```


```r
# samotna CV
for (gamma in gamma.cv) {
  K <- Kernel.RKHS(t.seq, gamma = gamma)
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'radial',
                    type = 'eps-regression',
                    epsilon = eps,
                    cost = C,
                    gamma = gamma)
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
  }
  
  # projdeme dimenze
  for(d.RKHS in dimensions) {
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) 
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                             alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    # projdeme folds
    for (index_cv in 1:k_cv) {
      # definice testovaci a trenovaci casti pro CV
      fold <- folds[[index_cv]]
      # rozdeleni na trenovaci a validacni data
      XX.train <- Lambda.RKHS[, fold]
      XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
      # pripravime si datovou tabulku pro ulozeni vysledku
      Res <- data.frame(model = c('SVM linear - RKHS', 
                                  'SVM poly - RKHS', 
                                  'SVM rbf - RKHS'), 
                        Err.test = NA)
      # projdeme jednotliva jadra
      for (kernel_number in 1:3) {
        kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
      
        data.RKHS.train <- as.data.frame(t(XX.train))
        data.RKHS.train$Y <- factor(Y.train[fold])
        
        data.RKHS.test <- as.data.frame(t(XX.test))
        data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
        
        # sestrojeni modelu
        clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                            type = 'C-classification',
                            scale = TRUE,
                            cost = C,
                            coef0 = coef0,
                            kernel = kernel_type)
        
        # presnost na validacnich datech
        predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
        presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
          prop.table() |> diag() |> sum()
        
        # ulozeni vysledku
        Res[kernel_number, 2] <- 1 - presnost.test
      }
      # presnosti vlozime na pozice pro dane d, gamma a fold
      CV.results$SVM.l[paste0('gamma:', round(gamma, 3)), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[1, 2]
      CV.results$SVM.p[paste0('gamma:', round(gamma, 3)), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[2, 2]
      CV.results$SVM.r[paste0('gamma:', round(gamma, 3)), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[3, 2]
    }
  }
}
```


```r
# spocitame prumerne presnosti pro jednotliva d pres folds
for (n_method in 1:length(CV.results)) {
  CV.results[[n_method]] <- apply(CV.results[[n_method]], c(1, 2), mean)
}

gamma.opt <- c(which.min(CV.results$SVM.l) %% length(gamma.cv), 
               which.min(CV.results$SVM.p) %% length(gamma.cv), 
               which.min(CV.results$SVM.r) %% length(gamma.cv))
gamma.opt[gamma.opt == 0] <- length(gamma.cv)
gamma.opt <- gamma.cv[gamma.opt]

d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
           which.min(t(CV.results$SVM.p)) %% length(dimensions), 
           which.min(t(CV.results$SVM.r)) %% length(dimensions))
d.opt[d.opt == 0] <- length(dimensions)
d.opt <- dimensions[d.opt]

err.opt.cv <- c(min(CV.results$SVM.l), 
                     min(CV.results$SVM.p),
                     min(CV.results$SVM.r))
df.RKHS.res <- data.frame(d = d.opt, gamma = gamma.opt, CV = err.opt.cv,
           Kernel = c('linear', 'polynomial', 'radial') |> factor(),
           row.names = c('linear', 'poly', 'radial'))
```


Table: (\#tab:unnamed-chunk-202)Souhrnné výsledky cross-validace pro metodu SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

          $\quad\quad\quad\quad\quad d$   $\quad\quad\quad\quad\quad\gamma$   $\widehat{Err}_{cross\_validace}$  Model                             
-------  ------------------------------  ----------------------------------  ----------------------------------  ----------------------------------
linear                               12                              1.0000                                   0  linear                            
poly                                 18                              1.0000                                   0  polynomial                        
radial                               17                              1.9307                                   0  radial                            

Vidíme, že nejlépe vychází hodnota parametru $d={}$ 12 a $\gamma={}$ 1 pro lineární jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 1, $d={}$ 18 a $\gamma={}$ 1 pro polynomiální jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 1 a $d={}$ 17 a $\gamma={}$ 1.9307 pro radiální jádro s hodnotou přesnosti 1.
Pro zajímavost si ještě vykresleme funkci validační chybovosti v závislosti na dimenzi $d$ a hodnotě hyperparametru $\gamma$.


```r
CV.results.plot <- data.frame(d = rep(dimensions |> rep(3), each = length(gamma.cv)), 
                              gamma = rep(gamma.cv, length(dimensions)) |> rep(3),
                               CV = c(c(CV.results$SVM.l), 
                                      c(CV.results$SVM.p), 
                                      c(CV.results$SVM.r)),
                               Kernel = rep(c('linear', 'polynomial', 'radial'), 
                                            each = length(dimensions) * 
                                              length(gamma.cv)) |> factor())
CV.results.plot |> 
  ggplot(aes(x = d, y = gamma, z = CV)) + 
  geom_contour_filled() +
  scale_y_continuous(trans='log10') +
  facet_wrap(~Kernel) +
  theme_bw() + 
  labs(x = expression(d),
       y = expression(gamma)) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_point(data = df.RKHS.res, aes(x = d, y = gamma),
             size = 5, pch = '+')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-203-1.png" alt="Závislost validační chybovosti na volbě hyperparametrů $d$ a $\gamma$, zvlášť pro všechna tři uvažovaná jádra v metodě SVM." width="672" />
<p class="caption">(\#fig:unnamed-chunk-203)Závislost validační chybovosti na volbě hyperparametrů $d$ a $\gamma$, zvlášť pro všechna tři uvažovaná jádra v metodě SVM.</p>
</div>

Na grafech výše vidíme, jak se měnila validační chybovost v závislosti na hodnotách hyperparametrů $d$ a $\gamma$.
Všimněme si zejména, že ve všech třech grafech pro jednotlivá jádra jsou patrné výrazné horizontální útvary.
Z toho můžeme usoudit významné teoretické i praktické zjištění -- uvažovaná klasifikační metoda (projekce na RKHS pomocí SVM + klasifikace SVM) je robustní na volbu hyperparametru $d$ (tj. při malé změně v hodnotě tohoto parametru nedojde k výraznému zhoršení validační chybovosti), zatímco při volbě hyperparametru $\gamma$ musíme být velmi obezřetní (i malá změna v jeho hodnotě může vést k velké změně validační chybovosti).
Toto chování je nejlépe patrné u Gaussova jádra.

Jelikož již máme nalezeny optimální hodnoty hyperparametrů, můžeme zkounstruovat finální modely a určit jejich úspěšnost klasifikace na testovacích datech.


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
```


```r
# pripravime si datovou tabulku pro ulozeni vysledku
Res <- data.frame(model = c('SVM linear - RKHS - radial', 
                            'SVM poly - RKHS - radial', 
                            'SVM rbf - RKHS - radial'), 
                  Err.train = NA,
                  Err.test = NA)

# projdeme jednotliva jadra
for (kernel_number in 1:3) {
  # spocitame matici K
  gamma <- gamma.opt[kernel_number] # hodnota gamma pomoci CV
  K <- Kernel.RKHS(t.seq, gamma = gamma)
  
  # urcime vlastni cisla a vektory
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  # urceni koeficientu alpha z SVM
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                       ncol = dim(data.RKHS)[2]) # prazdny objekt
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'radial',
                    type = 'eps-regression',
                    epsilon = eps,
                    cost = C,
                    gamma = gamma)
    # urceni alpha
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
  }
  # d
  d.RKHS <- d.opt[kernel_number]
  
  # urceni vektoru lambda
  Lambda.RKHS <- matrix(NA, 
                        ncol = dim(data.RKHS)[2], 
                        nrow = d.RKHS) # vytvoreni prazdneho objektu
  
  # vypocet reprezentace
  for(l in 1:dim(data.RKHS)[2]) {
    Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
  }
  
  # rozdeleni na trenovaci a testovaci data
  XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
  XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]

  kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]

  data.RKHS.train <- as.data.frame(t(XX.train))
  data.RKHS.train$Y <- factor(Y.train)
  
  data.RKHS.test <- as.data.frame(t(XX.test))
  data.RKHS.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                      type = 'C-classification',
                      scale = TRUE,
                      cost = C,
                      coef0 = coef0,
                      kernel = kernel_type)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  # ulozeni vysledku
  Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
}
```


Table: (\#tab:unnamed-chunk-206)Souhrnné výsledky metody SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

Model                                $\widehat{Err}_{train}\quad\quad\quad\quad\quad$         $\widehat{Err}_{test}\quad\quad\quad\quad\quad$       
----------------------------------  -------------------------------------------------------  -------------------------------------------------------
SVM linear - RKHS - radial                                                                0                                                   0.0308
SVM poly - RKHS - radial                                                                  0                                                   0.0769
SVM rbf - RKHS - radial                                                                   0                                                   0.0308

Přesnost metody SVM v kombinaci s projekcí na Reproducing Kernel Hilbert Space je tedy na trénovacích datech rovna 0 % pro lineární jádro, 0 % pro polynomiální jádro a 0 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 3.08 % pro lineární jádro, 7.69 % pro polynomiální jádro a 3.08 % pro radiální jádro.


```r
RESULTS <- rbind(RESULTS, Res)
```

###### Polynomiální jádro


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())

# jadro a jadrova matice ... polynomialni s parametrem p
Poly.kernel <- function(x, y, p) {
  return((1 + x * y)^p)
}

Kernel.RKHS <- function(x, p) {
  K <- matrix(NA, ncol = length(x), nrow = length(x))
  for(i in 1:nrow(K)) {
    for(j in 1:ncol(K)) {
      K[i, j] <- Poly.kernel(x = x[i], y = x[j], p)
    }
  }
  return(K)
}
```


```r
# rozdelime trenovaci data na k casti
folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()

# hodnoty hyperparametru, ktere budeme prochazet
dimensions <- 2:30 # rozumny rozsah hodnot d
poly.cv <- 2:5

# list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro dane
# v radcich budou hodnoty pro dane p a vrstvy odpovidaji folds
dim.names <- list(p = paste0('p:', poly.cv),
                  d = paste0('d:', dimensions),
                  CV = paste0('cv:', 1:k_cv))

CV.results <- list(
  SVM.l = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.p = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.r = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                dimnames = dim.names))
```


```r
# samotna CV
for (p in poly.cv) {
  K <- Kernel.RKHS(t.seq, p = p)
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'polynomial',
                    type = 'eps-regression',
                    coef0 = 1,
                    cost = C,
                    epsilon = eps,
                    degree = p)
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
  }
  
  # projdeme dimenze
  for(d.RKHS in dimensions) {
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) 
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                             alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    # projdeme folds
    for (index_cv in 1:k_cv) {
      # definice testovaci a trenovaci casti pro CV
      fold <- folds[[index_cv]]
      # rozdeleni na trenovaci a validacni data
      XX.train <- Lambda.RKHS[, fold]
      XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
      # pripravime si datovou tabulku pro ulozeni vysledku
      Res <- data.frame(model = c('SVM linear - RKHS', 
                                  'SVM poly - RKHS', 
                                  'SVM rbf - RKHS'), 
                        Err.test = NA)
      # projdeme jednotliva jadra
      for (kernel_number in 1:3) {
        kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
      
        data.RKHS.train <- as.data.frame(t(XX.train))
        data.RKHS.train$Y <- factor(Y.train[fold])
        
        data.RKHS.test <- as.data.frame(t(XX.test))
        data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
        
        # sestrojeni modelu
        clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                            type = 'C-classification',
                            scale = TRUE,
                            cost = C,                    
                            coef0 = 1,
                            gamma = 1,
                            kernel = kernel_type)
        
        # presnost na validacnich datech
        predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
        presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
          prop.table() |> diag() |> sum()
        
        # ulozeni vysledku
        Res[kernel_number, 2] <- 1 - presnost.test
      }
      # presnosti vlozime na pozice pro dane d, gamma a fold
      CV.results$SVM.l[paste0('p:', p), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[1, 2]
      CV.results$SVM.p[paste0('p:', p), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[2, 2]
      CV.results$SVM.r[paste0('p:', p), 
                       d.RKHS - min(dimensions) + 1, 
                       index_cv] <- Res[3, 2]
    }
  }
}
```


```r
# spocitame prumerne presnosti pro jednotliva d pres folds
for (n_method in 1:length(CV.results)) {
  CV.results[[n_method]] <- apply(CV.results[[n_method]], c(1, 2), mean)
}

poly.opt <- c(which.min(CV.results$SVM.l) %% length(poly.cv), 
               which.min(CV.results$SVM.p) %% length(poly.cv), 
               which.min(CV.results$SVM.r) %% length(poly.cv))
poly.opt[poly.opt == 0] <- length(poly.cv)
poly.opt <- poly.cv[poly.opt]

d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
           which.min(t(CV.results$SVM.p)) %% length(dimensions), 
           which.min(t(CV.results$SVM.r)) %% length(dimensions))
d.opt[d.opt == 0] <- length(dimensions)
d.opt <- dimensions[d.opt]

err.opt.cv <- c(min(CV.results$SVM.l), 
                     min(CV.results$SVM.p),
                     min(CV.results$SVM.r))
df.RKHS.res <- data.frame(d = d.opt, p = poly.opt, CV = err.opt.cv,
           Kernel = c('linear', 'polynomial', 'radial') |> factor(),
           row.names = c('linear', 'poly', 'radial'))
```


Table: (\#tab:unnamed-chunk-211)Souhrnné výsledky cross-validace pro metodu SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

          $\quad\quad\quad\quad\quad d$   $\quad\quad\quad\quad\quad p$   $\widehat{Err}_{cross\_validace}$  Model                             
-------  ------------------------------  ------------------------------  ----------------------------------  ----------------------------------
linear                               30                               3                              0.0077  linear                            
poly                                 17                               2                              0.0077  polynomial                        
radial                                3                               4                              0.0254  radial                            

Vidíme, že nejlépe vychází hodnota parametru $d={}$ 30 a $p={}$ 3 pro lineární jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.9923, $d={}$ 17 a $p={}$ 2 pro polynomiální jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.9923 a $d={}$ 3 a $p={}$ 4 pro radiální jádro s hodnotou přesnosti 0.9746.

Jelikož již máme nalezeny optimální hodnoty hyperparametrů, můžeme zkounstruovat finální modely a určit jejich úspěšnost klasifikace na testovacích datech.


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
```


```r
# pripravime si datovou tabulku pro ulozeni vysledku
Res <- data.frame(model = c('SVM linear - RKHS - poly', 
                            'SVM poly - RKHS - poly', 
                            'SVM rbf - RKHS - poly'), 
                  Err.train = NA,
                  Err.test = NA)

# projdeme jednotliva jadra
for (kernel_number in 1:3) {
  # spocitame matici K
  p <- poly.opt[kernel_number] # hodnota gamma pomoci CV
  K <- Kernel.RKHS(t.seq, p = p)
  
  # urcime vlastni cisla a vektory
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  # urceni koeficientu alpha z SVM
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                       ncol = dim(data.RKHS)[2]) # prazdny objekt
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'polynomial',
                    type = 'eps-regression',
                    epsilon = eps,
                    coef0 = 1,
                    cost = C,
                    gamma = 1,
                    degree = p)
    # urceni alpha
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
  }
  # d
  d.RKHS <- d.opt[kernel_number]
  
  # urceni vektoru lambda
  Lambda.RKHS <- matrix(NA, 
                        ncol = dim(data.RKHS)[2], 
                        nrow = d.RKHS) # vytvoreni prazdneho objektu
  
  # vypocet reprezentace
  for(l in 1:dim(data.RKHS)[2]) {
    Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
  }
  
  # rozdeleni na trenovaci a testovaci data
  XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
  XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]

  kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]

  data.RKHS.train <- as.data.frame(t(XX.train))
  data.RKHS.train$Y <- factor(Y.train)
  
  data.RKHS.test <- as.data.frame(t(XX.test))
  data.RKHS.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                      type = 'C-classification',
                      scale = TRUE,
                      coef0 = 1,
                      cost = C,
                      gamma = 1,
                      kernel = kernel_type)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  # ulozeni vysledku
  Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
}
```


Table: (\#tab:unnamed-chunk-214)Souhrnné výsledky metody SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

Model                                $\widehat{Err}_{train}\quad\quad\quad\quad\quad$         $\widehat{Err}_{test}\quad\quad\quad\quad\quad$       
----------------------------------  -------------------------------------------------------  -------------------------------------------------------
SVM linear - RKHS - poly                                                               0.00                                                   0.0308
SVM poly - RKHS - poly                                                                 0.00                                                   0.0154
SVM rbf - RKHS - poly                                                                  0.02                                                   0.0308

Přesnost metody SVM v kombinaci s projekcí na Reproducing Kernel Hilbert Space je tedy na trénovacích datech rovna 0 % pro lineární jádro, 0 % pro polynomiální jádro a 2 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 3.08 % pro lineární jádro, 1.54 % pro polynomiální jádro a 3.08 % pro radiální jádro.


```r
RESULTS <- rbind(RESULTS, Res)
```

###### Lineární jádro


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())

# jadro a jadrova matice ... polynomialni s parametrem p
Linear.kernel <- function(x, y) {
  return(x * y)
}

Kernel.RKHS <- function(x) {
  K <- matrix(NA, ncol = length(x), nrow = length(x))
  for(i in 1:nrow(K)) {
    for(j in 1:ncol(K)) {
      K[i, j] <- Linear.kernel(x = x[i], y = x[j])
    }
  }
  return(K)
}
```


```r
# rozdelime trenovaci data na k casti
folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()

# hodnoty hyperparametru, ktere budeme prochazet
dimensions <- 2:40 # rozumny rozsah hodnot d

# list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
# prazdna matice, do ktere vlozime jednotlive vysledky
# ve sloupcich budou hodnoty presnosti pro dane d
# v radcich budou hodnoty pro vrstvy odpovidaji folds
dim.names <- list(d = paste0('d:', dimensions),
                  CV = paste0('cv:', 1:k_cv))

CV.results <- list(
  SVM.l = array(NA, dim = c(length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.p = array(NA, dim = c(length(dimensions), k_cv),
                dimnames = dim.names),
  SVM.r = array(NA, dim = c(length(dimensions), k_cv),
                dimnames = dim.names))
```


```r
# samotna CV
K <- Kernel.RKHS(t.seq)
Eig <- eigen(K)
eig.vals <- Eig$values
eig.vectors <- Eig$vectors
alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 

# model
for(i in 1:dim(data.RKHS)[2]) {
  df.svm <- data.frame(x = t.seq,
                       y = data.RKHS[, i])
  svm.RKHS <- svm(y ~ x, data = df.svm, 
                  kernel = 'linear',
                  type = 'eps-regression',
                  epsilon = eps,                   
                  coef0 = 1,
                  gamma = 1,
                  cost = C)
  alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
}

# projdeme dimenze
for(d.RKHS in dimensions) {
  Lambda.RKHS <- matrix(NA, 
                        ncol = dim(data.RKHS)[2], 
                        nrow = d.RKHS) 
  # vypocet reprezentace
  for(l in 1:dim(data.RKHS)[2]) {
    Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                           alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
  }
  # projdeme folds
  for (index_cv in 1:k_cv) {
    # definice testovaci a trenovaci casti pro CV
    fold <- folds[[index_cv]]
    # rozdeleni na trenovaci a validacni data
    XX.train <- Lambda.RKHS[, fold]
    XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
    # pripravime si datovou tabulku pro ulozeni vysledku
    Res <- data.frame(model = c('SVM linear - RKHS', 
                                'SVM poly - RKHS', 
                                'SVM rbf - RKHS'), 
                      Err.test = NA)
    # projdeme jednotliva jadra
    for (kernel_number in 1:3) {
      kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
    
      data.RKHS.train <- as.data.frame(t(XX.train))
      data.RKHS.train$Y <- factor(Y.train[fold])
      
      data.RKHS.test <- as.data.frame(t(XX.test))
      data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
      
      # sestrojeni modelu
      clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                          type = 'C-classification',
                          scale = TRUE,
                          kernel = kernel_type,
                          epsilon = eps,                   
                          coef0 = 1,
                          gamma = 1,
                          cost = C)
      
      # presnost na validacnich datech
      predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
      presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
        prop.table() |> diag() |> sum()
      
      # ulozeni vysledku
      Res[kernel_number, 2] <- 1 - presnost.test
    }
    # presnosti vlozime na pozice pro dane d, gamma a fold
    CV.results$SVM.l[d.RKHS - min(dimensions) + 1, 
                     index_cv] <- Res[1, 2]
    CV.results$SVM.p[d.RKHS - min(dimensions) + 1, 
                     index_cv] <- Res[2, 2]
    CV.results$SVM.r[d.RKHS - min(dimensions) + 1, 
                     index_cv] <- Res[3, 2]
  }
}
```


```r
# spocitame prumerne presnosti pro jednotliva d pres folds
for (n_method in 1:length(CV.results)) {
  CV.results[[n_method]] <- apply(CV.results[[n_method]], 1, mean)
}

d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
           which.min(t(CV.results$SVM.p)) %% length(dimensions), 
           which.min(t(CV.results$SVM.r)) %% length(dimensions))
d.opt[d.opt == 0] <- length(dimensions)
d.opt <- dimensions[d.opt]

err.opt.cv <- c(min(CV.results$SVM.l), 
                     min(CV.results$SVM.p),
                     min(CV.results$SVM.r))
df.RKHS.res <- data.frame(d = d.opt, CV = err.opt.cv,
           Kernel = c('linear', 'polynomial', 'radial') |> factor(),
           row.names = c('linear', 'poly', 'radial'))
```


Table: (\#tab:unnamed-chunk-219)Souhrnné výsledky cross-validace pro metodu SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

          $\quad\quad\quad\quad\quad d$   $\widehat{Err}_{cross\_validace}$  Model                             
-------  ------------------------------  ----------------------------------  ----------------------------------
linear                               31                              0.0267  linear                            
poly                                 21                              0.0267  polynomial                        
radial                                6                              0.0592  radial                            

Vidíme, že nejlépe vychází hodnota parametru $d={}$ 31 pro lineární jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.9733, $d={}$ 21 pro polynomiální jádro s hodnotou přesnosti spočtenou pomocí 10-násobné CV 0.9733 a $d={}$ 6 pro radiální jádro s hodnotou přesnosti 0.9408.

Jelikož již máme nalezeny optimální hodnoty hyperparametrů, můžeme zkounstruovat finální modely a určit jejich úspěšnost klasifikace na testovacích datech.


```r
# odstranime posledni sloupec, ve kterem jsou hodnoty Y
data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
# pridame i testovaci data
data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
```


```r
# pripravime si datovou tabulku pro ulozeni vysledku
Res <- data.frame(model = c('SVM linear - RKHS - linear', 
                            'SVM poly - RKHS - linear', 
                            'SVM rbf - RKHS - linear'), 
                  Err.train = NA,
                  Err.test = NA)

# projdeme jednotliva jadra
for (kernel_number in 1:3) {
  # spocitame matici K
  K <- Kernel.RKHS(t.seq)
  
  # urcime vlastni cisla a vektory
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  # urceni koeficientu alpha z SVM
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                       ncol = dim(data.RKHS)[2]) # prazdny objekt
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'linear',
                    type = 'eps-regression',
                    epsilon = eps,                   
                    coef0 = 1,
                    gamma = 1,
                    cost = C)
    # urceni alpha
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
  }
  # d
  d.RKHS <- d.opt[kernel_number]
  
  # urceni vektoru lambda
  Lambda.RKHS <- matrix(NA, 
                        ncol = dim(data.RKHS)[2], 
                        nrow = d.RKHS) # vytvoreni prazdneho objektu
  
  # vypocet reprezentace
  for(l in 1:dim(data.RKHS)[2]) {
    Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
  }
  
  # rozdeleni na trenovaci a testovaci data
  XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
  XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]

  kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]

  data.RKHS.train <- as.data.frame(t(XX.train))
  data.RKHS.train$Y <- factor(Y.train)
  
  data.RKHS.test <- as.data.frame(t(XX.test))
  data.RKHS.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                      type = 'C-classification',
                      scale = TRUE,
                      kernel = kernel_type,
                      epsilon = eps,                   
                      coef0 = 1,
                      gamma = 1,
                      cost = C)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  # ulozeni vysledku
  Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
}
```


Table: (\#tab:unnamed-chunk-222)Souhrnné výsledky metody SVM v kombinaci s RKHS na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

Model                                $\widehat{Err}_{train}\quad\quad\quad\quad\quad$         $\widehat{Err}_{test}\quad\quad\quad\quad\quad$       
----------------------------------  -------------------------------------------------------  -------------------------------------------------------
SVM linear - RKHS - linear                                                             0.00                                                   0.0462
SVM poly - RKHS - linear                                                               0.00                                                   0.0308
SVM rbf - RKHS - linear                                                                0.02                                                   0.1231

Přesnost metody SVM v kombinaci s projekcí na Reproducing Kernel Hilbert Space je tedy na trénovacích datech rovna 0 % pro lineární jádro, 0 % pro polynomiální jádro a 2 % pro gaussovské jádro.
Na testovacích datech je potom přesnost metody 4.62 % pro lineární jádro, 3.08 % pro polynomiální jádro a 12.31 % pro radiální jádro.


```r
RESULTS <- rbind(RESULTS, Res)
```

## Tabulka výsledků

Z tabulky níže si všimněme zejména dvou podstatných věcí. Tou první je, že metody klasifikují data podstatně lépe než v situaci původních nederivovaných dat. U některých metod je zlepšení i v řádech desítek procent. Druhou podstatnou věcí je fakt, že nyní není takový výrazný rozdíl mezi výsledky jednotlivých metod.


Table: (\#tab:unnamed-chunk-224)Souhrnné výsledky použitých metod na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti a $\widehat{Err}_{test}$ testovací chybovosti.

Model                                $\widehat{Err}_{train}\quad\quad\quad\quad\quad$         $\widehat{Err}_{test}\quad\quad\quad\quad\quad$       
----------------------------------  -------------------------------------------------------  -------------------------------------------------------
KNN                                                                                  0.0133                                                   0.0769
LDA                                                                                  0.0400                                                   0.0923
QDA                                                                                  0.0067                                                   0.0154
LR functional                                                                        0.0000                                                   0.0769
LR score                                                                             0.0067                                                   0.0462
Tree - diskr.                                                                        0.0067                                                   0.0154
Tree - score                                                                         0.0067                                                   0.0615
Tree - Bbasis                                                                        0.0067                                                   0.0154
RForest - diskr                                                                      0.0000                                                   0.0462
RForest - score                                                                      0.0067                                                   0.0462
RForest - Bbasis                                                                     0.0000                                                   0.0308
SVM linear - diskr                                                                   0.0067                                                   0.0154
SVM poly - diskr                                                                     0.0067                                                   0.0615
SVM rbf - diskr                                                                      0.0067                                                   0.0154
SVM linear - PCA                                                                     0.0067                                                   0.0462
SVM poly - PCA                                                                       0.0067                                                   0.0154
SVM rbf - PCA                                                                        0.0067                                                   0.0308
SVM linear - Bbasis                                                                  0.0067                                                   0.0769
SVM poly - Bbasis                                                                    0.0067                                                   0.0769
SVM rbf - Bbasis                                                                     0.0067                                                   0.0923
SVM linear - projection                                                              0.0200                                                   0.0615
SVM poly - projection                                                                0.0267                                                   0.0615
SVM rbf - projection                                                                 0.0933                                                   0.1077
SVM linear - RKHS - radial                                                           0.0000                                                   0.0308
SVM poly - RKHS - radial                                                             0.0000                                                   0.0769
SVM rbf - RKHS - radial                                                              0.0000                                                   0.0308
SVM linear - RKHS - poly                                                             0.0000                                                   0.0308
SVM poly - RKHS - poly                                                               0.0000                                                   0.0154
SVM rbf - RKHS - poly                                                                0.0200                                                   0.0308
SVM linear - RKHS - linear                                                           0.0000                                                   0.0462
SVM poly - RKHS - linear                                                             0.0000                                                   0.0308
SVM rbf - RKHS - linear                                                              0.0200                                                   0.1231

## Simulační studie {#simstudyA3}

V celé předchozí části jsme se zabývali pouze jedním náhodně vygenerovaným souborem funkcí ze dvou klasifikačních tříd, který jsme následně opět náhodně rozdělili na testovací a trénovací část.
Poté jsme jednotlivé klasifikátory získané pomocí uvažovaných metod ohodnotili na základě testovací a trénovací chybovosti.

Jelikož se vygenerovaná data (a jejich rozdělení na dvě části) mohou při každém zopakování výrazně lišit, budou se i chybovosti jednotlivých klasifikačních algoritmů výrazně lišit.
Proto dělat jakékoli závěry o metodách a porovnávat je mezi sebou může být na základě jednoho vygenerovaného datového souboru velmi zavádějící.

Z tohoto důvodu se v této části zaměříme na opakování celého předchozího postupu pro různé vygenerované soubory.
Výsledky si budeme ukládat do tabulky a nakonec spočítáme průměrné charakteristiky modelů přes jednotlivá opakování.
Aby byly naše závěry dostatečně obecné, zvolíme počet opakování $n_{sim} = 100$.

### Simulace pro nederivovaná data

Nejprve se podívejme na simulaci původních, tedy nederivovaných, dat. Vzhledem k vysokému počtu uvažovaných klasifikačních metod proveďme simulace zvlášť se stejným nastavením generátoru pseudonáhodných čísel (`set.seed(42)`), poté můžeme výsledky simulací porovnat  mezi sebou.


```r
# nastaveni generatoru pseudonahodnych cisel
set.seed(42)

# pocet simulaci
n.sim <- 100

## list, do ktereho budeme ukladat hodnoty chybovosti
# ve sloupcich budou metody
# v radcich budou jednotliva opakovani
# list ma dve polozky ... train a test
methods <- c('KNN', 'LDA', 'QDA', 'LR_functional', 'LR_score', 'Tree_discr',
             'Tree_score', 'Tree_Bbasis', 'RF_discr', 'RF_score', 'RF_Bbasis', 
             'SVM linear - diskr', 'SVM poly - diskr', 'SVM rbf - diskr', 
             'SVM linear - PCA', 'SVM poly - PCA', 'SVM rbf - PCA', 
             'SVM linear - Bbasis', 'SVM poly - Bbasis', 'SVM rbf - Bbasis',
             'SVM linear - projection', 'SVM poly - projection', 
             'SVM rbf - projection', 'SVM linear - RKHS - radial', 
             'SVM poly - RKHS - radial', 'SVM rbf - RKHS - radial', 
             'SVM linear - RKHS - poly', 'SVM poly - RKHS - poly', 
             'SVM rbf - RKHS - poly', 'SVM linear - RKHS - linear', 
             'SVM poly - RKHS - linear', 'SVM rbf - RKHS - linear')

SIMULACE <- list(train = as.data.frame(matrix(NA, ncol = length(methods), 
                                             nrow = n.sim,
                                             dimnames = list(1:n.sim, methods))), 
                 test = as.data.frame(matrix(NA, ncol = length(methods), 
                                             nrow = n.sim,
                                             dimnames = list(1:n.sim, methods))))

# objekt na ulozeni optimalnich hodnot hyperparametru, ktere se urcuji pomoci CV
CV_RESULTS <- data.frame(KNN_K = rep(NA, n.sim), 
                         nharm = NA, 
                         LR_func_n_basis = NA,
                         SVM_d_Linear = NA,
                         SVM_d_Poly = NA,
                         SVM_d_Radial = NA, 
                         SVM_RKHS_radial_gamma1 = NA,
                         SVM_RKHS_radial_gamma2 = NA,
                         SVM_RKHS_radial_gamma3 = NA,
                         SVM_RKHS_radial_d1 = NA,
                         SVM_RKHS_radial_d2 = NA,
                         SVM_RKHS_radial_d3 = NA,
                         SVM_RKHS_poly_p1 = NA,
                         SVM_RKHS_poly_p2 = NA,
                         SVM_RKHS_poly_p3 = NA,
                         SVM_RKHS_poly_d1 = NA,
                         SVM_RKHS_poly_d2 = NA,
                         SVM_RKHS_poly_d3 = NA,
                         SVM_RKHS_linear_d1 = NA,
                         SVM_RKHS_linear_d2 = NA,
                         SVM_RKHS_linear_d3 = NA)
```

Nyní zopakujeme celou předchozí část 100-krát a hodnoty chybovostí si budeme ukládat to listu `SIMULACE`.
Do datové tabulky `CV_RESULTS` si potom budeme ukládat hodnoty optimálních hyperparametrů -- pro metodu $K$ nejbližších sousedů a pro SVM hodnotu dimenze $d$ v případě projekce na B-splinovou bázi. Také uložíme všechny hodnoty hyperparametrů pro metodu SVM + RKHS.


```r
# nastaveni generatoru pseudonahodnych cisel
set.seed(42)

## SIMULACE

for(sim in 1:n.sim) {
  # rozdeleni na testovaci a trenovaci cast
  split <- sample.split(XXfd$fdnames$reps, SplitRatio = 0.7)
  
  # vytvoreni vektoru 0 a 1, 0 pro < 20 a 1 pro > 20 
  Y <- ifelse(labels == 'large', 1, 0)
  
  X.train <- subset(XXfd, split == TRUE)
  X.test <- subset(XXfd, split == FALSE)
  
  Y.train <- subset(Y, split == TRUE)
  Y.test <- subset(Y, split == FALSE)
  
  x.train <- fdata(X.train)
  y.train <- as.numeric(factor(Y.train))
  
  ## 1) K nejbližších sousedů
  
  k_cv <- 10 # k-fold CV
  neighbours <- c(1:10) # pocet sousedu 
  
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)
  
  # kontrola, ze mame opravdu k = k_cv
  while (length(folds) != k_cv) {
    folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)
  }
  
  CV.results <- matrix(NA, nrow = length(neighbours), ncol = k_cv)
  
  for (index in 1:k_cv) {
    # definujeme danou indexovou mnozinu
    fold <- folds[[index]]
      
    x.train.cv <- subset(X.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
      fdata()
    y.train.cv <- subset(Y.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
      factor() |> as.numeric()
    
    x.test.cv <- subset(X.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
      fdata()
    y.test.cv <- subset(Y.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
      factor() |> as.numeric()
    
    # projdeme kazdou cast ... k-krat zopakujeme
    for(neighbour in neighbours) {
      # model pro konkretni volbu K
      neighb.model <- classif.knn(group = y.train.cv, 
                                fdataobj = x.train.cv, 
                                knn = neighbour) 
      # predikce na validacni casti
      model.neighb.predict <- predict(neighb.model, 
                                      new.fdataobj = x.test.cv)
      # presnost na validacni casti
      presnost <- table(y.test.cv, model.neighb.predict) |> 
        prop.table() |> diag() |> sum()
      
      # presnost vlozime na pozici pro dane K a fold
      CV.results[neighbour, index] <- presnost
    }
  }
  
  # spocitame prumerne presnosti pro jednotliva K pres folds
  CV.results <- apply(CV.results, 1, mean)
  K.opt <- which.max(CV.results)
  CV_RESULTS$KNN_K[sim] <- K.opt
  presnost.opt.cv <- max(CV.results)
  CV.results <- data.frame(K = neighbours, CV = CV.results)
  
  neighb.model <- classif.knn(group = y.train, fdataobj = x.train, knn = K.opt)
  
  # predikce
  model.neighb.predict <- predict(neighb.model, 
                                  new.fdataobj = fdata(X.test))
  
  presnost <- table(as.numeric(factor(Y.test)), model.neighb.predict) |>
    prop.table() |>
    diag() |>
    sum()
  
  RESULTS <- data.frame(model = 'KNN', 
                        Err.train = 1 - neighb.model$max.prob,
                        Err.test = 1 - presnost)
  
  ## 2) Lineární diskriminační analýza
  
  # analyza hlavnich komponent
  data.PCA <- pca.fd(X.train, nharm = 10) # nharm - maximalni pocet HK
  nharm <- which(cumsum(data.PCA$varprop) >= 0.9)[1] # urceni p
  CV_RESULTS$nharm[sim] <- nharm
  if(nharm == 1) nharm <- 2
  
  data.PCA <- pca.fd(X.train, nharm = nharm) 
  data.PCA.train <- as.data.frame(data.PCA$scores) # skore prvnich p HK
  data.PCA.train$Y <- factor(Y.train) # prislusnost do trid
  
  # vypocet skoru testovacich funkci
  scores <- matrix(NA, ncol = nharm, nrow = length(Y.test)) # prazdna matice 
  
  for(k in 1:dim(scores)[1]) {
    xfd = X.test[k] - data.PCA$meanfd[1] # k-te pozorovani - prumerna funkce
    scores[k, ] = inprod(xfd, data.PCA$harmonics) 
    # skalarni soucin rezidua a vlastnich funkci rho (funkcionalni hlavni komponenty)
  }
  
  data.PCA.test <- as.data.frame(scores)
  data.PCA.test$Y <- factor(Y.test)
  colnames(data.PCA.test) <- colnames(data.PCA.train) 
  
  # model
  clf.LDA <- lda(Y ~ ., data = data.PCA.train)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.LDA, newdata = data.PCA.train)
  presnost.train <- table(data.PCA.train$Y, predictions.train$class) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.LDA, newdata = data.PCA.test)
  presnost.test <- table(data.PCA.test$Y, predictions.test$class) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'LDA', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 3) Kvadratická diskriminační analýza
  
  # model
  clf.QDA <- qda(Y ~ ., data = data.PCA.train)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.QDA, newdata = data.PCA.train)
  presnost.train <- table(data.PCA.train$Y, predictions.train$class) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.QDA, newdata = data.PCA.test)
  presnost.test <- table(data.PCA.test$Y, predictions.test$class) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'QDA', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 4) Logistická regrese
  ### 4.1) Funkcionální logistická regrese
  
  # vytvorime vhodne objekty
  x.train <- fdata(X.train)
  y.train <- as.numeric(Y.train)
  
  # body, ve kterych jsou funkce vyhodnoceny
  tt <- x.train[["argvals"]]
  
  dataf <- as.data.frame(y.train) 
  colnames(dataf) <- "Y"
  # B-spline baze 
  nbasis.x <- 7
  basis1 <- create.bspline.basis(rangeval = range(tt), 
                                 nbasis = nbasis.x)
  
  ### 10-fold cross-validation
  n.basis.max <- 25
  n.basis <- 4:n.basis.max
  k_cv <- 10 # k-fold CV
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)
  # kontrola, ze mame opravdu k = k_cv
  while (length(folds) != k_cv) {
    folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)
  }
  ## prvky, ktere se behem cyklu nemeni
  # body, ve kterych jsou funkce vyhodnoceny
  tt <- x.train[["argvals"]]
  rangeval <- range(tt)
  # vztah
  f <- Y ~ x
  # baze pro x
  basis.x <- list("x" = basis1)
  
  CV.results <- matrix(NA, nrow = length(n.basis), ncol = k_cv, 
                       dimnames = list(n.basis, 1:k_cv))
  
  for (index in 1:k_cv) {
    # definujeme danou indexovou mnozinu
    fold <- folds[[index]]
      
    x.train.cv <- subset(X.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
      fdata()
    y.train.cv <- subset(Y.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
      as.numeric()
    
    x.test.cv <- subset(X.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
      fdata()
    y.test.cv <- subset(Y.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
      as.numeric()
    
    dataf <- as.data.frame(y.train.cv) 
    colnames(dataf) <- "Y"
    
    for (i in n.basis) {
      # baze pro bety
      basis2 <- create.bspline.basis(rangeval = rangeval, nbasis = i)
      
      basis.b <- list("x" = basis2)
      # vstupni data do modelu
      ldata <- list("df" = dataf, "x" = x.train.cv)
      # binomicky model ... model logisticke regrese
      model.glm <- fregre.glm(f, family = binomial(), data = ldata,
                              basis.x = basis.x, basis.b = basis.b)
      
      # presnost na validacni casti 
      newldata = list("df" = as.data.frame(y.test.cv), "x" = x.test.cv)
      predictions.valid <- predict(model.glm, newx = newldata)
      predictions.valid <- data.frame(Y.pred = ifelse(predictions.valid < 1/2, 0, 1))
      presnost.valid <- table(y.test.cv, predictions.valid$Y.pred) |>
        prop.table() |> diag() |> sum()
      
      # vlozime do matice
      CV.results[as.character(i), as.character(index)] <- presnost.valid
    } 
  }
  
  # spocitame prumerne presnosti pro jednotliva n pres folds
  CV.results <- apply(CV.results, 1, mean)
  n.basis.opt <- n.basis[which.max(CV.results)]
  CV_RESULTS$LR_func_n_basis[sim] <- n.basis.opt
  presnost.opt.cv <- max(CV.results)
  
  # optimalni model
  basis2 <- create.bspline.basis(rangeval = range(tt), nbasis = n.basis.opt)
  f <- Y ~ x
  # baze pro x a bety
  basis.x <- list("x" = basis1) 
  basis.b <- list("x" = basis2)
  # vstupni data do modelu
  dataf <- as.data.frame(y.train) 
  colnames(dataf) <- "Y"
  ldata <- list("df" = dataf, "x" = x.train)
  # binomicky model ... model logisticke regrese
  model.glm <- fregre.glm(f, family = binomial(), data = ldata,
                          basis.x = basis.x, basis.b = basis.b)
  
  # presnost na trenovacich datech
  predictions.train <- predict(model.glm, newx = ldata)
  predictions.train <- data.frame(Y.pred = ifelse(predictions.train < 1/2, 0, 1))
  presnost.train <- table(Y.train, predictions.train$Y.pred) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  newldata = list("df" = as.data.frame(Y.test), "x" = fdata(X.test))
  predictions.test <- predict(model.glm, newx = newldata)
  predictions.test <- data.frame(Y.pred = ifelse(predictions.test < 1/2, 0, 1))
  presnost.test <- table(Y.test, predictions.test$Y.pred) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'LR_functional', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 4.2) Logistická regrese s analýzou hlavních komponent
  
  # model
  clf.LR <- glm(Y ~  ., data = data.PCA.train, family = binomial)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.LR, newdata = data.PCA.train, type = 'response')
  predictions.train <- ifelse(predictions.train > 0.5, 1, 0)
  presnost.train <- table(data.PCA.train$Y, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.LR, newdata = data.PCA.test, type = 'response')
  predictions.test <- ifelse(predictions.test > 0.5, 1, 0)
  presnost.test <- table(data.PCA.test$Y, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'LR_score', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 5) Rozhodovací stromy
  ### 5.1) Diskretizace intervalu
  
  # posloupnost bodu, ve kterych funkce vyhodnotime
  t.seq <- seq(min(t), max(t), length = 101)
     
  grid.data <- eval.fd(fdobj = X.train, evalarg = t.seq)
  grid.data <- as.data.frame(t(grid.data)) # transpozice kvuli funkcim v radku
  grid.data$Y <- Y.train |> factor()
  
  grid.data.test <- eval.fd(fdobj = X.test, evalarg = t.seq)
  grid.data.test <- as.data.frame(t(grid.data.test))
  grid.data.test$Y <- Y.test |> factor()
  
  # sestrojeni modelu
  clf.tree <- train(Y ~ ., data = grid.data, 
                   method = "rpart", 
                   trControl = trainControl(method = "CV", number = 10),
                   metric = "Accuracy")
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.tree, newdata = grid.data)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.tree, newdata = grid.data.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'Tree_discr', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 5.2) Skóre hlavních komponent
  
  # sestrojeni modelu
  clf.tree.PCA <- train(Y ~ ., data = data.PCA.train, 
                   method = "rpart", 
                   trControl = trainControl(method = "CV", number = 10),
                   metric = "Accuracy")
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.tree.PCA, newdata = data.PCA.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.tree.PCA, newdata = data.PCA.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'Tree_score', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 5.3) Bázové koeficienty
  
  # trenovaci dataset
  data.Bbasis.train <- t(X.train$coefs) |> as.data.frame()
  data.Bbasis.train$Y <- factor(Y.train)
  
  # testovaci dataset
  data.Bbasis.test <- t(X.test$coefs) |> as.data.frame()
  data.Bbasis.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.tree.Bbasis <- train(Y ~ ., data = data.Bbasis.train, 
                   method = "rpart", 
                   trControl = trainControl(method = "CV", number = 10),
                   metric = "Accuracy")
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.tree.Bbasis, newdata = data.Bbasis.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.tree.Bbasis, newdata = data.Bbasis.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'Tree_Bbasis', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 6) Náhodné lesy
  
  ### 6.1) Diskretizace intervalu
  
  # sestrojeni modelu
  clf.RF <- randomForest(Y ~ ., data = grid.data, 
                         ntree = 500, # pocet stromu
                         importance = TRUE,
                         nodesize = 5)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.RF, newdata = grid.data)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.RF, newdata = grid.data.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'RF_discr', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 6.2) Skóre hlavních komponent
  
  # sestrojeni modelu
  clf.RF.PCA <- randomForest(Y ~ ., data = data.PCA.train, 
                             ntree = 500, # pocet stromu
                             importance = TRUE,
                             nodesize = 5)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.RF.PCA, newdata = data.PCA.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.RF.PCA, newdata = data.PCA.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'RF_score', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 6.3) Bázové koeficienty
  
  # sestrojeni modelu
  clf.RF.Bbasis <- randomForest(Y ~ ., data = data.Bbasis.train, 
                                ntree = 500, # pocet stromu
                                importance = TRUE,
                                nodesize = 5)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.RF.Bbasis, newdata = data.Bbasis.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.RF.Bbasis, newdata = data.Bbasis.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'RF_Bbasis', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 7) SVM
  
  ### 7.1) Diskretizace intervalu
  
  # rozdeleni na testovaci a trenovaci cast
  X.train_norm <- subset(XXfd_norm, split == TRUE)
  X.test_norm <- subset(XXfd_norm, split == FALSE)
  
  Y.train_norm <- subset(Y, split == TRUE)
  Y.test_norm <- subset(Y, split == FALSE)
  
  grid.data <- eval.fd(fdobj = X.train_norm, evalarg = t.seq)
  grid.data <- as.data.frame(t(grid.data)) 
  grid.data$Y <- Y.train_norm |> factor()
  
  grid.data.test <- eval.fd(fdobj = X.test_norm, evalarg = t.seq)
  grid.data.test <- as.data.frame(t(grid.data.test))
  grid.data.test$Y <- Y.test_norm |> factor()
  
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  # kontrola, ze mame opravdu k = k_cv
  while (length(folds) != k_cv) {
    folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  }
  
  # ktere hodnoty gamma chceme uvazovat
  gamma.cv <- 10^seq(-2, 2, length = 5)
  C.cv <- 10^seq(-3, 2, length = 5)
  p.cv <- 3
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
      cv_sample <- 1:dim(grid.data)[1] %in% fold
      
      data.grid.train.cv <- as.data.frame(grid.data[cv_sample, ])
      data.grid.test.cv <- as.data.frame(grid.data[!cv_sample, ])
      
      ## LINEARNI JADRO
      # sestrojeni modelu
      clf.SVM.l <- svm(Y ~ ., data = data.grid.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       kernel = 'linear')
      
      # presnost na validacnich datech
      predictions.test.l <- predict(clf.SVM.l, newdata = data.grid.test.cv)
      presnost.test.l <- table(data.grid.test.cv$Y, predictions.test.l) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane C a fold
      CV.results$SVM.l[(1:length(C.cv))[C.cv == C], 
                       index_cv] <- presnost.test.l
      
      ## POLYNOMIALNI JADRO
      for (p in p.cv) {
        # sestrojeni modelu
        clf.SVM.p <- svm(Y ~ ., data = data.grid.train.cv,
                         type = 'C-classification',
                         scale = TRUE,
                         cost = C,
                         coef0 = coef0,
                         degree = p,
                         kernel = 'polynomial')
        
        # presnost na validacnich datech
        predictions.test.p <- predict(clf.SVM.p, 
                                      newdata = data.grid.test.cv)
        presnost.test.p <- table(data.grid.test.cv$Y, predictions.test.p) |>
          prop.table() |> diag() |> sum()
        
        # presnosti vlozime na pozice pro dane C, p a fold
        CV.results$SVM.p[(1:length(C.cv))[C.cv == C], 
                         (1:length(p.cv))[p.cv == p],
                         index_cv] <- presnost.test.p
      }
          
      ## RADIALNI JADRO
      for (gamma in gamma.cv) {
        # sestrojeni modelu
        clf.SVM.r <- svm(Y ~ ., data = data.grid.train.cv,
                         type = 'C-classification',
                         scale = TRUE,
                         cost = C,
                         gamma = gamma,
                         kernel = 'radial')
        
        # presnost na validacnich datech
        predictions.test.r <- predict(clf.SVM.r, newdata = data.grid.test.cv)
        presnost.test.r <- table(data.grid.test.cv$Y, predictions.test.r) |>
          prop.table() |> diag() |> sum()
        
        # presnosti vlozime na pozice pro dane C, gamma a fold
        CV.results$SVM.r[(1:length(C.cv))[C.cv == C], 
                         (1:length(gamma.cv))[gamma.cv == gamma],
                         index_cv] <- presnost.test.r
      }
    }
  }
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
  
  # sestrojeni modelu
  clf.SVM.l <- svm(Y ~ ., data = grid.data,
                   type = 'C-classification',
                   scale = TRUE,
                   cost = C.opt[1],
                   kernel = 'linear')
  
  clf.SVM.p <- svm(Y ~ ., data = grid.data,
                   type = 'C-classification',
                   scale = TRUE,
                   cost = C.opt[2],
                   degree = p.opt,
                   coef0 = coef0,
                   kernel = 'polynomial')
  
  clf.SVM.r <- svm(Y ~ ., data = grid.data,
                   type = 'C-classification',
                   scale = TRUE, 
                   cost = C.opt[3],
                   gamma = gamma.opt,
                   kernel = 'radial')
  
  # presnost na trenovacich datech
  predictions.train.l <- predict(clf.SVM.l, newdata = grid.data)
  presnost.train.l <- table(Y.train, predictions.train.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.p <- predict(clf.SVM.p, newdata = grid.data)
  presnost.train.p <- table(Y.train, predictions.train.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.r <- predict(clf.SVM.r, newdata = grid.data)
  presnost.train.r <- table(Y.train, predictions.train.r) |>
    prop.table() |> diag() |> sum()
  
  # presnost na testovacich datech
  predictions.test.l <- predict(clf.SVM.l, newdata = grid.data.test)
  presnost.test.l <- table(Y.test, predictions.test.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.p <- predict(clf.SVM.p, newdata = grid.data.test)
  presnost.test.p <- table(Y.test, predictions.test.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.r <- predict(clf.SVM.r, newdata = grid.data.test)
  presnost.test.r <- table(Y.test, predictions.test.r) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = c('SVM linear - diskr', 
                              'SVM poly - diskr', 
                              'SVM rbf - diskr'), 
                    Err.train = 1 - c(presnost.train.l,
                                      presnost.train.p, presnost.train.r),
                    Err.test = 1 - c(presnost.test.l, 
                                     presnost.test.p, presnost.test.r))
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 7.2) Skóre hlavních komponent
  
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
                       coef0 = 1,
                       degree = p.opt,
                       kernel = 'polynomial')
  
  clf.SVM.r.PCA <- svm(Y ~ ., data = data.PCA.train,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C.opt[3],
                       gamma = gamma.opt,
                       kernel = 'radial')
  
  # presnost na trenovacich datech
  predictions.train.l <- predict(clf.SVM.l.PCA, newdata = data.PCA.train)
  presnost.train.l <- table(data.PCA.train$Y, predictions.train.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.p <- predict(clf.SVM.p.PCA, newdata = data.PCA.train)
  presnost.train.p <- table(data.PCA.train$Y, predictions.train.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.r <- predict(clf.SVM.r.PCA, newdata = data.PCA.train)
  presnost.train.r <- table(data.PCA.train$Y, predictions.train.r) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test.l <- predict(clf.SVM.l.PCA, newdata = data.PCA.test)
  presnost.test.l <- table(data.PCA.test$Y, predictions.test.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.p <- predict(clf.SVM.p.PCA, newdata = data.PCA.test)
  presnost.test.p <- table(data.PCA.test$Y, predictions.test.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.r <- predict(clf.SVM.r.PCA, newdata = data.PCA.test)
  presnost.test.r <- table(data.PCA.test$Y, predictions.test.r) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = c('SVM linear - PCA', 
                              'SVM poly - PCA', 
                              'SVM rbf - PCA'), 
                    Err.train = 1 - c(presnost.train.l, 
                                      presnost.train.p, presnost.train.r),
                    Err.test = 1 - c(presnost.test.l, 
                                     presnost.test.p, presnost.test.r))
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 7.3) Bázové koeficienty
  
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
      cv_sample <- 1:dim(grid.data)[1] %in% fold
      
      data.Bbasis.train.cv <- as.data.frame(data.Bbasis.train[cv_sample, ])
      data.Bbasis.test.cv <- as.data.frame(data.Bbasis.train[!cv_sample, ])
      
      ## LINEARNI JADRO
      # sestrojeni modelu
      clf.SVM.l <- svm(Y ~ ., data = data.Bbasis.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       kernel = 'linear')
      
      # presnost na validacnich datech
      predictions.test.l <- predict(clf.SVM.l, newdata = data.Bbasis.test.cv)
      presnost.test.l <- table(data.Bbasis.test.cv$Y, predictions.test.l) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane C a fold
      CV.results$SVM.l[(1:length(C.cv))[C.cv == C], 
                       index_cv] <- presnost.test.l
      
      ## POLYNOMIALNI JADRO
      for (p in p.cv) {
        # sestrojeni modelu
        clf.SVM.p <- svm(Y ~ ., data = data.Bbasis.train.cv,
                         type = 'C-classification',
                         scale = TRUE,
                         cost = C,
                         coef0 = coef0,
                         degree = p,
                         kernel = 'polynomial')
        
        # presnost na validacnich datech
        predictions.test.p <- predict(clf.SVM.p, 
                                      newdata = data.Bbasis.test.cv)
        presnost.test.p <- table(data.Bbasis.test.cv$Y, predictions.test.p) |>
          prop.table() |> diag() |> sum()
        
        # presnosti vlozime na pozice pro dane C, p a fold
        CV.results$SVM.p[(1:length(C.cv))[C.cv == C], 
                         (1:length(p.cv))[p.cv == p],
                         index_cv] <- presnost.test.p
      }
          
      ## RADIALNI JADRO
      for (gamma in gamma.cv) {
        # sestrojeni modelu
        clf.SVM.r <- svm(Y ~ ., data = data.Bbasis.train.cv,
                         type = 'C-classification',
                         scale = TRUE,
                         cost = C,
                         gamma = gamma,
                         kernel = 'radial')
        
        # presnost na validacnich datech
        predictions.test.r <- predict(clf.SVM.r, 
                                      newdata = data.Bbasis.test.cv)
        presnost.test.r <- table(data.Bbasis.test.cv$Y, predictions.test.r) |>
          prop.table() |> diag() |> sum()
        
        # presnosti vlozime na pozice pro dane C, gamma a fold
        CV.results$SVM.r[(1:length(C.cv))[C.cv == C], 
                         (1:length(gamma.cv))[gamma.cv == gamma],
                         index_cv] <- presnost.test.r
      }
    }
  }
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
  # sestrojeni modelu
  clf.SVM.l.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                          type = 'C-classification',
                          scale = TRUE,
                          cost = C.opt[1],
                          kernel = 'linear')
  
  clf.SVM.p.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                          type = 'C-classification',
                          scale = TRUE,
                          cost = C.opt[2],
                          coef0 = 1,
                          degree = p.opt,
                          kernel = 'polynomial')
  
  clf.SVM.r.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                          type = 'C-classification',
                          scale = TRUE,
                          cost = C.opt[3],
                          gamma = gamma.opt,
                          kernel = 'radial')
  
  # presnost na trenovacich datech
  predictions.train.l <- predict(clf.SVM.l.Bbasis, newdata = data.Bbasis.train)
  presnost.train.l <- table(Y.train, predictions.train.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.p <- predict(clf.SVM.p.Bbasis, newdata = data.Bbasis.train)
  presnost.train.p <- table(Y.train, predictions.train.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.r <- predict(clf.SVM.r.Bbasis, newdata = data.Bbasis.train)
  presnost.train.r <- table(Y.train, predictions.train.r) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test.l <- predict(clf.SVM.l.Bbasis, newdata = data.Bbasis.test)
  presnost.test.l <- table(Y.test, predictions.test.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.p <- predict(clf.SVM.p.Bbasis, newdata = data.Bbasis.test)
  presnost.test.p <- table(Y.test, predictions.test.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.r <- predict(clf.SVM.r.Bbasis, newdata = data.Bbasis.test)
  presnost.test.r <- table(Y.test, predictions.test.r) |>
    prop.table() |> diag() |> sum()

  Res <- data.frame(model = c('SVM linear - Bbasis', 
                              'SVM poly - Bbasis', 
                              'SVM rbf - Bbasis'), 
                    Err.train = 1 - c(presnost.train.l, 
                                      presnost.train.p, presnost.train.r),
                    Err.test = 1 - c(presnost.test.l, 
                                     presnost.test.p, presnost.test.r))

  RESULTS <- rbind(RESULTS, Res)
  
  ### 7.4) Projekce na B-splinovou bázi
  
  # hodnoty pro B-splinovou bazi
  rangeval <- range(t)
  norder <- 4
  n_basis_min <- norder
  n_basis_max <- 20
  dimensions <- n_basis_min:n_basis_max 
  folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  
  CV.results <- list(SVM.l = matrix(NA, nrow = length(dimensions), ncol = k_cv),
                     SVM.p = matrix(NA, nrow = length(dimensions), ncol = k_cv),
                     SVM.r = matrix(NA, nrow = length(dimensions), ncol = k_cv))
  
  for (d in dimensions) {
    bbasis <- create.bspline.basis(rangeval = rangeval, 
                                   nbasis = d)
    Projection <- project.basis(y = XX, argvals = t, basisobj = bbasis) 
    XX.train <- subset(t(Projection), split == TRUE)
    
    for (index_cv in 1:k_cv) {
      fold <- folds[[index_cv]]
      cv_sample <- 1:dim(XX.train)[1] %in% fold
      
      data.projection.train.cv <- as.data.frame(XX.train[cv_sample, ])
      data.projection.train.cv$Y <- factor(Y.train[cv_sample])
      data.projection.test.cv <- as.data.frame(XX.train[!cv_sample, ])
      Y.test.cv <- Y.train[!cv_sample]
      data.projection.test.cv$Y <- factor(Y.test.cv)
      # sestrojeni modelu
      clf.SVM.l.projection <- svm(Y ~ ., data = data.projection.train.cv,
                              type = 'C-classification',
                              scale = TRUE,
                              kernel = 'linear')
      
      clf.SVM.p.projection <- svm(Y ~ ., data = data.projection.train.cv,
                              type = 'C-classification',
                              scale = TRUE,
                              coef0 = 1,
                              kernel = 'polynomial')
      
      clf.SVM.r.projection <- svm(Y ~ ., data = data.projection.train.cv,
                              type = 'C-classification',
                              scale = TRUE,
                              kernel = 'radial')
      # presnost na validacnich datech
      ## linear kernel
      predictions.test.l <- predict(clf.SVM.l.projection,
                                    newdata = data.projection.test.cv)
      presnost.test.l <- table(Y.test.cv, predictions.test.l) |>
        prop.table() |> diag() |> sum()
      ## polynomial kernel
      predictions.test.p <- predict(clf.SVM.p.projection, 
                                    newdata = data.projection.test.cv)
      presnost.test.p <- table(Y.test.cv, predictions.test.p) |>
        prop.table() |> diag() |> sum()
      ## radial kernel
      predictions.test.r <- predict(clf.SVM.r.projection,
                                    newdata = data.projection.test.cv)
      presnost.test.r <- table(Y.test.cv, predictions.test.r) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane d a fold
      CV.results$SVM.l[d - min(dimensions) + 1, index_cv] <- presnost.test.l
      CV.results$SVM.p[d - min(dimensions) + 1, index_cv] <- presnost.test.p
      CV.results$SVM.r[d - min(dimensions) + 1, index_cv] <- presnost.test.r
    }
  }
  
  # spocitame prumerne presnosti pro jednotliva d pres folds
  for (n_method in 1:length(CV.results)) {
    CV.results[[n_method]] <- apply(CV.results[[n_method]], 1, mean)
  }
  
  d.opt <- c(which.max(CV.results$SVM.l) + n_basis_min - 1, 
             which.max(CV.results$SVM.p) + n_basis_min - 1, 
             which.max(CV.results$SVM.r) + n_basis_min - 1)
  
  # ulozime optimalni d do datove tabulky
  CV_RESULTS[sim, 4:6] <- d.opt
  
  # pripravime si datovou tabulku pro ulozeni vysledku
  Res <- data.frame(model = c('SVM linear - projection', 
                              'SVM poly - projection', 
                              'SVM rbf - projection'), 
                    Err.train = NA,
                    Err.test = NA)
  
  for (kernel_number in 1:3) {
    kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
    bbasis <- create.bspline.basis(rangeval = rangeval, 
                                   nbasis = d.opt[kernel_number])
    Projection <- project.basis(y = XX, argvals = t, basisobj = bbasis) 
    
    XX.train <- subset(t(Projection), split == TRUE)
    XX.test <- subset(t(Projection), split == FALSE)
    
    data.projection.train <- as.data.frame(XX.train)
    data.projection.train$Y <- factor(Y.train)
    
    data.projection.test <- as.data.frame(XX.test)
    data.projection.test$Y <- factor(Y.test)
    
    # sestrojeni modelu
    clf.SVM.projection <- svm(Y ~ ., data = data.projection.train,
                              type = 'C-classification',
                              scale = TRUE,
                              coef0 = 1,
                              kernel = kernel_type)
    
    # presnost na trenovacich datech
    predictions.train <- predict(clf.SVM.projection, newdata = data.projection.train)
    presnost.train <- table(Y.train, predictions.train) |>
      prop.table() |> diag() |> sum()
      
    # presnost na trenovacich datech
    predictions.test <- predict(clf.SVM.projection, newdata = data.projection.test)
    presnost.test <- table(Y.test, predictions.test) |>
      prop.table() |> diag() |> sum()
    
    # ulozeni vysledku
    Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
  }
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 7.5) SVM + RKHS
  
  C <- 1
  eps <- 0.01
  
  ### Gaussovo jadro
  
  # jadro a jadrova matice ... Gaussovske s parametrem gamma
  Gauss.kernel <- function(x, y, gamma) {
    return(exp(-gamma * norm(c(x - y) |> t(), type = 'F')))
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
  
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  
  # hodnoty hyperparametru, ktere budeme prochazet
  dimensions <- seq(2, 30, by = 2) # rozumny rozsah hodnot d
  gamma.cv <- 10^seq(-2, 2, length = 15)
  
  # list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
  # prazdna matice, do ktere vlozime jednotlive vysledky
  # ve sloupcich budou hodnoty presnosti pro dane
  # v radcich budou hodnoty pro danou gamma a vrstvy odpovidaji folds
  dim.names <- list(gamma = paste0('gamma:', round(gamma.cv, 3)),
                    d = paste0('d:', dimensions),
                    CV = paste0('cv:', 1:k_cv))
  
  CV.results <- list(
    SVM.l = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.p = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.r = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                  dimnames = dim.names))
  
  # samotna CV
  for (gamma in gamma.cv) {
    K <- Kernel.RKHS(t.seq, gamma = gamma)
    Eig <- eigen(K)
    eig.vals <- Eig$values
    eig.vectors <- Eig$vectors
    alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 
    
    # model
    for(i in 1:dim(data.RKHS)[2]) {
      df.svm <- data.frame(x = t.seq,
                           y = data.RKHS[, i])
      svm.RKHS <- svm(y ~ x, data = df.svm, 
                      kernel = 'radial',
                      type = 'eps-regression',
                      epsilon = eps, 
                      cost = C,
                      coef0 = 1,
                      gamma = gamma)
      alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
    }
    
    # projdeme dimenze
    for(d.RKHS in dimensions) {
      Lambda.RKHS <- matrix(NA, 
                            ncol = dim(data.RKHS)[2], 
                            nrow = d.RKHS) 
      # vypocet reprezentace
      for(l in 1:dim(data.RKHS)[2]) {
        Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                               alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
      }
      # projdeme folds
      for (index_cv in 1:k_cv) {
        # definice testovaci a trenovaci casti pro CV
        fold <- folds[[index_cv]]
        # rozdeleni na trenovaci a validacni data
        XX.train <- Lambda.RKHS[, fold]
        XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
        # pripravime si datovou tabulku pro ulozeni vysledku
        Res <- data.frame(model = c('SVM linear - RKHS', 
                                    'SVM poly - RKHS', 
                                    'SVM rbf - RKHS'), 
                          Err.test = NA)
        # projdeme jednotliva jadra
        for (kernel_number in 1:3) {
          kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
        
          data.RKHS.train <- as.data.frame(t(XX.train))
          data.RKHS.train$Y <- factor(Y.train[fold])
          
          data.RKHS.test <- as.data.frame(t(XX.test))
          data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
          
          # sestrojeni modelu
          clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                              type = 'C-classification',
                              scale = TRUE,
                              cost = C,
                              coef0 = 1,
                              kernel = kernel_type)
          
          # presnost na validacnich datech
          predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
          presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
            prop.table() |> diag() |> sum()
          
          # ulozeni vysledku
          Res[kernel_number, 2] <- 1 - presnost.test
        }
        # presnosti vlozime na pozice pro dane d, gamma a fold
        CV.results$SVM.l[paste0('gamma:', round(gamma, 3)), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[1, 2]
        CV.results$SVM.p[paste0('gamma:', round(gamma, 3)), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[2, 2]
        CV.results$SVM.r[paste0('gamma:', round(gamma, 3)), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[3, 2]
      }
    }
  }
  
  # spocitame prumerne presnosti pro jednotliva d pres folds
  for (n_method in 1:length(CV.results)) {
    CV.results[[n_method]] <- apply(CV.results[[n_method]], c(1, 2), mean)
  }
  
  gamma.opt <- c(which.min(CV.results$SVM.l) %% length(gamma.cv), 
                 which.min(CV.results$SVM.p) %% length(gamma.cv), 
                 which.min(CV.results$SVM.r) %% length(gamma.cv))
  gamma.opt[gamma.opt == 0] <- length(gamma.cv)
  gamma.opt <- gamma.cv[gamma.opt]
  
  d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
             which.min(t(CV.results$SVM.p)) %% length(dimensions), 
             which.min(t(CV.results$SVM.r)) %% length(dimensions))
  d.opt[d.opt == 0] <- length(dimensions)
  d.opt <- dimensions[d.opt]
  
  err.opt.cv <- c(min(CV.results$SVM.l), 
                       min(CV.results$SVM.p),
                       min(CV.results$SVM.r))
  df.RKHS.res <- data.frame(d = d.opt, gamma = gamma.opt, CV = err.opt.cv,
             Kernel = c('linear', 'polynomial', 'radial') |> factor(),
             row.names = c('linear', 'poly', 'radial'))
  
  CV_RESULTS[sim, 7:9] <- gamma.opt
  CV_RESULTS[sim, 10:12] <- d.opt
  
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  # pridame i testovaci data
  data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
  # pripravime si datovou tabulku pro ulozeni vysledku
  Res <- data.frame(model = c('SVM linear - RKHS - radial', 
                              'SVM poly - RKHS - radial', 
                              'SVM rbf - RKHS - radial'), 
                    Err.train = NA,
                    Err.test = NA)
  
  # projdeme jednotliva jadra
  for (kernel_number in 1:3) {
    # spocitame matici K
    gamma <- gamma.opt[kernel_number] # hodnota gamma pomoci CV
    K <- Kernel.RKHS(t.seq, gamma = gamma)
    
    # urcime vlastni cisla a vektory
    Eig <- eigen(K)
    eig.vals <- Eig$values
    eig.vectors <- Eig$vectors
    # urceni koeficientu alpha z SVM
    alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                         ncol = dim(data.RKHS)[2]) # prazdny objekt
    
    # model
    for(i in 1:dim(data.RKHS)[2]) {
      df.svm <- data.frame(x = t.seq,
                           y = data.RKHS[, i])
      svm.RKHS <- svm(y ~ x, data = df.svm, 
                      kernel = 'radial',
                      type = 'eps-regression',
                      epsilon = eps,
                      cost = C,
                      gamma = gamma)
      # urceni alpha
      alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
    }
    # d
    d.RKHS <- d.opt[kernel_number]
    
    # urceni vektoru lambda
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) # vytvoreni prazdneho objektu
    
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    
    # rozdeleni na trenovaci a testovaci data
    XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
    XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]
  
    kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
  
    data.RKHS.train <- as.data.frame(t(XX.train))
    data.RKHS.train$Y <- factor(Y.train)
    
    data.RKHS.test <- as.data.frame(t(XX.test))
    data.RKHS.test$Y <- factor(Y.test)
    
    # sestrojeni modelu
    clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                        type = 'C-classification',
                        scale = TRUE,
                        cost = C,
                        coef0 = 1,
                        kernel = kernel_type)
    
    # presnost na trenovacich datech
    predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
    presnost.train <- table(Y.train, predictions.train) |>
      prop.table() |> diag() |> sum()
      
    # presnost na testovacich datech
    predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
    presnost.test <- table(Y.test, predictions.test) |>
      prop.table() |> diag() |> sum()
    
    # ulozeni vysledku
    Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
  }
  
  RESULTS <- rbind(RESULTS, Res)

  ### Polynomialni jadro
  
  # jadro a jadrova matice ... polynomialni s parametrem p
  Poly.kernel <- function(x, y, p) {
    return((1 + x * y)^p)
  }
  
  Kernel.RKHS <- function(x, p) {
    K <- matrix(NA, ncol = length(x), nrow = length(x))
    for(i in 1:nrow(K)) {
      for(j in 1:ncol(K)) {
        K[i, j] <- Poly.kernel(x = x[i], y = x[j], p)
      }
    }
    return(K)
  }
  
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  
  # hodnoty hyperparametru, ktere budeme prochazet
  dimensions <- seq(2, 10, by = 1) # rozumny rozsah hodnot d
  poly.cv <- 2:5
  
  # list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
  # prazdna matice, do ktere vlozime jednotlive vysledky
  # ve sloupcich budou hodnoty presnosti pro dane
  # v radcich budou hodnoty pro dane p a vrstvy odpovidaji folds
  dim.names <- list(p = paste0('p:', poly.cv),
                    d = paste0('d:', dimensions),
                    CV = paste0('cv:', 1:k_cv))
  
  CV.results <- list(
    SVM.l = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.p = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.r = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                  dimnames = dim.names))
  
  # samotna CV
  for (p in poly.cv) {
    K <- Kernel.RKHS(t.seq, p = p)
    Eig <- eigen(K)
    eig.vals <- Eig$values
    eig.vectors <- Eig$vectors
    alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 
    
    # model
    for(i in 1:dim(data.RKHS)[2]) {
      df.svm <- data.frame(x = t.seq,
                           y = data.RKHS[, i])
      svm.RKHS <- svm(y ~ x, data = df.svm, 
                      kernel = 'polynomial',
                      type = 'eps-regression',
                      epsilon = eps,                   
                      coef0 = 1,
                      gamma = 1,
                      cost = C,
                      degree = p)
      alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
    }
    
    # projdeme dimenze
    for(d.RKHS in dimensions) {
      Lambda.RKHS <- matrix(NA, 
                            ncol = dim(data.RKHS)[2], 
                            nrow = d.RKHS) 
      # vypocet reprezentace
      for(l in 1:dim(data.RKHS)[2]) {
        Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                               alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
      }
      # projdeme folds
      for (index_cv in 1:k_cv) {
        # definice testovaci a trenovaci casti pro CV
        fold <- folds[[index_cv]]
        # rozdeleni na trenovaci a validacni data
        XX.train <- Lambda.RKHS[, fold]
        XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
        # pripravime si datovou tabulku pro ulozeni vysledku
        Res <- data.frame(model = c('SVM linear - RKHS', 
                                    'SVM poly - RKHS', 
                                    'SVM rbf - RKHS'), 
                          Err.test = NA)
        # projdeme jednotliva jadra
        for (kernel_number in 1:3) {
          kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
        
          data.RKHS.train <- as.data.frame(t(XX.train))
          data.RKHS.train$Y <- factor(Y.train[fold])
          
          data.RKHS.test <- as.data.frame(t(XX.test))
          data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
          
          # sestrojeni modelu
          clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                              type = 'C-classification',
                              scale = TRUE,
                              epsilon = eps,                   
                              coef0 = 1,
                              gamma = 1,
                              cost = C,
                              kernel = kernel_type)
          
          # presnost na validacnich datech
          predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
          presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
            prop.table() |> diag() |> sum()
          
          # ulozeni vysledku
          Res[kernel_number, 2] <- 1 - presnost.test
        }
        # presnosti vlozime na pozice pro dane d, gamma a fold
        CV.results$SVM.l[paste0('p:', p), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[1, 2]
        CV.results$SVM.p[paste0('p:', p), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[2, 2]
        CV.results$SVM.r[paste0('p:', p), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[3, 2]
      }
    }
  }
  
  # spocitame prumerne presnosti pro jednotliva d pres folds
  for (n_method in 1:length(CV.results)) {
    CV.results[[n_method]] <- apply(CV.results[[n_method]], c(1, 2), mean)
  }
  
  poly.opt <- c(which.min(CV.results$SVM.l) %% length(poly.cv), 
                 which.min(CV.results$SVM.p) %% length(poly.cv), 
                 which.min(CV.results$SVM.r) %% length(poly.cv))
  poly.opt[poly.opt == 0] <- length(poly.cv)
  poly.opt <- poly.cv[poly.opt]
  
  d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
             which.min(t(CV.results$SVM.p)) %% length(dimensions), 
             which.min(t(CV.results$SVM.r)) %% length(dimensions))
  d.opt[d.opt == 0] <- length(dimensions)
  d.opt <- dimensions[d.opt]
  
  err.opt.cv <- c(min(CV.results$SVM.l), 
                       min(CV.results$SVM.p),
                       min(CV.results$SVM.r))
  df.RKHS.res <- data.frame(d = d.opt, p = poly.opt, CV = err.opt.cv,
             Kernel = c('linear', 'polynomial', 'radial') |> factor(),
             row.names = c('linear', 'poly', 'radial'))
  
  CV_RESULTS[sim, 13:15] <- poly.opt
  CV_RESULTS[sim, 16:18] <- d.opt
  
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  # pridame i testovaci data
  data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
  
  # pripravime si datovou tabulku pro ulozeni vysledku
  Res <- data.frame(model = c('SVM linear - RKHS - poly', 
                              'SVM poly - RKHS - poly', 
                              'SVM rbf - RKHS - poly'), 
                    Err.train = NA,
                    Err.test = NA)
  
  # projdeme jednotliva jadra
  for (kernel_number in 1:3) {
    # spocitame matici K
    p <- poly.opt[kernel_number] # hodnota gamma pomoci CV
    K <- Kernel.RKHS(t.seq, p = p)
    
    # urcime vlastni cisla a vektory
    Eig <- eigen(K)
    eig.vals <- Eig$values
    eig.vectors <- Eig$vectors
    # urceni koeficientu alpha z SVM
    alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                         ncol = dim(data.RKHS)[2]) # prazdny objekt
    
    # model
    for(i in 1:dim(data.RKHS)[2]) {
      df.svm <- data.frame(x = t.seq,
                           y = data.RKHS[, i])
      svm.RKHS <- svm(y ~ x, data = df.svm, 
                      kernel = 'polynomial',
                      type = 'eps-regression',
                      epsilon = eps,                   
                      coef0 = 1,
                      gamma = 1,
                      cost = C,
                      degree = p)
      # urceni alpha
      alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
    }
    # d
    d.RKHS <- d.opt[kernel_number]
    
    # urceni vektoru lambda
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) # vytvoreni prazdneho objektu
    
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    
    # rozdeleni na trenovaci a testovaci data
    XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
    XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]
  
    kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
  
    data.RKHS.train <- as.data.frame(t(XX.train))
    data.RKHS.train$Y <- factor(Y.train)
    
    data.RKHS.test <- as.data.frame(t(XX.test))
    data.RKHS.test$Y <- factor(Y.test)
    
    # sestrojeni modelu
    clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                        type = 'C-classification',
                        scale = TRUE,
                        epsilon = eps,                   
                        coef0 = 1,
                        gamma = 1,
                        cost = C,
                        kernel = kernel_type)
    
    # presnost na trenovacich datech
    predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
    presnost.train <- table(Y.train, predictions.train) |>
      prop.table() |> diag() |> sum()
      
    # presnost na testovacich datech
    predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
    presnost.test <- table(Y.test, predictions.test) |>
      prop.table() |> diag() |> sum()
    
    # ulozeni vysledku
    Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
  }
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### Linearni jadro
  
  # jadro a jadrova matice ... polynomialni s parametrem p
  Linear.kernel <- function(x, y) {
    return(x * y)
  }
  
  Kernel.RKHS <- function(x) {
    K <- matrix(NA, ncol = length(x), nrow = length(x))
    for(i in 1:nrow(K)) {
      for(j in 1:ncol(K)) {
        K[i, j] <- Linear.kernel(x = x[i], y = x[j])
      }
    }
    return(K)
  }
  
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  
  # hodnoty hyperparametru, ktere budeme prochazet
  dimensions <- seq(3, 40, by = 2) # rozumny rozsah hodnot d
  
  # list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
  # prazdna matice, do ktere vlozime jednotlive vysledky
  # ve sloupcich budou hodnoty presnosti pro dane d
  # v radcich budou hodnoty pro vrstvy odpovidaji folds
  dim.names <- list(d = paste0('d:', dimensions),
                    CV = paste0('cv:', 1:k_cv))
  
  CV.results <- list(
    SVM.l = array(NA, dim = c(length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.p = array(NA, dim = c(length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.r = array(NA, dim = c(length(dimensions), k_cv),
                  dimnames = dim.names))
  
  # samotna CV
  K <- Kernel.RKHS(t.seq)
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'linear',
                    type = 'eps-regression',
                    cost = C,
                    epsilon = eps)
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
  }
  
  # projdeme dimenze
  for(d.RKHS in dimensions) {
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) 
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                             alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    # projdeme folds
    for (index_cv in 1:k_cv) {
      # definice testovaci a trenovaci casti pro CV
      fold <- folds[[index_cv]]
      # rozdeleni na trenovaci a validacni data
      XX.train <- Lambda.RKHS[, fold]
      XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
      # pripravime si datovou tabulku pro ulozeni vysledku
      Res <- data.frame(model = c('SVM linear - RKHS', 
                                  'SVM poly - RKHS', 
                                  'SVM rbf - RKHS'), 
                        Err.test = NA)
      # projdeme jednotliva jadra
      for (kernel_number in 1:3) {
        kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
      
        data.RKHS.train <- as.data.frame(t(XX.train))
        data.RKHS.train$Y <- factor(Y.train[fold])
        
        data.RKHS.test <- as.data.frame(t(XX.test))
        data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
        
        # sestrojeni modelu
        clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                            type = 'C-classification',
                            scale = TRUE,
                            cost = C,
                            coef0 = 1,
                            kernel = kernel_type)
        
        # presnost na validacnich datech
        predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
        presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
          prop.table() |> diag() |> sum()
        
        # ulozeni vysledku
        Res[kernel_number, 2] <- 1 - presnost.test
      }
      # presnosti vlozime na pozice pro dane d, gamma a fold
      CV.results$SVM.l[paste0('d:', d.RKHS), 
                       index_cv] <- Res[1, 2]
      CV.results$SVM.p[paste0('d:', d.RKHS), 
                       index_cv] <- Res[2, 2]
      CV.results$SVM.r[paste0('d:', d.RKHS), 
                       index_cv] <- Res[3, 2]
    }
  }
  
  # spocitame prumerne presnosti pro jednotliva d pres folds
  for (n_method in 1:length(CV.results)) {
    CV.results[[n_method]] <- apply(CV.results[[n_method]], 1, mean)
  }
  
  d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
             which.min(t(CV.results$SVM.p)) %% length(dimensions), 
             which.min(t(CV.results$SVM.r)) %% length(dimensions))
  d.opt[d.opt == 0] <- length(dimensions)
  d.opt <- dimensions[d.opt]
  
  err.opt.cv <- c(min(CV.results$SVM.l), 
                       min(CV.results$SVM.p),
                       min(CV.results$SVM.r))
  df.RKHS.res <- data.frame(d = d.opt, CV = err.opt.cv,
             Kernel = c('linear', 'polynomial', 'radial') |> factor(),
             row.names = c('linear', 'poly', 'radial'))
  
  CV_RESULTS[sim, 19:21] <- d.opt
  
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  # pridame i testovaci data
  data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
  
  # pripravime si datovou tabulku pro ulozeni vysledku
  Res <- data.frame(model = c('SVM linear - RKHS - linear', 
                              'SVM poly - RKHS - linear', 
                              'SVM rbf - RKHS - linear'), 
                    Err.train = NA,
                    Err.test = NA)
  
  # projdeme jednotliva jadra
  for (kernel_number in 1:3) {
    # spocitame matici K
    K <- Kernel.RKHS(t.seq)
    
    # urcime vlastni cisla a vektory
    Eig <- eigen(K)
    eig.vals <- Eig$values
    eig.vectors <- Eig$vectors
    # urceni koeficientu alpha z SVM
    alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                         ncol = dim(data.RKHS)[2]) # prazdny objekt
    
    # model
    for(i in 1:dim(data.RKHS)[2]) {
      df.svm <- data.frame(x = t.seq,
                           y = data.RKHS[, i])
      svm.RKHS <- svm(y ~ x, data = df.svm, 
                      kernel = 'linear',
                      type = 'eps-regression',
                      cost = C,
                      epsilon = eps)
      # urceni alpha
      alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
    }
    # d
    d.RKHS <- d.opt[kernel_number]
    
    # urceni vektoru lambda
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) # vytvoreni prazdneho objektu
    
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    
    # rozdeleni na trenovaci a testovaci data
    XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
    XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]
  
    kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
  
    data.RKHS.train <- as.data.frame(t(XX.train))
    data.RKHS.train$Y <- factor(Y.train)
    
    data.RKHS.test <- as.data.frame(t(XX.test))
    data.RKHS.test$Y <- factor(Y.test)
    
    # sestrojeni modelu
    clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                        type = 'C-classification',
                        scale = TRUE,
                        cost = C,
                        coef0 = 1,
                        kernel = kernel_type)
    
    # presnost na trenovacich datech
    predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
    presnost.train <- table(Y.train, predictions.train) |>
      prop.table() |> diag() |> sum()
      
    # presnost na testovacich datech
    predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
    presnost.test <- table(Y.test, predictions.test) |>
      prop.table() |> diag() |> sum()
    
    # ulozeni vysledku
    Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
  }
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## pridame vysledky do objektu SIMULACE
  
  SIMULACE$train[sim, ] <- RESULTS$Err.train
  SIMULACE$test[sim, ] <- RESULTS$Err.test
  
  cat('\r', sim)
}

# ulozime vysledne hodnoty 
save(SIMULACE, CV_RESULTS, file = 'RData/aplikace_03neder.RData')
```

Nyní spočítáme průměrné testovací a trénovací chybovosti pro jednotlivé klasifikační metody.


```r
# dame do vysledne tabulky

SIMULACE.df <- data.frame(Err.train = apply(SIMULACE$train, 2, mean),
                          Err.test = apply(SIMULACE$test, 2, mean),
                          SD.train = apply(SIMULACE$train, 2, sd),
                          SD.test = apply(SIMULACE$test, 2, sd))

# ulozime vysledne hodnoty 
save(SIMULACE.df, file = 'RData/aplikace_03neder_res.RData')
```

#### Výsledky




Table: (\#tab:unnamed-chunk-229)Souhrnné výsledky použitých metod na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti, $\widehat{Err}_{test}$ testovací chybovosti, $\widehat{SD}_{train}$ odhad směrodatné odchylky trénovacích chybovostí a $\widehat{SD}_{test}$ je odhad směrodatné odchylky testovacích chybovostí.

                              $\widehat{Err}_{train}$   $\widehat{Err}_{test}$   $\widehat{SD}_{train}$   $\widehat{SD}_{test}$
---------------------------  ------------------------  -----------------------  -----------------------  ----------------------
KNN                                            0.1730                   0.1777                   0.0265                  0.0504
LDA                                            0.2996                   0.3163                   0.0207                  0.0446
QDA                                            0.3017                   0.3182                   0.0217                  0.0441
LR_functional                                  0.0105                   0.0482                   0.0278                  0.0467
LR_score                                       0.2919                   0.3103                   0.0207                  0.0447
Tree_discr                                     0.1857                   0.2937                   0.0449                  0.0583
Tree_score                                     0.2470                   0.3360                   0.0451                  0.0574
Tree_Bbasis                                    0.1855                   0.2918                   0.0456                  0.0565
RF_discr                                       0.0121                   0.2072                   0.0091                  0.0502
RF_score                                       0.0360                   0.3100                   0.0110                  0.0507
RF_Bbasis                                      0.0121                   0.2066                   0.0087                  0.0453
SVM linear - diskr                             0.0036                   0.0162                   0.0058                  0.0227
SVM poly - diskr                               0.0113                   0.0455                   0.0129                  0.0234
SVM rbf - diskr                                0.0053                   0.0346                   0.0063                  0.0229
SVM linear - PCA                               0.2989                   0.3285                   0.0234                  0.0495
SVM poly - PCA                                 0.2833                   0.3515                   0.0329                  0.0452
SVM rbf - PCA                                  0.1553                   0.3469                   0.1128                  0.0479
SVM linear - Bbasis                            0.0113                   0.0277                   0.0074                  0.0194
SVM poly - Bbasis                              0.0137                   0.0495                   0.0104                  0.0287
SVM rbf - Bbasis                               0.0153                   0.0532                   0.0063                  0.0252
SVM linear - projection                        0.0312                   0.0425                   0.0103                  0.0244
SVM poly - projection                          0.0339                   0.0577                   0.0135                  0.0331
SVM rbf - projection                           0.1395                   0.2034                   0.0289                  0.0554
SVM linear - RKHS - radial                     0.0008                   0.0197                   0.0024                  0.0188
SVM poly - RKHS - radial                       0.0009                   0.0132                   0.0023                  0.0178
SVM rbf - RKHS - radial                        0.0036                   0.0195                   0.0047                  0.0163
SVM linear - RKHS - poly                       0.0543                   0.0832                   0.0138                  0.0321
SVM poly - RKHS - poly                         0.0305                   0.0889                   0.0152                  0.0305
SVM rbf - RKHS - poly                          0.0302                   0.0642                   0.0094                  0.0294
SVM linear - RKHS - linear                     0.0448                   0.0725                   0.0149                  0.0349
SVM poly - RKHS - linear                       0.0433                   0.0754                   0.0134                  0.0358
SVM rbf - RKHS - linear                        0.0803                   0.1188                   0.0217                  0.0434

V tabulce výše jsou uvedeny všechny vypočtené charakteristiky.
Jsou zde uvedeny také směrodatné odchylky, abychom mohli porovnat jakousi stálost či míru variability jednotlivých metod.


```r
wilcox.test(SIMULACE$test[, 'SVM poly - RKHS - radial'], SIMULACE$test[, 'SVM linear - diskr'], alternative = 't', paired = T)$p.value
```

```
## [1] 0.3874749
```

```r
wilcox.test(SIMULACE$test[, 'SVM rbf - RKHS - radial'], SIMULACE$test[, 'SVM linear - diskr'], alternative = 't', paired = T)$p.value
```

```
## [1] 0.201061
```

```r
wilcox.test(SIMULACE$test[, 'SVM linear - RKHS - radial'], SIMULACE$test[, 'SVM linear - diskr'], alternative = 't', paired = T)$p.value
```

```
## [1] 0.1067453
```

```r
wilcox.test(SIMULACE$test[, 'LR_functional'], SIMULACE$test[, 'SVM linear - diskr'], alternative = 'greater', paired = T)$p.value
```

```
## [1] 2.09865e-11
```

```r
wilcox.test(SIMULACE$test[, 'LR_functional'], SIMULACE$test[, 'SVM poly - RKHS - radial'], alternative = 'greater', paired = T)$p.value
```

```
## [1] 2.355944e-13
```

```r
wilcox.test(SIMULACE$test[, 'LR_functional'], SIMULACE$test[, 'SVM rbf - RKHS - radial'], alternative = 'greater', paired = T)$p.value
```

```
## [1] 1.332949e-10
```


Nakonec ještě můžeme graficky zobrazit vypočtené hodnoty ze simulace pro jednotlivé klasifikační metody pomocí krabicových diagramů, zvlášť pro testovací a trénovací chybovosti.


```r
# nastavime jinak nazvy klasifikacnich metod
methods_names <- c(
      '$K$ nejbližších sousedů',
      'Lineární diskriminační analýza',
      'Kvadratická diskriminační analýza',
      'Funkcionální logistická regrese',
      'Logistické regrese s fPCA',
      'Rozhodovací strom -- diskretizace',
      'Rozhodovací strom -- fPCA',
      'Rozhodovací strom -- bázové koeficienty',
      'Náhodný les -- diskretizace',
      'Náhodný les -- fPCA',
      'Náhodný les -- bázové koeficienty',
      'SVM (linear) -- diskretizace',
      'SVM (poly) -- diskretizace',
      'SVM (radial) -- diskretizace',
      'SVM (linear) -- fPCA',
      'SVM (poly) -- fPCA',
      'SVM (radial) -- fPCA',
      'SVM (linear) -- bázové koeficienty',
      'SVM (poly) -- bázové koeficienty',
      'SVM (radial) -- bázové koeficienty',
      'SVM (linear) -- projekce',
      'SVM (poly) -- projekce',
      'SVM (radial) -- projekce',
      'RKHS (radial SVR) $+$ SVM (linear)',
      'RKHS (radial SVR) $+$ SVM (poly)',
      'RKHS (radial SVR) $+$ SVM (radial)',
      'RKHS (poly SVR) $+$ SVM (linear)',
      'RKHS (poly SVR) $+$ SVM (poly)',
      'RKHS (poly SVR) $+$ SVM (radial)',
      'RKHS (linear SVR) $+$ SVM (linear)',
      'RKHS (linear SVR) $+$ SVM (poly)',
      'RKHS (linear SVR) $+$ SVM (radial)'
)


# barvy pro boxploty 
box_col <- c('#4dd2ff', '#0099cc', '#00ace6', '#00bfff',
             '#1ac5ff', rep('#33ccff', 3), rep('#0086b3', 3),
             rep('#ff3814', 3), rep('#ff6347', 3), rep('#ff7961', 3),
             rep('#ff4d2e', 3), rep('#fa2600', 9))

# box_col <- c('#CA0A0A', '#fa2600', '#fa2600', '#D15804',
#              '#D15804', rep('#D3006D', 3), rep('#BE090F', 3), c("#12DEE8", "#4ECBF3", "#127DE8", "#4C3CD3", "#4E65F3", "#4E9EF3", "#081D58") |> rep(each = 3))

# alpha pro boxploty
box_alpha <- c(0.9, 0.9, 0.8, 0.9, 0.8, 0.9, 0.8, 0.7, 0.9, 0.8, 0.7,
               0.9, 0.8, 0.7, 0.9, 0.8, 0.7, 0.9, 0.8, 0.7, 0.9, 0.8, 0.7,
               seq(0.9, 0.6, length = 9)) - 0.3
```



```r
# pro trenovaci data
SIMULACE$train |> 
  pivot_longer(cols = methods, names_to = 'method', values_to = 'Err') |>
  mutate(method = factor(method, levels = methods, labels = methods, ordered = TRUE)) |> 
  as.data.frame() |>
  ggplot(aes(x = method, y = Err, fill = method, colour = method, alpha = 0.3)) + 
  geom_boxplot(outlier.colour = "white", outlier.shape = 16, outlier.size = 0, 
               notch = FALSE, colour = 'black') + 
  theme_bw() + 
  labs(x = 'Klasifikační metoda',
       y = expression(widehat(Err)[train])) + 
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_jitter(position = position_jitter(0.15), alpha = 0.7, size = 1, pch = 21,
              colour = 'black') +
  stat_summary(fun = "mean", geom = "point", shape = '+',
               size = 4, color = "black", alpha = 0.9)+ 
  geom_hline(yintercept = min(SIMULACE.df$Err.train), 
             linetype = 'dashed', colour = 'grey') 
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-232-1.png" alt="Krabicové diagramy trénovacích chybovostí pro 100 simulací zvlášť pro jednotlivé klasifikační metody. Černými symboly $+$ jsou vyznačeny průměry." width="672" />
<p class="caption">(\#fig:unnamed-chunk-232)Krabicové diagramy trénovacích chybovostí pro 100 simulací zvlášť pro jednotlivé klasifikační metody. Černými symboly $+$ jsou vyznačeny průměry.</p>
</div>




```r
# pro testovaci data
SIMULACE$test |> 
  pivot_longer(cols = methods, names_to = 'method', values_to = 'Err') |>
  mutate(method = factor(method, levels = methods, labels = methods, ordered = TRUE)) |> 
  as.data.frame() |>
  ggplot(aes(x = method, y = Err, fill = method, colour = method, alpha = method)) + 
  geom_hline(yintercept = min(SIMULACE.df$Err.test), 
             linetype = 'dashed', colour = 'grey') +
  geom_boxplot(outlier.colour = "white", outlier.shape = 16, outlier.size = 0, 
               notch = FALSE, colour = 'black') + 
  theme_bw() + 
  labs(x = 'Klasifikační metoda',
       # y = "$\\widehat{\\textnormal{Err}}_{test}$"
       y = expression(widehat(Err)[test])
       ) + 
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  geom_jitter(position = position_jitter(0.15), alpha = 0.6, size = 0.9, pch = 21,
              colour = "black") +
  stat_summary(fun = "mean", geom = "point", shape = '+',
               size = 3, color = "black", alpha = 0.9)# +
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-234-1.png" alt="Krabicové diagramy testovacích chybovostí pro 100 simulací zvlášť pro jednotlivé klasifikační metody. Černými symboly $+$ jsou vyznačeny průměry." width="672" />
<p class="caption">(\#fig:unnamed-chunk-234)Krabicové diagramy testovacích chybovostí pro 100 simulací zvlášť pro jednotlivé klasifikační metody. Černými symboly $+$ jsou vyznačeny průměry.</p>
</div>

```r
  # scale_x_discrete(labels = methods_names) +
  # theme(plot.margin = unit(c(0.5, 0.5, 2, 2), "cm")) +
  # scale_fill_manual(values = box_col) +
  # scale_alpha_manual(values = box_alpha)

# ggsave("figures/kap7_tecator_box_test_neder.tex", device = tikz, width = 9, height = 7)
```



Nakonec se podívejme, jaké hodnoty hyperparametrů byly nejčastější volbou.


Table: (\#tab:unnamed-chunk-236)Mediány hodnot hyperparametrů pro vybrané metody, u nichž se určoval nějaký hyperparametr pomocí cross-validace.

                          Mediánová hodnota hyperparametru
-----------------------  ---------------------------------
KNN_K                                                    1
nharm                                                    1
LR_func_n_basis                                         19
SVM_d_Linear                                             6
SVM_d_Poly                                               6
SVM_d_Radial                                             6
SVM_RKHS_radial_gamma1                                   1
SVM_RKHS_radial_gamma2                                   1
SVM_RKHS_radial_gamma3                                   1
SVM_RKHS_radial_d1                                      17
SVM_RKHS_radial_d2                                      16
SVM_RKHS_radial_d3                                      20
SVM_RKHS_poly_p1                                         4
SVM_RKHS_poly_p2                                         4
SVM_RKHS_poly_p3                                         5
SVM_RKHS_poly_d1                                         8
SVM_RKHS_poly_d2                                         7
SVM_RKHS_poly_d3                                         7
SVM_RKHS_linear_d1                                      20
SVM_RKHS_linear_d2                                      17
SVM_RKHS_linear_d3                                      13


```r
CV_res <- CV_RESULTS |> 
  pivot_longer(cols = CV_RESULTS |> colnames(), names_to = 'method', values_to = 'hyperparameter') |>
  mutate(method = factor(method, 
                         levels = CV_RESULTS |> colnames(), 
                         labels = CV_RESULTS |> colnames(), ordered = TRUE)) |> 
  as.data.frame() 

CV_res |> 
  filter(method %in% c('KNN_K', 'nharm', 'LR_func_n_basis')) |>
  ggplot(aes(x = hyperparameter, #y = after_stat(density),
             fill = method, colour = method)) + 
  geom_histogram(binwidth = 1, alpha = 0.6) + 
  theme_bw() + 
  facet_grid(~method, scales = 'free') +
  labs(x = 'Hodnoty hyperparametru',
       y = 'Absolutní počet') + 
  theme(legend.position = 'none')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-237-1.png" alt="Histogramy hodnot hyperparametrů pro KNN, funkcionální logistickou regresi a také histogram pro počet hlavních komponent." width="672" />
<p class="caption">(\#fig:unnamed-chunk-237)Histogramy hodnot hyperparametrů pro KNN, funkcionální logistickou regresi a také histogram pro počet hlavních komponent.</p>
</div>




```r
CV_res |> 
  filter(method %in% c('SVM_d_Linear', 'SVM_d_Poly', 'SVM_d_Radial')) |>
  ggplot(aes(x = hyperparameter, #y = after_stat(density),
             fill = method, colour = method)) + 
  geom_histogram(binwidth = 1, alpha = 0.6) + 
  theme_bw() + 
  facet_grid(~method, scales = 'free') +
  labs(x = 'Hodnoty hyperparametru',
       y = 'Absolutní počet') + 
  theme(legend.position = 'none')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-239-1.png" alt="Histogramy hodnot hyperparametrů metody SVM s projekcí na B-splinovou bázi." width="672" />
<p class="caption">(\#fig:unnamed-chunk-239)Histogramy hodnot hyperparametrů metody SVM s projekcí na B-splinovou bázi.</p>
</div>




```r
CV_res |> 
  filter(method %in% c('SVM_RKHS_radial_gamma1', 'SVM_RKHS_radial_gamma2',
                       'SVM_RKHS_radial_gamma3', 'SVM_RKHS_radial_d1', 
                       'SVM_RKHS_radial_d2', 'SVM_RKHS_radial_d3')) |>
  ggplot(aes(x = hyperparameter, #y = after_stat(density),
             fill = method, colour = method)) + 
  geom_histogram(bins = 10, alpha = 0.6) + 
  theme_bw() + 
  facet_wrap(~method, scales = 'free', ncol = 3) +
  labs(x = 'Hodnoty hyperparametru',
       y = 'Absolutní počet') + 
  theme(legend.position = 'none')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-241-1.png" alt="Histogramy hodnot hyperparametrů metody RKHS + SVM s radiálním jádrem." width="672" />
<p class="caption">(\#fig:unnamed-chunk-241)Histogramy hodnot hyperparametrů metody RKHS + SVM s radiálním jádrem.</p>
</div>




```r
CV_res |> 
  filter(method %in% c('SVM_RKHS_poly_p1', 'SVM_RKHS_poly_p2',
                       'SVM_RKHS_poly_p3', 'SVM_RKHS_poly_d1',
                       'SVM_RKHS_poly_d2', 'SVM_RKHS_poly_d3')) |>
  ggplot(aes(x = hyperparameter, #y = after_stat(density),
             fill = method, colour = method)) + 
  geom_histogram(binwidth = 1, alpha = 0.6) + 
  theme_bw() + 
  facet_wrap(~method, scales = 'free', ncol = 3) +
  labs(x = 'Hodnoty hyperparametru',
       y = 'Absolutní počet') + 
  theme(legend.position = 'none')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-243-1.png" alt="Histogramy hodnot hyperparametrů metody RKHS + SVM s polynomiálním jádrem." width="672" />
<p class="caption">(\#fig:unnamed-chunk-243)Histogramy hodnot hyperparametrů metody RKHS + SVM s polynomiálním jádrem.</p>
</div>





```r
CV_res |> 
  filter(method %in% c('SVM_RKHS_linear_d1',
                       'SVM_RKHS_linear_d2', 'SVM_RKHS_linear_d3')) |>
  ggplot(aes(x = hyperparameter, #y = after_stat(density),
             fill = method, colour = method)) + 
  geom_histogram(binwidth = 5, alpha = 0.6) + 
  theme_bw() + 
  facet_wrap(~method, scales = 'free', ncol = 3) +
  labs(x = 'Hodnoty hyperparametru',
       y = 'Absolutní počet') + 
  theme(legend.position = 'none')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-245-1.png" alt="Histogramy hodnot hyperparametrů metody RKHS + SVM s lineárním jádrem." width="672" />
<p class="caption">(\#fig:unnamed-chunk-245)Histogramy hodnot hyperparametrů metody RKHS + SVM s lineárním jádrem.</p>
</div>



### Simulace pro derivovaná data

Nyní se konečně podívejme na simulaci derivovaných dat (určili jsme druhou derivaci křivek). Vzhledem k vysokému počtu uvažovaných klasifikačních metod proveďme simulace zvlášť se stejným nastavením generátoru pseudonáhodných čísel (`set.seed(42)`), abychom mohli výsledky simulací porovnat  mezi sebou.


```r
# nastaveni generatoru pseudonahodnych cisel
set.seed(42)

# pocet simulaci
n.sim <- 100

## list, do ktereho budeme ukladat hodnoty chybovosti
# ve sloupcich budou metody
# v radcich budou jednotliva opakovani
# list ma dve polozky ... train a test
methods <- c('KNN', 'LDA', 'QDA', 'LR_functional', 'LR_score', 'Tree_discr',
             'Tree_score', 'Tree_Bbasis', 'RF_discr', 'RF_score', 'RF_Bbasis', 
             'SVM linear - diskr', 'SVM poly - diskr', 'SVM rbf - diskr', 
             'SVM linear - PCA', 'SVM poly - PCA', 'SVM rbf - PCA', 
             'SVM linear - Bbasis', 'SVM poly - Bbasis', 'SVM rbf - Bbasis',
             'SVM linear - projection', 'SVM poly - projection', 
             'SVM rbf - projection', 'SVM linear - RKHS - radial', 
             'SVM poly - RKHS - radial', 'SVM rbf - RKHS - radial', 
             'SVM linear - RKHS - poly', 'SVM poly - RKHS - poly', 
             'SVM rbf - RKHS - poly', 'SVM linear - RKHS - linear', 
             'SVM poly - RKHS - linear', 'SVM rbf - RKHS - linear')

SIMULACE <- list(train = as.data.frame(matrix(NA, ncol = length(methods), 
                                             nrow = n.sim,
                                             dimnames = list(1:n.sim, methods))), 
                 test = as.data.frame(matrix(NA, ncol = length(methods), 
                                             nrow = n.sim,
                                             dimnames = list(1:n.sim, methods))))

# objekt na ulozeni optimalnich hodnot hyperparametru, ktere se urcuji pomoci CV
CV_RESULTS <- data.frame(KNN_K = rep(NA, n.sim), 
                         nharm = NA, 
                         LR_func_n_basis = NA,
                         SVM_d_Linear = NA,
                         SVM_d_Poly = NA,
                         SVM_d_Radial = NA, 
                         SVM_RKHS_radial_gamma1 = NA,
                         SVM_RKHS_radial_gamma2 = NA,
                         SVM_RKHS_radial_gamma3 = NA,
                         SVM_RKHS_radial_d1 = NA,
                         SVM_RKHS_radial_d2 = NA,
                         SVM_RKHS_radial_d3 = NA,
                         SVM_RKHS_poly_p1 = NA,
                         SVM_RKHS_poly_p2 = NA,
                         SVM_RKHS_poly_p3 = NA,
                         SVM_RKHS_poly_d1 = NA,
                         SVM_RKHS_poly_d2 = NA,
                         SVM_RKHS_poly_d3 = NA,
                         SVM_RKHS_linear_d1 = NA,
                         SVM_RKHS_linear_d2 = NA,
                         SVM_RKHS_linear_d3 = NA)
```

Nyní zopakujeme celou předchozí část 100-krát a hodnoty chybovostí si budeme ukládat to listu `SIMULACE`.
Do datové tabulky `CV_RESULTS` si potom budeme ukládat hodnoty optimálních hyperparametrů -- pro metodu $K$ nejbližších sousedů a pro SVM hodnotu dimenze $d$ v případě projekce na B-splinovou bázi. Také uložíme všechny hodnoty hyperparametrů pro metodu SVM + RKHS.


```r
# nastaveni generatoru pseudonahodnych cisel
set.seed(42)

## SIMULACE

for(sim in 1:n.sim) {
  # rozdeleni na testovaci a trenovaci cast
  split <- sample.split(XXder$fdnames$reps, SplitRatio = 0.7)
  
  # vytvoreni vektoru 0 a 1, 0 pro < 20 a 1 pro > 20 
  Y <- ifelse(labels == 'large', 1, 0)
  
  X.train <- subset(XXder, split == TRUE)
  X.test <- subset(XXder, split == FALSE)
  
  Y.train <- subset(Y, split == TRUE)
  Y.test <- subset(Y, split == FALSE)
  
  x.train <- fdata(X.train)
  y.train <- as.numeric(factor(Y.train))
  
  ## 1) K nejbližších sousedů
  
  k_cv <- 10 # k-fold CV
  neighbours <- c(1:10) # pocet sousedu 
  
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)
  
  # kontrola, ze mame opravdu k = k_cv
  while (length(folds) != k_cv) {
    folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)
  }
  
  CV.results <- matrix(NA, nrow = length(neighbours), ncol = k_cv)
  
  for (index in 1:k_cv) {
    # definujeme danou indexovou mnozinu
    fold <- folds[[index]]
      
    x.train.cv <- subset(X.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
      fdata()
    y.train.cv <- subset(Y.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
      factor() |> as.numeric()
    
    x.test.cv <- subset(X.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
      fdata()
    y.test.cv <- subset(Y.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
      factor() |> as.numeric()
    
    # projdeme kazdou cast ... k-krat zopakujeme
    for(neighbour in neighbours) {
      # model pro konkretni volbu K
      neighb.model <- classif.knn(group = y.train.cv, 
                                fdataobj = x.train.cv, 
                                knn = neighbour) 
      # predikce na validacni casti
      model.neighb.predict <- predict(neighb.model, 
                                      new.fdataobj = x.test.cv)
      # presnost na validacni casti
      presnost <- table(y.test.cv, model.neighb.predict) |> 
        prop.table() |> diag() |> sum()
      
      # presnost vlozime na pozici pro dane K a fold
      CV.results[neighbour, index] <- presnost
    }
  }
  
  # spocitame prumerne presnosti pro jednotliva K pres folds
  CV.results <- apply(CV.results, 1, mean)
  K.opt <- which.max(CV.results)
  CV_RESULTS$KNN_K[sim] <- K.opt
  presnost.opt.cv <- max(CV.results)
  CV.results <- data.frame(K = neighbours, CV = CV.results)
  
  neighb.model <- classif.knn(group = y.train, fdataobj = x.train, knn = K.opt)
  
  # predikce
  model.neighb.predict <- predict(neighb.model, 
                                  new.fdataobj = fdata(X.test))
  
  presnost <- table(as.numeric(factor(Y.test)), model.neighb.predict) |>
    prop.table() |>
    diag() |>
    sum()
  
  RESULTS <- data.frame(model = 'KNN', 
                        Err.train = 1 - neighb.model$max.prob,
                        Err.test = 1 - presnost)
  
  ## 2) Lineární diskriminační analýza
  
  # analyza hlavnich komponent
  data.PCA <- pca.fd(X.train, nharm = 10) # nharm - maximalni pocet HK
  nharm <- which(cumsum(data.PCA$varprop) >= 0.9)[1] # urceni p
  CV_RESULTS$nharm[sim] <- nharm
  if(nharm == 1) nharm <- 2
  
  data.PCA <- pca.fd(X.train, nharm = nharm) 
  data.PCA.train <- as.data.frame(data.PCA$scores) # skore prvnich p HK
  data.PCA.train$Y <- factor(Y.train) # prislusnost do trid
  
  # vypocet skoru testovacich funkci
  scores <- matrix(NA, ncol = nharm, nrow = length(Y.test)) # prazdna matice 
  
  for(k in 1:dim(scores)[1]) {
    xfd = X.test[k] - data.PCA$meanfd[1] # k-te pozorovani - prumerna funkce
    scores[k, ] = inprod(xfd, data.PCA$harmonics) 
    # skalarni soucin rezidua a vlastnich funkci rho (funkcionalni hlavni komponenty)
  }
  
  data.PCA.test <- as.data.frame(scores)
  data.PCA.test$Y <- factor(Y.test)
  colnames(data.PCA.test) <- colnames(data.PCA.train) 
  
  # model
  clf.LDA <- lda(Y ~ ., data = data.PCA.train)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.LDA, newdata = data.PCA.train)
  presnost.train <- table(data.PCA.train$Y, predictions.train$class) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.LDA, newdata = data.PCA.test)
  presnost.test <- table(data.PCA.test$Y, predictions.test$class) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'LDA', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 3) Kvadratická diskriminační analýza
  
  # model
  clf.QDA <- qda(Y ~ ., data = data.PCA.train)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.QDA, newdata = data.PCA.train)
  presnost.train <- table(data.PCA.train$Y, predictions.train$class) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.QDA, newdata = data.PCA.test)
  presnost.test <- table(data.PCA.test$Y, predictions.test$class) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'QDA', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 4) Logistická regrese
  ### 4.1) Funkcionální logistická regrese
  
  # vytvorime vhodne objekty
  x.train <- fdata(X.train)
  y.train <- as.numeric(Y.train)
  
  # body, ve kterych jsou funkce vyhodnoceny
  tt <- x.train[["argvals"]]
  
  dataf <- as.data.frame(y.train) 
  colnames(dataf) <- "Y"
  # B-spline baze 
  nbasis.x <- 7
  basis1 <- create.bspline.basis(rangeval = range(tt), 
                                 nbasis = nbasis.x)
  
  ### 10-fold cross-validation
  n.basis.max <- 25
  n.basis <- 4:n.basis.max
  k_cv <- 10 # k-fold CV
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)
  # kontrola, ze mame opravdu k = k_cv
  while (length(folds) != k_cv) {
    folds <- createMultiFolds(X.train$fdnames$reps, k = k_cv, time = 1)
  }
  ## prvky, ktere se behem cyklu nemeni
  # body, ve kterych jsou funkce vyhodnoceny
  tt <- x.train[["argvals"]]
  rangeval <- range(tt)
  # vztah
  f <- Y ~ x
  # baze pro x
  basis.x <- list("x" = basis1)
  
  CV.results <- matrix(NA, nrow = length(n.basis), ncol = k_cv, 
                       dimnames = list(n.basis, 1:k_cv))
  
  for (index in 1:k_cv) {
    # definujeme danou indexovou mnozinu
    fold <- folds[[index]]
      
    x.train.cv <- subset(X.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
      fdata()
    y.train.cv <- subset(Y.train, c(1:length(X.train$fdnames$reps)) %in% fold) |>
      as.numeric()
    
    x.test.cv <- subset(X.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
      fdata()
    y.test.cv <- subset(Y.train, !c(1:length(X.train$fdnames$reps)) %in% fold) |>
      as.numeric()
    
    dataf <- as.data.frame(y.train.cv) 
    colnames(dataf) <- "Y"
    
    for (i in n.basis) {
      # baze pro bety
      basis2 <- create.bspline.basis(rangeval = rangeval, nbasis = i)
      
      basis.b <- list("x" = basis2)
      # vstupni data do modelu
      ldata <- list("df" = dataf, "x" = x.train.cv)
      # binomicky model ... model logisticke regrese
      model.glm <- fregre.glm(f, family = binomial(), data = ldata,
                              basis.x = basis.x, basis.b = basis.b)
      
      # presnost na validacni casti 
      newldata = list("df" = as.data.frame(y.test.cv), "x" = x.test.cv)
      predictions.valid <- predict(model.glm, newx = newldata)
      predictions.valid <- data.frame(Y.pred = ifelse(predictions.valid < 1/2, 0, 1))
      presnost.valid <- table(y.test.cv, predictions.valid$Y.pred) |>
        prop.table() |> diag() |> sum()
      
      # vlozime do matice
      CV.results[as.character(i), as.character(index)] <- presnost.valid
    } 
  }
  
  # spocitame prumerne presnosti pro jednotliva n pres folds
  CV.results <- apply(CV.results, 1, mean)
  n.basis.opt <- n.basis[which.max(CV.results)]
  CV_RESULTS$LR_func_n_basis[sim] <- n.basis.opt
  presnost.opt.cv <- max(CV.results)
  
  # optimalni model
  basis2 <- create.bspline.basis(rangeval = range(tt), nbasis = n.basis.opt)
  f <- Y ~ x
  # baze pro x a bety
  basis.x <- list("x" = basis1) 
  basis.b <- list("x" = basis2)
  # vstupni data do modelu
  dataf <- as.data.frame(y.train) 
  colnames(dataf) <- "Y"
  ldata <- list("df" = dataf, "x" = x.train)
  # binomicky model ... model logisticke regrese
  model.glm <- fregre.glm(f, family = binomial(), data = ldata,
                          basis.x = basis.x, basis.b = basis.b)
  
  # presnost na trenovacich datech
  predictions.train <- predict(model.glm, newx = ldata)
  predictions.train <- data.frame(Y.pred = ifelse(predictions.train < 1/2, 0, 1))
  presnost.train <- table(Y.train, predictions.train$Y.pred) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  newldata = list("df" = as.data.frame(Y.test), "x" = fdata(X.test))
  predictions.test <- predict(model.glm, newx = newldata)
  predictions.test <- data.frame(Y.pred = ifelse(predictions.test < 1/2, 0, 1))
  presnost.test <- table(Y.test, predictions.test$Y.pred) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'LR_functional', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 4.2) Logistická regrese s analýzou hlavních komponent
  
  # model
  clf.LR <- glm(Y ~  ., data = data.PCA.train, family = binomial)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.LR, newdata = data.PCA.train, type = 'response')
  predictions.train <- ifelse(predictions.train > 0.5, 1, 0)
  presnost.train <- table(data.PCA.train$Y, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.LR, newdata = data.PCA.test, type = 'response')
  predictions.test <- ifelse(predictions.test > 0.5, 1, 0)
  presnost.test <- table(data.PCA.test$Y, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'LR_score', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 5) Rozhodovací stromy
  ### 5.1) Diskretizace intervalu
  
  # posloupnost bodu, ve kterych funkce vyhodnotime
  t.seq <- seq(min(t), max(t), length = 101)
     
  grid.data <- eval.fd(fdobj = X.train, evalarg = t.seq)
  grid.data <- as.data.frame(t(grid.data)) # transpozice kvuli funkcim v radku
  grid.data$Y <- Y.train |> factor()
  
  grid.data.test <- eval.fd(fdobj = X.test, evalarg = t.seq)
  grid.data.test <- as.data.frame(t(grid.data.test))
  grid.data.test$Y <- Y.test |> factor()
  
  # sestrojeni modelu
  clf.tree <- train(Y ~ ., data = grid.data, 
                   method = "rpart", 
                   trControl = trainControl(method = "CV", number = 10),
                   metric = "Accuracy")
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.tree, newdata = grid.data)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.tree, newdata = grid.data.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'Tree_discr', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 5.2) Skóre hlavních komponent
  
  # sestrojeni modelu
  clf.tree.PCA <- train(Y ~ ., data = data.PCA.train, 
                   method = "rpart", 
                   trControl = trainControl(method = "CV", number = 10),
                   metric = "Accuracy")
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.tree.PCA, newdata = data.PCA.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.tree.PCA, newdata = data.PCA.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'Tree_score', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 5.3) Bázové koeficienty
  
  # trenovaci dataset
  data.Bbasis.train <- t(X.train$coefs) |> as.data.frame()
  data.Bbasis.train$Y <- factor(Y.train)
  
  # testovaci dataset
  data.Bbasis.test <- t(X.test$coefs) |> as.data.frame()
  data.Bbasis.test$Y <- factor(Y.test)
  
  # sestrojeni modelu
  clf.tree.Bbasis <- train(Y ~ ., data = data.Bbasis.train, 
                   method = "rpart", 
                   trControl = trainControl(method = "CV", number = 10),
                   metric = "Accuracy")
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.tree.Bbasis, newdata = data.Bbasis.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.tree.Bbasis, newdata = data.Bbasis.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'Tree_Bbasis', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 6) Náhodné lesy
  
  ### 6.1) Diskretizace intervalu
  
  # sestrojeni modelu
  clf.RF <- randomForest(Y ~ ., data = grid.data, 
                         ntree = 500, # pocet stromu
                         importance = TRUE,
                         nodesize = 5)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.RF, newdata = grid.data)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.RF, newdata = grid.data.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'RF_discr', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 6.2) Skóre hlavních komponent
  
  # sestrojeni modelu
  clf.RF.PCA <- randomForest(Y ~ ., data = data.PCA.train, 
                             ntree = 500, # pocet stromu
                             importance = TRUE,
                             nodesize = 5)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.RF.PCA, newdata = data.PCA.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.RF.PCA, newdata = data.PCA.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'RF_score', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 6.3) Bázové koeficienty
  
  # sestrojeni modelu
  clf.RF.Bbasis <- randomForest(Y ~ ., data = data.Bbasis.train, 
                                ntree = 500, # pocet stromu
                                importance = TRUE,
                                nodesize = 5)
  
  # presnost na trenovacich datech
  predictions.train <- predict(clf.RF.Bbasis, newdata = data.Bbasis.train)
  presnost.train <- table(Y.train, predictions.train) |>
    prop.table() |> diag() |> sum()
    
  # presnost na trenovacich datech
  predictions.test <- predict(clf.RF.Bbasis, newdata = data.Bbasis.test)
  presnost.test <- table(Y.test, predictions.test) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = 'RF_Bbasis', 
                    Err.train = 1 - presnost.train,
                    Err.test = 1 - presnost.test)
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 7) SVM
  
  ### 7.1) Diskretizace intervalu
  
  # rozdeleni na testovaci a trenovaci cast
  X.train_norm <- subset(XXfd_norm_der, split == TRUE)
  X.test_norm <- subset(XXfd_norm_der, split == FALSE)
  
  Y.train_norm <- subset(Y, split == TRUE)
  Y.test_norm <- subset(Y, split == FALSE)
  
  grid.data <- eval.fd(fdobj = X.train_norm, evalarg = t.seq)
  grid.data <- as.data.frame(t(grid.data)) 
  grid.data$Y <- Y.train_norm |> factor()
  
  grid.data.test <- eval.fd(fdobj = X.test_norm, evalarg = t.seq)
  grid.data.test <- as.data.frame(t(grid.data.test))
  grid.data.test$Y <- Y.test_norm |> factor()
  
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  # kontrola, ze mame opravdu k = k_cv
  while (length(folds) != k_cv) {
    folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  }
  
  # ktere hodnoty gamma chceme uvazovat
  gamma.cv <- 10^seq(-2, 2, length = 5)
  C.cv <- 10^seq(-3, 2, length = 5)
  p.cv <- 3
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
      cv_sample <- 1:dim(grid.data)[1] %in% fold
      
      data.grid.train.cv <- as.data.frame(grid.data[cv_sample, ])
      data.grid.test.cv <- as.data.frame(grid.data[!cv_sample, ])
      
      ## LINEARNI JADRO
      # sestrojeni modelu
      clf.SVM.l <- svm(Y ~ ., data = data.grid.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       kernel = 'linear')
      
      # presnost na validacnich datech
      predictions.test.l <- predict(clf.SVM.l, newdata = data.grid.test.cv)
      presnost.test.l <- table(data.grid.test.cv$Y, predictions.test.l) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane C a fold
      CV.results$SVM.l[(1:length(C.cv))[C.cv == C], 
                       index_cv] <- presnost.test.l
      
      ## POLYNOMIALNI JADRO
      for (p in p.cv) {
        # sestrojeni modelu
        clf.SVM.p <- svm(Y ~ ., data = data.grid.train.cv,
                         type = 'C-classification',
                         scale = TRUE,
                         cost = C,
                         coef0 = coef0,
                         degree = p,
                         kernel = 'polynomial')
        
        # presnost na validacnich datech
        predictions.test.p <- predict(clf.SVM.p, 
                                      newdata = data.grid.test.cv)
        presnost.test.p <- table(data.grid.test.cv$Y, predictions.test.p) |>
          prop.table() |> diag() |> sum()
        
        # presnosti vlozime na pozice pro dane C, p a fold
        CV.results$SVM.p[(1:length(C.cv))[C.cv == C], 
                         (1:length(p.cv))[p.cv == p],
                         index_cv] <- presnost.test.p
      }
          
      ## RADIALNI JADRO
      for (gamma in gamma.cv) {
        # sestrojeni modelu
        clf.SVM.r <- svm(Y ~ ., data = data.grid.train.cv,
                         type = 'C-classification',
                         scale = TRUE,
                         cost = C,
                         gamma = gamma,
                         kernel = 'radial')
        
        # presnost na validacnich datech
        predictions.test.r <- predict(clf.SVM.r, newdata = data.grid.test.cv)
        presnost.test.r <- table(data.grid.test.cv$Y, predictions.test.r) |>
          prop.table() |> diag() |> sum()
        
        # presnosti vlozime na pozice pro dane C, gamma a fold
        CV.results$SVM.r[(1:length(C.cv))[C.cv == C], 
                         (1:length(gamma.cv))[gamma.cv == gamma],
                         index_cv] <- presnost.test.r
      }
    }
  }
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
  
  # sestrojeni modelu
  clf.SVM.l <- svm(Y ~ ., data = grid.data,
                   type = 'C-classification',
                   scale = TRUE,
                   cost = C.opt[1],
                   kernel = 'linear')
  
  clf.SVM.p <- svm(Y ~ ., data = grid.data,
                   type = 'C-classification',
                   scale = TRUE,
                   cost = C.opt[2],
                   degree = p.opt,
                   coef0 = coef0,
                   kernel = 'polynomial')
  
  clf.SVM.r <- svm(Y ~ ., data = grid.data,
                   type = 'C-classification',
                   scale = TRUE, 
                   cost = C.opt[3],
                   gamma = gamma.opt,
                   kernel = 'radial')
  
  # presnost na trenovacich datech
  predictions.train.l <- predict(clf.SVM.l, newdata = grid.data)
  presnost.train.l <- table(Y.train, predictions.train.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.p <- predict(clf.SVM.p, newdata = grid.data)
  presnost.train.p <- table(Y.train, predictions.train.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.r <- predict(clf.SVM.r, newdata = grid.data)
  presnost.train.r <- table(Y.train, predictions.train.r) |>
    prop.table() |> diag() |> sum()
  
  # presnost na testovacich datech
  predictions.test.l <- predict(clf.SVM.l, newdata = grid.data.test)
  presnost.test.l <- table(Y.test, predictions.test.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.p <- predict(clf.SVM.p, newdata = grid.data.test)
  presnost.test.p <- table(Y.test, predictions.test.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.r <- predict(clf.SVM.r, newdata = grid.data.test)
  presnost.test.r <- table(Y.test, predictions.test.r) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = c('SVM linear - diskr', 
                              'SVM poly - diskr', 
                              'SVM rbf - diskr'), 
                    Err.train = 1 - c(presnost.train.l,
                                      presnost.train.p, presnost.train.r),
                    Err.test = 1 - c(presnost.test.l, 
                                     presnost.test.p, presnost.test.r))
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 7.2) Skóre hlavních komponent
  
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
                       coef0 = coef0,
                       degree = p.opt,
                       kernel = 'polynomial')
  
  clf.SVM.r.PCA <- svm(Y ~ ., data = data.PCA.train,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C.opt[3],
                       gamma = gamma.opt,
                       kernel = 'radial')
  
  # presnost na trenovacich datech
  predictions.train.l <- predict(clf.SVM.l.PCA, newdata = data.PCA.train)
  presnost.train.l <- table(data.PCA.train$Y, predictions.train.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.p <- predict(clf.SVM.p.PCA, newdata = data.PCA.train)
  presnost.train.p <- table(data.PCA.train$Y, predictions.train.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.r <- predict(clf.SVM.r.PCA, newdata = data.PCA.train)
  presnost.train.r <- table(data.PCA.train$Y, predictions.train.r) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test.l <- predict(clf.SVM.l.PCA, newdata = data.PCA.test)
  presnost.test.l <- table(data.PCA.test$Y, predictions.test.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.p <- predict(clf.SVM.p.PCA, newdata = data.PCA.test)
  presnost.test.p <- table(data.PCA.test$Y, predictions.test.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.r <- predict(clf.SVM.r.PCA, newdata = data.PCA.test)
  presnost.test.r <- table(data.PCA.test$Y, predictions.test.r) |>
    prop.table() |> diag() |> sum()
  
  Res <- data.frame(model = c('SVM linear - PCA', 
                              'SVM poly - PCA', 
                              'SVM rbf - PCA'), 
                    Err.train = 1 - c(presnost.train.l, 
                                      presnost.train.p, presnost.train.r),
                    Err.test = 1 - c(presnost.test.l, 
                                     presnost.test.p, presnost.test.r))
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### 7.3) Bázové koeficienty
  
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
      cv_sample <- 1:dim(grid.data)[1] %in% fold
      
      data.Bbasis.train.cv <- as.data.frame(data.Bbasis.train[cv_sample, ])
      data.Bbasis.test.cv <- as.data.frame(data.Bbasis.train[!cv_sample, ])
      
      ## LINEARNI JADRO
      # sestrojeni modelu
      clf.SVM.l <- svm(Y ~ ., data = data.Bbasis.train.cv,
                       type = 'C-classification',
                       scale = TRUE,
                       cost = C,
                       kernel = 'linear')
      
      # presnost na validacnich datech
      predictions.test.l <- predict(clf.SVM.l, newdata = data.Bbasis.test.cv)
      presnost.test.l <- table(data.Bbasis.test.cv$Y, predictions.test.l) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane C a fold
      CV.results$SVM.l[(1:length(C.cv))[C.cv == C], 
                       index_cv] <- presnost.test.l
      
      ## POLYNOMIALNI JADRO
      for (p in p.cv) {
        # sestrojeni modelu
        clf.SVM.p <- svm(Y ~ ., data = data.Bbasis.train.cv,
                         type = 'C-classification',
                         scale = TRUE,
                         cost = C,
                         coef0 = coef0,
                         degree = p,
                         kernel = 'polynomial')
        
        # presnost na validacnich datech
        predictions.test.p <- predict(clf.SVM.p, 
                                      newdata = data.Bbasis.test.cv)
        presnost.test.p <- table(data.Bbasis.test.cv$Y, predictions.test.p) |>
          prop.table() |> diag() |> sum()
        
        # presnosti vlozime na pozice pro dane C, p a fold
        CV.results$SVM.p[(1:length(C.cv))[C.cv == C], 
                         (1:length(p.cv))[p.cv == p],
                         index_cv] <- presnost.test.p
      }
          
      ## RADIALNI JADRO
      for (gamma in gamma.cv) {
        # sestrojeni modelu
        clf.SVM.r <- svm(Y ~ ., data = data.Bbasis.train.cv,
                         type = 'C-classification',
                         scale = TRUE,
                         cost = C,
                         gamma = gamma,
                         kernel = 'radial')
        
        # presnost na validacnich datech
        predictions.test.r <- predict(clf.SVM.r, 
                                      newdata = data.Bbasis.test.cv)
        presnost.test.r <- table(data.Bbasis.test.cv$Y, predictions.test.r) |>
          prop.table() |> diag() |> sum()
        
        # presnosti vlozime na pozice pro dane C, gamma a fold
        CV.results$SVM.r[(1:length(C.cv))[C.cv == C], 
                         (1:length(gamma.cv))[gamma.cv == gamma],
                         index_cv] <- presnost.test.r
      }
    }
  }
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
  # sestrojeni modelu
  clf.SVM.l.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                          type = 'C-classification',
                          scale = TRUE,
                          cost = C.opt[1],
                          kernel = 'linear')
  
  clf.SVM.p.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                          type = 'C-classification',
                          scale = TRUE,
                          cost = C.opt[2],
                          degree = p.opt,
                          coef0 = coef0,
                          kernel = 'polynomial')
  
  clf.SVM.r.Bbasis <- svm(Y ~ ., data = data.Bbasis.train,
                          type = 'C-classification',
                          scale = TRUE,
                          cost = C.opt[3],
                          gamma = gamma.opt,
                          kernel = 'radial')
  
  # presnost na trenovacich datech
  predictions.train.l <- predict(clf.SVM.l.Bbasis, newdata = data.Bbasis.train)
  presnost.train.l <- table(Y.train, predictions.train.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.p <- predict(clf.SVM.p.Bbasis, newdata = data.Bbasis.train)
  presnost.train.p <- table(Y.train, predictions.train.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.train.r <- predict(clf.SVM.r.Bbasis, newdata = data.Bbasis.train)
  presnost.train.r <- table(Y.train, predictions.train.r) |>
    prop.table() |> diag() |> sum()
    
  # presnost na testovacich datech
  predictions.test.l <- predict(clf.SVM.l.Bbasis, newdata = data.Bbasis.test)
  presnost.test.l <- table(Y.test, predictions.test.l) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.p <- predict(clf.SVM.p.Bbasis, newdata = data.Bbasis.test)
  presnost.test.p <- table(Y.test, predictions.test.p) |>
    prop.table() |> diag() |> sum()
  
  predictions.test.r <- predict(clf.SVM.r.Bbasis, newdata = data.Bbasis.test)
  presnost.test.r <- table(Y.test, predictions.test.r) |>
    prop.table() |> diag() |> sum()

  Res <- data.frame(model = c('SVM linear - Bbasis', 
                              'SVM poly - Bbasis', 
                              'SVM rbf - Bbasis'), 
                    Err.train = 1 - c(presnost.train.l, 
                                      presnost.train.p, presnost.train.r),
                    Err.test = 1 - c(presnost.test.l, 
                                     presnost.test.p, presnost.test.r))

  RESULTS <- rbind(RESULTS, Res)
  
  ### 7.4) Projekce na B-splinovou bázi
  
  # hodnoty pro B-splinovou bazi
  rangeval <- range(t)
  norder <- 4
  n_basis_min <- norder
  n_basis_max <- 20
  dimensions <- n_basis_min:n_basis_max 
  folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  
  CV.results <- list(SVM.l = matrix(NA, nrow = length(dimensions), ncol = k_cv),
                     SVM.p = matrix(NA, nrow = length(dimensions), ncol = k_cv),
                     SVM.r = matrix(NA, nrow = length(dimensions), ncol = k_cv))
  
  for (d in dimensions) {
    bbasis <- create.bspline.basis(rangeval = rangeval, 
                                   nbasis = d)
    Projection <- project.basis(y = XX, argvals = t, basisobj = bbasis) 
    XX.train <- subset(t(Projection), split == TRUE)
    
    for (index_cv in 1:k_cv) {
      fold <- folds[[index_cv]]
      cv_sample <- 1:dim(XX.train)[1] %in% fold
      
      data.projection.train.cv <- as.data.frame(XX.train[cv_sample, ])
      data.projection.train.cv$Y <- factor(Y.train[cv_sample])
      data.projection.test.cv <- as.data.frame(XX.train[!cv_sample, ])
      Y.test.cv <- Y.train[!cv_sample]
      data.projection.test.cv$Y <- factor(Y.test.cv)
      # sestrojeni modelu
      clf.SVM.l.projection <- svm(Y ~ ., data = data.projection.train.cv,
                              type = 'C-classification',
                              scale = TRUE,
                              kernel = 'linear')
      
      clf.SVM.p.projection <- svm(Y ~ ., data = data.projection.train.cv,
                              type = 'C-classification',
                              scale = TRUE,
                              coef0 = coef0,
                              kernel = 'polynomial')
      
      clf.SVM.r.projection <- svm(Y ~ ., data = data.projection.train.cv,
                              type = 'C-classification',
                              scale = TRUE,
                              kernel = 'radial')
      # presnost na validacnich datech
      ## linear kernel
      predictions.test.l <- predict(clf.SVM.l.projection,
                                    newdata = data.projection.test.cv)
      presnost.test.l <- table(Y.test.cv, predictions.test.l) |>
        prop.table() |> diag() |> sum()
      ## polynomial kernel
      predictions.test.p <- predict(clf.SVM.p.projection, 
                                    newdata = data.projection.test.cv)
      presnost.test.p <- table(Y.test.cv, predictions.test.p) |>
        prop.table() |> diag() |> sum()
      ## radial kernel
      predictions.test.r <- predict(clf.SVM.r.projection,
                                    newdata = data.projection.test.cv)
      presnost.test.r <- table(Y.test.cv, predictions.test.r) |>
        prop.table() |> diag() |> sum()
      
      # presnosti vlozime na pozice pro dane d a fold
      CV.results$SVM.l[d - min(dimensions) + 1, index_cv] <- presnost.test.l
      CV.results$SVM.p[d - min(dimensions) + 1, index_cv] <- presnost.test.p
      CV.results$SVM.r[d - min(dimensions) + 1, index_cv] <- presnost.test.r
    }
  }
  
  # spocitame prumerne presnosti pro jednotliva d pres folds
  for (n_method in 1:length(CV.results)) {
    CV.results[[n_method]] <- apply(CV.results[[n_method]], 1, mean)
  }
  
  d.opt <- c(which.max(CV.results$SVM.l) + n_basis_min - 1, 
             which.max(CV.results$SVM.p) + n_basis_min - 1, 
             which.max(CV.results$SVM.r) + n_basis_min - 1)
  
  # ulozime optimalni d do datove tabulky
  CV_RESULTS[sim, 4:6] <- d.opt
  
  # pripravime si datovou tabulku pro ulozeni vysledku
  Res <- data.frame(model = c('SVM linear - projection', 
                              'SVM poly - projection', 
                              'SVM rbf - projection'), 
                    Err.train = NA,
                    Err.test = NA)
  
  for (kernel_number in 1:3) {
    kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
    bbasis <- create.bspline.basis(rangeval = rangeval, 
                                   nbasis = d.opt[kernel_number])
    Projection <- project.basis(y = XX, argvals = t, basisobj = bbasis) 
    
    XX.train <- subset(t(Projection), split == TRUE)
    XX.test <- subset(t(Projection), split == FALSE)
    
    data.projection.train <- as.data.frame(XX.train)
    data.projection.train$Y <- factor(Y.train)
    
    data.projection.test <- as.data.frame(XX.test)
    data.projection.test$Y <- factor(Y.test)
    
    # sestrojeni modelu
    clf.SVM.projection <- svm(Y ~ ., data = data.projection.train,
                              type = 'C-classification',
                              scale = TRUE,
                              coef0 = coef0,
                              kernel = kernel_type)
    
    # presnost na trenovacich datech
    predictions.train <- predict(clf.SVM.projection, newdata = data.projection.train)
    presnost.train <- table(Y.train, predictions.train) |>
      prop.table() |> diag() |> sum()
      
    # presnost na trenovacich datech
    predictions.test <- predict(clf.SVM.projection, newdata = data.projection.test)
    presnost.test <- table(Y.test, predictions.test) |>
      prop.table() |> diag() |> sum()
    
    # ulozeni vysledku
    Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
  }
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## 7.5) SVM + RKHS
  
  C <- 1
  eps <- 0.01
  
  ### Gaussovo jadro
  
  # jadro a jadrova matice ... Gaussovske s parametrem gamma
  Gauss.kernel <- function(x, y, gamma) {
    return(exp(-gamma * norm(c(x - y) |> t(), type = 'F')))
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
  
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  
  # hodnoty hyperparametru, ktere budeme prochazet
  dimensions <- seq(2, 30, by = 2) # rozumny rozsah hodnot d
  gamma.cv <- 10^seq(-2, 2, length = 15)
  
  # list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
  # prazdna matice, do ktere vlozime jednotlive vysledky
  # ve sloupcich budou hodnoty presnosti pro dane
  # v radcich budou hodnoty pro danou gamma a vrstvy odpovidaji folds
  dim.names <- list(gamma = paste0('gamma:', round(gamma.cv, 3)),
                    d = paste0('d:', dimensions),
                    CV = paste0('cv:', 1:k_cv))
  
  CV.results <- list(
    SVM.l = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.p = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.r = array(NA, dim = c(length(gamma.cv), length(dimensions), k_cv),
                  dimnames = dim.names))
  
  # samotna CV
  for (gamma in gamma.cv) {
    K <- Kernel.RKHS(t.seq, gamma = gamma)
    Eig <- eigen(K)
    eig.vals <- Eig$values
    eig.vectors <- Eig$vectors
    alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 
    
    # model
    for(i in 1:dim(data.RKHS)[2]) {
      df.svm <- data.frame(x = t.seq,
                           y = data.RKHS[, i])
      svm.RKHS <- svm(y ~ x, data = df.svm, 
                      kernel = 'radial',
                      type = 'eps-regression',
                      epsilon = eps, 
                      cost = C,
                      gamma = gamma)
      alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
    }
    
    # projdeme dimenze
    for(d.RKHS in dimensions) {
      Lambda.RKHS <- matrix(NA, 
                            ncol = dim(data.RKHS)[2], 
                            nrow = d.RKHS) 
      # vypocet reprezentace
      for(l in 1:dim(data.RKHS)[2]) {
        Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                               alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
      }
      # projdeme folds
      for (index_cv in 1:k_cv) {
        # definice testovaci a trenovaci casti pro CV
        fold <- folds[[index_cv]]
        # rozdeleni na trenovaci a validacni data
        XX.train <- Lambda.RKHS[, fold]
        XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
        # pripravime si datovou tabulku pro ulozeni vysledku
        Res <- data.frame(model = c('SVM linear - RKHS', 
                                    'SVM poly - RKHS', 
                                    'SVM rbf - RKHS'), 
                          Err.test = NA)
        # projdeme jednotliva jadra
        for (kernel_number in 1:3) {
          kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
        
          data.RKHS.train <- as.data.frame(t(XX.train))
          data.RKHS.train$Y <- factor(Y.train[fold])
          
          data.RKHS.test <- as.data.frame(t(XX.test))
          data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
          
          # sestrojeni modelu
          clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                              type = 'C-classification',
                              scale = TRUE,
                              cost = C,
                              coef0 = coef0,
                              kernel = kernel_type)
          
          # presnost na validacnich datech
          predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
          presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
            prop.table() |> diag() |> sum()
          
          # ulozeni vysledku
          Res[kernel_number, 2] <- 1 - presnost.test
        }
        # presnosti vlozime na pozice pro dane d, gamma a fold
        CV.results$SVM.l[paste0('gamma:', round(gamma, 3)), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[1, 2]
        CV.results$SVM.p[paste0('gamma:', round(gamma, 3)), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[2, 2]
        CV.results$SVM.r[paste0('gamma:', round(gamma, 3)), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[3, 2]
      }
    }
  }
  
  # spocitame prumerne presnosti pro jednotliva d pres folds
  for (n_method in 1:length(CV.results)) {
    CV.results[[n_method]] <- apply(CV.results[[n_method]], c(1, 2), mean)
  }
  
  gamma.opt <- c(which.min(CV.results$SVM.l) %% length(gamma.cv), 
                 which.min(CV.results$SVM.p) %% length(gamma.cv), 
                 which.min(CV.results$SVM.r) %% length(gamma.cv))
  gamma.opt[gamma.opt == 0] <- length(gamma.cv)
  gamma.opt <- gamma.cv[gamma.opt]
  
  d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
             which.min(t(CV.results$SVM.p)) %% length(dimensions), 
             which.min(t(CV.results$SVM.r)) %% length(dimensions))
  d.opt[d.opt == 0] <- length(dimensions)
  d.opt <- dimensions[d.opt]
  
  err.opt.cv <- c(min(CV.results$SVM.l), 
                       min(CV.results$SVM.p),
                       min(CV.results$SVM.r))
  df.RKHS.res <- data.frame(d = d.opt, gamma = gamma.opt, CV = err.opt.cv,
             Kernel = c('linear', 'polynomial', 'radial') |> factor(),
             row.names = c('linear', 'poly', 'radial'))
  
  CV_RESULTS[sim, 7:9] <- gamma.opt
  CV_RESULTS[sim, 10:12] <- d.opt
  
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  # pridame i testovaci data
  data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
  # pripravime si datovou tabulku pro ulozeni vysledku
  Res <- data.frame(model = c('SVM linear - RKHS - radial', 
                              'SVM poly - RKHS - radial', 
                              'SVM rbf - RKHS - radial'), 
                    Err.train = NA,
                    Err.test = NA)
  
  # projdeme jednotliva jadra
  for (kernel_number in 1:3) {
    # spocitame matici K
    gamma <- gamma.opt[kernel_number] # hodnota gamma pomoci CV
    K <- Kernel.RKHS(t.seq, gamma = gamma)
    
    # urcime vlastni cisla a vektory
    Eig <- eigen(K)
    eig.vals <- Eig$values
    eig.vectors <- Eig$vectors
    # urceni koeficientu alpha z SVM
    alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                         ncol = dim(data.RKHS)[2]) # prazdny objekt
    
    # model
    for(i in 1:dim(data.RKHS)[2]) {
      df.svm <- data.frame(x = t.seq,
                           y = data.RKHS[, i])
      svm.RKHS <- svm(y ~ x, data = df.svm, 
                      kernel = 'radial',
                      type = 'eps-regression',
                      epsilon = eps,
                      cost = C,
                      gamma = gamma)
      # urceni alpha
      alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
    }
    # d
    d.RKHS <- d.opt[kernel_number]
    
    # urceni vektoru lambda
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) # vytvoreni prazdneho objektu
    
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    
    # rozdeleni na trenovaci a testovaci data
    XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
    XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]
  
    kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
  
    data.RKHS.train <- as.data.frame(t(XX.train))
    data.RKHS.train$Y <- factor(Y.train)
    
    data.RKHS.test <- as.data.frame(t(XX.test))
    data.RKHS.test$Y <- factor(Y.test)
    
    # sestrojeni modelu
    clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                        type = 'C-classification',
                        scale = TRUE,
                        cost = C,
                        coef0 = coef0,
                        kernel = kernel_type)
    
    # presnost na trenovacich datech
    predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
    presnost.train <- table(Y.train, predictions.train) |>
      prop.table() |> diag() |> sum()
      
    # presnost na testovacich datech
    predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
    presnost.test <- table(Y.test, predictions.test) |>
      prop.table() |> diag() |> sum()
    
    # ulozeni vysledku
    Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
  }
  
  RESULTS <- rbind(RESULTS, Res)

  ### Polynomialni jadro
  
  # jadro a jadrova matice ... polynomialni s parametrem p
  Poly.kernel <- function(x, y, p) {
    return((1 + x * y)^p)
  }
  
  Kernel.RKHS <- function(x, p) {
    K <- matrix(NA, ncol = length(x), nrow = length(x))
    for(i in 1:nrow(K)) {
      for(j in 1:ncol(K)) {
        K[i, j] <- Poly.kernel(x = x[i], y = x[j], p)
      }
    }
    return(K)
  }
  
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  
  # hodnoty hyperparametru, ktere budeme prochazet
  dimensions <- seq(2, 10, by = 1) # rozumny rozsah hodnot d
  poly.cv <- 2:5
  
  # list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
  # prazdna matice, do ktere vlozime jednotlive vysledky
  # ve sloupcich budou hodnoty presnosti pro dane
  # v radcich budou hodnoty pro dane p a vrstvy odpovidaji folds
  dim.names <- list(p = paste0('p:', poly.cv),
                    d = paste0('d:', dimensions),
                    CV = paste0('cv:', 1:k_cv))
  
  CV.results <- list(
    SVM.l = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.p = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.r = array(NA, dim = c(length(poly.cv), length(dimensions), k_cv),
                  dimnames = dim.names))
  
  # samotna CV
  for (p in poly.cv) {
    K <- Kernel.RKHS(t.seq, p = p)
    Eig <- eigen(K)
    eig.vals <- Eig$values
    eig.vectors <- Eig$vectors
    alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 
    
    # model
    for(i in 1:dim(data.RKHS)[2]) {
      df.svm <- data.frame(x = t.seq,
                           y = data.RKHS[, i])
      svm.RKHS <- svm(y ~ x, data = df.svm, 
                      kernel = 'polynomial',
                      type = 'eps-regression',
                      epsilon = eps,                   
                      coef0 = 1,
                      gamma = 1,
                      cost = C,
                      degree = p)
      alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
    }
    
    # projdeme dimenze
    for(d.RKHS in dimensions) {
      Lambda.RKHS <- matrix(NA, 
                            ncol = dim(data.RKHS)[2], 
                            nrow = d.RKHS) 
      # vypocet reprezentace
      for(l in 1:dim(data.RKHS)[2]) {
        Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                               alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
      }
      # projdeme folds
      for (index_cv in 1:k_cv) {
        # definice testovaci a trenovaci casti pro CV
        fold <- folds[[index_cv]]
        # rozdeleni na trenovaci a validacni data
        XX.train <- Lambda.RKHS[, fold]
        XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
        # pripravime si datovou tabulku pro ulozeni vysledku
        Res <- data.frame(model = c('SVM linear - RKHS', 
                                    'SVM poly - RKHS', 
                                    'SVM rbf - RKHS'), 
                          Err.test = NA)
        # projdeme jednotliva jadra
        for (kernel_number in 1:3) {
          kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
        
          data.RKHS.train <- as.data.frame(t(XX.train))
          data.RKHS.train$Y <- factor(Y.train[fold])
          
          data.RKHS.test <- as.data.frame(t(XX.test))
          data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
          
          # sestrojeni modelu
          clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                              type = 'C-classification',
                              scale = TRUE,
                              epsilon = eps,                   
                              coef0 = 1,
                              gamma = 1,
                              cost = C,
                              kernel = kernel_type)
          
          # presnost na validacnich datech
          predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
          presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
            prop.table() |> diag() |> sum()
          
          # ulozeni vysledku
          Res[kernel_number, 2] <- 1 - presnost.test
        }
        # presnosti vlozime na pozice pro dane d, gamma a fold
        CV.results$SVM.l[paste0('p:', p), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[1, 2]
        CV.results$SVM.p[paste0('p:', p), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[2, 2]
        CV.results$SVM.r[paste0('p:', p), 
                         paste0('d:', d.RKHS), 
                         index_cv] <- Res[3, 2]
      }
    }
  }
  
  # spocitame prumerne presnosti pro jednotliva d pres folds
  for (n_method in 1:length(CV.results)) {
    CV.results[[n_method]] <- apply(CV.results[[n_method]], c(1, 2), mean)
  }
  
  poly.opt <- c(which.min(CV.results$SVM.l) %% length(poly.cv), 
                 which.min(CV.results$SVM.p) %% length(poly.cv), 
                 which.min(CV.results$SVM.r) %% length(poly.cv))
  poly.opt[poly.opt == 0] <- length(poly.cv)
  poly.opt <- poly.cv[poly.opt]
  
  d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
             which.min(t(CV.results$SVM.p)) %% length(dimensions), 
             which.min(t(CV.results$SVM.r)) %% length(dimensions))
  d.opt[d.opt == 0] <- length(dimensions)
  d.opt <- dimensions[d.opt]
  
  err.opt.cv <- c(min(CV.results$SVM.l), 
                       min(CV.results$SVM.p),
                       min(CV.results$SVM.r))
  df.RKHS.res <- data.frame(d = d.opt, p = poly.opt, CV = err.opt.cv,
             Kernel = c('linear', 'polynomial', 'radial') |> factor(),
             row.names = c('linear', 'poly', 'radial'))
  
  CV_RESULTS[sim, 13:15] <- poly.opt
  CV_RESULTS[sim, 16:18] <- d.opt
  
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  # pridame i testovaci data
  data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
  
  # pripravime si datovou tabulku pro ulozeni vysledku
  Res <- data.frame(model = c('SVM linear - RKHS - poly', 
                              'SVM poly - RKHS - poly', 
                              'SVM rbf - RKHS - poly'), 
                    Err.train = NA,
                    Err.test = NA)
  
  # projdeme jednotliva jadra
  for (kernel_number in 1:3) {
    # spocitame matici K
    p <- poly.opt[kernel_number] # hodnota gamma pomoci CV
    K <- Kernel.RKHS(t.seq, p = p)
    
    # urcime vlastni cisla a vektory
    Eig <- eigen(K)
    eig.vals <- Eig$values
    eig.vectors <- Eig$vectors
    # urceni koeficientu alpha z SVM
    alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                         ncol = dim(data.RKHS)[2]) # prazdny objekt
    
    # model
    for(i in 1:dim(data.RKHS)[2]) {
      df.svm <- data.frame(x = t.seq,
                           y = data.RKHS[, i])
      svm.RKHS <- svm(y ~ x, data = df.svm, 
                      kernel = 'polynomial',
                      type = 'eps-regression',
                      epsilon = eps,                   
                      coef0 = 1,
                      gamma = 1,
                      cost = C,
                      degree = p)
      # urceni alpha
      alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
    }
    # d
    d.RKHS <- d.opt[kernel_number]
    
    # urceni vektoru lambda
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) # vytvoreni prazdneho objektu
    
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    
    # rozdeleni na trenovaci a testovaci data
    XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
    XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]
  
    kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
  
    data.RKHS.train <- as.data.frame(t(XX.train))
    data.RKHS.train$Y <- factor(Y.train)
    
    data.RKHS.test <- as.data.frame(t(XX.test))
    data.RKHS.test$Y <- factor(Y.test)
    
    # sestrojeni modelu
    clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                        type = 'C-classification',
                        scale = TRUE,
                        epsilon = eps,                   
                        coef0 = 1,
                        gamma = 1,
                        cost = C,
                        kernel = kernel_type)
    
    # presnost na trenovacich datech
    predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
    presnost.train <- table(Y.train, predictions.train) |>
      prop.table() |> diag() |> sum()
      
    # presnost na testovacich datech
    predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
    presnost.test <- table(Y.test, predictions.test) |>
      prop.table() |> diag() |> sum()
    
    # ulozeni vysledku
    Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
  }
  
  RESULTS <- rbind(RESULTS, Res)
  
  ### Linearni jadro
  
  # jadro a jadrova matice ... polynomialni s parametrem p
  Linear.kernel <- function(x, y) {
    return(x * y)
  }
  
  Kernel.RKHS <- function(x) {
    K <- matrix(NA, ncol = length(x), nrow = length(x))
    for(i in 1:nrow(K)) {
      for(j in 1:ncol(K)) {
        K[i, j] <- Linear.kernel(x = x[i], y = x[j])
      }
    }
    return(K)
  }
  
  # rozdelime trenovaci data na k casti
  folds <- createMultiFolds(1:sum(split), k = k_cv, time = 1)
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  
  # hodnoty hyperparametru, ktere budeme prochazet
  dimensions <- seq(3, 40, by = 2) # rozumny rozsah hodnot d
  
  # list se tremi slozkami ... array pro jednotlive jadra -> linear, poly, radial
  # prazdna matice, do ktere vlozime jednotlive vysledky
  # ve sloupcich budou hodnoty presnosti pro dane d
  # v radcich budou hodnoty pro vrstvy odpovidaji folds
  dim.names <- list(d = paste0('d:', dimensions),
                    CV = paste0('cv:', 1:k_cv))
  
  CV.results <- list(
    SVM.l = array(NA, dim = c(length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.p = array(NA, dim = c(length(dimensions), k_cv),
                  dimnames = dim.names),
    SVM.r = array(NA, dim = c(length(dimensions), k_cv),
                  dimnames = dim.names))
  
  # samotna CV
  K <- Kernel.RKHS(t.seq)
  Eig <- eigen(K)
  eig.vals <- Eig$values
  eig.vectors <- Eig$vectors
  alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1], ncol = dim(data.RKHS)[2]) 
  
  # model
  for(i in 1:dim(data.RKHS)[2]) {
    df.svm <- data.frame(x = t.seq,
                         y = data.RKHS[, i])
    svm.RKHS <- svm(y ~ x, data = df.svm, 
                    kernel = 'linear',
                    type = 'eps-regression',
                    cost = C,
                    epsilon = eps)
    alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs 
  }
  
  # projdeme dimenze
  for(d.RKHS in dimensions) {
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) 
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% 
                             alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    # projdeme folds
    for (index_cv in 1:k_cv) {
      # definice testovaci a trenovaci casti pro CV
      fold <- folds[[index_cv]]
      # rozdeleni na trenovaci a validacni data
      XX.train <- Lambda.RKHS[, fold]
      XX.test <- Lambda.RKHS[, !(1:dim(Lambda.RKHS)[2] %in% fold)]
      # pripravime si datovou tabulku pro ulozeni vysledku
      Res <- data.frame(model = c('SVM linear - RKHS', 
                                  'SVM poly - RKHS', 
                                  'SVM rbf - RKHS'), 
                        Err.test = NA)
      # projdeme jednotliva jadra
      for (kernel_number in 1:3) {
        kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
      
        data.RKHS.train <- as.data.frame(t(XX.train))
        data.RKHS.train$Y <- factor(Y.train[fold])
        
        data.RKHS.test <- as.data.frame(t(XX.test))
        data.RKHS.test$Y <- factor(Y.train[!(1:dim(Lambda.RKHS)[2] %in% fold)])
        
        # sestrojeni modelu
        clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                            type = 'C-classification',
                            scale = TRUE,
                            cost = C,
                            coef0 = coef0,
                            kernel = kernel_type)
        
        # presnost na validacnich datech
        predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
        presnost.test <- table(data.RKHS.test$Y, predictions.test) |>
          prop.table() |> diag() |> sum()
        
        # ulozeni vysledku
        Res[kernel_number, 2] <- 1 - presnost.test
      }
      # presnosti vlozime na pozice pro dane d, gamma a fold
      CV.results$SVM.l[paste0('d:', d.RKHS), 
                       index_cv] <- Res[1, 2]
      CV.results$SVM.p[paste0('d:', d.RKHS), 
                       index_cv] <- Res[2, 2]
      CV.results$SVM.r[paste0('d:', d.RKHS), 
                       index_cv] <- Res[3, 2]
    }
  }
  
  # spocitame prumerne presnosti pro jednotliva d pres folds
  for (n_method in 1:length(CV.results)) {
    CV.results[[n_method]] <- apply(CV.results[[n_method]], 1, mean)
  }
  
  d.opt <- c(which.min(t(CV.results$SVM.l)) %% length(dimensions), 
             which.min(t(CV.results$SVM.p)) %% length(dimensions), 
             which.min(t(CV.results$SVM.r)) %% length(dimensions))
  d.opt[d.opt == 0] <- length(dimensions)
  d.opt <- dimensions[d.opt]
  
  err.opt.cv <- c(min(CV.results$SVM.l), 
                       min(CV.results$SVM.p),
                       min(CV.results$SVM.r))
  df.RKHS.res <- data.frame(d = d.opt, CV = err.opt.cv,
             Kernel = c('linear', 'polynomial', 'radial') |> factor(),
             row.names = c('linear', 'poly', 'radial'))
  
  CV_RESULTS[sim, 19:21] <- d.opt
  
  # odstranime posledni sloupec, ve kterem jsou hodnoty Y
  data.RKHS <- grid.data[, -dim(grid.data)[2]] |> t()
  # pridame i testovaci data
  data.RKHS <- cbind(data.RKHS, grid.data.test[, -dim(grid.data.test)[2]] |> t())
  
  # pripravime si datovou tabulku pro ulozeni vysledku
  Res <- data.frame(model = c('SVM linear - RKHS - linear', 
                              'SVM poly - RKHS - linear', 
                              'SVM rbf - RKHS - linear'), 
                    Err.train = NA,
                    Err.test = NA)
  
  # projdeme jednotliva jadra
  for (kernel_number in 1:3) {
    # spocitame matici K
    K <- Kernel.RKHS(t.seq)
    
    # urcime vlastni cisla a vektory
    Eig <- eigen(K)
    eig.vals <- Eig$values
    eig.vectors <- Eig$vectors
    # urceni koeficientu alpha z SVM
    alpha.RKHS <- matrix(0, nrow = dim(data.RKHS)[1],
                         ncol = dim(data.RKHS)[2]) # prazdny objekt
    
    # model
    for(i in 1:dim(data.RKHS)[2]) {
      df.svm <- data.frame(x = t.seq,
                           y = data.RKHS[, i])
      svm.RKHS <- svm(y ~ x, data = df.svm, 
                      kernel = 'linear',
                      type = 'eps-regression',
                      cost = C,
                      epsilon = eps)
      # urceni alpha
      alpha.RKHS[svm.RKHS$index, i] <- svm.RKHS$coefs # nahrazeni nul koeficienty
    }
    # d
    d.RKHS <- d.opt[kernel_number]
    
    # urceni vektoru lambda
    Lambda.RKHS <- matrix(NA, 
                          ncol = dim(data.RKHS)[2], 
                          nrow = d.RKHS) # vytvoreni prazdneho objektu
    
    # vypocet reprezentace
    for(l in 1:dim(data.RKHS)[2]) {
      Lambda.RKHS[, l] <- (t(eig.vectors[, 1:d.RKHS]) %*% alpha.RKHS[, l]) * eig.vals[1:d.RKHS]
    }
    
    # rozdeleni na trenovaci a testovaci data
    XX.train <- Lambda.RKHS[, 1:dim(grid.data)[1]]
    XX.test <- Lambda.RKHS[, (dim(grid.data)[1] + 1):dim(Lambda.RKHS)[2]]
  
    kernel_type <- c('linear', 'polynomial', 'radial')[kernel_number]
  
    data.RKHS.train <- as.data.frame(t(XX.train))
    data.RKHS.train$Y <- factor(Y.train)
    
    data.RKHS.test <- as.data.frame(t(XX.test))
    data.RKHS.test$Y <- factor(Y.test)
    
    # sestrojeni modelu
    clf.SVM.RKHS <- svm(Y ~ ., data = data.RKHS.train,
                        type = 'C-classification',
                        scale = TRUE,
                        cost = C,
                        coef0 = coef0,
                        kernel = kernel_type)
    
    # presnost na trenovacich datech
    predictions.train <- predict(clf.SVM.RKHS, newdata = data.RKHS.train)
    presnost.train <- table(Y.train, predictions.train) |>
      prop.table() |> diag() |> sum()
      
    # presnost na testovacich datech
    predictions.test <- predict(clf.SVM.RKHS, newdata = data.RKHS.test)
    presnost.test <- table(Y.test, predictions.test) |>
      prop.table() |> diag() |> sum()
    
    # ulozeni vysledku
    Res[kernel_number, c(2, 3)] <- 1 - c(presnost.train, presnost.test)
  }
  
  RESULTS <- rbind(RESULTS, Res)
  
  ## pridame vysledky do objektu SIMULACE
  
  SIMULACE$train[sim, ] <- RESULTS$Err.train
  SIMULACE$test[sim, ] <- RESULTS$Err.test
  
  cat('\r', sim)
}

# ulozime vysledne hodnoty 
save(SIMULACE, CV_RESULTS, file = 'RData/aplikace_03der.RData')
```

Nyní spočítáme průměrné testovací a trénovací chybovosti pro jednotlivé klasifikační metody.


```r
# dame do vysledne tabulky

SIMULACE.df <- data.frame(Err.train = apply(SIMULACE$train, 2, mean),
                          Err.test = apply(SIMULACE$test, 2, mean),
                          SD.train = apply(SIMULACE$train, 2, sd),
                          SD.test = apply(SIMULACE$test, 2, sd))

# ulozime vysledne hodnoty 
save(SIMULACE.df, file = 'RData/aplikace_03der_res.RData')
```

#### Výsledky




Table: (\#tab:unnamed-chunk-251)Souhrnné výsledky použitých metod na simulovaných datech. $\widehat{Err}_{train}$ značí odhad trénovací chybovosti, $\widehat{Err}_{test}$ testovací chybovosti, $\widehat{SD}_{train}$ odhad směrodatné odchylky trénovacích chybovostí a $\widehat{SD}_{test}$ je odhad směrodatné odchylky testovacích chybovostí.

                              $\widehat{Err}_{train}$   $\widehat{Err}_{test}$   $\widehat{SD}_{train}$   $\widehat{SD}_{test}$
---------------------------  ------------------------  -----------------------  -----------------------  ----------------------
KNN                                            0.0147                   0.0212                   0.0066                  0.0170
LDA                                            0.0558                   0.0626                   0.0089                  0.0287
QDA                                            0.0105                   0.0145                   0.0064                  0.0129
LR_functional                                  0.0009                   0.0405                   0.0032                  0.0327
LR_score                                       0.0081                   0.0145                   0.0059                  0.0157
Tree_discr                                     0.0109                   0.0263                   0.0454                  0.0603
Tree_score                                     0.0171                   0.0229                   0.0063                  0.0176
Tree_Bbasis                                    0.0107                   0.0251                   0.0455                  0.0595
RF_discr                                       0.0003                   0.0117                   0.0015                  0.0122
RF_score                                       0.0057                   0.0168                   0.0032                  0.0164
RF_Bbasis                                      0.0001                   0.0089                   0.0007                  0.0103
SVM linear - diskr                             0.0033                   0.0091                   0.0052                  0.0137
SVM poly - diskr                               0.0013                   0.0152                   0.0032                  0.0164
SVM rbf - diskr                                0.0025                   0.0148                   0.0040                  0.0135
SVM linear - PCA                               0.0097                   0.0197                   0.0062                  0.0203
SVM poly - PCA                                 0.0080                   0.0174                   0.0060                  0.0167
SVM rbf - PCA                                  0.0073                   0.0174                   0.0058                  0.0144
SVM linear - Bbasis                            0.0033                   0.0249                   0.0065                  0.0216
SVM poly - Bbasis                              0.0031                   0.0234                   0.0044                  0.0174
SVM rbf - Bbasis                               0.0033                   0.0220                   0.0068                  0.0197
SVM linear - projection                        0.0297                   0.0449                   0.0086                  0.0274
SVM poly - projection                          0.0339                   0.0560                   0.0142                  0.0393
SVM rbf - projection                           0.1454                   0.1954                   0.0306                  0.0605
SVM linear - RKHS - radial                     0.0007                   0.0238                   0.0020                  0.0152
SVM poly - RKHS - radial                       0.0024                   0.0238                   0.0037                  0.0170
SVM rbf - RKHS - radial                        0.0038                   0.0203                   0.0046                  0.0148
SVM linear - RKHS - poly                       0.0142                   0.0386                   0.0071                  0.0215
SVM poly - RKHS - poly                         0.0070                   0.0488                   0.0094                  0.0254
SVM rbf - RKHS - poly                          0.0127                   0.0535                   0.0102                  0.0232
SVM linear - RKHS - linear                     0.0063                   0.0442                   0.0089                  0.0219
SVM poly - RKHS - linear                       0.0035                   0.0397                   0.0058                  0.0250
SVM rbf - RKHS - linear                        0.0061                   0.0426                   0.0075                  0.0233

V tabulce výše jsou uvedeny všechny vypočtené charakteristiky.
Jsou zde uvedeny také směrodatné odchylky, abychom mohli porovnat jakousi stálost či míru variability jednotlivých metod.


```r
wilcox.test(SIMULACE$test[, 'RF_Bbasis'], SIMULACE$test[, 'RF_discr'], alternative = 'less', paired = T)$p.value
```

```
## [1] 0.0005059073
```

```r
wilcox.test(SIMULACE$test[, 'RF_Bbasis'], SIMULACE$test[, 'SVM linear - diskr'], alternative = 't', paired = T)$p.value
```

```
## [1] 0.8449667
```

Nakonec ještě můžeme graficky zobrazit vypočtené hodnoty ze simulace pro jednotlivé klasifikační metody pomocí krabicových diagramů, zvlášť pro testovací a trénovací chybovosti.


```r
# pro trenovaci data
SIMULACE$train |> 
  pivot_longer(cols = methods, names_to = 'method', values_to = 'Err') |>
  mutate(method = factor(method, levels = methods, labels = methods, ordered = TRUE)) |> 
  as.data.frame() |>
  ggplot(aes(x = method, y = Err, fill = method, colour = method, alpha = 0.3)) + 
  geom_boxplot(outlier.colour = "white", outlier.shape = 16, outlier.size = 0, 
               notch = FALSE, colour = 'black') + 
  theme_bw() + 
  labs(x = 'Klasifikační metoda',
       y = expression(widehat(Err)[train])) + 
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 40, hjust = 1)) +
  geom_jitter(position = position_jitter(0.15), alpha = 0.7, size = 1, pch = 21,
              colour = 'black') +
  stat_summary(fun = "mean", geom = "point", shape = '+',
               size = 4, color = "black", alpha = 0.9)+ 
  geom_hline(yintercept = min(SIMULACE.df$Err.train), 
             linetype = 'dashed', colour = 'grey')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-253-1.png" alt="Krabicové diagramy trénovacích chybovostí pro 100 simulací zvlášť pro jednotlivé klasifikační metody. Černými symboly $+$ jsou vyznačeny průměry." width="672" />
<p class="caption">(\#fig:unnamed-chunk-253)Krabicové diagramy trénovacích chybovostí pro 100 simulací zvlášť pro jednotlivé klasifikační metody. Černými symboly $+$ jsou vyznačeny průměry.</p>
</div>




```r
# pro testovaci data
SIMULACE$test |> 
  pivot_longer(cols = methods, names_to = 'method', values_to = 'Err') |>
  mutate(method = factor(method, levels = methods, labels = methods, ordered = TRUE)) |> 
  as.data.frame() |>
  ggplot(aes(x = method, y = Err, fill = method, colour = method, alpha = method)) + 
  geom_hline(yintercept = min(SIMULACE.df$Err.test), 
             linetype = 'dashed', colour = 'grey') +
  geom_boxplot(outlier.colour = "white", outlier.shape = 16, outlier.size = 0, 
               notch = FALSE, colour = 'black') + 
  theme_bw() + 
  labs(x = 'Klasifikační metoda',
       # y = "$\\widehat{\\textnormal{Err}}_{test}$"
       y = expression(widehat(Err)[train])) + 
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  geom_jitter(position = position_jitter(0.15), alpha = 0.6, size = 0.9, pch = 21,
              colour = "black") +
  stat_summary(fun = "mean", geom = "point", shape = '+',
               size = 3, color = "black", alpha = 0.9) #+
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-255-1.png" alt="Krabicové diagramy testovacích chybovostí pro 100 simulací zvlášť pro jednotlivé klasifikační metody. Černými symboly $+$ jsou vyznačeny průměry." width="672" />
<p class="caption">(\#fig:unnamed-chunk-255)Krabicové diagramy testovacích chybovostí pro 100 simulací zvlášť pro jednotlivé klasifikační metody. Černými symboly $+$ jsou vyznačeny průměry.</p>
</div>

```r
  # scale_x_discrete(labels = methods_names) +
  # theme(plot.margin = unit(c(0.5, 0.5, 2, 2), "cm")) +
  # coord_cartesian(ylim = c(0, 0.15)) +
  # scale_fill_manual(values = box_col) +
  # scale_alpha_manual(values = box_alpha)

# ggsave("figures/kap7_tecator_box_test_der.tex", device = tikz, width = 9, height = 7)
```



Nakonec se podívejme, jaké hodnoty hyperparametrů byly nejčastější volbou.


Table: (\#tab:unnamed-chunk-257)Mediány hodnot hyperparametrů pro vybrané metody, u nichž se určoval nějaký hyperparametr pomocí cross-validace.

                          Mediánová hodnota hyperparametru
-----------------------  ---------------------------------
KNN_K                                                  5.0
nharm                                                  2.0
LR_func_n_basis                                        6.0
SVM_d_Linear                                           6.0
SVM_d_Poly                                             6.0
SVM_d_Radial                                           6.0
SVM_RKHS_radial_gamma1                                 0.5
SVM_RKHS_radial_gamma2                                 0.3
SVM_RKHS_radial_gamma3                                 0.3
SVM_RKHS_radial_d1                                    14.0
SVM_RKHS_radial_d2                                    12.0
SVM_RKHS_radial_d3                                    10.0
SVM_RKHS_poly_p1                                       4.0
SVM_RKHS_poly_p2                                       4.0
SVM_RKHS_poly_p3                                       4.0
SVM_RKHS_poly_d1                                       5.0
SVM_RKHS_poly_d2                                       5.5
SVM_RKHS_poly_d3                                       4.0
SVM_RKHS_linear_d1                                    21.0
SVM_RKHS_linear_d2                                    21.0
SVM_RKHS_linear_d3                                    24.0


```r
CV_res <- CV_RESULTS |> 
  pivot_longer(cols = CV_RESULTS |> colnames(), names_to = 'method', values_to = 'hyperparameter') |>
  mutate(method = factor(method, 
                         levels = CV_RESULTS |> colnames(), 
                         labels = CV_RESULTS |> colnames(), ordered = TRUE)) |> 
  as.data.frame() 

CV_res |> 
  filter(method %in% c('KNN_K', 'nharm', 'LR_func_n_basis')) |>
  ggplot(aes(x = hyperparameter, #y = after_stat(density),
             fill = method, colour = method)) + 
  geom_histogram(binwidth = 1, alpha = 0.6) + 
  theme_bw() + 
  facet_grid(~method, scales = 'free') +
  labs(x = 'Hodnoty hyperparametru',
       y = 'Absolutní počet') + 
  theme(legend.position = 'none')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-258-1.png" alt="Histogramy hodnot hyperparametrů pro KNN, funkcionální logistickou regresi a také histogram pro počet hlavních komponent." width="672" />
<p class="caption">(\#fig:unnamed-chunk-258)Histogramy hodnot hyperparametrů pro KNN, funkcionální logistickou regresi a také histogram pro počet hlavních komponent.</p>
</div>




```r
CV_res |> 
  filter(method %in% c('SVM_d_Linear', 'SVM_d_Poly', 'SVM_d_Radial')) |>
  ggplot(aes(x = hyperparameter, #y = after_stat(density),
             fill = method, colour = method)) + 
  geom_histogram(binwidth = 1, alpha = 0.6) + 
  theme_bw() + 
  facet_grid(~method, scales = 'free') +
  labs(x = 'Hodnoty hyperparametru',
       y = 'Absolutní počet') + 
  theme(legend.position = 'none')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-260-1.png" alt="Histogramy hodnot hyperparametrů metody SVM s projekcí na B-splinovou bázi." width="672" />
<p class="caption">(\#fig:unnamed-chunk-260)Histogramy hodnot hyperparametrů metody SVM s projekcí na B-splinovou bázi.</p>
</div>




```r
CV_res |> 
  filter(method %in% c('SVM_RKHS_radial_gamma1', 'SVM_RKHS_radial_gamma2',
                       'SVM_RKHS_radial_gamma3', 'SVM_RKHS_radial_d1', 
                       'SVM_RKHS_radial_d2', 'SVM_RKHS_radial_d3')) |>
  ggplot(aes(x = hyperparameter, #y = after_stat(density),
             fill = method, colour = method)) + 
  geom_histogram(bins = 10, alpha = 0.6) + 
  theme_bw() + 
  facet_wrap(~method, scales = 'free', ncol = 3) +
  labs(x = 'Hodnoty hyperparametru',
       y = 'Absolutní počet') + 
  theme(legend.position = 'none')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-262-1.png" alt="Histogramy hodnot hyperparametrů metody RKHS + SVM s radiálním jádrem." width="672" />
<p class="caption">(\#fig:unnamed-chunk-262)Histogramy hodnot hyperparametrů metody RKHS + SVM s radiálním jádrem.</p>
</div>




```r
CV_res |> 
  filter(method %in% c('SVM_RKHS_poly_p1', 'SVM_RKHS_poly_p2',
                       'SVM_RKHS_poly_p3', 'SVM_RKHS_poly_d1',
                       'SVM_RKHS_poly_d2', 'SVM_RKHS_poly_d3')) |>
  ggplot(aes(x = hyperparameter, #y = after_stat(density),
             fill = method, colour = method)) + 
  geom_histogram(binwidth = 1, alpha = 0.6) + 
  theme_bw() + 
  facet_wrap(~method, scales = 'free', ncol = 3) +
  labs(x = 'Hodnoty hyperparametru',
       y = 'Absolutní počet') + 
  theme(legend.position = 'none')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-264-1.png" alt="Histogramy hodnot hyperparametrů metody RKHS + SVM s polynomiálním jádrem." width="672" />
<p class="caption">(\#fig:unnamed-chunk-264)Histogramy hodnot hyperparametrů metody RKHS + SVM s polynomiálním jádrem.</p>
</div>





```r
CV_res |> 
  filter(method %in% c('SVM_RKHS_linear_d1',
                       'SVM_RKHS_linear_d2', 'SVM_RKHS_linear_d3')) |>
  ggplot(aes(x = hyperparameter, #y = after_stat(density),
             fill = method, colour = method)) + 
  geom_histogram(binwidth = 5, alpha = 0.6) + 
  theme_bw() + 
  facet_wrap(~method, scales = 'free', ncol = 3) +
  labs(x = 'Hodnoty hyperparametru',
       y = 'Absolutní počet') + 
  theme(legend.position = 'none')
```

<div class="figure">
<img src="11-Application_3_files/figure-html/unnamed-chunk-266-1.png" alt="Histogramy hodnot hyperparametrů metody RKHS + SVM s lineárním jádrem." width="672" />
<p class="caption">(\#fig:unnamed-chunk-266)Histogramy hodnot hyperparametrů metody RKHS + SVM s lineárním jádrem.</p>
</div>


