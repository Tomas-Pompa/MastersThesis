## Klasifikace funkcionálních dat

Tento dokument slouží jako podpůrný materiál pro diplomovou práci 

**Metoda podpůrných vektorů pro funkcionální data**,

jejíž oficiální zadání zní následovně:

*Zadání*: V mnoha aplikacích naměřená data reprezentují hodnoty nějaké funkce. Proto je
výhodné, pokud to situace dovolí, pracovat s nimi ve funkcionální podobě, tj. jako s prvky
nekonečně rozměrného prostoru. Práce bude volně navazovat na předchozí studentovu bakalářskou
práci. Cílem práce je zobecnit úvahy metod strojového učení na situaci funkcionálních
dat, dále popsat vlastnosti takových přístupů. Získané výsledky budou demonstrovány na
simulovaných nebo reálných datech.

### Support vector machines pro funkcionální data

Cílem dokumentu bude aplikovat poznatky o metodě podpůrných vektorů (SVM) pro mnohorozměrná data na data funkcionálního typu, tedy nekonečně-rozměrné objekty.
K tomu využijeme jednak převod (redukci) objektů z nekonečné dimenze na objekty konečné dimenze a následným využitím známých postupů a také modifikaci SVM přímo pro funkcionální data, k čemuž využijeme poznatky o Hilbertových prostorech a skalárním součinu.

Dalším cílem bude porovnání jednotlivých metod pro klasifikaci funkcionálních dat na reálných a simulovaných datech. Bylo by dobré vymyslet nějakou zajímavou simulační studii, která bude demonstrovat různá chování uvažovaných metod.

Mezi uvažované klasifikační metody patří:

- [] $K$ nejbližších sousedů (KNN),

- [] logistická regrese (jak obyčejná (LR) tak její funkcionální modifikace (LR_fda)),

- [] lineární (LDA) a kvadratická (QDA) diskriminační analýza,

- [] rozhodovací stromy (DT),

- [] náhodné lesy (RF) a 

- [] Support Vector Machines.

Postupně jednotlivé metody projdeme, nejprve na simulovaných datech, a následně budeme konstruovat metodu podpůrných vektorů pro funkcionální data (SVM_fda).

Základním balíčkem v `R` pro práci s funkcionálními objekty je `fda`. Dalšími užitečnými balíčky budou `MASS`, `e1071`, `fda.usc`, `refund` a další.