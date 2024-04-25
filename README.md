## Klasifikace funkcionálních dat

Tento dokument slouží jako podpůrný materiál pro diplomovou práci:

**Metoda podpůrných vektorů pro funkcionální data**,

jejíž oficiální zadání zní následovně.

*Zadání*: V mnoha aplikacích naměřená data reprezentují hodnoty nějaké funkce. Proto je
výhodné, pokud to situace dovolí, pracovat s nimi ve funkcionální podobě, tj. jako s prvky
nekonečně rozměrného prostoru. Práce bude volně navazovat na předchozí studentovu bakalářskou
práci. Cílem práce je zobecnit úvahy metod strojového učení na situaci funkcionálních
dat, dále popsat vlastnosti takových přístupů. Získané výsledky budou demonstrovány na
simulovaných nebo reálných datech.

### Support vector machines pro funkcionální data

Cílem dokumentu bude aplikovat poznatky o metodě podpůrných vektorů (SVM) pro mnohorozměrná data na data funkcionálního typu, tedy nekonečně-rozměrné objekty.
K tomu využijeme zejména převod (redukci) objektů z nekonečné dimenze na objekty konečné dimenze a následným využitím známých postupů z konečných rozměrů. Ukážeme několik možných přístupů.

Dalším cílem bude porovnání jednotlivých metod pro klasifikaci funkcionálních dat na reálných a simulovaných datech. Zaměříme se primárně na simulovaná data a kromě porovnání metod mezi sebou na základě simulační studie se také podíváme na závislost úspěšnosti klasifikace uvažovaných metod na parametrech, které využíváme při generování (bude nás zajímat rozptyl kolem generujících křivek a také rozptyl vertikálního posunutí). Dále nás také bude zajímat závislost chabovosti klasifikačních metod na diskretizaci intervalu, což je jedna z možností, jak aplikovat konečně-rozměrné metody na funkcionální data.

Mezi uvažované klasifikační metody patří:

  - $K$ nejbližších sousedů (KNN),

  - logistická regrese (jak obyčejná (LR) tak její funkcionální modifikace (LR_fda)),

  - lineární (LDA) a kvadratická (QDA) diskriminační analýza,

  - rozhodovací stromy (DT),

  - náhodné lesy (RF) a 

  - Support Vector Machines: zde budeme uvažovat mnoho variant, všechny z nich jsou přitom postaveny na principu filtrace (redukce dimenze).

Postupně jednotlivé metody projdeme, nejprve na simulovaných datech, a následně budeme konstruovat metodu podpůrných vektorů pro funkcionální data.

Základním balíčkem v `R` pro práci s funkcionálními objekty je `fda`. Dalšími užitečnými balíčky budou `MASS`, `e1071`, `fda.usc`, `refund` a další.

V aplikační části dokumentu se podíváme na tři datové soubory -- `growth`, `phoneme` a `tecator`. Poslední kapitola pak obsahuje zdrojový kód k obrázkům které jsou součástí diplomové práce. Výsledky jsou prezentovány jak graficky, tak číselně, podrobné komentáře lze najít právě v diplomové práci.
